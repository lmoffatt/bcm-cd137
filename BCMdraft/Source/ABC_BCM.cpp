#include <deque>
#include "Include/ABC_BCM.h"
#include "Include/MatrixInverse.h"
#include <limits>
#include <chrono>
#include <random>


std::vector<double>& operator+=(std::vector<double>& x, const std::vector<double>& y)
{
  for (std::size_t i=0; i<x.size(); i++)
    x[i]+=y[i];
  return x;
}

std::vector<double> operator+(const std::vector<double>& x, const std::vector<double>& y)
{
  std::vector<double> z(x);
  z+=y;
  return z;

}



void ABC_BCM::init()
{
  trnas_=setTransformations();
  form_=setRangeFormat();

  externalVariables_=new MultipleVariables(this,"External Variables",externalStateVariables());
  internalVariables_=new MultipleVariables(this,"Internal Variables",internalStateVariables());
  observVariables_=new MultipleVariables(this,"Observed Variables",observedStateVariables());
  paramVariables_=new MultipleVariables(this,"Parameter Variables",parameterVariables());

}



std::vector<double> operator*(const std::vector<double>& x, double h)
{
  std::vector<double> z(x.size());
  for (std::size_t i=0; i<x.size(); i++)
    z[i]=x[i]*h;
  return z;

}


const Transformation* ABC_BCM::getTransform(int index)const
{
  auto  it=trnas_.find(index);
  if (it!=trnas_.end())
    return it->second;
  else
    return nullptr;
}

const RangeFormat* ABC_BCM::getRangeFormat(int index)const
{
  auto  it=this->form_.find(index);
  if (it!=form_.end())
    return it->second;
  else
    return nullptr;
}


const MultipleVariables* ABC_BCM::getInternalStateVariables()const{
  return internalVariables_;
}
const MultipleVariables* ABC_BCM::getParameterVariables()const{
  return paramVariables_;
}
const MultipleVariables* ABC_BCM::getExternalStateVariables()const{
  return externalVariables_;
}
const MultipleVariables* ABC_BCM::getObservedStateVariables()const{
  return observVariables_;
}



std::vector<double> ABC_BCM::RungeKutta4(const std::vector<double>& modelParameters,
                                         const std::vector<double>& externalModelState,
                                         const std::vector<double>& y_,
                                         double h)const
{

  std::vector<double> k1=internalStatesDerivative(modelParameters,
                                                  externalModelState,
                                                  y_)*h;

  std::vector<double> k2=internalStatesDerivative(modelParameters,
                                                  externalModelState,
                                                  y_+k1*0.5)*h;

  std::vector<double> k3=internalStatesDerivative(modelParameters,
                                                  externalModelState,
                                                  y_+k2*0.5)*h;

  std::vector<double> k4=internalStatesDerivative(modelParameters,
                                                  externalModelState,
                                                  y_+k3)*h;

  std::vector<double> delta=(k1+k2*2.0+k3*2.0+k4)*(1.0/6.0);

  return y_+delta;

}

const unsigned long ABC_BCM::seed = std::chrono::high_resolution_clock::now()
    .time_since_epoch().count();

double randNormal(double mean,double std)
{
  return randNormal()*std+mean;
}


double randNormal()
{
 // static auto const seed=std::random_device()();
  static std::mt19937 random_engine(ABC_BCM::seed);
  static std::normal_distribution<> normal;
  return normal(random_engine);

/*

  //Box-Muller method http://en.wikipedia.org/wiki/Normal_distribution#Generating_values_from_normal_distribution
  double U=(1.0*rand())/RAND_MAX;
  double V=(1.0*rand())/RAND_MAX;
  const double PI=3.1415926;

  double r=sqrt(-2*log(U))*cos(2*PI*V);
  return r;
*/
}





std::vector<double> ABC_BCM::yfit(const VariablesValue& par,
                                  const ExperimentSchedule& exper,
                                  double dt)const
{
  std::size_t nData=exper.measuredVariables().size();
  std::vector<double> intState=
      internalStatesStart(par.values(),exper.ExternalState().values());

  std::vector<double> intStateNext;
  double t=0;
  auto intObs=getIntegrMeasures(par,exper);
  auto itI=intObs.begin();
  double nextt=itI->first.first;
  double nextend=itI->first.second;


  //     time of expiration                             sumy
  std::map<double,std::vector<std::pair<const SingleVariable*,double> > > runningMeasures; //stores integrating measures

  std::map<int,std::map<double,double> > mapYfit;  // stores the yfit
  if (nextt<=t)  // process itI
    {
      if (nextend>nextt) // observation that needs integration
        {
          for (const SingleVariable* elem:itI->second)
            runningMeasures[nextend].
                push_back({elem,0.0});  // adds to the runningMeasures,nothing more
        }
      else  // means that it is an instantaneous observation
        {
          for (const SingleVariable* elem:itI->second)
            {

              mapYfit[elem->index()][nextt]=
                  observedVariablesEquations(par.values(),
                                             exper.ExternalState().values(),
                                             intState,
                                             elem->index(),
                                             0.0,
                                             0.0,
                                             false);
            }
        }
      ++itI;
      nextt=itI->first.first;
      nextend=itI->first.second;

    }
  std::size_t i=0;
  while (t<=exper.tmax())
    {
      intStateNext=RungeKutta4(par.values(),
                               exper.ExternalState().values(),
                               intState,
                               dt);

      if (isNaN(intStateNext))
        {
          // std::cerr<<intState<<"\n";

          // std::cerr<<intStateNext<<"\n";
          return std::vector<double>(nData,std::numeric_limits<double>::quiet_NaN());
        }


      intState=intStateNext;
      t=dt*i;
      if ((nextt<=t) &&(itI!=intObs.end()))  // process itI
        {
          if (nextend>nextt) // observation that needs integration
            {
              for (const SingleVariable* elem:itI->second)
                runningMeasures[nextend].
                    push_back({elem,0.0});  // adds to the runningMeasures,nothing more
            }
          else  // means that it is an instantaneous observation
            {
              for (const SingleVariable* elem:itI->second)
                {

                  mapYfit[elem->index()][nextt]=observedVariablesEquations(par.values(),
                                                                           exper.ExternalState().values(),
                                                                           intState,
                                                                           elem->index(),0,0,false);
                }
            }
          ++itI;
          if (itI!=intObs.end())
          {
              nextt=itI->first.first;
              nextend=itI->first.second;
          }

        }

      for (auto& it:runningMeasures)// add to the integrating variables
        {
          bool finalcalc=it.first<=t;
          for (auto& elem:it.second)
            {
              elem.second=observedVariablesEquations(par.values(),
                                                     exper.ExternalState().values(),
                                                     intState,
                                                     elem.first->index(),
                                                     elem.second,dt,finalcalc);
            }
          if (finalcalc)
            {
              for (auto& elem:it.second)
                {
                  mapYfit[elem.first->index()][it.first]=elem.second;
                }
            }
        }
      // check if the integrating varibales had finished

      while (!runningMeasures.empty()&& runningMeasures.begin()->first<=t)
        {

          runningMeasures.erase(runningMeasures.begin());
        }
      ++i;
    }
  // convert from mapY to yfit
  std::vector<double> yfit;
  for (auto it=mapYfit.begin(); it!=mapYfit.end();++it)
    {
      std::map<double,double>* mdd=&it->second;
      for (auto it2=mdd->begin(); it2!=mdd->end(); ++it2)
        {
          yfit.push_back(it2->second);
        }

    }
  return yfit;
}


std::map<std::pair<double, double>, std::vector<const SingleVariable *> >
ABC_BCM::getIntegrMeasures(const VariablesValue& par,
                           const ExperimentSchedule& anExp)const
{
  std::map<std::pair<double,double>,std::vector<const SingleVariable *>  > res;

  const MultipleVariables& vars=anExp.measuredVariables();
  const std::vector<double>& times=anExp.measuredTimes();

  const SingleVariable* currvar=vars.var(0);
  double integrationTime=observedVariableIntegrationTime(par.values(),
                                                         anExp.ExternalState().values(),
                                                         currvar->index());



  for (std::size_t i=0; i<vars.size(); ++i)
    {

      if (currvar!=vars.var(i))
        {
          currvar=vars.var(i);
          integrationTime=observedVariableIntegrationTime(par.values(),
                                                          anExp.ExternalState().values(),
                                                          currvar->index());
        }

      double t=times[i];

      res[{t-integrationTime,t}].push_back(currvar);


    }
  return res;
}

std::vector<double> ABC_BCM::Tyfit(const VariablesValue& par,
                                   const std::vector<const ExperimentDistribution *> &exp,
                                   double dt) const
{
  std::vector<double> Ty;
  for (std::size_t i=0; i<exp.size(); i++)
    {
      const ExperimentDistribution* e=exp[i];
      std::vector<double> yf=
          e->schedule().measuredVariables().transform(yfit(par,e->schedule(),dt));
      Ty.insert(Ty.end(),yf.begin(),yf.end());
    }
  return Ty;
}




std::vector<double> ABC_BCM::yfit(const VariablesValue &par,
                                  const std::vector<const ExperimentDistribution*>& exp,
                                  double dt) const
{
  std::vector<std::vector<double>> y(exp.size());
  for (std::size_t i=0; i<exp.size(); i++)
    {
      std::vector<double> yf=yfit(par,exp[i]->schedule(),dt);
      y[i]=yf;
    }
  std::vector<double> yr;
  for (std::size_t i=0; i<exp.size(); i++)
    {
      yr.insert(yr.end(),y[i].begin(),y[i].end());
    }

  return yr;
}
std::vector<double> ABC_BCM::y(const std::vector<const ExperimentDistribution *> exp) const
{
  std::vector<double> yr;

  for (std::size_t i=0; i<exp.size(); i++)
    {
      auto it=yr.end();
      const ExperimentDistribution* e=exp[i];
      std::vector<double> yp=e->y();
      yr.insert(yr.end(),yp.begin(),yp.end());
    }
  return yr;
}

std::vector<double> ABC_BCM::Ty(const std::vector<const ExperimentDistribution *> exp) const
{
  std::vector<double> Tyr;
  for (std::size_t i=0; i<exp.size(); i++)
    {
      auto it=Tyr.end();
      const ExperimentDistribution* e=exp[i];
      Tyr.insert(it,exp[i]->Ty().begin(),e->Ty().end());
    }
  return Tyr;
}



std::vector<double> ABC_BCM::w(const std::vector<const ExperimentDistribution*> exp) const
{
  std::vector<double> w;
  for (std::size_t i=0; i<exp.size(); i++)
    {
      const ExperimentDistribution* e=exp[i];
      w.insert(w.end(),e->w().begin(),e->w().end());
    }
  return w;
}

std::vector<double> ABC_BCM::yerr(const std::vector<const ExperimentDistribution *> exp) const
{
  std::vector<double> se;
  for (std::size_t i=0; i<exp.size(); i++)
    {
      se.insert(se.end(),exp[i]->yerr().begin(),exp[i]->yerr().end());
    }
  return se;
}

std::vector<const ExperimentDistribution*> ABC_BCM::simulate(
    const VariablesValue& par,
    const std::vector<const ExperimentDistribution*>& exper,
    double dt,
    double measureNoisefactor) const
{
  std::vector<const ExperimentDistribution*> res;

  for (std::size_t i=0; i<exper.size(); i++)
    {
      const ExperimentDistribution* e=exper[i];
      ExperimentDistribution* simE= new ExperimentDistribution(*e);
      if (measureNoisefactor>0)
        {
          std::vector<double> Tyf=
              e->schedule().measuredVariables().transform(yfit(par,e->schedule(),dt));
          std::vector<double> se=e->yerr();

          for (std::size_t j=0; j<Tyf.size(); ++j)
            {
              Tyf[j]=randNormal(Tyf[j],se[j]*measureNoisefactor);
            }
          simE->setTmeanValues(Tyf);
        }
      else
        {
          std::vector<double> yf=yfit(par,e->schedule(),dt);
          simE->setMeanValues(yf);

        }
      res.push_back(simE);

    }
  return res;
}



double ABC_BCM::SSdata(const std::vector<double>& yf,
                       const std::vector<double>& y,
                       const std::vector<double>& w) const{
  double result=0;
  for (std::size_t i=0;i<y.size();++i)
    {
      double e=(yf[i]-y[i]);
      result+=e*e*w[i];
    }
  return result;
}


double ABC_BCM::SSpar(const VariablesValue& par) const{
  double result=0;
  std::vector<double> tpar=par.tvalues();
  std::vector<double> tprior=getPrior()->center().tvalues();
  for (std::size_t i=0;i<par.size();++i)
    {
      double e=(tpar[i]-tprior[i])/getPrior()->tstd()[i];
      result+=e*e;
    }
  return result;


}






void ABC_BCM::push_backExperiment(const std::string& experimentIndentifier,
                                  const std::vector<VarValue>& externalState,
                                  const std::vector<VariableTimeDistribution>& data)
{
  expM_[experimentIndentifier]=new ExperimentDistribution(this,
                                                          experimentIndentifier,
                                                          externalState,
                                                          data);

}

void ABC_BCM::push_backParameters(std::string priorIdentifier,
                                  std::vector<VariableDistribution> paramters)
{
  const MultipleVariables*  variables=getParameterVariables();
  std::vector<double> tmean(variables->size());
  std::vector<double> tstd(variables->size());

  for (std::size_t i=0; i<paramters.size(); ++i)
    {
      VariableDistribution *vd=&paramters[i];
      const RangeFormat* F=getRangeFormat(vd->FormatIndex);
      const SingleVariable* v=variables->var(vd->index);
      auto meanstd=F->getTransformedMoments(v->transfromation(),vd->first,vd->second);
      tmean[vd->index]=meanstd.first;
      tstd[vd->index]=meanstd.second;
    }
  parM_[priorIdentifier]=new VariablesNormalDistribution(variables,tmean,tstd);
}

void ABC_BCM::setPrior(std::string priorIdentifier)
{
  prior_=priorIdentifier;
}





const ExperimentDistribution* ABC_BCM::getExperiment(std::string name)const{
  auto it=expM_.find(name);
  if (it!=expM_.end())
    return it->second;
  else
    return nullptr;

}

std::vector<const ExperimentDistribution*> ABC_BCM::getExperiments()const
{
  std::vector<const ExperimentDistribution*> res;
  for (auto it:expM_)
    res.push_back(it.second);
  return res;
}

std::vector<const ExperimentDistribution*> ABC_BCM::getExperiments(std::vector<std::string> nameList)const
{
  std::vector<const ExperimentDistribution*> res;
  for (auto name:nameList)
    {
      auto it=expM_.find(name);
      if (it!=expM_.end())
        res.push_back(it->second);
    }
  return res;
}







const VariablesNormalDistribution* ABC_BCM::getPrior()const
{
  auto it=parM_.find(prior_);
  if (it!=parM_.end())
    return it->second;
  else
    return nullptr;
}





VariablesNormalDistribution ABC_BCM::getJJwbeta(const VariablesValue &par,
                                                double beta,
                                                const std::vector<const ExperimentDistribution*> &exp,
                                                double dx,
                                                double dt)
{
  std::size_t npar=par.size();



  std::vector<double> wD=w(exp);

  std::vector<double>wPar=pow(getPrior()->tstd(),-2);
  std::size_t ndata=wD.size();




  std::vector<std::vector<double> > JJw(npar,std::vector<double>(npar));

  for (std::size_t i=0; i<npar; ++i)
    for (std::size_t j=0; j<npar; ++j)
      {
        if (i==j)
          JJw[i][j]=wPar[i];
        else
          JJw[i][j]=0;
      }

  if (beta>0)
    {
      std::vector<double> y=this->yfit(par,exp,dt);
      std::vector<std::vector<double> > Jacobian(npar,std::vector<double>(ndata));


      for (std::size_t i=0; i<npar; i++)
        {
          VariablesValue p1(par);
          double tx=p1.getTvalue(i);
          double ddx=dx;
          p1.setTvalue(i,tx+ddx);
          std::vector<double> y1=yfit(p1,exp,dt);
          if (isNaN(y1))
            {
              ddx=-ddx;
              p1.setTvalue(i,tx+ddx);
              y1=yfit(p1,exp,dt);
            }
          while (isNaN(y1))
            {
              ddx/=2;
              p1.setTvalue(i,tx+ddx);
              y1=yfit(p1,exp,dt);
            }


          for (std::size_t j=0; j<ndata; j++)
            {
              Jacobian[i][j]=(y1[j]-y[j])/ddx;
            }
        }
      for (std::size_t i=0; i<npar; ++i)
        for (std::size_t j=0; j<npar; ++j)
          {
            for (std::size_t n=0; n<ndata; ++n)
              {
                JJw[i][j]+=Jacobian[i][n]*Jacobian[j][n]*wD[n]*beta;
              }
          }
    }
  std::vector<std::vector<double> > cov=inv(JJw);
  VariablesNormalDistribution result(par,cov);

  return result;

}



double ABC_BCM::posteriorLogLikelihood(const VariablesValue &p,
                                       double beta,
                                       const std::vector<const ExperimentDistribution*> &exp,
                                       double dt)
{
  const double PI=3.1415926;

  std::vector<double> d=this->y(exp);
  std::vector<double> wD=this->w(exp);

  std::vector<double>pr=getPrior()->center().tvalues();

  std::vector<double>wPar=pow(getPrior()->tstd(),-2);


  double logLikelihood=0;
  double logPrior=0;

  std::size_t ndata=d.size();


  std::vector<double> y=this->yfit(p,exp,dt);



  for (std::size_t i=0; i<ndata;i++)
    {
      logLikelihood+= -0.5*pow(y[i]-d[i],2)*wD[i]+0.5 *log(wD[i]/2/PI);
    };
  for (std::size_t i=0;i<p.size();i++)
    {
      logPrior+= -0.5*pow(p.getTvalue(i)-pr[i],2)*wPar[i]+0.5 *log(wPar[i]/2/PI);
    };

  return beta*logLikelihood+logPrior;
}

double ABC_BCM::logLikelihood(const VariablesValue& p,
                              const std::vector<const ExperimentDistribution *> &exp,
                              double dt)
{
  const double PI=3.1415926;

  std::vector<double> d=this->y(exp);
  std::vector<double> wD=this->w(exp);

  double logL=0;

  std::size_t ndata=d.size();


  std::vector<double> y=this->yfit(p,exp,dt);

  for (std::size_t i=0; i<ndata;i++)
    {
      logL+= -0.5*pow(y[i]-d[i],2)*wD[i]+0.5 *log(wD[i]/2/PI);
    };

  return logL;
}





double ABC_BCM::SumWeighedSquare(const VariablesValue& p,
                                 const std::vector<const ExperimentDistribution*> &exp,
                                 double dt)
{

  std::vector<double> d=this->y(exp);
  std::vector<double> wD=this->w(exp);

  double SSw=0;

  std::size_t ndata=d.size();


  std::vector<double> y=this->yfit(p,exp,dt);

  for (std::size_t i=0; i<ndata;i++)
    {
      SSw+= pow(y[i]-d[i],2)*wD[i];
    };

  return SSw;
}



double ABC_BCM::SumWeighedSquareParameters(const VariablesValue &p)
{

  std::vector<double>pr=getPrior()->center().tvalues();

  std::vector<double>wPar=pow(getPrior()->tstd(),-2);


  double SSpar=0;




  for (std::size_t i=0;i<p.size();i++)
    {
      SSpar+= pow(p.getTvalue(i)-pr[i],2)*wPar[i];
    };

  return SSpar;
}

