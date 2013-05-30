#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include <map>
#include "BayesIteration.h"
#include "MatrixInverse.h"


std::ostream& operator<<(std::ostream& s, const Evidence& e)
{
  s<<"\n  Thermodynamic Integration \t"<<e.logZmelt;
  s<<"\n Harmonic mean logEvidence \t"<<e.HMA;
  s<<"\t Harmonic mean logEvidence coefficient of variation\t"<<e.HMAcv<<"\n";
  s<<"\t Gaussian fit logEvidence\t"<<e.logEvidence;
  s<<"\t Gaussian fit logPostLikelihood \t"<<e.logPostLik<<"\n";
  s<<"\t log det post Covariance\t"<<e.logdetCovPost<<"\n";
  s<<"\t Gaussian fit min square estimate\t"<<e.ss;
  s<<"\t Gaussian fit min square standard deviation\t"<<e.ssStd<<std::endl;


}


std::ostream& BayesIteration::
run(std::ostream& f,
    const std::vector<std::pair <double,Parameters> >& pars) const
{
  return m_->run(f,pars);
}


std::ostream& BayesIteration::
run(std::ostream& f,const std::vector<Parameters >& pars,double quantil0, double quantil1)const
{
  return m_->run(f,pars,quantil0,quantil1);

}



std::vector<double> ABC_data::getDataWeigth()
{
  std::vector<double> w=getDataStandardError();
  for (std::size_t i=0; i<w.size(); i++)
    {
      w[i]=1.0/w[i]/w[i];
    }
  std::vector<double> wW(w);
  return wW;
}

ABC_data::~ABC_data(){
}

ABC_model::~ABC_model(){}

BayesIteration::~BayesIteration(){
}



Parameters BayesIteration::Posterior()const
{
  return posterior_;
}

Parameters BayesIteration::Prior(std::size_t n)const
{
  return priors_[n];

}


BayesIteration::BayesIteration(const ABC_model* m,
                               Parameters prior,
                               const ABC_data* d,
                               const std::string& filename):
  m_(m),
  priors_(1,prior),
  data_(1,d),
  posterior_(),
  numSeeds_(20),
  filename_(filename)
{
  //getPosterior();
}


BayesIteration& BayesIteration::addNewData(ABC_data* d)
{
  data_.push_back(d);
  priors_.push_back(posterior_);
  getPosterior();
}



BayesIteration::BayesIteration(const BayesIteration& other):
  m_(other.m_),
  data_(other.data_),
  priors_(other.priors_),
  posterior_(other.posterior_)
{}


/*   void swap(BayesIteration& one, BayesIteration& other);

    BayesIteration& operator=(const BayesIteration& other);

    BayesIteration();

    ~BayesIteration(){}
*/
// void reset(const SimParameters& sp,const Treatment& tr);



ABC_model& BayesIteration::setData(const ABC_data& experimentalData)
{}


BayesIteration::BayesIteration(){}




std::vector<double> BayesIteration::yfit (const std::vector<double>& param){}
std::vector<double> BayesIteration::yfit(const Parameters& parameters)const
{

  std::vector<double> simulatedResults=m_->yfit(parameters);
  std::vector<double> param=parameters.pMeans();

  simulatedResults.insert(simulatedResults.end(),param.begin(),param.end());
  return simulatedResults;
}


std::vector<double> BayesIteration::getData()const
{
  std::vector<double>data =data_.back()->getData();
  std::vector<double> param=priors_.back().pMeans();

  data.insert(data.end(),param.begin(),param.end());
  std::vector<double> d(data);
  return d;


}

std::vector<double> BayesIteration::getDataStandardError()const
{

  std::vector<double>se =data_.back()->getDataStandardError();
  std::vector<double> param=priors_.back().pStds();

  se.insert(se.end(),param.begin(),param.end());
  std::vector<double> o(se);
  return o;

}



BayesIteration& BayesIteration::getImportancePosterior(const Parameters& startingPoint,
                                                       std::size_t numIterations)
{
  std::vector<double> data=getData();
  std::vector<double> w=getDataWeigth();

  Parameters p=priors_.back();


  LevenbergMarquardt LM(this,
                                  data,
                                  startingPoint,
                                  w,
                                  numIterations);
  LM.optimize();
  std::ofstream f;

  f.open(filename_.c_str(),std::ios_base::app);

  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  f<<"----------The prior is-----------------------------\n";
  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";

  put(f,p);


  /*
    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    f<<"----------Start from the following point------------------\n";
    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    put(f,startingPoint);

*/
  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  f<<"----------Result of Levenberg Marquardt------------------\n";
  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  f<<LM;
  f<<"SS \t"<<LM.SS()<<"\n";
  f<<"Evidence \t"<<LM.getEvidence()<<"\n";
  f<<"Posterior Likelihoo \t"<<LM.getLogPostLik()<<"\n";
  f<<"logDetPriorCov \t"<<LM.logDetPriorCov()<<"\n";
  f<<"logDetPostrCov \t"<<LM.logDetPostCov()<<"\n";
  f<<"logDetPostrStd \t"<<LM.logDetPostStd()<<"\n";
  f<<"SSdata \t"<<LM.SSdata()<<"\n";





  put(f,LM.OptimParameters());
  LM.OptimParameters().putCorrelationCoefficient(f);
  f<<"SS \t"<<LM.SS()<<"\n";



  f<<"Suma de cuadrados"<<LM.SS();
  f<<"\t numero iteraciones"<<LM.numIter()<<"\tNumero evaluaciones"<<LM.numEval()<<"\n";
  f<<"Evidence \t"<<LM.getEvidence()<<"\n";
  f<<"Posterior Likelihoo \t"<<LM.getLogPostLik()<<"\n";
  f<<"logDetPriorCov \t"<<LM.logDetPriorCov()<<"\n";
  f<<"logDetPostrCov \t"<<LM.logDetPostCov()<<"\n";
  f<<"logDetPostrStd \t"<<LM.logDetPostStd()<<"\n";
  f<<"SSdata \t"<<LM.SSdata()<<"\n";

  f<<"SS \t"<<LM.SS()<<std::endl;



  if (LM.SS()<100)
    {
      put(f,LM.OptimParameters());
      f<<"--------------------------------------------------"
         "---------------------------------------------------\n";
      f<<"----------Perform the Importance Sampling-----------\n";
      f<<"--------------------------------------------------"
         "---------------------------------------------------\n";
      f.close();
      Evidence ImpSam0=anaImportanceSampling(LM,numIterations);

      Evidence ImpSam1=anaImportanceSampling(LM,numIterations);
      ;

      f.open(filename_.c_str(),std::ios_base::app);


      f<<"Suma de cuadrados"<<LM.SS();
      f<<"\t numero iteraciones"<<LM.numIter()<<"\tNumero evaluaciones"<<LM.numEval()<<"\n";
      f<<"Evidence \t"<<LM.getEvidence()<<"\n";
      f<<"Importance Sampling sample 0"<<ImpSam0<<std::endl;
      f<<"MCMC direct sample 1"<<ImpSam1<<std::endl;
      f<<"MCMC parameters sample 0\n"<<ImpSam0.p;
      f<<"MCMC parameters sample 10\n"<<ImpSam1.p;
      /*
            f<<"MCMC mean curves sample 0\n";
            run(f,mcmc0);
            f<<"MCMC mean curves sample 1\n";
            run(f,mcmc1);

            f<<"Correlation coefficients 0\n";
            ImpSam0.p.putCorrelationCoefficient(f);

            f<<"Correlation coefficients 1\n";
            dmcmc1.p.putCorrelationCoefficient(f);

*/
    }



  f.close();
  return *this;
}



BayesIteration& BayesIteration::getPosterior(const Parameters& startingPoint,
                                             std::size_t numIterations)
{
  std::vector<double> data=getData();
  std::vector<double> w=getDataWeigth();

  Parameters p=priors_.back();


  LevenbergMarquardt LM(this,
                                  data,
                                  startingPoint,
                                  w,
                                  numIterations);
  LM.optimize();
  std::ofstream f;

  f.open(filename_.c_str(),std::ios_base::app);

  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  f<<"----------The prior is-----------------------------\n";
  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";

  put(f,p);


  /*
    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    f<<"----------Start from the following point------------------\n";
    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    put(f,startingPoint);

*/
  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  f<<"----------Result of Levenberg Marquardt------------------\n";
  f<<"--------------------------------------------------"
     "---------------------------------------------------\n";
  f<<LM;
  f<<"SS \t"<<LM.SS()<<"\n";
  f<<"Evidence \t"<<LM.getEvidence()<<"\n";
  f<<"Posterior Likelihoo \t"<<LM.getLogPostLik()<<"\n";
  f<<"logDetPriorCov \t"<<LM.logDetPriorCov()<<"\n";
  f<<"logDetPostrCov \t"<<LM.logDetPostCov()<<"\n";
  f<<"logDetPostrStd \t"<<LM.logDetPostStd()<<"\n";
  f<<"SSdata \t"<<LM.SSdata()<<"\n";





  put(f,LM.OptimParameters());
  LM.OptimParameters().putCorrelationCoefficient(f);
  f<<"SS \t"<<LM.SS()<<"\n";



  f<<"Suma de cuadrados"<<LM.SS();
  f<<"\t numero iteraciones"<<LM.numIter()<<"\tNumero evaluaciones"<<LM.numEval()<<"\n";
  f<<"Evidence \t"<<LM.getEvidence()<<"\n";
  f<<"Posterior Likelihoo \t"<<LM.getLogPostLik()<<"\n";
  f<<"logDetPriorCov \t"<<LM.logDetPriorCov()<<"\n";
  f<<"logDetPostrCov \t"<<LM.logDetPostCov()<<"\n";
  f<<"logDetPostrStd \t"<<LM.logDetPostStd()<<"\n";
  f<<"SSdata \t"<<LM.SSdata()<<"\n";

  f<<"SS \t"<<LM.SS()<<std::endl;



  if (LM.SS()<100)
    {
      put(f,LM.OptimParameters());
      f<<"--------------------------------------------------"
         "---------------------------------------------------\n";
      f<<"----------Perform the MCM--------------------\n";
      f<<"--------------------------------------------------"
         "---------------------------------------------------\n";
      f.close();
      std::vector<std::pair <double, Parameters> >
          mcmc0=runMCMC(LM,numIterations);

      std::vector<std::pair <double, Parameters> >
          mcmc1=runMCMC(LM,numIterations);


      Evidence dmcmc0=anaDiMCMC(mcmc0,LM);

      Evidence dmcmc1=anaDiMCMC(mcmc1,LM);

      f.open(filename_.c_str(),std::ios_base::app);
      f<<"--------------------------------------------------"
         "---------------------------------------------------\n";
      f<<"----------Result of MCMC-------------------------\n";
      f<<"--------------------------------------------------"
         "---------------------------------------------------\n";


      f<<"Suma de cuadrados"<<LM.SS();
      f<<"\t numero iteraciones"<<LM.numIter()<<"\tNumero evaluaciones"<<LM.numEval()<<"\n";
      f<<"Evidence \t"<<LM.getEvidence()<<"\n";
      f<<"MCMC direct sample 0"<<dmcmc0<<std::endl;
      f<<"MCMC direct sample 1"<<dmcmc1<<std::endl;
      f<<"MCMC parameters sample 0\n"<<dmcmc0.p;
      f<<"MCMC parameters sample 10\n"<<dmcmc1.p;

      f<<"MCMC mean curves sample 0\n";
      run(f,mcmc0);
      f<<"MCMC mean curves sample 1\n";
      run(f,mcmc1);

      f<<"Correlation coefficients 0\n";
      dmcmc0.p.putCorrelationCoefficient(f);

      f<<"Correlation coefficients 1\n";
      dmcmc1.p.putCorrelationCoefficient(f);


    }



  f.close();
  return *this;
}

BayesIteration& BayesIteration::getPosterior(const Parameters& startingPoint,
                                             double factor,
                                             std::size_t numSeeds,
                                             double probParChange)
{
  std::vector<double> data=getData();
  std::vector<double> w=getDataWeigth();

  Parameters p=priors_.back();

  std::size_t numIterations=1000;

  std::ofstream f;

  f.open(filename_.c_str(),std::ios_base::app);


  for (std::size_t i=0;i<numSeeds;i++)
    {

      Parameters seed;
      seed=startingPoint.randomSample(p,factor,probParChange);
      LevenbergMarquardt LM(this,
                                      data,
                                      seed,
                                      w,
                                      numIterations);
      LM.optimize();




      f<<"--------------------------------------------------"
         "---------------------------------------------------\n";
      f<<"----------Result of Levenberg Marquardt------------------\n";
      f<<"--------------------------------------------------"
         "---------------------------------------------------\n";
      if (LM.SS()<200)
        {f<<"--------------------------------------------------"
            "---------------------------------------------------\n";
          f<<"----------Start from the following point------------------\n";
          f<<"--------------------------------------------------"
             "---------------------------------------------------\n";
          put(f,seed);


          f<<LM;
        }
      f<<"Suma de cuadrados"<<LM.SS();
      f<<"\t numero iteraciones"<<LM.numIter()<<"\tNumero evaluaciones"<<LM.numEval()<<"\n";
      f<<"Evidence \t"<<LM.getEvidence()<<"\n";
      f<<"Posterior Likelihoo \t"<<LM.getLogPostLik()<<"\n";
      f<<"logDetPriorCov \t"<<LM.logDetPriorCov()<<"\n";
      f<<"logDetPostrCov \t"<<LM.logDetPostCov()<<"\n";
      f<<"logDetPostrStd \t"<<LM.logDetPostStd()<<"\n";
      f<<"SSdata \t"<<LM.SSdata()<<"\n";

      f<<"SS \t"<<LM.SS()<<std::endl;


      if (LM.SS()<100)
        {
          put(f,LM.OptimParameters());
          f<<"--------------------------------------------------"
             "---------------------------------------------------\n";
          f<<"----------Perform the MCM--------------------\n";
          f<<"--------------------------------------------------"
             "---------------------------------------------------\n";
          f.close();
          std::vector<std::pair <double, Parameters> >
              mcmc0=runMCMC(LM,numIterations);

          std::vector<std::pair <double, Parameters> >
              mcmc1=runMCMC(LM,numIterations);


          Evidence dmcmc0=anaDiMCMC(mcmc0,LM);

          Evidence dmcmc1=anaDiMCMC(mcmc1,LM);

          f.open(filename_.c_str(),std::ios_base::app);
          f<<"--------------------------------------------------"
             "---------------------------------------------------\n";
          f<<"----------Result of MCMC-------------------------\n";
          f<<"--------------------------------------------------"
             "---------------------------------------------------\n";


          f<<"SS \t"<<LM.SS()<<"\n";

          f<<"Suma de cuadrados"<<LM.SS();
          f<<"\t numero iteraciones"<<LM.numIter()<<"\tNumero evaluaciones"<<LM.numEval()<<"\n";
          f<<"Evidence \t"<<LM.getEvidence()<<"\n";
          f<<"MCMC direct sample 0"<<dmcmc0<<std::endl;
          f<<"MCMC direct sample 1"<<dmcmc1<<std::endl;
          f<<"MCMC parameters sample 0\n"<<dmcmc0.p;
          f<<"MCMC parameters sample 10\n"<<dmcmc1.p;

          f<<"MCMC mean curves sample 0\n";
          run(f,mcmc0);
          f<<"MCMC mean curves sample 1\n";
          run(f,mcmc1);



        }
      f.close();

    }
  return *this;

}
  Parameters BayesIteration::getFIMbeta(const Parameters& par, double beta)
  {
    std::size_t npar=par.size();

    double dx=1e-7;



    std::vector<double> logL=partialPosteriorLogLikelihood(par,beta);
    std::size_t ntot=logL.size();

    std::vector<std::vector<double> > Jacobian(npar,std::vector<double>(ntot));
    std::vector<std::vector<double> > FIM(npar,std::vector<double>(npar));

    std::vector<double> gradient(npar,0);



    for (std::size_t i=0; i<npar; i++)
      {
        Parameters p1(par);
        p1[i]+=dx;
        std::vector<double> logL1=partialPosteriorLogLikelihood(p1,beta);
        for (std::size_t j=0; j<ntot; j++)
          {
            Jacobian[i][j]=(logL1[j]-logL[j])/dx;
            gradient[i]+=Jacobian[i][j];
          }
        gradient[i]/=ntot;
        for (std::size_t j=0; j<ntot; j++)
          {
            //    Jacobian[i][j]-=gradient[i];
          }
      }
    FIM=prodTr(Jacobian,Jacobian);
    Parameters result(par);
    std::vector<std::vector<double> > cov=inv(FIM);
    result.setCovariance(FIM);

    return result;

  }

  Parameters BayesIteration::getJJwbeta(const Parameters& par, double beta)
  {
    std::size_t npar=par.size();

    double dx=1e-5;

    const double PI=3.1415926;

    std::vector<double> w=this->getDataWeigth();

    std::size_t ntot=w.size();


    std::size_t ndata=ntot-npar;

    std::vector<double> y=this->yfit(par);

    std::vector<std::vector<double> > Jacobian(npar,std::vector<double>(ntot));
    std::vector<std::vector<double> > JJw(npar,std::vector<double>(npar));

    for (std::size_t i=0; i<npar; i++)
      {
        Parameters p1(par);
        double ddx=dx;
        p1[i]+=ddx;
        std::vector<double> y1=yfit(p1);
        if (isNaN(y1))
          {
            ddx=-ddx;
            p1[i]=par[i]+ddx;
            y1=yfit(p1);
          }
        while (isNaN(y1))
          {
            ddx/=2;
            p1[i]=par[i]+ddx;
            y1=yfit(p1);
          }


        for (std::size_t j=0; j<ntot; j++)
          {
            Jacobian[i][j]=(y1[j]-y[j])/ddx;
          }
      }
    for (std::size_t i=0; i<npar; ++i)
      for (std::size_t j=0; j<npar; ++j)
        {
          JJw[i][j]=0;
          for (std::size_t n=0; n<ndata; ++n)
            {
              JJw[i][j]+=Jacobian[i][n]*Jacobian[j][n]*w[n]*beta;
            }
          for (std::size_t n=ndata; n<ntot; ++n)
            {
              JJw[i][j]+=Jacobian[i][n]*Jacobian[j][n]*w[n];
            }
        }
    Parameters result(par);
    std::vector<std::vector<double> > cov=inv(JJw);
    result.setCovariance(cov);

    return result;

  }



  BayesIteration& BayesIteration::getPosterior()
  {
    std::vector<double> data=getData();
    std::vector<double> w=getDataWeigth();
    std::vector<LevenbergMarquardt> LMs;
    std::vector<Parameters> Ps;


    Parameters p=priors_.back();


    std::size_t factor=2;
    std::size_t numIterations=30;
    std::size_t numSeeds=10;
    std::map<double,Parameters> seeds=getRandomParameters(numSeeds*factor,1);

    std::map<double,Parameters> friuts;

    LevenbergMarquardt LM(this,
                                    data,
                                    p,
                                    w,
                                    numIterations);
    LM.optimize();
    LMs.push_back(LM);
    Ps.push_back(LM.OptimParameters());

    friuts[LM.SS()]=LM.OptimParameters();
    std::ofstream f;

    f.open(filename_.c_str(),std::ios_base::app);

    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    f<<"----------The prior is-----------------------------\n";
    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";

    put(f,p);



    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    f<<"----------Start from the center-------------------\n";
    f<<"--------------------------------------------------"
       "---------------------------------------------------\n";
    f<<LM;
    put(f,Ps.back());
    f.close();


    std::map<double,Parameters>::iterator it=seeds.begin();

    for (std::size_t i=0; i<numSeeds; i++)
      {
        Parameters initParam=(*it).second;
        double ss=(*it).first;
        ++it;
        LevenbergMarquardt LM(this,
                                        data,
                                        initParam,
                                        w,
                                        numIterations);
        LM.optimize();
        friuts[LM.SS()]=LM.OptimParameters();

        LMs.push_back(LM);
        Ps.push_back(LM.OptimParameters());
        f.open(filename_.c_str(),std::ios_base::app);
        f<<"--------------------------------------------------"
           "---------------------------------------------------\n";
        f<<"----------Start from the perisphery-------------------\n";
        f<<"--------------------------------------------------"
           "---------------------------------------------------\n";
        f<<LM;
        put(f,Ps.back());
        f.close();

      }
    double errorFactor=0.2;
    double errorShrinkFactor=3;

    std::size_t numCycle=3;
    for (size_t iCycle=0;
         iCycle<numCycle;
         iCycle++)
      {
        seeds=friuts;
        it=seeds.begin();
        size_t j=0;
        size_t jfactor=5;
        Parameters initParam=(*it).second;
        errorFactor/=errorShrinkFactor;



        for (std::size_t i=0; i<numSeeds; i++)
          {
            ++j;
            if(j==jfactor)
              {
                ++it;
                initParam=(*it).second;
                j=0;
              }
            initParam=initParam.randomSample(errorFactor);

            LevenbergMarquardt LM(this,
                                            data,
                                            initParam,
                                            w,
                                            numIterations);
            LM.optimize();
            friuts[LM.SS()]=LM.OptimParameters();

            LMs.push_back(LM);
            Ps.push_back(LM.OptimParameters());
            f.open(filename_.c_str(),std::ios_base::app);
            f<<"--------------------------------------------------"
               "---------------------------------------------------\n";
            f<<"----------Start from the perisphery-------------------\n";
            f<<"--------------------------------------------------"
               "---------------------------------------------------\n";
            f<<LM;
            put(f,Ps.back());
            f.close();

          }


      }

    return *this;
  }


  std::ostream& BayesIteration::put(std::ostream& s,const Parameters& parameters)const
  {
    m_->put(s,parameters);
    return s;
  }


  void BayesIteration::setFilename(const std::string filename)
  {
    filename_=filename;
  }


  std::pair<double,double> BayesIteration::posteriorLogLikelihood(const Parameters& p)
  {
    const double PI=3.1415926;

    std::vector<double> d=this->getData();
    std::vector<double> w=this->getDataWeigth();

    double logLikelihood=0;
    double logPrior=0;

    std::size_t ntot=d.size();
    std::size_t npar=p.size();

    std::size_t ndata=ntot-npar;

    std::vector<double> y=this->yfit(p);



    for (std::size_t i=0; i<ndata;i++)
      {
        logLikelihood+= -0.5*pow(y[i]-d[i],2)*w[i]+0.5 *log(w[i]/2/PI);
      };
    for (std::size_t i=ndata;i<ntot;i++)
      {
        logPrior+= -0.5*pow(y[i]-d[i],2)*w[i]+0.5 *log(w[i]/2/PI);
      };

    return std::pair<double,double>(logLikelihood,logPrior);
  }


  double BayesIteration::logLikelihood(const Parameters& p)

  {
    const double PI=3.1415926;

    std::vector<double> d=this->getData();
    std::vector<double> w=this->getDataWeigth();

    double logLikelihoodVal=0;

    std::size_t ntot=d.size();
    std::size_t npar=p.size();

    std::size_t ndata=ntot-npar;

    std::vector<double> y=this->yfit(p);



    for (std::size_t i=0; i<ndata;i++)
      {
        logLikelihoodVal+= -0.5*pow(y[i]-d[i],2)*w[i]+0.5 *log(w[i]/2/PI);
      };

    return logLikelihoodVal;
  }



  double BayesIteration::posteriorLogLikelihood(const Parameters& p,double beta)
  {
    const double PI=3.1415926;

    std::vector<double> d=this->getData();
    std::vector<double> w=this->getDataWeigth();

    double logLikelihood=0;
    double logPrior=0;

    std::size_t ntot=d.size();
    std::size_t npar=p.size();

    std::size_t ndata=ntot-npar;

    std::vector<double> y=this->yfit(p);



    for (std::size_t i=0; i<ndata;i++)
      {
        logLikelihood+= -0.5*pow(y[i]-d[i],2)*w[i]+0.5 *log(w[i]/2/PI);
      };
    for (std::size_t i=ndata;i<ntot;i++)
      {
        logPrior+= -0.5*pow(y[i]-d[i],2)*w[i]+0.5 *log(w[i]/2/PI);
      };

    return beta*logLikelihood+logPrior;
  }



  std::vector<double> BayesIteration::partialPosteriorLogLikelihood(const Parameters& p, double beta)
  {
    const double PI=3.1415926;

    std::vector<double> d=this->getData();
    std::vector<double> w=this->getDataWeigth();

    std::size_t ntot=d.size();
    std::size_t npar=p.size();

    std::vector<double> logresult(ntot,0);

    std::size_t ndata=ntot-npar;

    std::vector<double> y=this->yfit(p);


    for (std::size_t i=0; i<ndata;i++)
      {
        logresult[i]= beta*(-0.5*pow(y[i]-d[i],2)*w[i]+0.5 *log(w[i]/2.0/PI));
      };
    for (std::size_t i=ndata;i<ntot;i++)
      {
        logresult[i]= -0.5*pow(y[i]-d[i],2)*w[i]+0.5 *log(w[i]/2.0/PI);
      };

    return logresult;
  }


  double BayesIteration::SumWeighedSquare(const Parameters& p)
  {
    std::vector<double> d=this->getData();
    std::vector<double> w=this->getDataWeigth();

    std::vector<double> y=this->yfit(p);
    double SSW=0;
    for (std::size_t i=0;i<d.size();i++)
      {
        SSW+=pow(y[i]-d[i],2)*w[i];
      };
    return SSW;
  }
  double BayesIteration::SumWeighedSquareParameters(const Parameters& p)
  {
    std::vector<double> d=this->getData();
    std::vector<double> w=this->getDataWeigth();

    double SSW=0;
    std::size_t ntot=d.size();
    std::size_t npar=p.size();

    for (std::size_t i=ntot-npar;i<ntot;i++)
      {
        SSW+=pow(p[i+npar-ntot]-d[i],2)*w[i];
      };
    return SSW;
  }


  double BayesIteration::logDetCovData()
  {

    std::vector<double> w=this->getDataWeigth();

    double logdetcov=0;
    std::size_t ntot=w.size();
    std::size_t npar=priors_.back().size();

    for (std::size_t i=0;i<ntot-npar;i++)
      {
        logdetcov-=log(w[i]);
      };
    return logdetcov;
  }

  Parameters BayesIteration::getHessian(const Parameters& MAP,double eps)
  {
    std::vector<std::vector<double> > Hessian(MAP.size(),std::vector<double> (MAP.size(),0));
    double ss=0.5*SumWeighedSquare(MAP);

    for (std::size_t i=0; i<MAP.size(); i++)
      {
        Parameters fip(MAP);
        fip[i]=fip[i]+eps;
        Parameters fin(MAP);
        fin[i]=fin[i]-eps;
        double ssip=0.5*SumWeighedSquare(fip);
        double ssin=0.5*SumWeighedSquare(fin);

        Hessian[i][i]=(ssip+ssin-2.0*ss)/eps/eps;
        for (std::size_t j=i+1; j<MAP.size(); j++)
          {
            Parameters fipjp(MAP);
            fipjp[i]=fipjp[i]+eps;
            Parameters fipjn(fipjp);
            fipjp[j]=fipjp[j]+eps;
            fipjn[j]=fipjn[j]-eps;

            Parameters finjp(MAP);
            finjp[i]=finjp[i]-eps;
            Parameters finjn(finjp);
            finjp[j]=finjp[j]+eps;
            finjn[j]=finjn[j]-eps;

            double ssipjp=0.5*SumWeighedSquare(fipjp);
            double ssinjp=0.5*SumWeighedSquare(finjp);
            double ssipjn=0.5*SumWeighedSquare(fipjn);
            double ssinjn=0.5*SumWeighedSquare(finjn);

            Hessian[i][j]=(ssipjp+ssinjn-ssipjn-ssinjp)/eps/eps/4.0;
            Hessian[j][i]=Hessian[i][j];
          }

      }
    Parameters result(MAP);
    result.setCovariance(inv(Hessian));
    return result;
  }


  Parameters BayesIteration::getHessianInterpol(const Parameters& MAP, double mindSS, double maxdSS)
  {
    std::vector<std::vector<double> > Hessian(MAP.size(),std::vector<double> (MAP.size(),0));
    double ss=0.5*SumWeighedSquare(MAP);

    for (std::size_t i=0; i<MAP.size(); i++)
      {
        double ep=1e-3;
        double h;

        Parameters fip(MAP);
        Parameters fin(MAP);
        double ssip;
        double ssin;

        fip[i]=MAP[i]+ep;
        fin[i]=MAP[i]-ep;
        ssip=0.5*SumWeighedSquare(fip);
        ssin=0.5*SumWeighedSquare(fin);
        h=(ssip+ssin-2.0*ss);
        while (h<mindSS)
          {

          }




        Hessian[i][i]=h;
        for (std::size_t j=i+1; j<MAP.size(); j++)
          {
            Parameters fipjp(MAP);
            fipjp[i]=fipjp[i]+ep;
            Parameters fipjn(fipjp);
            fipjp[j]=fipjp[j]+ep;
            fipjn[j]=fipjn[j]-ep;

            Parameters finjp(MAP);
            finjp[i]=finjp[i]-ep;
            Parameters finjn(finjp);
            finjp[j]=finjp[j]+ep;
            finjn[j]=finjn[j]-ep;

            double ssipjp=0.5*SumWeighedSquare(fipjp);
            double ssinjp=0.5*SumWeighedSquare(finjp);
            double ssipjn=0.5*SumWeighedSquare(fipjn);
            double ssinjn=0.5*SumWeighedSquare(finjn);

            Hessian[i][j]=(ssipjp+ssinjn-ssipjn-ssinjp)/ep/ep/4.0;
            Hessian[j][i]=Hessian[i][j];
          }

      }
    Parameters result(MAP);
    result.setCovariance(inv(Hessian));
    return result;
  }


  std::pair<double,Parameters> BayesIteration::anaLsMCMC(std::vector<std::pair<double,Parameters> > mcmc)
  {
    // fit the data to a parabola
    /*
     ss=ss0+(x-m)Cov(x-m)
     ss(k)=ss0+0.5*m(i)C(ij)m(j) -m(i)C(ij) xk(j)+0.5*C(ij) xk(i)xk(j)

     Xb =y
     b= inv(XT X) XT y

     X=[x(i)x(j),x(i),1]
      Xy
*/
    std::size_t npar=mcmc[0].second.size();
    std::size_t nruns=mcmc.size();

    std::size_t nspar=2*npar+(npar*(npar-1))/2+1;

    std::vector< std::vector <double> > X(mcmc.size(),std::vector<double> (nspar));

    std::vector<std::vector< double > > y(nruns,std::vector<double> (1));

    std::vector<std::vector< double > > b(nspar,std::vector<double> (1));
    for (std::size_t i=0; i< nruns; i++)
      {
        y[i][0]=mcmc[i].first;
        std::size_t jj=0;
        X[i][jj]=1;
        for (std::size_t j=0;j<npar;j++)
          {
            jj++;
            X[i][jj]=mcmc[i].second[j];
          }
        for (std::size_t j=0;j<npar;j++)
          {
            jj++;
            X[i][jj]=0.5*mcmc[i].second[j]*mcmc[i].second[j];

            for (std::size_t k=j+1; k<npar;k++)
              {
                jj++;
                X[i][jj]=mcmc[i].second[j]*mcmc[i].second[k];
              }
          }
      }
    std::vector< std::vector <double> > SS=trProd(X,X);
    std::vector< std::vector <double> > Sxy=trProd(X,y);

    b=prod(inv(SS),Sxy);

    std::vector< std::vector<double> > mH(npar,std::vector<double> (1));
    std::vector<std::vector<double> > H(npar,std::vector<double>(npar,0));
    double ss0mCm;
    std::size_t jj=0;
    ss0mCm=b[jj][0];
    for (std::size_t j=0;j<npar;j++)
      {
        jj++;
        mH[j][0]=b[jj][0];
      }
    for (std::size_t j=0;j<npar;j++)
      {
        jj++;
        H[j][j]=b[jj][0];

        for (std::size_t k=j+1; k<npar;k++)
          {
            jj++;
            H[k][j]=b[jj][0];
            H[j][k]=b[jj][0];
          }
      }
    std::vector< std::vector<double> > C=inv(H);
    std::vector< std::vector<double> > m=prod(C,mH);

    double ss0=ss0mCm-0.5*trProd(m,mH)[0][0];
    std::pair<double,Parameters> result(ss0,mcmc[0].second);
    result.second.setpMeans(m);
    result.second.setCovariance(H);
    return result;


  }


  Evidence
      BayesIteration::anaDiMCMC(std::vector<std::pair<double,Parameters> > mcmc,
                                const LevenbergMarquardt& LM)
  {
    const double PI=3.1415926;

    std::size_t npar=mcmc[0].second.size();
    std::size_t ndata=getData().size();
    std::size_t nruns=mcmc.size();
    std::vector<double> m(npar,0);
    std::vector<std::vector<double> > C(npar,std::vector<double>(npar,0));


    double HMA=0;
    double HMAcv=0;

    for(std::size_t im=0; im<mcmc.size(); ++im)
      {
        double lik=-0.5*mcmc[im].first+0.5*SumWeighedSquareParameters(mcmc[im].second)
            -0.5*logDetCovData()-0.5*(ndata-npar)*log(2*PI);
        HMA+=exp(-lik);
        HMAcv+=exp(-2*lik);

      }
    HMA/=nruns;
    HMAcv=std::sqrt(HMAcv/nruns-HMA*HMA);
    HMAcv=HMAcv/HMA;

    HMA=-log(HMA);
    for(std::size_t im=0; im<mcmc.size(); ++im)
      {
        for (std::size_t i=0; i<npar; i++)
          m[i]+=mcmc[im].second[i]/nruns;
      }
    for(std::size_t im=0; im<nruns; ++im)
      {
        for (std::size_t i=0; i<npar; i++)
          for (std::size_t j=0; j<npar; j++)
            C[i][j]+=(mcmc[im].second[i]-m[i])*
                (mcmc[im].second[j]-m[j])/
                nruns;
      }

    std::vector< std::vector<double> > H=inv(C);

    double ss0=0;
    double ssstd=0;
    for (std::size_t ii=0; ii< nruns; ii++)
      {
        double ss=mcmc[ii].first;
        for (std::size_t i=0;i<npar;i++)
          {
            ss-=pow(mcmc[ii].second[i]-m[i],2)*H[i][i];
            for (std::size_t j=i+1;j<npar;j++)
              {
                ss-=2*(mcmc[ii].second[i]-m[i])*(mcmc[ii].second[j]-m[i])*H[i][j];

              }
          }
        ss0+=ss;
        ssstd+=ss*ss;
      }
    ss0/=nruns;
    ssstd=std::sqrt(ssstd/nruns-ss0*ss0);


    Evidence result;
    result.p=mcmc[0].second;
    result.p.setpMeans(m);
    result.p.setCovariance(C);
    result.ss=ss0;
    result.ssStd=ssstd;
    result.HMA=HMA;
    result.HMAcv=HMAcv;

    result.logdetCovPost=log(det(C));

    result.logEvidence=-0.5*ss0+0.5*result.logdetCovPost-0.5*logDetCovData()-
        0.5*LM.logDetPriorCov()-0.5*(ndata-npar)*log(2*PI);

    result.logPostLik=-0.5*ss0-0.5*LM.logDetPriorCov()
        -0.5*logDetCovData()-0.5*ndata*log(2*PI);


    return result;
  }




  std::vector<std::pair<double, Parameters> >
      BayesIteration::runMCMC(LevenbergMarquardt& LM,
                              std::size_t num)


  {
    Parameters maximumPostLik=LM.OptimParameters();
    double ss0=SumWeighedSquare(maximumPostLik);
    std::vector<std::pair <double,Parameters> > parvec;

    std::size_t npar=maximumPostLik.size();

    std::ofstream f;

    f.open(filename_.c_str(),std::ios_base::app);

    f<<"------------------------------------------------------\n";
    f<<"-----------Markov chain Monte Carlo-------------------\n";
    f<<"------------------------------------------------------\n";


    // lets find the right factor
    std::size_t numsteps=500;

    bool priorMode=false;
    bool discardCOV=false;

    double factor=1.0/maximumPostLik.size();
    Parameters par;
    if (priorMode)
      par=maximumPostLik.randomSample(this->Prior(),factor);
    else
      par=maximumPostLik.randomSample(factor,discardCOV);

    double ss=SumWeighedSquare(par)-ss0;
    if (ss!=ss)
      {
        LM.reOptimize();
        maximumPostLik=LM.OptimParameters();
        ss0=SumWeighedSquare(maximumPostLik);
        if (priorMode)
          par=maximumPostLik.randomSample(this->Prior(),factor);
        else
          par=maximumPostLik.randomSample(factor,discardCOV);

        ss=SumWeighedSquare(par)-ss0;

      }
    while ((ss)>4)
      {
        factor/=sqrt(2);
        if (priorMode)
          par=maximumPostLik.randomSample(this->Prior(),factor);
        else
          par=maximumPostLik.randomSample(factor,discardCOV);
        ss=SumWeighedSquare(par)-ss0;
      }
    while ((ss)<0.5)
      {
        factor*=sqrt(2);
        if (priorMode)
          par=maximumPostLik.randomSample(this->Prior(),factor);
        else
          par=maximumPostLik.randomSample(factor,discardCOV);
        ss=SumWeighedSquare(par)-ss0;
      }
    double sump=0;
    for (std::size_t i=0; i<numsteps; i++)
      {
        if (priorMode)
          par=maximumPostLik.randomSample(this->Prior(),factor);
        else
          par=maximumPostLik.randomSample(factor,discardCOV);
        sump+=exp(-0.5*(SumWeighedSquare(par)-ss0));
      }
    double logmeanp=-2*log(sump/numsteps);
    const double dssN=-2*log(0.23);
    factor=std::sqrt(dssN/logmeanp)*factor;
    sump=0;

    for (std::size_t i=0; i<numsteps; i++)
      {
        if (priorMode)
          par=maximumPostLik.randomSample(this->Prior(),factor);
        else
          par=maximumPostLik.randomSample(factor,discardCOV);
        ss=SumWeighedSquare(par);
        sump+=exp(-0.5*(ss-ss0));
      }
    logmeanp=-2*log(sump/numsteps);
    factor=std::sqrt(dssN/logmeanp)*factor;
    sump=0;

    for (std::size_t i=0; i<numsteps; i++)
      {
        if (priorMode)
          par=maximumPostLik.randomSample(this->Prior(),factor);
        else
          par=maximumPostLik.randomSample(factor,discardCOV);
        ss=SumWeighedSquare(par);
        sump+=exp(-0.5*(ss-ss0));

      }
    logmeanp=-2*log(sump/numsteps);

    f<<"valor medio de ss diference ="<<logmeanp<<"\tfactor= "<<factor<<"\n";
    std::cerr<<"valor medio de ss diference ="<<logmeanp<<"\tfactor= "<<factor<<"\n";

    double lik=exp(-0.5*(ss-ss0));

    double sspar,ssdata;

    for (std::size_t i=0; i<num*1.1; i++ )
      {

        std::size_t j=0;
        std::size_t nsubsteps=par.size();
        for (std::size_t ii=0; ii<par.size(); ii++)
          {
            Parameters nextPar;
            if (priorMode)
              nextPar=par.randomSample(this->Prior(),factor);
            else
              nextPar=par.randomSample(factor,discardCOV);

            double  nextss=SumWeighedSquare(nextPar);
            double nextsspar=SumWeighedSquareParameters(nextPar);

            double nextLik=exp(-0.5*(nextss-ss0));

            double r=(1.0*rand())/RAND_MAX;

            if (nextLik/lik>r)
              {
                par=nextPar;
                ss=nextss;
                sspar=nextsspar;

                lik=nextLik;
                j++;
              }
          }
        double ratio=(1.0*j)/nsubsteps;
        if ((i/20)*20==i)
          f<<"ind \tsaltos ef\t ratio\tfactor\tss\tss param\tss data"<<std::endl;
        f<<i<<"\t"<<j<<"\t"<<ratio<<"\t"<<factor<<"\t"<<ss<<"\t"<<sspar<<"\t"<<ss-sspar<<std::endl;
        if ((i/20)*20==i)
          std::cerr<<"ind \tsaltos ef\t ratio\tfactor\tss\tss param\tss data"<<std::endl;
        std::cerr<<i<<"\t"<<j<<"\t"<<ratio<<"\t"<<factor<<"\t"<<ss<<"\t"<<sspar<<"\t"<<ss-sspar<<std::endl;

        if (ratio<0.10)
          factor/=sqrt(2);
        else if (ratio>0.50)
          factor*=sqrt(2);

        if (i>num*0.1)
          {
            parvec.push_back(std::pair<double,Parameters>(ss,par));

          }


      }

    return parvec;
  }



  MCMCrun
      BayesIteration::initMCMC(const Parameters& start,
                               double beta,
                               std::size_t num)
  {
    Parameters parCov=getJJwbeta(start,beta);
    double postLogLik0=posteriorLogLikelihood(parCov,beta);
    // lets find the right factor
    std::size_t numsteps=100;
    MCMCrun mcmc;

    mcmc.factor=1.0/parCov.size();
    mcmc.beta=beta;
    Parameters par=parCov.randomSample(mcmc.factor);
    double postLogLik=posteriorLogLikelihood(par,beta)-postLogLik0;
    while(postLogLik!=postLogLik)
      {
        par=parCov.randomSample(mcmc.factor);
        postLogLik=posteriorLogLikelihood(par,beta)-postLogLik0;
      }

    std::ofstream f;

    f.open(filename_.c_str(),std::ios_base::app);


    while ((std::abs(postLogLik)>2)||(postLogLik!=postLogLik))
      {
        mcmc.factor/=sqrt(2);
        par=parCov.randomSample(mcmc.factor);
        postLogLik=posteriorLogLikelihood(par,beta)-postLogLik0;

      }
    while (std::abs(postLogLik)<0.25)
      {
        mcmc.factor*=sqrt(2);
        par=parCov.randomSample(mcmc.factor);
        postLogLik=posteriorLogLikelihood(par,beta)-postLogLik0;
        while(postLogLik!=postLogLik)
          {
            par=parCov.randomSample(mcmc.factor);
            postLogLik=posteriorLogLikelihood(par,beta)-postLogLik0;
          }
      }
    double sump=0;
    for (std::size_t i=0; i<numsteps; i++)
      {
        par=parCov.randomSample(mcmc.factor);
        sump+=exp((posteriorLogLikelihood(par,beta)-postLogLik0));
      }
    double logmeanp=-2*log(sump/numsteps);
    const double dssN=-2*log(0.23);
    if (!(logmeanp!=logmeanp))
      mcmc.factor=std::sqrt(dssN/logmeanp)*mcmc.factor;
    else
      mcmc.factor/=2;
    sump=0;

    for (std::size_t i=0; i<numsteps; i++)
      {
        par=parCov.randomSample(mcmc.factor);
        postLogLik=posteriorLogLikelihood(par,beta);
        sump+=exp(postLogLik-postLogLik0);
      }
    logmeanp=-2*log(sump/numsteps);
    if (!(logmeanp!=logmeanp))
      mcmc.factor=std::sqrt(dssN/logmeanp)*mcmc.factor;
    else
      mcmc.factor/=2;
    mcmc.sumLogL=0;
    mcmc.sumLogLh=0;
    mcmc.nh=0;
    mcmc.sumLogL2=0;
    mcmc.sumLogL2h=0;
    mcmc.sumP=0;

    double lik=exp(postLogLik-postLogLik0);


    for (std::size_t i=0; i<num*1.1; i++ )
      {

        std::size_t j=0;
        std::size_t nsubsteps=par.size();
        for (std::size_t ii=0; ii<par.size(); ii++)
          {
            Parameters nextPar;
            nextPar=par.randomSample(mcmc.factor);

            double  nextpostLogL=posteriorLogLikelihood(nextPar,beta);

            double nextLik=0;
            if (!nextpostLogL!=nextpostLogL)
              nextLik=exp(nextpostLogL-postLogLik0);

            double r=(1.0*rand())/RAND_MAX;

            if (nextLik/lik>r)
              {
                par=nextPar;
                postLogLik=nextpostLogL;

                lik=nextLik;
                j++;
              }
          }
        double ratio=(1.0*j)/nsubsteps;

        if ((ratio<0.10)&&(mcmc.factor>1.5e-2))
          mcmc.factor/=sqrt(2);
        else if (ratio>0.30)
          mcmc.factor*=sqrt(2);

        if (i>num*0.01)
          {
            double logL=logLikelihood(par);
            double ss=SumWeighedSquare(par);
            double sspar=SumWeighedSquareParameters(par);

            mcmc.logL.push_back(logL);
            mcmc.p.push_back(par);
            f<<i<<"\t"<<beta<<"\t"<<mcmc.factor<<"\t"<<ratio<<"\t"<<logL<<"\t"<<ss<<"\t"<<ss-sspar<<"\n";
            std::cerr<<i<<"\t"<<beta<<"\t"<<mcmc.factor<<"\t"<<ratio<<"\t"<<logL<<"\t"<<ss<<"\t"<<ss-sspar<<"\n";
            mcmc.sumLogL+=logL;
            mcmc.sumLogL2=logL*logL;
            mcmc.sumP+=exp(logL);

            if (i>num*0.6)
              {
                mcmc.sumLogLh+=logL;
                mcmc.nh++;

              }
          }


      }

    mcmc.meanLogL=mcmc.sumLogL/mcmc.logL.size();
    mcmc.meanLogLh=mcmc.sumLogLh/mcmc.nh;
    mcmc.stdLogL=sqrt(mcmc.sumLogL2/mcmc.logL.size()-mcmc.meanLogL*mcmc.meanLogL);
    mcmc.stdLogLh=sqrt(mcmc.sumLogL2h/mcmc.nh-mcmc.meanLogLh*mcmc.meanLogLh);

    mcmc.HMA=log(mcmc.logL.size()/mcmc.sumP);


    return mcmc;

  }



  MCMCrun
      BayesIteration::runMCMC(const Parameters& start,
                              double beta,
                              double factor,
                              std::size_t num)
  {



    // lets find the right factor

    MCMCrun mcmc;
    mcmc.beta=beta;
    mcmc.factor=factor;
    mcmc.sumLogL=0;
    mcmc.sumLogLh=0;
    mcmc.nh=0;
    mcmc.sumLogL2=0;
    mcmc.sumLogL2h=0;
    mcmc.sumP=0;

    Parameters par=getJJwbeta(start,beta);

    double postLogLik0=posteriorLogLikelihood(par,beta);
    double postLogLik=postLogLik0;




    std::ofstream f;

    f.open(filename_.c_str(),std::ios_base::app);





    double lik=exp(postLogLik-postLogLik0);

    for (std::size_t i=0; i<num*1.1; i++ )
      {

        std::size_t j=0;
        std::size_t nsubsteps=par.size()*3;
        for (std::size_t ii=0; ii<par.size(); ii++)
          {
            Parameters nextPar;
            nextPar=par.randomSample(mcmc.factor);

            double  nextpostLogL=posteriorLogLikelihood(nextPar,beta);

            double nextLik=0;
            if (!nextpostLogL!=nextpostLogL)
              nextLik=exp(nextpostLogL-postLogLik0);

            double r=(1.0*rand())/RAND_MAX;

            if (nextLik/lik>r)
              {
                par=nextPar;
                postLogLik=nextpostLogL;

                lik=nextLik;
                j++;
              }
          }
        double ratio=(1.0*j)/nsubsteps;

        if (false)
          {
            if (ratio<0.10)
              mcmc.factor/=sqrt(2);
            else if (ratio>0.30)
              mcmc.factor*=sqrt(2);
          }

        //         if (i>num*0.1)
        if (i>0)
          {

            double logL=logLikelihood(par);

            double ss=SumWeighedSquare(par);
            double sspar=SumWeighedSquareParameters(par);

            mcmc.logL.push_back(logL);
            mcmc.p.push_back(par);
            f<<i<<"\t"<<beta<<"\t"<<mcmc.factor<<"\t"<<ratio<<"\t"<<logL<<"\t"<<ss<<"\t"<<ss-sspar<<"\n";
            std::cerr<<i<<"\t"<<beta<<"\t"<<mcmc.factor<<"\t"<<ratio<<"\t"<<logL<<"\t"<<ss<<"\t"<<ss-sspar<<"\n";
            mcmc.sumLogL+=logL;
            mcmc.sumLogL2=logL*logL;
            mcmc.sumP+=exp(logL);
            if(i>num*0.6)
              {
                mcmc.sumLogLh+=logL;
                mcmc.sumLogL2h=logL*logL;
                mcmc.nh++;
              }

          }


      }

    mcmc.meanLogL=mcmc.sumLogL/mcmc.logL.size();
    mcmc.meanLogLh=mcmc.sumLogLh/mcmc.nh;
    mcmc.stdLogL=sqrt(mcmc.sumLogL2/mcmc.logL.size()-mcmc.meanLogL*mcmc.meanLogL);
    mcmc.stdLogLh=sqrt(mcmc.sumLogL2h/mcmc.nh-mcmc.meanLogLh*mcmc.meanLogLh);

    mcmc.HMA=log(mcmc.logL.size()/mcmc.sumP);

    return mcmc;

  }





  MCMCrun
      BayesIteration::runMCMC(BayesIteration finalModel,
                              const Parameters& start,
                              double beta,
                              double factor,
                              std::size_t num)
  {



    // lets find the right factor

    MCMCrun mcmc;
    mcmc.beta=beta;
    mcmc.factor=factor;
    mcmc.sumLogL=0;
    mcmc.sumLogLh=0;
    mcmc.nh=0;
    mcmc.sumLogL2=0;
    mcmc.sumLogL2h=0;
    mcmc.sumP=0;

    Parameters parCov0=getJJwbeta(start,1.0);
    Parameters parCov1=finalModel.getJJwbeta(start,1.0);

    Parameters par=start;


    par.setCovariance(inv(sum(mult(inv(parCov0.getCovariance()),1.0-beta),mult(inv(parCov1.getCovariance()),beta))));

    double postLogLik0=posteriorLogLikelihood(par,1.0)*(1.0-beta)+
        finalModel.posteriorLogLikelihood(par,1.0)*(beta);

    double postLogLik=postLogLik0;


    std::ofstream f;

    f.open(filename_.c_str(),std::ios_base::app);

    double lik=exp(postLogLik-postLogLik0);

    f<<"i"<<"\t"<<"beta"<<"\t"<<"factor"<<"\t"<<"ratio"<<"\t"<<"logL1-logL0"<<"\t"<<"logL1"<<"\t";
    f<<"logL0"<<"\t"<<"ss1"<<"\t"<<"ss0"<<"\n";
    std::cerr<<"i"<<"\t"<<"beta"<<"\t"<<"factor"<<"\t"<<"ratio"<<"\t"<<"logL1-logL0"<<"\t"<<"logL1"<<"\t";
    std::cerr<<"logL0"<<"\t"<<"ss1"<<"\t"<<"ss0"<<"\n";

    for (std::size_t i=0; i<num*1.1; i++ )
      {

        std::size_t j=0;
        std::size_t nsubsteps=par.size()*3;
        for (std::size_t ii=0; ii<par.size(); ii++)
          {
            Parameters nextPar;
            nextPar=par.randomSample(mcmc.factor);

            double  nextpostLogL=posteriorLogLikelihood(nextPar,1.0)*(1.0-beta)+
                finalModel.posteriorLogLikelihood(nextPar,1.0)*beta;

            double nextLik=0;
            if (!nextpostLogL!=nextpostLogL)
              nextLik=exp(nextpostLogL-postLogLik0);

            double r=(1.0*rand())/RAND_MAX;

            if (nextLik/lik>r)
              {
                par=nextPar;
                postLogLik=nextpostLogL;

                lik=nextLik;
                j++;
              }
          }
        double ratio=(1.0*j)/nsubsteps;

        if (false)
          {
            if (ratio<0.10)
              mcmc.factor/=sqrt(2);
            else if (ratio>0.30)
              mcmc.factor*=sqrt(2);
          }

        if (i>num*0.1)
          {

            double logL0=posteriorLogLikelihood(par,1.0);

            double logL1=finalModel.posteriorLogLikelihood(par,1.0);

            double logL=logL1-logL0;
            double ss0=SumWeighedSquare(par);
            double ss1=finalModel.SumWeighedSquare(par);

            mcmc.logL.push_back(logL);
            mcmc.p.push_back(par);
            f<<i<<"\t"<<beta<<"\t"<<mcmc.factor<<"\t"<<ratio<<"\t"<<logL1-logL0<<"\t"<<logL1<<"\t";
            f<<logL0<<"\t"<<ss1<<"\t"<<ss0<<"\n";
            std::cerr<<i<<"\t"<<beta<<"\t"<<mcmc.factor<<"\t"<<ratio<<"\t"<<logL1-logL0<<"\t"<<logL1<<"\t";
            std::cerr<<logL0<<"\t"<<ss1<<"\t"<<ss0<<"\n";
            mcmc.sumLogL+=logL;
            mcmc.sumLogL2=logL*logL;
            mcmc.sumP+=exp(logL);

            if (i>num*0.6)
              {
                mcmc.sumLogLh+=logL;
                mcmc.sumLogL2h=logL*logL;
                mcmc.nh++;
              }




          }


      }

    mcmc.meanLogL=mcmc.sumLogL/mcmc.logL.size();
    mcmc.meanLogLh=mcmc.sumLogLh/mcmc.nh;
    mcmc.stdLogL=sqrt(mcmc.sumLogL2/mcmc.logL.size()-mcmc.meanLogL*mcmc.meanLogL);
    mcmc.stdLogLh=sqrt(mcmc.sumLogL2h/mcmc.nh-mcmc.meanLogLh*mcmc.meanLogLh);


    return mcmc;

  }



  void BayesIteration::runMCMC(const Parameters &start, double factor, std::size_t num)
  {
    std::ofstream f;

    filename_="rMCMC"+filename_;

    f.open(filename_.c_str(),std::ios_base::app);

    f<<"------------------------------------------------------\n";
    f<<"-----------Monte Carlo Markov Chain-------------------\n";
    f<<"------------------------------------------------------\n";
    f.close();

    MCMCrun mcmc=runMCMC(start,1,factor,num);
    f.open(filename_.c_str(),std::ios_base::app);
    f<<mcmc.p;
    f.close();
    quantilesMCMC(mcmc);
    runQuantilesMCMC(mcmc);
  }


  void
      BayesIteration::quantilesMCMC(MCMCrun mcmcrun)
  {
    std::ofstream f;

    f.open(filename_.c_str(),std::ios_base::app);
    f<<"\n MCMC Analisis\n";
    f<<"Parameters  quantiles\n";
    f<<"Parameter \tq0.025\tq0.25\tq0.5\tq0.75\tq0.975 \n";
    Parameters p=mcmcrun.p.back();

    for (std::size_t i=0; i<p.size(); i++)
      {
        std::vector<double> pvalues;
        for (std::size_t j=0; j<mcmcrun.p.size(); j++)
          {
            pvalues.push_back(pow(10,mcmcrun.p[j][i]));
          }
        sort(pvalues);

        f<<p.indexToName(i)<<"\t"<<getQuantil(pvalues,0.025);
        f<<"\t"<<getQuantil(pvalues,0.25);
        f<<"\t"<<getQuantil(pvalues,0.5);
        f<<"\t"<<getQuantil(pvalues,0.75);
        f<<"\t"<<getQuantil(pvalues,0.975);
        f<<"\n";

      }
    f<<"End quantiles\n";


  }



  void
      BayesIteration::runQuantilesMCMC(MCMCrun mcmcrun)
  {
    std::ofstream f;

    f.open(filename_.c_str(),std::ios_base::app);
    f<<"\n MCMC Analisis\n";
    f<<"Run  quantiles\n";
    f<<"Run \tq0.025\tq0.975 \n";

    m_->run(f,mcmcrun.p,0.025,0.975);
    f<<"Run end\n";

    f<<"Run \tq0.25\tq0.75 \n";

    m_->run(f,mcmcrun.p,0.25,0.75);
    f<<"Run end\n";


    f<<"Run \tq0.5\tq0.5 \n";

    m_->run(f,mcmcrun.p,0.5,0.5);

    f<<"Run end\n";

  }







  Evidence BayesIteration::runThermodynamicIntegration(const Parameters& start, double factor,
                                                       std::size_t num,
                                                       std::size_t numBeta,
                                                       double alfa)
  {

    // melting

    // uso notacion de Syst. Biol. 55(2):195-207, 2006: Computing Bayes factors using thermodynamic integration

    // metodo de cambio de temperatura de
    // Estimating Bayes factors via thermodynamic integration and population MCMC,
    // "Ben Calderhead a,∗ , Mark Girolami"
    // Computational Statistics and Data Analysis 53 (2009) 4028–4045


    //  el tiempo que pasamos en cada temperatura es proporcional a
    // la raiz cuadrada de la esperanza del cuadrado de  la logLikelihood.

    // para la jumping distribution usamos la estimacion local del hessiano a atraves de la varianza del score.

    Evidence result;
    std::ofstream f;

    filename_="TI"+filename_;

    f.open(filename_.c_str(),std::ios_base::app);

    f<<"------------------------------------------------------\n";
    f<<"-----------Thermodynamic Integration-------------------\n";
    f<<"------------------------------------------------------\n";
    f.close();


    // first run a standard MCMC from the maximum


    double beta=1.0;

    MCMCrun mcmc0;
    MCMCrun mcmc;
    result.logZmelt=0;
    result.logZanneal=0;
    result.logZmelth=0;
    result.logZannealh=0;

    mcmc0=runMCMC(start,beta,factor,num*(1+numBeta/10));
    result.HMA=mcmc0.HMA;
    for (std::size_t i=1; i<numBeta; i++)
      {
        beta=1.0*pow(1.0*(numBeta-i)/numBeta,alfa);
        mcmc=runMCMC(mcmc0.p.back(),beta,factor,num);
        f.open(filename_.c_str(),std::ios_base::app);
        f<<"beta\t"<<beta<<"\tmeanLogL\t"<<mcmc.meanLogL<<"\tmeanLogL half\t"<<mcmc.meanLogLh<<"\n";
        result.logZmelt+=(mcmc0.meanLogL+mcmc.meanLogL)*(mcmc0.beta-mcmc.beta)/2.0;
        result.logZmelth+=(mcmc0.meanLogLh+mcmc.meanLogLh)*(mcmc0.beta-mcmc.beta)/2.0;

        mcmc0=mcmc;
      }

    result.logZmelt+=(mcmc.meanLogL)*(mcmc.beta)/2.0;
    result.logZmelth+=(mcmc.meanLogLh)*(mcmc.beta)/2.0;
    f.open(filename_.c_str(),std::ios_base::app);

    f<<"Thermodynamic integration melt logZ \t"<<result.logZmelt<<"\n";
    f<<"Thermodynamic integration melt last half logZ \t"<<result.logZmelth<<"\n";
    f<<"Harmonic mean  \t"<<result.HMA<<"\n";
    std::cerr<<"Thermodynamic integration melt logZ \t"<<result.logZmelt<<"\n";
    std::cerr<<"Thermodynamic integration last half melt logZ \t"<<result.logZmelth<<"\n";
    std::cerr<<"Harmonic mean  \t"<<result.HMA<<"\n";

    f.close();
    result.logZanneal+=(mcmc0.meanLogL)*(mcmc0.beta)/2.0;
    result.logZannealh+=(mcmc0.meanLogLh)*(mcmc0.beta)/2.0;

    for (std::size_t i=1; i<=numBeta; i++)
      {
        beta=1.0*pow((1.0*i)/numBeta,alfa);
        mcmc=runMCMC(mcmc0.p.back(),beta,factor,num);
        result.logZanneal+=(mcmc0.meanLogL+mcmc.meanLogL)*(mcmc.beta-mcmc0.beta)/2.0;
        result.logZannealh+=(mcmc0.meanLogLh+mcmc.meanLogLh)*(mcmc.beta-mcmc0.beta)/2.0;
        mcmc0=mcmc;
      }



    f.open(filename_.c_str(),std::ios_base::app);

    f<<"Thermodynamic integration anneal logZ \t"<<result.logZanneal<<"\n";
    std::cerr<<"Thermodynamic integration anneal logZ \t"<<result.logZanneal<<"\n";
    f<<"Thermodynamic integration last half anneal logZ \t"<<result.logZannealh<<"\n";
    std::cerr<<"Thermodynamic integration last half anneal logZ \t"<<result.logZannealh<<"\n";
    f.close();
    return result;
  }





  Evidence BayesIteration::runThermodynamicIntegration(BayesIteration& finalModel,
                                                       const Parameters& start,
                                                       double factor,
                                                       std::size_t num,
                                                       std::size_t numBeta,
                                                       double alfa)
  {

    // melting

    // uso notacion de Syst. Biol. 55(2):195-207, 2006: Computing Bayes factors using thermodynamic integration

    // metodo de cambio de temperatura de
    // Estimating Bayes factors via thermodynamic integration and population MCMC,
    // "Ben Calderhead a,∗ , Mark Girolami"
    // Computational Statistics and Data Analysis 53 (2009) 4028–4045


    //  el tiempo que pasamos en cada temperatura es proporcional a
    // la raiz cuadrada de la esperanza del cuadrado de  la logLikelihood.

    // para la jumping distribution usamos la estimacion local del hessiano a atraves de la varianza del score.

    Evidence result;
    std::ofstream f;
    filename_="SM"+filename_;

    f.open(filename_.c_str(),std::ios_base::app);

    f<<"------------------------------------------------------\n";
    f<<"-----------Thermodynamic Integration-------------------\n";
    f<<"-----------Model Switch--------------------------------\n";
    f<<"------------------------------------------------------\n";
    f.close();


    // first run a standard MCMC from the initial model


    double beta=0;

    MCMCrun mcmc0;
    MCMCrun mcmc;
    result.logZmelt=0;
    result.logZanneal=0;
    result.logZmelth=0;
    result.logZannealh=0;
    beta=1.0*pow(1.0/numBeta,alfa);

    mcmc0=runMCMC(finalModel,start,beta,factor,num*(1+numBeta/10));
    result.logZanneal+=mcmc0.meanLogL*mcmc0.beta;
    result.logZannealh+=mcmc0.meanLogLh*mcmc0.beta;
    result.HMA=mcmc0.HMA;
    for (std::size_t i=2; i<numBeta; i++)
      {
        beta=1.0*pow((1.0*i)/numBeta,alfa);
        mcmc=runMCMC(finalModel,mcmc0.p.back(),beta,factor,num);
        result.logZanneal+=(mcmc0.meanLogL+mcmc.meanLogL)*(mcmc.beta-mcmc0.beta)/2.0;
        result.logZannealh+=(mcmc0.meanLogLh+mcmc.meanLogLh)*(mcmc.beta-mcmc0.beta)/2.0;
        mcmc0=mcmc;
      }
    result.logZanneal+=mcmc0.meanLogL*(1.0-mcmc0.beta);
    result.logZannealh+=mcmc0.meanLogLh*(1.0-mcmc0.beta);

    f.open(filename_.c_str(),std::ios_base::app);

    f<<"Thermodynamic integration from M0 to M1 logZ \t"<<result.logZanneal<<"\n";
    std::cerr<<"Thermodynamic integration from M0 to M1 logZ \t"<<result.logZanneal<<"\n";
    f<<"Thermodynamic integration last half from M0 to M1 logZ \t"<<result.logZannealh<<"\n";
    std::cerr<<"Thermodynamic integration last half from M0 to M1 logZ \t"<<result.logZannealh<<"\n";

    f.close();
    result.logZmelt+=mcmc0.meanLogL*(1.0-mcmc0.beta);
    result.logZmelth+=mcmc0.meanLogLh*(1.0-mcmc0.beta);

    for (std::size_t i=2; i<numBeta; i++)
      {
        beta=1.0*pow(1.0*(numBeta-i)/numBeta,alfa);
        mcmc=runMCMC(finalModel,mcmc0.p.back(),beta,factor,num);
        result.logZmelt+=(mcmc0.meanLogL+mcmc.meanLogL)*(mcmc.beta-mcmc0.beta)/2.0;
        result.logZmelth+=(mcmc0.meanLogLh+mcmc.meanLogLh)*(mcmc.beta-mcmc0.beta)/2.0;
        mcmc0=mcmc;
      }
    result.logZmelt+=mcmc0.meanLogL*mcmc0.beta;
    result.logZmelth+=mcmc0.meanLogLh*mcmc0.beta;
    f.open(filename_.c_str(),std::ios_base::app);
    f<<"Thermodynamic integration from M1 to M0 logZ \t"<<result.logZmelt<<"\n";
    f<<"Thermodynamic integration last half from M1 to M0 logZ \t"<<result.logZmelth<<"\n";
    f<<"Harmonic mean  \t"<<result.HMA<<"\n";
    std::cerr<<"Thermodynamic integration from M1 to M0 logZ \t"<<result.logZmelt<<"\n";
    std::cerr<<"Thermodynamic integration last half from M1 to M0 logZ \t"<<result.logZmelth<<"\n";
    std::cerr<<"Harmonic mean  \t"<<result.HMA<<"\n";



    f.close();

    return result;
  }






  Parameters BayesIteration::runImportanceSampling(Parameters& gDist,
                                                   std::size_t num,
                                                   double& factor)
  {
    double ss0=SumWeighedSquare(gDist);

    std::vector<Parameters> pars (num);
    std::vector<double> dSSs (num);
    std::vector<double> logImportances(num);
    std::vector<double> chi2s(num,0);
    for (std::size_t i=0; i<num; i++ )
      {
        std::pair<Parameters,double> nextPar=gDist.randomSampleChi2(factor);

        pars[i]=nextPar.first;
        dSSs[i]=ss0-SumWeighedSquare(nextPar.first);
        logImportances[i]=dSSs[i]+nextPar.second;
        chi2s[i]=gDist.chi2Distance(nextPar.first);
      }


    double Sumw=0;
    double sumSS=0;
    double wcv=0;
    std::size_t npar=gDist.size();
    std::vector<double> m(npar,0);
    std::vector<std::vector<double> > C(npar,std::vector<double>(npar,0));

    for(std::size_t im=0; im<num; ++im)
      {
        if (!(logImportances[im]!=logImportances[im]))
          {
            double w=exp(0.5*logImportances[im]);
            Sumw+=w;
            sumSS+=dSSs[im]*w;
            wcv+=w*w;
            for (std::size_t i=0; i<npar; i++)
              m[i]+=pars[im][i]*w;

          }
      }
    double SSm=sumSS/Sumw;
    for (std::size_t i=0; i<npar; i++)
      m[i]/=Sumw;
    double wm=Sumw/num;
    wcv=std::sqrt(wcv/num-wm*wm);
    wcv=wcv/wm;

    for(std::size_t im=0; im<num; ++im)
      {
        double w=exp(logImportances[im]);
        for (std::size_t i=0; i<npar; i++)
          for (std::size_t j=0; j<npar; j++)
            C[i][j]+=(pars[im][i]-m[i])*
                (pars[im][j]-m[j])*w;
      }

    for (std::size_t i=0; i<npar; i++)
      for (std::size_t j=0; j<npar; j++)
        C[i][j]/=Sumw;


    Parameters result(gDist);
    result.setpMeans(m);
    result.setCovariance(C);

    return result;
  }

  Evidence
      BayesIteration::anaImportanceSampling(LevenbergMarquardt& LM,
                                            std::size_t num)
  {
    const double PI=3.1415926;
    std::size_t ndata=getData().size();
    Parameters MAP=LM.OptimParameters();
    std::size_t npar=MAP.size();
    std::vector<double> m(npar,0);
    std::vector<std::vector<double> > C(npar,std::vector<double>(npar,0));

    double ss0=SumWeighedSquare(MAP);

    std::vector<Parameters> pars (num);
    std::vector<double> dSSs (num);
    std::vector<double> logImportances(num);
    double factor=1.0;

    std::ofstream f;

    f.open(filename_.c_str(),std::ios_base::app);

    f<<"------------------------------------------------------\n";
    f<<"-----------Importance Sampling based on MAP-------------------\n";
    f<<"------------------------------------------------------\n";


    factor=0.2;
    Parameters prun=runImportanceSampling(MAP,num*npar*0.1,factor);

    factor=1;
    Parameters prun2=runImportanceSampling(MAP,num*npar*0.1,factor);

    std::vector<double> chi2s(num,0);
    for (std::size_t i=0; i<num; i++ )
      {
        std::pair<Parameters,double> nextPar=MAP.randomSampleChi2(factor);

        pars[i]=nextPar.first;
        dSSs[i]=ss0-SumWeighedSquare(nextPar.first);
        logImportances[i]=dSSs[i]+nextPar.second;
        chi2s[i]=MAP.chi2Distance(nextPar.first);
      }


    double Sumw=0;
    double sumSS=0;
    double wcv=0;

    for(std::size_t im=0; im<num; ++im)
      {
        if (!(logImportances[im]!=logImportances[im]))
          {
            double w=exp(0.5*logImportances[im]);
            Sumw+=w;
            sumSS+=dSSs[im]*w;
            wcv+=w*w;
            for (std::size_t i=0; i<npar; i++)
              m[i]+=pars[im][i]*w;

          }
      }
    double SSm=sumSS/Sumw;
    for (std::size_t i=0; i<npar; i++)
      m[i]/=Sumw;
    double wm=Sumw/num;
    wcv=std::sqrt(wcv/num-wm*wm);
    wcv=wcv/wm;

    for(std::size_t im=0; im<num; ++im)
      {
        double w=exp(logImportances[im]);
        for (std::size_t i=0; i<npar; i++)
          for (std::size_t j=0; j<npar; j++)
            C[i][j]+=(pars[im][i]-m[i])*
                (pars[im][j]-m[j])*w;
      }

    for (std::size_t i=0; i<npar; i++)
      for (std::size_t j=0; j<npar; j++)
        C[i][j]/=Sumw;

    Evidence result;
    result.p=MAP;
    result.p.setpMeans(m);
    result.p.setCovariance(C);
    double logdetCovMAP=log(det(MAP.getCovariance()));
    double logdetCovImp=logdetCovMAP+log(wm);
    result.ImportanceEvidence=-0.5*ss0+
        0.5*logdetCovImp-0.5*logDetCovData()-
        0.5*LM.logDetPriorCov()-0.5*(ndata-npar)*log(2*PI)+log(wm);
    result.ImportanceNumber=wcv;

    result.logdetCovPost=log(det(C));

    result.logEvidence=-0.5*ss0+0.5*result.logdetCovPost-0.5*logDetCovData()-
        0.5*LM.logDetPriorCov()-0.5*(ndata-npar)*log(2*PI);

    result.logPostLik=-0.5*ss0-0.5*LM.logDetPriorCov()
        -0.5*logDetCovData()-0.5*ndata*log(2*PI);


    return result;
  }











  Parameters BayesIteration::getEvidence(const Parameters& maximumPostLik, std::size_t num)
  {
    double ss0=SumWeighedSquare(maximumPostLik);
    std::vector<Parameters> parvec;
    std::vector<double> ssvec;
    std::vector<double> m=maximumPostLik.pMeans();
    std::vector<std::vector<double> > covinv=inv(maximumPostLik.getCovariance());
    std::size_t npar=maximumPostLik.size();
    double detCovinv=det(covinv);
    double detdiam=std::pow(detCovinv,1.0/npar);

    std::ofstream f;

    f.open(filename_.c_str(),std::ios_base::app);

    f<<"------------------------------------------------------\n";
    f<<"-----------Markov chain Monte Carlo-------------------\n";
    f<<"------------------------------------------------------\n";


    double sumw=0;
    double sumP=0;
    //   try a Metropolis


    // lets find the right factor
    std::size_t numsteps=100;

    double factor=1.0/maximumPostLik.size();
    Parameters par=maximumPostLik.randomSample(factor);
    double ss=SumWeighedSquare(par)-ss0;
    while ((ss)>4)
      {
        factor/=sqrt(2);
        par=maximumPostLik.randomSample(factor);
        ss=SumWeighedSquare(par)-ss0;
      }
    while ((ss)<0.5)
      {
        factor*=sqrt(2);
        par=maximumPostLik.randomSample(factor);
        ss=SumWeighedSquare(par)-ss0;
      }
    double sumdss=0;
    for (std::size_t i=0; i<numsteps; i++)
      {
        par=maximumPostLik.randomSample(factor);
        sumdss+=SumWeighedSquare(par)-ss0;
      }
    double meandss=sumdss/numsteps;
    factor=std::sqrt(-2*log(0.23)/meandss)*factor;
    sumdss=0;
    for (std::size_t i=0; i<numsteps; i++)
      {
        par=maximumPostLik.randomSample(factor);
        ss=SumWeighedSquare(par);
        sumdss+=ss-ss0;

      }
    meandss=sumdss/numsteps;

    f<<"valor medio de ss diference ="<<meandss<<"\tfactor= "<<factor<<"\n";

    double lik=exp(-0.5*(ss-ss0));

    for (std::size_t i=0; i<num*1.1; i++ )
      {

        std::size_t j=0;
        for (std::size_t ii=0; ii<par.size(); ii++)
          {
            Parameters nextPar=par.randomSample(factor);

            double  nextss=SumWeighedSquare(nextPar);
            double nextLik=exp(-0.5*(nextss-ss0));

            double r=(1.0*rand())/RAND_MAX;

            if (nextLik/lik>r)
              {
                par=nextPar;
                ss=nextss;
                lik=nextLik;
                j++;
              }
          }
        if (i>num*0.1)
          {
            parvec.push_back(par);
            ssvec.push_back(ss);
            f<<"muestra nro "<<i<<"\tsaltos efectivos\t"<<j<<"\tss\t"<<ss<<std::endl;

          }
        std::cerr<<i<<" "<<std::endl;

      }

    // now calculate expected mean and covariance
    Parameters result(maximumPostLik);
    std::vector<double> mout(m.size(),0);
    std::vector<std::vector<double> > covout(m.size(),std::vector<double>(m.size(),0));


    double sss=0;
    for(std::size_t im=0; im<parvec.size(); ++im)
      {
        sss+=ssvec[im]/parvec.size();
        for (std::size_t i=0; i<parvec[im].size(); i++)
          mout[i]+=parvec[im][i]/parvec.size();
      }
    for(std::size_t im=0; im<parvec.size(); ++im)
      {
        for (std::size_t i=0; i<parvec[im].size(); i++)
          for (std::size_t j=0; j<parvec[im].size(); j++)
            covout[i][j]+=(parvec[im][i]-mout[i])*(parvec[im][j]-mout[j])/
                parvec.size();
      }





    result.setpMeans(mout);
    result.setCovariance(covout);

    return result;




  }




  std::map<double,Parameters> BayesIteration::getRandomParameters(const Parameters& per,
                                                                  std::size_t num,
                                                                  double factor)
  {
    std::map<double,Parameters> myMap;
    for (std::size_t i=0; i<num; i++)
      {
        Parameters p=per.randomSample(factor);
        double ss=SumWeighedSquare(p);

        std::cout<<ss<<"\n";
        //    std::cout<<p;
        if (!(ss!=ss))
          myMap[ss]=p;
      }


    for (std::map<double,Parameters>::iterator it=myMap.begin();it!=myMap.end();++it)
      {
        double ss=(*it).first;
        std::cout<<(*it).first<<"\n";

      }


    return myMap;
  }

  std::map<double,Parameters> BayesIteration::getRandomParameters(std::size_t num,
                                                                  double factor)
  {
    std::map<double,Parameters> myMap;
    for (std::size_t i=0; i<num; i++)
      {
        Parameters p=priors_.back().randomSample(factor);
        double ss=SumWeighedSquare(p);

        std::cout<<ss<<"\n";
        //    std::cout<<p;
        if (!(ss!=ss))
          myMap[ss]=p;
      }


    for (std::map<double,Parameters>::iterator it=myMap.begin();it!=myMap.end();++it)
      {
        double ss=(*it).first;
        std::cout<<(*it).first<<"\n";

      }


    return myMap;
  }



