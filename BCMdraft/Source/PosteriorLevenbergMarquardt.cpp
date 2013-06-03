#include <cmath>
#include<limits>
#include <sstream>
#include <thread>
#include <fstream>
#include "Include/PosteriorLevenbergMarquardt.h"
#include "Include/MatrixInverse.h"


VariablesNormalDistribution PosteriorLevenbergMarquardt::OptimParameters()const
{
  return VariablesNormalDistribution(iter_.currParam_,this->JTWJinv_);
}


double PosteriorLevenbergMarquardt::elapsedTime()const
{
  std::chrono::nanoseconds d =iter_.t1_-t0_;
  return 1E-9*d.count();
}




std::size_t PosteriorLevenbergMarquardt::numEval()const
{
  return iter_.numEval();
}
std::size_t PosteriorLevenbergMarquardt::numIter()const
{
  return iter_.numIter();
}
double PosteriorLevenbergMarquardt::SS()const
{
  return iter_.SS();
}
std::vector<double> PosteriorLevenbergMarquardt::Gradient()const
{
  return iter_.Gradient();
}


PosteriorLevenbergMarquardt::
PosteriorLevenbergMarquardt(ABC_BCM* m,
                            const VariablesValue &initial,
                            const std::vector<const ExperimentDistribution *>& expset,
                            double dt):
  t0_(std::chrono::steady_clock::now()),
  initialParam_(initial),
  iter_(m,initial,expset,dt),
  convCriter_(convergence()),
  G_(std::vector<double>(initial.size())),
  JTWJ_(std::vector< std::vector<double> > (initial.size(),std::vector<double>(initial.size()))),
  JTWJinv_(std::vector< std::vector<double> >(initial.size(),std::vector<double>(initial.size())))
{}

PosteriorLevenbergMarquardt& PosteriorLevenbergMarquardt::optimize(std::size_t numIterations,
                                                                   double maxTime)
{
  std::cout<<iter_.header();
  std::stringstream ss1;
  ss1<<"optimize";
  ss1<<ABC_BCM::seed;
  ss1<<".txt";
  std::string filename=ss1.str();
  std::ofstream ff;
  ff.open(filename.c_str(),std::fstream::out|std::fstream::app);
  ff<<iter_.header()<<std::endl;
  convCriter_.addIter(numIterations);
  convCriter_.addTime(maxTime);
  while (!convCriter_.meetConvergenceCriteria(iter_))
    iter_.iterate();



  for (std::size_t i=0; i<iter_.nPar_; ++i)
    for (std::size_t j=0; j<iter_.nPar_; ++j)
      {
        if (i==j)
          JTWJ_[i][j]=iter_.wPar_[i];
        else
          JTWJ_[i][j]=0;
        for (std::size_t n=0; n<iter_.nData_; ++n)
          {
            JTWJ_[i][j]+=iter_.J_[n][i]*iter_.J_[n][j]*iter_.w_[n];
          }
      }

  JTWJinv_=inv(JTWJ_);
   /*    evidence_=0;
    logDetDataCov_=0;
    logDetPriorCov_=0;
    LogPostLik_=0;
    logDetPostStd_=0;
    for (std::size_t i=0; i<nData_-nPar_;i++)
    {
        logDetDataCov_-=log(w_[i]);
    }
    for (std::size_t i=nData_-nPar_; i<nData_;i++)
    {
        logDetPriorCov_-=log(w_[i]);
    }
    for (std::size_t i=0;i<currParam_.size(); i++)
    {
        logDetPostStd_-=log(JTWJinv_[i][i]);
    }

    const double PI=3.1415926;
    LogPostLik_+= -0.5*currSS_-0.5*logDetDataCov_-0.5*logDetPriorCov_-0.5*nData_*log(2*PI);
    logDetPostCov_=log(det(JTWJinv_));
    evidence_=LogPostLik_+0.5*logDetPostCov_+0.5*nPar_*log(2*PI);
    */

  std::cout<<report();
  return *this;
}

double PosteriorLevenbergMarquardt::testedValue::elapsedTime()const
{

  std::chrono::nanoseconds d =t1_-t0_;
  return 1E-9*d.count();
}

PosteriorLevenbergMarquardt::testedValue::
testedValue(ABC_BCM* m,
            VariablesValue initial,
            const std::vector<const ExperimentDistribution *>& expset,
            double dt):
  bcm_(m),
  exp(expset),
  t0_(std::chrono::steady_clock::now()),
  t1_(),
  dt_(dt),
  data_(m->Ty(expset)),
  w_(m->w(expset)),
  par_(m->getPrior()->center().tvalues()),
  wPar_(m->getPrior()->w()),
  nPar_(initial.size()),
  nData_(data_.size()),
  nIter_(0),
  nFeval_(0),
  currParam_(initial),
  newParam_(initial),
  newParam0_(initial),
  currYfit_(m->Tyfit(initial,expset,dt_)),
  currSS_(m->SSdata(currYfit_,data_,w_)+m->SSpar(initial)),
  J_(std::vector<std::vector<double>>(nData_,std::vector<double>(nPar_))),
  G_(std::vector<double>(nPar_)),
  JTWJ_(std::vector< std::vector<double> > (nPar_,std::vector<double>(nPar_))),
  JTWJinv_(std::vector< std::vector<double> >(nPar_,std::vector<double>(nPar_))),
  d_(std::vector<double>(nPar_)),
  ParamChange_(std::numeric_limits<double>::infinity()),
  SSChange_(std::numeric_limits<double>::infinity()),
  NormGrad_(std::numeric_limits<double>::infinity())
{
  nFeval_++;
}




void PosteriorLevenbergMarquardt::testedValue::iterate()
{
  computeJacobian();
  computeSearchDirection();
  updateLanda();
  t1_=std::chrono::steady_clock::now();
  std::stringstream ss1;
  ss1<<"optimize";
  ss1<<ABC_BCM::seed;
  ss1<<".txt";
  std::string filename=ss1.str();
  std::ofstream ff;
  ff.open(filename.c_str(),std::fstream::out|std::fstream::app);




  std::cerr<<nIter_<<"\t"<<currSS_<<"\t"<<elapsedTime()<<"\t\t"<<landa_<<"\t";
  std::cerr<<ParamChange_<<"\t"<<SSChange_<<"\t"<<NormGrad_<<"\n";

  ff<<nIter_<<"\t"<<currSS_<<"\t"<<elapsedTime()<<"\t\t"<<landa_<<"\t";
  ff<<ParamChange_<<"\t"<<SSChange_<<"\t"<<NormGrad_<<"\n";


  nIter_++;
}


void PosteriorLevenbergMarquardt::testedValue::computeJacobian()
{


  std::vector<double> curryf=bcm_->Tyfit(currParam_,exp,dt_);

#pragma omp parallel for
  for (std::size_t i=0; i<nPar_; i++)
    {
      VariablesValue x(currParam_);
      double tx=x.getTvalue(i);
       double ddx=dx_;
      x.setTvalue(i,tx+ddx);
      std::vector<double> yf=this->bcm_->Tyfit(x,exp,dt_);
      while (isNaN(yf))
        {
          ddx/=2;
          x.setTvalue(i,tx+ddx);
          yf=this->bcm_->Tyfit(x,exp,dt_);
        }


      for (std::size_t n=0;n<nData_;++n)
        {
          double jj=(yf[n]-curryf[n])/ddx;
          J_[n][i]=jj;
        }
    }


}



std::string PosteriorLevenbergMarquardt::testedValue::header()const
{
  std::stringstream ss;
  ss<<"nIter"<<"\t"<<"currSS"<<"\t\t"<<"time"<<"\t"<<"landa"<<"\t";
  ss<<"ParamChange"<<"\t"<<"SSChange"<<"\t"<<"NormGrad"<<"\n";
  return ss.str();
}

std::string PosteriorLevenbergMarquardt::report()const{

  std::stringstream output;

  output<<"Convergence critera: \t";
  if (convCriter_.surpassIter_)
    output<<iter_.nIter_<<" surpass number of iterations="<<convCriter_.maxIter_<<"\t";
  if (convCriter_.surpassFeval_)
    output<<iter_.nFeval_<<" surpass number of function evaluations="<<convCriter_.maxFeval_<<"\t";
  if (convCriter_.surpassLanda_)
    output<<iter_.landa_<<" surpass maximum landa value="<<convCriter_.maxLanda_<<"\t";
  if (convCriter_.smallParamChange_)
    output<<iter_.ParamChange_<<" surpass minium parameter change value="<<convCriter_.minParamChange_<<"\t";
  if (convCriter_.smallSSChange_)
    output<<iter_.SSChange_<<" surpass minium SS change value="<<convCriter_.minSSChange_<<"\t";
  if (convCriter_.smallGradient_)
    output<<iter_.NormGrad_<<" surpass minium gradient norm value="<<convCriter_.minGradient_<<"\t";
  if (iter_.currSS_!=iter_.currSS_)
    output<<iter_.currSS_<<" invalid value of the square sum"<<"\t";




  return output.str();


}





void PosteriorLevenbergMarquardt::testedValue::computeSearchDirection()
{
  for (std::size_t i=0; i<nPar_; ++i)
    for (std::size_t j=0; j<nPar_; ++j)
      {
        if (i==j)
          JTWJ_[i][j]=wPar_[i];
        else
          JTWJ_[i][j]=0;

        for (std::size_t n=0; n<nData_; ++n)
          {
            JTWJ_[i][j]+=J_[n][i]*J_[n][j]*w_[n];
          }
      }



  for (std::size_t i=0; i<nPar_; ++i)
    {
      JTWJ_[i][i]*=1+landa_;
    }

  JTWJinv_=inv(JTWJ_);



  for (std::size_t i=0; i<nPar_; ++i)
    {
      G_[i]=(par_[i]-currParam_.getTvalue(i))*wPar_[i];
      for (std::size_t n=0; n<nData_;++n)
        {
          G_[i]+=(data_[n]-currYfit_[n])*w_[n]*J_[n][i];
        }

    }
  for (std::size_t i=0; i<nPar_; ++i)
    {
      d_[i]=0;
      for (std::size_t j=0; j<nPar_;++j)
        {
          d_[i]+=JTWJinv_[i][j]*G_[j];
        }

    }
  newParam_=currParam_;
  newParam_.setTransfValues(currParam_.tvalues()+d_);

  newYfit_=bcm_->Tyfit(newParam_,exp,dt_);
  nFeval_++;
  if (newYfit_.size()==0)
    newSSW_=std::numeric_limits<double>::quiet_NaN();
  else
    {

      newSSW_=0;
      for (std::size_t n=0; n<nData_; n++)
        {
          newSSW_+=(newYfit_[n]-data_[n])*(newYfit_[n]-data_[n])*w_[n];
        }
      for (std::size_t n=0; n<nPar_; n++)
        {
          newSSW_+=std::pow((newParam_.getTvalue(n)-par_[n]),2)*wPar_[n];
        }
    }
}


void PosteriorLevenbergMarquardt::testedValue::updateLanda()
{
  std::size_t ifevalLoop=0;
  if ((newSSW_>=currSS_)||(newSSW_!=newSSW_))
    {
      while(((newSSW_>=currSS_)&&(ifevalLoop<maxFeval_))||(newSSW_!=newSSW_))
        {
          if (landa_*v_>=maxLanda_) break;
          landa0_=landa_;
          landa_=landa0_*v_;
          newSSW0_=newSSW_;
          newParam0_=newParam_;
          computeSearchDirection();
          ifevalLoop++;
          std::cerr<<"landa "<<landa_<<"\n";
        }

    }
  else
    {

      landa0_=landa_;
      landa_=landa_/v_;
      newSSW0_=newSSW_;
      newParam0_=newParam_;
      newYfit0_=newYfit_;
      computeSearchDirection();
      ifevalLoop++;
      while((newSSW_<newSSW0_)&&(!newSSW_!=newSSW_)&&(landa_>0.5))
        {
          landa0_=landa_;
          landa_=landa_/v_;
          newSSW0_=newSSW_;
          newParam0_=newParam_;
          newYfit0_=newYfit_;
          computeSearchDirection();
          ifevalLoop++;
        }
      if ((newSSW_>=newSSW0_)||(newSSW_!=newSSW_))
        {
          landa_=landa0_;
          newParam_=newParam0_;
          newSSW_=newSSW0_;
          newYfit_=newYfit0_;
        }
    }
  if (newSSW_>=currSS_)
    {
      ParamChange_=0;
      SSChange_=0;

    }
  else
    {
      ParamChange_=0;
      for (std::size_t i=0; i<nPar_; ++i)
        ParamChange_+=std::pow(currParam_.getTvalue(i)-newParam_.getTvalue(i),2);
      ParamChange_=sqrt(ParamChange_);
      SSChange_=currSS_-newSSW_;
      currParam_=newParam_;
      currSS_=newSSW_;
      currYfit_=newYfit_;
    }
  NormGrad_=0;
  for (std::size_t i=0; i<nPar_; ++i)
    NormGrad_+=G_[i]*G_[i];
  NormGrad_=sqrt(NormGrad_);
}


bool PosteriorLevenbergMarquardt::convergence::meetConvergenceCriteria(testedValue &val)
{
  surpassIter_=(bool(val.nIter_>=maxIter_)&&(val.currSS_>100))||(val.nIter_>=maxIter_*1);
  surpassFeval_=val.nFeval_>=maxFeval_;
  surpassLanda_=val.landa_>=maxLanda_;
  smallParamChange_=val.ParamChange_<minParamChange_;
  smallSSChange_=val.SSChange_<minSSChange_;
  smallGradient_=val.NormGrad_<minGradient_;
  invalidSS_=val.currSS_!=val.currSS_;


  return surpassIter_||
      surpassFeval_||
      smallParamChange_||
      smallSSChange_||
      smallGradient_||
      surpassLanda_||
      invalidSS_;
}



std::size_t PosteriorLevenbergMarquardt::testedValue::numEval()const
{
  return nFeval_;
}
std::size_t PosteriorLevenbergMarquardt::testedValue::numIter()const
{
  return nIter_;
}
double PosteriorLevenbergMarquardt::testedValue::SS()const
{
  return currSS_;
}
std::vector<double> PosteriorLevenbergMarquardt::testedValue::Gradient()const
{
  return G_;
}
std::string PosteriorLevenbergMarquardt::convergence::report(const testedValue& iter)const{
  std::stringstream output;

  output<<"Convergence critera: \t";
  if (surpassIter_)
    output<<iter.nIter_<<" surpass number of iterations="<<maxIter_<<"\t";
  if (surpassFeval_)
    output<<iter.nFeval_<<" surpass number of function evaluations="<<maxFeval_<<"\t";
  if (surpassLanda_)
    output<<iter.landa_<<" surpass maximum landa value="<<maxLanda_<<"\t";
  if (smallParamChange_)
    output<<iter.ParamChange_<<" surpass minium parameter change value="<<minParamChange_<<"\t";
  if (smallSSChange_)
    output<<iter.SSChange_<<" surpass minium SS change value="<<minSSChange_<<"\t";
  if (smallGradient_)
    output<<iter.NormGrad_<<" surpass minium gradient norm value="<<minGradient_<<"\t";
  if (iter.currSS_!=iter.currSS_)
    output<<iter.currSS_<<" invalid value of the square sum"<<"\t";




  return output.str();



}


std::ostream& operator<<(std::ostream& s, PosteriorLevenbergMarquardt& LM)
{
  /*
    s<<"LevenbergMarquardtParameters\n";
    s<<"Number of iterations \t"<<LM.numIter()<<"\t";
    s<<"Sum of squares"<<LM.currSS_<<"\n";
    s<<"Report \t"<<LM.report()<<"\n";
    s<<"Landa value \t"<<LM.landa_<<"\t";
    s<<"Parameter change \t"<<LM.minParamChange_<<"\t";
    s<<"SS change \t"<<LM.minSSChange_<<"\t";
    s<<"Gradient norm \t"<<LM.NormGrad_<<"\n";
    s<<"Evidence \t"<<LM.evidence_<<"\n";
    s<<"End\n";
*/
}



