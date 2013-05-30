#ifndef POSTERIORLEVENBERGMARQUARDT_H
#define POSTERIORLEVENBERGMARQUARDT_H
#include <vector>
#include <iostream>
#include "ABC_BCM.h"
#include <chrono>


class PosteriorLevenbergMarquardt
{
public:
  VariablesNormalDistribution OptimParameters()const;
  double elapsedTime()const;
  std::size_t numEval()const;
  std::size_t numIter()const;
  double SS()const;
  std::vector<double> Gradient()const;

   PosteriorLevenbergMarquardt(ABC_BCM* m,
                              const VariablesValue& initial,
                              const std::vector<const ExperimentDistribution*>& expset,
                               double dt);

    PosteriorLevenbergMarquardt& optimize(std::size_t numIterations, double maxTime=0);


/*
    double getEvidence()const;

    double getLogPostLik()const;

    double logDetPriorCov()const;
    double logDetPostCov()const;
    double logDetPostStd()const;
    double SSdata()const;

*/
    /*
    PosteriorLevenbergMarquardt(const PosteriorLevenbergMarquardt& other);

    friend void swap(PosteriorLevenbergMarquardt& one, PosteriorLevenbergMarquardt& other);

    PosteriorLevenbergMarquardt& operator=(const PosteriorLevenbergMarquardt& other);

    PosteriorLevenbergMarquardt();

    ~PosteriorLevenbergMarquardt(){}
    */
    std::string report()const;

   // void reset(const SimParameters& sp,const Treatment& tr);


    friend std::ostream& operator<<(std::ostream& s, PosteriorLevenbergMarquardt& LM);

    class convergence;
    class testedValue
    {
    public:
      double elapsedTime()const;
      std::size_t numEval()const;
      std::size_t numIter()const;
      double SS()const;
      std::vector<double> Gradient()const;


      testedValue(ABC_BCM* m,
                  VariablesValue initial,
                  const std::vector<const ExperimentDistribution*>& expset,
                  double dt);

      std::string header()const;
      void iterate();
      void computeJacobian();
      void computeSearchDirection();
      void updateLanda();
    friend class  PosteriorLevenbergMarquardt::convergence;
    friend class PosteriorLevenbergMarquardt;
    private:
     ABC_BCM* bcm_;
     const std::vector<const ExperimentDistribution*> exp;
     std::chrono::steady_clock::time_point t0_;
     std::chrono::steady_clock::time_point t1_;

     double dx_=1e-8;
     double dt_;

     std::vector<double> data_;
     std::vector<double> w_;

     std::vector<double> par_;
     std::vector<double> wPar_;


     std::size_t nPar_;
     std::size_t nData_;

      double landa_=1E2;
      double landa0_;
      double v_=3;
      double maxLanda_=1E10;

      std::size_t nIter_;
      std::size_t nFeval_;
      std::size_t maxFeval_=30;

      VariablesValue currParam_;
      VariablesValue newParam_;
      VariablesValue newParam0_;
      std::vector<double> dPar_;


      std::vector<double> currYfit_;
      std::vector<double> newYfit_;
      std::vector<double> newYfit0_;

      double currSS_;
      double newSSW_;
      double newSSW0_;



      std::vector< std::vector< double> > J_;
      std::vector<double> G_;
      std::vector< std::vector<double> > JTWJ_;
      std::vector< std::vector<double> > JTWJinv_;
      std::vector<double> d_;

      double evidence_;
      double LogPostLik_;
      double logDetDataCov_;
      double logDetPriorCov_;
      double logDetPostCov_;
      double logDetPostStd_;
      double SSdata_;

      double ParamChange_;
      double SSChange_;
      double NormGrad_;


    };

    class convergence
    {
    public:
       convergence(std::size_t maxIter,
                   double maxTime,
                 std::size_t maxFeval,
                 double minParamChange,
                 double minSSChange,
                 double minGradient,
                 double maxLanda):
         maxIter_(maxIter),
         maxTime_(maxTime),
         maxFeval_(maxFeval),
         minParamChange_(minParamChange),
         minSSChange_(minSSChange),
         minGradient_(minGradient),
         maxLanda_(maxLanda){}
       convergence(){}

      bool meetConvergenceCriteria(testedValue& val);
      std::string report(const testedValue& iter)const;
      void addIter(std::size_t maxNumberIterations){maxIter_+=maxNumberIterations;}
      void addTime(double maxElapsedTime){maxTime_+=maxElapsedTime;}

      friend class PosteriorLevenbergMarquardt;
    private:
      std::size_t maxIter_=1;
      std::size_t maxFeval_=10000;
      double maxTime_=1000;

      double minParamChange_=1e-6;
      double minSSChange_=1e-6;
      double minGradient_=1e-3;
      double maxLanda_=1e7;

      bool surpassIter_=false;
      bool surpassFeval_=false;
      bool surpassLanda_=false;
      bool smallParamChange_=false;
      bool smallSSChange_=false;
      bool smallGradient_=false;
      bool invalidSS_=false;

    };

private:
    std::chrono::steady_clock::time_point t0_;
    VariablesValue initialParam_;
    testedValue iter_;
    convergence convCriter_;
    std::vector<double> G_;
    std::vector< std::vector<double> > JTWJ_;
    std::vector< std::vector<double> > JTWJinv_;

};







#endif // POSTERIORLEVENBERGMARQUARDT_H
