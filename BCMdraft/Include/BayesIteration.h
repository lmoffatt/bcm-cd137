#ifndef BAYESITERATION_H
#define BAYESITERATION_H
#include <vector>
#include <string>
#include <map>
#include "LevenbergMarquardtParameters.h"


class ABC_data
{
public:
    virtual std::vector<double> getData()const=0;
    virtual std::vector<double> getDataStandardError()const=0;
    virtual std::vector<double> getDataWeigth();


    virtual ~ABC_data();
};


class ABC_model:public ABC_function
{
public:

    virtual ABC_model& setData(const ABC_data& experimentalData)=0;
    virtual ~ABC_model();
    virtual std::ostream&
    run(std::ostream& f,const std::vector<std::pair <double,Parameters> >& pars)const=0;
  virtual std::ostream&
  run(std::ostream& f,const std::vector<Parameters >& pars,double quantil0, double quantil1)const=0;

};



struct Evidence
{
    Parameters p;
    double ss;
    double ssStd;
    double HMA;
    double HMAcv;
    double logPostLik;
    double logdetCovPost;
    double logEvidence;
    double ImportanceEvidence;
    double ImportanceNumber;
    double logZmelt;
    double logZanneal;
    double logZmelth;
    double logZannealh;

};



struct MCMCrun
{
  std::vector<Parameters> p;
  std::vector<double> logL;
  double beta;
  double factor;
  double sumLogL;
  double sumLogLh;
  std::size_t nh;
  double sumLogL2;
  double sumLogL2h;
  double meanLogL;
  double meanLogLh;

  double stdLogL;
  double stdLogLh;
  double sumP;
  double HMA;
};



std::ostream& operator<<(std::ostream& s, const Evidence& e);




class BayesIteration:public ABC_model, public ABC_data
{
public:

    Parameters Posterior()const;

    Parameters Prior(std::size_t n=0)const;


    BayesIteration(const ABC_model* f,
                   Parameters prior,
                   const ABC_data* d,
                   const std::string &filename);


    std::map<double,Parameters> getRandomParameters(std::size_t num,double factor);

    std::map<double,Parameters> getRandomParameters(const Parameters& per,std::size_t num, double factor);


    double SumWeighedSquare(const Parameters& p);
    double SumWeighedSquareParameters(const Parameters& p);

    double posteriorLogLikelihood(const Parameters& p, double beta);
    double logLikelihood(const Parameters& p);

    std::pair<double,double> posteriorLogLikelihood(const Parameters& p);


    std::vector<double> partialPosteriorLogLikelihood(const Parameters& p, double beta);


    BayesIteration& addNewData(ABC_data* d);

    virtual ~BayesIteration();

    virtual void setFilename(const std::string filename);

    virtual std::vector<double> yfit (const std::vector<double>& param);
    virtual std::vector<double> yfit(const Parameters& parameters)const;

    virtual std::vector<double> getData()const;
    virtual std::vector<double> getDataStandardError()const;

    virtual double logDetCovData();


    virtual ABC_model& setData(const ABC_data& experimentalData);

    std::ostream& run(std::ostream& f,
                      const std::vector<std::pair <double,Parameters> >& pars)const;

    virtual std::ostream&
    run(std::ostream& f,const std::vector<Parameters >& pars,double quantil0, double quantil1)const;


    BayesIteration(const BayesIteration& other);
    /*

    friend void swap(BayesIteration& one, BayesIteration& other);

    BayesIteration& operator=(const BayesIteration& other);
*/
    BayesIteration();

    BayesIteration& getPosterior();

    BayesIteration& getPosterior(const Parameters& startingPoint, std::size_t numIterations);

    BayesIteration& getImportancePosterior(const Parameters& startingPoint, std::size_t numIterations);


    BayesIteration& getPosterior(const Parameters& startingPoint,double factor, std::size_t numSeeds,double probParChange);

    Parameters getEvidence(const Parameters& maximumPostLik, std::size_t num);

    Parameters getFIMbeta(const Parameters&par, double beta=1.0);
    Parameters getJJwbeta(const Parameters&par, double beta=1.0);

    std::vector<double> yfit(const BayesIteration& finalModel,const Parameters& parameters,double beta)const;

    std::pair<double,Parameters>
    anaLsMCMC(std::vector<std::pair<double,Parameters> > mcmc);


    Evidence
    anaDiMCMC(std::vector<std::pair<double,Parameters> > mcmc,
              const LevenbergMarquardtParameters& LM);

    std::vector<std::pair<double, Parameters> > runMCMC(LevenbergMarquardtParameters& LM,
                                                        std::size_t num);



    MCMCrun
    runMCMC(const Parameters& start,
            double beta, double factor,
            std::size_t num);


    void
    runMCMC(const Parameters& start,
            double factor,
            std::size_t num);

    void
    quantilesMCMC(MCMCrun mcmcrun);

    void
    runQuantilesMCMC(MCMCrun mcmcrun);

    MCMCrun
    runMCMC(BayesIteration finalModel,
         const Parameters& start,
                             double beta, double factor,
                             std::size_t num);


    MCMCrun
    initMCMC(const Parameters& start,
            double beta,
            std::size_t num);

    Evidence
    anaImportanceSampling(LevenbergMarquardtParameters& LM,
                                          std::size_t num);



    Parameters runImportanceSampling(Parameters& MAP,
                                                        std::size_t num,
                               double &factor);


    Evidence runThermodynamicIntegration(const Parameters &start,
                                         double factor,
                                     std::size_t num,
                                     std::size_t numBeta, double alfa);

    Evidence runThermodynamicIntegration(BayesIteration& finalModel,
                                         const Parameters &start, double factor,
                                     std::size_t num,
                                     std::size_t numBeta, double alfa);


    Parameters getHessian(const Parameters& MAP, double eps=1e-3);


    Parameters getHessianInterpol(const Parameters& MAP, double minP, double maxP);

    virtual std::ostream& put(std::ostream& s,const Parameters& parameters)const;



private:

    const ABC_model* m_;

    std::vector<const ABC_data*> data_;

    std::vector<Parameters> priors_;

    Parameters posterior_;
    std::size_t numSeeds_;

    std::vector<LevenbergMarquardtParameters> LM_;

    std::string filename_;


    };

#endif // BAYESITERATION_H
