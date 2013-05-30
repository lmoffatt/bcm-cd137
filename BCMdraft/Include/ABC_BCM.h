#ifndef ABC_BCM_H
#define ABC_BCM_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <tuple>
#include <chrono>

#include "Include/Experiment.h"
#include "Include/Variables.h"
#include "Include/VariablesNormalDistribution.h"
//#include "Include/PosteriorLevenbergMarquardt.h"

class PosteriorLevenbergMarquardt;

/*
BCM
inputs:
model
arquitecture
1. initial equations
2. dynamic equations
3. observational equations
parameters
4. Parameters distributions

Experiment
5. Experimental conditions

6. Experimental observations


output:
1. e data
2. e parameter


LM
3.lista de maximos parciales con sus hessianos

MCMC
4. distribuciones a posteriori

Termodynamic integration
5. evidencia del modelo
*/
double randNormal();
double randNormal(double mean,double std);



class ABC_BCM
{
public:


  virtual std::string identifierName()const=0;
  ~ABC_BCM(){}
  virtual std::vector<double> internalStatesDerivative(const std::vector<double>& modelParameters,
                                                       const std::vector<double>& externalModelState,
                                                       const std::vector<double>& internalModelState)const=0;

  virtual std::vector<double> internalStatesStart(const std::vector<double>& log10modelParameters,
                                                  const std::vector<double>& externalModelState)const=0;



  virtual double observedVariablesEquations(const std::vector<double>& modelParameters,
                                            const std::vector<double>& externalModelState,
                                            const std::vector<double>& internalModelState,
                                            int observedIndex,
                                            double currentsum,
                                            double dt,
                                            bool finalcalculation)const=0;


  virtual double observedVariableIntegrationTime(const std::vector<double>& p,
                                                 const std::vector<double>& e,
                                                 int observedIndex)const=0;


  virtual std::vector<VariableInformation> internalStateVariables()const=0;
  virtual std::vector<VariableInformation> parameterVariables()const=0;
  virtual std::vector<VariableInformation> externalStateVariables()const=0;
  virtual std::vector<VariableInformation> observedStateVariables()const=0;

  virtual std::map<int,Transformation*> setTransformations()const=0;
  virtual std::map<int,RangeFormat*> setRangeFormat()const=0;


  virtual const Transformation *getTransform(int index)const;

  virtual const RangeFormat *getRangeFormat(int index)const;

  virtual const MultipleVariables* getInternalStateVariables()const;
  virtual const MultipleVariables* getParameterVariables()const;
  virtual const MultipleVariables* getExternalStateVariables()const;
  virtual const MultipleVariables* getObservedStateVariables()const;



  /* virtual std::vector<Experiment> getExperiments(std::vector<std::string> names)const;
*/



  virtual std::vector<double> RungeKutta4(const std::vector<double>& modelParameters,
                                          const std::vector<double>& externalModelState,
                                          const std::vector<double>& internalModelState,
                                          double h)const;




  virtual void push_backExperiment(const std::string& experimentIndentifier,
                                   const std::vector<VarValue> &externalState,
                                   const std::vector<VariableTimeDistribution> &data);


  virtual void push_backParameters(std::string priorIdentifier,
                                   std::vector< VariableDistribution > paramters);

  virtual void setPrior(std::string priorIdentifier);

  virtual const ExperimentDistribution *getExperiment(std::string name)const;

  virtual std::vector<const ExperimentDistribution *> getExperiments()const;

  virtual std::vector<const ExperimentDistribution*> getExperiments(std::vector<std::string> nameList)const;



  virtual const VariablesNormalDistribution *getPrior()const;

  virtual std::vector<double> yerr(const std::vector<const ExperimentDistribution*> exp) const;


  virtual std::vector<double> y(const std::vector<const ExperimentDistribution *> exp) const;

  virtual std::vector<double> Ty(const std::vector<const ExperimentDistribution *> exp) const;

  virtual std::vector<double> w(const std::vector<const ExperimentDistribution *> exp) const;

  virtual std::vector<double> yfit(const VariablesValue& parameters,
                                   const ExperimentSchedule& exper, double dt) const;

  virtual std::vector<const ExperimentDistribution *> simulate(
      const VariablesValue& parameters,
      const std::vector<const ExperimentDistribution *> &exper,
      double dt,
      double measureNoisefactor) const;



  virtual std::vector<double> yfit(const VariablesValue& par,
                                   const std::vector<const ExperimentDistribution *> &exp,
                                   double dt) const;

  virtual std::vector<double> Tyfit(const VariablesValue& par,
                                    const std::vector<const ExperimentDistribution *> &exp,
                                    double dt) const;
  virtual double SSdata(const std::vector<double>& yf,
                        const std::vector<double>& y,
                        const std::vector<double>& w) const;


  virtual double SSpar(const VariablesValue &par) const;



  /*virtual std::vector<std::vector<double>> getJdata(const Parameters& par,
                                                       const std::vector<Experiment >& exp,
                                                       double delta,
                                                       double dt)const;



  virtual std::map<double,std::vector<PosteriorLevenbergMarquardt> >
  findMAPs(const std::vector<Experiment> exp,double searchRadius, std::size_t numSeeds, std::size_t numIter);




*/

  virtual VariablesNormalDistribution getJJwbeta(const VariablesValue& par,
                                                 double beta,
                                                 const std::vector<const ExperimentDistribution *> &exp,
                                                 double dx,
                                                 double dt);


  virtual double posteriorLogLikelihood(const VariablesValue& p,
                                        double beta,
                                        const std::vector<const ExperimentDistribution *> &exp,
                                        double dt);


  virtual double logLikelihood(const VariablesValue& p,
                               const std::vector<const ExperimentDistribution*> &exp,
                               double dt);


  virtual double SumWeighedSquare(const VariablesValue& p,
                                  const std::vector<const ExperimentDistribution *> &exp,
                                  double dt);

  virtual double SumWeighedSquareParameters(const VariablesValue& p);


  virtual std::map<std::pair<double, double>, std::vector<const SingleVariable *> >
  getIntegrMeasures(const VariablesValue &par,
                    const ExperimentSchedule &anExp)const;


  virtual void init();
static const unsigned long seed;
private:
  MultipleVariables* externalVariables_;
  MultipleVariables* internalVariables_;
  MultipleVariables* observVariables_;
  MultipleVariables* paramVariables_;

  std::map<int,Transformation*> trnas_;
  std::map<int,RangeFormat*> form_;

  std::map<std::string,ExperimentDistribution*> expM_;
  std::map<std::string,VariablesNormalDistribution*> parM_;
  std::string prior_;


};


std::vector<double>& operator+=(std::vector<double>& x, const std::vector<double>& y);


std::vector<double> operator+(const std::vector<double>& x, const std::vector<double>& y);


std::vector<double> operator*(const std::vector<double>& x, double h);


inline std::ostream& operator<<(std::ostream& s, std::vector<std::string> v)
{
  for (std::string n:v)
    s<<n<<"\t";
  return s;
}


class ABC_BCM_Table
{

  virtual ABC_BCM* getBCM(const std::string& name)=0;

};




/*

class ABC_Model
{

  virtual std::vector<double> yfit(std::vector<double>& )=0;
  virtual std::vector<double> ymeasured()=0;
  virtual std::vector<double> yerr()=0;
  virtual double SSwdata()=0;
  virtual double SSwpar()=0;
};

class ABC_ODEModel:public ABC_Model
{
  std::set<Variable> internal_variables() const=0;
  std::set<Variable> observable_variables()const=0;
  std::set<Variable> external_variables()const=0;

  virtual std::vector<double> Start()=0;

  virtual std::vector<double> Derivative(double t,std::vector<double> y)=0;

  virtual double Observation(double t,std::string obsvarname)=0;




  Parameters parameters();

};



class Experiment
{
  double externalVariable(std::string name,double time)=0;

};

class Result
{
  std::tuple<time,observation,stderror> measure(std::size_t i);
  std::size_t n()const;

};


*/
#endif // ABC_BCM_H
