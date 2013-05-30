#ifndef EXPERIMENT_H
#define EXPERIMENT_H
#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <set>

#include "Include/Variables.h"
#include "Include/VariablesNormalDistribution.h"



class ABC_BCM;
class ExperimentSchedule
{
public:
  std::string nameId()const;
  const VariablesValue& ExternalState()const;
  const MultipleVariables& measuredVariables()const;
  const std::vector<double>& measuredTimes()const;

  double tmax()const;
  ExperimentSchedule(const ABC_BCM *model,
                     const std::string& experimentIndentifier,
                     const std::vector<VarValue> &externalState,
                     const std::vector<VariableTimeDistribution>& data);

private:
  VariablesValue getExternalState(std::vector<VarValue> externalSate);
  MultipleVariables getMeasuredValues(const std::string &experimentIdentifier,
                                      std::vector<VariableTimeDistribution> data);

  std::vector<double> getMeasuredTimes(std::vector<VariableTimeDistribution> data);

  std::string name_;
  const ABC_BCM* model_;
  VariablesValue externalVar_;
  MultipleVariables measuredVars_;
  std::vector<double> measuredTimes_;

  double tmax_;
};

class ExperimentDistribution
{
public:
  const std::vector<double> &w()const;
  const std::vector<double> &y()const;
  const std::vector<double> &Ty()const;
  const std::vector<double> &yerr()const;
  const ExperimentSchedule& schedule()const;
  const VariablesNormalDistribution& distr()const;

  virtual ExperimentDistribution& setTmeanValues(const std::vector<double> tmean);
  virtual ExperimentDistribution& setMeanValues(const std::vector<double> mean);


  ExperimentDistribution(const ABC_BCM *model,
                         const std::string& experimentIndentifier,
                         const std::vector<VarValue>& externalState,
                         const std::vector<VariableTimeDistribution> &data);


private:
  VariablesNormalDistribution getYvalues(const ABC_BCM* model,
                                         const ExperimentSchedule& sche,
                                         std::vector<VariableTimeDistribution> data);
  ExperimentSchedule sche_;
  VariablesNormalDistribution distr_;
};

class ExperimentValue
{
public:
  const VariablesValue& y()const;
  const ExperimentSchedule* schedule()const;
};

std::ostream& operator<<(std::ostream& s, const ExperimentDistribution& exp);

std::ostream& operator<<(std::ostream& s, const std::vector<const ExperimentDistribution*>& exp);


/*

class Measure
{
   int i_;
  double val_;
  double err_;
public:
  Measure(int indexObser ,double v,double err):
    i_(indexObser),
    val_(v),
    err_(err){}
  int index()const{return i_;}
  double value()const{return val_;}
  double se()const{return err_;}
};




class Experiment{
      std::string id;
      MultipleVariables eState;
      std::map<double,std::vector<SingleVariableNormalDistribution> > obs_;
      std::set<int> measVar_;
      std::map<int,std::map<double,Measure>> di_;
      std::vector<double> y_;
      std::vector<double> yerr_;
      std::vector<double> w_;
      double tmax_;
public:
     Experiment(const std::string& experimentIndentifier,
                std::vector<std::pair<int,double>> externalState,
                std::vector<std::tuple<int,double,double,double>> data);





     void setyfit(const std::vector<double>& newy);

     Experiment()=default;
     Experiment(const Experiment& other)=default;
     Experiment& operator=(const Experiment& other)=default;

     Experiment(Experiment&& other)=default;
     Experiment& operator=(Experiment&& other)=default;

     const std::set<int> &getMeasVarSet()const;

     std::string name()const;
     const VariablesValue& ExternalState()const;

     const std::map<double,SingleVariableNormalDistribution>& timeMeasure(int index)const;

     std::size_t size()const;
     const std::vector<double>& y()const;
     const std::vector<double>& yerr()const;
     const std::vector<double>& w()const;


     double tmax()const;


};
*/
#endif // EXPERIMENT_H
