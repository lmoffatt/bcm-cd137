#ifndef BCM_CD137_H
#define BCM_CD137_H

#include "Include/ABC_BCM.h"

class BCM_CD137:public ABC_BCM
{
public:

  virtual std::string identifierName()const;
  virtual std::vector<double> internalStatesDerivative(const std::vector<double>& modelParameters,
                                                       const std::vector<double>& externalModelState,
                                                       const std::vector<double>& internalModelState)const;

  virtual std::vector<double> internalStatesStart(const std::vector<double>& log10modelParameters,
                                                  const std::vector<double>& externalModelState)const;



  virtual double observedVariablesEquations(const std::vector<double>& modelParameters,
                                            const std::vector<double>& externalModelState,
                                            const std::vector<double>& internalModelState,
                                            int observedIndex,
                                            double currentsum,
                                            double dt,
                                            bool finalcalculation)const;



  virtual double observedVariableIntegrationTime(const std::vector<double>& p,
                                                 const std::vector<double>& e,
                                                 int observedIndex)const;




  virtual std::vector<VariableInformation> internalStateVariables()const;
  virtual std::vector<VariableInformation> parameterVariables()const;
  virtual std::vector<VariableInformation> externalStateVariables()const;
  virtual std::vector<VariableInformation> observedStateVariables()const;

  virtual std::map<int,Transformation*> setTransformations()const;
  virtual std::map<int,RangeFormat*> setRangeFormat()const;


  BCM_CD137();
};


#endif // BCM_CD137_H
