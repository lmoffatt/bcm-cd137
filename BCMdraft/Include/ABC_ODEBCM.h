#ifndef ABC_ODEBCM_H
#define ABC_ODEBCM_H

#include <vector>
#include <string>
#include <map>







class ABC_ODEBCM
{
public:


private:

  virtual void pushExperiment(const std::vector<std::map<double, double> >& externalStateEvolutions,
                             const std::vector<std::map<double, double> >& observedVariablevolution);

  virtual void setPriors(const std::vector<std::pair<double,double>>& mean);











};

#endif // ABC_ODEBCM_H
