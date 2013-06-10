#ifndef MCMC_H
#define MCMC_H

#include "Include/ABC_BCM.h"






struct MCMCrun
{
  std::vector<VariablesValue> p;
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


MCMCrun
runMCMC(std::string filename,
    ABC_BCM* m,
        const VariablesValue& start,
        double beta,
        double factor,
        std::size_t num,
        const std::vector<const ExperimentDistribution *> &exper,
        double delta,
        double dt);








#endif // MCMC_H
