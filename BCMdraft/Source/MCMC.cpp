#include <fstream>
#include <iostream>
#include "Include/MCMC.h"



MCMCrun
runMCMC(std::string filename,
    ABC_BCM* m,
        const VariablesValue &start,
        double beta,
        double factor,
        std::size_t num,
        const std::vector<const ExperimentDistribution *> &exper,
        double delta,
        double dt)
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

  VariablesNormalDistribution par=m->getJJwbeta(start,beta,exper,delta,dt);
  VariablesNormalDistribution parcenter=par;

  double postLogLik0=m->posteriorLogLikelihood(parcenter.center(),beta,exper,dt);
  double postLogLik=postLogLik0;

  double lik=exp(postLogLik-postLogLik0);
  std::ofstream f;

  f.open(filename.c_str(),std::ios_base::app);


  for (std::size_t i=0; i<num*1.1; i++ )
    {

      std::size_t j=0;
      std::size_t nsubsteps=par.size()*3;
      for (std::size_t ii=0; ii<par.size(); ii++)
        {
          VariablesValue nextPar=parcenter.randomSample(mcmc.factor);

          double  nextpostLogL=m->posteriorLogLikelihood(nextPar,beta,exper,dt);

          double nextLik=0;
          if (!nextpostLogL!=nextpostLogL)
            nextLik=exp(nextpostLogL-postLogLik0);

          double r=(1.0*rand())/RAND_MAX;

          if (nextLik/lik>r)
            {
              parcenter.setMeanValues(nextPar.values());
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

          double logL=m->logLikelihood(parcenter.center(),exper,dt);

          double ss=m->SumWeighedSquare(parcenter.center(),exper,dt);
          double sspar=m->SumWeighedSquareParameters(parcenter.center());
          mcmc.logL.push_back(logL);

          mcmc.p.push_back(parcenter.center());
          f<<i<<"\t"<<beta<<"\t"<<mcmc.factor<<"\t"<<ratio<<"\t"<<logL<<"\t"<<ss+sspar<<"\t"<<ss<<"\n";
          std::cerr<<i<<"\t"<<beta<<"\t"<<mcmc.factor<<"\t"<<ratio<<"\t"<<logL<<"\t"<<ss+sspar<<"\t"<<ss<<"\n";


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


