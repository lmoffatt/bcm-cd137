#include <iostream>
#include <numeric>
#include <fstream>
#include <sstream>

#include "Include/BCM_CD137.h"
#include "Include/MatrixInverse.h"
#include "Include/PosteriorLevenbergMarquardt.h"
#include "Include/MCMC.h"

int main(int argc, char **argv)
{
  BCM_CD137 m;

  if (argc>1)
    {
      std::string  command=argv[1];
      if (command=="optimize")
        {
          auto e=m.getExperiments({"media","mtb","CD137block"});
          std::size_t nseeds=1;
          std::size_t niter=100;
          double radius=0.5;
          if (argc>2)
            {
              std::string siter=argv[2];
              std::stringstream ss(siter);
              ss>>niter;
            }
          if (argc>3)
            {
              std::string seg=argv[3];
              std::stringstream ss(seg);
              ss>>radius;
            }



          for(std::size_t i=0; i<nseeds;i++)
            {
              PosteriorLevenbergMarquardt LM(
                    &m,m.getPrior()->randomSample(radius),e,0.01);
              LM.optimize(niter,10000);

              std::stringstream ss1;
              ss1<<LM.SS();
              ss1<<"optimize";
              ss1<<ABC_BCM::seed;
              ss1<<".txt";
              std::string filename=ss1.str();
              std::ofstream ff;
              ff.open(filename.c_str()  , std::fstream::out);



              ff<<LM;
              ff<<LM.OptimParameters();
              ff<<m.simulate(LM.OptimParameters().center(),e,0.01,0);



              ff.flush();

            }


        }

      if (command=="reoptimize")
        {
          std::string fnpar="testParameters.txt";
          if (argc>2)
            {
              fnpar=argv[2];
            }
          std::fstream f;
          f.open(fnpar);
          VariablesNormalDistribution testP(m.getPrior()->variables());
          f>>testP;
          auto e=m.getExperiments({"media","mtb","CD137block"});

          std::size_t nseeds=1;
          std::size_t niter=100;
          double radius=0.5;
          if (argc>3)
            {
              std::string siter=argv[3];
              std::stringstream ss(siter);
              ss>>niter;
            }
          if (argc>4)
            {
              std::string seg=argv[4];
              std::stringstream ss(seg);
              ss>>radius;
            }



          for(std::size_t i=0; i<nseeds;i++)
            {
              PosteriorLevenbergMarquardt LM(
                    &m,testP.randomSample(radius),e,0.01);
              LM.optimize(niter,10000);

              std::stringstream ss1;
              ss1<<LM.SS();
              ss1<<"reoptimize";
              ss1<<ABC_BCM::seed;
              ss1<<".txt";
              std::string filename=ss1.str();
              std::ofstream ff;
              ff.open(filename.c_str()  , std::fstream::out);



              ff<<LM;
              ff<<LM.OptimParameters();
              ff<<m.simulate(LM.OptimParameters().center(),e,0.01,0);



              ff.flush();

            }


        }
      if (command=="mcmc")

        {
          double beta=1;
          double factor=0.001;
          double delta=1e-7;
          double dt=0.01;
          std::string fnpar="testParameters.txt";
          if (argc>2)
            {
              fnpar=argv[2];
            }




          std::fstream f;
          f.open(fnpar);
          VariablesNormalDistribution start(m.getPrior()->variables());
          f>>start;
          auto e=m.getExperiments({"media","mtb","CD137block"});

          std::stringstream ss1;
          ss1<<m.SumWeighedSquare(start.center(),e,dt)+m.SumWeighedSquareParameters(start.center());
          ss1<<"MCMC";
          ss1<<ABC_BCM::seed;
          ss1<<".txt";
          std::string filename=ss1.str();

          std::size_t niter=100;
          if (argc>3)
            {
              std::string siter=argv[3];
              std::stringstream ss(siter);
              ss>>niter;
            }
          if (argc>4)
            {
              std::string seg=argv[4];
              std::stringstream ss(seg);
              ss>>factor;
            }


          MCMCrun mcmc=
              runMCMC(filename,&m,start.center(),beta,factor,niter,e,delta,dt);





        }
      else if (command=="simulate"){
          std::string fnpar="testParameters.txt";
          double dt=0.01;

          if (argc>2)
            {
              fnpar=argv[2];
            }
          std::fstream f;
          f.open(fnpar);
          VariablesNormalDistribution testP(m.getPrior()->variables());
          auto e=m.getExperiments({"media","mtb","CD137block"});
          f>>testP;
          auto sim=m.simulate(testP.center(),e,dt,0.0);


          std::string fnout="simulate"+fnpar;

          std::ofstream fout;
          fout.open(fnout);
          fout<<sim;
          for (std::size_t i=0; i< sim.size(); i++)
            delete sim[i];





        }
    }
  return 0;

}





