#include <iostream>
#include <numeric>
#include <fstream>
#include <sstream>

#include "Include/BCM_CD137.h"
#include "Include/MatrixInverse.h"
#include "Include/PosteriorLevenbergMarquardt.h"

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

            if (argc>2)
            {
                std::string seg=argv[2];
                std::stringstream ss(seg);
                ss>>nseeds;
            }

            for(std::size_t i=0; i<nseeds;i++)
            {
                PosteriorLevenbergMarquardt LM(&m,m.getPrior()->randomSample(0.1),e,0.01);
                LM.optimize(10,1000);

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
    }
    return 0;

}





