#include <limits>
#include <chrono>
#include "Include/RandGenerator.h"

const unsigned long RandomGenerator::seed = std::chrono::high_resolution_clock::now()
    .time_since_epoch().count();



std::mt19937 RandomGenerator::random_engine_=std::mt19937(seed);


double RandomGenerator::randNormal()
{
 return normal_(random_engine_);
}

double RandomGenerator::randNormal(double mean,double std)
{
  return mean+randNormal()*std;
}

double RandomGenerator::rand()
{
 return uniform_(random_engine_);
}
