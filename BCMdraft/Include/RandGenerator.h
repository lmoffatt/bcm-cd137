#ifndef RANDGENERATOR_H
#define RANDGENERATOR_H
#include <random>


class RandomGenerator
{
public:
double randNormal();
double randNormal(double mean,double std);
double rand();
static const unsigned long seed;
private:
 static std::mt19937 random_engine_;
 std::normal_distribution<> normal_;
 std::uniform_real_distribution<> uniform_;
};

#endif // RANDGENERATOR_H
