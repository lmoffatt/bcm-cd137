#ifndef VARIABLESNORMALDISTRIBUTION_H
#define VARIABLESNORMALDISTRIBUTION_H


#include <ostream>
#include <vector>
#include <cmath>
#include "Include/Variables.h"

double randNormal(double mean,double stddev);
double randNormal();


class RangeFormat
{
public:
  virtual std::pair<double,double> getTransformedMoments(const Transformation*,double first,double second)const=0;

};

class Range_1S:public RangeFormat
{
public:
  std::pair<double,double> getTransformedMoments(const Transformation* T,double first,double second)const;
};

class Mean_SD:public RangeFormat
{
public:
  std::pair<double,double> getTransformedMoments(const Transformation* T,double first,double second)const;
};

class Mean_DB:public RangeFormat
{
public:
  std::pair<double,double> getTransformedMoments(const Transformation* T,double first,double second)const;
};






class VariablesNormalDistribution
{
public:
  virtual VariablesValue randomSample()const;

  virtual VariablesValue randomSample(double factor)const;

  virtual const VariablesValue &center()const;


  virtual std::size_t size() const;

  //virtual double probability(const Variables& inferiorCorner,const Variables& superiorCorner);

  // virtual double mean(ABC_function* f);

  //  virtual double variance(ABC_function* f);

  virtual const std::vector<double>& tstd()const;

  virtual const std::vector<double>& w()const;


  const MultipleVariables* variables()const;

  virtual  const std::vector< std::vector <double> >& getCovariance()const;

  virtual VariablesNormalDistribution& setTmeanValues(const std::vector<double>& tmean);
  virtual VariablesNormalDistribution& setMeanValues(const std::vector<double>& mean);
  virtual VariablesNormalDistribution& setStdValues(const std::vector<double>& tstd);


  VariablesNormalDistribution(const MultipleVariables *var,
                              const std::vector<double>& tmean,
                              const std::vector< std::vector <double> >& cov );

  VariablesNormalDistribution(const MultipleVariables* var,
                              const std::vector<double>& tmean,
                              const std::vector< double >& tstd);

  VariablesNormalDistribution(const VariablesValue& par,
                              const std::vector< std::vector <double> >& cov );


  VariablesNormalDistribution(const MultipleVariables* var);


public:
  friend std::istream& operator>>(std::istream& s, VariablesNormalDistribution& v);

  VariablesValue var_;
  std::vector< std::vector <double> > cov_;
  std::vector<double> tStd_;  // not in dB
  std::vector< std::vector <double> > cho_;
  std::vector<double> w_;
};


std::ostream& operator<<(std::ostream& s, const VariablesNormalDistribution& v);

  std::istream& operator>>(std::istream& s, VariablesNormalDistribution& v);



#endif // VARIABLESNORMALDISTRIBUTION_H
