#ifndef SAMPLEDDISTRIBUTION_H
#define SAMPLEDDISTRIBUTION_H


#include "Include/Variables.h"


class SampledDistribution
{
public:
  virtual VariablesValue randomSample()const;

  virtual double logProbability(const VariablesValue& sample)const;

  const MultipleVariables* variables()const;

  virtual std::size_t size() const;




  // virtual double mean(ABC_function* f);

  //  virtual double variance(ABC_function* f);

  virtual const std::vector<double>& tstd()const;

  virtual const std::vector<double>& w()const;

  virtual double KLdivergence()const;

  virtual  const std::vector< std::vector <double> >& getCovariance()const;

  SampledDistribution(const MultipleVariables *var,
                              const std::vector<double>& tmean,
                              const std::vector< std::vector <double> >& cov );

  SampledDistribution(const MultipleVariables* var,
                              const std::vector<double>& tmean,
                              const std::vector< double >& tstd);

  SampledDistribution(const VariablesValue& par,
                              const std::vector< std::vector <double> >& cov );


  SampledDistribution(const MultipleVariables* var);


public:
  friend std::istream& operator>>(std::istream& s, SampledDistribution& v);

  VariablesValue var_;
  std::vector< std::vector <double> > cov_;
  std::vector<double> tStd_;  // not in dB
  std::vector< std::vector <double> > cho_;
  std::vector<double> w_;
};


std::ostream& operator<<(std::ostream& s, const SampledDistribution& v);

  std::istream& operator>>(std::istream& s, SampledDistribution& v);



#endif // SAMPLEDDISTRIBUTION_H
