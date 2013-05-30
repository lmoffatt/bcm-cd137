#include <cmath>
#include <limits>
#include <random>
#include <chrono>
#include <sstream>

#include "Include/VariablesNormalDistribution.h"
#include "Include/MatrixInverse.h"
#include "Include/ABC_BCM.h"


std::pair<double,double> Range_1S::getTransformedMoments(const Transformation* T,double first,double second)const
{
  double tfirst=T->eval(first);
  double tsecond=T->eval(second);
  return {0.5*(tfirst+tsecond),0.5*(tsecond-tfirst)};
}

std::pair<double,double> Mean_SD::getTransformedMoments(const Transformation* T,double first,double second)const
{
  double tmax=T->eval(first+second);
  double tmin=T->eval(first-second);
  if (std::isfinite(tmax)&&std::isfinite(tmin))
    return {0.5*(tmax+tmin),0.5*(tmax-tmin)};
else
{
double tmean=T->eval(first);
if (std::isfinite(tmax))
return {tmean,tmax-tmean};
else
return{tmean,tmean-tmin};
}
}
std::pair<double,double> Mean_DB::getTransformedMoments(const Transformation* T,double first,double second)const
{
  return {T->eval(first),10.0*second};

}


VariablesNormalDistribution& VariablesNormalDistribution::setTmeanValues(const std::vector<double>& tmean)
{
  var_.setTransfValues(tmean);
  return *this;
}

VariablesNormalDistribution& VariablesNormalDistribution::setMeanValues(const std::vector<double>& mean)
{
  var_.setValues(mean);
  return *this;
}

VariablesNormalDistribution& VariablesNormalDistribution::setStdValues(const std::vector<double>& tstd)
{
  tStd_=tstd;
  cov_.clear();
  cho_.clear();
}


const VariablesValue &VariablesNormalDistribution::center()const
{
  return var_;
}

VariablesValue VariablesNormalDistribution::randomSample(double factor)const
{
  double chi2z=0;
  double sumz=0;
  std::vector<double> val(size());
  if (cho_.empty())
    {
      for (std::size_t i=0; i<size();i++)

        {
          val[i]=randNormal(var_.tvalues()[i],tStd_[i]*factor);

        }
    }
  else
    {
      std::vector<double> z(size());
      for (std::size_t i=0; i<size();i++)
        {
          z[i]=randNormal()*factor;
          chi2z+=z[i]*z[i];
          sumz+=z[i];

        }      for (std::size_t i=0; i<size();i++)
        {

          val[i]=var_.tvalues()[i];
          for (std::size_t j=0; j<i+1;j++)
            {
              val[i]+=cho_[i][j]*z[j];
            }

        }

    }

  VariablesValue s(variables());
  return s.setTransfValues(val);

}

const MultipleVariables* VariablesNormalDistribution::variables()const
{
  return var_.variables();
}






VariablesValue VariablesNormalDistribution::randomSample()const
{
  return randomSample(1.0);
}



const std::vector<double> &VariablesNormalDistribution::tstd()const
{
  return tStd_;
}

const std::vector<double> &VariablesNormalDistribution::w()const
{
  return w_;
}

const std::vector<std::vector<double> > &VariablesNormalDistribution::getCovariance()const
{
  return cov_;
}

std::size_t VariablesNormalDistribution::size() const
{
  return var_.size();
}

VariablesNormalDistribution::VariablesNormalDistribution(
    const MultipleVariables* var,
    const std::vector<double>& tmean,
    const std::vector< std::vector <double> >& cov ):
  var_((VariablesValue(var)).setTransfValues(tmean)),
  cov_(cov),
  tStd_(sqrt(diag(cov))),
  cho_(chol(cov)),
  w_(pow(tStd_,-2))
{}

VariablesNormalDistribution::VariablesNormalDistribution(const MultipleVariables *var,
                                                         const std::vector<double>& tmean,
                                                         const std::vector <double> & tstd ):
  var_((VariablesValue(var)).setTransfValues(tmean)),
  cov_(std::vector<std::vector<double> > ()),
  tStd_(tstd),
  cho_(),
  w_(pow(tStd_,-2))
{}


VariablesNormalDistribution::VariablesNormalDistribution(const VariablesValue &par,
                                                         const std::vector< std::vector <double> >& cov ):
  var_(par),
  cov_(cov),
  tStd_(sqrt(diag(cov))),
  cho_(chol(cov)),
  w_(pow(tStd_,-2))
{}


std::ostream& operator<<(std::ostream& s, const VariablesNormalDistribution& v)
{
  s<<"Variables_Normal_Distriubion \n";
  s<<"Begin \n";

  s<<"Bayesian_Model \n";
  s<<v.variables()->model()->identifierName()<<"\n";
  s<<"name\t meanValue\t x-se  \t x+se \t transf \n";
  for (std::size_t i=0; i<v.size(); i++)
    {

      s<<v.variables()->var(i)->name()<<"\t";
      s<<v.center().values()[i]<<"\t";
      s<<v.variables()->var(i)->transfromation()->inverse(
           v.center().tvalues()[i]-v.tstd()[i])<<"\t";
      s<<v.variables()->var(i)->transfromation()->inverse(
           v.center().tvalues()[i]+v.tstd()[i])<<"\t";
      s<<v.variables()->var(i)->transfromation()->myClass()<<std::endl;

    }
  s<<"End\n";

  return s;

}

std::istream& operator>>(std::istream& s, VariablesNormalDistribution& v)
{
  std::string line;
  std::getline(s,line);//<<"Variables_Normal_Distriubion \n";
  std::getline(s,line);//  s<<"Begin \n";

  std::getline(s,line);//  s<<"Bayesian_Model \n";
  std::getline(s,line);//  s<<v.variables()->model()->identifierName()<<"\n";
  std::getline(s,line);//  s<<"name\t meanValue\t x-se  \t x+se \t transf \n";
  std::vector<double> meanValue(v.size());
  std::vector<double> std(v.size(),0);

  for (std::size_t i=0; i<v.size(); i++)
    {
      std::getline(s,line);
      std::stringstream ss(line);
      std::string name;
      ss>>name;
      double mean;
      ss>>mean;
      double x_se;
      ss>> x_se;
      double xpse;
      ss>>xpse;
      std::string trans;
      ss>>trans;

      for (std::size_t i=0; i<v.size(); ++i)
        {
          if (v.variables()->var(i)->name()==name)
            {
              meanValue[i]=mean;
              if (v.variables()->var(i)->transfromation()->myClass()==trans)
                {
                  double tmean=v.variables()->var(i)->transfromation()->eval(mean);
                  std[i]=v.variables()->var(i)->transfromation()->eval(x_se)-tmean;
                }

            }
        }
    }
  std::getline(s,line); //  s<<"End\n";
  v.setMeanValues(meanValue);
  v.setStdValues(std);

  return s;

}


