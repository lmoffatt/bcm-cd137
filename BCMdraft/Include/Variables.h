#ifndef VARIABLES_H
#define VARIABLES_H
#include <vector>
#include <string>
#include <cmath>


inline std::vector<double> pow(const std::vector<double>& x,int n)
{
  std::vector<double> res(x.size());
  for (unsigned i=0; i<res.size(); i++)
    res[i]=pow(x[i],n);
  return res;
}
enum TRANSF{LINEAR,LOG10,LOG10RATIO,LOG10PERCENT};
enum FORMAT{RANGE_1S,MEAN_SD,MEAN_DB};

struct VariableInformation
{
  int index;
  std::string shortName;
  std::string units;
  int TransformationIndex;
 };

struct VarValue
{
  int index;
 double value;
};

struct VariableDistribution
{
  int index;
 double first;
 double second;
 int FormatIndex;
};
struct VariableTimeDistribution
{
  int index;
  double time;
 double first;
 double second;
 int FormatIndex;
};


class Transformation
{
public:
  virtual std::string myClass()const=0;
  virtual double eval(double x)const=0;
  virtual double inverse(double x)const=0;

};

class Log10Tranformation:public Transformation
{
  std::string myClass()const;
  double eval(double x)const;
  double inverse(double x)const;

  static std::string ClassName();
};

class Log10RatioTranformation:public Transformation
{
  std::string myClass()const;
  double eval(double x)const;
  double inverse(double x)const;

  static std::string ClassName();
};
class Log10PercentTranformation:public Transformation
{
  std::string myClass()const;
  double eval(double x)const;
  double inverse(double x)const;

  static std::string ClassName();
};


class LinearTranformation:public Transformation
{
  std::string myClass()const;
  double eval(double x)const;
  double inverse(double x)const;

  static std::string ClassName();
};

class MultipleVariables;
class ABC_BCM;

class SingleVariable
{
public:
  const MultipleVariables& VariableSet()const;
  int index()const;
  const std::string& name()const;
  const std::string& unit() const;
  const Transformation *transfromation()const;

  SingleVariable(const MultipleVariables* parent,
                 int myindex,
                 const std::string& myname,
                 const std::string& myunit,
                 const Transformation* myTransformation);



  SingleVariable(const SingleVariable&)=delete;
  SingleVariable& operator=(const SingleVariable&)=delete;

  SingleVariable()=delete;


private:
  const MultipleVariables* parent_;
  int index_;
   std::string name_;
   std::string unit_;
  const Transformation* transfromation_;
};


class MultipleVariables
{
public:
  std::string variablesSetName()const;

  const std::vector<const SingleVariable *> &var()const;
  const SingleVariable *var(unsigned i)const;

  std::vector<double> transform(const std::vector<double> val)const;

  std::vector<double> trans_inverse(const std::vector<double> &val)const;

  std::size_t size()const;

  ~MultipleVariables();

  const ABC_BCM* model()const;

  MultipleVariables(const ABC_BCM* model,
                   const std::string name,
                    const std::vector<const SingleVariable *> &vars);
  MultipleVariables(const ABC_BCM* model,
                    const std::string name,
                    const std::vector <VariableInformation> & initialList);


private:
  std::vector<const SingleVariable*> buildList(
      const std::vector<VariableInformation >& initialList);
  const ABC_BCM* model_;
  std::string variableName_;
   std::vector<const SingleVariable*> var_;
   bool owner_;
};


class VariablesValue
{
public:
  const MultipleVariables* variables()const;

  const std::vector<double> &values()const;

  const std::vector<double> &tvalues()const;

  std::size_t size()const;

  VariablesValue& setValues(const std::vector<double>& newValues);
  VariablesValue& setTransfValues(const std::vector<double> &newTransValues);

  double getTvalue(unsigned i)const;
  VariablesValue& setTvalue(unsigned i, double x);

  VariablesValue& addTransfDx(std::size_t i, double dx);


  VariablesValue(const MultipleVariables *var);
private:
  const MultipleVariables* variables_;
  std::vector<double> values_;
  std::vector<double> tvalues_;

};



std::ostream& operator <<(std::ostream& s, const VariablesValue& val);
















#endif // VARIABLES_H
