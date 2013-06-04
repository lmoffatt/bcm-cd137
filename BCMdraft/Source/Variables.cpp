#include "Include/Variables.h"
#include "Include/ABC_BCM.h"



std::string Log10Tranformation::myClass()const
{
  return ClassName();
}
double Log10Tranformation::eval(double x)const
{
  return std::log10(x);
}
double Log10Tranformation::inverse(double x)const
{
  return std::pow(10.0,x);
}

std::string Log10Tranformation::ClassName()
{
  return "Log10";
}


std::string Log10RatioTranformation::myClass()const
{
  return ClassName();
}
double Log10RatioTranformation::eval(double x)const
{
  return std::log10(x/(1.0-x));
}
double Log10RatioTranformation::inverse(double x)const
{
  double r= std::pow(10.0,x);
  return r/(r+1.0);
}

std::string Log10RatioTranformation::ClassName()
{
  return "Log10Ratio";
}



std::string Log10PercentTranformation::myClass()const
{
  return ClassName();
}
double Log10PercentTranformation::eval(double x)const
{
  return std::log10(x/(100.0-x));
}
double Log10PercentTranformation::inverse(double x)const
{
  double r= std::pow(10.0,x);
  return r/(r+1.0)*100.0;
}

std::string Log10PercentTranformation::ClassName()
{
  return "Log10Percent";
}


std::string LinearTranformation::myClass()const
{
  return ClassName();
}
double LinearTranformation::eval(double x)const
{
  return x;
}
double LinearTranformation::inverse(double x)const
{
  return x;
}

std::string LinearTranformation::ClassName()
{
  return "Linear";
}


const ABC_BCM* MultipleVariables::model()const
{
  return model_;
}

int SingleVariable::index()const
{
  return index_;
}

const std::string& SingleVariable::name()const{
  return name_;
}
const std::string& SingleVariable::unit() const{
  return unit_;
}
const Transformation* SingleVariable::transfromation()const{
  return transfromation_;

}

SingleVariable::SingleVariable(const MultipleVariables *parent,
                               int myindex,
                               const std::string& myname,
                               const std::string& myunit,
                               const Transformation *myTransformation):
  parent_(parent),
  index_(myindex),
  name_(myname),
  unit_(myunit),
  transfromation_(myTransformation)
{}




std::string MultipleVariables::variablesSetName()const{
  return variableName_;

}

const std::vector<const SingleVariable*>& MultipleVariables::var()const{
  return var_;

}
const SingleVariable* MultipleVariables::var(unsigned i)const{
  return var_[i];
}

std::vector<double> MultipleVariables::transform(const std::vector<double> val)const
{
  std::vector<double> res(var_.size());
  for (std::size_t i=0; i<res.size(); ++i)
    {
      res[i]=var(i)->transfromation()->eval(val[i]);
    }
  return res;
}

std::vector<double> MultipleVariables::trans_inverse(const std::vector<double> &val)const{
  std::vector<double> res(var_.size());
  for (std::size_t i=0; i<res.size(); ++i)
    {
      res[i]=var(i)->transfromation()->inverse(val[i]);
    }
  return res;
}
std::size_t MultipleVariables::size()const{

  return var_.size();
}

MultipleVariables::MultipleVariables(const ABC_BCM *model,
                                     const std::string name,
                                     const std::vector<const SingleVariable*>& vars):
  model_(model),
  variableName_(name),
  var_(vars),
  owner_(false){}

MultipleVariables::MultipleVariables(const ABC_BCM* model,
                                     const std::string name,
                                     const std::vector<VariableInformation >& initialList):
  model_(model),
  variableName_(name),
  var_(buildList(initialList)),
  owner_(true){}

std::vector<const SingleVariable*> MultipleVariables::buildList(const std::vector<VariableInformation> &initialList)
{
  std::vector<const SingleVariable*> var(initialList.size());
  for (std::size_t i=0; i<initialList.size(); ++i)
    {
      const VariableInformation* va=&initialList[i];
      const Transformation* tr=model_->getTransform(va->TransformationIndex);
      SingleVariable* v=new  SingleVariable(this,va->index,va->shortName,va->units,tr);
      var[va->index]=v;
   }
  return var;
}

MultipleVariables::~MultipleVariables(){
  if (owner_)
    {
      for (const SingleVariable* ptr: var_)
          delete ptr;
    }
}



std::size_t VariablesValue::size()const
{
  return variables()->size();
}

VariablesValue::VariablesValue(const MultipleVariables* var):
  variables_(var),
  values_(),
  tvalues_()
{}


VariablesValue& VariablesValue::setValues(const std::vector<double> &newValues)
{
  values_=newValues;
  tvalues_=variables()->transform(newValues);
  return *this;
}

VariablesValue& VariablesValue::setTransfValues(const std::vector<double>& newTransValues)
{
  values_=variables()->trans_inverse(newTransValues);
  tvalues_=newTransValues;
  return *this;
}


VariablesValue& VariablesValue::addTransfDx(std::size_t i, double dx)
{
  tvalues_[i]+=dx;
  values_[i]=variables()->var(i)->transfromation()->inverse(tvalues_[i]);
return *this;
}


double VariablesValue::getTvalue(unsigned i)const
{
  return tvalues_[i];
}
VariablesValue& VariablesValue::setTvalue(unsigned i, double x)
{
  tvalues_[i]=x;
  values_[i]=variables()->var(i)->transfromation()->inverse(tvalues_[i]);
return *this;
}


const std::vector<double>& VariablesValue::values()const
{
  return values_;
}


const MultipleVariables *VariablesValue::variables()const
{
  return variables_;
}

const std::vector<double>& VariablesValue::tvalues()const
{
  return tvalues_;
}

std::ostream& operator <<(std::ostream& s, const VariablesValue& val)
{
  s<<val.variables()->variablesSetName()<<"\n";
  s<<"name \t value \t  tvalue \t transf \n";
  for (std::size_t i=0; i<val.size(); ++i)
    {
      s<<val.variables()->var(i)->name()<<"\t";
      s<<val.values()[i]<<"\t";
      s<<val.tvalues()[i]<<"\t";
      s<<val.variables()->var(i)->transfromation()->myClass()<<"\n";
    }
  s<<"\n";
  return s;
}

