#include <cmath>
#include <tuple>

#include "Include/Experiment.h"
#include "Include/ABC_BCM.h"
#include "Include/MatrixInverse.h"

const VariablesValue& ExperimentSchedule::ExternalState()const
{
  return externalVar_;
}

const MultipleVariables& ExperimentSchedule::measuredVariables()const
{
  return measuredVars_;
}
const std::vector<double>& ExperimentSchedule::measuredTimes()const
{
  return measuredTimes_;
}
double ExperimentSchedule::tmax()const
{
  return tmax_;
}

ExperimentSchedule::ExperimentSchedule(const ABC_BCM* model,
                                       const std::string& experimentIndentifier,
                                       const std::vector<VarValue> &externalState,
                                       const std::vector<VariableTimeDistribution> &data):
  name_(experimentIndentifier),
  model_(model),
  externalVar_(getExternalState(externalState)),
  measuredVars_(getMeasuredValues(experimentIndentifier,data)),
  measuredTimes_(getMeasuredTimes(data)),
  tmax_(max(measuredTimes_)){}

VariablesValue ExperimentSchedule::getExternalState(
    std::vector<VarValue> externalSate)
{
  const MultipleVariables* mv=model_->getExternalStateVariables();
  VariablesValue extState(mv);
  std::vector<double> val(externalSate.size());
  for (std::size_t i=0; i<externalSate.size(); i++)
    {
      VarValue *va=&externalSate[i];
      val[va->index]=va->value;
    }
  extState.setValues(val);
  return extState;
}

MultipleVariables ExperimentSchedule::getMeasuredValues(const std::string& experimentIdentifier,
                                                        std::vector<VariableTimeDistribution> data)
{
  const MultipleVariables* obserVar=model_->getObservedStateVariables();
  std::multimap<int,const SingleVariable*> vars;
  std::vector<const SingleVariable*> var;
  for (std::size_t i=0; i<data.size(); i++)
    {

      VariableTimeDistribution *vt=&data[i];
      vars.insert({vt->index,obserVar->var(vt->index)});

    }
  for (auto it=vars.begin(); it!=vars.end(); ++it)
    {
      var.push_back(it->second);
    }
  MultipleVariables res(model_,experimentIdentifier+".measures",var);
  return res;

}

std::string ExperimentSchedule::nameId()const
{
  return name_;
}

std::vector<double> ExperimentSchedule::getMeasuredTimes(
    std::vector<VariableTimeDistribution> data)
{
  std::map<int,std::set<double> > times;
  std::vector<double> time;
  for (std::size_t i=0; i<data.size(); i++)
    {
      VariableTimeDistribution* va=&data[i];
      times[va->index].insert(va->time);
    }
  for (auto it=times.begin(); it!=times.end(); ++it)
    {
      for (auto it2=it->second.begin(); it2!=it->second.end(); ++it2)
        time.push_back(*it2);
    }
  return time;
}



const ExperimentSchedule& ExperimentDistribution::schedule()const
{
  return sche_;
}
const VariablesNormalDistribution& ExperimentDistribution::distr()const
{
  return distr_;
}

ExperimentDistribution::ExperimentDistribution(const ABC_BCM* model,
                                               const std::string& experimentIndentifier,
                                               const std::vector<VarValue> &externalState,
                                               const std::vector<VariableTimeDistribution>& data):
  sche_(ExperimentSchedule(model,experimentIndentifier,externalState,data)),
  distr_(getYvalues(model,sche_,data)){}



ExperimentDistribution& ExperimentDistribution::setTmeanValues(const std::vector<double> tmean)
{
  distr_.setTmeanValues(tmean);
  return *this;
}

ExperimentDistribution& ExperimentDistribution::setMeanValues(const std::vector<double> mean)
{
  distr_.setMeanValues(mean);
  return *this;
}




const std::vector<double>& ExperimentDistribution::w()const{
  return distr().w();
}
const std::vector<double>& ExperimentDistribution::y()const
{
  return this->distr().center().values();
}
const std::vector<double> &ExperimentDistribution::Ty()const{
  return this->distr().center().tvalues();

}
const std::vector<double> &ExperimentDistribution::yerr()const
{
  return distr().tstd();
}



VariablesNormalDistribution
ExperimentDistribution::getYvalues(const ABC_BCM *model,
                                   const ExperimentSchedule& sche,
                                   std::vector<VariableTimeDistribution> data)
{
  std::map<int,std::map<double,std::pair<double,double> > > timeyyerr;
  std::vector<double> yv;
  std::vector<double> yerrv;

  for (std::size_t i=0; i<data.size(); i++)
    {
      VariableTimeDistribution* vt=&data[i];
      const RangeFormat* f=model->getRangeFormat(vt->FormatIndex);
      const SingleVariable* v=model->getObservedStateVariables()->var(vt->index);
      const Transformation* T=v->transfromation();
      auto p=f->getTransformedMoments(T,vt->first,vt->second);
          timeyyerr[vt->index][vt->time]=p;
    }
  for (auto it=timeyyerr.begin(); it!=timeyyerr.end(); ++it)
    {
      for (auto it2=it->second.begin(); it2!=it->second.end(); ++it2)
        {
          yv.push_back(it2->second.first);
          yerrv.push_back(it2->second.second);

        }
    }
  VariablesNormalDistribution dist(&sche.measuredVariables(),yv,yerrv);
  return dist;

}


std::ostream& operator<<(std::ostream& s, const std::vector<const ExperimentDistribution*>& exp)
{
  for (std::size_t i=0;i<exp.size(); ++i)
    s<<*exp[i]<<"\n";
  return s;
}

std::ostream& operator<<(std::ostream& s, const ExperimentDistribution& exp)
{
  s<<exp.schedule().nameId()<<std::endl;
  s<<exp.schedule().ExternalState();
  s<<"Observed variables\n";
  s<<"name \t time \t mean \t m-se  \t  m+se \t transf\n";
  for (std::size_t i=0; i<exp.distr().size(); ++i)
    {
      s<<exp.distr().variables()->var(i)->name()<<"\t";
      s<<exp.schedule().measuredTimes()[i]<<"\t";
      s<<exp.y()[i]<<"\t";
      s<<exp.distr().variables()->var(i)->transfromation()->inverse(exp.Ty()[i]-exp.yerr()[i])<<"\t";
      s<<exp.distr().variables()->var(i)->transfromation()->inverse(exp.Ty()[i]+exp.yerr()[i])<<"\t";
      s<<exp.distr().variables()->var(i)->transfromation()->myClass()<<"\n";

    }

  return s;
}





/*

const VariablesValue &Experiment::ExternalState()const
{
  return this->eState;
}

const std::map<double, std::vector<SingleVariableNormalDistribution> > &Experiment::timeMeasures()
const{
  return obs_;
}

void Experiment::setyfit(const std::vector<double> &newy)
{
  y_=newy;
}


Experiment::Experiment(const std::string& experimentIndentifier,
                       std::vector<std::pair<int,double>> externalState,
                       std::vector<std::tuple<int,double,double,double>> data):
  id(experimentIndentifier),
  y_(),yerr_(),w_()
{
  for (std::size_t i=0; i<externalState.size(); i++)
    {
      int j=externalState[i].first;
      if (j<eState.size())
        {
          eState[j]=externalState[i].second;
        }
      else
        {
          eState.resize(j+1,0);
          eState[j]=externalState[i].second;
        }
    }
  for (std::size_t i=0; i<data.size(); i++)
    {

      Measure m(std::get<0>(data[i]),std::get<2>(data[i]),std::get<3>(data[i]));
      obs_[std::get<1>(data[i])].push_back(m);
       di_[std::get<0>(data[i])].insert(std::pair<double,Measure>(std::get<1>(data[i]),m));
    }

  for (auto& it :di_)
    {
      measVar_.insert(it.first);
      for (auto& t: it.second)
        {
          y_.push_back(t.second.value());
          yerr_.push_back(t.second.se());
          w_.push_back(1.0/std::pow(t.second.se(),2));

        }
    }

  tmax_=obs_.rbegin()->first;
}

std::size_t Experiment::size()const
{
  return y_.size();
}
const std::vector<double>& Experiment::y()const
{
  return y_;
}
const std::vector<double>& Experiment::yerr()const
{
  return yerr_;
}
const std::vector<double> &Experiment::w()const
{
  return w_;
}

double Experiment::tmax()const
{
  return tmax_;
}

const std::map<double, SingleVariableNormalDistribution> &Experiment::timeMeasure(int index)const
{
  return di_.at(index);
}


const std::set<int>& Experiment::getMeasVarSet()const
{
  return this->measVar_;
}

*/

