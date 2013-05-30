#ifndef RUNGEKUTTA4_H
#define RUNGEKUTTA4_H
#include <vector>


class ABC_ODE
{
public:
    virtual std::vector<double> Derivative(double t,std::vector<double> y)=0;

};


class RungeKutta4
{

public:
    RungeKutta4(ABC_ODE* ode,std::vector<double> initialY,double initialTime=0);

    const std::vector<double>& next(double h);

private:
    ABC_ODE* ode_;
    std::vector<double> y_;
    double t_;

};

std::vector<double>& operator+=(std::vector<double>& x, const std::vector<double>& y);


std::vector<double> operator+(const std::vector<double>& x, const std::vector<double>& y);


std::vector<double> operator*(const std::vector<double>& x, double h);


#endif // RUNGEKUTTA4_H
