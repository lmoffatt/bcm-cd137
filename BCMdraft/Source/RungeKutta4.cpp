#include "Include/RungeKutta4.h"


RungeKutta4::RungeKutta4(ABC_ODE* ode,std::vector<double> initialY,double initialTime):
    ode_(ode),
   y_(initialY),
    t_(initialTime)

{}


std::vector<double>& operator+=(std::vector<double>& x, const std::vector<double>& y)
{
        for (std::size_t i=0; i<x.size(); i++)
            x[i]+=y[i];
        return x;
}



std::vector<double> operator+(const std::vector<double>& x, const std::vector<double>& y)
{
    std::vector<double> z(x);
    z+=y;
    return z;

}

std::vector<double> operator*(const std::vector<double>& x, double h)
{
    std::vector<double> z(x.size());
    for (std::size_t i=0; i<x.size(); i++)
        z[i]=x[i]*h;
    return z;

}



const std::vector<double>& RungeKutta4::next(double h)
{

    std::vector<double> k1=ode_->Derivative(t_,y_)*h;

    std::vector<double> k2=ode_->Derivative(t_+h*0.5,y_+k1*0.5)*h;

    std::vector<double> k3=ode_->Derivative(t_+h*0.5,y_+k2*0.5)*h;

    std::vector<double> k4=ode_->Derivative(t_+h,y_+k3)*h;

    std::vector<double> delta=(k1+k2*2.0+k3*2.0+k4)*(1.0/6.0);

    y_+=delta;
    t_+=h;

    return y_;

}


