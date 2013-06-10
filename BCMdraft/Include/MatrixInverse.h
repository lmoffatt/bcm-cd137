#ifndef MATRIXINVERSE_H
#define MATRIXINVERSE_H
#include <vector>
#include <iostream>

std::vector< std::vector< double> >
inv(const std::vector< std::vector< double> >& matrix);

double det(const std::vector< std::vector< double> >& matrix);


double diagProduct(const std::vector< std::vector< double> >& matrix);


std::vector<double> operator-(const std::vector<double>& x,const std::vector<double>& y );


std::vector< std::vector< double> >
chol(const std::vector< std::vector< double> >& matrix);

std::vector< double>
sqrt(const std::vector<double> &avector);

std::vector< double>
diag(const std::vector< std::vector< double> >& matrix);

double
max(const std::vector<double>& avector);

std::vector< std::vector< double> >
prod(const std::vector< std::vector< double> >& A,
     const std::vector< std::vector< double> >& B);

std::vector< std::vector< double> >
trProd(const std::vector< std::vector< double> >& A,
     const std::vector< std::vector< double> >& B);

std::vector< std::vector< double> >
prodTr(const std::vector< std::vector< double> >& A,
     const std::vector< std::vector< double> >& B);



std::vector< std::vector< double> >
sum(const std::vector< std::vector< double> >& A,
     const std::vector< std::vector< double> >& B);


std::vector< std::vector< double> >
mult(const std::vector< std::vector< double> >& A,
    double b);



bool isNaN(const std::vector<double>);


std::ostream& operator<<(
    std::ostream& s,const std::vector< std::vector< double> >& matrix);


std::ostream& operator<<(
    std::ostream& s,const std::vector<  double> & aVector);


double getQuantil(std::vector<double> sortedVector,double p);

void sort(std::vector<double>& v);

double erfinv(double y);

double normal_cdf(double mean,double s, double x);

double normal_cdf(double x);

double normal_pdf(double x);

double normal_pdf( double x,double mean,double std);



double normal_cdf_inv(double mean,double s, double p);

double normal_cdf_inv(double p);



#endif // MATRIXINVERSE_H
