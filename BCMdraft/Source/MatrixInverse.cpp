#include "Include/MatrixInverse.h"
#include <cmath>
#include <algorithm>

extern "C" void dgetrf_(int *M,
                        int* N,
                        double *A,
                        int* LDA,
                        int* IPIV,
                        int * INFO );

extern "C" void dgetri_(int* n,
                        double *B,
                        int* dla,
                        int* ipiv,
                        double* work1,
                        int* lwork,
                        int* info);


std::vector< double>
sqrt(const std::vector< double>& avector)
{
  std::vector<double> res(avector.size());
  for (unsigned i=0; i<avector.size(); ++i)
    res[i]=std::sqrt(avector[i]);
  return res;
}

std::vector< double>
diag(const std::vector< std::vector< double> >& matrix){
  std::vector<double> res(std::min(matrix.size(),matrix[0].size()));
for (unsigned i=0; i< res.size(); ++i)
  res[i]=matrix[i][i];
return res;
}


double
max(const std::vector<double>& avector)
{
  double themax=avector[0];
  for (std::size_t i=1; i<avector.size(); ++i)
    if (avector[i]>themax)
      themax=avector[i];
  return themax;
}






std::vector< std::vector< double> >
inv(const std::vector< std::vector< double> >& matrix)
{
  double *A, *work, *work1;
  int info=0;
  //  char msg[101];
  int *ipiv;
  int lwork;
  int n =matrix.size();
  int m=n;
  int dla=n;
  //A=new double[n*n];
  A= new double[n*n]; //more efficient code
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      *(A+i+n*j) = matrix[i][j];

  ipiv = new int[n];
  dgetrf_(&n, &m, A, &dla,ipiv,&info);
  lwork= -1;
  work1 = new double[2];
  dgetri_(&n,A,&dla,ipiv,work1,&lwork,&info);
  lwork = (int)(work1[0]);
  work = new double [2*lwork];
  dgetri_(&n,A,&dla,ipiv,work,&lwork,&info);

  std::vector< std::vector<double> > result(n,std::vector<double>(n));
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      result[i][j]=*(A+i+n*j);
  delete [] ipiv;
  delete [] work;
  delete [] work1;
  delete [] A;

  if (info>0)
    std::cerr<<"\n singular matrix\n";

  return result;
}


bool isNaN(const std::vector<double> x)
{
  for (std::size_t i=0; i<x.size();i++)
    if (x[i]!=x[i])
      return true;
  return false;
}


std::vector< std::vector< double> >
prod(const std::vector< std::vector< double> >& A,
     const std::vector< std::vector< double> >& B)
{
  std::size_t nrow=A.size();
  std::size_t ncol=B[0].size();
  std::vector< std::vector<double> > result(nrow,std::vector<double>(ncol));

  for (std::size_t i=0; i<nrow;i++)
    for(std::size_t j=0; j<ncol; j++)
      {
        result[i][j]=0;
        for (std::size_t k=0; k<B.size(); k++)
          result[i][j]+=A[i][k]*B[k][j];
      }
  return result;

}

std::vector< std::vector< double> >
trProd(const std::vector< std::vector< double> >& A,
       const std::vector< std::vector< double> >& B)
{
  std::size_t nrow=A[0].size();
  std::size_t ncol=B[0].size();
  std::vector< std::vector<double> > result(nrow,std::vector<double>(ncol));

  for (std::size_t i=0; i<nrow;i++)
    for(std::size_t j=0; j<ncol; j++)
      {
        result[i][j]=0;
        for (std::size_t k=0; k<B.size(); k++)
          result[i][j]+=A[k][i]*B[k][j];
      }
  return result;
}

std::vector< std::vector< double> >
prodTr(const std::vector< std::vector< double> >& A,
       const std::vector< std::vector< double> >& B)
{
  std::size_t nrow=A.size();
  std::size_t ncol=B.size();
  std::vector< std::vector<double> > result(nrow,std::vector<double>(ncol));

  for (std::size_t i=0; i<nrow;i++)
    for(std::size_t j=0; j<ncol; j++)
      {
        result[i][j]=0;
        for (std::size_t k=0; k<B[0].size(); k++)
          result[i][j]+=A[i][k]*B[j][k];
      }
  return result;
}




std::vector< std::vector< double> >
sum(const std::vector< std::vector< double> >& A,
    const std::vector< std::vector< double> >& B)
{
  std::vector< std::vector<double> > result(A);
  std::size_t nrow=A.size();
  std::size_t ncol=A[0].size();

  for (std::size_t i=0; i<nrow;i++)
    for(std::size_t j=0; j<ncol; j++)
      {
        result[i][j]+=B[i][j];
      }
  return result;
}


double diagProduct(const std::vector< std::vector< double> >& matrix)
{
  double p=1;

  for (std::size_t i=0; i<(std::min(matrix.size(),matrix[0].size())); ++i)
    p*=matrix[i][i];
  return p;
}



std::vector< std::vector< double> >
mult(const std::vector< std::vector< double> >& A,
     double x)
{
  std::vector< std::vector<double> > result(A);
  std::size_t nrow=A.size();
  std::size_t ncol=A[0].size();

  for (std::size_t i=0; i<nrow;i++)
    for(std::size_t j=0; j<ncol; j++)
      {
        result[i][j]*=x;
      }
  return result;
}




std::ostream& operator<<(
    std::ostream& s,const std::vector< std::vector< double> >& matrix)
{
  s<<"\n";
  for (std::size_t i=0; i<matrix.size();++i)
    {
      for (std::size_t j=0; j<matrix[i].size();j++)
        s<<matrix[i][j]<<"\t";
      s<<"\n";
    }
  return s;
}
std::ostream& operator<<(
    std::ostream& s,const std::vector<  double> & aVector){
  for (std::size_t j=0; j<aVector.size();j++)
    s<<aVector[j]<<"\t";
  s<<"\n";
  return s;
}

double getQuantil(std::vector<double> sortedVector,double p)
{
  double r=sortedVector.size()*p;
  std::size_t n=floor(r);

  if ((r>n)&&(n+1<sortedVector.size()))
    return sortedVector[n]*(n+1-r)+sortedVector[n+1]*(r-n);
  else
    return sortedVector[n];
}

void sort(std::vector<double>& v)
{
  std::sort(v.begin(),v.end());
}


double erfinv(double y)
{

  //rational approx coefficients
  const  double a[]= { 0.886226899, -1.645349621,  0.914624893, -0.140543331};
  const  double  b[] = {-2.118377725,  1.442710462, -0.329097515,  0.012229801};
  const  double c[]= {-1.970840454, -1.624906493,  3.429567803,  1.641345311};
  const  double d[]={ 3.543889200,  1.637067800};

  const double y0 = 0.7;

  double x;

  if ((y<-1.0 ||y> 1.0))
    return std::numeric_limits<double>::quiet_NaN();

  if (std::abs(y)==1.0)
    x= -y*std::log(0.0);
  else if (y<-y0)
    {
      double z = std::sqrt(-std::log((1.0 + y)/2.0));
      x=  -(((c[3]*z+c[2])*z+c[1])*z+c[0])/((d[1]*z+d[0])*z+1.0);
    }
  else
    {
      if (y<y0)
        {
          double z= y*y;
          x = y*(((a[3]*z+a[2])*z+a[1])*z+a[0])/((((b[3]*z+b[3])*z+b[1])*z+b[0])*z+1.0);
        }
      else
        {
          double z = sqrt(-std::log((1.0-y)/2.0));
          x = (((c[3]*z+c[2])*z+c[1])*z+c[0])/((d[1]*z+d[0])*z+1.0);
        }
      //polish x to full accuracy
      x = x - (std::erf(x) - y) / (2.0/sqrt(M_PI) * std::exp(-x*x));
      x = x - (std::erf(x) - y) / (2.0/sqrt(M_PI) * std::exp(-x*x));
    }
  return x;
}

std::vector<double> operator-(const std::vector<double>& x,const std::vector<double>& y )
{
  std::vector<double> res(x.size());
  for (std::size_t i=0; i<x.size(); i++)
    res[i]=x[i]-y[i];
  return res;

}



double normal_cdf_inv(double p)
{
  const double r2=std::sqrt(2);
  return r2*erfinv(2.0*p-1.0);
}

double normal_cdf(double x)
{
  const double r2=std::sqrt(2);
  return 0.5*erf(1.0+x/r2);
}


double normal_cdf(double mean, double s, double x)
{
  return normal_cdf((x-mean)/s);
}


double normal_cdf_inv(double mean,double s, double p)
{
  return normal_cdf_inv(p)*s+mean;
}


double normal_pdf(double x)
{
  const double sqrt2Pi=1.0/std::sqrt(2.0*M_PI);
  return  sqrt2Pi*std::exp(-x*x/2.0);
}



double normal_pdf( double x,double mean,double std)
{
  return normal_pdf((x-mean)/std);
}
