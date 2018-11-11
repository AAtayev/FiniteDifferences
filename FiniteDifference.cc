#include<iostream>
#include<vector>
#include <string>
#include <algorithm>
#include "SparseMatrix.hh"
#include "FiniteDifference.hh"

FiniteDifference::FiniteDifference()
{}

FiniteDifference::FiniteDifference(double h, double L, double alpha, double beta, double gamma)
{
  h_ = h;
  L_ = L;
  alpha_ = alpha;
  beta_ = beta;
  gamma_ = gamma;
}

FiniteDifference::FiniteDifference(const FiniteDifference& scheme)
{
  h_ = (*this).h_;
  L_ = (*this).L_;
  alpha_ = (*this).alpha_;
  beta_ = (*this).beta_;
  gamma_ = (*this).gamma_;
}

FiniteDifference::~FiniteDifference()
{}

SparseMatrix FiniteDifference::constructLAD(double alpha, double beta)
{
  int N = 1/h_;
  double D = 2*alpha_/(double) h_*h_;
  double UD = -(alpha_ + h_*beta_)/(double) (h_*h_);
  double LD = -(alpha_ - h_*beta_)/(double) (h_*h_);
  SparseMatrix A = SparseMatrix(N,N);
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      if(j == i - 1)
      {
        A.addEntry(i,j, UD);
      }
      else if(j == i)
      {
        A.addEntry(i,j,D);
      }
      else if (j == i - 1)
      {
        A.addEntry(i,j,LD);
      }
    }
  }
  return A;
}

void FiniteDifference::linearAdvectionDiffusion(double alpha, double beta)
{
  int N = 1/h_;
  std::vector<double> x_0(N);
  std::vector<double> b(N,0);
  b[N-1] = (alpha_ - h_*beta_)/(double) (h_*h_);
  SparseMatrix A = constructLAD(alpha, beta);
  A.GaussSeidel(x_0,10e-6,1000,b,"Test.dat");
}
