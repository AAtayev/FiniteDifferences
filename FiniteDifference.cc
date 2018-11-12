#include<iostream>
#include<vector>
#include <string>
#include <algorithm>
#include "SparseMatrix.hh"
#include "FiniteDifference.hh"

FiniteDifference::FiniteDifference()
{}

FiniteDifference::FiniteDifference(double L, double alpha, double beta, double gamma)
{
  L_ = L;
  alpha_ = alpha;
  beta_ = beta;
  gamma_ = gamma;
}

FiniteDifference::FiniteDifference(const FiniteDifference& scheme)
{
  L_ = (*this).L_;
  alpha_ = (*this).alpha_;
  beta_ = (*this).beta_;
  gamma_ = (*this).gamma_;
}

FiniteDifference::~FiniteDifference()
{}

double FiniteDifference::getAlpha()
{
  return alpha_;
}

double FiniteDifference::getBeta()
{
  return beta_;
}

double FiniteDifference::getGamma()
{
  return gamma_;
}

SparseMatrix FiniteDifference::getA()
{
  return A_;
}

SparseMatrix FiniteDifference::constructMatrix(double h)
{
  int N = 1/h;
  // std::cout << "Size N = " << N << '\n';
  SparseMatrix A = SparseMatrix(N,N);
  double D = 2*alpha_/(double) (h*h);
  double UD = -(alpha_ + h*beta_)/(double) (h*h);
  double LD = -(alpha_ - h*beta_)/(double) (h*h);
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      if(j == i + 1)
      {
        A.addEntry(i, j, UD);
      }
      else if(j == i)
      {
        A.addEntry(i, j, D);
      }
      else if (j == i - 1)
      {
        A.addEntry(i, j, LD);
      }
    }
  }
  // std::cout << A.getEntry(9,0) << '\n';
  // A.printMatrix();
  A_ = A;
  return A;
}

void runScheme(double h, double L, double alpha, double beta, double gamma)
{
  int N = 1/h;
  std::vector<double> x_0(N);
  std::vector<double> b(N);
  b[N-1] += (alpha - h*beta)/(double) (h*h);
  FiniteDifference Fin = FiniteDifference(L, alpha, beta, gamma);
  // Fin.constructMatrix(h);
  SparseMatrix A = Fin.constructMatrix(h);
  A.printMatrix();
  A.GaussSeidel(x_0,10e-6,1000,b,"Test.dat");
}
