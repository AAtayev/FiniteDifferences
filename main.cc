#include<iostream>
#include<vector>
#include <string>
#include <algorithm>
#include "SparseMatrix.hh"
#include "FiniteDifference.hh"

int main ()
{
  double h = 0.1;
  double L = 1;
  double alpha = 10;
  double beta = 1;
  double gamma = 0;
  FiniteDifference LAD = FiniteDifference(h,L,alpha,beta,gamma);
  LAD.linearAdvectionDiffusion(alpha,beta);
  return 0;
}
