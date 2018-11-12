#include<iostream>
#include<vector>
#include <string>
#include <algorithm>
#include "SparseMatrix.hh"
#include "FiniteDifference.hh"

int main ()
{
  // SparseMatrix A = SparseMatrix(2,2);
  // A.printMatrix();
  double h = 0.1;
  double L = 1;
  double alpha = 1;
  double beta = 0;
  double gamma = 0;
  // FiniteDifference LAD = FiniteDifference(L,alpha,beta,gamma);
  // LAD.constructMatrix(h);
  runScheme(h, L, alpha, beta, gamma);
  return 0;
}
