#include<iostream>
#include<vector>
#include <string>
#include <algorithm>
#include <cmath>
// #include "SparseMatrix.hh"
#include "FiniteDifference.hh"

int main ()
{
  std::vector<int> Jvec = {10, 50, 100, 200, 300, 400, 500, 600, 700};
  std::vector<double> Alpha = {1000,0.25,1};
  std::vector<double> Beta = {1,0.5,21};
  // SparseMatrix A = SparseMatrix(2,2);
  // A.printMatrix();
  // int J = 9;
  double L = 1;
  // double alpha = 1;
  // double beta = 1;
  double gamma = 0;
  double b_0 = 0;
  double b_L = 1;
  // FiniteDifference LAD = FiniteDifference(L,alpha,beta,gamma);
  // LAD.constructMatrix(h);
  // runScheme(J, L, alpha, beta, gamma, b_0, b_L);
  for (int i = 0; i < Alpha.size(); ++i)
  {
    runMultiple(Jvec, L, Alpha[i], Beta[i], gamma, b_0, b_L, "IterResidual", "Error");
  }
  return 0;
}
