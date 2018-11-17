#include<iostream>
#include<vector>
#include <string>
#include <algorithm>
#include <cmath>
// #include "SparseMatrix.hh"
#include "FiniteDifference.hh"

int main ()
{
  std::vector<int> Jvec = {9, 19, 39, 79, 159, 319, 639};
  // std::vector<int> Jvec = {9};
  std::vector<double> Alpha = {1,1000,0.25,1};
  std::vector<double> Beta = {0,1,0.5,21};
  double L = 1;
  double gamma = 0;
  double b_0 = 0;
  double b_L = 1;
  for (int i = 0; i < Alpha.size(); ++i)
  {
    runMultiple(Jvec, L, Alpha[i], Beta[i], gamma, b_0, b_L, "IterResidual", "Error");
  }
  return 0;
}
