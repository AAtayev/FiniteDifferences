#include<iostream>
#include<vector>
#include <string>
#include <algorithm>
#include <cmath>
// #include "SparseMatrix.hh"
#include "FiniteDifference.hh"

int main ()
{
  std::vector<int> Jvec = {9, 19, 39, 79, 159, 319, 639}; // Set number of discretisation points to run over
  std::vector<double> Alpha = {1,1000,0.25,1}; // Set alpha values to run over
  std::vector<double> Beta = {0,1,0.5,21}; // Set beta values to run over
  double L = 1; // Set interval length
  double gamma = 0; // Set gamma = 0 by default
  double b_0 = 0; // Boundary value at x = 0
  double b_L = 1; // Boundary value at x = L
  for (int i = 0; i < Alpha.size(); ++i) // Run over all alpha and beta values
  {
    runMultiple(Jvec, L, Alpha[i], Beta[i], gamma, b_0, b_L, "IterResidual", "Error");
  }
  return 0;
}
