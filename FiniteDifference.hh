#ifndef CLASS_FINITEDIFFERENCE
#define CLASS_FINITEDIFFERENCE

#include<iostream>
#include<vector>
#include <string>
#include <algorithm>
#include "SparseMatrix.hh"

class FiniteDifference
{
public:
  FiniteDifference();
  FiniteDifference(double L, double alpha, double beta, double gamma);
  FiniteDifference(const FiniteDifference& scheme);
  ~FiniteDifference();

  double getAlpha();
  double getBeta();
  double getGamma();
  SparseMatrix getA();

  SparseMatrix constructMatrix(double h);

private:
  double L_, alpha_, beta_, gamma_;
  SparseMatrix A_;
};

void runScheme(double h, double L, double alpha, double beta, double gamma);

#endif
