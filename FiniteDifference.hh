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
  FiniteDifference(double h, double L, double alpha, double beta, double gamma);
  FiniteDifference(const FiniteDifference& scheme);
  ~FiniteDifference();

  SparseMatrix constructLAD(double alpha, double beta);

  void linearAdvectionDiffusion(double alpha, double beta);
  void linearDiffusion();
  
private:
  double L_;
  double h_;
  double alpha_, beta_, gamma_;
};

#endif
