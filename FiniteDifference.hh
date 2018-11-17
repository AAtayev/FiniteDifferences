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
  FiniteDifference(); // Default Constructor
  FiniteDifference(int J, double L, double alpha, double beta, double gamma, double b_0, double b_L);
  FiniteDifference(const FiniteDifference& scheme); // Copy Constructor
  ~FiniteDifference(); // Destructor

  int getJ();
  double getL();
  double getAlpha();
  double getBeta();
  double getGamma();
  double getb_0();
  double getb_L();


  SparseMatrix constructMatrix(); // Constructs the "differential operator" matrix A

private:
  int J_; // The number of points in the discretisation
  double L_; // Length of interval
  double b_0_, b_L_; // Boundary conditions
  double alpha_, beta_, gamma_; // Equation parameters
};

// runScheme runs the finite-difference method for given parameter values
std::vector<double> runScheme(int J, double L, double alpha, double beta, double gamma, double b_0, double b_L, std::string fileName);
std::vector<double> analyticalSol(int J, double L, double alpha, double beta, double gamma);
void runMultiple(std::vector<int> Jvec, double L, double alpha, double beta, double gamma, double b_0, double b_L, std::string fileName, std::string errorFileName);
#endif
