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

  int getJ(); // retrives private variable J_
  double getL(); // retrives private variable L_
  double getAlpha(); // retrives private variable alpha_
  double getBeta(); // retrives private variable beta_
  double getGamma(); // retrives private variable gamma_
  double getb_0(); // retrives private variable b_0_
  double getb_L(); // retrives private variable b_L_


  SparseMatrix constructMatrix(); // Constructs the "differential operator" matrix A

private:
  int J_; // The number of points in the discretisation
  double L_; // Length of interval
  double b_0_, b_L_; // Boundary conditions
  double alpha_, beta_, gamma_; // Equation parameters
};

// runScheme: runs the finite-difference method for given parameter values
std::vector<double> runScheme(int J, double L, double alpha, double beta, double gamma, double b_0, double b_L, std::string fileName);

// analyticalSol: for given values J, L, alpha, beta and gamma, constructs the analytical solution on a vector of size J (i.e. r_h u in the notes).
//                We assume that b_0 = 0 and b_L = 1 here for simplicity, but this can also be changed if required.
std::vector<double> analyticalSol(int J, double L, double alpha, double beta, double gamma);

// runMultiple: for the given inputs, runMultiple performs runScheme for all the given discretisation points J in Jvec,
//              and further prints the error between the analytical solution and the numerical solution on a file.
void runMultiple(std::vector<int> Jvec, double L, double alpha, double beta, double gamma, double b_0, double b_L, std::string fileName, std::string errorFileName);
#endif
