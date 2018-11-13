#include<iostream>
#include <fstream>
#include<vector>
#include <string>
#include <algorithm>
#include <cmath>
#include "SparseMatrix.hh"
#include "FiniteDifference.hh"

FiniteDifference::FiniteDifference()
{}

FiniteDifference::FiniteDifference(int J, double L, double alpha, double beta, double gamma, double b_0, double b_L)
{
  J_ = J;
  L_ = L;
  alpha_ = alpha;
  beta_ = beta;
  gamma_ = gamma;
  b_0_ = b_0;
  b_L_ = b_L;
}

FiniteDifference::FiniteDifference(const FiniteDifference& scheme)
{
  J_ = (*this).J_;
  L_ = (*this).L_;
  alpha_ = (*this).alpha_;
  beta_ = (*this).beta_;
  gamma_ = (*this).gamma_;
  b_0_ = (*this).b_0_;
  b_L_ = (*this).b_L_;
}

FiniteDifference::~FiniteDifference()
{}

int FiniteDifference::getJ()
{
  return J_;
}

double FiniteDifference::getL()
{
  return L_;
}

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

double FiniteDifference::getb_0()
{
  return b_0_;
}

double FiniteDifference::getb_L()
{
  return b_L_;
}
SparseMatrix FiniteDifference::constructMatrix()
{
  double h = L_/(double) (J_ + 1); // Setting the mesh size
  SparseMatrix A = SparseMatrix(J_,J_);
  double D = 2*alpha_/(double) (h*h) + gamma_; // Diagonal terms of A
  double UD = -(alpha_ - h*beta_/2.0)/(double) (h*h); // Upper-diagonal terms of A
  double LD = -(alpha_ + h*beta_/2.0)/(double) (h*h); // Lower-diagonal terms of A
  for (int i = 0; i < J_; ++i)
  {
    for (int j = 0; j < J_; ++j)
    {
      if(j == i + 1) // Upper-diagonal entries
      {
        A.addEntry(i, j, UD);
      }
      else if(j == i) // Diagonal entries
      {
        A.addEntry(i, j, D);
      }
      else if (j == i - 1) // Lower-diagonal entries
      {
        A.addEntry(i, j, LD);
      }
    }
  }
  return A;
}

// runScheme: takes in variables required to perform a finite difference method and equation parameters
// along with boundary conditions and runs the Gauss-Seidel algorithm on the constructed system AU=f,
// where U = x_0 and f = b below.
std::vector<double> runScheme(int J, double L, double alpha, double beta, double gamma, double b_0, double b_L, std::string fileName)
{
  FiniteDifference Fin = FiniteDifference(J, L, alpha, beta, gamma, b_0, b_L); // Initialises the Finite Difference method
  double h = L/(double)(J + 1); // Setting mesh-size
  std::vector<double> x_0(J); // Setting guessed solution at the zero vector, but this can change if needed
  std::vector<double> b(J); // Initialing the vector f = b which is obtained from the finite difference method
  b[0] = Fin.getb_0()*(alpha + h*beta/2.0)/(double) (h*h); // This is equal to f_1
  b[J-1] = Fin.getb_L()*(alpha - h*beta/2.0)/(double) (h*h); // This is equal to f_J
  // Fin.constructMatrix(h);
  SparseMatrix A = Fin.constructMatrix(); // Construct the "differential operator"(or difference method) matrix A
  A.GaussSeidel(x_0,1e-6,10000000,b, fileName); // Perform Gauss-Seidel on the system AU = f, i.e. Ax = b here
  // A.printMatrix();
  // for (int i = 0; i < J; ++i) // This for loop prints the solution obtained from the Gauss-Seidel method
  // {
  //   std::cout.width(20);
  //   std::cout << std::left << x_0[i] << std::endl;
  // }
  // for (int i = 0; i < J; ++i)
  // {
  //   std::cout.width(20);
  //   std::cout << std::left << x_0[i] << std::endl;
  // }
  return x_0;
}

std::vector<double> analyticalSol(int J, double L, double alpha, double beta, double gamma)
{
  std::vector<double> solution(J);
  double h = L/(double)(J + 1);
  if (alpha != 0 && gamma == 0)
  {
    for(int i = 0; i < J; ++i)
    {
      solution[i] = (1 - exp((beta/alpha)*(i+1)*h))/(1 - exp((beta/alpha)*L));
    }
  }
  else if (alpha != 0 && beta == 0)
  {
    for(int i = 0; i < J; ++i)
    {
      solution[i] = (sinh((gamma/alpha)*(i+1)*h))/(sinh((gamma/alpha)*L));
    }
  }
  return solution;
}

void runMultiple(std::vector<int> Jvec, double L, double alpha, double beta, double gamma, double b_0, double b_L, std::string fileName, std::string errorFileName)
{
  double Pe = fabs(beta)*L/(double) (2*alpha);
  std::ofstream myFile;
  myFile.open(errorFileName + "_Pe_" + std::to_string(Pe) + ".dat", std::ios::out); // Creates and opens a file in the name of myName
  if (!myFile.good())
  {
    throw std::invalid_argument("Failed to open file");
  }
  myFile.width(20);
  myFile << std::left << "Number of points J" << "Error" << std::endl;
  for (int j: Jvec)
  {
    std::vector<double> x_sol = runScheme(j, L, alpha, beta, gamma, b_0, b_L, fileName + std::to_string(j) + "_Pe_" + std::to_string(Pe) + ".dat");
    std::vector<double> x_an = analyticalSol(j, L, alpha, beta, gamma);
    myFile.width(20);
    myFile << std::left << j << inftyNorm(v_minus_w(x_sol,x_an)) << std::endl;
    // for (int i = 0; i < j; ++i)
    // {
    //   std::cout.width(20);
    //   std::cout << std::left << x_sol[i] << x_an[i] <<std::endl;
    // }
  }
  myFile.close();
}
