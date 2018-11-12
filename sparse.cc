#include <cassert>

#include "sparse.hh"

int SparseMatrix::findEntry(int i, int j) const {
  for ( int k = 0; k<index_[i].size(); ++k ) {
	if ( index_[i][k] == j ) {
	  return k;
	}
  }
  return -1;
}

double SparseMatrix::getEntry(int i, int j) const {
  int idx = findEntry(i,j);
  if (idx  == -1 ) {
	return 0.0;
  } else {
	return data_[i][idx];
  }
}

void SparseMatrix::addEntry(int i, int j, double v) {
  int idx = findEntry(i,j);
  if (idx  == -1 ) {
	index_[i].push_back(j);
	data_[i].push_back(v);
  } else {
	data_[i][idx] = v;
  }
}

Vector SparseMatrix::multiply(const Vector& other) const {
  Vector x(N);
  for ( int i = 0; i<N; ++i ) {
	x[i] = multiplyLine(other, i);
  }
  return x;
}

double SparseMatrix::multiplyLine(const Vector& x, int i) const
{
  double ret = 0.0;
  for ( int k = 0; k<index_[i].size(); ++k ) {
	ret += data_[i][k]*x[index_[i][k]];
  }
  return ret;
}

Vector SparseMatrix::GaussSeidel(const Vector& b, const Vector& x0, double tolerance, int maxIter) const {
  assert(b.size() == N);
  const SparseMatrix& A = (*this);
  Vector x = x0;
  for ( int iter=0; iter < maxIter; ++iter ) {
	for ( int i = 0; i<N; ++i ) {
	  const double a_ii = getEntry(i,i);
	  x[i] += (b[i] - multiplyLine(x, i))/a_ii;
	}
	if ( (iter % 100) == 0 ) {
	  const double resi = (b - A*x).maxNorm();
	  std::cout << "Iteration " << iter << ": residual=" << resi << std::endl;
	  if ( resi < tolerance ) {
		return x;
	  }
	}
  }
  std::cerr << "Tolerance not reached after " << maxIter << " iterations. Aborting!" << std::endl;
  return x;
}

Vector SparseMatrix::SteepestDescent(const Vector& b, const Vector& x0, double tolerance, int maxIter) const {
  assert(b.size() == N);
  const SparseMatrix& A = (*this);
  Vector x = x0;
  for ( int iter=0; iter < maxIter; ++iter ) {
	Vector r = b - A*x;
	if ( (iter % 100) == 0 ) {
	  const double resi = r.maxNorm();
	  std::cout << "Iteration " << iter << ": residual=" << resi << std::endl;
	  if ( resi < tolerance ) {
		return x;
	  }
	}
	double alpha = (r*r)/(A*r*r);
	x = x + r*alpha;
  }
  std::cerr << "Tolerance not reached after " << maxIter << " iterations. Aborting!" << std::endl;
  return x;
}

Vector SparseMatrix::ConjugateGradient(const Vector& b, const Vector& x0, double tolerance, int maxIter) const {
  assert(b.size() == N);
  const SparseMatrix& A = (*this);
  Vector x = x0;
  Vector r = b - A*x;
  Vector p = r;

  for ( int iter=0; iter < maxIter; ++iter ) {
	const double App = A*p*p;
	const double alpha = (p*r)/App;
	x = x + p*alpha;
	r = b - A*x;
	const double beta = -(A*p*r)/App;
	p = r + p*beta;
	if ( (iter % 10) == 0 ) {
	  const double resi = r.maxNorm();
	  std::cout << "Iteration " << iter << ": residual=" << resi << std::endl;
	  if ( resi < tolerance ) {
		return x;
	  }
	}
  }
  std::cerr << "Tolerance not reached after " << maxIter << " iterations. Aborting!" << std::endl;
  return x;
}

void SparseMatrix::printMatrix()
{
  for (unsigned int i = 0; i < N; ++i)
  {
    for (unsigned int j = 0; j < N; ++j)
    {
      std::cout.width(5);
      std::cout << getEntry(i,j);
    }
    std::cout << std::endl;
  }
}
