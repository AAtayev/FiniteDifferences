#include <cassert>
#include <cmath>

#include "vector.hh"

std::ostream& operator<<(std::ostream& s, const Vector& v) {
  for (double l : v ) {
	s << l << std::endl;
  }
  return s;
}

double Vector::norm() const {
  double ret = 0.0;
  for ( int i=0; i<size(); ++i ) {
	ret += (*this)[i] * (*this)[i];
  }
  return std::sqrt(ret);
}

double Vector::maxNorm() const {
  double ret = 0.0;
  for ( int i=0; i<size(); ++i ) {
	const double v = std::abs((*this)[i]);
	if (v > ret) ret = v;
  }
  return ret;
}

Vector Vector::operator+(const Vector& other) const {
  Vector ret(size());
  for ( int i=0; i<size(); ++i ) {
	ret[i] = (*this)[i] + other[i];
  }
  return ret;
}

Vector Vector::operator-(const Vector& other) const {
  Vector ret(size());
  for ( int i=0; i<size(); ++i ) {
	ret[i] = (*this)[i] - other[i];
  }
  return ret;
}

double Vector::operator*(const Vector& other) const {
  double ret = 0.0;
  for ( int i=0; i<size(); ++i ) {
	ret += (*this)[i]*other[i];
  }
  return ret;
}

Vector Vector::operator*(double v) const {
  Vector ret(*this);
  for ( int i=0; i<size(); ++i ) {
	ret[i] = (*this)[i]*v;
  }
  return ret;
}
