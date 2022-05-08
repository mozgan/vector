#include "vector.h"

#include <limits.h>
#include <math.h>
#include <string.h>

#include <iomanip>
#include <iostream>
#include <numeric>

#include "random.h"

using namespace std;

namespace algebra {
namespace linear {

/* *** *** *** *** *** *** *** Constructors *** *** *** *** *** *** *** */
/* *** *** *** default constructor *** *** *** */
Vector::Vector() noexcept : size_(0), values_(nullptr) {}

/* *** *** *** constructor with a given size *** *** *** */
Vector::Vector(const size_t size) noexcept
    : size_(size), values_(new double[size_]) {
  memset(values_, 0, sizeof(double) * size_);
}

/* *** *** *** copy constructor *** *** *** */
Vector::Vector(const Vector& src) noexcept : Vector(src.size()) {
  memcpy(values_, src.values_, sizeof(double) * size());
}

/* *** *** *** move constructor *** *** *** */
Vector::Vector(Vector&& src) noexcept : Vector(0) { swap(*this, src); }

/* *** *** *** *** *** *** *** *** Operators *** *** *** *** *** *** *** *** */
/* *** *** *** copy-assignment operator *** *** *** */
Vector& Vector::operator=(const Vector& rhs) noexcept {
  if (this != &rhs) {
    Vector temp(rhs.size());
    swap(*this, temp);
    memcpy(values_, rhs.values_, sizeof(double) * size());
  }

  return *this;
}

/* *** *** *** move-assignment operator *** *** *** */
Vector& Vector::operator=(Vector&& rhs) noexcept {
  if (this != &rhs) {
    Vector empty(0);
    swap(*this, rhs);
    swap(rhs, empty);
  }

  return *this;
}

/* *** *** *** operator[]: setter *** *** *** */
double& Vector::operator[](const size_t i) {
  if (i >= size()) throw std::out_of_range("Vector::operator[]: out of range");

  return values_[i];
}

/* *** *** *** operator[]: getter *** *** *** */
const double Vector::operator[](const size_t i) const {
  if (i >= size()) throw std::out_of_range("Vector::operator[]: out of range");

  return values_[i];
}

/* *** *** *** operator-(): negation *** *** *** */
Vector Vector::operator-() const noexcept {
  Vector result{*this};
  for (size_t i = 0; i < size_; ++i) result[i] = -result[i];
  return result;
}

/* *** *** *** operator+=: vector addition and assignment *** *** *** */
Vector& Vector::operator+=(const Vector& rhs) {
  if (size() != rhs.size())
    throw logic_error("Vector::operator+=: Missmatch array dimensions");

  for (size_t i = 0; i < size_; ++i) values_[i] += rhs[i];
  return *this;
}

/* *** *** *** operator-=: vector subtraction and assignment *** *** *** */
Vector& Vector::operator-=(const Vector& rhs) {
  if (size() != rhs.size())
    throw logic_error("Vector::operator-=: Missmatch array dimensions");

  for (size_t i = 0; i < size_; ++i) values_[i] -= rhs[i];
  return *this;
}

/* *** *** operator*=: scalar-vector multiplication and assignment *** *** */
Vector& Vector::operator*=(const double scalar) noexcept {
  for (size_t i = 0; i < size(); ++i) values_[i] *= scalar;
  return *this;
}

/* *** operator/=: division of a vector by a scalar and assignment *** */
Vector& Vector::operator/=(const double scalar) {
  if (scalar == 0)
    throw logic_error("Vector::operator/=: Divison by zero not allowed.");

  for (size_t i = 0; i < size(); ++i) values_[i] /= scalar;
  return *this;
}

/* *** *** *** *** *** *** *** Destructor *** *** *** *** *** *** *** */
Vector::~Vector() { delete[] values_; }

/* *** *** *** *** *** *** *** Friend Functions *** *** *** *** *** *** *** */
/* *** *** *** swap() *** *** *** */
void swap(Vector& lhs, Vector& rhs) noexcept {
  ::swap(lhs.size_, rhs.size_);
  ::swap(lhs.values_, rhs.values_);
}

/* *** *** *** operator<< *** *** *** */
ostream& operator<<(ostream& os, const Vector& src) {
  for (size_t i = 0; i < src.size(); ++i) os << setw(5) << src[i] << endl;
  return os;
}

/* *** *** *** *** *** *** *** Member Functions *** *** *** *** *** *** *** */
/* *** *** *** size() *** *** *** */
const size_t Vector::size() const noexcept { return size_; }

/* *** *** *** *** *** *** Non-member Functions *** *** *** *** *** *** */
/* *** *** *** operator+ *** *** *** */
Vector operator+(const Vector& v, const Vector& u) {
  if (v.size() != u.size())
    throw logic_error("operator+: Missmatch array dimensions.");

  Vector result(v);
  return (result += u);
}

/* *** *** *** operator- *** *** *** */
Vector operator-(const Vector& v, const Vector& u) {
  if (v.size() != u.size())
    throw logic_error("operator-: Missmatch array dimensions.");

  Vector result(v);
  return (result -= u);
}

/* *** *** *** operator* *** *** *** */
Vector operator*(const double scalar, const Vector& v) {
  Vector u{v};
  return (u *= scalar);
}

/* *** *** *** operator/ *** *** *** */
Vector operator/(const Vector& v, const double scalar) {
  if (scalar == 0) throw logic_error("operator/: Divison by zero not allowed.");

  Vector u{v};
  return (u /= scalar);
}

/* *** *** *** operator== *** *** *** */
const bool operator==(const Vector& lhs, const Vector& rhs) {
  if (lhs.size() != rhs.size()) return false;

  for (size_t i = 0; i < lhs.size(); ++i)
    if (lhs[i] != rhs[i]) return false;

  return true;
}

/* *** *** *** operator!= *** *** *** */
const bool operator!=(const Vector& lhs, const Vector& rhs) {
  return !(lhs == rhs);
}

/* *** *** *** dot() (dot product (aka. inner product)) *** *** *** */
const double dot(const Vector& v, const Vector& u) {
  if (v.size() != u.size())
    throw logic_error("dot(): Missmatch array dimensions.");

  double acc{0};
  for (size_t i = 0; i < v.size(); ++i) acc += (v[i] * u[i]);
  return acc;
}

/* *** *** *** l1_norm() (Manhattan distance) *** *** *** */
const double l1_norm(const Vector& v) {
  double acc{0};
  for (size_t i = 0; i < v.size(); ++i) acc += abs(v[i]);
  return acc;
}

/* *** *** *** l2_norm() (Euclidean norm or Euclidean distance) *** *** *** */
const double l2_norm(const Vector& v) {
  double acc{0};
  for (size_t i = 0; i < v.size(); ++i) acc += (v[i] * v[i]);
  return sqrtl(acc);
}

/* *** *** *** p_norm() *** *** *** */
const double p_norm(const Vector& v) {
  double acc{0};
  for (size_t i = 0; i < v.size(); ++i) acc += pow(abs(v[i]), v.size());
  return pow(acc, 1.0L / v.size());
}

/* *** *** *** max_norm() *** *** *** */
const double max_norm(const Vector& v) {
  double max = v[0];
  for (size_t i = 1; i < v.size(); ++i)
    if (max < v[i]) max = v[i];
  return max;
}

/* *** *** *** orthogonal() *** *** *** */
const bool orthogonal(const Vector& v, const Vector& u) {
  return (dot(v, u) == 0 ? true : false);
}

/* *** *** *** radian() *** *** *** */
const double radian(const Vector& v, const Vector& u) {
  return acosl(dot(v, u) / (l2_norm(v) * l2_norm(u)));
}

/* *** *** *** x_radian() *** *** *** */
const double x_radian(const Vector& v) {
  Vector x{v.size()};
  x[0] = 1;
  x[1] = 0;

  return radian(v, x);
}

/* *** *** *** degree() *** *** *** */
const double degree(const Vector& v, const Vector& u) {
  return (radian(v, u) * 180) / M_PI;
}

/* *** *** *** x_degree() *** *** *** */
const double x_degree(const Vector& v) { return (x_radian(v) * 180) / M_PI; }

/* *** *** *** unit() *** *** *** */
const Vector unit(const Vector& v) {
  Vector u{v.size()};
  const double magnitude = l2_norm(v);
  for (size_t i = 0; i < v.size(); ++i) u[i] = v[i] / magnitude;
  return u;
}

}  // namespace linear
}  // namespace algebra
