#pragma once

#include <ostream>

namespace algebra {
namespace linear {

/**
 * @brief Vector class and its operations.
 */
class Vector {
 public:
  /* *** *** *** *** *** *** *** Constructors *** *** *** *** *** *** *** */
  /**
   * @brief C-tor.
   *
   * Create an empty Vector with size of zero.
   */
  Vector() noexcept;

  /**
   * @brief C-tor.
   *
   * @param size The size of vector.
   *
   * Each element of Vector is 0.
   */
  Vector(const std::size_t size) noexcept;

  /**
   * @brief Copy constructor.
   *
   * @param src The source Vector to have a copy.
   *
   * Usage: Vector v{src};
   */
  Vector(const Vector& src) noexcept;

  /**
   * @brief Move constructor.
   *
   * @param src The source Vector to move.
   *
   * Usage: Vector v{move(src)};
   */
  Vector(Vector&& src) noexcept;

  /* *** *** *** *** *** *** *** Operators *** *** *** *** *** *** *** */
  /**
   * @brief Copy-assignment operator.
   *
   * @param rhs The source Vector to make a copy.
   *
   * @return The new Vector.
   *
   * Usage: lhs = rhs.
   * Make a copy of rhs and assignmen to lhs.
   */
  Vector& operator=(const Vector& rhs) noexcept;

  /**
   * @brief Move-assignment operator.
   *
   * @param rhs The source Vector to move.
   *
   * @return The new Vector.
   *
   * Usage: lhs = move(rhs). Move rhs to lhs and empty the rhs.
   */
  Vector& operator=(Vector&& rhs) noexcept;

  /**
   * @brief Brackets operator [] overloading to set a single value at the place.
   *
   * @param i The index of value to change.
   *
   * @return The reference of indexed value.
   *
   * @throws std::out_of_range
   *
   * Usage: v[i] = a. Set the value ar index i of v to a.
   */
  double& operator[](const std::size_t i);

  /**
   * @brief Brackets operator [] overloading to get a single value from the
   * place.
   *
   * @param i The index of value to get.
   *
   * @return The value.
   *
   * @throws std::out_of_range
   *
   * Usage: v[i]. Get the value at index i of v.
   */
  const double operator[](const std::size_t i) const;

  /**
   * @brief Negation operator - overloading.
   *
   * @return The negated version of Vector.
   */
  Vector operator-() const noexcept;

  /**
   * @brief Addition assignment operator += overloading.
   *
   * @param rhs The Vector to add.
   *
   * @return The accumulated Vector.
   *
   * @throws invalid_argument
   *
   * Usage: v +=u. The summation of vectors v and u assigned to v.
   */
  Vector& operator+=(const Vector& rhs);

  /**
   * @brief Subtraction assignment operator -= overloading.
   *
   * @param rhs The Vector to subtract.
   *
   * @return The subtracted Vector.
   *
   * @throws invalid_argument
   *
   * Usage: v -= u. The subtraction of u from v assigned to v.
   *
   */
  Vector& operator-=(const Vector& rhs);

  /**
   * @brief Multiplication assignment operator *= overloading.
   *
   * @param scalar The scalar to multiply the vector.
   *
   * @return The multiplied by scalar vector.
   *
   * Usage: v *= a. The vector multiplied by scalar "a" assigned to v.
   */
  Vector& operator*=(const double scalar) noexcept;

  /**
   * @brief Division assignment operator /= overloading.
   *
   * @param scalar The scalar to divide the vector.
   *
   * @return The divided by scalar vector.
   *
   * @throws logic_error
   *
   * Usage: v /= a. The vector divided by scalar "a" assigned to v.
   */
  Vector& operator/=(const double scalar);

  /* *** *** *** *** *** *** *** Deconstructor *** *** *** *** *** *** *** */
  /**
   * @brief D-tor.
   */
  virtual ~Vector();

  /* *** *** *** *** *** *** Friend Functions *** *** *** *** *** *** */
  /**
   * @brief std::swap() overloading.
   *
   * @param lhs The vector on the left hand side.
   * @param rhs The vector on the right hand side.
   *
   * Usage: swap(v, u). The elements and size of vectors swapped each other.
   */
  friend void swap(Vector& lhs, Vector& rhs) noexcept;

  /**
   * @brief operator<< overloading.
   *
   * @param os The output stream.
   * @param src The source vector to put the os.
   *
   * Usage: stream << v.
   */
  friend std::ostream& operator<<(std::ostream& os, const Vector& src);

  /* *** *** *** *** *** *** *** Member Functions *** *** *** *** *** *** *** */
  /**
   * @brief Get the size of vector.
   *
   * @return The size of vector.
   *
   * Usage: v.size() returns the size of v.
   */
  const std::size_t size() const noexcept;

 private:
  /**
   * @brief The size of Vector.
   */
  std::size_t size_;

  /**
   * @brief The values of Vector.
   */
  double* values_;
};

/**
 * @brief Addition operator + overloading.
 *
 * @param v The first vector.
 * @param u The second vector.
 *
 * @return The summation of two vectors.
 *
 * @throws invalid_argument
 *
 * Usage: m = u + v. The summation of u and v assigned to m.
 */
Vector operator+(const Vector& v, const Vector& u);

// operator-: vector subtraction
/**
 * @brief Subtraction operator - overloading.
 *
 * @param v The first vector.
 * @param u The second vector.
 *
 * @return The subtraction of two vectors.
 *
 * @throws invalid_argument
 *
 * Usage: m = u - v. The subtraction v from u assigned to m.
 */
Vector operator-(const Vector& v, const Vector& u);

/**
 * @brief Scalar vector multiplication operator * overloading.
 *
 * @param scalar The scalar.
 * @param v The vector.
 *
 * @return The vector multiplied by scalar.
 *
 * Usage: u = a * v. The vector v multiplied by "a" and assigned to u.
 */
Vector operator*(const double scalar, const Vector& v);

/**
 * @brief Scalar vector division operator / overloading.
 *
 * @param v The vector.
 * @param scalar The scalar.
 *
 * @return The vector divided by scalar.
 *
 * @throws logic_error
 *
 * Usage: u = v / a. The vector divided by "a" and assigned to u.
 */
Vector operator/(const Vector& v, const double scalar);

/**
 * @brief Equality operator == overloading.
 *
 * @param lhs The vector on the left hand side.
 * @param rhs The vector on the right hand side.
 *
 * @return True if the both vectors are equal, otherwise, false.
 *
 * Usage: u == v. True if u and v are equal, false, otherwise.
 */
const bool operator==(const Vector& lhs, const Vector& rhs);

/**
 * @brief Inequality operator != overloading.
 *
 * @param lhs The vector on the left hand side.
 * @param rhs The vector on the right hand side.
 *
 * @return True if the vectors are not equal, otherwise, true.
 *
 * Usage: u != v. True if u an d are not equal, false, otherwise.
 */
const bool operator!=(const Vector& lhs, const Vector& rhs);

/**
 * @brief Calculate dot production (aka. inner product) of two vectors.
 *
 * @param v The first vector.
 * @param u The second vector.
 *
 * @return doube Dot production of two vectors.
 *
 * @throws invalid_argument
 *
 * Usage: dot(u, v) gives the inner product of u and v.
 */
const double dot(const Vector& v, const Vector& u);

/**
 * @brief \f$l_1-norm\f$ (also called Manhattan distance).
 *
 * @param v The vector.
 *
 * @return The sum of the absolute values of the vector.
 *
 * \f$l_1-norm\f$ calculates the sum of the absolute values of the vector
 * such that \f$||v||_1 = \sum_{i = 1}^n |x_i|\f$.
 */
const double l1_norm(const Vector& v);

/**
 * @brief \f$l_2-norm\f$ (also called Euclidean norm or Euclidean distance).
 *
 * @param v The vector.
 *
 * @return The square root of the sum of the squared vector values (also the
 * length or amplitude of vector).
 *
 * \f$l_2-norm\f$ calculates the square root of the sum of the squared vector
 * values such that \f$||v||_2 = \sum_{i = 1}^n x_i^2\f$.
 */
const double l2_norm(const Vector& v);

/**
 * @brief \f$p-norm\f$.
 *
 * @param v The vector.
 *
 * @return The \f$p^{th}\f$ root of the sum of the \f$p^{th}\f$ power of the
 * vector values.
 *
 * \f$p-norm\f$ calculates the \f$p^{th}\f$ root of the sum of the \f$p^{th}\f$
 * power of the vector values such that \f$||v||_p = (\sum_{i = 1}^n
 * |x_i|^p)^{1/p}\f$.
 */
const double p_norm(const Vector& v);

/**
 * @brief \f$max-norm\f$ (also called \f$l-infinity\f$).
 *
 * @param v The vector.
 *
 * @return The maximum vector of values.
 *
 * The \f$max-norm\f$ calculates the maximum vector values such that
 * \f$||v||_p = \lim_{p \to \infty} (\sum_{i = 1}^n |x_i|^p)^{1/p} = \sup_{1
 * \leq i \leq n} |x_i|\f$.
 */
const double max_norm(const Vector& v);

/**
 * @brief Test the orthogonality of two vectors.
 *
 * @param v The first vector.
 * @param u The second vector.
 *
 * @return True if two vect are orthogonal, oth, false.
 *
 * Two vectors are orthogonal (perpendicular) iff the dot (inner) product of
 * them equals to zero.
 */
const bool orthogonal(const Vector& v, const Vector& u);

/**
 * @brief Find the angle in radians between two vectors.
 *
 * @param v The first vector.
 * @param u The second vector.
 *
 * @return The angle in radians between \f$0\f$ and \f$\pi\f$.
 *
 * We find this angle by using cosine formula:
 * \f$ u \cdot v = ||u||_2 * ||v||_2 * sin(\theta) \f$.
 */
const double radian(const Vector& v, const Vector& u);

/**
 * @brief Find the angle in radians between the given vector and \f$x\f$-axis.
 *
 * @param v The vector.
 *
 * @return The angle in radians between \f$0\f$ and \f$\pi\f$.
 */
const double x_radian(const Vector& v);

/**
 * @brief Find the angle in degrees between two vectors.
 *
 * @param v The first vector.
 * @param u The second vector.
 *
 * @return The angle in degrees between \f$0\f$ and \f$180\f$.
 */
const double degree(const Vector& v, const Vector& u);

/**
 * @brief Find the angle in degrees between the given vector and \f$x\f$-axis.
 *
 * @param v The vector.
 *
 * @return The angle in degrees between \f$0\f$ and \f$180\f$.
 */
const double x_degree(const Vector& v);

/**
 * @brief Find the unit vector of the given vector.
 *
 * @param v The vector.
 *
 * @return The unit vector.
 *
 * The unit vector of vector \f$v\f$ is calculated \f$ \hat{v} = \frac{v}{|v|}
 * \f$, where \f$ |v| \f$ is the magnitude (\f$l_2-norm\f$) of vector
 * \f$v\f$.
 */
const Vector unit(const Vector& v);

}  // namespace linear
}  // namespace algebra
