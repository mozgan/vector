#include "vector.h"

#include <gtest/gtest.h>
#include <limits.h>
#include <math.h>

#include <iostream>

#include "random.h"

using namespace std;
using namespace algebra::linear;

void fill_array(double* arr, size_t size, int min = INT_MIN,
                int max = INT_MAX) {
  Random r(min, max);
  for (size_t i = 0; i < size; ++i) arr[i] = r.generate();
}

TEST(vector_test, constructors) {
  Vector v;
  ASSERT_EQ(v.size(), 0);

  Vector u(19);
  ASSERT_EQ(u.size(), 19);
}

TEST(vector_test, bracket_operators) {
  const double n{989};

  // getter and setter
  Vector v{5};
  v[2] = n;
  ASSERT_EQ(v[2], n);

  // exception in setter
  EXPECT_THROW(v[v.size()] = 10, out_of_range);

  // exception in getter
  EXPECT_THROW(v[v.size() * 2], out_of_range);
}

TEST(vector_test, swap) {
  const size_t n{20};
  double arr[n];
  fill_array(arr, n, 10, 50);

  Vector v(n);
  for (size_t i = 0; i < n; ++i) v[i] = arr[i];

  Vector u{n};
  for (size_t i = 0; i < n; ++i) EXPECT_NE(u[i], arr[i]);

  swap(u, v);

  for (size_t i = 0; i < n; ++i) EXPECT_EQ(u[i], arr[i]);
  for (size_t i = 0; i < n; ++i) EXPECT_NE(v[i], arr[i]);
}

TEST(vector_test, copy_constructor) {
  const size_t n{37};
  Vector v{n};

  Random r(0, 100);
  for (size_t i = 0; i < n; ++i) v[i] = r.generate();

  Vector u{v};
  ASSERT_EQ(u.size(), v.size());
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(u[i], v[i]);
}

TEST(vector_test, move_constructor) {
  const size_t n{100};
  double arr[n];
  fill_array(arr, n, 0, 100);

  Vector v{n};
  for (size_t i = 0; i < n; ++i) v[i] = arr[i];

  Vector u{move(v)};
  ASSERT_EQ(v.size(), 0);
  ASSERT_EQ(u.size(), n);
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(u[i], arr[i]);
}

TEST(vector_test, copy_assignment) {
  const size_t n{93};
  Vector v{n};

  Random r(100, 1000);
  for (size_t i = 0; i < n; ++i) v[i] = r.generate();

  Vector u = v;
  ASSERT_EQ(u.size(), v.size());
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(u[i], v[i]);
}

TEST(vector_test, move_assignment) {
  const size_t n{149};
  double arr[n];
  fill_array(arr, n, 200, 1000);

  Vector v{n};
  for (size_t i = 0; i < n; ++i) v[i] = arr[i];

  Vector u = move(v);
  ASSERT_EQ(v.size(), 0);
  ASSERT_EQ(u.size(), n);
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(u[i], arr[i]);
}

TEST(vector_test, negation_operator) {
  const size_t n{23};
  Random r(0, 100);
  Vector v{n};

  for (size_t i = 0; i < n; ++i) v[i] = r.generate();

  Vector u = -v;
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(u[i], -v[i]);

  Vector k = -(u);
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(k[i], v[i]);
}

TEST(vector_test, vector_addition) {
  const size_t n{10};
  double arr[n];
  fill_array(arr, n, 50, 100);

  Vector v(n), u(n);
  for (size_t i = 0; i < n; ++i) v[i] = u[i] = arr[i];

  // operator+=
  v += u;
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(v[i], 2 * u[i]);

  // operator+
  Vector m = v + u;
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(m[i], 3 * u[i]);

  // exceptions
  Vector k(2 * n);
  EXPECT_THROW(k += v, logic_error);

  Vector z;
  EXPECT_THROW(z = k + v, logic_error);
}

TEST(vector_test, vector_subtraction) {
  const size_t n{20};
  double arr[n];
  fill_array(arr, n, 50, 100);

  Vector v(n), u(n);
  for (size_t i = 0; i < n; ++i) v[i] = u[i] = arr[i];

  // operator-=
  v -= u;
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(v[i], 0);

  // operator-
  Vector m = v - u;
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(m[i], -u[i]);

  // exceptions
  Vector k(2 * n);
  EXPECT_THROW(k -= v, logic_error);

  Vector z;
  EXPECT_THROW(z = k - v, logic_error);
}

TEST(vector_test, scalar_multiplication) {
  const size_t n{20};
  const double scalar{4.59};
  double arr[n];
  fill_array(arr, n, 50, 100);

  Vector v(n);
  for (size_t i = 0; i < n; ++i) v[i] = arr[i];

  v *= scalar;
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(v[i], arr[i] * scalar);

  Vector u = scalar * v;
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(u[i], v[i] * scalar);
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(u[i], arr[i] * scalar * scalar);
}

TEST(vector_test, division_by_scalar) {
  const size_t n{10};
  const double scalar{0.009};
  double arr[n];
  fill_array(arr, n, 50, 100);

  Vector v(n);
  for (size_t i = 0; i < n; ++i) v[i] = arr[i];

  v /= scalar;
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(v[i], arr[i] / scalar);

  Vector u = v / scalar;
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(u[i], v[i] / scalar);
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(u[i], (arr[i] / scalar) / scalar);
}

TEST(vector_test, equality) {
  const size_t n{210};
  Random r(1'000, 10'000);

  Vector v{n}, u{n};
  for (size_t i = 0; i < n; ++i) v[i] = u[i] = r.generate();

  EXPECT_EQ(v, u);
  EXPECT_TRUE(v == u);

  v[0] = -1;
  u[0] = -2;
  EXPECT_NE(v, u);
  EXPECT_FALSE(v == u);
}

TEST(vector_test, dot_product) {
  const size_t n{397};

  Vector v{n}, u{n};
  EXPECT_EQ(dot(v, u), 0);

  for (size_t i = 0; i < n; ++i) v[i] = i;
  EXPECT_EQ(dot(v, u), 0);

  double sum{0};
  for (size_t i = 0; i < n; ++i) {
    sum += i * i;
    u[i] = i;
  }
  EXPECT_EQ(dot(u, v), sum);
}

TEST(vector_test, properties) {
  const size_t n{10};
  const double c{43}, d{98};
  Random r{0, 100};

  Vector u{n}, v{n}, empty{n};
  for (size_t i = 0; i < n; ++i) {
    v[i] = r.generate();
    u[i] = r.generate();
  }

  // u + (-u) = [0...0]
  EXPECT_EQ(u + (-u), empty);

  // (c + d) v = cv + cv
  EXPECT_EQ((c + d) * v, c * v + d * v);

  // c(u + v) = cu + cv
  EXPECT_EQ(c * (u + v), c * u + c * v);
}

TEST(vector_test, l1_norm) {
  const size_t n{11};
  double scalar{-100};
  Random r;
  Vector v{n};

  EXPECT_EQ(0, l1_norm(v));
  v[0] = scalar;
  EXPECT_EQ(abs(scalar), l1_norm(v));

  double acc{0};
  for (size_t i = 0; i < n; ++i) {
    double number = (double)r.generate();
    acc += abs(number);
    v[i] = number;
  }

  EXPECT_EQ(acc, l1_norm(v));
}

TEST(vector_test, l2_norm) {
  const size_t n{25};
  double scalar{-193};
  Random r;
  Vector v{n};

  EXPECT_EQ(0, l2_norm(v));
  v[0] = scalar;
  EXPECT_EQ(abs(scalar), l2_norm(v));

  double acc{0};
  for (size_t i = 0; i < n; ++i) {
    double number = (double)r.generate();
    acc += number * number;
    v[i] = number;
  }

  EXPECT_DOUBLE_EQ(sqrt(acc), l2_norm(v));
}

TEST(vector_test, p_norm) {
  const size_t n{10};
  double scalar{-10'000'000'000};
  Random r;
  Vector v{n};

  EXPECT_EQ(0, p_norm(v));
  v[0] = scalar;
  EXPECT_EQ(-scalar, p_norm(v));

  double acc{0};
  for (size_t i = 0; i < n; ++i) {
    double number = (double)r.generate();
    acc += pow(abs(number), n);
    v[i] = number;
  }

  EXPECT_DOUBLE_EQ(pow(acc, 1.0L / n), p_norm(v));
}

TEST(vector_test, max_norm) {
  const size_t n{1'000};
  const double max{1'111'111.111'111'111};
  Random r(-1'000, 1'000);

  Vector v{n};
  for (size_t i = 0; i < n - 1; ++i) v[i] = r.generate();
  v[n - 1] = max;

  EXPECT_DOUBLE_EQ(max_norm(v), max);
}

TEST(vector_test, orthogonality) {
  Vector v{2}, u{2};

  v[0] = 5;
  v[1] = 4;
  u[0] = 8;
  u[1] = -10;
  EXPECT_TRUE(orthogonal(u, v));

  u[0] = 10;
  u[1] = -10;
  EXPECT_FALSE(orthogonal(u, v));

  Vector m{v.size() * 2};
  EXPECT_THROW(orthogonal(m, v), logic_error);
}

TEST(vector_test, angle_of_two_vectors) {
  const size_t n{3};
  Vector v{n}, u{n};

  v[0] = 1;
  v[1] = 2;
  v[2] = 3;

  u[0] = 3;
  u[1] = -2;
  u[2] = 1;

  // check angle: https://www.omnicalculator.com/math/angle-between-two-vectors
  EXPECT_TRUE(abs(radian(u, v) - 1.42745) < 0.0001);
  EXPECT_TRUE(abs(degree(u, v) - 81.7868) < 0.0001);

  // degree: 45
  v[0] = 1;
  v[1] = v[2] = 0;
  u[0] = 1;
  u[1] = 1;
  u[2] = 0;

  EXPECT_DOUBLE_EQ(45, degree(u, v));
}

TEST(vector_test, angle_to_x_axis) {
  const size_t n{2};
  Vector v{n};

  // degree: 30
  v[0] = sqrt(3);
  v[1] = 1;
  EXPECT_DOUBLE_EQ(30, x_degree(v));

  // degree: 60
  v[0] = 1;
  v[1] = sqrt(3);
  EXPECT_DOUBLE_EQ(60, x_degree(v));

  // degree: 90
  v[0] = 0;
  v[1] = 1;
  EXPECT_DOUBLE_EQ(90, x_degree(v));

  // degree: 180
  v[0] = -1;
  v[1] = 0;
  EXPECT_DOUBLE_EQ(180, x_degree(v));

  // degree: 270 -> 90
  v[0] = 0;
  v[1] = -1;
  EXPECT_DOUBLE_EQ(90, x_degree(v));

  // degree: 360 -> 0
  v[0] = 1;
  v[1] = 0;
  EXPECT_DOUBLE_EQ(0, x_degree(v));
}

TEST(vector_test, angle_of_perpendiculars) {
  const size_t n{2};
  Vector v{n}, u{n};

  v[0] = 1.0L / sqrt(2);
  v[1] = 1.0L / sqrt(2);

  u[0] = -v[0];
  u[1] = v[1];

  // angle of perpendiculars is pi/2 in radians, 90 in degrees
  EXPECT_DOUBLE_EQ(M_PI / 2, radian(u, v));
  EXPECT_DOUBLE_EQ(90, degree(u, v));
}

TEST(vector_test, unit_vector) {
  const size_t n{3};
  Vector v{n};

  v[0] = 12;
  v[1] = -3;
  v[2] = -4;

  // A unit vector u of a vector v is a vector whose length equals one. Then the
  // dot product of u equals 1.
  Vector u{unit(v)};
  EXPECT_DOUBLE_EQ(dot(u, u), 1);
}

TEST(vector_test, sin_cos_unit_vector) {
  const size_t n{2};
  const double theta{59};
  Vector v{n};

  // u = [cos(Θ) sin(Θ)]' is a unit vector.
  v[0] = cos(theta);
  v[1] = sin(theta);

  EXPECT_DOUBLE_EQ(dot(v, v), 1);
}
