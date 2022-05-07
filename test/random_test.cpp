#include "random.h"

#include <gtest/gtest.h>
#include <limits.h>

#include <iostream>

TEST(random_test, no_arg) {
  Random r;
  int number = r.generate();

  ASSERT_TRUE(INT_MIN <= number);
  ASSERT_TRUE(INT_MAX >= number);
}

TEST(randomTest, with_arg) {
  int min = 10, max = 20;

  Random r(min, max);
  int number = r.generate();

  ASSERT_TRUE((min <= number) <= max);
}

TEST(randomTest, binary) {
  int min = 0, max = 1;

  Random r(min, max);
  int number = r.generate();

  ASSERT_TRUE((number == 0) || (number == 1));
}
