#pragma once

#include <random>

class Random {
 public:
  Random();
  Random(int min, int max);

  virtual ~Random();

  int generate();

 private:
  int min_, max_;
};
