#include "random.h"

#include <limits.h>

Random::Random() : min_(INT_MIN), max_(INT_MAX) {}

Random::Random(int min, int max) : min_(min), max_(max) {}

Random::~Random() {}

int Random::generate() {
  // only used once to initialise (seed) engine
  std::random_device rd;
  // random-number engine used (Mersenne-Twister in this case)
  std::mt19937 rng(rd());
  // guaranteed unbiased
  std::uniform_int_distribution<int> uni(min_, max_);

  return uni(rng);
}
