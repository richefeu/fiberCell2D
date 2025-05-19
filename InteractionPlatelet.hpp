#pragma once

#include "vec2.hpp"

struct Interaction {
  size_t i, j;
  size_t ki, kj;
  double fn, ft, delta_t;

  Interaction();
  Interaction(size_t I, size_t J, size_t KI, size_t KJ);
};
