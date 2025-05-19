#include "PeriodicCell.hpp"

PeriodicCell::PeriodicCell()
    : h(), h0(), vh(), ah() {}

PeriodicCell::PeriodicCell(double a1x, double a1y, double a2x, double a2y)
    : h(), h0(), vh(), ah() {
  Define(a1x, a1y, a2x, a2y);
}

void PeriodicCell::Define(double a1x, double a1y, double a2x, double a2y) {
  h.xx = a1x;
  h.xy = a2x;
  h.yx = a1y;
  h.yy = a2y;
}