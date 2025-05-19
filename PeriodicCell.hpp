#pragma once

#include "mat4.hpp"

struct PeriodicCell {
  //double hxx, hxy;
  //double hyx, hyy;
  mat4r h;

  //double hxx_0, hxy_0;
  //double hyx_0, hyy_0;
  mat4r h0;

  //double vhxx, vhxy;
  //double vhyx, vhyy;
  mat4r vh;

  //double ahxx, ahxy;
  //double ahyx, ahyy;
  mat4r ah;

  double mass;

  PeriodicCell();
  PeriodicCell(double a1x, double a1y, double a2x, double a2y);
  void Define(double a1x, double a1y, double a2x, double a2y);
};
