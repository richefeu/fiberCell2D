#pragma once

#include <functional>

#include "InteractionPlatelet.hpp"
#include "PeriodicCell.hpp"
#include "Platelet.hpp"
#include "mat4.hpp"

const bool ForceDriven = true;
const bool VelocityDriven = false;

// Loading applied to collective degrees of freedom
struct Loading {
  // 1 for pressure and 0 for velocity
  // bool xxDrive, xyDrive;
  // bool yxDrive, yyDrive;
  mat4<bool> Drive;

  // imposed external stress
  // double Sigxx, Sigxy;
  // double Sigyx, Sigyy;
  mat4r Sig;

  // imposed velocities (TODO rename vhxx...)
  // double vxx, vxy;
  // double vyx, vyy;
  mat4r vh;

  // for path load
  int ipath = 0;
  int Nincr, Ncount;
  char StoredCommand[256];

  // This function will be set to a lambda (c++11)
  std::function<void(Loading &load, PeriodicCell &cell)> ServoFunction; // paramettre -> FiberCell2DSimulation & MD
  std::function<void(Loading &load, PeriodicCell &cell, double &t)> CycleServoFunction;
  std::function<void(Loading &load, double &Sigxx, double &Sigxy, double &Sigyx, double &Sigyy)> IncreServoFunction;
  Loading();
  void BiaxialCompression(double pressure, double velocity);
  void IsostaticCompression(double pressure);
  void SimpleShear(double pressure, double gammaDot);
  void VelocityControl(double Vxx, double Vxy, double Vyx, double Vyy);
  void SimpleShearCycle(double PressureAmp, double gammaDotAmp, double T0);
  void ProportionalBC(double dSigma, double alpha);
  void IncrementalCompression(double dSigma, int Nincr);
  void CyclicSimpleShear(double Pressure, double gammaDotAmp, double T0);
  void VelocityBiaxialCompression(double Vxx, double Vyy);
};
