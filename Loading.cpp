#include "Loading.hpp"
#include <cmath>
#include <cstdio>
#include <iostream>

Loading::Loading() {}

void Loading::BiaxialCompression(double pressure, double velocity) {
  sprintf(StoredCommand, "BiaxialCompression %g %g", pressure, velocity);
  Drive.xx = true;
  Drive.yy = false;
  Drive.xy = Drive.yx = false;
  Sig.xx = pressure;
  Sig.xy = Sig.yx = Sig.yy = 0.0;
  vh.yy = -velocity;
  vh.xy = vh.yx = 0.0;
  vh.xx = 0.0; // free in fact
  ServoFunction = nullptr;
  CycleServoFunction = nullptr;
  IncreServoFunction = nullptr;
}

void Loading::IsostaticCompression(double pressure) {
  sprintf(StoredCommand, "IsostaticCompression %g", pressure);
  Drive.xx = Drive.yy = ForceDriven;
  Drive.xy = Drive.yx = VelocityDriven;
  Sig.xx = Sig.yy = pressure;
  Sig.xy = Sig.yx = 0.0;
  vh.xx = vh.yy = 0.0; // free in fact
  vh.xy = vh.yx = 0.0;
  ServoFunction = nullptr;
  CycleServoFunction = nullptr;
  IncreServoFunction = nullptr;
}

void Loading::SimpleShear(double pressure, double gammaDot) {
  sprintf(StoredCommand, "SimpleShear %g %g", pressure, gammaDot);
  Drive.xx = Drive.xy = Drive.yx = VelocityDriven;
  Drive.yy = ForceDriven;
  Sig.xx = Sig.xy = Sig.yx = 0.0;
  Sig.yy = pressure;
  vh.xx = vh.yx = vh.yy = 0.0;
  vh.xy = 0.0; // will be driven by the servoFunction
  ServoFunction = [gammaDot](Loading &load, PeriodicCell &cell) -> void { load.vh.xy = gammaDot * cell.h.yy; };
  CycleServoFunction = nullptr;
  IncreServoFunction = nullptr;
}

void Loading::VelocityControl(double Vxx, double Vxy, double Vyx, double Vyy) {
  sprintf(StoredCommand, "VelocityControl %g %g %g %g", Vxx, Vxy, Vyx, Vyy);
  Drive.xx = Drive.xy = Drive.yx = Drive.yy = VelocityDriven;
  Sig.xx = Sig.xy = Sig.yx = Sig.yy = 0.0;
  vh.xx = Vxx;
  vh.xy = Vxy;
  vh.yx = Vyx;
  vh.yy = Vyy;
  ServoFunction = nullptr;
  CycleServoFunction = nullptr;
  IncreServoFunction = nullptr;
}

void Loading::CyclicSimpleShear(double Pressure, double gammaDotAmp, double T0) {
  sprintf(StoredCommand, "CyclicSimpleShear %g %g %g", Pressure, gammaDotAmp, T0);
  Drive.xx = Drive.xy = Drive.yx = VelocityDriven;
  Drive.yy = ForceDriven;
  Sig.xx = Sig.xy = Sig.yx = 0.0;
  vh.xx = vh.yx = vh.yy = 0.0;
  Sig.yy = Pressure;
  vh.xy = 0.0; // will be driven by the servoFunction
  CycleServoFunction = [gammaDotAmp, T0](Loading &load, PeriodicCell &cell, double &t) -> void {
    load.vh.xy = gammaDotAmp * cell.h.yy * cos(M_PI / 2 - 2 * M_PI * t / T0);
  };
  ServoFunction = nullptr;
  IncreServoFunction = nullptr;
}

void Loading::ProportionalBC(double dSigma, double alpha) {
  sprintf(StoredCommand, "ProportionalBiaxialCompression %g %g ", dSigma, alpha);
  Drive.xx = Drive.yy = ForceDriven;
  Drive.xy = Drive.yx = VelocityDriven;
  Sig.xy = Sig.yx = 0.0;
  vh.xx = vh.yy = 0.0; // free in fact
  vh.xy = vh.yx = 0.0;
  IncreServoFunction = [dSigma, alpha](Loading &load, double &Sigxx, double & /*Sigxy*/, double & /*Sigyx*/,
                                       double &Sigyy) -> void {
    if (load.Sig.xx == 0 && load.Sig.yy == 0) {
      std::cout << load.Sig.xx << " " << load.Sig.yy << std::endl;
      load.Sig.xx = Sigxx;
      load.Sig.yy = Sigyy;
    }
    load.Sig.xx += dSigma;
    load.Sig.yy += dSigma * alpha;
  };
  ServoFunction = nullptr;
  CycleServoFunction = nullptr;
}

void Loading::VelocityBiaxialCompression(double Vxx, double Vyy) {
  sprintf(StoredCommand, "VelocityBiaxialCompression %g %g ", Vxx, Vyy);
  Drive.xx = VelocityDriven;
  Drive.yy = VelocityDriven;
  Drive.xy = Drive.yx = ForceDriven;
  Sig.xy = Sig.yx = 0.0;
  vh.xx = Vxx;
  vh.yy = Vyy;
  ServoFunction = nullptr;
  CycleServoFunction = nullptr;
  IncreServoFunction = nullptr;
}

void Loading::IncrementalCompression(double dSigma, int Nincr) {
  sprintf(StoredCommand, "IncrementalCompression %g %d %d ", dSigma, Nincr, 0);
  Drive.xx = Drive.yy = ForceDriven;
  Drive.xy = Drive.yx = VelocityDriven;
  //   Sigxx += dSigma;
  //   Sigyy += dSigma;
  Sig.xy = Sig.yx = 0.0;
  vh.xx = vh.yy = 0.0; // free in fact
  vh.xy = vh.yx = 0.0;
  IncreServoFunction = [dSigma](Loading &load, double &Sigxx, double & /*Sigxy*/, double & /*Sigyx*/,
                                double &Sigyy) -> void {
    if (load.Sig.xx == dSigma && load.Sig.yy == dSigma) {
      std::cout << load.Sig.xx << " " << load.Sig.yy << std::endl;
      load.Sig.xx = Sigxx;
      load.Sig.yy = Sigyy;
    }
    if (load.Ncount > load.Nincr) {
      load.Sig.xx += dSigma;
      load.Sig.yy += dSigma;
      load.Ncount = 0;
    }
    load.Ncount += 1;
  };
  ServoFunction = nullptr;
  CycleServoFunction = nullptr;
}
