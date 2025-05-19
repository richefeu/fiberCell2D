#pragma once

#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "InteractionPlatelet.hpp"
#include "Loading.hpp"
#include "PeriodicCell.hpp"
#include "Platelet.hpp"

#include "fileTool.hpp"

// FIXME: rename FiberCell2DSimulation

class FiberCell2DSimulation {
public:
  std::vector<Platelet> Particles;
  std::vector<Interaction> Interactions;
  Loading Load;
  PeriodicCell Cell;
  mat4r Sig;
  double L0;
  double t0;
  double t;
  double tmax;
  double dt;
  double interVerletC, interVerlet;
  double dVerlet;
  double interOutC, interOut;
  double interConfC, interConf;
  int autoVerlet;
  double kn;
  double kt;
  double kb;
  double kc_ff_n;
  double kc_ff_t;
  double kc_fg_n;
  double kc_fg_t;
  double kc_gg_n;
  double kc_gg_t;
  double alpha;
  // ff = fibre-fibre, fg = fibre-grain, gg = grain-grain
  double nu_ff; 
  double nu_fg;
  double nu_gg;
  double density;
  double Surf_p;
  double Cell_Damping;
  int iconf;
  size_t Nfibre, Ngrain;
  std::string result_folder;
  std::string command;

  // Functions
  FiberCell2DSimulation();
  void check();
  void integrate();
  void CundallDamping();
  void ModularTransformation();
  void updatePlateletMoments();
  void accelerations();
  void particleContactForces();
  void FiberInternalForces();
  void ResetCloseList(double dmax);
  void saveConf(int i);
  void loadConf(const char *name);
};

