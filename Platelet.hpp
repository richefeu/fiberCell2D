#pragma once

#include <cmath>
#include <vector>

#include "vec2.hpp"
#include "mat4.hpp"

// FIXME: rename FibreGrain (?)
struct Platelet {
  size_t n;               // number of nodes
  std::vector<vec2r> pos; // Positions (reduced coordinantes)
  std::vector<vec2r> vel; // Velocities (reduced coordinantes)
  std::vector<vec2r> acc; // Accelerations (reduced coordinantes)

  // displacement since last neighbour list update
  std::vector<vec2r> dpos; // update FIXME: rename cumulatedDisp (for example)
  std::vector<vec2r> prev_pos;

  // rotations
  std::vector<double> rot;  // Angular positions
  std::vector<double> vrot; // Angular velocities
  std::vector<double> arot; // Angular accelerations

  // ___
  std::vector<double> l0;     // initial distances between nodes (not reduced length)
  std::vector<double> dt;     // ___
  std::vector<double> dtheta; // ___

  double radius;  // same radius for all nodes
  double inertia; // inertia of a SINGLE node
  double mass;    // mass of a SINGLE node so that the Platelet mass is n*mass
  double surf0;   // Particle initial surface

  std::vector<vec2r> force;   // resultant force at each node
  std::vector<double> moment; // resultant moment at each node

  double Mxx, Mxy, Myx, Myy; // Tensorial moment TODO: use mat4r

  // Ctors
  Platelet();
  Platelet(size_t N);
  Platelet(vec2r pos1, vec2r pos2, size_t Ns, double R);
  Platelet(vec2r pos, double R);
  
  void getShape(mat4r h, std::vector<vec2r> & coords, double angleStepPI = 0.05);
};
