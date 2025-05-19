#include "Platelet.hpp"

Platelet::Platelet() {}

Platelet::Platelet(size_t N)
    : n(N), pos(N), vel(N), acc(N), rot(N), vrot(N), arot(N), l0(N - 1), dt(N - 1), dtheta(N - 1), radius(0.0),
      inertia(0.0), mass(1.0), surf0(0.), force(N), moment(N) {
  for (size_t i = 0; i < N - 1; i++) {
    dt[i] = 0.0;
    dtheta[i] = 0.0;
  }
}

Platelet::Platelet(vec2r pos1, vec2r pos2, size_t N, double R)
    : n(N), pos(N), vel(N), acc(N), rot(N), vrot(N), arot(N), l0(N - 1), dt(N - 1), dtheta(N - 1), radius(R),
      inertia(0.0), mass(1.0), surf0(0.), force(N), moment(N) {
  pos[0] = pos1;
  vec2r l = (pos2 - pos1) / (n - 1);
  double L0 = norm(l);
  // surf0 = M_PI * R*R + 2*R * Norm(pos2-pos1);
  surf0 = M_PI * R * R;
  for (size_t i = 0; i < n - 1; i++) {
    pos[i + 1] = pos[i] + l;
    l0[i] = L0;
    dt[i] = 0.0;
    dtheta[i] = 0.0;
  }
}

Platelet::Platelet(vec2r pos1, double R)
    : n(1), pos(n), vel(n), acc(n), rot(n), vrot(n), arot(n), l0(n - 1), dt(n - 1), dtheta(n - 1), radius(R),
      inertia(0.0), mass(1.0), force(n), moment(n) {
  pos[0] = pos1;
  surf0 = M_PI * R * R;
}

void Platelet::getShape(mat4r h, std::vector<vec2r> &coords, double angleStepPI) {

  if (n < 2) {
    vec2r center = h * pos[0];
    for (double angle = 0.0; angle < 2.0 * M_PI; angle += angleStepPI * M_PI) {
      coords.emplace_back(center.x + radius * cos(angle), center.y + radius * sin(angle));
    }
    return;
  }

  // from here n >= 2 ==================================

  // compute disk centers in real coordinates
  std::vector<vec2r> nodes;
  for (size_t i = 0; i < n; i++) {
    nodes.push_back(h * pos[i]);
  }

  // compute the angles between 'bars'
  std::vector<double> angles(n, 0.0);
  for (size_t i = 1; i < n - 1; i++) {
    vec2r prev = nodes[i] - nodes[i - 1];
    vec2r next = nodes[i + 1] - nodes[i];
    angles[i] = angleBetweenVectors(prev, next);
  }

  // compute the unit vector perpendicular to the 'bars'
  std::vector<vec2r> T(n - 1);
  for (size_t i = 0; i < n - 1; i++) {
    vec2r u = nodes[i + 1] - nodes[i];
    u.normalize();
    T[i].set(-u.y, u.x);
  }

  // start-cap
  double incl = inclinationX(nodes[1] - nodes[0]);
  for (double a = incl + 0.5 * M_PI; a < incl + 1.5 * M_PI; a += angleStepPI * M_PI) {
    coords.emplace_back(nodes[0].x + radius * cos(a), nodes[0].y + radius * sin(a));
  }

  // ==> aller
  vec2r newNode;
  newNode = nodes[0] - radius * T[0];
  coords.push_back(newNode);

  for (size_t i = 1; i < n - 1; i++) {
    if (angles[i] < 1e-6) {
      vec2r TT = T[i-1] + T[i];
      TT.normalize();
      newNode = nodes[i] - radius * TT;
      coords.push_back(newNode);
    } else {
      // TODO
      /*
      double angle0 = inclinationX(T[i-1]);
      double angle1 = inclinationX(T[i]);
      double da = (angle1-angle0) / 5.0;
      for (double a = angle0; a <= angle1; a += da) {
        coords.emplace_back(nodes[i].x + radius * cos(a), nodes[i].y + radius * sin(a));
      }
      */
      
      vec2r TT = T[i-1] + T[i];
      TT.normalize();
      newNode = nodes[i] - radius * TT;
      coords.push_back(newNode);
    }
  }

  newNode = nodes[n - 1] - radius * T[n - 2];
  coords.push_back(newNode);

  // end-cap
  incl = inclinationX(nodes[n - 1] - nodes[n - 2]);
  for (double a = incl - 0.5 * M_PI; a < incl + 0.5 * M_PI; a += angleStepPI * M_PI) {
    coords.emplace_back(nodes[n - 1].x + radius * cos(a), nodes[n - 1].y + radius * sin(a));
  }

  // <== retour
  for (size_t i = n-2; i >= 1; i--) {
    if (angles[i] < 1e-6) {
      vec2r TT = T[i-1] + T[i];
      TT.normalize();
      newNode = nodes[i] + radius * TT;
      coords.push_back(newNode);
    } else {
      // TODO
      vec2r TT = T[i-1] + T[i];
      TT.normalize();
      newNode = nodes[i] + radius * TT;
      coords.push_back(newNode);
      
      /*
      double angle0 = inclinationX(T[i]);
      double angle1 = inclinationX(T[i-1]);
      double da = (angle1-angle0) / 5.0;
      for (double a = angle0; a <= angle1; a += da) {
        coords.emplace_back(nodes[i].x + radius * cos(a), nodes[i].y + radius * sin(a));
      }
      */
      
    }
  }

  return;
}
