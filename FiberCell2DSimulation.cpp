#include "FiberCell2DSimulation.hpp"

FiberCell2DSimulation::FiberCell2DSimulation() {
  t = 0.0;
  tmax = 40.0;
  dt = 1e-4;
  interVerletC = 0.0;
  interVerlet = 0.01;
  interOutC = 0.0;
  interOut = 0.02;
  interConfC = 0.0;
  interConf = 0.25;
  autoVerlet = 1;
  kn = 1e3;
  kt = 1e3;
  kb = 1e3;
  // Rigidité de contact ff fibre-fibre, fg fibre-grain, gg grain-grain
  kc_ff_n = 1.e6;
  kc_ff_t = 1.e6;
  kc_fg_n = 1.e6;
  kc_fg_t = 1.e6;
  kc_gg_n = 1.e6;
  kc_gg_t = 1.e6;
  density = 2.7;
  // Coeficient de friction ff fibre-fibre, fg fibre-grain, gg grain-grain
  nu_ff = 0.2;
  nu_fg = 0.2;
  nu_gg = 0.2;
  // Damping Coefficient
  alpha = 0.1;
  Cell_Damping = 0.0;
  iconf = 1;
  result_folder = "./RESULTS";
}

void FiberCell2DSimulation::check() {
  double m = Particles[0].mass;
  double omega0 = sqrt(kn / m);
  std::cout << "dt = " << dt << ", dt/dt0 = " << dt / (1.0 / omega0) << std::endl;
}

// The integration loop (Velocity Verlet)
void FiberCell2DSimulation::integrate() {

  double dt_2 = 0.5 * dt;
  double dt2_2 = 0.5 * dt * dt;

  std::cout << "Beginning iterations..." << std::endl;

  std::ofstream timeStream(result_folder + std::string("/timeStream.txt"));
  std::ofstream Energ(result_folder + std::string("/Energy.txt"));
  Energ << " # 1  2       3        4      5 " << std::endl;
  Energ << " # i  En[i]   Et[i]    Eb[i]  En[i]+Et[i]+Eb[i] " << std::endl;
  timeStream << "#1   2       3        4       5      6      7          8         9 10  11   12       13        14 "
             << std::endl;
  timeStream << "#t   t_exec  Sigxx    Sigxy   Sigyx  Sigyy  Cell_Surf  Compacity  q p  q/p  Ep  Eq En_total Et_total  "
                "Eb_total "
             << std::endl;
  double Compacity = Surf_p / (Cell.h.xx * Cell.h.yy - Cell.h.xy * Cell.h.yx);
  timeStream << t << " " << clock() / 1000.0 << " " << Sig.xx << " " << Sig.xy << " " << Sig.yx << " " << Sig.yy << " "
             << Cell.h.xx * Cell.h.yy - Cell.h.xy * Cell.h.yx << " " << Compacity << " ";
  double En_f = 0.0;
  double Et_f = 0.0;
  double Eb_f = 0.0;

  std::vector<double> En;
  std::vector<double> Et;
  std::vector<double> Eb;
  double En_total = 0.0;
  double Et_total = 0.0;
  double Eb_total = 0.0;

  for (size_t i = 0; i < Particles.size(); i++) {
    if (Particles[i].n > 2) {
      for (size_t p = 0; p < Particles[i].n - 1; p++) {
        vec2r n = Particles[i].pos[p + 1] - Particles[i].pos[p];
        vec2r nReal(Cell.h.xx * n.x + Cell.h.xy * n.y, Cell.h.yx * n.x + Cell.h.yy * n.y);
        // vec2r branch = nReal;
        //  normal
        double l = nReal.normalize();
        double dn = l - Particles[i].l0[p];
        En_f += 0.5 * kn * dn * dn;
        // tangent
        Et_f += 0.5 * kt * Particles[i].dt[p] * Particles[i].dt[p];
        // moment
        Eb_f += 0.5 * kb * Particles[i].dtheta[p] * Particles[i].dtheta[p];
      }
      En.push_back(En_f);
      Et.push_back(Et_f);
      Eb.push_back(Eb_f);

      En_total += En_f;
      Et_total += Et_f;
      Eb_total += Eb_f;
      Energ << i << " " << En[i] << " " << Et[i] << " " << Eb[i] << " " << En[i] + Et[i] + Eb[i] << std::endl;
    }
  }

  double b = -Sig.yy - Sig.xx;
  double c = -Sig.xy * Sig.yx + Sig.xx * Sig.yy;
  double delta = b * b - 4.0 * c;

  double S1 = 0.5 * (-b - sqrt(delta));
  double S2 = 0.5 * (-b + sqrt(delta));

  if (S1 < S2) {
    double swap = S2;
    S2 = S1;
    S1 = swap;
  }
  double q = 0.5 * (S1 - S2);
  double p = 0.5 * (S1 + S2);
  timeStream << q << " " << p << " " << q / p << " " << En_total << " " << Et_total << " " << Eb_total << std::endl;

  for (size_t i = 0; i < Particles.size(); ++i) {
    Particles[i].dpos.resize(Particles[i].n);
    Particles[i].prev_pos.resize(Particles[i].n);

    for (size_t p = 0; p < Particles[i].n; ++p) {
      vec2r realPos(Cell.h.xx * Particles[i].pos[p].x + Cell.h.xy * Particles[i].pos[p].y,
                    Cell.h.yx * Particles[i].pos[p].x + Cell.h.yy * Particles[i].pos[p].y);

      Particles[i].dpos[p] = vec2r(0.0, 0.0);
      Particles[i].prev_pos[p] = Particles[i].pos[p];
    }
  }

  // Limit_Liste();
  ResetCloseList(dVerlet);

  bool updateVerlet;
  double threshold = dVerlet * 0.5;

  while (t < tmax) {

    ModularTransformation();

    // Input parameters for the automatic updating

    if (autoVerlet) {
      updateVerlet = false;

      // displacements for each grain
      for (size_t i = 0; i < Particles.size(); ++i) {
        for (size_t p = 0; p < Particles[i].n; ++p) {

          vec2r delta = Particles[i].pos[p] - Particles[i].prev_pos[p];
          delta.x -= floor(delta.x + 0.5);
          delta.y -= floor(delta.y + 0.5);

          // vec2r delta_real(Cell.h.xx * delta.x + Cell.h.xy * delta.y, Cell.h.yx * delta.x + Cell.h.yy * delta.y);
          vec2r delta_real = Cell.h * delta;
          Particles[i].dpos[p] += delta_real;
          // Particles[i].dpos[p].set(Cell.hxx * delta.x + Cell.hxy * delta.y, Cell.hyx * delta.x + Cell.hyy * delta.y);

          if (norm(Particles[i].dpos[p]) > threshold) {
            updateVerlet = true;
            break; // if it finds one, it breaks the inner loop
          }

          Particles[i].prev_pos[p] = Particles[i].pos[p]; // realPos;
        }
        if (updateVerlet) { // if one has been found, it breaks the outer loop
          break;
        }
      }
    }

    // if the limit is surpased the neighbour list is updated
    if ((autoVerlet != 0 && updateVerlet == true) || (autoVerlet == 0 && interVerletC >= interVerlet)) {

      // std::cout << "Actualizando lista de vecinos..." << std::endl;

      // Limit_Liste();

      ResetCloseList(dVerlet);
      // ResetCloseList_Image(dVerlet);
      interVerletC = 0.0;

      // std::cout << "Lista de vecinos actualizada." << std::endl;

      // Reset of the acumulated desplacement and save the position as prev_position
      if (autoVerlet) {
        for (size_t i = 0; i < Particles.size(); ++i) {
          for (size_t p = 0; p < Particles[i].n; ++p) {
            Particles[i].dpos[p].set(0.0, 0.0);
          }
        }
      }
    }

    if (Load.ServoFunction != nullptr) {
      Load.ServoFunction(Load, Cell);
    }
    if (Load.CycleServoFunction != nullptr) {
      Load.CycleServoFunction(Load, Cell, t);
    }
    if (Load.IncreServoFunction != nullptr) {
      Load.IncreServoFunction(Load, Sig.xx, Sig.xy, Sig.yx, Sig.yy);
    }
    // PlateletXp.clear();
    // PlateletXm.clear();
    // PlateletYm.clear();
    // PlateletYp.clear();
    // PlateletAngGb.clear();
    // PlateletAngGup.clear();
    // PlateletAngDup.clear();
    // PlateletAngDb.clear();
    // LimitsG.clear();
    // LimitsF.clear();

    vec2r C;

    for (size_t i = 0; i < Particles.size(); i++) {
      C.reset();
      for (size_t p = 0; p < Particles[i].n; p++) {

        Particles[i].pos[p] += dt * Particles[i].vel[p] + dt2_2 * Particles[i].acc[p];
        C += Particles[i].pos[p];

        Particles[i].vel[p] += dt_2 * Particles[i].acc[p];
      }
      C /= (double)(Particles[i].n);

      for (size_t p = 0; p < Particles[i].n; p++) {
        Particles[i].rot[p] += dt * Particles[i].vrot[p] + dt2_2 * Particles[i].arot[p];
        Particles[i].vrot[p] += dt_2 * Particles[i].arot[p];
      }

      // Periodicity in position
      // Remark: the reduced velocities do not need to be corrected
      //         since they are periodic
      while (C.x < 0.0) {
        C.x += 1.0;
        for (size_t p = 0; p < Particles[i].n; p++) {
          Particles[i].pos[p].x += 1.0;
        }
      }
      while (C.x > 1.0) {
        C.x -= 1.0;
        for (size_t p = 0; p < Particles[i].n; p++) {
          Particles[i].pos[p].x -= 1.0;
        }
      }
      while (C.y < 0.0) {
        C.y += 1.0;
        for (size_t p = 0; p < Particles[i].n; p++) {
          Particles[i].pos[p].y += 1.0;
        }
      }
      while (C.y > 1.0) {
        C.y -= 1.0;
        for (size_t p = 0; p < Particles[i].n; p++) {
          Particles[i].pos[p].y -= 1.0;
        }
      }

      /*
      double Lmim = 0.2;
      double Lmax = 0.8;
      if (C.x < Lmim) {
        PlateletXm.push_back(i);
        Platelet P = Particles[i];
        for (size_t ip = 0; ip < Particles[i].n; ip++) {
          P.pos[ip].x = 1 + Particles[i].pos[ip].x;
        }
        if (Particles[i].n > 1)
          LimitsF.push_back(P);
        else if (Particles[i].n < 2)
          LimitsG.push_back(P);
      }
      if (C.x > Lmax) {
        PlateletXp.push_back(i);
        if (Particles[i].n > 1)
          LimitsF.push_back(Particles[i]);
        else if (Particles[i].n < 2)
          LimitsG.push_back(Particles[i]);
      }
      if (C.y < Lmim) {
        PlateletYm.push_back(i);

        Platelet P = Particles[i];
        for (size_t ip = 0; ip < Particles[i].n; ip++) {
          P.pos[ip].y = 1 + Particles[i].pos[ip].y;
        }
        if (Particles[i].n > 1)
          LimitsF.push_back(P);
        else if (Particles[i].n < 2)
          LimitsG.push_back(P);
      }
      if (C.y > Lmax) {
        PlateletYp.push_back(i);
        if (Particles[i].n > 1)
          LimitsF.push_back(Particles[i]);
        else if (Particles[i].n < 2)
          LimitsG.push_back(Particles[i]);
      }
      if (C.x < Lmim && C.y < Lmim) {
        PlateletAngGb.push_back(i);
        Platelet P = Particles[i];
        for (size_t ip = 0; ip < Particles[i].n; ip++) {
          P.pos[ip].x = 1 + Particles[i].pos[ip].x;
          P.pos[ip].y = 1 + Particles[i].pos[ip].y;
        }
        if (Particles[i].n > 1)
          LimitsF.push_back(P);
        else if (Particles[i].n < 2)
          LimitsG.push_back(P);
      }
      if (C.x < Lmim && C.y > Lmax) {
        PlateletAngGup.push_back(i);
        if (Particles[i].n > 1)
          LimitsF.push_back(Particles[i]);
        else if (Particles[i].n < 2)
          LimitsG.push_back(Particles[i]);
      }
      if (C.x > Lmax && C.y > Lmax) {
        PlateletAngDup.push_back(i);

        if (Particles[i].n > 1)
          LimitsF.push_back(Particles[i]);
        else if (Particles[i].n < 2)
          LimitsG.push_back(Particles[i]);
      }
      if (C.x > Lmax && C.y < Lmim) {

        PlateletAngDb.push_back(i);
        Platelet P = Particles[i];
        for (size_t ip = 0; ip < Particles[i].n; ip++) {
          P.pos[ip].x = -1 + Particles[i].pos[ip].x;
          P.pos[ip].y = 1 + Particles[i].pos[ip].y;
        }
        if (Particles[i].n > 1)
          LimitsF.push_back(P);
        else if (Particles[i].n < 2)
          LimitsG.push_back(P);
      }
      */
    }

    if (Load.Drive.xx == ForceDriven) {
      Cell.h.xx += dt * Cell.vh.xx + dt2_2 * Cell.ah.xx;
      Cell.vh.xx += dt_2 * Cell.ah.xx;
    } else {
      Cell.h.xx += dt * Load.vh.xx;
      Cell.vh.xx = Load.vh.xx;
      Cell.ah.xx = 0.0;
    }

    if (Load.Drive.xy == ForceDriven) {
      Cell.h.xy += dt * Cell.vh.xy + dt2_2 * Cell.ah.xy;
      Cell.vh.xy += dt_2 * Cell.ah.xy;
    } else {
      Cell.h.xy += dt * Load.vh.xy;
      Cell.vh.xy = Load.vh.xy;
      Cell.ah.xy = 0.0;
    }

    if (Load.Drive.yx == ForceDriven) {
      Cell.h.yx += dt * Cell.vh.yx + dt2_2 * Cell.ah.yx;
      Cell.vh.yx += dt_2 * Cell.ah.yx;
    } else {
      Cell.h.yx += dt * Load.vh.yx;
      Cell.vh.yx = Load.vh.yx;
      Cell.ah.yx = 0.0;
    }

    if (Load.Drive.yy == ForceDriven) {
      Cell.h.yy += dt * Cell.vh.yy + dt2_2 * Cell.ah.yy;
      Cell.vh.yy += dt_2 * Cell.ah.yy;
    } else {
      Cell.h.yy += dt * Load.vh.yy;
      Cell.vh.yy = Load.vh.yy;
      Cell.ah.yy = 0.0;
    }

    accelerations();

    vec2r vmean;
    size_t nbs = 0;
    for (size_t i = 0; i < Particles.size(); i++) {
      nbs += Particles[i].n;
      for (size_t p = 0; p < Particles[i].n; p++) {

        Particles[i].vel[p] += dt_2 * Particles[i].acc[p];
        vmean += Particles[i].vel[p];
        Particles[i].vrot[p] += dt_2 * Particles[i].arot[p];
      }
    }
    vmean /= (double)(nbs);
    for (size_t i = 0; i < Particles.size(); i++) {
      for (size_t p = 0; p < Particles[i].n; p++) {
        Particles[i].vel[p] -= vmean;
      }
    }

    if (Load.Drive.xx == ForceDriven) {
      Cell.vh.xx += dt_2 * Cell.ah.xx;
    }
    if (Load.Drive.xy == ForceDriven) {
      Cell.vh.xy += dt_2 * Cell.ah.xy;
    }
    if (Load.Drive.yx == ForceDriven) {
      Cell.vh.yx += dt_2 * Cell.ah.yx;
    }
    if (Load.Drive.yy == ForceDriven) {
      Cell.vh.yy += dt_2 * Cell.ah.yy;
    }

    // ----

    /*
    if (interVtkC >= interVtk) {
      std::cout << "ivtk = " << ivtk << ", " << 100.0 * (t / tmax) << "%" << std::endl;

      if (ImageCell == 1)
        SaveParticles_GrainVtk_Limite(Cell, LimitsG, ivtk, result_folder.c_str(), 24, LimitsG.size());
      if (ImageCell == 1)
        SaveParticlesVtk_Limite(Cell, LimitsF, ivtk, result_folder.c_str(), 24, LimitsF.size());

      SaveParticlesVtk(Cell, Particles, ivtk, result_folder.c_str(), 24, Nfibre);
      SaveParticles_GrainVtk(Cell, Particles, ivtk, result_folder.c_str(), 24, Ngrain);
      SaveCellVtk(Cell, ivtk, result_folder.c_str());
      ivtk++;
      interVtkC = 0.0;
    }
*/

    if (interOutC >= interOut) {
      // 1    2     3     4     5     6  7 8 9 10
      // time Sigxx Sigxy Sigyx Sigyy V  e q p q/p
      std::vector<double> En;
      std::vector<double> Et;
      std::vector<double> Eb;
      double En_total = 0.;
      double Et_total = 0.;
      double Eb_total = 0.;
      for (size_t i = 0; i < Particles.size(); i++) {
        if (Particles[i].n > 2) {
          double En_f = 0.;
          double Et_f = 0.;
          double Eb_f = 0.;

          for (size_t p = 0; p < Particles[i].n - 1; p++) {
            vec2r n = Particles[i].pos[p + 1] - Particles[i].pos[p];
            vec2r nReal(Cell.h.xx * n.x + Cell.h.xy * n.y, Cell.h.yx * n.x + Cell.h.yy * n.y);
            // vec2r branch = nReal;
            //  normal
            double l = nReal.normalize();
            double dn = l - Particles[i].l0[p];
            En_f += 0.5 * kn * dn * dn;
            // tangent
            Et_f += 0.5 * kt * Particles[i].dt[p] * Particles[i].dt[p];
            // moment
            Eb_f += 0.5 * kb * Particles[i].dtheta[p] * Particles[i].dtheta[p];
          }
          En.push_back(En_f);
          Et.push_back(Et_f);
          Eb.push_back(Eb_f);

          En_total += En_f;
          Et_total += Et_f;
          Eb_total += Eb_f;
          Energ << i << " " << En[i] << " " << Et[i] << " " << Eb[i] << " " << En[i] + Et[i] + Eb[i] << std::endl;
        }
      }
      // SaveParticlesVtk_Energy(Cell, Particles, ivtk - 1, result_folder.c_str(), 24, Nfibre, En, Et, Eb);

      double Compacity = Surf_p / (Cell.h.xx * Cell.h.yy - Cell.h.xy * Cell.h.yx);
      timeStream << t << " " << clock() / 1000.
                 << " "
                 // timeStream << t << " " << clock() / (double)CLOCKS_PER_SEC << " "
                 << Sig.xx << " " << Sig.xy << " " << Sig.yx << " " << Sig.yy << " "
                 << Cell.h.xx * Cell.h.yy - Cell.h.xy * Cell.h.yx << " " << Compacity << " ";
      double b = -Sig.yy - Sig.xx;
      double c = -Sig.xy * Sig.yx + Sig.xx * Sig.yy;
      double delta = b * b - 4.0 * c;

      double S1 = 0.5 * (-b - sqrt(delta));
      double S2 = 0.5 * (-b + sqrt(delta));

      if (S1 < S2) {
        double swap = S2;
        S2 = S1;
        S1 = swap;
      }
      double q = 0.5 * (S1 - S2);
      double p = 0.5 * (S1 + S2);
      timeStream << q << " " << p << " " << q / p << " " << En_total << " " << Et_total << " " << Eb_total << std::endl;

      interOutC = 0.0;
    }

    if (interConfC >= interConf) {
      std::cout << "Save Conf" << iconf << std::endl;
      saveConf(iconf);
      iconf++;
      interConfC = 0.0;
    }

    interConfC += dt;
    interOutC += dt;
    // interVtkC += dt;
    interVerletC += dt;
    t += dt;
  }
  return;
}

void FiberCell2DSimulation::updatePlateletMoments() {
  for (size_t i = 0; i < Particles.size(); i++) {
    Particles[i].Mxx = 0.0;
    Particles[i].Mxy = 0.0;
    Particles[i].Myx = 0.0;
    Particles[i].Myy = 0.0;
  }

  size_t i, j, pi, pj;
  for (size_t k = 0; k < Interactions.size(); k++) {
    i = Interactions[k].i;
    j = Interactions[k].j;
    pi = Interactions[k].ki;
    pj = Interactions[k].kj;

    vec2r sij = Particles[j].pos[pj] - Particles[i].pos[pi];
    // sij.x -= floor(sij.x + 0.5);
    // sij.y -= floor(sij.y + 0.5);
    vec2r branch(Cell.h.xx * sij.x + Cell.h.xy * sij.y, Cell.h.yx * sij.x + Cell.h.yy * sij.y);
    vec2r posReal(Cell.h.xx * Particles[i].pos[pi].x + Cell.h.xy * Particles[i].pos[pi].y,
                  Cell.h.yx * Particles[i].pos[pi].x + Cell.h.yy * Particles[i].pos[pi].y);
    vec2r n = branch;
    n.normalize();
    branch *= 0.5;
    vec2r pos = (posReal + branch);
    // vec2r posi = pos - Particles[i].pos[0];
    // vec2r posj = pos - Particles[j].pos[0];
    vec2r f = Interactions[k].fn * n;

    Particles[i].Mxx += f.x * pos.x;
    Particles[i].Mxy += f.x * pos.y;
    Particles[i].Myx += f.y * pos.x;
    Particles[i].Myy += f.y * pos.y;

    Particles[j].Mxx += -f.x * pos.x;
    Particles[j].Mxy += -f.x * pos.y;
    Particles[j].Myx += -f.y * pos.x;
    Particles[j].Myy += -f.y * pos.y;
  }
}

void FiberCell2DSimulation::ModularTransformation() {
  // this is ok only for shearing (should be embbeded in Loading???)
  if (Cell.h.xy <= Cell.h.xx) {
    return;
  }

  // Save the current periodic cell
  mat4r prevh = Cell.h;

  // Make the modular transformation of the cell
  // double dx = prevh.xy;
  // Cell.h.xy -= dx;
  Cell.h.xy = 0.0;

  mat4r newhinv = Cell.h.get_inverse();

  vec2r C;
  for (size_t i = 0; i < Particles.size(); i++) {

    // barycenter of the fiber (reduced coordinates)
    C.reset();
    for (size_t p = 0; p < Particles[i].n; p++) {
      C += prevh * Particles[i].pos[p];
    }
    C /= (double)(Particles[i].n);

    //vec2r newRealCenterPos = (prevh * C);

    if (C.x > prevh.xx) {
      for (size_t p = 0; p < Particles[i].n; p++) {
        vec2r prevRealPos = prevh * Particles[i].pos[p];
        prevRealPos.x -= prevh.xx;
        Particles[i].pos[p] = newhinv * prevRealPos;
      }
    } else {
      for (size_t p = 0; p < Particles[i].n; p++) {
        vec2r prevRealPos = prevh * Particles[i].pos[p];
        Particles[i].pos[p] = newhinv * prevRealPos;
      }
    }

    /*
    if (C.y > 1.0 - C.x) { // upper part
      for (size_t p = 0; p < Particles[i].n; p++) {
        vec2r prevRealPos = prevh * Particles[i].pos[p];
        prevRealPos.x -= prevh.xx;
        Particles[i].pos[p] = newhinv * prevRealPos;
      }
    } else { // lower part
      for (size_t p = 0; p < Particles[i].n; p++) {
        vec2r prevRealPos = prevh * Particles[i].pos[p];
        Particles[i].pos[p] = newhinv * prevRealPos;
      }
    }
    */
  }
}

void FiberCell2DSimulation::CundallDamping() {
  if (Load.Drive.xx == ForceDriven) {
    if (Cell.ah.xx * Cell.vh.xx >= 0.0) {
      Cell.ah.xx *= (1.0 - Cell_Damping);
    } else {
      Cell.ah.xx *= (1.0 + Cell_Damping);
    }
  }

  if (Load.Drive.yy == ForceDriven) {
    if (Cell.ah.yy * Cell.vh.yy >= 0.0) {
      Cell.ah.yy *= (1.0 - Cell_Damping);
    } else {
      Cell.ah.yy *= (1.0 + Cell_Damping);
    }
  }
}

void FiberCell2DSimulation::ResetCloseList(double dmax) {

  // First we need to store...
  std::vector<Interaction> Ibak;
  Interaction I;
  for (size_t k = 0; k < Interactions.size(); k++) {
    I = Interactions[k];
    Ibak.push_back(I);
  }

  // ... Now we rebuild the list
  for (size_t i = 0; i < Particles.size(); i++) {
    for (size_t j = i + 1; j < Particles.size(); j++) {

      for (size_t ki = 0; ki < Particles[i].n; ki++) {
        for (size_t kj = 0; kj < Particles[j].n; kj++) {

          vec2r sij = Particles[j].pos[kj] - Particles[i].pos[ki];
          sij.x -= floor(sij.x + 0.5);
          sij.y -= floor(sij.y + 0.5);

          vec2r branch(Cell.h.xx * sij.x + Cell.h.xy * sij.y, Cell.h.yx * sij.x + Cell.h.yy * sij.y);
          vec2r n = branch;
          double Lij = n.normalize();
          double dn = Lij - Particles[i].radius - Particles[j].radius;
          if (dn <= dmax) {
            Interactions.push_back(Interaction(i, j, ki, kj));
          }

        } // pj
      } // pi

    } // j
  } // i

  // And finally, we retrieve the data
  // retrieve previous contacts or bonds
  size_t k, kold = 0;
  for (k = 0; k < Interactions.size(); ++k) {
    while (kold < Ibak.size() && Ibak[kold].i < Interactions[k].i) {
      ++kold;
    }
    if (kold == Ibak.size()) {
      break;
    }

    while (kold < Ibak.size() && Ibak[kold].i == Interactions[k].i && Ibak[kold].j < Interactions[k].j) {
      ++kold;
    }
    if (kold == Ibak.size()) {
      break;
    }

    while (kold < Ibak.size() && Ibak[kold].i == Interactions[k].i && Ibak[kold].j == Interactions[k].j &&
           Ibak[kold].ki < Interactions[k].ki) {
      ++kold;
    }
    if (kold == Ibak.size()) {
      break;
    }

    while (kold < Ibak.size() && Ibak[kold].i == Interactions[k].i && Ibak[kold].j == Interactions[k].j &&
           Ibak[kold].ki == Interactions[k].ki && Ibak[kold].kj < Interactions[k].kj) {
      ++kold;
    }
    if (kold == Ibak.size()) {
      break;
    }

    if (Ibak[kold].i == Interactions[k].i && Ibak[kold].j == Interactions[k].j && Ibak[kold].ki == Interactions[k].ki &&
        Ibak[kold].kj == Interactions[k].kj) {
      Interactions[k] = Ibak[kold];
      ++kold;
    }
  }
}

void FiberCell2DSimulation::particleContactForces() {

  size_t i, j, ki, kj;

  for (size_t k = 0; k < Interactions.size(); k++) {
    i = Interactions[k].i;
    j = Interactions[k].j;
    ki = Interactions[k].ki;
    kj = Interactions[k].kj;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vec2r sij = Particles[j].pos[kj] - Particles[i].pos[ki];
    sij.x -= floor(sij.x + 0.5);
    sij.y -= floor(sij.y + 0.5);

    vec2r vij = Particles[j].vel[kj] - Particles[i].vel[ki];
    vec2r branch(Cell.h.xx * sij.x + Cell.h.xy * sij.y, Cell.h.yx * sij.x + Cell.h.yy * sij.y);
    vec2r Vij(Cell.h.xx * vij.x + Cell.h.xy * vij.y, Cell.h.yx * vij.x + Cell.h.yy * vij.y);

    double Rij = (Particles[i].radius + Particles[j].radius);

    double dn = norm(branch) - Rij;
    double kc_n = 0.0, kc_t = 0.0, nu = 0.0;

    if (Particles[i].n == 1 && Particles[j].n == 1) {
      kc_n = kc_gg_n;
      kc_t = kc_gg_t;
      nu = nu_gg;
    }

    if (Particles[i].n == 1 && Particles[j].n > 1) {
      kc_n = kc_fg_n;
      kc_t = kc_fg_t;
      nu = nu_fg;
    }
    if (Particles[i].n > 1 && Particles[j].n == 1) {
      kc_n = kc_fg_n;
      kc_t = kc_fg_t;
      nu = nu_fg;
    }

    if (Particles[i].n > 1 && Particles[j].n > 1) {
      kc_n = kc_ff_n;
      kc_t = kc_ff_t;
      nu = nu_ff;
    }
    if (dn < 0.0) {
      // std::cout << "i= " << i << " j = " << j << " Ri = " << Particles[i].radius << " Rj = " << Particles[j].radius
      // << std::endl;

      // std::cout<< "dn= "<< dn << " Rij = "<< Rij << " Branch= " << Norm(branch) << std::endl;
      // Normal force (elastic + viscuous)
      vec2r n = branch.normalized();
      vec2r vt(-n.y, n.x);
      double m1 = Particles[i].mass;
      double m2 = Particles[j].mass;
      double m_eff = (m1 * m2) / (m1 + m2);
      double eta = alpha * sqrt(m_eff * kc_n);
      vec2r Vnij = Vij - (Vij * vt) * vt;
      double fn = -kc_n * dn - eta * (Vnij * n);
      // double fn = -kc_n * dn;
      if (fn < 0.)
        fn = 0;
      double Vrotij = Particles[i].radius * Particles[i].vrot[ki] - Particles[j].radius * Particles[j].vrot[kj];
      vec2r Vsij = Vij - (Vij * n) * n + Vrotij * vt;
      // vec2r Vsij = Vij - (Vij * n) * n;
      vec2r delta_t = Vsij * dt;
      Interactions[k].delta_t -= delta_t * vt;
      Interactions[k].fn = fn;
      // double eta_t = alpha*sqrt(m_eff * kc_t);
      // double Ft = kc_t * Interactions[k].delta_t + vsij * eta_t*vt;
      // Interactions[k].ft = fmin(fabs(Ft), fabs(nu*fn));

      // Tangential force (friction)
      double eta_t = alpha * sqrt(m_eff * kc_t);
      double Ft = kc_t * Interactions[k].delta_t - eta_t * (Vsij * vt);
      double sig = (fabs(Vsij * vt) > 1e-20) ? -1 * (Vsij * vt) / fabs(Vsij * vt) : 0.0;
      double FT = sig * fmin(fabs(Ft), nu * fn);
      Interactions[k].ft = FT;

      // Resultant force
      vec2r f = Interactions[k].fn * n + Interactions[k].ft * vt;

      // #pragma omp critical
      {
        Particles[i].force[ki] -= f;
        Particles[j].force[kj] += f;
        Particles[i].moment[ki] += Interactions[k].ft * Particles[i].radius;
        Particles[j].moment[kj] -= Interactions[k].ft * Particles[j].radius;
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // Internal stress
        Sig.xx += f.x * branch.x;
        Sig.xy += f.x * branch.y;
        Sig.yx += f.y * branch.x;
        Sig.yy += f.y * branch.y;
      }
    }
    if (dn >= 0.0) {
      Interactions[k].delta_t = 0.0;
    }
  } // Loop over interactions k
}

void FiberCell2DSimulation::accelerations() {
  // Set forces and moments to zero
  for (size_t i = 0; i < Particles.size(); i++) {
    for (size_t p = 0; p < Particles[i].n; p++) {
      Particles[i].force[p].reset();
      Particles[i].moment[p] = 0.0;
    }
  }

  Sig.reset();

  // Compute the forces between ff/fg/gg
  particleContactForces();

  // Compute the forces related to the flexibility of the fibers
  FiberInternalForces();

  double invV = 1.0 / (Cell.h.xx * Cell.h.yy - Cell.h.yx * Cell.h.xy);
  Sig.xx *= invV;
  Sig.xy *= invV;
  Sig.yx *= invV;
  Sig.yy *= invV;

  // Finally compute the accelerations (translation and rotation)
  for (size_t i = 0; i < Particles.size(); i++) {
    for (size_t p = 0; p < Particles[i].n; p++) {
      Particles[i].arot[p] = Particles[i].moment[p] / Particles[i].inertia;
      Particles[i].acc[p] = Particles[i].force[p] / Particles[i].mass;
      // The following 4 lines may actually be removed
      // Particles[i].acc[p].x -= 2.0 * (Cell.vhxx * Particles[i].vel[p].x + Cell.vhxy * Particles[i].vel[p].y);
      // Particles[i].acc[p].y -= 2.0 * (Cell.vhyx * Particles[i].vel[p].x + Cell.vhyy * Particles[i].vel[p].y);
      // Be carreful: if the system is globally damped (Cundall damping), ah may not be okay
      // Particles[i].acc[p].x -= Cell.ahxx * Particles[i].pos[p].x + Cell.ahxy * Particles[i].pos[p].y;
      // Particles[i].acc[p].y -= Cell.ahyx * Particles[i].pos[p].x + Cell.ahyy * Particles[i].pos[p].y;
      // Compute inverse of the matrix h
      double invDet = 1.0 / (Cell.h.xx * Cell.h.yy - Cell.h.yx * Cell.h.xy);
      double hxxinv = invDet * Cell.h.yy;
      double hxyinv = -invDet * Cell.h.xy;
      double hyxinv = -invDet * Cell.h.yx;
      double hyyinv = invDet * Cell.h.xx;

      vec2r acc = Particles[i].acc[p];
      Particles[i].acc[p].x = hxxinv * acc.x + hxyinv * acc.y;
      Particles[i].acc[p].y = hyxinv * acc.x + hyyinv * acc.y;
    }
  }

  if (Load.Drive.xx == ForceDriven) {
    Cell.ah.xx = ((Sig.xx - Load.Sig.xx) * Cell.h.yy - (Sig.yx - Load.Sig.yx) * Cell.h.xy) / Cell.mass;
  } else {
    Cell.ah.xx = 0.0;
  }

  if (Load.Drive.xy == ForceDriven) {
    Cell.ah.xy = ((Sig.xy - Load.Sig.xy) * Cell.h.yy - (Sig.yy - Load.Sig.yy) * Cell.h.xy) / Cell.mass;
  } else {
    Cell.ah.xy = 0.0;
  }

  if (Load.Drive.yx == ForceDriven) {
    Cell.ah.yx = ((Sig.yx - Load.Sig.yx) * Cell.h.xx - (Sig.xx - Load.Sig.xx) * Cell.h.yx) / Cell.mass;
  } else {
    Cell.ah.yx = 0.0;
  }

  if (Load.Drive.yy == ForceDriven) {
    Cell.ah.yy = ((Sig.yy - Load.Sig.yy) * Cell.h.xx - (Sig.xy - Load.Sig.xy) * Cell.h.yx) / Cell.mass;
  } else {
    Cell.ah.yy = 0.0;
  }

  // Globally damp the system
  if (Cell_Damping != 0.0) {
    CundallDamping();
  }
}
void FiberCell2DSimulation::FiberInternalForces() {
  for (size_t i = 0; i < Particles.size(); i++) {
    if (Particles[i].n > 1) {
      for (size_t p = 0; p < Particles[i].n - 1; p++) {
        vec2r n = Particles[i].pos[p + 1] - Particles[i].pos[p];
        vec2r nReal(Cell.h.xx * n.x + Cell.h.xy * n.y, Cell.h.yx * n.x + Cell.h.yy * n.y);
        vec2r branch = nReal;

        // normal
        double l = nReal.normalize();
        double dn = l - Particles[i].l0[p];
        vec2r fn = kn * dn * nReal;

        // tangent
        vec2r tReal(-nReal.y, nReal.x);
        vec2r vt = (Particles[i].vel[p + 1] - Particles[i].vel[p]);
        vec2r vtReal(Cell.h.xx * vt.x + Cell.h.xy * vt.y + Cell.vh.xx * n.x + Cell.vh.xy * n.y,
                     Cell.h.yx * vt.x + Cell.h.yy * vt.y + Cell.vh.yx * n.x + Cell.vh.yy * n.y);

        vtReal -= 0.5 * Particles[i].l0[p] * (Particles[i].vrot[p + 1] + Particles[i].vrot[p]) *
                  tReal; // absolument necessaire
        Particles[i].dt[p] += vtReal * tReal * dt;
        double ft = kt * Particles[i].dt[p];

        // moment
        Particles[i].dtheta[p] += (Particles[i].vrot[p + 1] - Particles[i].vrot[p]) * dt;
        double mom = kb * Particles[i].dtheta[p];

        vec2r f = fn + ft * tReal;
        Particles[i].force[p] += f;
        Particles[i].force[p + 1] -= f;
        l *= 0.5;
        Particles[i].moment[p] += l * ft + mom;
        Particles[i].moment[p + 1] += l * ft - mom;

        // Internal stress
        Sig.xx -= f.x * branch.x;
        Sig.xy -= f.x * branch.y;
        Sig.yx -= f.y * branch.x;
        Sig.yy -= f.y * branch.y;
      }
    }

  } // end for i
}

void FiberCell2DSimulation::saveConf(int i) {
  char fname[256];
  sprintf(fname, "%s/conf%d", result_folder.c_str(), i);
  std::ofstream conf(fname);

  conf << "claySOMEFv2 09-02-2015" << std::endl; // format: progName version-date
  conf << "result_folder " << result_folder << std::endl;
  conf << "t " << t << std::endl;
  // conf << "t_exec " << clock() / 1000.0 << std::endl;
  conf << "tmax " << tmax << std::endl;
  conf << "dt " << dt << std::endl;
  conf << "interVerlet " << interVerlet << std::endl;
  conf << "interOut " << interOut << std::endl;
  conf << "interConf " << interConf << std::endl;
  conf << "dVerlet " << dVerlet << std::endl;
  conf << "autoVerlet " << autoVerlet << std::endl;
  conf << "kn " << kn << std::endl;
  conf << "kt " << kt << std::endl;
  conf << "kb " << kb << std::endl;
  conf << "kc_ff_n " << kc_ff_n << std::endl;
  conf << "kc_ff_t " << kc_ff_t << std::endl;
  conf << "kc_fg_n " << kc_fg_n << std::endl;
  conf << "kc_fg_t " << kc_fg_t << std::endl;
  conf << "kc_gg_n " << kc_gg_n << std::endl;
  conf << "kc_gg_t " << kc_gg_t << std::endl;
  conf << "nu_ff " << nu_ff << std::endl;
  conf << "nu_fg " << nu_fg << std::endl;
  conf << "nu_gg " << nu_gg << std::endl;
  conf << "alpha " << alpha << std::endl;
  conf << "density " << density << std::endl;
  conf << "Cell_Damping " << Cell_Damping << std::endl;
  conf << "iconf " << iconf << std::endl;
  conf << std::scientific << std::setprecision(5);
  conf << "Cell " << Cell.h << std::endl;
  conf << "vCell " << Cell.vh << std::endl;
  conf << "aCell " << Cell.ah << std::endl;
  conf.unsetf(std::ios_base::floatfield); // because conf << std::defaultfloat; doesn't work with gcc-4.9
  conf << "massCell " << Cell.mass << std::endl;
  conf << "Load " << Load.StoredCommand << std::endl;
  conf << "Stress " << Sig << std::endl;
  conf << std::scientific << std::setprecision(5);
  conf << "Platelets " << Particles.size() << std::endl;
  for (size_t i = 0; i < Particles.size(); i++) {
    conf << Particles[i].n << " " << Particles[i].radius << " " << Particles[i].inertia << " " << Particles[i].mass
         << std::endl;
    if (Particles[i].n > 1) {
      for (size_t p = 0; p < Particles[i].n - 1; p++) {
        conf << Particles[i].l0[p] << " " << Particles[i].dt[p] << " " << Particles[i].dtheta[p] << " ";
      }
      conf << std::endl;
    }
    for (size_t p = 0; p < Particles[i].n; p++) {
      conf << Particles[i].pos[p] << " " << Particles[i].vel[p] << " " << Particles[i].acc[p] << " "
           << Particles[i].rot[p] << " " << Particles[i].vrot[p] << " " << Particles[i].arot[p] << " " << std::endl;
    }
  }

  size_t nI = 0;
  for (size_t i = 0; i < Interactions.size(); i++) {
    if (fabs(Interactions[i].fn) < 1.e-20) {
      continue;
    }
    nI++;
  }
  conf << "Interactions " << nI << std::endl;
  for (size_t i = 0; i < Interactions.size(); i++) {
    if (fabs(Interactions[i].fn) < 1.e-20) {
      continue;
    }
    conf << Interactions[i].i << " " << Interactions[i].j << " " << Interactions[i].ki << " " << Interactions[i].kj
         << " " << Interactions[i].fn << " " << Interactions[i].ft << " " << Interactions[i].delta_t << " "
         << std::endl;
  }
}

void FiberCell2DSimulation::loadConf(const char *name) {
  std::ifstream conf(name);
  if (!conf.is_open()) {
    std::cerr << "Cannot read " << name << std::endl;
  }
  // Check header
  std::string prog;
  conf >> prog;
  if (prog != "claySOMEFv2") {
    std::cerr << "This is not file for claySOMEFv2 executable!" << std::endl;
  }
  std::string date;
  conf >> date;
  if (date != "09-02-2015") {
    std::cerr << "The version-date should be 09-02-2015!" << std::endl;
  }
  std::string token;
  conf >> token;
  while (conf.good()) {
    if (token == "result_folder")
      conf >> result_folder;
    else if (token == "t")
      conf >> t;
    else if (token == "t_exec") { // TODO: remove
      double t_exec;
      conf >> t_exec;
      std::cerr << "'t_exec' is deprecated\n";
    } else if (token == "dcut") { // TODO: remove
      double dcut;
      conf >> dcut;
      std::cerr << "'dcut' is deprecated\n";
    } else if (token == "ImageCell") { // TODO: remove
      int ImageCell;
      conf >> ImageCell;
      std::cerr << "'ImageCell' is deprecated\n";
    } else if (token == "tmax") {
      conf >> tmax;
    } else if (token == "dt") {
      conf >> dt;
    } else if (token == "interVerlet") {
      conf >> interVerlet;
    } else if (token == "interOut") {
      conf >> interOut;
    } else if (token == "interConf") {
      conf >> interConf;
    } else if (token == "dVerlet") {
      conf >> dVerlet;
    } else if (token == "autoVerlet") {
      conf >> autoVerlet;
    } else if (token == "kn") {
      conf >> kn;
    } else if (token == "kt") {
      conf >> kt;
    } else if (token == "kb") {
      conf >> kb;
    } else if (token == "kc_ff_n") {
      conf >> kc_ff_n;
    } else if (token == "kc_ff_t") {
      conf >> kc_ff_t;
    } else if (token == "kc_fg_n") {
      conf >> kc_fg_n;
    } else if (token == "kc_fg_t") {
      conf >> kc_fg_t;
    } else if (token == "kc_gg_n") {
      conf >> kc_gg_n;
    } else if (token == "kc_gg_t") {
      conf >> kc_gg_t;
    } else if (token == "nu_ff") {
      conf >> nu_ff;
    } else if (token == "nu_fg") {
      conf >> nu_fg;
    } else if (token == "nu_gg") {
      conf >> nu_gg;
    } else if (token == "alpha") {
      conf >> alpha;
    } else if (token == "density") {
      conf >> density;
    } else if (token == "Cell_Damping") {
      conf >> Cell_Damping;
    } else if (token == "iconf") {
      conf >> iconf;
    } else if (token == "Cell") {
      conf >> Cell.h;
      Cell.h0 = Cell.h;
    } else if (token == "vCell")
      conf >> Cell.vh;
    else if (token == "aCell")
      conf >> Cell.ah;
    else if (token == "massCell")
      conf >> Cell.mass;
    // else if (token == "ImageCell")
    //  conf >> ImageCell;
    else if (token == "fibres") {
      size_t Ns;
      conf >> Nfibre >> Ns;
      for (size_t i = 0; i < Nfibre; i++) {
        vec2r pos1, pos2;
        double R;
        conf >> pos1.x >> pos1.y >> pos2.x >> pos2.y >> R;
        Particles.push_back(Platelet(pos1, pos2, Ns, R));
      }
    } else if (token == "grains") {
      conf >> Ngrain;
      for (size_t i = 0; i < Ngrain; i++) {
        vec2r pos;
        double R;
        conf >> pos.x >> pos.y >> R;
        Particles.push_back(Platelet(pos, R));
      }
      // reduced values
      double invDet = 1.0 / (Cell.h.xx * Cell.h.yy - Cell.h.yx * Cell.h.xy);
      double hxxinv = invDet * Cell.h.yy;
      double hxyinv = -invDet * Cell.h.xy;
      double hyxinv = -invDet * Cell.h.yx;
      double hyyinv = invDet * Cell.h.xx;
      double TotalMass = 0.0;
      double PI = 3.141592653589793;
      Surf_p = 0.0;
      for (size_t i = 0; i < Particles.size(); i++) {
        double R2 = Particles[i].radius * Particles[i].radius;
        Particles[i].mass = PI * density * R2; // a mettree à jour avec la densité
        Particles[i].inertia = 0.5 * Particles[i].mass * R2;
        if (Particles[i].n > 1)
          Particles[i].mass = PI * density * R2 / 2.;
        TotalMass += double(Particles[i].n) * Particles[i].mass;
        // Surf_p += Particles[i].surf0;
        Surf_p += double(Particles[i].n) * Particles[i].surf0;
        for (size_t p = 0; p < Particles[i].n; p++) {
          vec2r pos = Particles[i].pos[p];
          Particles[i].pos[p].x = hxxinv * pos.x + hxyinv * pos.y;
          Particles[i].pos[p].y = hyxinv * pos.x + hyyinv * pos.y;
          Particles[i].rot[p] = 0.0;
          Particles[i].vel[p].reset();
          Particles[i].acc[p].reset();
          Particles[i].vrot[p] = 0.0;
          Particles[i].arot[p] = 0.0;
        }
      }
      if (Cell.mass < 1e-20) {
        Cell.mass = TotalMass;
      }
      std::cout << "With command 'grains' " << Particles.size() << " particles have been loaded" << std::endl;
    } else if (token == "Load") {
      // std::string command;
      conf >> command;
      if (command == "BiaxialCompression") {
        double pressure, velocity;
        conf >> pressure >> velocity;
        Load.BiaxialCompression(pressure, velocity);
      } else if (command == "VelocityBiaxialCompression") {
        double Vxx, Vyy;
        conf >> Vxx >> Vyy;
        Load.VelocityBiaxialCompression(Vxx, Vyy);
      } else if (command == "IsostaticCompression") {
        double pressure;
        conf >> pressure;
        Load.IsostaticCompression(pressure);
      } else if (command == "VelocityControl") {
        double vhxx, vhxy, vhyx, vhyy;
        conf >> vhxx >> vhxy >> vhyx >> vhyy;
        Load.VelocityControl(vhxx, vhxy, vhyx, vhyy);
      } else if (command == "SimpleShear") {
        double pressure, gammaDot;
        conf >> pressure >> gammaDot;
        Load.SimpleShear(pressure, gammaDot);
      } else if (command == "CyclicSimpleShear") {
        double pressure, gammaDot, T0;
        conf >> pressure >> gammaDot >> T0;
        Load.CyclicSimpleShear(pressure, gammaDot, T0); // 1/T0 Clycle loding frequency
      } else if (command == "ProportionalBC") {
        double dSigma, alpha;
        conf >> dSigma;
        conf >> alpha;
        Load.ProportionalBC(dSigma, alpha);
      } else if (command == "IncrementalCompression") {
        double dSigma;
        int Nincr /*, Ncount*/;
        conf >> dSigma;
        conf >> Nincr;
        Load.Nincr = Nincr;
        Load.IncrementalCompression(dSigma, Nincr);
      } else {
        std::cerr << "Unknown command for loading: " << command << std::endl;
      }
    } else if (token == "Stress") {
      conf >> Sig;
    } else if (token == "Platelets") {
      size_t nb;
      conf >> nb;
      Particles.clear();
      Nfibre = Ngrain = 0;
      Surf_p = 0.;
      for (size_t i = 0; i < nb; i++) {
        size_t nn;
        conf >> nn;
        Platelet P(nn);
        conf >> P.radius >> P.inertia >> P.mass;
        Surf_p += nn * M_PI * P.radius * P.radius;
        if (nn > 1) {
          Nfibre += 1;
          for (size_t p = 0; p < nn - 1; p++) {
            conf >> P.l0[p] >> P.dt[p] >> P.dtheta[p];
          }
        }
        if (nn < 2) {
          Ngrain += 1;
        }

        for (size_t p = 0; p < nn; p++) {
          conf >> P.pos[p] >> P.vel[p] >> P.acc[p] >> P.rot[p] >> P.vrot[p] >> P.arot[p];
        }
        Particles.push_back(P);
      }

    } else if (token == "Interactions") {
      Interactions.clear();
      size_t nb;
      conf >> nb;
      Interactions.clear();
      Interaction I;
      for (size_t i = 0; i < nb; i++) {
        conf >> I.i >> I.j >> I.ki >> I.kj >> I.fn >> I.ft >> I.delta_t;
        Interactions.push_back(I);
      }
    } else {
      std::cerr << "Unknown token: " << token << std::endl;
    }

    conf >> token;
  }
}

/*
void FiberCell2DSimulation::Limit_Liste() {
  vec2r C;
  PlateletXp.clear();
  PlateletXm.clear();
  PlateletYm.clear();
  PlateletYp.clear();
  PlateletAngGb.clear();
  PlateletAngGup.clear();
  PlateletAngDup.clear();
  PlateletAngDb.clear();
  LimitsG.clear();
  LimitsF.clear();
  for (size_t i = 0; i < Particles.size(); i++) {
    C.reset();

    for (size_t p = 0; p < Particles[i].n; p++) {
      C += Particles[i].pos[p];
    }
    C /= (double)(Particles[i].n);
    while (C.x < 0.0) {
      C.x += 1.;
      for (size_t p = 0; p < Particles[i].n; p++)
        Particles[i].pos[p].x += 1.0;
    }
    while (C.x > 1.0) {
      C.x -= 1.;
      for (size_t p = 0; p < Particles[i].n; p++)
        Particles[i].pos[p].x -= 1.0;
    }
    while (C.y < 0.0) {
      C.y += 1.;
      for (size_t p = 0; p < Particles[i].n; p++)
        Particles[i].pos[p].y += 1.0;
    }
    while (C.y > 1.0) {
      C.y -= 1.;
      for (size_t p = 0; p < Particles[i].n; p++)
        Particles[i].pos[p].y -= 1.0;
    }
    double Lmim = 0.2;
    double Lmax = 0.8;

    if (C.x < Lmim) {
      PlateletXm.push_back(i);
      Platelet P = Particles[i];
      if (ImageCell == 1) {
        for (size_t ip = 0; ip < Particles[i].n; ip++) {
          P.pos[ip].x = 1 + Particles[i].pos[ip].x;
        }
        if (Particles[i].n > 1)
          LimitsF.push_back(P);
        else if (Particles[i].n < 2)
          LimitsG.push_back(P);
      }
    }
    if (C.x > Lmax) {
      PlateletXp.push_back(i);
      if (ImageCell == 1) {
        if (Particles[i].n > 1)
          LimitsF.push_back(Particles[i]);
        else if (Particles[i].n < 2)
          LimitsG.push_back(Particles[i]);
      }
    }
    if (C.y < Lmim) {
      PlateletYm.push_back(i);

      Platelet P = Particles[i];
      if (ImageCell == 1) {
        for (size_t ip = 0; ip < Particles[i].n; ip++) {
          P.pos[ip].y = 1 + Particles[i].pos[ip].y;
        }
        if (Particles[i].n > 1)
          LimitsF.push_back(P);
        else if (Particles[i].n < 2)
          LimitsG.push_back(P);
      }
    }
    if (C.y > Lmax) {
      PlateletYp.push_back(i);
      if (ImageCell == 1) {
        if (Particles[i].n > 1)
          LimitsF.push_back(Particles[i]);
        else if (Particles[i].n < 2)
          LimitsG.push_back(Particles[i]);
      }
    }
    if (C.x < Lmim && C.y < Lmim) {
      PlateletAngGb.push_back(i);
      Platelet P = Particles[i];
      if (ImageCell == 1) {
        for (size_t ip = 0; ip < Particles[i].n; ip++) {
          P.pos[ip].x = 1 + Particles[i].pos[ip].x;
          P.pos[ip].y = 1 + Particles[i].pos[ip].y;
        }
        if (Particles[i].n > 1)
          LimitsF.push_back(P);
        else if (Particles[i].n < 2)
          LimitsG.push_back(P);
      }
    }
    if (C.x < Lmim && C.y > Lmax) {
      PlateletAngGup.push_back(i);
      if (ImageCell == 1) {
        if (Particles[i].n > 1)
          LimitsF.push_back(Particles[i]);
        else if (Particles[i].n < 2)
          LimitsG.push_back(Particles[i]);
      }
    }
    if (C.x > Lmax && C.y > Lmax) {
      PlateletAngDup.push_back(i);
      if (ImageCell == 1) {
        if (Particles[i].n > 1)
          LimitsF.push_back(Particles[i]);
        else if (Particles[i].n < 2)
          LimitsG.push_back(Particles[i]);
      }
    }
    if (C.x > Lmax && C.y < Lmim) {
      PlateletAngDb.push_back(i);
      Platelet P = Particles[i];
      if (ImageCell == 1) {
        for (size_t ip = 0; ip < Particles[i].n; ip++) {
          P.pos[ip].x = -1 + Particles[i].pos[ip].x;
          P.pos[ip].y = 1 + Particles[i].pos[ip].y;
        }
        if (Particles[i].n > 1)
          LimitsF.push_back(P);
        else if (Particles[i].n < 2)
          LimitsG.push_back(P);
      }
    }
  }
}
*/