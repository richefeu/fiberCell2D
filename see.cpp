#include "see.hpp"

void printHelp() {
  std::cout << std::endl;
  std::cout << "+         load next configuration file" << std::endl;
  std::cout << "-         load previous configuration file" << std::endl;
  std::cout << "=         fit the view" << std::endl;
  std::cout << "a/A       particle transparency" << std::endl;
  std::cout << "b/B       ghost particle transparency" << std::endl;
  std::cout << "c         show/hide periodic cell" << std::endl;
  std::cout << "C         show/hide contacts" << std::endl;
  std::cout << "f         show/hide normal forces" << std::endl;
  std::cout << "g         show/hide ghost particles" << std::endl;
  std::cout << "h         print this help" << std::endl;
  std::cout << "i         print information" << std::endl;
  std::cout << "n         go to file (see terminal to enter the file number)" << std::endl;
  std::cout << "o         show/hide particle orientations" << std::endl;
  std::cout << "p         show/hide particles" << std::endl;
  std::cout << "P         switch pipe display" << std::endl;
  std::cout << "q         quit" << std::endl;
  std::cout << "s/S       tune vector sizes" << std::endl;
  std::cout << "w/W       tune displayed ghost width" << std::endl;
  // std::cout << "x         xxxx" << std::endl;
  std::cout << std::endl;
  std::cout << "0         particles colored with light gray" << std::endl;
  std::cout << "1         particles colored with pressure" << std::endl;
  std::cout << std::endl;
}

void printInfo() {
  int V = glutGet(GLUT_VERSION);
  int V1 = V / 10000;
  int V2 = (V - V1 * 10000) / 100;
  int V3 = V - V1 * 10000 - V2 * 100;
  std::cout << "glut version " << V1 << "." << V2 << "." << V3 << "\n";
}

void keyboard(unsigned char Key, int /*x*/, int /*y*/) {
  switch (Key) {

  case '0': {
    color_option = 0;
  } break;

  case '1': { // particle pressures
    color_option = 1;
  } break;

  case '2': {
    // colorTable.setMinMax(pmin, pmax);
    colorTable.setTableID(2);
    colorTable.Rebuild();
    color_option = 2;
  } break;

  case 'a': {
    alpha_particles = std::max(0.0f, alpha_particles - 0.05f);
  } break;
  case 'A': {
    alpha_particles = std::min(1.0f, alpha_particles + 0.05f);
  } break;

  case 'b': {
    alpha_ghosts = std::max(0.0f, alpha_ghosts - 0.05f);
  } break;
  case 'B': {
    alpha_ghosts = std::min(1.0f, alpha_ghosts + 0.05f);
  } break;

  case 'c': {
    show_cell = 1 - show_cell;
  } break;

  case 'C': {
    show_contacts = 1 - show_contacts;
  } break;

  case 'f': {
    show_forces = 1 - show_forces;
  } break;

  case 'g': {
    show_ghosts = 1 - show_ghosts;
  } break;

  case 'h': {
    printHelp();
  } break;

  case 'i': {
    printInfo();
  } break;

  case 'n': {
    std::cout << "Go to file number ";
    int conNumTry;
    std::cin >> conNumTry;
    try_to_readConf(conNumTry, Conf, confNum);
  } break;

  case 'o': {
    showOrientations = 1 - showOrientations;
  } break;

  case 'p': {
    show_particles = 1 - show_particles;
  } break;

    // case 'P': {
    //  show_pipe = 1 - show_pipe;
    //} break;

  case 'q': {
    exit(0);
  } break;

  case 'S': {
    vScale *= 1.05;
  } break;

  case 's': {
    vScale *= 0.95;
    if (vScale < 0.0)
      vScale = 1.0;
  } break;

  case 'w': {
    ghost_width = std::max(0.0, ghost_width - 0.05);
  } break;
  case 'W': {
    ghost_width = std::min(0.5, ghost_width + 0.05);
  } break;

  case '-': {
    if (confNum > 0) {
      try_to_readConf(confNum - 1, Conf, confNum);
    }
    updateTextLine();
  } break;

  case '+': {
    try_to_readConf(confNum + 1, Conf, confNum);
    updateTextLine();
  } break;

  case '=': {
    fit_view();
    reshape(width, height);
  } break;
  };

  glutPostRedisplay();
}

void updateTextLine() { textZone.addLine("conf %d,  t %0.4g s", confNum, Conf.t); }

void mouse(int button, int state, int x, int y) {

  if (state == GLUT_UP) {
    mouse_mode = NOTHING;
    display();
  } else if (state == GLUT_DOWN) {
    mouse_start[0] = x;
    mouse_start[1] = y;
    switch (button) {
    case GLUT_LEFT_BUTTON:
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
        mouse_mode = PAN;
      else
        mouse_mode = ROTATION;
      break;
    case GLUT_MIDDLE_BUTTON:
      mouse_mode = ZOOM;
      break;
    }
  }
}

void motion(int x, int y) {

  if (mouse_mode == NOTHING)
    return;

  double dx = (double)(x - mouse_start[0]) / (double)width;
  double dy = (double)(y - mouse_start[1]) / (double)height;

  switch (mouse_mode) {

  case ZOOM: {
    double ddy = (worldBox.max.y - worldBox.min.y) * dy;
    double ddx = (worldBox.max.x - worldBox.min.x) * dy;
    worldBox.min.x -= ddx;
    worldBox.max.x += ddx;
    worldBox.min.y -= ddy;
    worldBox.max.y += ddy;
  } break;

  case PAN: {
    double ddx = (worldBox.max.x - worldBox.min.x) * dx;
    double ddy = (worldBox.max.y - worldBox.min.y) * dy;
    worldBox.min.x -= ddx;
    worldBox.max.x -= ddx;
    worldBox.min.y += ddy;
    worldBox.max.y += ddy;
  } break;

  default:
    break;
  }
  mouse_start[0] = x;
  mouse_start[1] = y;

  reshape(width, height);
  display();
}

void display() {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);

  if (show_cell) {
    drawBox();
  }

  if (show_particles == 1) {
    drawParticles();
  }

  if (show_contacts == 1) {
    drawContacts();
  }

  if (show_forces == 1) {
    drawForces();
  }

  if (show_ghosts == 1) {
    drawGhosts();
  }

  textZone.draw();

  glFlush();
  glutSwapBuffers();
}

void fit_view() {
  //
  // 3 x ------- x 2
  //   |         |
  // 0 x ------- x 1
  double x0 = 0.0;
  double y0 = 0.0;
  double x1 = Conf.Cell.h.xy;
  double y1 = Conf.Cell.h.yy;
  double x2 = Conf.Cell.h.xy + Conf.Cell.h.xx;
  double y2 = Conf.Cell.h.yy + Conf.Cell.h.yx;
  double x3 = Conf.Cell.h.xx;
  double y3 = Conf.Cell.h.yx;
  worldBox.min.x = std::min(std::min(std::min(x0, x1), x2), x3);
  worldBox.min.y = std::min(std::min(std::min(y0, y1), y2), y3);
  worldBox.max.x = std::max(std::max(std::max(x0, x1), x2), x3);
  worldBox.max.y = std::max(std::max(std::max(y0, y1), y2), y3);
  reshape(width, height);
}

void reshape(int w, int h) {
  width = w;
  height = h;

  double left = worldBox.min.x;
  double right = worldBox.max.x;
  double bottom = worldBox.min.y;
  double top = worldBox.max.y;
  double worldW = right - left;
  double worldH = top - bottom;
  double dW = 0.1 * worldW;
  double dH = 0.1 * worldH;
  left -= dW;
  right += dW;
  top += dH;
  bottom -= dH;
  worldW = right - left;
  worldH = top - bottom;

  if (worldW > worldH) {
    worldH = worldW * ((GLfloat)height / (GLfloat)width);
    top = 0.5 * (bottom + top + worldH);
    bottom = top - worldH;
  } else {
    worldW = worldH * ((GLfloat)width / (GLfloat)height);
    right = 0.5 * (left + right + worldW);
    left = right - worldW;
  }

  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(left, right, bottom, top);

  glutPostRedisplay();
}

// Draw the periodic cell
void drawBox() {
  glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
  glLineWidth(1.0f);

  glBegin(GL_LINE_LOOP);
  glVertex2f(0.0f, 0.0f);
  glVertex2f(Conf.Cell.h.xy, Conf.Cell.h.yy);
  glVertex2f(Conf.Cell.h.xy + Conf.Cell.h.xx, Conf.Cell.h.yy + Conf.Cell.h.yx);
  glVertex2f(Conf.Cell.h.xx, Conf.Cell.h.yx);
  glEnd();
}

void setColor(int i, GLfloat alpha) {
  switch (color_option) {

  case 0: {
    glColor4f(0.8f, 0.8f, 0.9f, alpha);
  } break;

  case 1: { // pressure
    // colorRGBA col;
    // colorTable.getRGB(pdata[i].p, &col);
    // glColor4f(col.r / 255.0, col.g / 255.0, col.b / 255.0, 1.0f);
  } break;

  case 2: {
    /*
    colorRGBA col;
    colorTable.getRGB(Conf.grains[i].zncc2, &col);
    glColor4f(col.r / 255.0, col.g / 255.0, col.b / 255.0, 1.0f);
    */
  } break;

  default: {
    glColor4f(0.8f, 0.8f, 0.9f, alpha);
  } break;
  }
}

void add_ghost_pos(int i, double mn, double mx, std::vector<vec2r> &lst) {
  /*
  lst.clear();
  vec2r pos = Conf.Particles[i].pos;
  if (pos.x > mn && pos.x < mx && pos.y > mn && pos.y < mx) {
    return;
  }

  vec2r ghostPos;

  if (pos.x <= mn) {
    ghostPos.set(pos.x + 1.0, pos.y);
    lst.push_back(ghostPos);

    if (pos.y <= mn) {
      ghostPos.set(pos.x + 1.0, pos.y + 1.0);
      lst.push_back(ghostPos);
    }
    if (pos.y >= mx) {
      ghostPos.set(pos.x + 1.0, pos.y - 1.0);
      lst.push_back(ghostPos);
    }
  }

  if (pos.x >= mx) {
    ghostPos.set(pos.x - 1.0, pos.y);
    lst.push_back(ghostPos);

    if (pos.y <= mn) {
      ghostPos.set(pos.x - 1.0, pos.y + 1.0);
      lst.push_back(ghostPos);
    }
    if (pos.y >= mx) {
      ghostPos.set(pos.x - 1.0, pos.y - 1.0);
      lst.push_back(ghostPos);
    }
  }

  if (pos.y <= mn) {
    ghostPos.set(pos.x, pos.y + 1.0);
    lst.push_back(ghostPos);

    if (pos.x <= mn) {
      ghostPos.set(pos.x + 1.0, pos.y + 1.0);
      lst.push_back(ghostPos);
    }
    if (pos.x >= mx) {
      ghostPos.set(pos.x - 1.0, pos.y + 1.0);
      lst.push_back(ghostPos);
    }
  }

  if (pos.y >= mx) {
    ghostPos.set(pos.x, pos.y - 1.0);
    lst.push_back(ghostPos);

    if (pos.x <= mn) {
      ghostPos.set(pos.x + 1.0, pos.y - 1.0);
      lst.push_back(ghostPos);
    }
    if (pos.x >= mx) {
      ghostPos.set(pos.x - 1.0, pos.y - 1.0);
      lst.push_back(ghostPos);
    }
  }
  */
}

void loop2triangles(std::vector<vec2r> &coords, std::vector<int> &indx) {
  indx.clear();
  if (coords.size() < 3) {
    return; 
  }

  size_t iend = coords.size() - 1;
  size_t imid = 0;
  size_t ibeg = 1;

  indx.push_back(ibeg);
  indx.push_back(imid);
  indx.push_back(iend);

  int toggle = 0;
  while (iend != imid && ibeg != imid) {
    if (toggle == 0) {
      imid = ibeg;
      ibeg++;
      if (ibeg != imid && imid != iend) {
        indx.push_back(ibeg);
        indx.push_back(imid);
        indx.push_back(iend);
      }
    } else {
      imid = iend;
      iend--;
      if (ibeg != imid && imid != iend) {
        indx.push_back(ibeg);
        indx.push_back(imid);
        indx.push_back(iend);
      }
    }
    toggle = 1 - toggle;
  }
}

void drawParticles() {
  glLineWidth(1.0f);

  for (size_t i = 0; i < Conf.Particles.size(); ++i) {

    /// on changera ça...
    //mat4r H;
    //H.xx = Conf.Cell.hxx;
    //H.xy = Conf.Cell.hxy;
    //H.yx = Conf.Cell.hyx;
    //H.yy = Conf.Cell.hyy;

    vec2r pos = Conf.Cell.h * Conf.Particles[i].pos[0];

    double R = Conf.Particles[i].radius;

    setColor(i, alpha_particles);

    std::vector<vec2r> coords;
    Conf.Particles[i].getShape(Conf.Cell.h, coords, 0.05);

    if (Conf.Particles[i].n < 2) { // FIXME: pour le moment on ne colorie pas les fibres
      glBegin(GL_POLYGON);
      for (size_t v = 0; v < coords.size(); v++) {
        glVertex2f(coords[v].x, coords[v].y);
      }
      glEnd();
      
    }
    
    glColor4f(0.0f, 0.0f, 0.0f, alpha_particles);
    glBegin(GL_LINE_LOOP);
    for (size_t v = 0; v < coords.size(); v++) {
      glVertex2f(coords[v].x, coords[v].y);
    }
    glEnd();

    if (showOrientations) {
      if (Conf.Particles[i].n < 2) {
        double rot = Conf.Particles[i].rot[0];
        glBegin(GL_LINES);
        glVertex2f(pos.x, pos.y);
        glVertex2f(pos.x + R * cos(rot), pos.y + R * sin(rot));
        glEnd();
      } else {
        vec2r center = Conf.Cell.h * Conf.Particles[i].pos[0];
        glBegin(GL_LINE_STRIP);
        glVertex2f(center.x, center.y);
        for (size_t c = 0; c < Conf.Particles[i].pos.size(); c++) {
          center = Conf.Cell.h * Conf.Particles[i].pos[c];
          glVertex2f(center.x, center.y);
        }
        glEnd();
      }
    }
  }
}

void drawGhosts() {
  /*
  // if (mouse_mode != NOTHING && box.Particles.size() > 2000) return;

  std::vector<vec2r> lst_pos; // list of reduced positions of ghost particles
  double mn = ghost_width;
  double mx = 1.0 - ghost_width;
  // GLColorRGBA color;
  glLineWidth(1.0f);
  for (size_t i = 0; i < Conf.Particles.size(); ++i) {
    add_ghost_pos(i, mn, mx, lst_pos);
    double R = Conf.Particles[i].radius;
    for (size_t ig = 0; ig < lst_pos.size(); ig++) {

      vec2r pos = Conf.Cell.h * lst_pos[ig];

      setColor(i, alpha_ghosts);
      glBegin(GL_POLYGON);
      for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
        glVertex2f(pos.x + R * cos(angle), pos.y + R * sin(angle));
      }
      glEnd();

      glColor4f(0.0f, 0.0f, 0.0f, alpha_particles);
      glBegin(GL_LINE_LOOP);
      for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
        glVertex2f(pos.x + R * cos(angle), pos.y + R * sin(angle));
      }
      glEnd();

      if (showOrientations) {
        double rot = Conf.Particles[i].rot;
        glBegin(GL_LINES);
        glVertex2f(pos.x, pos.y);
        glVertex2f(pos.x + R * cos(rot), pos.y + R * sin(rot));
        glEnd();
      }
    }
  }
  */
}

void drawContacts() {
  /*
  glLineWidth(1.5f);

  // grain-grain
  glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
  glBegin(GL_LINES);
  for (size_t k = 0; k < Conf.Interactions.size(); ++k) {
    size_t i = Conf.Interactions[k].i;
    size_t j = Conf.Interactions[k].j;
    vec2r posi = Conf.Cell.h * Conf.Particles[i].pos;
    vec2r sij = Conf.Particles[j].pos - Conf.Particles[i].pos;
    sij.x -= floor(sij.x + 0.5);
    sij.y -= floor(sij.y + 0.5);
    vec2r posj = posi + Conf.Cell.h * sij;
    glVertex2f(posi.x, posi.y);
    glVertex2f(posj.x, posj.y);
  }
  glEnd();

  // with pipe
  glColor4f(0.0f, 1.0f, 0.0f, 1.0f);
  glBegin(GL_LINES);
  for (size_t k = 0; k < Conf.InteractionsPipe.size(); ++k) {
    size_t i = Conf.InteractionsPipe[k].i;
    size_t inode = Conf.InteractionsPipe[k].inode;
    vec2r posi = Conf.Cell.h * Conf.Particles[i].pos;
    vec2r sij = Conf.pipe.pos[inode] - Conf.Particles[i].pos;
    double proj = -sij * Conf.pipe.u[inode];
    if (proj > 0.0) {
      vec2r n = Conf.pipe.u[inode].quarterLeftTurned();
      sij = (Conf.Particles[i].radius + Conf.pipe.node_radius) * n;
    } else {
      sij = Conf.Cell.h * sij;
    }
    vec2r posj = posi + sij;
    glVertex2f(posi.x, posi.y);
    glVertex2f(posj.x, posj.y);
  }
  glEnd();
  */
}

void drawForces() {
  /*
  // grain-grain
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);

  for (size_t k = 0; k < Conf.Interactions.size(); ++k) {
    size_t i = Conf.Interactions[k].i;
    size_t j = Conf.Interactions[k].j;
    vec2r posi = Conf.Cell.h * Conf.Particles[i].pos;
    vec2r sij = Conf.Particles[j].pos - Conf.Particles[i].pos;
    sij.x -= floor(sij.x + 0.5);
    sij.y -= floor(sij.y + 0.5);
    vec2r posj = posi + Conf.Cell.h * sij;

    // Calculate the width of the rectangle
    GLfloat width = Conf.radiusMean * (Conf.Interactions[k].fn / Conf.fnMax);

    // Calculate the direction vector and the perpendicular vector
    vec2r dir = posj - posi;
    vec2r perp(-dir.y, dir.x);
    perp.normalize();
    perp *= 0.5 * width;

    // Calculate the four corners of the rectangle
    vec2r p1 = posi + perp;
    vec2r p2 = posi - perp;
    vec2r p3 = posj - perp;
    vec2r p4 = posj + perp;

    // Draw the filled rectangle
    glBegin(GL_QUADS);
    glVertex2f(p1.x, p1.y);
    glVertex2f(p2.x, p2.y);
    glVertex2f(p3.x, p3.y);
    glVertex2f(p4.x, p4.y);
    glEnd();
  }

  // with pipe
  // glColor4f(0.0f, 1.0f, 0.0f, 1.0f);
  glBegin(GL_LINES);
  for (size_t k = 0; k < Conf.InteractionsPipe.size(); ++k) {
    size_t i = Conf.InteractionsPipe[k].i;
    size_t inode = Conf.InteractionsPipe[k].inode;
    vec2r posi = Conf.Cell.h * Conf.Particles[i].pos;
    vec2r sij = Conf.pipe.pos[inode] - Conf.Particles[i].pos;
    double proj = -sij * Conf.pipe.u[inode];
    if (proj > 0.0) {
      vec2r n = Conf.pipe.u[inode].quarterLeftTurned();
      sij = (Conf.Particles[i].radius + Conf.pipe.node_radius) * n;
    } else {
      sij = Conf.Cell.h * sij;
    }
    vec2r posj = posi + sij;

    // Calculate the width of the rectangle
    GLfloat width = Conf.radiusMean * (Conf.InteractionsPipe[k].fn / Conf.fnMax);

    // Calculate the direction vector and the perpendicular vector
    vec2r dir = posj - posi;
    vec2r perp(-dir.y, dir.x);
    perp.normalize();
    perp *= 0.5 * width;

    // Calculate the four corners of the rectangle
    vec2r p1 = posi + perp;
    vec2r p2 = posi - perp;
    vec2r p3 = posj - perp;
    vec2r p4 = posj + perp;

    // Draw the filled rectangle
    glBegin(GL_QUADS);
    glVertex2f(p1.x, p1.y);
    glVertex2f(p2.x, p2.y);
    glVertex2f(p3.x, p3.y);
    glVertex2f(p4.x, p4.y);
    glEnd();
  }
  glEnd();
  */
}

bool try_to_readConf(int num, FiberCell2DSimulation &CF, int &OKNum) {
  char file_name[256];
  snprintf(file_name, 256, "conf%d", num);
  if (fileTool::fileExists(file_name)) {
    std::cout << file_name << std::endl;
    OKNum = num;
    CF.loadConf(file_name);
    // CF.pipe.updateShape(CF.Cell.h);
    // CF.accelerations();
    // computeParticleData();
    return true;
  } else {
    std::cout << file_name << " does not exist" << std::endl;
  }
  return false;
}

void menu(int num) {
  switch (num) {

  case 0: {
    exit(0);
  } break;

    /*
    case 10: { // None
      show_pipe = 0;
    } break;
    case 11: { // Pipe alone
      show_pipe = 1;
      show_pipe_nodes = 0;
    } break;
    case 12: { // Show nodes
      show_pipe = 1;
      show_pipe_nodes = 1;
    } break;
    case 13: { // Show loading
      show_pipe = 1;
      show_pipe_nodes = 0;
    } break;
    case 14: { // Show internal stress
      show_pipe = 1;
      show_pipe_nodes = 0;
    } break;
  */
    // case xx: {} break;
  };

  glutPostRedisplay();
}

void buildMenu() {
  int submenu1 = glutCreateMenu(menu);
  glutAddMenuEntry("None", 10);
  glutAddMenuEntry("Pipe alone", 11);
  glutAddMenuEntry("Show nodes", 12);
  glutAddMenuEntry("Show loading", 13);
  glutAddMenuEntry("Show internal stress", 14);

  int submenu2 = glutCreateMenu(menu); // Force Colors
  glutAddMenuEntry("None", 200);
  glutAddMenuEntry("Magnitude", 201);

  glutCreateMenu(menu); // Main menu
  glutAddSubMenu("Pipe Display", submenu1);
  glutAddSubMenu("Force Colors", submenu2);
  glutAddMenuEntry("Quit", 0);
}

// =====================================================================
// Main function
// =====================================================================

int main(int argc, char *argv[]) {

  if (argc == 1) {
    confNum = 0;
  } else if (argc == 2) {
    confNum = atoi(argv[1]);
  }

  std::cout << "Current Configuration: ";
  try_to_readConf(confNum, Conf, confNum);

  mouse_mode = NOTHING;

  glEnable(GL_BLEND);
  glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // ==== Init GLUT and create window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA);
  int X0 = (glutGet(GLUT_SCREEN_WIDTH) - width) / 2;
  int Y0 = (glutGet(GLUT_SCREEN_HEIGHT) - height) / 2;
  glutInitWindowPosition(X0, Y0);
  glutInitWindowSize(width, height);

  main_window = glutCreateWindow("CONF VISUALIZER");

  // ==== Register callbacks
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);

  // ==== Menu
  buildMenu();
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  // ==== Other initialisations
  glText::init();
  updateTextLine();

  // ==== Enter GLUT event processing cycle
  fit_view();
  glutMainLoop();
  return 0;
}
