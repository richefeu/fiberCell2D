using namespace std;
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

// Input parameters
float R_fg; // = surface occupied by fibers / (surface occupied by fibers and grains) = surface ratio of fibers
int Nd_f; // number of disks per fiber
float eta; // overlap between the discretisation disks of the fibers
int Lb; // the box size is Lb times the length of a fiber
int echantillon_i; // index of the generated sample

// Variables computed from the input parameters
int ng; // number of grains
int nf; // number of fibers
double L_max; // full length of the fibers
double L; // length of the fibers (between the centers of the end grains)
const double R_f=1e-3; // radius of fibers is fixed to 1
const float C_fg = 0.50; // target compacity is fixed to 0.5

// miscellaneous variables
const double PI = acos(-1.0);
const int nth = 10; // for the representation of the fibers in paraview, nb of points -1 per half-circle (end of the fiber)
const double dth = PI/nth; // delta_theta for the representation of fibers in paraview
const float b_xmin = 0.0; // x coordinate of the lower side of the ximulation box
float b_xmax; // x coordinate of the upper side of the simulation box
const float b_ymin = 0.0; // y coordinate of the left side of the ximulation box
float b_ymax; // y coordinate of the right side of the simulation box

struct grain // definition of a grain
{
    double R; // radius
    double x,y; // coordinates of the center
};
std::vector <grain> vgrains; // vector containing the grains
//change para mostrar resultado
std::ostream& operator<<(std::ostream& os, const grain& g) {
    os << "{ Radius: " << g.R << ", x: " << g.x << ", y: " << g.y << " }";
    return os;
}

struct fibre // definition of a fiber
{
  std::vector <double> pos_c_1;
  std::vector <double> pos_c_2;
  double R; // radius
  int Nd; // number of discretisation disks
  double alpha; // orientation of the fiber
};
std::vector <fibre> vfibres; // vector containing the fibers

struct p_fibres // point of a fiber for representation in paraview
{
  std::vector <double> pos_1;
  std::vector <double> pos_2;
};
std::vector <p_fibres> vp_f; // vector containing the coordinates of the points for drawing the fibers in paraview
p_fibres p_f;

// Function declaration

// Creating output files
void output_grains(const string& nomFichier, const string& cheminDestination); // save grains to a vtk file
void output_fibres(const string& nomFichier, const string& cheminDestination); // save fibers to a vtk file
void output_fibres_up(int numfile); // save fibers to a vtk file
void output_fibres_down(int numfile); // save fibers to a vtk file
void output_fibres_left(int numfile); // save fibers to a vtk file
void output_fibres_right(int numfile); // save fibers to a vtk file
void output_fibres_down_left(int numfile); // save fibers to a vtk file
void output_fibres_down_right(int numfile); // save fibers to a vtk file
void output_fibres_up_left(int numfile); // save fibers to a vtk file
void output_fibres_up_right(int numfile); // save fibers to a vtk file
void output_fibres_grains(const string& nomFichier, const string& cheminDestination); // creates the input file for the DEM simulation

// Checking for collisions
bool doIntersect(double, double, double, double, double, double, double, double, float, float);
bool collisions_fibres_grains(int, double, double, double, double);

// random generator
std::random_device rd;
double seed = rd();
std::mt19937 gen(seed); 

// Main program
int main(int argc, char* argv[])
{
  int i;
  if (argc < 6)
    {
      std::cerr << "Usage: " << argv[0] << "<R_fg> <Nd_f> <eta> <Lb> <sample_index>" << std::endl;
      return 1;
    }

  R_fg = std::stof(argv[1]);
  Nd_f = std::stoi(argv[2]);
  eta = std::stof(argv[3]);
  Lb = std::stoi(argv[4]);
  echantillon_i = std::stoi(argv[5]);

  // length of the fibers
  L_max = 2*(1-eta)*Nd_f*R_f;
  L = 2*(1-eta)*(Nd_f-1)*R_f;

  b_xmax = Lb*L_max;
  b_ymax=b_xmax;
  cout << "b_xmax : " << b_xmax << " ; b_ymax : " << b_ymax << endl;

  float S_b = b_xmax*b_ymax; // surface of simulation box
  double S_1f = L*2*R_f + (PI*R_f*R_f); // surface of a fiber
  double inv_S_1f = 1/S_1f;
  //double S_1g = PI*1.2*R_f*1.2*R_f; // average surface of a grain (if random uniform distribution of the radii)
  double avR = (sqrt(2)*R_f+R_f)/2;
  double S_1g = PI*avR*avR;
  cout << "av : " << avR << endl;
  double inv_S_1g = 1/S_1g;

  // compute the number of fibers and grains requested in the sample
  if (R_fg == 0) //if there is no fibers, the compacity is going to be filled just by grains.
    {
      nf = 0;
      ng = C_fg*S_b*inv_S_1g;
    }
  else
    {
      double inv_R_fg = 1/R_fg;
      nf = C_fg*R_fg*S_b*inv_S_1f;
      ng = nf*S_1f*inv_S_1g*(1-R_fg)*inv_R_fg;
    }
  cout << "ng : " << ng << " ; nf : " << nf << endl;

  /////////////////////
  // FIBRE INSERTION //
  /////////////////////

  // uniform random distribution of fiber positions and orientations
  std::uniform_real_distribution<float> distribution_pos_f_x(b_xmin, b_xmax);
  std::uniform_real_distribution<float> distribution_pos_f_y(b_ymin, b_ymax);
  std::uniform_real_distribution<float> distribution_alpha(0, 2*PI);

  int att_f = 0; // nb of fiber insertion attempts
  while (vfibres.size()<nf) // while the nf fibers have not yet been positioned in the box
    {
      att_f++;
      int compteur_f = 0; // counting the collisions

      // generate the position and orientation of the candidate fiber
      float x_f = distribution_pos_f_x(gen);
      float y_f = distribution_pos_f_y(gen);
      float alpha = distribution_alpha(gen);

      // coordinates of the centers of the ends of the fiber
      double x_c_1 = x_f - L/2*cos(alpha);
      double y_c_1 = y_f - L/2*sin(alpha);
      double x_c_2 = x_f + L/2*cos(alpha);
      double y_c_2 = y_f + L/2*sin(alpha);

      fibre f = {{x_c_1, y_c_1},{x_c_2,y_c_2}, R_f, Nd_f, alpha}; // the candidate fiber

      //////////////////////////////////////////
      // Testing if collisions between fibers //
      //////////////////////////////////////////

      for (i=0; i<vfibres.size(); i++) // loop on fibers already accepted in the box
      {
        // 1. testing if fibers cross each other along their main segments (not involving the ends) //

        // main box
        if (doIntersect(f.pos_c_1[0], f.pos_c_1[1], f.pos_c_2[0], f.pos_c_2[1],
                        vfibres[i].pos_c_1[0], vfibres[i].pos_c_1[1], vfibres[i].pos_c_2[0], vfibres[i].pos_c_2[1],
                        f.alpha, vfibres[i].alpha)){
	          //crossing detected
            compteur_f++;
            break;
          }

        // lower periodic box
        if (doIntersect(f.pos_c_1[0], f.pos_c_1[1]-b_ymax, f.pos_c_2[0], f.pos_c_2[1]-b_ymax,
                        vfibres[i].pos_c_1[0], vfibres[i].pos_c_1[1], vfibres[i].pos_c_2[0], vfibres[i].pos_c_2[1],
                        f.alpha, vfibres[i].alpha)){
            compteur_f++;
            break;
	         }

        // upper periodic box
        if (doIntersect(f.pos_c_1[0], f.pos_c_1[1]+b_ymax, f.pos_c_2[0], f.pos_c_2[1]+b_ymax,
                        vfibres[i].pos_c_1[0], vfibres[i].pos_c_1[1], vfibres[i].pos_c_2[0], vfibres[i].pos_c_2[1],
                        f.alpha, vfibres[i].alpha)){
            compteur_f++;
            break;
	         }

	       // left periodic box
         if (doIntersect(f.pos_c_1[0]-b_xmax, f.pos_c_1[1], f.pos_c_2[0]-b_xmax, f.pos_c_2[1],
                         vfibres[i].pos_c_1[0], vfibres[i].pos_c_1[1], vfibres[i].pos_c_2[0], vfibres[i].pos_c_2[1],
                         f.alpha, vfibres[i].alpha)){
            compteur_f++;
            break;
           }

	       // right periodic box
	       if (doIntersect(f.pos_c_1[0]+b_xmax, f.pos_c_1[1], f.pos_c_2[0]+b_xmax, f.pos_c_2[1],
                         vfibres[i].pos_c_1[0], vfibres[i].pos_c_1[1], vfibres[i].pos_c_2[0], vfibres[i].pos_c_2[1],
                         f.alpha, vfibres[i].alpha)){
	           compteur_f++;
             break;
          }

	       // top right periodic corner
         if (doIntersect(f.pos_c_1[0]+b_xmax, f.pos_c_1[1]+b_ymax, f.pos_c_2[0]+b_xmax, f.pos_c_2[1]+b_ymax,
                         vfibres[i].pos_c_1[0], vfibres[i].pos_c_1[1], vfibres[i].pos_c_2[0], vfibres[i].pos_c_2[1],
                         f.alpha, vfibres[i].alpha)){
    	      compteur_f++;
	          break;
	       }

	       // top left periodic corner
	       if (doIntersect(f.pos_c_1[0]-b_xmax, f.pos_c_1[1]+b_ymax, f.pos_c_2[0]-b_xmax, f.pos_c_2[1]+b_ymax,
                         vfibres[i].pos_c_1[0], vfibres[i].pos_c_1[1], vfibres[i].pos_c_2[0], vfibres[i].pos_c_2[1],
                         f.alpha, vfibres[i].alpha)){
    	      compteur_f++;
	          break;
	       }


	       // bottom right periodic corner
	       if (doIntersect(f.pos_c_1[0]+b_xmax, f.pos_c_1[1]-b_ymax, f.pos_c_2[0]+b_xmax, f.pos_c_2[1]-b_ymax,
                         vfibres[i].pos_c_1[0], vfibres[i].pos_c_1[1], vfibres[i].pos_c_2[0], vfibres[i].pos_c_2[1],
                         f.alpha, vfibres[i].alpha)){
    	      compteur_f++;
	          break;
	       }

	       // bottom left periodic corner
	       if (doIntersect(f.pos_c_1[0]-b_xmax, f.pos_c_1[1]-b_ymax, f.pos_c_2[0]-b_xmax, f.pos_c_2[1]-b_ymax,
                         vfibres[i].pos_c_1[0], vfibres[i].pos_c_1[1], vfibres[i].pos_c_2[0], vfibres[i].pos_c_2[1],
                         f.alpha, vfibres[i].alpha)){
    	      compteur_f++;
	          break;
	       }

	        // 3. testing if the ends of two fibers collide

	        // main box
          // compute the distance between end 1 of the candidate fiber and end 1 of fiber i
          float xji_f1f1 = f.pos_c_1[0]-vfibres[i].pos_c_1[0];
          float yji_f1f1 = f.pos_c_1[1]-vfibres[i].pos_c_1[1];
          float CjCi_f1f1 = sqrt(xji_f1f1*xji_f1f1+yji_f1f1*yji_f1f1);
          // compute the distance between end 1 of the candidate fiber and end 2 of fiber i
          float xji_f1f2 = f.pos_c_1[0]-vfibres[i].pos_c_2[0];
          float yji_f1f2 = f.pos_c_1[1]-vfibres[i].pos_c_2[1];
          float CjCi_f1f2 = sqrt(xji_f1f2*xji_f1f2+yji_f1f2*yji_f1f2);
          // compute the distance between end 2 of the candidate fiber and end 1 of fiber i
          float xji_f2f1 = f.pos_c_2[0]-vfibres[i].pos_c_1[0];
          float yji_f2f1 = f.pos_c_2[1]-vfibres[i].pos_c_1[1];
          float CjCi_f2f1 = sqrt(xji_f2f1*xji_f2f1+yji_f2f1*yji_f2f1);
          // compute the distance between end 2 of the candidate fiber and end 2 of fiber i
          float xji_f2f2 = f.pos_c_2[0]-vfibres[i].pos_c_2[0];
          float yji_f2f2 = f.pos_c_2[1]-vfibres[i].pos_c_2[1];
          float CjCi_f2f2 = sqrt(xji_f2f2*xji_f2f2+yji_f2f2*yji_f2f2);

          // compute the corresponding contact distances
          float dn_f1f1 = CjCi_f1f1 - (f.R+vfibres[i].R);
          float dn_f1f2 = CjCi_f1f2 - (f.R+vfibres[i].R);
          float dn_f2f1 = CjCi_f2f1 - (f.R+vfibres[i].R);
          float dn_f2f2 = CjCi_f2f2 - (f.R+vfibres[i].R);

          if(dn_f1f1<0.0 || dn_f1f2<0.0 || dn_f2f1<0.0 || dn_f2f2<0.0){
             // overlap detected
             compteur_f++;
	           break;
	        }

	         // lower periodic box
           xji_f1f1 = f.pos_c_1[0]-vfibres[i].pos_c_1[0];
           yji_f1f1 = (f.pos_c_1[1]-b_ymax)-vfibres[i].pos_c_1[1];
           xji_f1f2 = f.pos_c_1[0]-vfibres[i].pos_c_2[0];
           yji_f1f2 = (f.pos_c_1[1]-b_ymax)-vfibres[i].pos_c_2[1];
           xji_f2f1 = f.pos_c_2[0]-vfibres[i].pos_c_1[0];
           yji_f2f1 = (f.pos_c_2[1]-b_ymax)-vfibres[i].pos_c_1[1];
           xji_f2f2 = f.pos_c_2[0]-vfibres[i].pos_c_2[0];
           yji_f2f2 = (f.pos_c_2[1]-b_ymax)-vfibres[i].pos_c_2[1];

           CjCi_f1f1 = sqrt(xji_f1f1*xji_f1f1+yji_f1f1*yji_f1f1);
           CjCi_f1f2 = sqrt(xji_f1f2*xji_f1f2+yji_f1f2*yji_f1f2);
           CjCi_f2f1 = sqrt(xji_f2f1*xji_f2f1+yji_f2f1*yji_f2f1);
           CjCi_f2f2 = sqrt(xji_f2f2*xji_f2f2+yji_f2f2*yji_f2f2);

           dn_f1f1 = CjCi_f1f1 - (f.R+vfibres[i].R);
           dn_f1f2 = CjCi_f1f2 - (f.R+vfibres[i].R);
           dn_f2f1 = CjCi_f2f1 - (f.R+vfibres[i].R);
           dn_f2f2 = CjCi_f2f2 - (f.R+vfibres[i].R);

           if(dn_f1f1<0.0 || dn_f1f2<0.0 || dn_f2f1<0.0 || dn_f2f2<0.0) {
	            compteur_f++;
	            break;
	         }

	         // upper periodic box
           xji_f1f1 = f.pos_c_1[0]-vfibres[i].pos_c_1[0];
           yji_f1f1 = (f.pos_c_1[1]+b_ymax)-vfibres[i].pos_c_1[1];
           xji_f1f2 = f.pos_c_1[0]-vfibres[i].pos_c_2[0];
           yji_f1f2 = (f.pos_c_1[1]+b_ymax)-vfibres[i].pos_c_2[1];
           xji_f2f1 = f.pos_c_2[0]-vfibres[i].pos_c_1[0];
           yji_f2f1 = (f.pos_c_2[1]+b_ymax)-vfibres[i].pos_c_1[1];
           xji_f2f2 = f.pos_c_2[0]-vfibres[i].pos_c_2[0];
           yji_f2f2 = (f.pos_c_2[1]+b_ymax)-vfibres[i].pos_c_2[1];

           CjCi_f1f1 = sqrt(xji_f1f1*xji_f1f1+yji_f1f1*yji_f1f1);
           CjCi_f1f2 = sqrt(xji_f1f2*xji_f1f2+yji_f1f2*yji_f1f2);
           CjCi_f2f1 = sqrt(xji_f2f1*xji_f2f1+yji_f2f1*yji_f2f1);
           CjCi_f2f2 = sqrt(xji_f2f2*xji_f2f2+yji_f2f2*yji_f2f2);

           dn_f1f1 = CjCi_f1f1 - (f.R+vfibres[i].R);
           dn_f1f2 = CjCi_f1f2 - (f.R+vfibres[i].R);
           dn_f2f1 = CjCi_f2f1 - (f.R+vfibres[i].R);
           dn_f2f2 = CjCi_f2f2 - (f.R+vfibres[i].R);

           if(dn_f1f1<0.0 || dn_f1f2<0.0 || dn_f2f1<0.0 || dn_f2f2<0.0){
	            compteur_f++;
	            break;
          }

          // left periodic box
          xji_f1f1 = (f.pos_c_1[0]-b_xmax)-vfibres[i].pos_c_1[0];
          yji_f1f1 = f.pos_c_1[1]-vfibres[i].pos_c_1[1];
          xji_f1f2 = (f.pos_c_1[0]-b_xmax)-vfibres[i].pos_c_2[0];
          yji_f1f2 = f.pos_c_1[1]-vfibres[i].pos_c_2[1];
          xji_f2f1 = (f.pos_c_2[0]-b_xmax)-vfibres[i].pos_c_1[0];
          yji_f2f1 = f.pos_c_2[1]-vfibres[i].pos_c_1[1];
          xji_f2f2 = (f.pos_c_2[0]-b_xmax)-vfibres[i].pos_c_2[0];
          yji_f2f2 = f.pos_c_2[1]-vfibres[i].pos_c_2[1];

          CjCi_f1f1 = sqrt(xji_f1f1*xji_f1f1+yji_f1f1*yji_f1f1);
          CjCi_f1f2 = sqrt(xji_f1f2*xji_f1f2+yji_f1f2*yji_f1f2);
          CjCi_f2f1 = sqrt(xji_f2f1*xji_f2f1+yji_f2f1*yji_f2f1);
          CjCi_f2f2 = sqrt(xji_f2f2*xji_f2f2+yji_f2f2*yji_f2f2);

          dn_f1f1 = CjCi_f1f1 - (f.R+vfibres[i].R);
          dn_f1f2 = CjCi_f1f2 - (f.R+vfibres[i].R);
          dn_f2f1 = CjCi_f2f1 - (f.R+vfibres[i].R);
          dn_f2f2 = CjCi_f2f2 - (f.R+vfibres[i].R);

          if(dn_f1f1<0.0 || dn_f1f2<0.0 || dn_f2f1<0.0 || dn_f2f2<0.0){
            compteur_f++;
            break;
          }

          // right periodic box
          xji_f1f1 = (f.pos_c_1[0]+b_xmax)-vfibres[i].pos_c_1[0];
          yji_f1f1 = f.pos_c_1[1]-vfibres[i].pos_c_1[1];
          xji_f1f2 = (f.pos_c_1[0]+b_xmax)-vfibres[i].pos_c_2[0];
          yji_f1f2 = f.pos_c_1[1]-vfibres[i].pos_c_2[1];
          xji_f2f1 = (f.pos_c_2[0]+b_xmax)-vfibres[i].pos_c_1[0];
          yji_f2f1 = f.pos_c_2[1]-vfibres[i].pos_c_1[1];
          xji_f2f2 = (f.pos_c_2[0]+b_xmax)-vfibres[i].pos_c_2[0];
          yji_f2f2 = f.pos_c_2[1]-vfibres[i].pos_c_2[1];

          CjCi_f1f1 = sqrt(xji_f1f1*xji_f1f1+yji_f1f1*yji_f1f1);
          CjCi_f1f2 = sqrt(xji_f1f2*xji_f1f2+yji_f1f2*yji_f1f2);
          CjCi_f2f1 = sqrt(xji_f2f1*xji_f2f1+yji_f2f1*yji_f2f1);
          CjCi_f2f2 = sqrt(xji_f2f2*xji_f2f2+yji_f2f2*yji_f2f2);

          dn_f1f1 = CjCi_f1f1 - (f.R+vfibres[i].R);
          dn_f1f2 = CjCi_f1f2 - (f.R+vfibres[i].R);
          dn_f2f1 = CjCi_f2f1 - (f.R+vfibres[i].R);
          dn_f2f2 = CjCi_f2f2 - (f.R+vfibres[i].R);

          if(dn_f1f1<0.0 || dn_f1f2<0.0 || dn_f2f1<0.0 || dn_f2f2<0.0){
            compteur_f++;
            break;
          }

          // upper left periodic corner
          xji_f1f1 = (f.pos_c_1[0]-b_xmax)-vfibres[i].pos_c_1[0];
          yji_f1f1 = (f.pos_c_1[1]+b_ymax)-vfibres[i].pos_c_1[1];
          xji_f1f2 = (f.pos_c_1[0]-b_xmax)-vfibres[i].pos_c_2[0];
          yji_f1f2 = (f.pos_c_1[1]+b_ymax)-vfibres[i].pos_c_2[1];
          xji_f2f1 = (f.pos_c_2[0]-b_xmax)-vfibres[i].pos_c_1[0];
          yji_f2f1 = (f.pos_c_2[1]+b_ymax)-vfibres[i].pos_c_1[1];
          xji_f2f2 = (f.pos_c_2[0]-b_xmax)-vfibres[i].pos_c_2[0];
          yji_f2f2 = (f.pos_c_2[1]+b_ymax)-vfibres[i].pos_c_2[1];

          CjCi_f1f1 = sqrt(xji_f1f1*xji_f1f1+yji_f1f1*yji_f1f1);
          CjCi_f1f2 = sqrt(xji_f1f2*xji_f1f2+yji_f1f2*yji_f1f2);
          CjCi_f2f1 = sqrt(xji_f2f1*xji_f2f1+yji_f2f1*yji_f2f1);
          CjCi_f2f2 = sqrt(xji_f2f2*xji_f2f2+yji_f2f2*yji_f2f2);

          dn_f1f1 = CjCi_f1f1 - (f.R+vfibres[i].R);
          dn_f1f2 = CjCi_f1f2 - (f.R+vfibres[i].R);
          dn_f2f1 = CjCi_f2f1 - (f.R+vfibres[i].R);
          dn_f2f2 = CjCi_f2f2 - (f.R+vfibres[i].R);

          if(dn_f1f1<0.0 || dn_f1f2<0.0 || dn_f2f1<0.0 || dn_f2f2<0.0){
            compteur_f++;
            break;
          }

          // lower left periodic corner
          xji_f1f1 = (f.pos_c_1[0]-b_xmax)-vfibres[i].pos_c_1[0];
          yji_f1f1 = (f.pos_c_1[1]-b_ymax)-vfibres[i].pos_c_1[1];
          xji_f1f2 = (f.pos_c_1[0]-b_xmax)-vfibres[i].pos_c_2[0];
          yji_f1f2 = (f.pos_c_1[1]-b_ymax)-vfibres[i].pos_c_2[1];
          xji_f2f1 = (f.pos_c_2[0]-b_xmax)-vfibres[i].pos_c_1[0];
          yji_f2f1 = (f.pos_c_2[1]-b_ymax)-vfibres[i].pos_c_1[1];
          xji_f2f2 = (f.pos_c_2[0]-b_xmax)-vfibres[i].pos_c_2[0];
          yji_f2f2 = (f.pos_c_2[1]-b_ymax)-vfibres[i].pos_c_2[1];

          CjCi_f1f1 = sqrt(xji_f1f1*xji_f1f1+yji_f1f1*yji_f1f1);
          CjCi_f1f2 = sqrt(xji_f1f2*xji_f1f2+yji_f1f2*yji_f1f2);
          CjCi_f2f1 = sqrt(xji_f2f1*xji_f2f1+yji_f2f1*yji_f2f1);
          CjCi_f2f2 = sqrt(xji_f2f2*xji_f2f2+yji_f2f2*yji_f2f2);

          dn_f1f1 = CjCi_f1f1 - (f.R+vfibres[i].R);
          dn_f1f2 = CjCi_f1f2 - (f.R+vfibres[i].R);
          dn_f2f1 = CjCi_f2f1 - (f.R+vfibres[i].R);
          dn_f2f2 = CjCi_f2f2 - (f.R+vfibres[i].R);

          if(dn_f1f1<0.0 || dn_f1f2<0.0 || dn_f2f1<0.0 || dn_f2f2<0.0){
            compteur_f++;
            break;
          }

          // upper right periodic corner
          xji_f1f1 = (f.pos_c_1[0]+b_xmax)-vfibres[i].pos_c_1[0];
          yji_f1f1 = (f.pos_c_1[1]+b_ymax)-vfibres[i].pos_c_1[1];
          xji_f1f2 = (f.pos_c_1[0]+b_xmax)-vfibres[i].pos_c_2[0];
          yji_f1f2 = (f.pos_c_1[1]+b_ymax)-vfibres[i].pos_c_2[1];
          xji_f2f1 = (f.pos_c_2[0]+b_xmax)-vfibres[i].pos_c_1[0];
          yji_f2f1 = (f.pos_c_2[1]+b_ymax)-vfibres[i].pos_c_1[1];
          xji_f2f2 = (f.pos_c_2[0]+b_xmax)-vfibres[i].pos_c_2[0];
          yji_f2f2 = (f.pos_c_2[1]+b_ymax)-vfibres[i].pos_c_2[1];

          CjCi_f1f1 = sqrt(xji_f1f1*xji_f1f1+yji_f1f1*yji_f1f1);
          CjCi_f1f2 = sqrt(xji_f1f2*xji_f1f2+yji_f1f2*yji_f1f2);
          CjCi_f2f1 = sqrt(xji_f2f1*xji_f2f1+yji_f2f1*yji_f2f1);
          CjCi_f2f2 = sqrt(xji_f2f2*xji_f2f2+yji_f2f2*yji_f2f2);

          dn_f1f1 = CjCi_f1f1 - (f.R+vfibres[i].R);
          dn_f1f2 = CjCi_f1f2 - (f.R+vfibres[i].R);
          dn_f2f1 = CjCi_f2f1 - (f.R+vfibres[i].R);
          dn_f2f2 = CjCi_f2f2 - (f.R+vfibres[i].R);

          if(dn_f1f1<0.0 || dn_f1f2<0.0 || dn_f2f1<0.0 || dn_f2f2<0.0){
            compteur_f++;
            break;
          }

          // lower right periodic corner
          xji_f1f1 = (f.pos_c_1[0]+b_xmax)-vfibres[i].pos_c_1[0];
          yji_f1f1 = (f.pos_c_1[1]-b_ymax)-vfibres[i].pos_c_1[1];
          xji_f1f2 = (f.pos_c_1[0]+b_xmax)-vfibres[i].pos_c_2[0];
          yji_f1f2 = (f.pos_c_1[1]-b_ymax)-vfibres[i].pos_c_2[1];
          xji_f2f1 = (f.pos_c_2[0]+b_xmax)-vfibres[i].pos_c_1[0];
          yji_f2f1 = (f.pos_c_2[1]-b_ymax)-vfibres[i].pos_c_1[1];
          xji_f2f2 = (f.pos_c_2[0]+b_xmax)-vfibres[i].pos_c_2[0];
          yji_f2f2 = (f.pos_c_2[1]-b_ymax)-vfibres[i].pos_c_2[1];

          CjCi_f1f1 = sqrt(xji_f1f1*xji_f1f1+yji_f1f1*yji_f1f1);
          CjCi_f1f2 = sqrt(xji_f1f2*xji_f1f2+yji_f1f2*yji_f1f2);
          CjCi_f2f1 = sqrt(xji_f2f1*xji_f2f1+yji_f2f1*yji_f2f1);
          CjCi_f2f2 = sqrt(xji_f2f2*xji_f2f2+yji_f2f2*yji_f2f2);

          dn_f1f1 = CjCi_f1f1 - (f.R+vfibres[i].R);
          dn_f1f2 = CjCi_f1f2 - (f.R+vfibres[i].R);
          dn_f2f1 = CjCi_f2f1 - (f.R+vfibres[i].R);
          dn_f2f2 = CjCi_f2f2 - (f.R+vfibres[i].R);

          if(dn_f1f1<0.0 || dn_f1f2<0.0 || dn_f2f1<0.0 || dn_f2f2<0.0) {
            compteur_f++;
            break;
          }

	  // 2. testing if the end of the candidate fiber collides with the length of the other fiber
    // for the detail of these calculations, see document on Contact detection by JYD

	  // main box
	  double x12, y12;
	  x12 = vfibres[i].pos_c_2[0] - vfibres[i].pos_c_1[0]; // vector C1C2 of fiber i
	  y12 = vfibres[i].pos_c_2[1] - vfibres[i].pos_c_1[1]; //

	  double milieu_f_x = vfibres[i].pos_c_1[0]+x12/2; // center of fiber i
	  double milieu_f_y = vfibres[i].pos_c_2[1]-y12/2; //

    // vector between end1 of candidate fiber and center of fiber i
    double milieu_ff1_x = f.pos_c_1[0] - milieu_f_x;
	  double milieu_ff1_y = f.pos_c_1[1] - milieu_f_y;
    // vector between end2 of candidate fiber and center of fiber i
	  double milieu_ff2_x = f.pos_c_2[0] - milieu_f_x;
	  double milieu_ff2_y = f.pos_c_2[1] - milieu_f_y;

	  std::vector <double> v_milieu_ff1 = {milieu_ff1_x, milieu_ff1_y};
	  std::vector <double> v_milieu_ff2 = {milieu_ff2_x, milieu_ff2_y};

    // normalized vector along fiber i
    std::vector <double> u = {cos(vfibres[i].alpha), sin(vfibres[i].alpha)};

    // dot product
	  double a_f1 = v_milieu_ff1[0]*u[0] + v_milieu_ff1[1]*u[1];
	  double a_f2 = v_milieu_ff2[0]*u[0] + v_milieu_ff2[1]*u[1];
	  double dn, abs_b_f1, abs_b_f2;

    // Case 2.1 ---> does end1 of fiber f overlap along the other fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f1 = v_milieu_ff1[0]*v[0] + v_milieu_ff1[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - vfibres[i].R-f.R;
      if (dn<0.0)
      {
        // overlap
        compteur_f++;
        break;
      }
    }

    // Case 2.2 ---> does end2 of fiber f overlap along the other fiber ?
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f2 = v_milieu_ff2[0]*v[0] + v_milieu_ff2[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - vfibres[i].R-f.R;
      if (dn<0.0)
      {
        // overlap
        compteur_f++;
        break;
      }
    }

    // lower periodic box
    milieu_ff1_x = f.pos_c_1[0] - milieu_f_x;
    milieu_ff1_y = (f.pos_c_1[1]-b_ymax) - milieu_f_y;
    milieu_ff2_x = f.pos_c_2[0] - milieu_f_x;
    milieu_ff2_y = (f.pos_c_2[1]-b_ymax) - milieu_f_y;

    v_milieu_ff1 = {milieu_ff1_x, milieu_ff1_y};
    v_milieu_ff2 = {milieu_ff2_x, milieu_ff2_y};

    a_f1 = v_milieu_ff1[0]*u[0] + v_milieu_ff1[1]*u[1];
    a_f2 = v_milieu_ff2[0]*u[0] + v_milieu_ff2[1]*u[1];

    // Case 2.1 --->does end1 of fiber f overlap along the fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f1 = v_milieu_ff1[0]*v[0] + v_milieu_ff1[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - vfibres[i].R-f.R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Case 2.2 --->does end2 of fiber f overlap along the fiber ?
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f2 = v_milieu_ff2[0]*v[0] + v_milieu_ff2[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - vfibres[i].R-f.R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // upper periodic box
    milieu_ff1_x = f.pos_c_1[0] - milieu_f_x; // Calculs pour déterminer la distance entre l'extrémité 1 de la fibre candidate (f1) et le milieu la fibre f
    milieu_ff1_y = (f.pos_c_1[1]+b_ymax) - milieu_f_y; //
    milieu_ff2_x = f.pos_c_2[0] - milieu_f_x; // Calculs pour déterminer la distance entre l'extrémité 2 de la fibre candidate (f2) et le milieu la fibre f
    milieu_ff2_y = (f.pos_c_2[1]+b_ymax) - milieu_f_y; //

    v_milieu_ff1 = {milieu_ff1_x, milieu_ff1_y};
    v_milieu_ff2 = {milieu_ff2_x, milieu_ff2_y};

    a_f1 = v_milieu_ff1[0]*u[0] + v_milieu_ff1[1]*u[1];
    a_f2 = v_milieu_ff2[0]*u[0] + v_milieu_ff2[1]*u[1];

    // Case 2.1 --->does end1 of fiber f overlap along the fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f1 = v_milieu_ff1[0]*v[0] + v_milieu_ff1[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - vfibres[i].R-f.R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Case 2.2 --->does end2 of fiber f overlap along the fiber ?
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f2 = v_milieu_ff2[0]*v[0] + v_milieu_ff2[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - vfibres[i].R-f.R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // left periodic box
    milieu_ff1_x = (f.pos_c_1[0]-b_xmax) - milieu_f_x;
    milieu_ff1_y = f.pos_c_1[1] - milieu_f_y;
    milieu_ff2_x = (f.pos_c_2[0]-b_xmax) - milieu_f_x;
    milieu_ff2_y = f.pos_c_2[1] - milieu_f_y;

    v_milieu_ff1 = {milieu_ff1_x, milieu_ff1_y};
    v_milieu_ff2 = {milieu_ff2_x, milieu_ff2_y};

    a_f1 = v_milieu_ff1[0]*u[0] + v_milieu_ff1[1]*u[1];
    a_f2 = v_milieu_ff2[0]*u[0] + v_milieu_ff2[1]*u[1];

    // Case 2.1 --->does end1 of fiber f overlap along the fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f1 = v_milieu_ff1[0]*v[0] + v_milieu_ff1[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - vfibres[i].R-f.R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Case 2.2 --->does end2 of fiber f overlap along the fiber ?
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f2 = v_milieu_ff2[0]*v[0] + v_milieu_ff2[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - vfibres[i].R-f.R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // right periodic box
    milieu_ff1_x = (f.pos_c_1[0]+b_xmax) - milieu_f_x;
    milieu_ff1_y = f.pos_c_1[1] - milieu_f_y;
    milieu_ff2_x = (f.pos_c_2[0]+b_xmax) - milieu_f_x;
    milieu_ff2_y = f.pos_c_2[1] - milieu_f_y;

    v_milieu_ff1 = {milieu_ff1_x, milieu_ff1_y};
    v_milieu_ff2 = {milieu_ff2_x, milieu_ff2_y};

    a_f1 = v_milieu_ff1[0]*u[0] + v_milieu_ff1[1]*u[1];
    a_f2 = v_milieu_ff2[0]*u[0] + v_milieu_ff2[1]*u[1];

    // Case 2.1 --->does end1 of fiber f overlap along the fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f1 = v_milieu_ff1[0]*v[0] + v_milieu_ff1[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - vfibres[i].R-f.R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Case 2.2 --->does end2 of fiber f overlap along the fiber ?
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f2 = v_milieu_ff2[0]*v[0] + v_milieu_ff2[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - vfibres[i].R-f.R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // top left periodic corner
    milieu_ff1_x = (f.pos_c_1[0]-b_xmax) - milieu_f_x;
    milieu_ff1_y = f.pos_c_1[1]+b_ymax - milieu_f_y;
    milieu_ff2_x = (f.pos_c_2[0]-b_xmax) - milieu_f_x;
    milieu_ff2_y = f.pos_c_2[1]+b_ymax - milieu_f_y;

    v_milieu_ff1 = {milieu_ff1_x, milieu_ff1_y};
    v_milieu_ff2 = {milieu_ff2_x, milieu_ff2_y};

    a_f1 = v_milieu_ff1[0]*u[0] + v_milieu_ff1[1]*u[1];
    a_f2 = v_milieu_ff2[0]*u[0] + v_milieu_ff2[1]*u[1];

    // Case 2.1 --->does end1 of fiber f overlap along the fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f1 = v_milieu_ff1[0]*v[0] + v_milieu_ff1[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - vfibres[i].R-f.R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Case 2.2 --->does end2 of fiber f overlap along the fiber ?
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f2 = v_milieu_ff2[0]*v[0] + v_milieu_ff2[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - vfibres[i].R-f.R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // top right periodic corner
    milieu_ff1_x = (f.pos_c_1[0]+b_xmax) - milieu_f_x;
    milieu_ff1_y = f.pos_c_1[1]+b_ymax - milieu_f_y;
    milieu_ff2_x = (f.pos_c_2[0]+b_xmax) - milieu_f_x;
    milieu_ff2_y = f.pos_c_2[1]+b_ymax - milieu_f_y;

    v_milieu_ff1 = {milieu_ff1_x, milieu_ff1_y};
    v_milieu_ff2 = {milieu_ff2_x, milieu_ff2_y};

    a_f1 = v_milieu_ff1[0]*u[0] + v_milieu_ff1[1]*u[1];
    a_f2 = v_milieu_ff2[0]*u[0] + v_milieu_ff2[1]*u[1];

    // Case 2.1 --->does end1 of fiber f overlap along the fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f1 = v_milieu_ff1[0]*v[0] + v_milieu_ff1[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - vfibres[i].R-f.R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Case 2.2 --->does end2 of fiber f overlap along the fiber ?
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f2 = v_milieu_ff2[0]*v[0] + v_milieu_ff2[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - vfibres[i].R-f.R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // lower left periodic corner
    milieu_ff1_x = (f.pos_c_1[0]-b_xmax) - milieu_f_x;
    milieu_ff1_y = f.pos_c_1[1]-b_ymax - milieu_f_y;
    milieu_ff2_x = (f.pos_c_2[0]-b_xmax) - milieu_f_x;
    milieu_ff2_y = f.pos_c_2[1]-b_ymax - milieu_f_y;

    v_milieu_ff1 = {milieu_ff1_x, milieu_ff1_y};
    v_milieu_ff2 = {milieu_ff2_x, milieu_ff2_y};

    a_f1 = v_milieu_ff1[0]*u[0] + v_milieu_ff1[1]*u[1];
    a_f2 = v_milieu_ff2[0]*u[0] + v_milieu_ff2[1]*u[1];

    // Case 2.1 --->does end1 of fiber f overlap along the fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f1 = v_milieu_ff1[0]*v[0] + v_milieu_ff1[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - vfibres[i].R-f.R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Case 2.2 --->does end2 of fiber f overlap along the fiber ?
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f2 = v_milieu_ff2[0]*v[0] + v_milieu_ff2[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - vfibres[i].R-f.R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // lower right periodic corner
    milieu_ff1_x = (f.pos_c_1[0]+b_xmax) - milieu_f_x;
    milieu_ff1_y = f.pos_c_1[1]-b_ymax - milieu_f_y;
    milieu_ff2_x = (f.pos_c_2[0]+b_xmax) - milieu_f_x;
    milieu_ff2_y = f.pos_c_2[1]-b_ymax - milieu_f_y;

    v_milieu_ff1 = {milieu_ff1_x, milieu_ff1_y};
    v_milieu_ff2 = {milieu_ff2_x, milieu_ff2_y};

    a_f1 = v_milieu_ff1[0]*u[0] + v_milieu_ff1[1]*u[1];
    a_f2 = v_milieu_ff2[0]*u[0] + v_milieu_ff2[1]*u[1];

    // Case 2.1 --->does end1 of fiber f overlap along the fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f1 = v_milieu_ff1[0]*v[0] + v_milieu_ff1[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - vfibres[i].R-f.R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Case 2.2 ---> does end2 of fiber f overlap along the fiber ?
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
      double b_f2 = v_milieu_ff2[0]*v[0] + v_milieu_ff2[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - vfibres[i].R-f.R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // do any ends of fiber f (already in the box) overlap along the candidate fiber ?

    // main box
    x12 = f.pos_c_2[0] - f.pos_c_1[0]; // C2C1 segment of the candidate fiber
    y12 = f.pos_c_2[1] - f.pos_c_1[1]; //

    milieu_f_x = f.pos_c_1[0]+x12/2; // coordinates of the center of the candidate fiber
    milieu_f_y = f.pos_c_2[1]-y12/2; //

    // dist between center of candidate fiber and end1 of fiber i
    double milieu_f1f_x = vfibres[i].pos_c_1[0] - milieu_f_x;
    double milieu_f1f_y = vfibres[i].pos_c_1[1] - milieu_f_y;
    // dist between center of candidate fiber and end2 of fiber i
    double milieu_f2f_x = vfibres[i].pos_c_2[0] - milieu_f_x;
    double milieu_f2f_y = vfibres[i].pos_c_2[1] - milieu_f_y;

    std::vector <double> v_milieu_f1f = {milieu_f1f_x, milieu_f1f_y};
    std::vector <double> v_milieu_f2f = {milieu_f2f_x, milieu_f2f_y};

    // norm vector along candidate fiber
    u = {cos(f.alpha), sin(f.alpha)};

    // dot product
    a_f1 = v_milieu_f1f[0]*u[0] + v_milieu_f1f[1]*u[1];
    a_f2 = v_milieu_f2f[0]*u[0] + v_milieu_f2f[1]*u[1];

    // Case 2.1 ---> does end1 of fiber i overlap along candidate fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f1 = v_milieu_f1f[0]*v[0] + v_milieu_f1f[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Case 2.2 ---> does end2 of fiber i overlap along candidate fiber ?
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f2 = v_milieu_f2f[0]*v[0] + v_milieu_f2f[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // lower perdiodic box
    /// x12 = f.pos_c_2[0] - f.pos_c_1[0];
    /// y12 = (f.pos_c_2[1]-b_ymax) - (f.pos_c_1[1]-b_ymax);

    /// milieu_f_x = f.pos_c_1[0]+x12/2;
    /// milieu_f_y = (f.pos_c_2[1]-b_ymax)-y12/2; //

    milieu_f1f_x = vfibres[i].pos_c_1[0] - milieu_f_x;
    milieu_f1f_y = vfibres[i].pos_c_1[1] - (milieu_f_y - b_ymax);
    milieu_f2f_x = vfibres[i].pos_c_2[0] - milieu_f_x;
    milieu_f2f_y = vfibres[i].pos_c_2[1] -  (milieu_f_y - b_ymax);

    v_milieu_f1f = {milieu_f1f_x, milieu_f1f_y};
    v_milieu_f2f = {milieu_f2f_x, milieu_f2f_y};

    u = {cos(f.alpha), sin(f.alpha)};

    a_f1 = v_milieu_f1f[0]*u[0] + v_milieu_f1f[1]*u[1];
    a_f2 = v_milieu_f2f[0]*u[0] + v_milieu_f2f[1]*u[1];

    // Case 2.1 ---> does end1 of fiber i overlap along candidate fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f1 = v_milieu_f1f[0]*v[0] + v_milieu_f1f[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Cas 2.2 ---> extrémité 2 de la fibre f recouvre la fibre sur sa longueur
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f2 = v_milieu_f2f[0]*v[0] + v_milieu_f2f[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Upper periodic box
    /// x12 = f.pos_c_2[0] - f.pos_c_1[0];
    /// y12 = (f.pos_c_2[1]+b_ymax) - (f.pos_c_1[1]+b_ymax);

    /// milieu_f_x = f.pos_c_1[0]+x12/2;
    /// milieu_f_y = (f.pos_c_2[1]+b_ymax)-y12/2;

    milieu_f1f_x = vfibres[i].pos_c_1[0] - milieu_f_x;
    milieu_f1f_y = vfibres[i].pos_c_1[1] - (milieu_f_y + b_ymax);
    milieu_f2f_x = vfibres[i].pos_c_2[0] - milieu_f_x;
    milieu_f2f_y = vfibres[i].pos_c_2[1] - (milieu_f_y + b_ymax);

    v_milieu_f1f = {milieu_f1f_x, milieu_f1f_y};
    v_milieu_f2f = {milieu_f2f_x, milieu_f2f_y};

    u = {cos(f.alpha), sin(f.alpha)};

    a_f1 = v_milieu_f1f[0]*u[0] + v_milieu_f1f[1]*u[1];
    a_f2 = v_milieu_f2f[0]*u[0] + v_milieu_f2f[1]*u[1];

    // Case 2.1 ---> does end1 of fiber i overlap along candidate fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f1 = v_milieu_f1f[0]*v[0] + v_milieu_f1f[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Cas 2.2 ---> extrémité 2 de la fibre f recouvre la fibre sur sa longueur
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f2 = v_milieu_f2f[0]*v[0] + v_milieu_f2f[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // left periodic box
    /// x12 = (f.pos_c_2[0]-b_xmax) - (f.pos_c_1[0]-b_xmax);
    /// y12 = f.pos_c_2[1] - f.pos_c_1[1];

    /// milieu_f_x = (f.pos_c_1[0]-b_xmax)+x12/2;
    /// milieu_f_y = f.pos_c_2[1]-y12/2;

    milieu_f1f_x = vfibres[i].pos_c_1[0] - (milieu_f_x - b_xmax);
    milieu_f1f_y = vfibres[i].pos_c_1[1] - milieu_f_y;
    milieu_f2f_x = vfibres[i].pos_c_2[0] - (milieu_f_x - b_xmax);
    milieu_f2f_y = vfibres[i].pos_c_2[1] - milieu_f_y;

    v_milieu_f1f = {milieu_f1f_x, milieu_f1f_y};
    v_milieu_f2f = {milieu_f2f_x, milieu_f2f_y};

    u = {cos(f.alpha), sin(f.alpha)};

    a_f1 = v_milieu_f1f[0]*u[0] + v_milieu_f1f[1]*u[1];
    a_f2 = v_milieu_f2f[0]*u[0] + v_milieu_f2f[1]*u[1];

    // Case 2.1 ---> does end1 of fiber i overlap along candidate fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f1 = v_milieu_f1f[0]*v[0] + v_milieu_f1f[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Cas 2.2 ---> extrémité 2 de la fibre f recouvre la fibre sur sa longueur
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f2 = v_milieu_f2f[0]*v[0] + v_milieu_f2f[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // right periodic box
    /// x12 = (f.pos_c_2[0]+b_xmax) - (f.pos_c_1[0]+b_xmax);
    /// y12 = f.pos_c_2[1] - f.pos_c_1[1];

    /// milieu_f_x = (f.pos_c_1[0]+b_xmax)+x12/2;
    /// milieu_f_y = f.pos_c_2[1]-y12/2;

    milieu_f1f_x = vfibres[i].pos_c_1[0] - (milieu_f_x + b_xmax);
    milieu_f1f_y = vfibres[i].pos_c_1[1] - milieu_f_y;
    milieu_f2f_x = vfibres[i].pos_c_2[0] - (milieu_f_x + b_xmax);
    milieu_f2f_y = vfibres[i].pos_c_2[1] - milieu_f_y;

    v_milieu_f1f = {milieu_f1f_x, milieu_f1f_y};
    v_milieu_f2f = {milieu_f2f_x, milieu_f2f_y};

    u = {cos(f.alpha), sin(f.alpha)};

    a_f1 = v_milieu_f1f[0]*u[0] + v_milieu_f1f[1]*u[1];
    a_f2 = v_milieu_f2f[0]*u[0] + v_milieu_f2f[1]*u[1];

    // Case 2.1 ---> does end1 of fiber i overlap along candidate fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f1 = v_milieu_f1f[0]*v[0] + v_milieu_f1f[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Cas 2.2 ---> extrémité 2 de la fibre f recouvre la fibre sur sa longueur
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f2 = v_milieu_f2f[0]*v[0] + v_milieu_f2f[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }


    // top right periodic corner
    /// x12 = (f.pos_c_2[0]+b_xmax) - (f.pos_c_1[0]+b_xmax);
    /// y12 = (f.pos_c_2[1]+b_ymax) - (f.pos_c_1[1]+b_ymax);

    /// milieu_f_x = (f.pos_c_1[0]+b_xmax)+x12/2;
    /// milieu_f_y = (f.pos_c_2[1]+b_ymax)-y12/2;

    milieu_f1f_x = vfibres[i].pos_c_1[0] - (milieu_f_x + b_xmax);
    milieu_f1f_y = vfibres[i].pos_c_1[1] - (milieu_f_y + b_ymax);
    milieu_f2f_x = vfibres[i].pos_c_2[0] - (milieu_f_x + b_xmax);
    milieu_f2f_y = vfibres[i].pos_c_2[1] - (milieu_f_y + b_ymax);

    v_milieu_f1f = {milieu_f1f_x, milieu_f1f_y};
    v_milieu_f2f = {milieu_f2f_x, milieu_f2f_y};

    u = {cos(f.alpha), sin(f.alpha)};

    a_f1 = v_milieu_f1f[0]*u[0] + v_milieu_f1f[1]*u[1];
    a_f2 = v_milieu_f2f[0]*u[0] + v_milieu_f2f[1]*u[1];

    // Case 2.1 ---> does end1 of fiber i overlap along candidate fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f1 = v_milieu_f1f[0]*v[0] + v_milieu_f1f[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Cas 2.2 ---> extrémité 2 de la fibre f recouvre la fibre sur sa longueur
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f2 = v_milieu_f2f[0]*v[0] + v_milieu_f2f[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // top left periodic corner
    /// x12 = (f.pos_c_2[0]-b_xmax) - (f.pos_c_1[0]-b_xmax);
    /// y12 = (f.pos_c_2[1]+b_ymax) - (f.pos_c_1[1]+b_ymax);

    /// milieu_f_x = (f.pos_c_1[0]-b_xmax)+x12/2;
    /// milieu_f_y = (f.pos_c_2[1]+b_ymax)-y12/2;

    milieu_f1f_x = vfibres[i].pos_c_1[0] - (milieu_f_x - b_xmax);
    milieu_f1f_y = vfibres[i].pos_c_1[1] - (milieu_f_y + b_ymax);
    milieu_f2f_x = vfibres[i].pos_c_2[0] - (milieu_f_x - b_xmax);
    milieu_f2f_y = vfibres[i].pos_c_2[1] - (milieu_f_y + b_ymax);

    v_milieu_f1f = {milieu_f1f_x, milieu_f1f_y};
    v_milieu_f2f = {milieu_f2f_x, milieu_f2f_y};

    u = {cos(f.alpha), sin(f.alpha)};

    a_f1 = v_milieu_f1f[0]*u[0] + v_milieu_f1f[1]*u[1];
    a_f2 = v_milieu_f2f[0]*u[0] + v_milieu_f2f[1]*u[1];

    // Case 2.1 ---> does end1 of fiber i overlap along candidate fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f1 = v_milieu_f1f[0]*v[0] + v_milieu_f1f[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Cas 2.2 ---> extrémité 2 de la fibre f recouvre la fibre sur sa longueur
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f2 = v_milieu_f2f[0]*v[0] + v_milieu_f2f[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // lower right periodic corner
    /// x12 = (f.pos_c_2[0]+b_xmax) - (f.pos_c_1[0]+b_xmax);
    /// y12 = (f.pos_c_2[1]-b_ymax) - (f.pos_c_1[1]-b_ymax);

    /// milieu_f_x = (f.pos_c_1[0]+b_xmax)+x12/2;
    /// milieu_f_y = (f.pos_c_2[1]-b_ymax)-y12/2;

    milieu_f1f_x = vfibres[i].pos_c_1[0] - (milieu_f_x + b_xmax);
    milieu_f1f_y = vfibres[i].pos_c_1[1] - (milieu_f_y - b_ymax);
    milieu_f2f_x = vfibres[i].pos_c_2[0] - (milieu_f_x + b_xmax);
    milieu_f2f_y = vfibres[i].pos_c_2[1] - (milieu_f_y - b_ymax);

    v_milieu_f1f = {milieu_f1f_x, milieu_f1f_y};
    v_milieu_f2f = {milieu_f2f_x, milieu_f2f_y};

    u = {cos(f.alpha), sin(f.alpha)};

    a_f1 = v_milieu_f1f[0]*u[0] + v_milieu_f1f[1]*u[1];
    a_f2 = v_milieu_f2f[0]*u[0] + v_milieu_f2f[1]*u[1];

    // Case 2.1 ---> does end1 of fiber i overlap along candidate fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f1 = v_milieu_f1f[0]*v[0] + v_milieu_f1f[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Cas 2.2 ---> extrémité 2 de la fibre f recouvre la fibre sur sa longueur
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f2 = v_milieu_f2f[0]*v[0] + v_milieu_f2f[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // lower left periodic corner
    /// x12 = (f.pos_c_2[0]-b_xmax) - (f.pos_c_1[0]-b_xmax);
    /// y12 = (f.pos_c_2[1]-b_ymax) - (f.pos_c_1[1]-b_ymax);

    /// milieu_f_x = (f.pos_c_1[0]-b_xmax)+x12/2;
    /// milieu_f_y = (f.pos_c_2[1]-b_ymax)-y12/2;

    milieu_f1f_x = vfibres[i].pos_c_1[0] - (milieu_f_x - b_xmax);
    milieu_f1f_y = vfibres[i].pos_c_1[1] - (milieu_f_y - b_ymax);
    milieu_f2f_x = vfibres[i].pos_c_2[0] - (milieu_f_x - b_xmax);
    milieu_f2f_y = vfibres[i].pos_c_2[1] - (milieu_f_y - b_ymax);

    v_milieu_f1f = {milieu_f1f_x, milieu_f1f_y};
    v_milieu_f2f = {milieu_f2f_x, milieu_f2f_y};

    u = {cos(f.alpha), sin(f.alpha)};

    a_f1 = v_milieu_f1f[0]*u[0] + v_milieu_f1f[1]*u[1];
    a_f2 = v_milieu_f2f[0]*u[0] + v_milieu_f2f[1]*u[1];

    // Case 2.1 ---> does end1 of fiber i overlap along candidate fiber ?
    if (a_f1>-L/2 && a_f1<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f1 = v_milieu_f1f[0]*v[0] + v_milieu_f1f[1]*v[1];

      if (b_f1<0)
      {
        abs_b_f1 = -b_f1;
      }
      else
      {
        abs_b_f1 = b_f1;
      }

      dn = abs_b_f1 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }

    // Cas 2.2 ---> extrémité 2 de la fibre f recouvre la fibre sur sa longueur
    if (a_f2>-L/2 && a_f2<L/2)
    {
      std::vector <double> v = {-sin(f.alpha), cos(f.alpha)};
      double b_f2 = v_milieu_f2f[0]*v[0] + v_milieu_f2f[1]*v[1];

      if (b_f2<0)
      {
        abs_b_f2 = -b_f2;
      }
      else
      {
        abs_b_f2 = b_f2;
      }

      dn = abs_b_f2 - f.R-vfibres[i].R;
      if (dn<0.0) // contact
      {
        compteur_f++;
        break;
      }
    }
  }

  if (compteur_f == 0)
  {
    // no contact, ACCEPT new fiber
    vfibres.push_back(f);  // store the fiber
    // compute the coordinates of the points for the description of the fiber as a polygon (paraview)
    double x12, y12;
    for (int th=0; th<nth+1; th++)
    {
      double xp_1, yp_1, xp_2, yp_2;
      x12 = f.pos_c_2[0] - f.pos_c_1[0];
      y12 = f.pos_c_2[1] - f.pos_c_1[1];
      xp_1 = f.pos_c_1[0] + f.R*cos(th*dth+alpha+PI/2);
      yp_1 = f.pos_c_1[1] + f.R*sin(th*dth+alpha+PI/2);
      xp_2 = f.pos_c_2[0] + f.R*cos(th*dth+alpha-PI/2);
      yp_2 = f.pos_c_2[1] + f.R*sin(th*dth+alpha-PI/2);
      p_f.pos_1.push_back(xp_1);
      p_f.pos_1.push_back(yp_1);
      p_f.pos_2.push_back(xp_2);
      p_f.pos_2.push_back(yp_2);
    }
    vp_f.push_back(p_f);
    p_f.pos_1.clear();
    p_f.pos_2.clear();
  }
} // end loop for inserting fibers


  ////////////
  // GRAINS //
  ////////////

/// @@@@@@@@@

  float R_min = R_f;
  float R_max = sqrt(2)*R_f;

  // Draw in random uniform laws the sizes and positions of the grains
  std::uniform_real_distribution<float> distribution_pos_x(b_xmin, b_xmax);
  std::uniform_real_distribution<float> distribution_pos_y(b_ymin, b_ymax);
  std::uniform_real_distribution<float> distribution_R(R_min, R_max);
//radius order  
  std::vector<float> radii;
for (int i = 0; i < ng; ++i) {
    radii.push_back(distribution_R(gen));
}
std::sort(radii.begin(), radii.end(), std::greater<float>());
//just for showing the radius
//for (int i = 0; i < ng; ++i) {
  //  cout << "grain radius: " << radii[i] << endl;
//}

// for consider the radius array
int radiusIndex = 0;
  int att_g = 0;
  while (vgrains.size()<ng) // While all the grains have been put in the box
  {
    att_g++;
    int compteur = 0;
    // properties of the candidate grain
    float x_al = distribution_pos_x(gen);
    float y_al = distribution_pos_y(gen);
    //float R_al = distribution_R(gen);
    float R_al = radii[radiusIndex];
    grain g = {R_al, x_al, y_al};
    //showing the structure g
    //std::cout << "trying to insert grain " << vgrains.size() + 1 << ": " << g << std::endl;
    
    // Testing if grains collide
    for (i=0; i<vgrains.size(); i++)
    {
      // main box
      float xji = g.x-vgrains[i].x;
      float yji = g.y-vgrains[i].y;
      float CjCi = sqrt(xji*xji+yji*yji);
      float dn = CjCi - (g.R+vgrains[i].R);
      if(dn<0.) // Contact
      {
        compteur++;
        break;
      }

      // lower periodic box
      xji = g.x-vgrains[i].x;
      yji = (g.y-b_ymax)-vgrains[i].y;
      CjCi = sqrt(xji*xji+yji*yji);
      dn = CjCi - (g.R+vgrains[i].R);
      if(dn<0.) // Contact
      {
        compteur++;
        break;
      }

      // upper periodic box
      xji = g.x-vgrains[i].x;
      yji = (g.y+b_ymax)-vgrains[i].y;
      CjCi = sqrt(xji*xji+yji*yji);
      dn = CjCi - (g.R+vgrains[i].R);
      if(dn<0.) // Contact
      {
        compteur++;
        break;
      }

      // left periodic box
      yji = g.y-vgrains[i].y;
      xji = (g.x-b_xmax)-vgrains[i].x;
      CjCi = sqrt(xji*xji+yji*yji);
      dn = CjCi - (g.R+vgrains[i].R);
      if(dn<0.) // Contact
      {
        compteur++;
        break;
      }

      // right periodic box
      yji = g.y-vgrains[i].y;
      xji = (g.x+b_xmax)-vgrains[i].x;
      CjCi = sqrt(xji*xji+yji*yji);
      dn = CjCi - (g.R+vgrains[i].R);
      if(dn<0.) // Contact
      {
        compteur++;
        break;
      }

      // upper left periodic corner
      yji = (g.y+b_ymax)-vgrains[i].y;
      xji = (g.x-b_xmax)-vgrains[i].x;
      CjCi = sqrt(xji*xji+yji*yji);
      dn = CjCi - (g.R+vgrains[i].R);
      if(dn<0.) // Contact
      {
        compteur++;
        break;
      }

      // upper right periodic corner
      yji = (g.y+b_ymax)-vgrains[i].y;
      xji = (g.x+b_xmax)-vgrains[i].x;
      CjCi = sqrt(xji*xji+yji*yji);
      dn = CjCi - (g.R+vgrains[i].R);
      if(dn<0.) // Contact
      {
        compteur++;
        break;
      }

      // lower right periodic corner
      yji = (g.y-b_ymax)-vgrains[i].y;
      xji = (g.x+b_xmax)-vgrains[i].x;
      CjCi = sqrt(xji*xji+yji*yji);
      dn = CjCi - (g.R+vgrains[i].R);
      if(dn<0.) // Contact
      {
        compteur++;
        break;
      }

      // lower left periodic corner
      yji = (g.y-b_ymax)-vgrains[i].y;
      xji = (g.x-b_xmax)-vgrains[i].x;
      CjCi = sqrt(xji*xji+yji*yji);
      dn = CjCi - (g.R+vgrains[i].R);
      if(dn<0.) // Contact
      {
        compteur++;
        break;
      }

    }

    // testing if the candidate grain overlaps with fibers

    // main box
    for (i=0; i<vfibres.size(); i++)
    {
      double x12, y12;
      x12 = vfibres[i].pos_c_2[0] - vfibres[i].pos_c_1[0];
      y12 = vfibres[i].pos_c_2[1] - vfibres[i].pos_c_1[1];

      double milieu_f_x = vfibres[i].pos_c_1[0]+x12/2;
      double milieu_f_y = vfibres[i].pos_c_2[1]-y12/2;
      double milieu_fg_x = g.x - milieu_f_x;
      double milieu_fg_y = g.y - milieu_f_y;
      std::vector <double> v_milieu_fg = {milieu_fg_x, milieu_fg_y};

      std::vector <double> u = {cos(vfibres[i].alpha), sin(vfibres[i].alpha)};
      std::vector <double> v_gf;

      double a = v_milieu_fg[0]*u[0] + v_milieu_fg[1]*u[1];
      double dn, abs_b, v_gf_x, v_gf_y;

      if (a<-L/2)
      {
        v_gf_x = vfibres[i].pos_c_1[0] - g.x;
        v_gf_y = vfibres[i].pos_c_1[1] - g.y;
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1])-vfibres[i].R-g.R;
        if (dn<0.)
        {
          compteur++;
          break;
        }
      }

      if (a>-L/2 && a<L/2)
      {
        std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
        double b = v_milieu_fg[0]*v[0] + v_milieu_fg[1]*v[1];

        if (b<0)
        {
          abs_b = -b;
        }
        else
        {
          abs_b = b;
        }

        dn = abs_b - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }

      if(a>L/2)
      {
        v_gf_x = vfibres[i].pos_c_2[0]-g.x;
        v_gf_y = vfibres[i].pos_c_2[1]-g.y;
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1]) - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }

      // lower periodic box
      milieu_fg_x = g.x - milieu_f_x;
      milieu_fg_y = (g.y-b_ymax) - milieu_f_y;
      v_milieu_fg = {milieu_fg_x, milieu_fg_y};

      a = v_milieu_fg[0]*u[0] + v_milieu_fg[1]*u[1];

      if (a<-L/2)
      {
        v_gf_x = vfibres[i].pos_c_1[0] - g.x;
        v_gf_y = vfibres[i].pos_c_1[1] - (g.y-b_ymax);
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1])-vfibres[i].R-g.R;
        if (dn<0.)
        {
          compteur++;
          break;
        }
      }

      if (a>-L/2 && a<L/2)
      {
        std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
        double b = v_milieu_fg[0]*v[0] + v_milieu_fg[1]*v[1];

        if (b<0)
        {
          abs_b = -b;
        }
        else
        {
          abs_b = b;
        }

        dn = abs_b - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }

      if(a>L/2)
      {
        v_gf_x = vfibres[i].pos_c_2[0]-g.x;
        v_gf_y = vfibres[i].pos_c_2[1]-(g.y-b_ymax);
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1]) - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }

      // upper periodic box
      milieu_fg_x = g.x - milieu_f_x;
      milieu_fg_y = (g.y+b_ymax) - milieu_f_y;
      v_milieu_fg = {milieu_fg_x, milieu_fg_y};

      a = v_milieu_fg[0]*u[0] + v_milieu_fg[1]*u[1];

      if (a<-L/2)
      {
        v_gf_x = vfibres[i].pos_c_1[0] - g.x;
        v_gf_y = vfibres[i].pos_c_1[1] - (g.y+b_ymax);
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1])-vfibres[i].R-g.R;
        if (dn<0.)
        {
          compteur++;
          break;
        }
      }

      if (a>-L/2 && a<L/2)
      {
        std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
        double b = v_milieu_fg[0]*v[0] + v_milieu_fg[1]*v[1];

        if (b<0)
        {
          abs_b = -b;
        }
        else
        {
          abs_b = b;
        }

        dn = abs_b - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }

      if(a>L/2)
      {
        v_gf_x = vfibres[i].pos_c_2[0]-g.x;
        v_gf_y = vfibres[i].pos_c_2[1]-(g.y+b_ymax);
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1]) - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }


      // left periodic box
      milieu_fg_x = (g.x-b_xmax) - milieu_f_x;
      milieu_fg_y = g.y - milieu_f_y;
      v_milieu_fg = {milieu_fg_x, milieu_fg_y};

      a = v_milieu_fg[0]*u[0] + v_milieu_fg[1]*u[1];

      if (a<-L/2)
      {
        v_gf_x = vfibres[i].pos_c_1[0] - (g.x-b_xmax);
        v_gf_y = vfibres[i].pos_c_1[1] - g.y;
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1])-vfibres[i].R-g.R;
        if (dn<0.)
        {
          compteur++;
          break;
        }
      }

      if (a>-L/2 && a<L/2)
      {
        std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
        double b = v_milieu_fg[0]*v[0] + v_milieu_fg[1]*v[1];

        if (b<0)
        {
          abs_b = -b;
        }
        else
        {
          abs_b = b;
        }

        dn = abs_b - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }

      if(a>L/2)
      {
        v_gf_x = vfibres[i].pos_c_2[0]-(g.x-b_xmax);
        v_gf_y = vfibres[i].pos_c_2[1]-g.y;
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1]) - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }


      // right periodic box
      milieu_fg_x = (g.x+b_xmax) - milieu_f_x;
      milieu_fg_y = g.y - milieu_f_y;
      v_milieu_fg = {milieu_fg_x, milieu_fg_y};

      a = v_milieu_fg[0]*u[0] + v_milieu_fg[1]*u[1];

      if (a<-L/2)
      {
        v_gf_x = vfibres[i].pos_c_1[0] - (g.x+b_xmax);
        v_gf_y = vfibres[i].pos_c_1[1] - g.y;
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1])-vfibres[i].R-g.R;
        if (dn<0.)
        {
          compteur++;
          break;
        }
      }

      if (a>-L/2 && a<L/2)
      {
        std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
        double b = v_milieu_fg[0]*v[0] + v_milieu_fg[1]*v[1];

        if (b<0)
        {
          abs_b = -b;
        }
        else
        {
          abs_b = b;
        }

        dn = abs_b - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }

      if(a>L/2)
      {
        v_gf_x = vfibres[i].pos_c_2[0]-(g.x+b_xmax);
        v_gf_y = vfibres[i].pos_c_2[1]-g.y;
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1]) - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }

      // upper left periodic corner
      milieu_fg_x = (g.x-b_xmax) - milieu_f_x;
      milieu_fg_y = (g.y+b_ymax) - milieu_f_y;
      v_milieu_fg = {milieu_fg_x, milieu_fg_y};

      a = v_milieu_fg[0]*u[0] + v_milieu_fg[1]*u[1];

      if (a<-L/2)
      {
        v_gf_x = vfibres[i].pos_c_1[0] - (g.x-b_xmax);
        v_gf_y = vfibres[i].pos_c_1[1] - (g.y+b_ymax);
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1])-vfibres[i].R-g.R;
        if (dn<0.)
        {
          compteur++;
          break;
        }
      }

      if (a>-L/2 && a<L/2)
      {
        std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
        double b = v_milieu_fg[0]*v[0] + v_milieu_fg[1]*v[1];

        if (b<0)
        {
          abs_b = -b;
        }
        else
        {
          abs_b = b;
        }

        dn = abs_b - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }

      if(a>L/2)
      {
        v_gf_x = vfibres[i].pos_c_2[0]-(g.x-b_xmax);
        v_gf_y = vfibres[i].pos_c_2[1]-(g.y+b_ymax);
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1]) - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }

      // upper right periodic corner
      milieu_fg_x = (g.x+b_xmax) - milieu_f_x;
      milieu_fg_y = (g.y+b_ymax) - milieu_f_y;
      v_milieu_fg = {milieu_fg_x, milieu_fg_y};

      a = v_milieu_fg[0]*u[0] + v_milieu_fg[1]*u[1];

      if (a<-L/2)
      {
        v_gf_x = vfibres[i].pos_c_1[0] - (g.x+b_xmax);
        v_gf_y = vfibres[i].pos_c_1[1] - (g.y+b_ymax);
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1])-vfibres[i].R-g.R;
        if (dn<0.)
        {
          compteur++;
          break;
        }
      }

      if (a>-L/2 && a<L/2)
      {
        std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
        double b = v_milieu_fg[0]*v[0] + v_milieu_fg[1]*v[1];

        if (b<0)
        {
          abs_b = -b;
        }
        else
        {
          abs_b = b;
        }

        dn = abs_b - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }

      if(a>L/2)
      {
        v_gf_x = vfibres[i].pos_c_2[0]-(g.x+b_xmax);
        v_gf_y = vfibres[i].pos_c_2[1]-(g.y+b_ymax);
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1]) - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }

      // lower left periodic corner
      milieu_fg_x = (g.x-b_xmax) - milieu_f_x;
      milieu_fg_y = (g.y-b_ymax) - milieu_f_y;
      v_milieu_fg = {milieu_fg_x, milieu_fg_y};

      a = v_milieu_fg[0]*u[0] + v_milieu_fg[1]*u[1];

      if (a<-L/2)
      {
        v_gf_x = vfibres[i].pos_c_1[0] - (g.x-b_xmax);
        v_gf_y = vfibres[i].pos_c_1[1] - (g.y-b_ymax);
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1])-vfibres[i].R-g.R;
        if (dn<0.)
        {
          compteur++;
          break;
        }
      }

      if (a>-L/2 && a<L/2)
      {
        std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
        double b = v_milieu_fg[0]*v[0] + v_milieu_fg[1]*v[1];

        if (b<0)
        {
          abs_b = -b;
        }
        else
        {
          abs_b = b;
        }

        dn = abs_b - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }

      if(a>L/2)
      {
        v_gf_x = vfibres[i].pos_c_2[0]-(g.x-b_xmax);
        v_gf_y = vfibres[i].pos_c_2[1]-(g.y-b_ymax);
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1]) - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }

      // lower right periodic corner
      milieu_fg_x = (g.x+b_xmax) - milieu_f_x;
      milieu_fg_y = (g.y-b_ymax) - milieu_f_y;
      v_milieu_fg = {milieu_fg_x, milieu_fg_y};

      a = v_milieu_fg[0]*u[0] + v_milieu_fg[1]*u[1];

      if (a<-L/2)
      {
        v_gf_x = vfibres[i].pos_c_1[0] - (g.x+b_xmax);
        v_gf_y = vfibres[i].pos_c_1[1] - (g.y-b_ymax);
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1])-vfibres[i].R-g.R;
        if (dn<0.)
        {
          compteur++;
          break;
        }
      }

      if (a>-L/2 && a<L/2)
      {
        std::vector <double> v = {-sin(vfibres[i].alpha), cos(vfibres[i].alpha)};
        double b = v_milieu_fg[0]*v[0] + v_milieu_fg[1]*v[1];

        if (b<0)
        {
          abs_b = -b;
        }
        else
        {
          abs_b = b;
        }

        dn = abs_b - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }

      if(a>L/2)
      {
        v_gf_x = vfibres[i].pos_c_2[0]-(g.x+b_xmax);
        v_gf_y = vfibres[i].pos_c_2[1]-(g.y-b_ymax);
        v_gf = {v_gf_x, v_gf_y};

        dn = sqrt(v_gf[0]*v_gf[0] + v_gf[1]*v_gf[1]) - vfibres[i].R-g.R;
        if (dn<0)
        {
          compteur++;
          break;
        }
      }

    }

    if (compteur==0)
    {
      vgrains.push_back(g);
      radiusIndex++;
    }
  }


  cout << "Number de fibres : " << vfibres.size() << " ; " << "Number of grains : " << vgrains.size() << endl;

  char chemin[1024];
      snprintf(chemin, sizeof(chemin), "./Ensemble_echantillons/Rfg_%02d_Ndf_%02d_eta_%02d_Lb_%02d/Echantillon_%d", 
         (int)(R_fg * 10), (int)Nd_f, (int)(eta * 10), (int)Lb, echantillon_i);

  // Construire et exécuter les commandes
  char command1[256];
  snprintf(command1, sizeof(command1), "mkdir -p \"%s\"", chemin);
  //std::cout << "Commande 1: " << command1 << endl;

  char chemin_simu[1024];
  snprintf(chemin_simu, sizeof(chemin), "./Ensemble_echantillons/Rfg_%02d_Ndf_%02d_eta_%02d_Lb_%02d/Echantillon_%d/Simulation", (int)(R_fg * 10), (int)Nd_f, (int)(eta * 10), (int)Lb, echantillon_i);

  // Construire et exécuter les commandes
  char command2[256];
  snprintf(command2, sizeof(command2), "mkdir -p \"%s\"", chemin_simu);
  //std::cout << "Commande 1: " << command1 << endl;


  // Exécuter les commandes
  if (system(command1) == 0)
    {
      //std::cout << "Dossier créé avec succès : " << chemin << std::endl;
    }
  else
    {
      std::cerr << "Erreur lors de la création du dossier : " << chemin << std::endl;
    }

    output_grains("grains.vtk", chemin);
    output_fibres("fibres.vtk", chemin);
    output_fibres_down(echantillon_i);
    output_fibres_down_right(echantillon_i);
    output_fibres_down_left(echantillon_i);
    output_fibres_up(echantillon_i);
    output_fibres_up_right(echantillon_i);
    output_fibres_up_left(echantillon_i);
    output_fibres_right(echantillon_i);
    output_fibres_left(echantillon_i);
    output_fibres_grains("fibres_grains.txt", chemin);


  float sum_R=0;
  double  S_g=0;
  for (i=0; i<vgrains.size(); i++)
    {
      sum_R += vgrains[i].R;
      S_g += PI*vgrains[i].R*vgrains[i].R;
    }

  double R_moy = sum_R/vgrains.size();
  double R_moy_th = (R_min+R_max)/2;
  double S_f = vfibres.size()*S_1f;

  double R_fg_calc = S_f/(S_f+S_g);
  double C_fg_calc = (S_f+S_g)/(S_b);

  cout << "Number of fiber insertion attemps : " << att_f << " " << "Number of grain insertion attempts : " << att_g << endl;
  cout << "Surface ratio of the generated system : " << R_fg_calc << " ; Compacity of the generated system : " << C_fg_calc << endl;
  cout <<  "Average radius of the grains : " << R_moy  << endl;

  // Nom du fichier de sortie
  std::string filename_1 = "calculated_values.dat";

  string cheminStr(chemin);
  string cheminComplet_1 = cheminStr + "/" + filename_1;

  std::ofstream outfile;

  outfile.open(cheminComplet_1);
  if (outfile.is_open())
    {
      // Ecrire les valeurs calculées dans le fichier
      outfile << "Rfg Ns eta Lb nf ng nb_essais_fibres nb_essais_grains nb_essais_total Moy_rayons_grains_calculee Rfg_calculee Cfg_calculee" << endl;
      outfile << R_fg << " " << Nd_f << " " << eta << " " << std::stoi(argv[4]) << " " << nf << " " << ng << " " << att_f << " " << att_g << " " << att_f+att_g << " " << R_moy << " " << R_fg_calc << " " << C_fg_calc << endl;

      // Fermer le fichier
      outfile.close();
    }
  else
    {
      std::cerr << "Impossible d'ouvrir le fichier " << filename_1 << " pour écriture." << std::endl;
    }

  std::string filename_2 = "graine.dat";

  string cheminComplet_2 = cheminStr + "/" + filename_2;


  outfile.open(cheminComplet_2);
  if (outfile.is_open())
    {
      // Ecrire les valeurs calculées dans le fichier
      //outfile << R_moy_th << " " << R_fg << " " << C_fg << endl;
      outfile << seed << endl;

      // Fermer le fichier
      outfile.close();
    }
  else
    {
      std::cerr << "Impossible d'ouvrir le fichier " << filename_2 << " pour écriture." << std::endl;
    }


  // // Construction des chemins
  // char chemin1[1024];
  // snprintf(chemin1, sizeof(chemin1), "/home/pierre/Bureau/Programmes/Ensemble_echantillons/Rfg_%.2f_Ndf_%d_eta_%.2f_Lb_%d", R_fg, Nd_f, eta, Lb);

  // char chemin2[1024];
  // snprintf(chemin2, sizeof(chemin2), "%s/Echantillon_%d", chemin1, echantillon_i);

  // // // Afficher les chemins pour vérifier
  // // std::cout << "Chemin 1: " << chemin1 << std::endl;
  // // std::cout << "Chemin 2: " << chemin2 << std::endl;


  // // Construire et exécuter les commandes
  // char command1[256];
  // snprintf(command1, sizeof(command1), "mkdir -p \"%s\"", chemin1);
  // //std::cout << "Commande 1: " << command1 << endl;

  // char command2[256];
  // snprintf(command2, sizeof(command2), "mkdir -p \"%s\"", chemin2);
  // //std::cout << "Commande 2: " << command2 << std::endl;



  // // Exécuter les commandes
  // if (system(command1) == 0) {
  //   //std::cout << "Dossier créé avec succès : " << chemin1 << std::endl;
  // } else {
  //   std::cerr << "Erreur lors de la création du dossier : " << chemin1 << std::endl;
  // }

  // if (system(command2) == 0) {
  //   //std::cout << "Dossier créé avec succès : " << chemin2 << std::endl;
  // } else {
  //   std::cerr << "Erreur lors de la création du dossier : " << chemin2 << std::endl;
  // }


  // // fichiers à déplacer
  // const char* fichierSource1 = "/home/pierre/Bureau/Programmes/fibres0001.vtk";
  // char fichierDestination1[1024];
  // snprintf(fichierDestination1, sizeof(fichierDestination1), "%s/fibres0001.vtk", chemin2);


  // try
  //   {
  //     // Vérifier si le dossier de destination existe
  //     if (!fs::exists(chemin2))
  // 	{
  // 	  std::cerr << "Erreur : le dossier de destination n'existe pas." << std::endl;
  // 	  return 1;
  // 	}

  //     // Déplacer le fichier
  //     fs::rename(fichierSource1, fichierDestination1);
  //     //std::cout << "Fichier déplacé avec succès vers : " << fichierDestination1 << std::endl;
  //   }
  // catch (const fs::filesystem_error& e)
  //   {
  //     std::cerr << "Erreur lors du déplacement du fichier : " << e.what() << std::endl;
  //   }


  // const char* fichierSource2 = "/home/pierre/Bureau/Programmes/grains0001.vtk";
  // char fichierDestination2[1024];
  // snprintf(fichierDestination2, sizeof(fichierDestination2), "%s/grains0001.vtk", chemin2);


  // try
  //   {
  //     // Vérifier si le dossier de destination existe
  //     if (!fs::exists(chemin2))
  // 	{
  // 	  std::cerr << "Erreur : le dossier de destination n'existe pas." << std::endl;
  // 	  return 1;
  // 	}

  //     // Déplacer le fichier
  //     fs::rename(fichierSource2, fichierDestination2);
  //     //std::cout << "Fichier déplacé avec succès vers : " << fichierDestination2 << std::endl;
  //   }
  // catch (const fs::filesystem_error& e)
  //   {
  //     std::cerr << "Erreur lors du déplacement du fichier : " << e.what() << std::endl;
  //   }


  // const char* fichierSource3 = "/home/pierre/Bureau/Programmes/fibres_grains0001.txt";
  // char fichierDestination3[1024];
  // snprintf(fichierDestination3, sizeof(fichierDestination3), "%s/fibres_grains0001.txt", chemin2);


  // try
  //   {
  //     // Vérifier si le dossier de destination existe
  //     if (!fs::exists(chemin2))
  // 	{
  // 	  std::cerr << "Erreur : le dossier de destination n'existe pas." << std::endl;
  // 	  return 1;
  // 	}

  //     // Déplacer le fichier
  //     fs::rename(fichierSource3, fichierDestination3);
  //     //std::cout << "Fichier déplacé avec succès vers : " << fichierDestination3 << std::endl;
  //   }
  // catch (const fs::filesystem_error& e)
  //   {
  //     std::cerr << "Erreur lors du déplacement du fichier : " << e.what() << std::endl;
  //   }


  // const char* fichierSource4 = "/home/pierre/Bureau/Programmes/graine.dat";
  // char fichierDestination4[1024];
  // snprintf(fichierDestination4, sizeof(fichierDestination4), "%s/graine.dat", chemin2);


  // try
  //   {
  //     // Vérifier si le dossier de destination existe
  //     if (!fs::exists(chemin2))
  // 	{
  // 	  std::cerr << "Erreur : le dossier de destination n'existe pas." << std::endl;
  // 	  return 1;
  // 	}

  //     // Déplacer le fichier
  //     fs::rename(fichierSource4, fichierDestination4);
  //     //std::cout << "Fichier déplacé avec succès vers : " << fichierDestination4 << std::endl;
  //   }
  // catch (const fs::filesystem_error& e)
  //   {
  //     std::cerr << "Erreur lors du déplacement du fichier : " << e.what() << std::endl;
  //   }

  return 0;

}
////////////////////
//  Functions
///////////////////

// Save grains vtk file
void output_grains(const string& nomFichier, const string& cheminDestination)
{
  int i;

  ofstream fog(cheminDestination + "/" + nomFichier);
  if (fog)
  {
    fog.precision(5); fog << scientific;
    fog << "# vtk DataFile Version 3.0" << endl;
    fog << "My grains" << endl;
    fog << "ASCII" << endl;
    fog << "DATASET UNSTRUCTURED_GRID" << endl;
    fog << "POINTS " << vgrains.size() << " float" << endl;
    for(i=0;i<vgrains.size();i++) fog << vgrains[i].x <<" "<< vgrains[i].y << " 0.0" << endl;
    fog << "POINT_DATA " << ng << endl;
    fog << "VECTORS Radius float" << endl;
    for(i=0;i<vgrains.size();i++) fog << vgrains[i].R << " 0.0 0.0" << endl;
  }
}

// // Save grains vtk file in the periodic boxes
// void output_grains_right(int numfile)
// {
//     int i;
//     char fname[25]; // file name
//     sprintf(fname,"grains_right%04d.vtk",numfile);

//     ofstream fog(fname, ios::out);
// 	if (fog)
// 	{
//         fog.precision(5); fog << scientific;
//         fog << "# vtk DataFile Version 3.0" << endl;
//         fog << "My grains" << endl;
//         fog << "ASCII" << endl;
//         fog << "DATASET UNSTRUCTURED_GRID" << endl;
//         fog << "POINTS " << vgrains.size() << " float" << endl;
//         for(i=0;i<vgrains.size();i++) fog << vgrains[i].x + b_xmax <<" "<< vgrains[i].y << " 0.0" << endl;
//         fog << "POINT_DATA " << ng << endl;
//         fog << "VECTORS Radius float" << endl;
//         for(i=0;i<vgrains.size();i++) fog << vgrains[i].R << " 0.0 0.0" << endl;
// 	}
// }

// // Save grains vtk file
// void output_grains_left(int numfile)
// {
//     int i;
//     char fname[25]; // file name
//     sprintf(fname,"grains_left%04d.vtk",numfile);

//     ofstream fog(fname, ios::out);
// 	if (fog)
// 	{
//         fog.precision(5); fog << scientific;
//         fog << "# vtk DataFile Version 3.0" << endl;
//         fog << "My grains" << endl;
//         fog << "ASCII" << endl;
//         fog << "DATASET UNSTRUCTURED_GRID" << endl;
//         fog << "POINTS " << vgrains.size() << " float" << endl;
//         for(i=0;i<vgrains.size();i++) fog << vgrains[i].x - b_xmax <<" "<< vgrains[i].y << " 0.0" << endl;
//         fog << "POINT_DATA " << ng << endl;
//         fog << "VECTORS Radius float" << endl;
//         for(i=0;i<vgrains.size();i++) fog << vgrains[i].R << " 0.0 0.0" << endl;
// 	}
// }

// // Save grains vtk file
// void output_grains_up(int numfile)
// {
//     int i;
//     char fname[25]; // file name
//     sprintf(fname,"grains_up%04d.vtk",numfile);

//     ofstream fog(fname, ios::out);
// 	if (fog)
// 	{
//         fog.precision(5); fog << scientific;
//         fog << "# vtk DataFile Version 3.0" << endl;
//         fog << "My grains" << endl;
//         fog << "ASCII" << endl;
//         fog << "DATASET UNSTRUCTURED_GRID" << endl;
//         fog << "POINTS " << vgrains.size() << " float" << endl;
//         for(i=0;i<vgrains.size();i++) fog << vgrains[i].x <<" "<< vgrains[i].y + b_ymax << " 0.0" << endl;
//         fog << "POINT_DATA " << ng << endl;
//         fog << "VECTORS Radius float" << endl;
//         for(i=0;i<vgrains.size();i++) fog << vgrains[i].R << " 0.0 0.0" << endl;
// 	}
// }

// // Save grains vtk file
// void output_grains_down(int numfile)
// {
//     int i;
//     char fname[25]; // file name
//     sprintf(fname,"grains_down%04d.vtk",numfile);

//     ofstream fog(fname, ios::out);
// 	if (fog)
// 	{
//         fog.precision(5); fog << scientific;
//         fog << "# vtk DataFile Version 3.0" << endl;
//         fog << "My grains" << endl;
//         fog << "ASCII" << endl;
//         fog << "DATASET UNSTRUCTURED_GRID" << endl;
//         fog << "POINTS " << vgrains.size() << " float" << endl;
//         for(i=0;i<vgrains.size();i++) fog << vgrains[i].x <<" "<< vgrains[i].y - b_ymax << " 0.0" << endl;
//         fog << "POINT_DATA " << ng << endl;
//         fog << "VECTORS Radius float" << endl;
//         for(i=0;i<vgrains.size();i++) fog << vgrains[i].R << " 0.0 0.0" << endl;
// 	}
// }

// // Save grains vtk file
// void output_grains_down_right(int numfile)
// {
//     int i;
//     char fname[30]; // file name
//     sprintf(fname,"grains_down_r%04d.vtk",numfile);

//     ofstream fog(fname, ios::out);
// 	if (fog)
// 	{
//         fog.precision(5); fog << scientific;
//         fog << "# vtk DataFile Version 3.0" << endl;
//         fog << "My grains" << endl;
//         fog << "ASCII" << endl;
//         fog << "DATASET UNSTRUCTURED_GRID" << endl;
//         fog << "POINTS " << vgrains.size() << " float" << endl;
//         for(i=0;i<vgrains.size();i++) fog << vgrains[i].x + b_xmax <<" "<< vgrains[i].y - b_ymax << " 0.0" << endl;
//         fog << "POINT_DATA " << ng << endl;
//         fog << "VECTORS Radius float" << endl;
//         for(i=0;i<vgrains.size();i++) fog << vgrains[i].R << " 0.0 0.0" << endl;
// 	}
// }

// // Save grains vtk file
// void output_grains_down_left(int numfile)
// {
//     int i;
//     char fname[25]; // file name
//     sprintf(fname,"grains_down_left%04d.vtk",numfile);

//     ofstream fog(fname, ios::out);
// 	if (fog)
// 	{
//         fog.precision(5); fog << scientific;
//         fog << "# vtk DataFile Version 3.0" << endl;
//         fog << "My grains" << endl;
//         fog << "ASCII" << endl;
//         fog << "DATASET UNSTRUCTURED_GRID" << endl;
//         fog << "POINTS " << vgrains.size() << " float" << endl;
//         for(i=0;i<vgrains.size();i++) fog << vgrains[i].x - b_xmax <<" "<< vgrains[i].y - b_ymax << " 0.0" << endl;
//         fog << "POINT_DATA " << ng << endl;
//         fog << "VECTORS Radius float" << endl;
//         for(i=0;i<vgrains.size();i++) fog << vgrains[i].R << " 0.0 0.0" << endl;
// 	}
// }

// // Save grains vtk file
// void output_grains_up_left(int numfile)
// {
//     int i;
//     char fname[25]; // file name
//     sprintf(fname,"grains_up_left%04d.vtk",numfile);

//     ofstream fog(fname, ios::out);
// 	if (fog)
// 	{
//         fog.precision(5); fog << scientific;
//         fog << "# vtk DataFile Version 3.0" << endl;
//         fog << "My grains" << endl;
//         fog << "ASCII" << endl;
//         fog << "DATASET UNSTRUCTURED_GRID" << endl;
//         fog << "POINTS " << vgrains.size() << " float" << endl;
//         for(i=0;i<vgrains.size();i++) fog << vgrains[i].x - b_xmax <<" "<< vgrains[i].y + b_ymax << " 0.0" << endl;
//         fog << "POINT_DATA " << ng << endl;
//         fog << "VECTORS Radius float" << endl;
//         for(i=0;i<vgrains.size();i++) fog << vgrains[i].R << " 0.0 0.0" << endl;
// 	}
// }

// // Save grains vtk file
// void output_grains_up_right(int numfile)
// {
//     int i;
//     char fname[25]; // file name
//     sprintf(fname,"grains_up_right%04d.vtk",numfile);

//     ofstream fog(fname, ios::out);
// 	if (fog)
// 	{
//         fog.precision(5); fog << scientific;
//         fog << "# vtk DataFile Version 3.0" << endl;
//         fog << "My grains" << endl;
//         fog << "ASCII" << endl;
//         fog << "DATASET UNSTRUCTURED_GRID" << endl;
//         fog << "POINTS " << vgrains.size() << " float" << endl;
//         for(i=0;i<vgrains.size();i++) fog << vgrains[i].x + b_xmax <<" "<< vgrains[i].y + b_ymax << " 0.0" << endl;
//         fog << "POINT_DATA " << ng << endl;
//         fog << "VECTORS Radius float" << endl;
//         for(i=0;i<vgrains.size();i++) fog << vgrains[i].R << " 0.0 0.0" << endl;
// 	}
// }


// function that tests if two segments cross each other.
bool doIntersect(double x1, double y1, double x2, double y2,
  double x3, double y3, double x4, double y4, float o1, float o2) {


    double c_d1 = (y2-y1)/(x2-x1); // coefficient directeur du segment 1 (la fibre candidate)
    double b1 = y1 - c_d1*x1;
    double c_d2 = (y4-y3)/(x4-x3); // coefficient directeur du segment 2 (la fibre i)
    double b2 = y3 - c_d2*x3;

    double x12 = x2-x1;
    double x34 = x4-x3;
    double abs_x12, abs_x34;

    if (x12<0)
    {
      abs_x12 = -x12;
    }
    else
    {
      abs_x12 = x12;
    }

    if (x34<0)
    {
      abs_x34 = -x34;
    }
    else
    {
      abs_x34 = x34;
    }


    if (c_d1==c_d2) // segements are parallel, they do not cross
    {
      return false;
    }

    else // tous les autres cas
    {
      double x_inter = (b2-b1)/(c_d1-c_d2); // calcule la coordonée en x du point d'intersection

      // Dans la suite, on teste si le point d'intersection est situé au niveau des segments où s'il est en en dehors, dans leur prolongement.
      // Pour ce faire, on détermine si x_inter se situe entre les coordonées en x des extrémités des deux segments.
      // Il y a donc plusieurs cas à tester selon l'orientation des fibres

      // Note : je peux réduire la fonction en utilisant la valeur absolue...
      if ((o1>=0 && o1<PI/2) || (o1>3*PI/2 && o1<=2*PI))
      {
        if ((o2>=0 && o2<PI/2) || (o2>3*PI/2 && o2<=2*PI))
        {
          if ((x_inter>=x1 && x_inter<=x2) && (x_inter>=x3 && x_inter<=x4))
          {
            return true; // intersectent
          }
          else
          {
            return false;
          }
        }
        if (o2>PI/2 && o2<3*PI/2)
        {
          if ((x_inter>=x1 && x_inter<=x2) && (x_inter>=x4 && x_inter<=x3))
          {
            return true;
          }
          else
          {
            return false;
          }
        }
        else
        {
          return false;
        }
      }

      if (o1>PI/2 && o1<3*PI/2)
      {
        if ((o2>=0 && o2<PI/2) || (o2>3*PI/2 && o2<=2*PI))
        {
          if ((x_inter>=x2 && x_inter<=x1) && (x_inter>=x3 && x_inter<=x4))
          {
            return true;
          }
          else
          {
            return false;
          }
        }
        if (o2>PI/2 && o2<3*PI/2)
        {
          if ((x_inter>=x2 && x_inter<=x1) && (x_inter>=x4 && x_inter<=x3))
          {
            return true;
          }
          else
          {
            return false;
          }
        }
        else
        {
          return false;
        }
      }
      else
      {
        return false;
      }
    }
  }

// Save fibres vtk file
void output_fibres(const string& nomFichier, const string& cheminDestination)
{
  int i;
  // char fname[25]; // file name
  // sprintf(fname,"fibres%04d.vtk",numfile);

  ofstream fog(cheminDestination + "/" + nomFichier);
  if (fog)
  {
    fog.precision(5); fog << scientific;
    fog << "# vtk DataFile Version 3.0" << endl;
    fog << "My fibers" << endl;
    fog << "ASCII" << endl;
    fog << "DATASET POLYDATA" << endl;
    fog << "POINTS " << 2*(nth+1)*nf + 4 << " float" << endl;
    // coordinates of the simulation box
    fog << "0.0 0.0 0.0" << endl;
    fog << b_xmax << " 0.0 0.0" << endl;
    fog << b_xmax << " " << b_ymax << " 0.0" << endl;
    fog << 0.0 << " " << b_ymax << " 0.0" << endl;
    for(i=0;i<nf;i++)
    {
      for (int j=0; j<nth+1; j++)
      {
        fog << vp_f[i].pos_1[2*j] << " " << vp_f[i].pos_1[2*(j+0.5)] << " 0.0" << endl;
      }
      for (int k=0; k<nth+1; k++)
      {
        fog << vp_f[i].pos_2[2*k] << " " << vp_f[i].pos_2[2*(k+0.5)] << " 0.0" << endl;
      }
    }
    fog << "POLYGONS " << nf + 1 << " " << 2*(nth+1)*nf + 4 + (nf+1) << endl;
    fog << "4 0 1 2 3";
    fog << endl;
    for (i=0; i<nf; i++)
    {
      fog << 2*(nth+1) << " ";
      for (int l=0; l<2*(nth+1); l++) fog << i*(nth+1)*2 + l + 4 << " ";
      fog << endl;
    }
  }
}

// Save fibres vtk file (periodic boxes)
void output_fibres_right(int numfile)
{
  int i;
  char fname[25]; // file name
  sprintf(fname,"fibres_right%04d.vtk",numfile);

  ofstream fog(fname, ios::out);
  if (fog)
  {
    fog.precision(5); fog << scientific;
    fog << "# vtk DataFile Version 3.0" << endl;
    fog << "My fibers" << endl;
    fog << "ASCII" << endl;
    fog << "DATASET POLYDATA" << endl;
    fog << "POINTS " << 2*(nth+1)*nf << " float" << endl;
    for(i=0;i<nf;i++)
    {
      for (int j=0; j<nth+1; j++)
      {
        fog << vp_f[i].pos_1[2*j] + b_xmax << " " << vp_f[i].pos_1[2*(j+0.5)] << " 0.0" << endl;
      }
      for (int k=0; k<nth+1; k++)
      {
        fog << vp_f[i].pos_2[2*k] + b_xmax << " " << vp_f[i].pos_2[2*(k+0.5)] << " 0.0" << endl;
      }
    }
    fog << "POLYGONS " << nf << " " << 2*(nth+1)*nf + nf << endl;
    for (i=0; i<nf; i++)
    {
      fog << 2*(nth+1) << " ";
      for (int l=0; l<2*(nth+1); l++) fog << i*(nth+1)*2 + l << " ";
      fog << endl;
    }
  }
}

// Save fibres vtk file
void output_fibres_left(int numfile)
{
  int i;
  char fname[25]; // file name
  sprintf(fname,"fibres_left%04d.vtk",numfile);

  ofstream fog(fname, ios::out);
  if (fog)
  {
    fog.precision(5); fog << scientific;
    fog << "# vtk DataFile Version 3.0" << endl;
    fog << "My fibers" << endl;
    fog << "ASCII" << endl;
    fog << "DATASET POLYDATA" << endl;
    fog << "POINTS " << 2*(nth+1)*nf << " float" << endl;
    for(i=0;i<nf;i++)
    {
      for (int j=0; j<nth+1; j++)
      {
        fog << vp_f[i].pos_1[2*j] - b_xmax << " " << vp_f[i].pos_1[2*(j+0.5)] << " 0.0" << endl;
      }
      for (int k=0; k<nth+1; k++)
      {
        fog << vp_f[i].pos_2[2*k] - b_xmax << " " << vp_f[i].pos_2[2*(k+0.5)] << " 0.0" << endl;
      }
    }
    fog << "POLYGONS " << nf << " " << 2*(nth+1)*nf + nf << endl;
    for (i=0; i<nf; i++)
    {
      fog << 2*(nth+1) << " ";
      for (int l=0; l<2*(nth+1); l++) fog << i*(nth+1)*2 + l << " ";
      fog << endl;
    }
  }
}

// Save fibres vtk file
void output_fibres_up(int numfile)
{
  int i;
  char fname[25]; // file name
  sprintf(fname,"fibres_up%04d.vtk",numfile);

  ofstream fog(fname, ios::out);
  if (fog)
  {
    fog.precision(5); fog << scientific;
    fog << "# vtk DataFile Version 3.0" << endl;
    fog << "My fibers" << endl;
    fog << "ASCII" << endl;
    fog << "DATASET POLYDATA" << endl;
    fog << "POINTS " << 2*(nth+1)*nf << " float" << endl;
    for(i=0;i<nf;i++)
    {
      for (int j=0; j<nth+1; j++)
      {
        fog << vp_f[i].pos_1[2*j] << " " << vp_f[i].pos_1[2*(j+0.5)] + b_ymax << " 0.0" << endl;
      }
      for (int k=0; k<nth+1; k++)
      {
        fog << vp_f[i].pos_2[2*k] << " " << vp_f[i].pos_2[2*(k+0.5)] + b_ymax << " 0.0" << endl;
      }
    }
    fog << "POLYGONS " << nf << " " << 2*(nth+1)*nf + nf << endl;
    for (i=0; i<nf; i++)
    {
      fog << 2*(nth+1) << " ";
      for (int l=0; l<2*(nth+1); l++) fog << i*(nth+1)*2 + l << " ";
      fog << endl;
    }
  }
}

// Save fibres vtk file
void output_fibres_down(int numfile)
{
  int i;
  char fname[25]; // file name
  sprintf(fname,"fibres_down%04d.vtk",numfile);

  ofstream fog(fname, ios::out);
  if (fog)
  {
    fog.precision(5); fog << scientific;
    fog << "# vtk DataFile Version 3.0" << endl;
    fog << "My fibers" << endl;
    fog << "ASCII" << endl;
    fog << "DATASET POLYDATA" << endl;
    fog << "POINTS " << 2*(nth+1)*nf << " float" << endl;
    for(i=0;i<nf;i++)
    {
      for (int j=0; j<nth+1; j++)
      {
        fog << vp_f[i].pos_1[2*j] << " " << vp_f[i].pos_1[2*(j+0.5)] - b_ymax << " 0.0" << endl;
      }
      for (int k=0; k<nth+1; k++)
      {
        fog << vp_f[i].pos_2[2*k] << " " << vp_f[i].pos_2[2*(k+0.5)] - b_ymax << " 0.0" << endl;
      }
    }
    fog << "POLYGONS " << nf << " " << 2*(nth+1)*nf + nf << endl;
    for (i=0; i<nf; i++)
    {
      fog << 2*(nth+1) << " ";
      for (int l=0; l<2*(nth+1); l++) fog << i*(nth+1)*2 + l << " ";
      fog << endl;
    }
  }
}

void output_fibres_down_right(int numfile)
{
  int i;
  char fname[30]; // file name
  sprintf(fname,"fibres_down_right%04d.vtk",numfile);

  ofstream fog(fname, ios::out);
  if (fog)
  {
    fog.precision(5); fog << scientific;
    fog << "# vtk DataFile Version 3.0" << endl;
    fog << "My fibers" << endl;
    fog << "ASCII" << endl;
    fog << "DATASET POLYDATA" << endl;
    fog << "POINTS " << 2*(nth+1)*nf << " float" << endl;
    for(i=0;i<nf;i++)
    {
      for (int j=0; j<nth+1; j++)
      {
        fog << vp_f[i].pos_1[2*j] + b_xmax << " " << vp_f[i].pos_1[2*(j+0.5)] - b_ymax << " 0.0" << endl;
      }
      for (int k=0; k<nth+1; k++)
      {
        fog << vp_f[i].pos_2[2*k] + b_xmax << " " << vp_f[i].pos_2[2*(k+0.5)] - b_ymax << " 0.0" << endl;
      }
    }
    fog << "POLYGONS " << nf << " " << 2*(nth+1)*nf + nf << endl;
    for (i=0; i<nf; i++)
    {
      fog << 2*(nth+1) << " ";
      for (int l=0; l<2*(nth+1); l++) fog << i*(nth+1)*2 + l << " ";
      fog << endl;
    }
  }
}

void output_fibres_down_left(int numfile)
{
  int i;
  char fname[25]; // file name
  sprintf(fname,"fibres_down_left%04d.vtk",numfile);

  ofstream fog(fname, ios::out);
  if (fog)
  {
    fog.precision(5); fog << scientific;
    fog << "# vtk DataFile Version 3.0" << endl;
    fog << "My fibers" << endl;
    fog << "ASCII" << endl;
    fog << "DATASET POLYDATA" << endl;
    fog << "POINTS " << 2*(nth+1)*nf << " float" << endl;
    for(i=0;i<nf;i++)
    {
      for (int j=0; j<nth+1; j++)
      {
        fog << vp_f[i].pos_1[2*j] - b_xmax << " " << vp_f[i].pos_1[2*(j+0.5)] - b_ymax << " 0.0" << endl;
      }
      for (int k=0; k<nth+1; k++)
      {
        fog << vp_f[i].pos_2[2*k] - b_xmax << " " << vp_f[i].pos_2[2*(k+0.5)] - b_ymax << " 0.0" << endl;
      }
    }
    fog << "POLYGONS " << nf << " " << 2*(nth+1)*nf + nf << endl;
    for (i=0; i<nf; i++)
    {
      fog << 2*(nth+1) << " ";
      for (int l=0; l<2*(nth+1); l++) fog << i*(nth+1)*2 + l << " ";
      fog << endl;
    }
  }
}

void output_fibres_up_right(int numfile)
{
  int i;
  char fname[25]; // file name
  sprintf(fname,"fibres_up_right%04d.vtk",numfile);

  ofstream fog(fname, ios::out);
  if (fog)
  {
    fog.precision(5); fog << scientific;
    fog << "# vtk DataFile Version 3.0" << endl;
    fog << "My fibers" << endl;
    fog << "ASCII" << endl;
    fog << "DATASET POLYDATA" << endl;
    fog << "POINTS " << 2*(nth+1)*nf << " float" << endl;
    for(i=0;i<nf;i++)
    {
      for (int j=0; j<nth+1; j++)
      {
        fog << vp_f[i].pos_1[2*j] + b_xmax << " " << vp_f[i].pos_1[2*(j+0.5)] + b_ymax << " 0.0" << endl;
      }
      for (int k=0; k<nth+1; k++)
      {
        fog << vp_f[i].pos_2[2*k] + b_xmax << " " << vp_f[i].pos_2[2*(k+0.5)] + b_ymax << " 0.0" << endl;
      }
    }
    fog << "POLYGONS " << nf << " " << 2*(nth+1)*nf + nf << endl;
    for (i=0; i<nf; i++)
    {
      fog << 2*(nth+1) << " ";
      for (int l=0; l<2*(nth+1); l++) fog << i*(nth+1)*2 + l << " ";
      fog << endl;
    }
  }
}

void output_fibres_up_left(int numfile)
{
  int i;
  char fname[25]; // file name
  sprintf(fname,"fibres_up_left%04d.vtk",numfile);

  ofstream fog(fname, ios::out);
  if (fog)
  {
    fog.precision(5); fog << scientific;
    fog << "# vtk DataFile Version 3.0" << endl;
    fog << "My fibers" << endl;
    fog << "ASCII" << endl;
    fog << "DATASET POLYDATA" << endl;
    fog << "POINTS " << 2*(nth+1)*nf << " float" << endl;
    for(i=0;i<nf;i++)
    {
      for (int j=0; j<nth+1; j++)
      {
        fog << vp_f[i].pos_1[2*j] - b_xmax << " " << vp_f[i].pos_1[2*(j+0.5)] + b_ymax << " 0.0" << endl;
      }
      for (int k=0; k<nth+1; k++)
      {
        fog << vp_f[i].pos_2[2*k] - b_xmax << " " << vp_f[i].pos_2[2*(k+0.5)] + b_ymax << " 0.0" << endl;
      }
    }
    fog << "POLYGONS " << nf << " " << 2*(nth+1)*nf + nf << endl;
    for (i=0; i<nf; i++)
    {
      fog << 2*(nth+1) << " ";
      for (int l=0; l<2*(nth+1); l++) fog << i*(nth+1)*2 + l << " ";
      fog << endl;
    }
  }
}

// Save fibres et grains file
void output_fibres_grains(const string& nomFichier, const string& cheminDestination)
{
  int i;
  // char fname[25]; // file name
  // sprintf(fname,"fibres_grains%04d.txt",numfile);
  char fname2[1024];
  sprintf(fname2, "./Ensemble_echantillons/Rfg_%02d_Ndf_%02d_eta_%02d_Lb_%02d/Echantillon_%d/Simulation", (int)(R_fg * 10), (int)Nd_f, (int)(eta * 10), (int)Lb, echantillon_i);


  ofstream fog(cheminDestination + "/" + nomFichier);
  if (fog)
  {
    fog.precision(5); fog << scientific;
    fog << "ClaySOMEFv2 09-02-2015" << endl;
    fog << "result_folder " << fname2 << endl;
    //fog << "result_folder ./Test/sous_test" << endl;
    fog << "tmax 60" << endl;
    fog << "dt 2e-6" << endl;
    fog << "interVtk 1" << endl;
    fog << "interVerlet 0.01" << endl;
    fog << "interOut 0.01" << endl;
    fog << "interConf 1" << endl;
    fog << "dVerlet 1e-3" << endl;
    fog << "autoVerlet 1" << endl;
    fog << "dcut 4" << endl;
    fog << "kn 1000" << endl;
    fog << "kt 1000" << endl;
    fog << "kb 1000" << endl;
    fog << "kc_ff_n 1e6" << endl;
    fog << "kc_ff_t 1e6" << endl;
    fog << "kc_fg_n 1e6" << endl;
    fog << "kc_fg_t 1e6" << endl;
    fog << "kc_gg_n 1e6" << endl;
    fog << "kc_gg_t 1e6" << endl;
    fog << "nu_ff 0.0" << endl;
    fog << "nu_fg 0.0" << endl;
    fog << "nu_gg 0.0" << endl;
    fog << "density 2600" << endl;
    fog << "alpha 0.1" << endl;
    fog << "Cell_Damping 0.0" <<endl;
    fog << "ivtk 1" << endl;
    fog << "iconf 1" << endl;
    fog << "Cell " << Lb*L_max << " 0 0 " << Lb*L_max << endl;
    fog << "vCell 0 0 0 0" << endl;
    fog << "aCell 0 0 0 0" << endl;
    fog << "ImageCell 1" << endl;
    fog << "Load IsostaticCompression 3" << endl;
    fog << "fibres " << vfibres.size() << " " << Nd_f << endl;
    for(i=0;i<vfibres.size();i++)
    {
      fog << vfibres[i].pos_c_1[0] << " " << vfibres[i].pos_c_1[1] << " " << vfibres[i].pos_c_2[0] << " " << vfibres[i].pos_c_2[1] << " " << vfibres[i].R << endl;
    }
    fog << "grains " << vgrains.size() << endl;
    for (i=0; i<vgrains.size(); i++)
    {
      fog << vgrains[i].x  << " " << vgrains[i].y << " " << vgrains[i].R << endl;
    }
  }
}

// // Save data
// void output_data(int numfile)
// {
//     int i;
//     char fname[25]; // file name
//     sprintf(fname,"data.txt",numfile);

//     ofstream fog(fname, ios::out);
// 	if (fog)
// 	{
//         fog.precision(5); fog << scientific;
// 	fog << "fibres " << vfibres.size() << endl;
//         for(i=0;i<nf;i++)
// 	  {
// 	    fog << vfibres[i].pos_c_1[0] << " " << vfibres[i].pos_c_1[1] << " " << vfibres[i].pos_c_2[0] << " " << vfibres[i].pos_c_2[1] << " " << vfibres[i].R << endl;
// 	  }
//         fog << "Grains " << vgrains.size() << endl;
// 	for (i=0; i<vgrains.size(); i++)
// 	  {
// 	    fog << vgrains[i].x  << " " << vgrains[i].y << " " << vgrains[i].R << endl;
// 	  }
// 	}
// }

// void output_postscript(int numfile)
// {
// int i;
//     char fname[25]; // file name
//     sprintf(fname,"fibres%04d.ps",numfile);

//     ofstream fog(fname, ios::out);
// 	if (fog)
// 	{
//         fog.precision(5); fog << scientific;
//         fog << "# entête poscript" << endl;
//         for(i=0;i<nf;i++)
// 	  {
// 	    for (int j=0; j<nth+1; j++)
// 	      {
// 		fog << vp_f[i].pos_1[2*j] << " " << vp_f[i].pos_1[2*(j+0.5)] << " 0.0" << endl;
// 	      }
// 	    for (int k=0; k<nth+1; k++)
// 	      {
// 		fog << vp_f[i].pos_2[2*k] << " " << vp_f[i].pos_2[2*(k+0.5)] << " 0.0" << endl;
// 	      }
// 	  }
//         fog << "POLYGONS " << nf << " " << 2*(nth+1)*nf + nf << endl;
// 	for (i=0; i<nf; i++)
// 	  {
// 	    fog << 2*(nth+1) << " ";
// 	    for (int l=0; l<2*(nth+1); l++) fog << i*(nth+1)*2 + l << " ";
// 	    fog << endl;
// 	  }
// 	}
// }
