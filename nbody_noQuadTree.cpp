#include <math.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <iostream>

using namespace std;

static double grav_const = 6.67430;
static double repulsion = 0.0000005;
static double theta = 0.5;
static double L = 1;

// Define a particle in 2D with a certain mass,
// a position (xpos,ypos) for each coordinate and
// two velocity components (xvel,yvel)
struct Particle {
  double xpos;
  double ypos;
  double xvel;
  double yvel;
  double mass;

  Particle() {
    xpos = 0.0;
    ypos = 0.0;
    xvel = 0.0;
    yvel = 0.0;
    mass = 0.0;
  }

  Particle(double _xpos, double _ypos, double _xvel, double _yvel,
           double _mass) {
    xpos = _xpos;
    ypos = _ypos;
    xvel = _xvel;
    yvel = _yvel;
    mass = _mass;
  }
};

struct vector2D {
  double x;
  double y;
  vector2D (double _x = 0, double _y = 0){
    x = _x;
    y = _y;
  }

  vector2D operator+(const vector2D& other) const {
    return vector2D(x + other.x, y + other.y);
  }

  vector2D operator*(const int other) const {
    return vector2D( other * x, other * y);
  }

  vector2D operator*(const double other) const {
    return vector2D( other * x, other * y);
  }

};

void print2txt(Particle **p, int n) {
  printf("\nid; xPos; yPos; xVel; yVel; Mass;\n");
  for (int i = 0; i < n; i++) {
    printf("%06d; %.3f; %.3f; %.3f; %.3f; %.3f;\n", i, p[i]->xpos, p[i]->ypos,
           p[i]->xvel, p[i]->yvel, p[i]->mass);
  }
}

void print2svg(char *fname, Particle **p, int n) {
  FILE *fp;
  fp = fopen(fname, "w");
  fprintf(fp, "<?xml version=\"1.0\" standalone=\"no\"?>\n\
 	   <!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n\
 	   \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n\
 	   <svg width=\"500mm\" height=\"500mm\" viewBox=\"0 0 500 500\"\n\
 	   xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n");
  for (int i = 0; i < n; i++) {
    fprintf(fp, "<circle\n \
 		style = \"fill:#ff0000;stroke:none;stroke-width:0.5;\"\n \
 		cx = \"%.5f\"\n \
 	    	cy = \"%.5f\"\n \
 	     	r = \"1\" />\n",
            400 * (p[i]->xpos + L/2) + 50, 400 * (p[i]->ypos + L/2) + 50);
  }
  fprintf(fp, "</svg>\n");
  fclose(fp);
}

int main(int argc, char **argv) {
  int n, nt;
  char fname[64];
  double dt = 1e-6;

  double starttime = omp_get_wtime();

  // Read input arguments from console
  if (argc > 2) {
    n = atoi(argv[1]);
    nt = atoi(argv[2]);
    if (argc > 3) {
      dt = atof(argv[3]);
    }
  } else {
    printf("Usage: %s n nt (dt) \n \
          (mandatory) number of particles = n \n \
          (mandatory) number of time steps = nt\n\
          (optional) time step size = dt [Default 1e-6]\n\
          ",
           argv[0]);
    printf("========================================================\n");
    return -1;
  }

  srand(time(NULL));
  // Generate n particles with (random) coordinates,
  // velocities and masses
  double xpos, ypos, xvel, yvel, mass;
  double dist, force, xupdate, yupdate;

  Particle **p = new Particle *[n];
  for (int i = 0; i < n; i++) {
    xpos = ((double)rand() / RAND_MAX - 0.5)*L;
    ypos = ((double)rand() / RAND_MAX - 0.5)*L;

    // give it a little spin ;)
    xvel = -100*ypos; //(double)rand() / RAND_MAX ;
    yvel = 100*xpos; //(double)rand() / RAND_MAX ;
    mass = 1; //(double)rand() / RAND_MAX; mass set to 1 because the shown particles have the same size
    p[i] = new Particle(xpos, ypos, xvel, yvel, mass);
  }

  // Iterate over time steps

  int time_average = 0;

  for (int timestep = 0; timestep < nt; timestep++) {
    
    starttime = omp_get_wtime();
  #pragma omp parallel for shared(p) num_threads(12)
    for (int i = 0; i < n; i++){
      vector2D F;
      double temp = 0;
      double xdiff = 0;
      double ydiff = 0;
      double d = 0;
      for (int j = 0; j < n; j++)
      {
        xdiff = p[i]->xpos - p[j]->xpos;
        ydiff = p[i]->ypos - p[j]->ypos;
        d = xdiff*xdiff + ydiff*ydiff;
        if ((d > 1e-7)) {
            temp = (grav_const) / pow(d, 1.5) - repulsion / pow(d, 3);
        }
        else{
            temp = 0;
        }
        F = F + vector2D(-xdiff * temp, -ydiff * temp);
      }
      p[i]->xvel += dt * F.x;
      p[i]->yvel += dt * F.y;
    }

    time_average += (omp_get_wtime() - starttime)*1e6;

    // clearup space for the next Tree

    // Update coordinates of all particles based on
    // computed velocities
    for (int i = 0; i < n; i++) {
      // if the particle is far out, just ignore it
      double outside_factor = 1;
      double restitution_factor = .5;
      double reflection_scale = 0.001;
      p[i]->xpos += dt * p[i]->xvel;
      p[i]->ypos += dt * p[i]->yvel;
      if (p[i]->xpos > outside_factor*L/2 || p[i]->xpos < -outside_factor*L/2){
        //continue;
        p[i]->xvel *= -restitution_factor;
        if (p[i] -> xpos > 0) {
          p[i]->xpos = L/2 - (L*reflection_scale);
        }else {
          p[i]->xpos = -L/2 + (L*reflection_scale);
        }
      }
      if (p[i]->ypos > outside_factor*L/2 || p[i]->ypos < -outside_factor*L/2) {
        p[i]->yvel *= -restitution_factor;
        if (p[i] -> ypos > 0) {
          p[i]->ypos = L/2 - (L*reflection_scale);
        }else {
          p[i]->ypos = -L/2 + (L*reflection_scale);
        }
      }
    }

    // Save particles as a vector graphic in svg format
    sprintf(fname, "frames/frame_%05d.svg", timestep);
    print2svg(fname, p, n);

    // Some textual output for debugging
    // print2txt(p, n);
  }
  cout << "Average Time to Loop Over Tree: " << time_average/nt <<"ns" << endl;

  // Remove particles
  for (int i = 0; i < n; i++) {
    delete p[i];
  }
  delete[] p;
  return 0;
}
