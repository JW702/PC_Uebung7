// Use "convert -delay 1 -loop 0 *.svg result.gif" to generate an animated gif
// from svg-files
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

static int global_num_particles = 0;

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

//QuadTree for Barnes Hutt Algorithm

struct QuadTree {
  int level;

  // Center of the Tree/node
  double xcenter;
  double ycenter;

  // side length 
  double l;

  Particle **particles = nullptr;
  int num_particles = 0;

  QuadTree* ne = nullptr;
  QuadTree* nw = nullptr;
  QuadTree* sw = nullptr;
  QuadTree* se = nullptr;

  // contructor
  QuadTree(double _xcenter, double _ycenter, double _l, Particle **p, int n, int depth = 0){
    xcenter = _xcenter;
    ycenter = _ycenter;
    l = _l;
    num_particles = n;
    particles = p;
    // if 1 or less particles are in the node, do not subdivide further
    if(n <= 1 || depth > 2*log2(global_num_particles)){
      return;
    }


    int n_ne = 0;
    int n_nw = 0;
    int n_sw = 0;
    int n_se = 0;

    // sets the gap which is already sorted
    int k = 0;

    for(int i = 0; i < n; i++){
      // check if p[i] is in the north eastern node of the tree
      if(p[i]->xpos >= xcenter && p[i]->ypos <= ycenter){

        // put the particle in the first unsorted spot
        Particle *temp = p[n_ne + k];
        p[n_ne + k] = p[i];
        p[i] = temp;
        n_ne++;
      }
    }

    k += n_ne;
    for(int i = k; i < n; i++){
      // check if p[i] is in the north western node of the tree
      if(p[i]->xpos <= xcenter && p[i]->ypos <= ycenter){

        // put the particle in the first unsorted spot
        Particle *temp = p[k + n_nw];
        p[k + n_nw] = p[i];
        p[i] = temp;
        n_nw++;
      }
    }

    k += n_nw;
    for(int i = k; i < n; i++){
      // check if p[i] is in the south western node of the tree
      if(p[i]->xpos <= xcenter && p[i]->ypos >= ycenter){

        // put the particle in the first unsorted spot
        Particle *temp = p[k + n_sw];
        p[k + n_sw] = p[i];
        p[i] = temp;
        n_sw++;
      }
    }

    // the remaining points must be in SE
    n_se = n - k - n_sw;
    ne = new QuadTree(xcenter + l/4, ycenter - l/4, l/2, p, n_ne ,depth + 1);
    nw = new QuadTree(xcenter - l/4, ycenter - l/4, l/2, p + n_ne, n_nw, depth + 1);
    sw = new QuadTree(xcenter - l/4, ycenter + l/4, l/2, p + n_ne + n_nw, n_sw, depth + 1);
    se = new QuadTree(xcenter + l/4, ycenter + l/4, l/2, p + n_ne + n_nw + n_sw, n_se, depth + 1);
  }

  bool isLeaf() const {
    return (nw == nullptr && ne == nullptr && sw == nullptr && se == nullptr);
  }

  void clear() {
        // recursively clear the children
        if (nw) { nw->clear(); delete nw; nw = nullptr; }
        if (ne) { ne->clear(); delete ne; ne = nullptr; }
        if (sw) { sw->clear(); delete sw; sw = nullptr; }
        if (se) { se->clear(); delete se; se = nullptr; }

        particles = nullptr;
    }

  ~QuadTree(){
    clear();
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

// calculating the Force acting on a particle 
vector2D ForceToQuadTree(double xPos, double yPos, QuadTree* QT){
  // empty leaf -> no force from it
  if(QT->num_particles <= 0){
    return vector2D(0,0);
  }

  // determine if it is a cluster (avoid the sqrt)
  double xdiff = xPos - QT->xcenter;
  double ydiff = yPos - QT->ycenter;
  double d = xdiff*xdiff + ydiff*ydiff;
  double l = QT->l;

  double force = 0;
  // smallest possible node, contains only one particle -> direct calculation is as fast as clustering
  if(QT->isLeaf()){
    xdiff = xPos - QT->particles[0]->xpos;
    ydiff = yPos - QT->particles[0]->ypos;
    d = xdiff*xdiff + ydiff*ydiff;
    if ((d > 1e-8)) {
      force = (grav_const * 1) / pow(d, 1.5) - repulsion / pow(d, 3);
      return  vector2D(-xdiff * force, -ydiff * force);
    }
  }
  else if (xPos >= 2*L || xPos <= -2*L || yPos >=  2*L || yPos <= -2*L){
    // if the particle is far out, just ignore it
    return vector2D(0,0);
  }
  else if (l*l <= d * theta*theta){
    // do the cluster stuff (calculate force)
    if ((d > 1e-8)) {
      force = (grav_const * QT->num_particles) / pow(d, 1.5) - repulsion / pow(d, 3);
      return vector2D(-xdiff * force, -ydiff * force);
    }
  }
  else {
    vector2D ForceNE, ForceNW, ForceSW, ForceSE, Force;
    ForceNE = ForceToQuadTree(xPos, yPos, QT->ne);
    ForceNW = ForceToQuadTree(xPos, yPos, QT->nw);
    ForceSW = ForceToQuadTree(xPos, yPos, QT->sw);
    ForceSE = ForceToQuadTree(xPos, yPos, QT->se);
    return vector2D(ForceNE + ForceNW + ForceSW + ForceSE);    
  }
  return vector2D(0,0);
}

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
    global_num_particles = n;
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
    // regenerating the QuadTree foreach timestep
    QuadTree *QTree = new QuadTree(0, 0, 2*L, p, n);
  //  if (timestep % 100 == 0){
    //  cout << " Generating the Tree for step " << timestep << " took: " << (omp_get_wtime() - starttime)*1e6 << "ns" << endl;
    //}
     
    
    //

    // This part sucks
    starttime = omp_get_wtime();
#pragma omp parallel for shared(p) num_threads(1)
    for (int i = 0; i < n; i++){
      vector2D F;
      // RK4? beacause Euler only converges in O(dt^2) and not O(dt^4) now just midpoint which is O(dt^3)
      {
        F = ForceToQuadTree(p[i]->xpos, p[i]->ypos, QTree);
      }
      p[i]->xvel += dt * F.x;
      p[i]->yvel += dt * F.y;
    }
   // if (timestep % 100 == 0){

     // cout << " Looping over the Tree for step " << timestep << " took: " << (omp_get_wtime() - starttime)*1e6 << "ns" <<endl;
    //}
    time_average += (omp_get_wtime() - starttime)*1e6;

    // clearup space for the next Tree
    QTree->clear();

    // Update coordinates of all particles based on
    // computed velocities
    for (int i = 0; i < n; i++) {
      // if the particle is far out, just ignore it
      double outsite_factor = 1;
      double bounce_scale = .1;
      double reflection_scale = 0.001;
      p[i]->xpos += dt * p[i]->xvel;
      p[i]->ypos += dt * p[i]->yvel;
      if (p[i]->xpos > outsite_factor*L/2 || p[i]->xpos < -outsite_factor*L/2){
        //continue;
        p[i]->xvel *= -bounce_scale;
        if (p[i] -> xpos > 0) {
          p[i]->xpos = L/2 - (L*reflection_scale);
        }else {
          p[i]->xpos = -L/2 + (L*reflection_scale);
        }
      }
      if (p[i]->ypos > outsite_factor*L/2 || p[i]->ypos < -outsite_factor*L/2) {
        p[i]->yvel *= -bounce_scale;
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
