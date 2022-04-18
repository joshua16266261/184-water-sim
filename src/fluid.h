#ifndef FLUID_H
#define FLUID_H

#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "CGL/vector3D.h"
#include "collision/collisionObject.h"
#include "particle.h"

using namespace CGL;
using namespace std;

struct FluidParameters {
  FluidParameters() {}
  /*
    relaxation = epsilon (see page 2 of Macklin and Muller)
    delta_q is a point some fixed distance inside the smoothing kernel radius (see page 3 of Macklin and Muller)
    k is a small positive constant used in calculating an artificial pressure (see page 3 of Macklin and Muller)
    n is some constant used in calculating an artificial pressure (see page 3 of Macklin and Muller)
    c is some constant used in applying XSPH viscosity (see page 3 of Macklin and Muller)
    total_time is the number of seconds that this simulation will run for
    fps is the number of frames per second
    h is the max distance that 2 particles can be considered neighbors
    */
  FluidParameters( double density, double relaxation, Vector3D delta_q, float total_time, int fps, int solverIters, float hfloat, float h = 1.0, float k = 0.1, float n = 4, float c = 0.01)
      : density(density), relaxation(relaxation), delta_q(delta_q), total_time(total_time), fps(fps), solverIters(solverIters), h(h), k(k), n(n), c(c) {}
  ~FluidParameters() {}

  // Simulation parameters
    float total_time;
    int fps;
    int solverIters;
    float h;

  // Fluid parameters
  double density;
  double relaxation;
  Vector3D delta_q;
  float k;
  float n;
  float c;
};

struct Fluid {
  Fluid() {}
  Fluid(float length, float width, float height, int num_length_particles, int num_width_particles, int num_height_particles);
  ~Fluid();

  void buildGrid();

  void simulate(FluidParameters *fp,
                vector<Vector3D> external_accelerations,
                vector<CollisionObject *> *collision_objects);

  void build_spatial_map(float h);
  string hash_position(Vector3D pos, float h);

  // Fluid properties
  double length;
  double width;
  double height;
  int num_length_particles;
  int num_width_particles;
  int num_height_particles;

  // Fluid components
  vector<Particle> particles;

  // Spatial hashing
  unordered_map<string, vector<Particle *> *> map;
};

#endif /* FLUID_H */
