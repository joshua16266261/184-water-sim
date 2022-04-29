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
  FluidParameters(double relaxation, float total_time, int fps, int solverIters, float h = 1.0, float k = 0.1, float n = 4, float c = 0.01, float delta_q = 0.2, double density = 1000, float cr = 0.45, float vorticity_eps = 0.01)
      : density(density), relaxation(relaxation), delta_q(delta_q), total_time(total_time), fps(fps), solverIters(solverIters), h(h), k(k), n(n), c(c), cr(cr), vorticity_eps(vorticity_eps) {}
  ~FluidParameters() {}

  // Simulation parameters
    float total_time;
    int fps;
    int solverIters;
    float h;

  // Fluid parameters
  double density;
  double relaxation;
  float delta_q;
  float k;
  float n;
  float c;
	float cr;
	float vorticity_eps;
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
	
	void set_neighbors(Particle *p, float h);
	float get_avg_spacing();
	
	float poly6_kernel(Vector3D r, float h) {
		if (r.norm() > h) {
			return 0;
		}
		return 315.0 / (64 * PI * pow(h, 9)) * pow(pow(h, 2) - r.norm2(), 3);
	}
	
	Vector3D grad_spiky_kernel(Vector3D r, float h) {
		if (r.norm() > h) {
			return Vector3D();
		}
		return -r * 45.0 / (PI * pow(h, 6) * r.norm()) * pow(h - r.norm(), 2);
	}
	
	float M4_kernel(Vector3D r, float h) {
		if (r.norm() > 2 * h) {
			return 0;
		}
		float q = r.norm() / h;
		if (q < 1) {
			return 10.0 / (7 * PI * pow(h, 2)) * (1 - 1.5 * pow(q, 2) + 0.75 * pow(q, 3));
		} else if (q < 2) {
			return 10.0 / (28 * PI * pow(h, 2)) * pow(2 - q, 3);
		}
	}

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
