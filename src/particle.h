#ifndef PARTICLE_H
#define PARTICLE_H

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "CGL/vector3D.h"

#include <vector>

using namespace CGL;
using namespace std;

struct Particle {
  Particle(Vector3D position)
      : position(position),
        last_position(position),
		velocity(Vector3D()),
		temp_velocity(Vector3D()),
		delta_p(Vector3D()),
        neighbors(new vector<Particle*>()),
		omega(Vector3D()){}

  // dynamic values
  Vector3D position;
  Vector3D last_position;
	Vector3D velocity;
	Vector3D temp_velocity;
	Vector3D delta_p;
  vector<Particle*> *neighbors;
	double lambda; // As calculated in line 10 of Algorithm 1
	Vector3D omega;
	Vector3D forces;
	double rho;
	Vector3D normal;
};


#endif /* PARTICLE_H */
