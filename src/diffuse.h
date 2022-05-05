#ifndef DIFFUSE_H
#define DIFFUSE_H

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "CGL/vector3D.h"
#include "particle.h"

#include <vector>

using namespace CGL;
using namespace std;

enum particle_type { SPRAY = 0, FOAM = 1, BUBBLE = 2 };

struct DiffuseParticle {
    //TODO: Mass might be wack
  DiffuseParticle(Vector3D position, Vector3D velocity, particle_type type, double ttl = 1)
      : position(position),
        last_position(position),
		velocity(Vector3D()),
		last_velocity(Vector3D()),
        type(type),
        ttl(ttl),
        neighbors(new vector<Particle*>()){}

  // dynamic values
  Vector3D position;
  Vector3D last_position;
  Vector3D velocity;
	Vector3D last_velocity;
  vector<Particle*> *neighbors;
	particle_type type;
    double ttl;
    double mass;
};

#endif /* DIFFUSE_H */
