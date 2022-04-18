#ifndef PARTICLE_H
#define PARTICLE_H

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "CGL/vector3D.h"

#include <vector>

using namespace CGL;

struct Particle {
  Particle(Vector3D position)
      : position(position),
        last_position(position),
        neighbors = new vector<Particle*> {}

  Vector3D velocity(double delta_t) {
    return (position - last_position) / delta_t;
  }

  // dynamic values
  Vector3D position;
  Vector3D last_position;
  vector<Particle*> *neighbors;
};

#endif /* PARTICLE_H */
