#ifndef PARTICLE_H
#define PARTICLE_H

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "CGL/vector3D.h"

using namespace CGL;

struct Particle {
  Particle(Vector3D position, float mass)
      : position(position),
        last_position(position),
        mass(mass) {}

  Vector3D velocity(double delta_t) {
    return (position - last_position) / delta_t;
  }

  // static values
  float mass;

  // dynamic values
  Vector3D position;
  Vector3D last_position;
  Vector3D forces;

};

#endif /* PARTICLE_H */
