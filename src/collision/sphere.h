#ifndef COLLISIONOBJECT_SPHERE_H
#define COLLISIONOBJECT_SPHERE_H

#include "../fluid.h"
#include "collisionObject.h"
#include "../particle.h"

using namespace CGL;
using namespace std;

struct Sphere : public CollisionObject {
public:
  Sphere(const Vector3D &origin, double radius, double friction, int num_lat = 40, int num_lon = 40)
      : origin(origin), radius(radius), radius2(radius * radius),
        friction(friction) {}

    void collide(Particle &p);
  

private:
  Vector3D origin;
  double radius;
  double radius2;

  double friction;
  
};

#endif /* COLLISIONOBJECT_SPHERE_H */


