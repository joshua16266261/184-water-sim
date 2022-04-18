#ifndef COLLISIONOBJECT_PLANE_H
#define COLLISIONOBJECT_PLANE_H

#include "../fluid.h"
#include "collisionObject.h"

using namespace CGL;
using namespace std;

struct Plane : public CollisionObject {
public:
  Plane(const Vector3D &point, const Vector3D &normal, double friction)
      : point(point), normal(normal.unit()), friction(friction) {}

  void collide(Particle &p);
  Vector3D point;
  Vector3D normal;

  double friction;
};

#endif /* COLLISIONOBJECT_PLANE_H */
