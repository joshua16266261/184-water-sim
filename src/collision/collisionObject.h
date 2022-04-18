#ifndef COLLISIONOBJECT
#define COLLISIONOBJECT

#include "../particle.h"

using namespace CGL;
using namespace std;

class CollisionObject {
public:
  virtual void collide(Particle &p) = 0;

private:
  double friction;
};

#endif /* COLLISIONOBJECT */
