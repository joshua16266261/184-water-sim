#ifndef COLLISIONOBJECT
#define COLLISIONOBJECT

#include <nanogui/nanogui.h>

#include "../particle.h"

using namespace CGL;
using namespace std;
using namespace nanogui;

class CollisionObject {
public:
  virtual void collide(Particle &p, float cr, float delta_t) = 0;

private:
  double friction;
};

#endif /* COLLISIONOBJECT */
