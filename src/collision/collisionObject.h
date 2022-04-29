#ifndef COLLISIONOBJECT
#define COLLISIONOBJECT

#include <nanogui/nanogui.h>

#include "../particle.h"
#include "../diffuse.h"

using namespace CGL;
using namespace std;
using namespace nanogui;

class CollisionObject {
public:
  virtual void collide(Particle &p, double cr, double delta_t) = 0;
  virtual void collide(DiffuseParticle &p) = 0;

private:
  double friction;
};

#endif /* COLLISIONOBJECT */
