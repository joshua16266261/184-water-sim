#ifndef COLLISIONOBJECT
#define COLLISIONOBJECT

using namespace CGL;
using namespace std;
using namespace nanogui;

class CollisionObject {
public:
  virtual void collide(Particle &p) = 0;

private:
  double friction;
};

#endif /* COLLISIONOBJECT */
