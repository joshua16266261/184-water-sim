#include "../fluid.h"
#include "sphere.h"
#include "../particle.h"

using namespace CGL;

void Sphere::collide(Particle &p) {
  // TODO (Part 3): Handle collisions with spheres.
	Vector3D path = p.position - origin;
	if (path.norm() < radius) {
		path.normalize();
		path *= radius;
		path += origin;
		Vector3D corr = path - p.last_position;
		p.position = p.last_position + (1 - friction) * corr;
    }
}

