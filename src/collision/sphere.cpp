#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../misc/sphere_drawing.h"
#include "../particle.h"
#include "sphere.h"

using namespace nanogui;
using namespace CGL;

void Sphere::collide(Particle &p, double cr, double delta_t) {
	//TODO: Implement sphere collision
//	Vector3D path = pm.position - origin;
//	if (path.norm() < radius) {
//		path.normalize();
//		path *= radius;
//		path += origin;
//		Vector3D corr = path - pm.last_position;
//		pm.position = pm.last_position + (1 - friction) * corr;
//	}
}


