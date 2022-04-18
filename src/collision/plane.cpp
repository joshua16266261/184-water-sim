#include "iostream"

#include "../fluid.h"
#include "plane.h"
#include "../particle.h"

using namespace std;
using namespace CGL;

#define SURFACE_OFFSET 0.0001

void Plane::collide(Particle &p) {
  // TODO (Part 3): Handle collisions with planes.
	if (dot(p.position - point, normal) * dot(p.last_position - point, normal) <= 0) {
		Vector3D dir = dot(point - p.position, normal) / dot(normal, normal) * normal;
		Vector3D tang = p.position + dir;
		int sign = 1;
		if (dot(p.last_position - point, normal) < 0) {
			sign = -1;
		}
		Vector3D corr = tang - p.last_position + sign * SURFACE_OFFSET * normal;
		p.position = p.last_position + (1 - friction) * corr;
	}
}

