#include "iostream"
#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../clothSimulator.h"
#include "plane.h"

using namespace std;
using namespace CGL;

#define SURFACE_OFFSET 0.0001

void Plane::collide(Particle &p, float cr, float delta_t) {
	//TODO: Implement plane collision
//	if (dot(pm.position - point, normal) * dot(pm.last_position - point, normal) <= 0) {
//		Vector3D dir = dot(point - pm.position, normal) / dot(normal, normal) * normal;
//		Vector3D tang = pm.position + dir;
//		int sign = 1;
//		if (dot(pm.last_position - point, normal) < 0) {
//			sign = -1;
//		}
//		Vector3D corr = tang - pm.last_position + sign * SURFACE_OFFSET * normal;
//		pm.position = pm.last_position + (1 - friction) * corr;
//	}
	
	Vector3D pos = p.position + p.delta_p;
	if (dot(pos - point, normal) < 0) {
		Vector3D dir = dot(point - pos, normal) / dot(normal, normal) * normal;
		Vector3D tang = p.position + dir;
		float penetration = (pos - tang).norm();
		p.velocity -= (1 + cr * penetration / (delta_t * p.velocity.norm())) * dot(p.velocity, normal) * normal;
		p.delta_p = tang - p.position + SURFACE_OFFSET * normal;
	}
	
}


