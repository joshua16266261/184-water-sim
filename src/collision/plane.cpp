#include "iostream"
#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../clothSimulator.h"
#include "plane.h"

using namespace std;
using namespace CGL;

#define SURFACE_OFFSET 0.0001

void Plane::collide(Particle &p, double cr, double delta_t) {
	//TODO: Revert to Proj 4 collisions
//	Vector3D pos = p.position + p.delta_p;
//	if (dot(pos - point, normal) * dot(p.last_position - point, normal) <= 0) {
//		Vector3D dir = dot(point - pos, normal) / dot(normal, normal) * normal;
//		Vector3D tang = pos + dir;
//		int sign = 1;
//		if (dot(p.last_position - point, normal) < 0) {
//			sign = -1;
//		}
//		Vector3D corr = tang - p.position + sign * SURFACE_OFFSET * normal;
////		p.position = p.last_position + (1 - friction) * corr;
//		p.delta_p = corr;
//	}
	
//	Vector3D pos = p.position + p.delta_p;
//	if (dot(pos - point, normal) < 0) {
//		Vector3D dir = dot(point - pos, normal) / dot(normal, normal) * normal;
//		Vector3D tang = pos + dir;
////		std::cout << tang << '\n';
//		Vector3D corr = tang - p.position + SURFACE_OFFSET * normal;
////		p.position = p.last_position + (1 - friction) * corr;
//		p.delta_p = corr;
//	}
	
	Vector3D pos = p.position + p.delta_p;
	if (dot(pos - point, normal) < 0) {
		Vector3D dir = dot(point - pos, normal) / dot(normal, normal) * normal;
		Vector3D tang = p.position + dir;
		double penetration = (pos - tang).norm();
		p.temp_velocity -= (1 + cr * penetration / (delta_t * p.temp_velocity.norm())) * dot(p.temp_velocity, normal) * normal;
		p.delta_p = tang - p.position + SURFACE_OFFSET * normal;
	}
	
}


