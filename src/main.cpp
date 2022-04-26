#include "fluid.h"
#include "particle.h"
#include "marchingCube.h"

using namespace std;

int main(int argc, char **argv) {
	// Just testing build
	Fluid *f = new Fluid(50, 50, 50, 20, 20, 20);
	f->buildGrid();
	/*
	for (auto p = begin(f->particles); p != end(f->particles); p++) {
		std::cout << p->position << '\n';
	}
	*/

	// Call the constructor to creat marchingCube object
	Vector3D bDim = Vector3D(50, 50., 50.);
	Vector3D unitDim = Vector3D(1., 1., 1.);
	float h = 1.;
	float search_radius = 1.;
	float particle_mass = 10.;
	float density = 1.;
	float isovalue = 0.1;

	marchingCube* m = new marchingCube(bDim, unitDim, f->particles, f->map, h, search_radius, particle_mass, density, isovalue);
	m->main_March("test.obj");

  return 0;
}
