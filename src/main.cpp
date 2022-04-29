#include "fluid.h"
#include "particle.h"
#include "marchingCube.h"

using namespace std;

int main(int argc, char **argv) {
	// Just testing build
	Fluid *f = new Fluid(4, 4, 4, 40, 40, 40);
	//Fluid* f = new Fluid(10, 10, 10, 10, 10, 10);
	f->buildGrid();
	f->build_spatial_map(1.);

	/*
	for (auto p = begin(f->particles); p != end(f->particles); p++) {
		std::cout << p->position << '\n';
	}
	*/

	// Call the constructor to creat marchingCube object
	Vector3D bDim = Vector3D(4., 4., 4.);
	Vector3D partDim = Vector3D(40., 40., 40.);
	float h = 1.;
	float search_radius = 1.0;
	float particle_mass = 1.;
	float density = 1.;
	float isovalue = 0.1;

	marchingCube* m = new marchingCube(bDim, partDim, f->particles, f->map, h, search_radius, particle_mass, density, isovalue);
	m->main_March("test.obj");

  return 0;
}
