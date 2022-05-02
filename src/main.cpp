#include "fluid.h"
#include "particle.h"
#include "marchingCube.h"


#include <fstream>

using namespace std;

vector<Particle> readtxt(string filename) {

	ifstream infile(filename);
	float a, b, c;
	char c1, c2;
	char p1, p2;
	vector<Particle> par_pos = vector<Particle>();
	while (infile >> p1 >> a >> c1 >> b >> c2 >> c >> p2)
	{
		Vector3D pos = Vector3D(a, b, c);
		par_pos.emplace_back(Particle(pos));
	}
	return par_pos;
}

int main(int argc, char **argv) {
	// Just testing build
	Fluid *f = new Fluid(2, 2, 2, 40, 40, 40);
	//Fluid* f = new Fluid(10, 10, 10, 10, 10, 10);
	f->buildGrid();
	f->build_spatial_map(1.);



	/*
	for (auto p = begin(f->particles); p != end(f->particles); p++) {
		std::cout << p->position << '\n';
	}
	*/

	// Call the constructor to creat marchingCube object
	
	Vector3D bDim = Vector3D(2., 2., 2.);
	Vector3D partDim = Vector3D(40., 40., 40.);
	float h = 1.;
	float search_radius = .05;
	float particle_mass = 1.;
	float density = 1.;
	float step_size_multiplier = 1.;
	float isovalue = 1710.;
	//vector<Particle> par_pos = readtxt("fluid_particles.txt");

	marchingCube* m = new marchingCube(bDim, partDim, f->particles, f->map, h, search_radius, particle_mass, density, isovalue, step_size_multiplier);
	m->main_March("test.obj");
	
	
	//vector<Particle> par_pos = readtxt("pos.txt");
	//cout << par_pos.size() << "\n";
	return 0;
}


