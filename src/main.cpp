#include "fluid.h"
#include "particle.h"
#include "collision/plane.h"
#include "marchingCube.h"
#include <vector>

#include <fstream>

using namespace std;

/*
void write_pos_to_file(Fluid* f, string filename) {
	string s = "";
	for (auto p = begin(f->particles); p != end(f->particles); p++) {
		s.append("(" + to_string(p->position[0]) + "," + to_string(p->position[1]) + "," + to_string(p->position[2]) + ")" + '\n');
	}
	ofstream file(filename);
	file << s;
	file.close();
}



int main(int argc, char** argv) {
	// Best params so far
	Fluid *f = new Fluid(4, 4, 4, 40, 40, 40);
	f->build_spatial_map(1.);

	FluidParameters *fp = new FluidParameters(EPS_F, 2, 60, 5);
	fp->h = 0.15;
	fp->delta_q = 0.02 * fp->h * Vector3D(1, 1, 1);
	fp->density = 1000;
	fp->k = 0.001;
	fp->c = 0.00025;
	fp->vorticity_eps = 0.005;
	fp->relaxation = 5000;
	fp->total_time = 1.;
	fp->fps = 60.;
	// Fluid* f = new Fluid(4, 4, 4, 40, 40, 40);

	//#pragma omp parallel for
	//for (auto p = begin(f->particles); p != end(f->particles); p++) {
	//	p->position += Vector3D(0, 0, 4);
	//}

	// FluidParameters* fp = new FluidParameters(EPS_F, 3, 60, 5);

	Vector3D g = Vector3D(0, 0, -9.81);
	vector<Vector3D> accel = vector<Vector3D>{ g };

	// They're doing this because they want to drop it but
	// we only set our box to be positive dimensions so rip lol
	Plane* floor = new Plane(Vector3D(0, 0, -2), Vector3D(0, 0, 1), 0);
	Plane* left_wall = new Plane(Vector3D(-2, 0, 0), Vector3D(1, 0, 0), 0);
	Plane* right_wall = new Plane(Vector3D(6, 0, 0), Vector3D(-1, 0, 0), 0);
	Plane* front_wall = new Plane(Vector3D(0, -2, 0), Vector3D(0, 1, 0), 0);
	Plane* back_wall = new Plane(Vector3D(0, 6, 0), Vector3D(0, -1, 0), 0);
	vector<CollisionObject*> collision = vector<CollisionObject*>{ floor, left_wall, right_wall, front_wall, back_wall };

	
	const string h = "--h";
	const string delta_q = "--q";
	const string k = "--k";
	const string c = "--c";
	const string vort = "--vort";
	const string relax = "--relax";
	const string time = "--time";
	const string fps = "--fps";

	for (int i = 1; i < argc - 1; i++) {
		std::cout << argv[i] << '\n';
		if (argv[i] == h) {
			fp->h = stod(argv[i + 1]);
		}
		else if (argv[i] == delta_q) {
			// --q param must come after --h param
			fp->delta_q = stod(argv[i + 1]) * fp->h * Vector3D(1, 0, 0);
		}
		else if (argv[i] == k) {
			fp->k = stod(argv[i + 1]);
		}
		else if (argv[i] == c) {
			fp->c = stod(argv[i + 1]);
		}
		else if (argv[i] == vort) {
			fp->vorticity_eps = stod(argv[i + 1]);
		}
		else if (argv[i] == relax) {
			fp->relaxation = stod(argv[i + 1]);
		}
		else if (argv[i] == time) {
			fp->total_time = stod(argv[i + 1]);
		}
		else if (argv[i] == fps) {
			fp->fps = stod(argv[i + 1]);
		}
	}
	

	//	fp->h = 0.15;
	//	fp->delta_q = 0.1 * fp->h * Vector3D(1, 0, 0);
	//  fp->density = 1000;
	//	fp->k = 0.001;
	//	fp->c = 0.0005;
	//	fp->c = 0.00075;
	//	fp->vorticity_eps = 0.0005;
	//	fp->relaxation = 4000;

	for (int frame = 0; frame < fp->total_time * fp->fps; frame++) {
		
		
		
		std::cout << frame << '\n';
		f->simulate(fp, accel, &collision);
		cout << "Simulated frame " + to_string(frame) << endl;

		write_pos_to_file(f, "floor " + to_string(frame) + ".txt");

		// Fluid* f = new Fluid(4, 4, 4, 40, 40, 40);


		Vector3D bDim = Vector3D(4., 4., 4.);
		Vector3D partDim = Vector3D(40., 40., 40.);
		float search_radius = .01;
		float particle_mass = 1.;
		float step_size_multiplier = 0.5;
		float isovalue = 0.001;
		
		marchingCube* m = new marchingCube(bDim, partDim, f->particles, f->map, fp->h, search_radius,
			particle_mass, fp->density, isovalue, step_size_multiplier);

		m->main_March("Frame-" + to_string(frame) + ".obj");
		cout << "Generated frame #" + to_string(frame) << endl;
	}

	return 0;
}
*/



/////////////////////////////
// CODE FOR TESTING IT OUT //
/////////////////////////////
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

int main(int argc, char** argv) {
	// Just testing build
	Fluid* f = new Fluid(2, 2, 2, 40, 40, 40);
	f->buildGrid();
	f->build_spatial_map(1.);



	
	//for (auto p = begin(f->particles); p != end(f->particles); p++) {
	//	std::cout << p->position << '\n';
	//}
	

	// Call the constructor to creat marchingCube object

	Vector3D bDim = Vector3D(2., 2., 2.);
	Vector3D partDim = Vector3D(40., 40., 40.);
	float h = 1.;
	float search_radius = .03;
	float particle_mass = 1.;
	float density = 1000.;
	float step_size_multiplier = 0.5;
	float isovalue = 0.1;
	//vector<Particle> par_pos = readtxt("fluid_particles.txt");

	marchingCube* m = new marchingCube(bDim, partDim, f->particles, f->map, h, search_radius, particle_mass, density, isovalue, step_size_multiplier);
	m->main_March("test.obj");


	//vector<Particle> par_pos = readtxt("pos.txt");
	//cout << par_pos.size() << "\n";
	return 0;
}
