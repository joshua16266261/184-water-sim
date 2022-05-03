#include "fluid.h"
#include "parent_fluid.h"
#include "particle.h"
#include "collision/plane.h"
#include <vector>

#include <fstream>

using namespace std;

void write_fluid_pos_to_file(Fluid *f, string filename) {
	string s = "";
	for (auto p = begin(f->particles); p != end(f->particles); p++) {
		s.append("(" + to_string(p->position[0]) + "," + to_string(p->position[1]) + "," + to_string(p->position[2]) + ")" + '\n');
	}
	ofstream file(filename);
	file << s;
	file.close();
}

void write_diffuse_pos_to_file_helper(ParentFluid *pf, string filename, particle_type type) {
	string s = "";
	for (auto p = begin(*pf->diffuse_particles); p != end(*pf->diffuse_particles); p++) {
		if ((*p)->type == type) {
			s.append("(" + to_string((*p)->position[0]) + "," + to_string((*p)->position[1]) + "," + to_string((*p)->position[2]) + ")" + '\n');
		}
	}
	ofstream file(filename);
	file << s;
	file.close();
}

void write_diffuse_pos_to_file(ParentFluid *pf, string filename) {
	write_diffuse_pos_to_file_helper(pf, "foam_" + filename, FOAM);
	write_diffuse_pos_to_file_helper(pf, "bub_" + filename, BUBBLE);
	write_diffuse_pos_to_file_helper(pf, "spray_" + filename, SPRAY);
}

int main(int argc, char **argv) {
	Fluid *f = new Fluid(4, 4, 4, 40, 40, 40);

	#pragma omp parallel for
	for (auto p = begin(f->particles); p != end(f->particles); p++) {
		p->position += Vector3D(0, 0, 4);
	}

//	FluidParameters *fp = new FluidParameters(EPS_F, 2.5, 60, 5);
	FluidParameters *fp = new FluidParameters(EPS_F, 7, 60, 5);

	Vector3D g = Vector3D(0, 0, -9.81);
	vector<Vector3D> accel = vector<Vector3D>{g};

	Plane *floor = new Plane(Vector3D(0, 0, -2), Vector3D(0, 0, 1), 0);
	Plane *left_wall = new Plane(Vector3D(-2, 0, 0), Vector3D(1, 0, 0), 0);
	Plane *right_wall = new Plane(Vector3D(6, 0, 0), Vector3D(-1, 0, 0), 0);
	Plane *front_wall = new Plane(Vector3D(0, -2, 0), Vector3D(0, 1, 0), 0);
	Plane *back_wall = new Plane(Vector3D(0, 6, 0), Vector3D(0, -1, 0), 0);
	vector<CollisionObject *> collision = vector<CollisionObject *>{floor, left_wall, right_wall, front_wall, back_wall};


//	const string h = "--h";
//	const string delta_q = "--q";
//	const string k = "--k";
//	const string c = "--c";
//	const string vort = "--vort";
//	const string relax = "--relax";
	const string time = "--time";
	const string fps = "--fps";

	for (int i = 1; i < argc - 1; i++) {
//		if (argv[i] == h) {
//			fp->h = stod(argv[i + 1]);
//		} else if (argv[i] == delta_q) {
//			// --q param must come after --h param
//			fp->delta_q = stod(argv[i + 1]) * fp->h * Vector3D(1, 0, 0);
//		} else if (argv[i] == k) {
//			fp->k = stod(argv[i + 1]);
//		} else if (argv[i] == c) {
//			fp->c = stod(argv[i + 1]);
//		} else if (argv[i] == vort) {
//			fp->vorticity_eps = stod(argv[i + 1]);
//		} else if (argv[i] == relax) {
//			fp->relaxation = stod(argv[i + 1]);
//		} else if (argv[i] == time) {
//			fp->total_time = stod(argv[i + 1]);
//		} else if (argv[i] == fps) {
//			fp->fps = stod(argv[i + 1]);
//		}
		if (argv[i] == time) {
			fp->total_time = stod(argv[i + 1]);
		} else if (argv[i] == fps) {
			fp->fps = stod(argv[i + 1]);
		}
	}

	fp->h = 0.15;
	fp->delta_q = 0.1 * fp->h * Vector3D(1, 0, 0);
	fp->density = 1000;
	fp->k = 0.001;
	fp->c = 0.0008;
	fp->vorticity_eps = 0.0002;
	fp->relaxation = 1600;

	DiffuseParameters *dp = new DiffuseParameters(fp->h, 0.5, 0.5, 1000, 1000, 13);
	ParentFluid *pf = new ParentFluid(4, 4, 4, 40, 60, 2.5, dp);
	pf->fp = fp;

	for (int frame = 0; frame < fp->total_time * fp->fps; frame++) {
		std::cout << frame << '\n';
		f->simulate(fp, accel, &collision);
		pf->fluid = f;
		pf->simulate_step(accel, &collision);
		write_fluid_pos_to_file(pf->fluid, "fluid_" + to_string(frame) + ".txt");
		write_diffuse_pos_to_file(pf, "diffuse_" + to_string(frame) + ".txt");
	}

//	#pragma omp parallel for
//	for (auto p = begin(f->particles); p != end(f->particles); p++) {
//		p->position /= 4;
//		p->position += Vector3D(0.5, 0.5, 0.5);
//	}
//
//	write_pos_to_file(f, "fluid_particles.txt");

  return 0;
}
