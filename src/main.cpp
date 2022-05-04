#include "fluid.h"
#include "parent_fluid.h"
#include "particle.h"
#include "collision/plane.h"
#include <vector>

#include <fstream>

using namespace std;

void write_fluid_pos_to_file(Fluid *f, string filename) {
	// Write all fluid particle positions to a text file
	string s = "";
	for (auto p = begin(f->particles); p != end(f->particles); p++) {
		s.append("(" + to_string(p->position[0]) + "," + to_string(p->position[1]) + "," + to_string(p->position[2]) + ")" + '\n');
		// Divide by 4 and offset by 0.5 to constrain to 2x2x2 box with nonnegative coordinates
		// s.append("(" + to_string(p->position[0] / 4 + 0.5) + "," + to_string(p->position[1] / 4 + 0.5) + "," + to_string(p->position[2] / 4 + 0.5) + ")" + '\n');
	}
	ofstream file(filename);
	file << s;
	file.close();
}

void write_diffuse_pos_to_file_helper(ParentFluid *pf, string filename, particle_type type) {
	// Write all diffuse particle positions of particle_type type to a text file
	string s = "";
	for (auto p = begin(*pf->diffuse_particles); p != end(*pf->diffuse_particles); p++) {
		if ((*p)->type == type) {
			s.append("(" + to_string((*p)->position[0]) + "," + to_string((*p)->position[1]) + "," + to_string((*p)->position[2]) + ")" + '\n');
			// Divide by 4 and offset by 0.5 to constrain to 2x2x2 box with nonnegative coordinates
			// s.append("(" + to_string((*p)->position[0] / 4 + 0.5) + "," + to_string((*p)->position[1] / 4 + 0.5) + "," + to_string((*p)->position[2] / 4 + 0.5) + ")" + '\n');
		}
	}
	ofstream file(filename);
	file << s;
	file.close();
}

void write_diffuse_pos_to_file(ParentFluid *pf, string filename) {
	// Write all diffuse particle positions to text files based on their type
	write_diffuse_pos_to_file_helper(pf, "foam_" + filename, FOAM);
	write_diffuse_pos_to_file_helper(pf, "bub_" + filename, BUBBLE);
	write_diffuse_pos_to_file_helper(pf, "spray_" + filename, SPRAY);
}

int main(int argc, char **argv) {
	// DO NOT TOUCH: Initialize fluid particles
//	Fluid *f = new Fluid(4, 4, 4, 40, 40, 40);
	Fluid *f = new Fluid(4, 4, 1, 40, 40, 10);

	// Falling cube: shift starting position of fluid up
//	#pragma omp parallel for
//	for (auto p = begin(f->particles); p != end(f->particles); p++) {
//		p->position += Vector3D(0, 0, 4);
//	}

	// Default: 7 seconds at 60 fps
	FluidParameters *fp = new FluidParameters(EPS_F, 7, 60, 5);

	// Gravity
	Vector3D g = Vector3D(0, 0, -9.81);
	vector<Vector3D> accel = vector<Vector3D>{g};

	// Walls and floor
//	Plane *floor = new Plane(Vector3D(0, 0, -2), Vector3D(0, 0, 1), 0);
//	Plane *left_wall = new Plane(Vector3D(-2, 0, 0), Vector3D(1, 0, 0), 0);
//	Plane *right_wall = new Plane(Vector3D(6, 0, 0), Vector3D(-1, 0, 0), 0);
//	Plane *front_wall = new Plane(Vector3D(0, -2, 0), Vector3D(0, 1, 0), 0);
//	Plane *back_wall = new Plane(Vector3D(0, 6, 0), Vector3D(0, -1, 0), 0);
	Plane *floor = new Plane(Vector3D(0, 0, -0.1), Vector3D(0, 0, 1), 0);
	Plane *left_wall = new Plane(Vector3D(-0.1, 0, 0), Vector3D(1, 0, 0), 0);
	Plane *right_wall = new Plane(Vector3D(4.1, 0, 0), Vector3D(-1, 0, 0), 0);
	Plane *front_wall = new Plane(Vector3D(0, -0.1, 0), Vector3D(0, 1, 0), 0);
	Plane *back_wall = new Plane(Vector3D(0, 4.1, 0), Vector3D(0, -1, 0), 0);
	vector<CollisionObject *> collision = vector<CollisionObject *>{floor, left_wall, right_wall, front_wall, back_wall};
	
	// Command line args for video time and fps
	const string time = "--time";
	const string fps = "--fps";

	for (int i = 1; i < argc - 1; i++) {
		if (argv[i] == time) {
			fp->total_time = stod(argv[i + 1]);
		} else if (argv[i] == fps) {
			fp->fps = stod(argv[i + 1]);
		}
	}

	// DO NOT TOUCH: Fluid simulation parameters
	fp->h = 0.15;
	fp->delta_q = 0.1 * fp->h * Vector3D(1, 0, 0);
	fp->density = 1000;
	fp->k = 0.001;
	fp->c = 0.0008;
	fp->vorticity_eps = 0.0002;
	fp->relaxation = 1600;

	// CURRENTLY TUNING: Diffuse fluid parameters
	DiffuseParameters *dp = new DiffuseParameters(fp->h, 0.5, 0.5, 200, 60, 13);
	dp->t_k_min = 5;
	dp->t_k_max = 60;
	ParentFluid *pf = new ParentFluid(4, 4, 4, 40, 60, 2.5, dp);
	pf->fp = fp;
	
	double freq = 3.0 / fp->fps;

	// Simulate all frames
	for (int frame = 0; frame < fp->total_time * fp->fps; frame++) {
		std::cout << frame << '\n';
		
		// Simulate 1 step
		f->simulate(fp, accel, &collision);
		pf->fluid = f;
		pf->simulate_step(accel, &collision);
		
		// Write positions to text files
		write_fluid_pos_to_file(pf->fluid, "fluid_" + to_string(frame) + ".txt");
		write_diffuse_pos_to_file(pf, "diffuse_" + to_string(frame) + ".txt");
		
		if (frame > 2 * fp->fps) {
			left_wall->point.x = sin(freq * (frame - 2 * fp->fps));
		}
	}

	return 0;
}
