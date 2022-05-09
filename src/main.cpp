#include "fluid.h"
#include "parent_fluid.h"
#include "particle.h"
#include "collision/plane.h"
#include "marchingCube.h"
#include <vector>

#include <fstream>

using namespace std;

void write_pos_to_file(Fluid* f, string filename) {
	// FOR TESTING ONLY
	string s = "";
	for (auto p = begin(f->particles); p != end(f->particles); p++) {
		s.append("(" + to_string(p->position[0]) + "," + to_string(p->position[1]) + "," + to_string(p->position[2]) + ")" + '\n');
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

void write_ply(vector<Particle> particles, string filename) {
	ofstream file(filename);
	file << "ply" << endl;
	file << "format ascii 1.0" << endl;
	file << "element vertex " << to_string(particles.size()) << endl;
	file << "property float x" << endl;
	file << "property float y" << endl;
	file << "property float z" << endl;
	file << "end_header" << endl;
	
	for (auto p = begin(particles); p != end(particles); p++) {
		Vector3D pos = p->position;
		file << pos.x << " " << pos.y << " " << pos.z << endl;
	}
	
	file.close();
}


int main(int argc, char** argv) {
	// Default: 4x4x4 cube with 40x40x40 particles
	Fluid *f = new Fluid(4, 4, 4, 40, 40, 40);

	// Falling cube: shift starting cube to have corners (2, 2, 2) and (6, 6, 6)
	#pragma omp parallel for
	for (auto p = begin(f->particles); p != end(f->particles); p++) {
		p->position += Vector3D(2, 2, 2);
	}
	
	// Gravity
	Vector3D g = Vector3D(0, 0, -9.81);
	vector<Vector3D> accel = vector<Vector3D>{g};

	// Walls and floor
	Plane* floor = new Plane(Vector3D(0, 0, 0), Vector3D(0, 0, 1), 0);
	Plane* left_wall = new Plane(Vector3D(0, 0, 0), Vector3D(1, 0, 0), 0);
	Plane* right_wall = new Plane(Vector3D(8, 0, 0), Vector3D(-1, 0, 0), 0);
	Plane* front_wall = new Plane(Vector3D(0, 0, 0), Vector3D(0, 1, 0), 0);
	Plane* back_wall = new Plane(Vector3D(0, 8, 0), Vector3D(0, -1, 0), 0);
	vector<CollisionObject*> collision = vector<CollisionObject*>{floor, left_wall, right_wall, front_wall, back_wall};

	// Command line args for video time and output fps
	const string time = "--time";
	const string fps = "--fps";
	const string no_march = "--no-march";
	// Alternatively, set total_time and output_fps
	double total_time = 10;
	int output_fps = 30;
	bool march = true;
	for (int i = 1; i < argc - 1; i++) {
		if (argv[i] == time) {
			total_time = stod(argv[i + 1]);
		} else if (argv[i] == fps) {
			output_fps = stoi(argv[i + 1]);
		} else if (argv[i] == no_march) {
			march = false;
		}
	}
	
	cout << "Total time: " + to_string(total_time) << endl;
	cout << "Output fps: " + to_string(output_fps) << endl;
	cout << "Marching: " + to_string(march) << endl;
	
	int downsample_rate = 60 / output_fps;
	
	/* DO NOT TOUCH ANYTHING BELOW THIS POINT */
	FluidParameters* fp = new FluidParameters(EPS_F, total_time, 60, 5);
	
	// Fluid simulation parameters
	fp->h = 0.15;
	fp->delta_q = 0.1 * fp->h * Vector3D(1, 0, 0);
	fp->density = 1000;
	fp->k = 0.001;
	fp->c = 0.0008;
	fp->vorticity_eps = 0.0002;
	fp->relaxation = 1600;
	
	// Marching cube parameters
	Vector3D bDim = Vector3D(2.5, 2.5, 2.5);
	Vector3D partDim = Vector3D(40., 40., 40.);
	float search_radius = .06;
	float particle_mass = 1.;
	float step_size_multiplier = 0.25;
	float isovalue = 0.01;
	
	// Diffuse particle simulation parameters for falling cube
	DiffuseParameters *dp = new DiffuseParameters(fp->h, 0.5, 0.5, 200, 60, 13);
	dp->t_k_min = 3;
	dp->t_k_max = 40;
	ParentFluid *pf = new ParentFluid(4, 4, 4, 10, 60, 2.5, dp);
	pf->fluid = f;
	pf->fp = fp;

	// Simulate all frames
	vector<Particle> divided_particles_4;
	
	for (int frame = 0; frame < fp->total_time * fp->fps; frame++) {
		cout << "Starting on frame #: " + to_string(frame) << endl;
		
		if (frame % downsample_rate == 0) {
			// Create a deep copy of all the particles and divide positions to keep everything within (0, 0, 0) and (2, 2, 2)
			divided_particles_4 = f->particles;
			#pragma omp parallel for
			for (auto p = begin(divided_particles_4); p != end(divided_particles_4); p++) {
				p->position /= 4;
			}
			cout << "Done Splitting on frame #: " + to_string(frame) << endl;

			if (march) {
				cout << "Marching on fluid" << endl;
				// Perform marching cubes and generate .obj file
				marchingCube* m = new marchingCube(bDim, partDim, divided_particles_4, f->map, fp->h, search_radius,
												   particle_mass, fp->density, isovalue, step_size_multiplier, 0.01);
				m->main_March("Frame-" + to_string(frame) + ".obj");
				delete m;
			}
			
			write_ply(divided_particles_4, "Fluid-" + to_string(frame) + ".ply");
			
			divided_particles_4.clear();
			for (auto p = begin(*pf->diffuse_particles); p != end(*pf->diffuse_particles); p++) {
				divided_particles_4.emplace_back(Particle((*p)->position / 4));
			}
			
			write_ply(divided_particles_4, "Diffuse-" + to_string(frame) + ".ply");
			
			cout << "" << endl;
			cout << "Generated frame #" + to_string(frame) << endl;
		} else {
			cout << "Frame #: " + to_string(frame) + " skipped." << endl;
		}

		// Simulate particle positions for next time step
		// Simulate 1 step
		f->simulate(fp, accel, &collision);
		pf->simulate_step(accel, &collision);
		
		cout << "Simulated frame " + to_string(frame) << endl;
		cout << " " << endl;
	}
	
	if (march) {
		remove("v.obj");
		remove("n.obj");
		remove("f.obj");
	}
	
	return 0;
}
