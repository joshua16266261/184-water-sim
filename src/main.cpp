#include "fluid.h"
#include "particle.h"
#include "collision/plane.h"
#include <vector>

#include <fstream>

using namespace std;

void write_pos_to_file(Fluid *f, string filename) {
	string s = "";
	for (auto p = begin(f->particles); p != end(f->particles); p++) {
		s.append("(" + to_string(p->position[0]) + "," + to_string(p->position[1]) + "," + to_string(p->position[2]) + ")" + '\n');
	}
	ofstream file(filename);
	file << s;
	file.close();
}

int main(int argc, char **argv) {
	Fluid *f = new Fluid(2, 2, 2, 20, 20, 20);
	f->frame = 0;
	
	FluidParameters *fp = new FluidParameters(EPS_F, 2, 60, 5);
	fp->h = 0.15;
	fp->delta_q = 0.02 * fp->h * Vector3D(1, 1, 1);
	fp->density = 1000;
	fp->k = 0.001;
	fp->c = 0.00025;
	fp->vorticity_eps = 0.005;
	fp->relaxation = 5000;
	
	Vector3D g = Vector3D(0, 0, -9.81);
	vector<Vector3D> accel = vector<Vector3D>{g};
	
	Plane *floor = new Plane(Vector3D(0, 0, -1), Vector3D(0, 0, 1), 0);
	Plane *left_wall = new Plane(Vector3D(-1, 0, 0), Vector3D(1, 0, 0), 0);
	Plane *right_wall = new Plane(Vector3D(3, 0, 0), Vector3D(-1, 0, 0), 0);
	Plane *front_wall = new Plane(Vector3D(0, -1, 0), Vector3D(0, 1, 0), 0);
	Plane *back_wall = new Plane(Vector3D(0, 3, 0), Vector3D(0, -1, 0), 0);
	vector<CollisionObject *> collision = vector<CollisionObject *>{floor, left_wall, right_wall, front_wall, back_wall};
	
//	int frame = 0;
	for (int frame = 0; frame < fp->total_time * fp->fps; frame++) {
		std::cout << frame << '\n';
		f->simulate(fp, accel, &collision);
		write_pos_to_file(f, "floor " + to_string(frame) + ".txt");
	}

  return 0;
}
