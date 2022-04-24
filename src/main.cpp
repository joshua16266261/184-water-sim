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
	Fluid *f = new Fluid(5, 5, 5, 10, 10, 10);
	
	FluidParameters *fp = new FluidParameters(EPS_F, 5, 20, 10);
	fp->h = 0.5;
	
	Vector3D g = Vector3D(0, 0, -9.81);
	vector<Vector3D> accel = vector<Vector3D>{g};
	
	Plane *floor = new Plane(Vector3D(0, 0, -1), Vector3D(0, 0, 1), 0);
	Plane *left_wall = new Plane(Vector3D(-1, 0, 0), Vector3D(1, 0, 0), 0);
	Plane *right_wall = new Plane(Vector3D(6, 0, 0), Vector3D(-1, 0, 0), 0);
	Plane *front_wall = new Plane(Vector3D(0, -1, 0), Vector3D(0, 1, 0), 0);
	Plane *back_wall = new Plane(Vector3D(0, 6, 0), Vector3D(0, -1, 0), 0);
	vector<CollisionObject *> collision = vector<CollisionObject *>{floor, left_wall, right_wall, front_wall, back_wall};
	
	
	for (int frame = 0; frame < fp->total_time * fp->fps; frame++) {
		std::cout << frame << '\n';
		f->simulate(fp, accel, &collision);
		write_pos_to_file(f, "floor " + to_string(frame) + ".txt");
	}

  return 0;
}
