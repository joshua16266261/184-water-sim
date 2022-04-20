#include "fluid.h"
#include "particle.h"
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
	Fluid *f = new Fluid(5, 5, 5, 5, 5, 5);
//	for (auto p = begin(f->particles); p != end(f->particles); p++) {
//		std::cout << p->position << '\n';
//	}
//	
//	std::cout << "----------------------------" << '\n';
//	std::cout << "----------------------------" << '\n';
//	std::cout << "----------------------------" << '\n';
//	std::cout << "----------------------------" << '\n';
//	std::cout << "----------------------------" << '\n';
//	std::cout << "----------------------------" << '\n';
//	std::cout << "----------------------------" << '\n';
//	std::cout << "----------------------------" << '\n';
//	std::cout << "----------------------------" << '\n';
//	std::cout << "----------------------------" << '\n';
//	std::cout << "----------------------------" << '\n';
//	std::cout << "----------------------------" << '\n';
//	std::cout << "----------------------------" << '\n';
//	std::cout << "----------------------------" << '\n';
	
	FluidParameters *fp = new FluidParameters(EPS_F, 1, 20, 10);
	
	Vector3D g = Vector3D(0, 0, -9.81);
	vector<Vector3D> accel = vector<Vector3D>{g};
	
	f->simulate(fp, accel, nullptr);
	
	write_pos_to_file(f, "output.txt");
	
	for (auto p = begin(f->particles); p != end(f->particles); p++) {
		std::cout << p->position << '\n';
	}

  return 0;
}
