#include "fluid.h"
#include "particle.h"
#include <vector>

using namespace std;

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
	
	for (auto p = begin(f->particles); p != end(f->particles); p++) {
		std::cout << p->position << '\n';
	}

  return 0;
}
