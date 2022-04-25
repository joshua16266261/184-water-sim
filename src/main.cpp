#include "fluid.h"
#include "particle.h"
#include "marchingCube.h"

using namespace std;

int main(int argc, char **argv) {
	// Just testing build
	Fluid *f = new Fluid(10, 10, 10, 10, 10, 10);
	f->buildGrid();
	for (auto p = begin(f->particles); p != end(f->particles); p++) {
		std::cout << p->position << '\n';
	}

	// Call the constructor to creat marchingCube object
	marchingCube* m = new marchingCube();
	m->main_March();

  return 0;
}
