#include <iostream>
#include <math.h>
#include <random>
#include <vector>
#include <string>
#include <fstream>

#pragma GCC target ("avx2")
#pragma GCC optimization ("O3")
#pragma GCC optimization ("unroll-loops")

#include "fluid.h"
#include "particle.h"
#include "collision/plane.h"
#include "collision/sphere.h"

#include "marchingCube.h"

using namespace std;

// Default Constructor just for build TEST not used
marchingCube::marchingCube(){}

// Marching Cube object constructor
marchingCube::marchingCube(Vector3D box_dimensions, Vector3D particle_dimensions, vector<Particle> particles,
						   unordered_map<string, vector<Particle*>*> hash_to_particles, float h,
						   float search_radius, float particle_mass, float density, float isovalue,
						   float step_size_multiplier, float the_box_hash_size) {
	
    init(box_dimensions, particle_dimensions, particles, hash_to_particles, h, search_radius, particle_mass, density, isovalue, step_size_multiplier, the_box_hash_size);
}


// Initializer function
void marchingCube::init(Vector3D box_dimensions, Vector3D particle_dimensions,
						vector<Particle> particles, unordered_map<string, vector<Particle*>*> hash_to_particles,
						float h, float search_radius, float particle_mass,
						float density, float isovalue, float step_size_multiplier, float the_box_hash_size) {
	
	m_box_dimensions = box_dimensions;
	m_particle_dimensions = particle_dimensions;
	m_particles = particles;
	m_h = h;
	m_isovalue = isovalue;
	m_search_radius = search_radius;
	m_particle_mass = particle_mass;
	m_density = density;
	m_step_size_multiplier = step_size_multiplier;
	// This variable set the size of the box hash used to speedup
	box_hash_size = the_box_hash_size;

	m_hash_to_particles.clear();

	// Build the Initial Hash Table
	for (auto p = begin(particles); p != end(particles); p++) {
		string key = hash_position(p->position, box_hash_size);
		if (m_hash_to_particles.count(key) > 0) {
			m_hash_to_particles[key]->emplace_back(&*p);
		} else {
			vector<Particle*>* v = new vector<Particle*>;
			m_hash_to_particles[key] = v;
			v->emplace_back(&*p);
		}
	}


	// Sets the unit lwh
	m_unit_dimensions = Vector3D(box_dimensions.x * m_step_size_multiplier / (m_particle_dimensions.x - 1),
								 box_dimensions.y * m_step_size_multiplier / (m_particle_dimensions.y - 1),
								 box_dimensions.z * m_step_size_multiplier / (m_particle_dimensions.z - 1));

	// TO CHANGE THE STEP SIZE WHILE KEEPING THE TOTAL LWH,
	// WE JUST NEED TO MULTIPLY x,y,z's numerator by a number (#)
	// and divide the depth,row,col < variable by a number (#)
	// for example, depth < m_part_dim.z / [(5)]
	//              z = depth * box_dim.z * [(5)] / (m_part_dim.z - 1)
	// where we want to increase the step size 5 times
	
	int max_depth = ceil(m_particle_dimensions.z / m_step_size_multiplier);
	int max_row = ceil(m_particle_dimensions.y / m_step_size_multiplier);
	int max_col = ceil(m_particle_dimensions.x / m_step_size_multiplier);
	
	isovalues.clear();

	// Creating the cubes takes the longest time
	// Everything else is fast
	// have to expand the box to fit the particles (like -2, -2, -2 to 2, 2, 2, or wherever the walls are etc.)
	Vector3D index;
//	for (int depth = 0; depth <= max_depth; depth++) {
//		index.z = depth * m_unit_dimensions.z - 0.25;
//		for (int row = 0; row <= max_row; row++) {
//			index.y = row * m_unit_dimensions.y - 0.25;
//			for (int col = 0; col <= max_col; col++) {
//				index.x = col * m_unit_dimensions.x - 0.25;
//				// Will create a empty cube pass in top left positional index of the cube
//				// Iterate over cube creation and then just emplace back after you fill it up
//				Cube marchCube = Cube();
//				createCube(marchCube, index);
//			}
//		}
//	}
	for (int col = 0; col <= max_col; col++) {
		index.x = col * m_unit_dimensions.x - 0.25;
		// Erase isovalues 2 rows to the left
		isovalues.erase(index.x - 2 * m_unit_dimensions.x);
		for (int row = 0; row <= max_row; row++) {
			index.y = row * m_unit_dimensions.y - 0.25;
			// Erase isovalues 2 cols below
			isovalues[index.x].erase(index.y - 2 * m_unit_dimensions.y);
			for (int depth = 0; depth <= max_depth; depth++) {
				index.z = depth * m_unit_dimensions.z - 0.25;
				// Erase isovalues 2 layers under
				isovalues[index.x][index.y].erase(index.z - 2 * m_unit_dimensions.z);
				// Will create a empty cube pass in top left positional index of the cube
				createCube(Cube(), index);
			}
		}
	}
}

void marchingCube::main_March(string filename) {
    // we have the marching cube vector now just iterate over the list and
    for (auto q = begin(cube_Vector); q != end(cube_Vector); q++) {
        Polygonise(*q, m_isovalue);
    }
    // will call the file.obj function in the main function btw and not here
    triToObj(filename);
	cout << "Written OBJ File" << endl;
}



// r is the vector from particle to position p
float marchingCube::isotropic_kernel(Vector3D r, float h) {
    float constant = 315 / (64 * PI * pow(h, 9));

    if (r.norm() >= 0 && r.norm() <= h) {
        return constant * pow(h * h - r.norm() * r.norm(), 3);
    } else {
        return 0.0;
    }
}

// Compute the density of a pos (will be used as the isovalue)
float marchingCube::isovalue(Vector3D pos, float h) {
    float rho = 0.0;
    // Hash the position vector
    // Find particles that are neighboring
    // Check for particles that are in specific radius
    
	// Need to check surrounding cube for edge case
//	vector<Vector3D> cube_off = { Vector3D(0, 0, 0), Vector3D(0, 0, 1), Vector3D(0, 0, -1), Vector3D(0, 1, 0),
//									Vector3D(0, 1, 1), Vector3D(0, 1, -1), Vector3D(0, -1, 0), Vector3D(0, -1, 1), Vector3D(0, -1, -1),
//								  Vector3D(1, 0, 0), Vector3D(1, 0, 1), Vector3D(1, 0, -1), Vector3D(1, 1, 0), Vector3D(1, 1, 1), Vector3D(1, 1, -1), Vector3D(1, -1, 0), Vector3D(1, -1, 1), Vector3D(1, -1, -1),
//								  Vector3D(-1, 0, 0), Vector3D(-1, 0, 1), Vector3D(-1, 0, -1), Vector3D(-1, 1, 0), Vector3D(-1, 1, 1), Vector3D(-1, 1, -1), Vector3D(-1, -1, 0), Vector3D(-1, -1, 1), Vector3D(-1, -1, -1)};
	Vector3D update_pos;
	for (int i = -1; i <= 1; i++) {
		update_pos.x = pos.x + box_hash_size * i;
		for (int j = -1; j <= 1; j++) {
			update_pos.y = pos.y + box_hash_size * j;
			for (int k = -1; k <= 1; k++) {
				update_pos.z = pos.z + box_hash_size * k;
				// Compute the Vortex of with the offset applied
				string vortex_key = hash_position(update_pos, box_hash_size);
				auto c = m_hash_to_particles[vortex_key];
				if (c) {
					for (auto q = begin(*c); q != end(*c); q++) {
						// Check if particle is inside the radius
						if ((pos - (*q)->position).norm() <= m_search_radius) {
							rho += m_particle_mass * isotropic_kernel(pos - (*q)->position, h);
						}
					}
				}
			}
		}
	}
//	for (auto off : cube_off) {
//		// Compute the Vortex of with the offset applied
//		Vector3D update_pos = Vector3D(pos.x + box_hash_size * off.x, pos.y + box_hash_size * off.y, pos.z + box_hash_size * off.z);
//		string vortex_key = hash_position(update_pos, box_hash_size);
//		auto c = m_hash_to_particles[vortex_key];
//		if (c == NULL) {
//			continue;
//		}
//		for (auto q = begin(*c); q != end(*c); q++) {
//			// Check if particle is inside the radius
//			if ((pos - (*q)->position).norm() <= m_search_radius) {
//				rho += m_particle_mass * isotropic_kernel(pos - (*q)->position, h);
//			}
//		}
//	}

    return rho;
}

string marchingCube::hash_position(Vector3D pos, float h) {
    // Hash a 3D position into a unique Vector3D identifier that represents membership in some 3D box volume.
    // Assume 30x30x30 cube and divide by h=5
    int x_box = floor(pos.x / h);
    int y_box = floor(pos.y / h);
    int z_box = floor(pos.z / h);

    string s = to_string(x_box);
    s.append(":");
    s.append(to_string(y_box));
    s.append(":");
    s.append(to_string(z_box));

    return s;
}



// Create the Marching Cube Grid (each sub-cube within the whole space)
// Inputs: the index of the cube in all directions (i.e. 1 in X, 4 in Y, and 2 in Z)
void marchingCube::createCube(Cube cube, Vector3D index) {
    // Cube placement in verticies is based on http://paulbourke.net/geometry/polygonise/polygonise1.gif
    // Configuration:
    // bottom-backward-left = 0, bottom-back-right = 1, bottom-front-right = 2,
    // bottom-front-left = 3, top-back-left = 4, top-back-right = 5,
    // top-front-right = 6, top-front-left = 7 
    // The top-back-left index we pass in is vertex 4 

    // Assign the verticies
    float unit_X = m_unit_dimensions.x;
    float unit_Y = m_unit_dimensions.y;
    float unit_Z = m_unit_dimensions.z;

    // The way I did this was that to go from one vertex to another of the cube you just add 1 unit length
    // Can do this visually, checking the picture aswell to do re-verify this works
    // NOTE TO SELF FOR SANITY: IF YOU USE A CLASS VAR NAME DO NOT RE-INITIALIZE IT OR ELSE IT WILL GIVE YOU 0  AAAAAAAAAAAAAAA
    // FIXED ACCORING TO https://i.stack.imgur.com/oLUUQ.png
    cube.vertices[0] = Vector3D(index.x, index.y, index.z - unit_Z);
    cube.vertices[1] = Vector3D(index.x + unit_X, index.y, index.z - unit_Z);
    cube.vertices[2] = Vector3D(index.x + unit_X, index.y + unit_Y, index.z - unit_Z);
    cube.vertices[3] = Vector3D(index.x, index.y + unit_Y, index.z - unit_Z);
    cube.vertices[4] = index;
    cube.vertices[5] = Vector3D(index.x + unit_X, index.y, index.z);
    cube.vertices[6] = Vector3D(index.x + unit_X, index.y + unit_Y, index.z);
    cube.vertices[7] = Vector3D(index.x, index.y + unit_Y, index.z);

    // Normals
    cube.normals[0] = cube.vertices[0] - (cube.vertices[0] - cube.vertices[6]).unit();
    cube.normals[1] = cube.vertices[1] - (cube.vertices[1] - cube.vertices[7]).unit();
    cube.normals[2] = cube.vertices[2] - (cube.vertices[2] - cube.vertices[4]).unit();
    cube.normals[3] = cube.vertices[3] - (cube.vertices[3] - cube.vertices[5]).unit();
    cube.normals[4] = cube.vertices[4] + (cube.vertices[4] - cube.vertices[2]).unit();
    cube.normals[5] = cube.vertices[5] + (cube.vertices[5] - cube.vertices[3]).unit();
    cube.normals[6] = cube.vertices[6] + (cube.vertices[6] - cube.vertices[0]).unit();
    cube.normals[7] = cube.vertices[7] + (cube.vertices[7] - cube.vertices[1]).unit();
	
    // Assign the isovalue for each cube vertex
	Vector3D vertex_pos;
    for (int i = 0; i < 8; i++) {
		vertex_pos = cube.vertices[i];
		if (isovalues.count(vertex_pos.x) > 0) {
			if (isovalues[vertex_pos.x].count(vertex_pos.y) > 0) {
				if (isovalues[vertex_pos.x][vertex_pos.y].count(vertex_pos.z) == 0) {
					isovalues[vertex_pos.x][vertex_pos.y][vertex_pos.z] = isovalue(vertex_pos, m_h);
				}
			} else {
				isovalues[vertex_pos.x][vertex_pos.y] = unordered_map<float, float>();
				isovalues[vertex_pos.x][vertex_pos.y][vertex_pos.z] = isovalue(vertex_pos, m_h);
			}
		} else {
			isovalues[vertex_pos.x] = unordered_map<float, unordered_map<float, float>>();
			isovalues[vertex_pos.x][vertex_pos.y] = unordered_map<float, float>();
			isovalues[vertex_pos.x][vertex_pos.y][vertex_pos.z] = isovalue(vertex_pos, m_h);
		}
		
		cube.isovalues[i] = isovalues[vertex_pos.x][vertex_pos.y][vertex_pos.z];
    }
	cube_Vector.emplace_back(cube);
}


/*
   Given a grid cell (cube in our case) and an isolevel, calculate the triangular
   facets required to represent the isosurface through the cell (or cube).
   Return the number of triangular facets, the array "triangles"
   will be loaded up with the vertices at most 5 triangular facets.
    0 will be returned if the grid cell is either totally above
   of totally below the isolevel.
*/
// http://paulbourke.net/geometry/polygonise/ is origin of code
int marchingCube::Polygonise(Cube cube, double isolevel) {
    int i, ntriang;
    int cubeindex;
    Vector3D vertlist[12];
    Vector3D normlist[12];

    /*
     Determine the index into the edge table which
     tells us which vertices are inside of the surface
    */
    cubeindex = 0;

    if (cube.isovalues[0] < isolevel) cubeindex |= 1;
    if (cube.isovalues[1] < isolevel) cubeindex |= 2;
    if (cube.isovalues[2] < isolevel) cubeindex |= 4;
    if (cube.isovalues[3] < isolevel) cubeindex |= 8;
    if (cube.isovalues[4] < isolevel) cubeindex |= 16;
    if (cube.isovalues[5] < isolevel) cubeindex |= 32;
    if (cube.isovalues[6] < isolevel) cubeindex |= 64;
    if (cube.isovalues[7] < isolevel) cubeindex |= 128;


    /* Cube is entirely in/out of the surface */
    if (edgeTable[cubeindex] == 0)
        return 0;

    /* Find the vertices where the surface intersects the cube */
    // FOR THIS FUNCTION I ADDED IN GETTING THE NORMALS TOO //
    if (edgeTable[cubeindex] & 1)
        vertlist[0] =
        VertexInterp(isolevel, cube.vertices[0], cube.vertices[1], cube.isovalues[0], cube.isovalues[1]);
        normlist[0] = VertexInterpNormals(isolevel, cube.normals[0], cube.normals[1], cube.isovalues[0], cube.isovalues[1]);

    if (edgeTable[cubeindex] & 2)
        vertlist[1] =
        VertexInterp(isolevel, cube.vertices[1], cube.vertices[2], cube.isovalues[1], cube.isovalues[2]);
        normlist[1] = VertexInterpNormals(isolevel, cube.normals[1], cube.normals[2], cube.isovalues[1], cube.isovalues[2]);

    if (edgeTable[cubeindex] & 4)
        vertlist[2] =
        VertexInterp(isolevel, cube.vertices[2], cube.vertices[3], cube.isovalues[2], cube.isovalues[3]);
        normlist[2] = VertexInterpNormals(isolevel, cube.normals[2], cube.normals[3], cube.isovalues[2], cube.isovalues[3]);

    if (edgeTable[cubeindex] & 8)
        vertlist[3] =
        VertexInterp(isolevel, cube.vertices[3], cube.vertices[0], cube.isovalues[3], cube.isovalues[0]);
        normlist[3] = VertexInterpNormals(isolevel, cube.normals[3], cube.normals[0], cube.isovalues[3], cube.isovalues[0]);

    if (edgeTable[cubeindex] & 16)
        vertlist[4] =
        VertexInterp(isolevel, cube.vertices[4], cube.vertices[5], cube.isovalues[4], cube.isovalues[5]);
        normlist[4] = VertexInterpNormals(isolevel, cube.normals[4], cube.normals[5], cube.isovalues[4], cube.isovalues[5]);

    if (edgeTable[cubeindex] & 32)
        vertlist[5] =
        VertexInterp(isolevel, cube.vertices[5], cube.vertices[6], cube.isovalues[5], cube.isovalues[6]);
        normlist[5] = VertexInterpNormals(isolevel, cube.normals[5], cube.normals[6], cube.isovalues[5], cube.isovalues[6]);

    if (edgeTable[cubeindex] & 64)
        vertlist[6] =
        VertexInterp(isolevel, cube.vertices[6], cube.vertices[7], cube.isovalues[6], cube.isovalues[7]);
        normlist[6] = VertexInterpNormals(isolevel, cube.normals[6], cube.normals[7], cube.isovalues[6], cube.isovalues[7]);


    if (edgeTable[cubeindex] & 128)
        vertlist[7] =
        VertexInterp(isolevel, cube.vertices[7], cube.vertices[4], cube.isovalues[7], cube.isovalues[4]);
        normlist[7] = VertexInterpNormals(isolevel, cube.normals[7], cube.normals[4], cube.isovalues[7], cube.isovalues[4]);


    if (edgeTable[cubeindex] & 256)
        vertlist[8] =
        VertexInterp(isolevel, cube.vertices[0], cube.vertices[4], cube.isovalues[0], cube.isovalues[4]);
        normlist[8] = VertexInterpNormals(isolevel, cube.normals[0], cube.normals[4], cube.isovalues[0], cube.isovalues[4]);


    if (edgeTable[cubeindex] & 512)
        vertlist[9] =
        VertexInterp(isolevel, cube.vertices[1], cube.vertices[5], cube.isovalues[1], cube.isovalues[5]);
        normlist[9] = VertexInterpNormals(isolevel, cube.normals[1], cube.normals[5], cube.isovalues[1], cube.isovalues[5]);


    if (edgeTable[cubeindex] & 1024)
        vertlist[10] =
        VertexInterp(isolevel, cube.vertices[2], cube.vertices[6], cube.isovalues[2], cube.isovalues[6]);
        normlist[10] = VertexInterpNormals(isolevel, cube.normals[2], cube.normals[6], cube.isovalues[2], cube.isovalues[6]);


    if (edgeTable[cubeindex] & 2048)
        vertlist[11] =
        VertexInterp(isolevel, cube.vertices[3], cube.vertices[7], cube.isovalues[3], cube.isovalues[7]);
        normlist[11] = VertexInterpNormals(isolevel, cube.normals[3], cube.normals[7], cube.isovalues[3], cube.isovalues[7]);


    /* Create the triangle */
    // Here we have a passed in list of triangles, 
    // but later on we should place the list of triangles we passed in onto a triangle vector
    ntriang = 0;
    
    for (i = 0; triTable[cubeindex][i] != -1; i += 3) {
        // We will initate a new triangle for each look and then will emplace it back every time
        newTriangle triangles = newTriangle();
        triangles.coordinates[0] = vertlist[triTable[cubeindex][i]];
        triangles.coordinates[1] = vertlist[triTable[cubeindex][i + 1]];
        triangles.coordinates[2] = vertlist[triTable[cubeindex][i + 2]];

        triangles.normal[0] = normlist[triTable[cubeindex][i]];
        triangles.normal[1] = normlist[triTable[cubeindex][i + 1]];
        triangles.normal[2] = normlist[triTable[cubeindex][i + 2]];

        // THIS IS THE VECTOR WE GET ALL OUR TRIANGLES WE WANT TO RASTERIZE INTO 3D FROM
        tri_Vector.emplace_back(triangles);
        ntriang++;
    }

    return ntriang;
}


/*
   Linearly interpolate the position where an isosurface cuts
   an edge between two vertices, each with their own scalar value
*/
// http://paulbourke.net/geometry/polygonise/ is origin of code. This code is modified

Vector3D VertexInterp(double isolevel, Vector3D p1, Vector3D p2, double valp1, double valp2) {
    double mu;
    Vector3D p;
    if (abs(isolevel - valp1) < 0.00001)
        return p1;
    if (abs(isolevel - valp2) < 0.00001)
        return p2;
    if (abs(valp1 - valp2) < 0.00001)
        return p1;
    mu = (isolevel - valp1) / (valp2 - valp1);
    p.x = p1.x + mu * (p2.x - p1.x);
    p.y = p1.y + mu * (p2.y - p1.y);
    p.z = p1.z + mu * (p2.z - p1.z);

    // We may also need to have it so we return p's norm, so maybe return something else instead of a 3D triangle
    return p;
}

Vector3D VertexInterpNormals(double isolevel, Vector3D n1, Vector3D n2, double valp1, double valp2) {
    double mu;
    Vector3D n;
    if (abs(isolevel - valp1) < 0.00001)
        return n1;
    if (abs(isolevel - valp2) < 0.00001)
        return n2;
    if (abs(valp1 - valp2) < 0.00001)
        return n1;
    mu = (isolevel - valp1) / (valp2 - valp1);
    n.x = n1.x + mu * (n2.x - n1.x);
    n.y = n1.y + mu * (n2.y - n1.y);
    n.z = n1.z + mu * (n2.z - n1.z);

    // We may also need to have it so we return p's norm, so maybe return something else instead of a 3D triangle
    return n;
}

void marchingCube::triToObj(string fName) {
    ofstream ofile;
    ofile.open(fName);
    if (ofile.is_open()) {
        cout << "Open Success" << "\n";
    }
	Vector3D v;
    for (auto& tri : tri_Vector) {
        // Add Vertices to the OBJ file
        for (int i = 0; i < 3; i++) {
            v = tri.coordinates[i];
            string pos = "v " + to_string(v.x) + " " + to_string(v.y) + " " + to_string(v.z);
            ofile << pos << "\n";
        }
    }

    for (auto& tri : tri_Vector) {
        // Add Normals to the OBJ file
        for (int i = 0; i < 3; i++) {
            v = tri.normal[i];
            string norm = "vn " + to_string(v.x) + " " + to_string(v.y) + " " + to_string(v.z);
            ofile << norm << "\n";
        }
    }

    // Add the faces
    int idx = 1;
    for (auto& tri : tri_Vector) {
        string face = "f " + to_string(idx) + "//" + to_string(idx) + " " + to_string(idx + 1) + "//" + to_string(idx + 1) + " " + to_string(idx + 2) + "//" + to_string(idx + 2);
        ofile << face << "\n";
        idx += 3;
    }
    ofile.close();
}
