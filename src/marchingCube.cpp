#include <iostream>
#include <math.h>
#include <random>
#include <vector>
#include <string>

#include "fluid.h"
#include "particle.h"
#include "collision/plane.h"
#include "collision/sphere.h"

#include "marchingCube.h"

// Main grid starts from Vector3D(0,0,0) and ends at Vector3D(box_dimensions.xyz)
// The list of particles just incase 
// The hashmap of string to list particles within the same box
// Cube dimensions we want: ex: box_dim = (9,9,9), cube_dim = (3,3,3) so (27 in total since 9/3 ^ 3)
// search_radius, mass, h, and other variables

void marchingCube::init(Vector3D box_dimensions, Vector3D unit_cube,
 vector<Particle> particles, unordered_map<string, vector<Particle*>*> hash_to_particles,
 float h, float search_radius, float particle_mass, 
 float density, float isovalue) {

// init variables to use in this files functions
m_box_dimensions = box_dimensions;
m_unit_cube = unit_cube;
m_particles = particles;
m_hash_to_particles = hash_to_particles;
m_h = h;
m_isovalue = isovalue;
m_search_radius = search_radius;
m_particle_mass = particle_mass;
m_density = density;

// RESET vector of triangles and then re-make it
// do we need to do this????
// delete(tri_vector);
// tri_Vector = new Vector<Triangle>();
///
vector<newTriangle> tri_Vector = vector<newTriangle>();

// Slice up the main box into cubes
// we need to use the create cube function to input 
// and then we just throw our cube into a vector of cubes

// The number of times we iterate through each dimension (l,w,h)
Vector3D iter_dimensions = Vector3D(ceil(box_dimensions.x / m_unit_cube.x), ceil(box_dimensions.y / m_unit_cube.y),
    ceil(box_dimensions.z / m_unit_cube.z));

for (int i = 0; i < iter_dimensions.x; i++) {
    for (int j = 0; j < iter_dimensions.y; j++) {
        for (int k = 0; k < iter_dimensions.z; k++) {
            Cube &marchCube = Cube();

            // Will create a empty cube pass in top left positional index of the cube 
            // Iterate over cube creation and then just emplace back after you fill it up
            Vector3D index = Vector3D(i, j, k);
            createCube(marchCube, index);
            cube_Vector.emplace_back(marchCube);

        }
    }
}
}

// r is the vector from particle to position p
float marchingCube::isotropic_kernel(Vector3D r, float h) {
    float constant = 315 / (64 * PI * pow(h, 9));
    if (r.norm() >= 0 && r.norm() <= h) {
        return constant * pow(h * h - r.norm() * r.norm());
    } else {
        return 0;
    }
}

// Compute the density of a pos (will be used as the isovalue)
float marchingCube::isovalue(Vector3D pos, float h){
    float rho = 0;
    // Hash the position vector
    // Find particles that are neighboring
    // Check for particles that are in specific radius
    string vortex_key = hash_position(pos, h);
    for (auto q = begin(*(hash_to_particles[vortex_key])); q != end(*(hash_to_particles[vortex_key])); q++) {
        if ((pos - (*q)->position).norm() <= search_radius) {
            rho += mass * isotropic_kernel(pos - (q*)->pos, h);
        }
    }
    return rho;
}

// Create the Marching Cube Grid (each sub-cube within the whole space)
// Inputs: the index of the cube in all directions (i.e. 1 in X, 4 in Y, and 2 in Z)
void marchingCube::createCube(Cube &cube, Vector3D index) {
    // Assign all the vertices within a cube to its approperiate (x, y, z) positions
    // Assign all approperiate isovalues to each vertex
    // Then in the init() function, push the cube into our cube_vector

    // Cube placement in verticies is based on http://paulbourke.net/geometry/polygonise/polygonise1.gif
    // 
    // Configuration:
    // bottom-backward-left = 0, bottom-back-right = 1, bottom-front-right = 2,
    // bottom-front-left = 3, top-back-left = 4, top-back-right = 5,
    // top-front-right = 6, top-front-left = 7 
    
    // The top-back-left index we pass in is vertex 4 
    // We may need to have edge cases where we surpass the original whole box with our marching cubes because we ceil in iter_dimensions?

    
    // Assign the verticies
    float unit_L = m_unit_cube.x;
    float unit_W = m_unit_cube.y;
    float unit_H = m_unit_cube.z;

    // The way I did this was that to go from one verticie to another of the cube you just add 1 unit length
    // Can do this visually, checking the picture aswell to do re-verify this works
    cube.vertices[0] = Vector3D(index.x, index.y, index.z + unit_H);
    cube.vertices[1] = Vector3D(index.x + unit_L, index.y, index.z + unit_H);
    cube.vertices[2] = Vector3D(index.x + unit_L, index.y + unit_W, index.z + unit_H);
    cube.vertices[3] = Vector3D(index.x, index.y + unit_W, index.z + unit_H);
    cube.vertices[4] = index;
    cube.vertices[5] = Vector3D(index.x + unit_L, index.y, index.z);
    cube.vertices[6] = Vector3D(index.x + unit_L, index.y + unit_W, index.z);
    cube.vertices[7] = Vector3D(index.x, index.y + unit_W, index.z);
    
    // Assign the isovalue for each cube verticie 
    for (int i = 0; i < 8; i++) {
        cube.isovalues[i] = isovalue(cube.vertices[i], m_h);
    }

    return;
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
// We just need to modify it to suit our code or modify our code to suit it

// We replace grid with cube and 
// The triangles is just a empty list of triangles you pass in to have it filled,
// THIS IS WHERE WE FILL IN THE TRIANGLES WHERE WE WILL USE IT TO FILL IN THE RASTERIZER
// This code is modified btw
int marchingCube::Polygonise(Cube cube, double isolevel, newTriangle* triangles) {
    int i, ntriang;
    int cubeindex;
    Vector3D vertlist[12];

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
    if (edgeTable[cubeindex] & 1)
        vertlist[0] =
        VertexInterp(isolevel, cube.vertices[0], cube.vertices[1], cube.isovalues[0], cube.isovalues[1]);
    if (edgeTable[cubeindex] & 2)
        vertlist[1] =
        VertexInterp(isolevel, cube.vertices[1], cube.vertices[2], cube.isovalues[1], cube.isovalues[2]);
    if (edgeTable[cubeindex] & 4)
        vertlist[2] =
        VertexInterp(isolevel, cube.vertices[2], cube.vertices[3], cube.isovalues[2], cube.isovalues[3]);
    if (edgeTable[cubeindex] & 8)
        vertlist[3] =
        VertexInterp(isolevel, cube.vertices[3], cube.vertices[0], cube.isovalues[3], cube.isovalues[0]);
    if (edgeTable[cubeindex] & 16)
        vertlist[4] =
        VertexInterp(isolevel, cube.vertices[4], cube.vertices[5], cube.isovalues[4], cube.isovalues[5]);
    if (edgeTable[cubeindex] & 32)
        vertlist[5] =
        VertexInterp(isolevel, cube.vertices[5], cube.vertices[6], cube.isovalues[5], cube.isovalues[6]);
    if (edgeTable[cubeindex] & 64)
        vertlist[6] =
        VertexInterp(isolevel, cube.vertices[6], cube.vertices[7], cube.isovalues[6], cube.isovalues[7]);
    if (edgeTable[cubeindex] & 128)
        vertlist[7] =
        VertexInterp(isolevel, cube.vertices[7], cube.vertices[4], cube.isovalues[7], cube.isovalues[4]);
    if (edgeTable[cubeindex] & 256)
        vertlist[8] =
        VertexInterp(isolevel, cube.vertices[0], cube.vertices[4], cube.isovalues[0], cube.isovalues[4]);
    if (edgeTable[cubeindex] & 512)
        vertlist[9] =
        VertexInterp(isolevel, cube.vertices[1], cube.vertices[5], cube.isovalues[1], cube.isovalues[5]);
    if (edgeTable[cubeindex] & 1024)
        vertlist[10] =
        VertexInterp(isolevel, cube.vertices[2], cube.vertices[6], cube.isovalues[2], cube.isovalues[6]);
    if (edgeTable[cubeindex] & 2048)
        vertlist[11] =
        VertexInterp(isolevel, cube.vertices[3], cube.vertices[7], cube.isovalues[3], cube.isovalues[7]);

    /* Create the triangle */
    // Here we have a passed in list of triangles, 
    // but later on we should place the list of triangles we passed in onto a triangle vector
    //
    ntriang = 0;
    for (i = 0; triTable[cubeindex][i] != -1; i += 3) {
        triangles[ntriang].coordinates[0] = vertlist[triTable[cubeindex][i]];
        triangles[ntriang].coordinates[1] = vertlist[triTable[cubeindex][i + 1]];
        triangles[ntriang].coordinates[2] = vertlist[triTable[cubeindex][i + 2]];

        // THIS IS THE VECTOR WE GET ALL OUR TRIANGLES WE WANT TO RASTERIZE INTO 3D FROM
        tri_Vector.emplace_back(triangles[ntriang]);
        ntriang++;
    }

    return ntriang;
}


/*
   Linearly interpolate the position where an isosurface cuts
   an edge between two vertices, each with their own scalar value
*/
// http://paulbourke.net/geometry/polygonise/ is origin of code. This code is modified

Vector3D marchingCube::VertexInterp(double isolevel, Vector3D p1, Vector3D p2, double valp1, double valp2) {
    double mu;
    Vector3D p;
    if (abs(isolevel - valp1) < 0.00001)
        return(p1);
    if (abs(isolevel - valp2) < 0.00001)
        return(p2);
    if (abs(valp1 - valp2) < 0.00001)
        return(p1);
    mu = (isolevel - valp1) / (valp2 - valp1);
    p.x = p1.x + mu * (p2.x - p1.x);
    p.y = p1.y + mu * (p2.y - p1.y);
    p.z = p1.z + mu * (p2.z - p1.z);

    // We may also need to have it so we return p's norm, so maybe return something else instead of a 3D triangle

    return p;
}

