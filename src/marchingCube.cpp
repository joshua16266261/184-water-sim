#include <iostream>
#include <math.h>
#include <random>
#include <vector>
#include <string>
#include <fstream>

#include "fluid.h"
#include "particle.h"
#include "collision/plane.h"
#include "collision/sphere.h"

#include "marchingCube.h"

using namespace std;


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

// CHANGE FOR TESTING PURPOSES
m_box_dimensions = box_dimensions;
m_unit_cube = unit_cube;
// m_particles is outside of the object for this function to work
// m_particles = vector<Particle>();
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
vector<Cube> cube_Vector = vector<Cube>();

// Slice up the main box into cubes
// we need to use the create cube function to input 
// and then we just throw our cube into a vector of cubes

//////////////////////
// TEST CODE//
/////////////////////
float x, y, z;
for (int depth = 0; depth < 10; depth++) {
    z = depth * 10 / (10 - 1);
    for (int row = 0; row < 10; row++) {
        y = row * 10 / (10 - 1);
        for (int col = 0; col < 10; col++) {
            x = col * 10 / (10 - 1);
            m_particles.emplace_back(Particle(Vector3D(x, y, z)));
        }
    }
}

// NOW m_particles holdsall the particles



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

void marchingCube::main_March() {
        // we have the marching cube vector now just iterate over the list and 
    for (const Cube &singleCube : cube_Vector) {

        // is this how we init a new array of triangles???
        newTriangle *triangleHolder;
        Polygonise(singleCube, m_isovalue);

    }
    // will call the file.obj function in the main function btw and not here
}



// r is the vector from particle to position p
float marchingCube::isotropic_kernel(Vector3D r, float h) {
    float constant = 315 / (64 * PI * pow(h, 9));
    if (r.norm() >= 0 && r.norm() <= h) {
        return constant * pow(h * h - r.norm() * r.norm(), 3);
    } else {
        return 0;
    }
}

// Compute the density of a pos (will be used as the isovalue)
float marchingCube::isovalue(Vector3D pos, float h) {
    float rho = 0;
    // Hash the position vector
    // Find particles that are neighboring
    // Check for particles that are in specific radius
    string vortex_key = hash_position(pos, h);
    for (auto q = begin(*(m_hash_to_particles[vortex_key])); q != end(*(m_hash_to_particles[vortex_key])); q++) {
        if ((pos - (*q)->position).norm() <= m_search_radius) {
            rho += m_particle_mass * isotropic_kernel(pos - (*q)->position, h);
        }
    }
    return rho;
}

string marchingCube::hash_position(Vector3D pos, float h) {
    // Hash a 3D position into a unique Vector3D identifier that represents membership in some 3D box volume.
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
    

    // The way we calculate the normals is to take a vertex and it's diagonal counterpart and subtract it to get the
    // 3d ray that is orthogonal/diagonal to the vertex and it's 3 edges
    // 
    //      x--------x
    //     /|       /|
    //    x--------x |
    //    | |      | |
    //    | x------|-x  
    //    x--------x
    // 
    // However, another way is to do something with the isovalues
    // but we can try it if this one does not work.
    // 
    cube.normals[0] = cube.vertices[0] - (cube.vertices[0] - cube.vertices[6]).unit();
    cube.normals[1] = cube.vertices[1] - (cube.vertices[1] - cube.vertices[7]).unit();
    cube.normals[2] = cube.vertices[2] - (cube.vertices[2] - cube.vertices[4]).unit();
    cube.normals[3] = cube.vertices[3] - (cube.vertices[3] - cube.vertices[5]).unit();
    cube.normals[4] = cube.vertices[4] - (cube.vertices[4] - cube.vertices[2]).unit();
    cube.normals[5] = cube.vertices[5] - (cube.vertices[5] - cube.vertices[3]).unit();
    cube.normals[6] = cube.vertices[6] - (cube.vertices[6] - cube.vertices[0]).unit();
    cube.normals[7] = cube.vertices[7] - (cube.vertices[7] - cube.vertices[1]).unit();

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
        normlist[10] = VertexInterpNormals(isolevel, cube.normals[3], cube.normals[7], cube.isovalues[3], cube.isovalues[7]);


    /* Create the triangle */
    // Here we have a passed in list of triangles, 
    // but later on we should place the list of triangles we passed in onto a triangle vector
    //
    ntriang = 0;
    for (i = 0; triTable[cubeindex][i] != -1; i += 3) {

        // We will initate a new triangle for each look and then will emplace it back every time
        newTriangle triangles = newTriangle();
        triangles.coordinates[0] = vertlist[triTable[cubeindex][i]];
        triangles.coordinates[1] = vertlist[triTable[cubeindex][i + 1]];
        triangles.coordinates[2] = vertlist[triTable[cubeindex][i + 2]];

        triangles.normal[0] = normlist[triTable[cubeindex][i + 2]];
        triangles.normal[1] = normlist[triTable[cubeindex][i + 2]];
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

Vector3D VertexInterpNormals(double isolevel, Vector3D n1, Vector3D n2, double valp1, double valp2) {
    double mu;
    Vector3D n;
    if (abs(isolevel - valp1) < 0.00001)
        return(n1);
    if (abs(isolevel - valp2) < 0.00001)
        return(n2);
    if (abs(valp1 - valp2) < 0.00001)
        return(n1);
    mu = (isolevel - valp1) / (valp2 - valp1);
    n.x = n1.x + mu * (n2.x - n1.x);
    n.y = n1.y + mu * (n2.y - n1.y);
    n.z = n1.z + mu * (n2.z - n1.z);

    // We may also need to have it so we return p's norm, so maybe return something else instead of a 3D triangle

    return n;
}

void marchingCube::triToObj(string fName) {
    ofstream file;
    file.open(fName + ".obj");
    for (auto& tri : tri_Vector) {
        // Add Vertices to the OBJ file
        for (int i = 0; i < 3; i++) {
            Vector3D vert = tri.coordinates[i];
            string pos = "v " + to_string(vert.x) + " " + to_string(vert.y) + " " + to_string(vert.z);
            file << pos;
        }
        
    }

    for (auto& tri : tri_Vector) {
        // Add Normals to the OBJ file
        for (int i = 0; i < 3; i++) {
            Vector3D normal = tri.normal[i];
            string norm = "vn " + to_string(normal.x) + " " + to_string(normal.y) + " " + to_string(normal.z);
            file << norm;
        }
    }

    //Add the faces
    int idx = 1;
    for (auto& tri : tri_Vector) {
        string face = "f " + to_string(idx) + "//" + to_string(idx) + " " + to_string(idx + 1) + "//" + to_string(idx + 1) + " " + to_string(idx + 2) + "//" + to_string(idx + 2);
        file << face;
        idx += 3;
    }
}



