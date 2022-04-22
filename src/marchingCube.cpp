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

// grid starts from (0,0,0) and ends at box_dimensions
// the list of particles
// the hashmap of string to list particles within the same box
// cube dimensions we want: ex: box_dim = (9,9,9), cube_dim = (3,3,3) so (27 in total since 9/3 ^ 3)
// search_radius, mass, h,  
void init(Vector3D box_dimensions, Vector3D singular_cube_dimensions,
&vector<Particle> particles, &unordered_map<string, &vector<Particle *> *> hash_to_particles,
 float h, float search_radius, float particle_mass, 
 float density) {

// init variables to use in this files functions
m_box_dimensions = box_dimensions;
m_singular_cube_dimensions = singular_cube_dimensions;
m_particles = particles;
m_hash_to_particles = hash_to_particles;
m_h = h;
m_search_radius = search_radius;
m_particle_mass = particle_mass;

// slice up the main box into cubes
// one way is to 
Vector3D 

}

// r is the vector from particle to position p
float isotropic_kernel(Vector3D r, float h) {
    float constant = 315 / (64 * PI * pow(h, 9));
    if (r.norm() >= 0 and r.norm() <= h) {
        return constant * pow(h * h - r.norm() * r.norm());
    } else {
        return 0;
    }
}

// Compute the density of a pos (will be used as the isovalue)
float isovalue(Vector3D pos, float h){
    float rho = 0;
    // Hash the position vector
    // Find particles that are neighboring
    // Check for particles that are in specific radius
    String vortex_key = hash_position(pos, h);
    for (auto q = begin(*(hash_to_particles[vortex_key])); q != end(*(hash_to_particles[vortex_key])); q++) {
        if ((pos - (*q)->position).norm() <= search_radius) {
            rho += mass * isotropic_kernel(pos - (q*)->pos, h);
        }
    }
    return rho;
}

// Create the Marching Cube Grid (each sub-cube within the whole space)
// Inputs: the index of the cube in all directions (i.e. 1 in X, 4 in Y, and 2 in Z)
void createCube(Cube &cube, Vector3D index) {
    // Assign all the vertices within a cube to its approperiate (x, y, z) positions
    // Assign all approperiate isovalues to each vertex
    for (int i = 0; i < 8; i++) {
        cube.vertices[i] = ...;
        cube.isovalues[i] = isovalue(cube.vertices[i]);
    }
    
}


