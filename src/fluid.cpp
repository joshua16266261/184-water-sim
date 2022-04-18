#include <iostream>
#include <math.h>
#include <random>
#include <vector>
#include <string>

#include "fluid.h"
#include "particle.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Fluid::Fluid(float length, float width, float height, int num_length_particles, int num_width_particles, int num_height_particles) {
  this->length = length;
  this->width = width;
  this->height = height;
  this->num_length_particles = num_length_particles;
  this->num_width_particles = num_width_particles;
  this->num_height_particles = num_height_particles;

  buildGrid();
}

Fluid::~Fluid() {
  particles.clear();
}

void Fluid::buildGrid() {
    // Fill the rectangular prism with particles
    // Depth major order -> row major order
	float x, y, z;
    for (int depth = 0; depth < num_height_particles; depth++) {
        z = depth * height / (num_height_particles - 1);
        for (int row = 0; row < num_width_particles; row++) {
            y = row * width / (num_width_particles - 1);
		    for (int col = 0; col < num_length_particles; col++) {
                x = col * length / (num_length_particles - 1);
                particles.emplace_back(Particle(Vector3D(x, y, z)));
            }
		}
	}
}

void Fluid::simulate(FluidParameters *fp,
                vector<Vector3D> external_accelerations,
                vector<CollisionObject *> *collision_objects) {
    double mass = length * width * height * fp->density / num_length_particles / num_width_particles / num_height_particles;
    double delta_t = 1.0f / fp->fps;

    // Compute total force acting on each point mass.
	Vector3D total_external_force = Vector3D();
	for (Vector3D a : external_accelerations) {
		total_external_force += mass * a;
	}
	
    // Apply external forces
	for (auto p = begin(particles); p != end(particles); p++) {
		Vector3D v = p->velocity(delta_t);
        v += delta_t * total_external_force;
        Vector3D pos = p->position;
        p->position = pos + delta_t * v;
        p->last_position = pos;
	}

    build_spatial_map(fp->h);

    // Find neighboring particles
    for (auto p = begin(particles); p != end(particles); p++) {
        string key = hash_position(p->position, fp->h);
        (p->neighbors)->clear();
		
        for (int i = -1; i < 2; i++) {
            for (int j = -1; j < 2; j++) {
                for (int k = -1; k < 2; k++) {
					int index1 = key.find(":");
                    int x = stoi(key.substr(0, index1), nullptr, 10);
                    int index2 = key.substr(index1+1, key.size()).find(":");
                    int y = stoi(key.substr(index1+1, index2), nullptr, 10);
                    int z = stoi(key.substr(index2+1, key.size()), nullptr, 10);
                    string neighbor_key = to_string(x + i);
                    neighbor_key.append(":");
                    neighbor_key.append(to_string(y + j));
                    neighbor_key.append(":");
                    neighbor_key.append(to_string(z + k));
                    if (map.count(neighbor_key) > 0) {
                        for (auto q = begin(*(map[neighbor_key])); q != end(*(map[neighbor_key])); q++) {
							if ((p->position - (*q)->position).norm() <= fp->h) {
								(p->neighbors)->emplace_back(*q);
							}
                        }
                    }
                }
            }
        }
    }

    for (int _ = 0; _ < fp->solverIters; _++) {
		for (auto p = begin(particles); p != end(particles); p++) {
			// Calculate lambda_i (line 13 of Algorithm 1)
			float rho_i = 0;
			for (auto q = begin(*(p->neighbors)); q != end(*(p->neighbors)); q++) {
				rho_i += pow(pow(fp->h, 2) - (p->position - (*q)->position).norm2(), 3);
				// TODO
			}
			rho_i *= mass * 315.0 / (64 * PI * pow(fp->h, 9));
		}
        
    }
	
	

//   // TODO (Part 2): Use Verlet integration to compute new point mass positions
// 	for (auto pm = begin(point_masses); pm != end(point_masses); pm++) {
// 		if (!pm->pinned) {
// 			Vector3D xt = pm->position;
// 			float d = cp->damping / 100;
// 			Vector3D x_prev = pm->last_position;
// 			Vector3D at = pm->forces / mass;
// 			pm->last_position = pm->position;
// 			pm->position = xt + (1 - d) * (xt - x_prev) + at * delta_t * delta_t;
// 		}
// 	}

//   // TODO (Part 4): Handle self-collisions.
// 	build_spatial_map();
// 	for (auto pm = begin(point_masses); pm != end(point_masses); pm++) {
// 		self_collide(*pm, simulation_steps);
// 	}


//   // TODO (Part 3): Handle collisions with other primitives.
// 	for (auto c = begin(*collision_objects); c != end(*collision_objects); c++) {
// 		for (auto pm = begin(point_masses); pm != end(point_masses); pm++) {
// 			(*c)->collide(*pm);
// 		}
// 	}


//   // TODO (Part 2): Constrain the changes to be such that the spring does not change
//   // in length more than 10% per timestep [Provot 1995].
// 	for (auto s = begin(springs); s != end(springs); s++) {
// 		Vector3D a_to_b = s->pm_b->position - s->pm_a->position;
// 		float diff = a_to_b.norm() - 1.1 * s->rest_length;
// 		if (diff > 0) {
// 			a_to_b.normalize();
// 			if (!s->pm_a->pinned && !s->pm_b->pinned) {
// 				a_to_b *= diff / 2;
// 				s->pm_a->position += a_to_b;
// 				s->pm_b->position -= a_to_b;
// 			} else if (!s->pm_a->pinned) {
// 				a_to_b *= diff;
// 				s->pm_a->position += a_to_b;
// 			} else if (!s->pm_b->pinned) {
// 				a_to_b *= diff;
// 				s->pm_b->position -= a_to_b;
// 			}
// 		}

// 	}

}

void Fluid::build_spatial_map(float h) {
    for (const auto &entry : map) {
        delete(entry.second);
    }
    map.clear();

    // Build a spatial map out of all of the particles
	for (auto p = begin(particles); p != end(particles); p++) {
		string key = hash_position(p->position, h);
		if (map.count(key) > 0) {
			map[key]->emplace_back(&*p);
		} else {
			vector<Particle*> *v = new vector<Particle*>;
			map[key] = v;
			v->emplace_back(&*p);
		}
	}
}

string Fluid::hash_position(Vector3D pos, float h) {
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

