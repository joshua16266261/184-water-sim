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
	float random_offset_bound = 0.001;
	
	particles.clear();
	
	float x, y, z;
    for (int depth = 0; depth < num_height_particles; depth++) {
        z = depth * height / (num_height_particles - 1);
        for (int row = 0; row < num_width_particles; row++) {
            y = row * width / (num_width_particles - 1);
		    for (int col = 0; col < num_length_particles; col++) {
                x = col * length / (num_length_particles - 1);
                particles.emplace_back(Particle(Vector3D(x + float(rand()) / RAND_MAX * 2 * random_offset_bound - random_offset_bound,
														 y + float(rand()) / RAND_MAX * 2 * random_offset_bound - random_offset_bound,
														 z + float(rand()) / RAND_MAX * 2 * random_offset_bound - random_offset_bound)));
            }
		}
	}
}

void Fluid::set_neighbors(Particle *p, float h) {
	string key = hash_position(p->position, h);
	(p->neighbors)->clear();
	
	for (int i = -1; i < 2; i++) {
		for (int j = -1; j < 2; j++) {
			for (int k = -1; k < 2; k++) {
				int index1 = key.find(":");
				int x = stoi(key.substr(0, index1), nullptr, 10);
				int index2 = key.substr(index1 + 1, key.size()).find(":");
				int y = stoi(key.substr(index1 + 1, index2), nullptr, 10);
				int z = stoi(key.substr(index1 + index2 + 2, key.size()), nullptr, 10);
				
				string neighbor_key = to_string(x + i);
				neighbor_key.append(":");
				neighbor_key.append(to_string(y + j));
				neighbor_key.append(":");
				neighbor_key.append(to_string(z + k));
				
				if (map.count(neighbor_key) > 0 && map[neighbor_key]->size() > 0) {
					for (auto q = begin(*(map[neighbor_key])); q != end(*(map[neighbor_key])); q++) {
						float dist = (p->position - (*q)->position).norm();
						if (dist <= h && dist > 0) {
							(p->neighbors)->emplace_back(*q);
						}
					}
				}
			}
		}
	}
}

float Fluid::get_avg_spacing() {
	//TODO: Set fp->h based on 3 times average particle spacing
//	build_spatial_map(<#float h#>)
	
	float spacing = 0;
	for (auto p = begin(particles); p != end(particles); p++) {
		
	}
	
	return 0;
}

void Fluid::simulate(FluidParameters *fp,
                vector<Vector3D> external_accelerations,
                vector<CollisionObject *> *collision_objects) {
    double mass = length * width * height * fp->density / num_length_particles / num_width_particles / num_height_particles;
    double delta_t = 1.0 / fp->fps;
	
    // Compute total force acting on each point mass.
	Vector3D total_external_force = Vector3D();
	for (Vector3D a : external_accelerations) {
		total_external_force += mass * a;
	}
	
    // Apply external forces
	#pragma omp parallel for
	for (auto p = begin(particles); p != end(particles); p++) {
		Vector3D v = p->velocity;
        v += delta_t * total_external_force / mass;
        Vector3D pos = p->position;
        p->position = pos + delta_t * v;
        p->last_position = pos;
	}

    build_spatial_map(fp->h);

    // Find neighboring particles
	#pragma omp parallel for
    for (auto p = begin(particles); p != end(particles); p++) {
		set_neighbors(&*p, fp->h);
    }
	
	// Line 8 of Algorithm 1
    for (int _ = 0; _ < fp->solverIters; _++) {
		// Line 9 of Algorithm 1
		#pragma omp parallel for
		for (auto p = begin(particles); p != end(particles); p++) {
			// Calculate lambda_i (line 10 of Algorithm 1)
			
			// Calculate C_i
			float rho_i = 0;
			for (auto q = begin(*(p->neighbors)); q != end(*(p->neighbors)); q++) {
				rho_i += poly6_kernel(p->position - (*q)->position, fp->h);
			}
			rho_i *= mass;
			
			float C_i = rho_i / fp->density - 1;
			Vector3D pi_term;
			float not_pi_term;
			for (auto pj = begin(*(p->neighbors)); pj != end(*(p->neighbors)); pj++) {
				Vector3D r = p->position - (*pj)->position;
				Vector3D grad_j = grad_spiky_kernel(r, fp->h);
				
				pi_term += grad_j;
				not_pi_term += grad_j.norm2();
			}
			float denom = 1.0 / fp->density * (pi_term.norm2() + not_pi_term);
			
			p->lambda = -C_i / (denom + fp->relaxation);
		}
		
		// Line 12 of Algorithm 1
		#pragma omp parallel for
		for (auto p = begin(particles); p != end(particles); p++) {
			// Calculate delta pi (line 13 of Algorithm 1)
			Vector3D delta_p = Vector3D();
			for (auto pj= begin(*(p->neighbors)); pj != end(*(p->neighbors)); pj++) {
				Vector3D r = p->position - (*pj)->position;
				Vector3D grad_j = grad_spiky_kernel(r, fp->h);
				
				float numer = M4_kernel(r, fp->h);
				float denom = M4_kernel(fp->delta_q, fp->h);
				
				float s_corr = -fp->k * pow(numer / denom, fp->n);
				
				delta_p += (p->lambda + (*pj)->lambda + s_corr) * grad_j;
			}
			delta_p *= 1.0 / fp->density;
			
			// Collision detection and response (line 14 of Algorithm 1)
			p->delta_p = delta_p;
			for (auto o = begin(*collision_objects); o != end(*collision_objects); o++) {
				p->temp_velocity = (p->position + delta_p - p->last_position) / delta_t;
				(*o)->collide(*p, fp->cr, delta_t);
			}
		}
		
		// Line 16 of Algorithm 1
		#pragma omp parallel for
		for (auto p = begin(particles); p != end(particles); p++) {
			// Update position (line 17 of Algorithm 1)
			p->position += p->delta_p;
		}
		
    }
	
	// Vorticity (line 22 of Algorithm 1)
	#pragma omp parallel for
	for (auto p = begin(particles); p != end(particles); p++) {
		p->omega = Vector3D();
		
		for (auto pj = begin(*p->neighbors); pj != end(*p->neighbors); pj++) {
			Vector3D v_ij = (*pj)->velocity - p->velocity;
			p->omega += cross(v_ij, -grad_spiky_kernel(p->position - (*pj)->position, fp->h));
			
			// Viscosity
			p->temp_velocity += fp->c * v_ij * M4_kernel((*pj)->position - p->position, fp->h);
		}
	}
	
	#pragma omp parallel for
	for (auto p = begin(particles); p != end(particles); p++) {
		Vector3D N = Vector3D();
		
		for (auto pj = begin(*p->neighbors); pj != end(*p->neighbors); pj++) {
			float d_omega = ((*pj)->omega - p->omega).norm();
			Vector3D r= (*pj)->position - p->position;
			N += Vector3D(d_omega / r.x, d_omega / r.y, d_omega / r.z);
		}
		
		if (N.norm() > EPS_F) {
			N.normalize();
		}
		
		p->temp_velocity += delta_t * fp->vorticity_eps * cross(N, p->omega);
	}
	
	#pragma omp parallel for
	for (auto p = begin(particles); p != end(particles); p++) {
		p->velocity = p->temp_velocity;
	}
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

