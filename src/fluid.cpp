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

Fluid::Fluid(double length, double width, double height, int num_length_particles, int num_width_particles, int num_height_particles) {
  this->length = length;
  this->width = width;
  this->height = height;
  this->num_length_particles = num_length_particles;
  this->num_width_particles = num_width_particles;
  this->num_height_particles = num_height_particles;
	
	particles.clear();
	buildGrid();
}

Fluid::~Fluid() {
  particles.clear();
}

void Fluid::buildGrid() {
    // Fill the rectangular prism with particles
    // Depth major order -> row major order
	double random_offset_bound = 0.001;
	
	particles.clear();
	
	double x, y, z;
    for (int depth = 0; depth < num_height_particles; depth++) {
        z = depth * height / (num_height_particles - 1);
        for (int row = 0; row < num_width_particles; row++) {
            y = row * width / (num_width_particles - 1);
		    for (int col = 0; col < num_length_particles; col++) {
                x = col * length / (num_length_particles - 1);
                particles.emplace_back(Particle(Vector3D(x + double(rand()) / RAND_MAX * 2 * random_offset_bound - random_offset_bound,
														 y + double(rand()) / RAND_MAX * 2 * random_offset_bound - random_offset_bound,
														 z + double(rand()) / RAND_MAX * 2 * random_offset_bound - random_offset_bound)));
            }
		}
	}
}

void Fluid::set_neighbors(Particle *p, double h) {
	// Get all particles distance h or less away and put them in p->neighbors
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
						double dist = (p->position - (*q)->position).norm();
						if (dist <= h*2 && dist > 0) {
							(p->neighbors)->emplace_back(*q);
						}
					}
				}
			}
		}
	}
}

void Fluid::calculate_lambda(Particle *p, double mass, double density, double h, double relaxation) {
	// Calculate p->lambda as given in Equation 9
	if (p->neighbors->size() == 0) {
		p->lambda = 0;
		return;
	}
	// Calculate C_i
	double rho_i = 0;
	for (auto q = begin(*(p->neighbors)); q != end(*(p->neighbors)); q++) {
		rho_i += poly6_kernel(p->position - (*q)->position, h);
	}
	rho_i *= mass;
	
	double C_i = rho_i / density - 1;
	
	Vector3D pi_term = Vector3D();
	double not_pi_term = 0;
	for (auto pj = begin(*(p->neighbors)); pj != end(*(p->neighbors)); pj++) {
		Vector3D r = p->position - (*pj)->position;
		Vector3D grad_j = grad_spiky_kernel(r, h);

		pi_term += grad_j;
		not_pi_term += grad_j.norm2();
	}
	double denom = 1.0 / pow(density, 2) * (pi_term.norm2() + not_pi_term);
	
	p->lambda = -C_i / (denom + relaxation);
}

void Fluid::calculate_delta_p(Particle *p, double h, Vector3D delta_q, double k, double n, double density) {
	// Calculate p->delta_p as given in Equation 14
	Vector3D delta_p = Vector3D();
	for (auto pj= begin(*(p->neighbors)); pj != end(*(p->neighbors)); pj++) {
		Vector3D r = p->position - (*pj)->position;
		Vector3D grad_j = grad_spiky_kernel(r, h);
		
		double numer = poly6_kernel(r, h);
		double denom = poly6_kernel(delta_q, h);
		double s_corr = -k * pow(numer / denom, n);
		
		delta_p += (p->lambda + (*pj)->lambda + s_corr) * grad_j;
	}
	delta_p *= 1.0 / density;
	p->delta_p = delta_p;
}

void Fluid::collision_detection(Particle *p, vector<CollisionObject *> *collision_objects, double cr, double delta_t) {
	// Set p->temp_velocity
	// If p collides with any objects, modify p->delta_p and p->temp_velocity
	for (auto o = begin(*collision_objects); o != end(*collision_objects); o++) {
		p->temp_velocity = (p->position + p->delta_p - p->last_position) / delta_t;
		(*o)->collide(*p, cr, delta_t);
	}
}

void Fluid::calculate_omega(Particle *p, double h) {
	// Calculate p->omega as given in Equation 15
	p->omega = Vector3D();

	for (auto pj = begin(*p->neighbors); pj != end(*p->neighbors); pj++) {
		Vector3D v_ij = (*pj)->velocity - p->velocity;
		//TODO: Remove negative sign of grad
		p->omega += cross(v_ij, -grad_spiky_kernel(p->position - (*pj)->position, h));
	}
}

void Fluid::vorticity(Particle *p, double h, double delta_t, double vorticity_eps, double mass) {
	// Calculate vorticity force on p as given in Equation 16 and apply it to p->temp_velocity
	Vector3D N = Vector3D();
	
	if (p->omega.norm() > EPS_F) {
		for (auto pj = begin(*p->neighbors); pj != end(*p->neighbors); pj++) {
			Vector3D r = (*pj)->position - p->position;
            if (r.norm() <= h) {
                double d_omega = (*pj)->omega.norm() - p->omega.norm();
                N += Vector3D(d_omega / r.x, d_omega / r.y, d_omega / r.z);
            }
		}

		if (N.norm() > EPS_F) {
			N.normalize();
		}
	}

	p->temp_velocity += delta_t * vorticity_eps * cross(N, p->omega) / mass;
}

void Fluid::viscosity(Particle *p, double c, double h) {
	// Update p->temp_velocity with viscosity as given in Equation 17
	for (auto pj = begin(*p->neighbors); pj != end(*p->neighbors); pj++) {
		Vector3D v_ij = (*pj)->velocity - p->velocity;
		p->temp_velocity += c * v_ij * poly6_kernel(p->position - (*pj)->position, h);
	}
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
			calculate_lambda(&*p, mass, fp->density, fp->h, fp->relaxation);
		}
		
		// Line 12 of Algorithm 1
		#pragma omp parallel for
		for (auto p = begin(particles); p != end(particles); p++) {
			// Calculate delta pi (line 13 of Algorithm 1)
			calculate_delta_p(&*p, fp->h, fp->delta_q, fp->k, fp->n, fp->density);
			
			// Collision detection and response (line 14 of Algorithm 1)
			collision_detection(&*p, collision_objects, fp->cr, delta_t);
		}
		
		// Line 16 of Algorithm 1
		#pragma omp parallel for
		for (auto p = begin(particles); p != end(particles); p++) {
			// Update position (line 17 of Algorithm 1)
			p->position += p->delta_p;
		}
    }
	
	#pragma omp parallel for
	for (auto p = begin(particles); p != end(particles); p++) {
		viscosity(&*p, fp->c, fp->h);
	}

	#pragma omp parallel for
	for (auto p = begin(particles); p != end(particles); p++) {
		calculate_omega(&*p, fp->h);
	}

	#pragma omp parallel for
	for (auto p = begin(particles); p != end(particles); p++) {
		vorticity(&*p, fp->h, delta_t, fp->vorticity_eps, mass);
	}
	
	#pragma omp parallel for
	for (auto p = begin(particles); p != end(particles); p++) {
		p->velocity = p->temp_velocity;
	}
}

void Fluid::build_spatial_map(double h) {
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

string Fluid::hash_position(Vector3D pos, double h) {
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

