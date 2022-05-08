#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "CGL/vector3D.h"
#include "collision/collisionObject.h"
#include "particle.h"
#include "fluid.h"
#include "diffuse.h"
#include "parent_fluid.h"

#include <stdio.h>

using namespace std;

ParentFluid::ParentFluid(double length, double width, double height, int particle_density, int fps, double total_time, DiffuseParameters *dp) {
    this->length = length;
    this->width = width;
    this->height = height;
    this->particle_density = particle_density;
    this->fps = fps;
    this->total_time = total_time;
    this->diffuse_particles = new vector<DiffuseParticle*>();
    this->dp = dp;
    this->vol_radius = 1.0 / particle_density;
	
    this->diffuse_particles->clear();
}

ParentFluid::~ParentFluid() {
    this->diffuse_particles->clear();
}

void ParentFluid::simulate_step(vector<Vector3D> external_accelerations, vector<CollisionObject *> *collision_objects) {
    Vector3D net_accel = Vector3D();
    for (auto p = begin(external_accelerations); p != end(external_accelerations); p++) {
        net_accel += *p;
    }
	
	// Build spatial map with width 2h for advection
	fluid->build_spatial_map(2 * dp->h);
    
    // Dissolution
    #pragma omp parallel for
    for (auto p = begin(*diffuse_particles); p != end(*diffuse_particles); ) {
        if ((*p)->type == FOAM) {
            bool dissolved = dissolve(*p);
            if (dissolved) {
                p = diffuse_particles->erase(p);
            } else {
				p++;
			}
        } else {
			p++;
		}
    }
    
    // Advection
    #pragma omp parallel for
    for (auto p = begin(*diffuse_particles); p != end(*diffuse_particles); ) {
        bool dead = advect(*p, net_accel);
        if (dead) {
            p = diffuse_particles->erase(p);
        } else {
            collide(*p, collision_objects);
            (*p)->last_position = (*p)->position;
            (*p)->last_velocity = (*p)->velocity;
            p++;
        }
    }
	
	// Rebuild spatial map with width h for generation
	fluid->build_spatial_map(dp->h);
	
	#pragma omp parallel for
	for (auto p = begin(fluid->particles); p != end(fluid->particles); p++) {
		fluid->set_neighbors(&*p, dp->h);
	}
    
	// Set rho for all fluid particles
	#pragma omp parallel for
	for (auto p = begin(fluid->particles); p != end(fluid->particles); p++) {
		double rho_i = 0;
		for (auto q = begin(*(p->neighbors)); q != end(*(p->neighbors)); q++) {
			rho_i += fluid->poly6_kernel(p->position - (*q)->position, dp->h);
		}
		p->rho = rho_i;
	}
	
	// Calculate normals for all fluid particles (0 if not surface particle)
	#pragma omp parallel for
	for (auto p = begin(fluid->particles); p != end(fluid->particles); p++) {
		set_normal(&*p);
	}
	
	// Calculate potentials and generate diffuse particles
    for (auto p = begin(fluid->particles); p != end(fluid->particles); p++) {
        double I_k = k_potential(&*p);
        double I_ta = ta_potential(&*p);
        double I_wc = wc_potential(&*p);
        int n_d = int(I_k * (dp->k_ta * I_ta + dp->k_wc * I_wc) * dp->delta_t);
        generate(&*p, n_d);
    }
}

bool ParentFluid::advect(DiffuseParticle *p, Vector3D external_accelerations) {
    if (p->type == SPRAY) {
        return advect_spray(p, external_accelerations);
    } else if (p->type == FOAM) {
        return advect_foam(p, external_accelerations);
    } else {
        return advect_bubbles(p, external_accelerations);
    }
}

bool ParentFluid::advect_foam(DiffuseParticle *p, Vector3D external_accelerations) {
	// Advect p as a foam particle
	// Return true if p has no neighbors and should be deleted
    Particle fake = Particle(p->last_position);
    fluid->set_neighbors(&fake, 2 * dp->h);
    if (fake.neighbors->size() == 0) {
        return true;
    }
	
    Vector3D num = Vector3D();
    double den = 0;
    for (auto q = begin(*(fake.neighbors)); q != end(*(fake.neighbors)); q++) {
		num += ((*q)->position - (*q)->last_position)/dp->delta_t * cubic_spline_kernel(fake.position - (*q)->last_position, dp->h);
        den += cubic_spline_kernel(fake.position - (*q)->last_position, dp->h);
	}
    Vector3D avg_velocity = num / den;
    p->position = p->last_position + dp->delta_t * avg_velocity;
    return false;
}

bool ParentFluid::advect_spray(DiffuseParticle *p, Vector3D external_accelerations) {
	// Advect p as a spray particle
	// Spray particles do not require neighbors so we never delete them
    p->velocity = p->last_velocity + dp->delta_t * external_accelerations;
    p->position = p->last_position + dp->delta_t * p->velocity;
    return false;
}

bool ParentFluid::advect_bubbles(DiffuseParticle *p, Vector3D external_accelerations) {
	// Advect p as a bubble particle
	// Return true if p has no neighbors and should be deleted
    Particle fake = Particle(p->last_position);
    fluid->set_neighbors(&fake, 2 * dp->h);
    if (fake.neighbors->size() == 0) {
        return true;
    }
	
    Vector3D num = Vector3D();
    double den = 0;
    for (auto q = begin(*(fake.neighbors)); q != end(*(fake.neighbors)); q++) {
		num += ((*q)->position - (*q)->last_position)/dp->delta_t * cubic_spline_kernel(fake.position - (*q)->last_position, dp->h);
        den += cubic_spline_kernel(fake.position - (*q)->last_position, dp->h);
	}
    Vector3D avg_velocity = num / den;
    p->velocity = p->last_velocity + dp->delta_t * (-dp->k_b * Vector3D(0, 0, -9.81) + dp->k_d * (avg_velocity - p->last_velocity) / dp->delta_t);
    p->position = p->last_position + dp->delta_t * p->velocity;
    return false;
}

bool ParentFluid::dissolve(DiffuseParticle *p) {
	// Return true if p has exceeded liftime
	// Otherwise, decrease p->ttl and return false
    if (p->ttl <= 0) {
        return true;
    }
    p->ttl -= dp->delta_t;
    return false;
}

void ParentFluid::collide(DiffuseParticle *p, vector<CollisionObject *> *collision_objects) {
    for (auto c = begin(*collision_objects); c != end(*collision_objects); c++) {
        (*c)->collide(*p);
        p->velocity = (p->position - p->last_position) / dp->delta_t;
    }
}

double ParentFluid::k_potential(Particle *p) {
    // Assuming mass is 1
    double kin = 0.5 * p->velocity.norm2();
    return clamp2(kin, dp->t_k_min, dp->t_k_max);
}

double ParentFluid::ta_potential(Particle *p) {
	// Calculate trapped air potential
    double v_diff = 0;
    for (auto q = begin(*p->neighbors); q != end(*p->neighbors); q++) {
        Vector3D vij = p->velocity - (*q)->velocity;
		Vector3D xij = p->position - (*q)->position;
		if (vij.norm() > EPS_F && xij.norm() > EPS_F) {
			Vector3D vij_hat = vij / vij.norm();
			Vector3D xij_hat = xij / xij.norm();
			v_diff += vij.norm() * (1 - dot(vij_hat, xij_hat)) * symm_kernel(xij, dp->h);
		}
    }
    return clamp2(v_diff, dp->t_ta_min, dp->t_ta_max);
}

void ParentFluid::set_normal(Particle *p) {
	// Calculate normal vector of p
	// See equations 15, 16 of https://matthias-research.github.io/pages/publications/sca03.pdf
	Vector3D n = Vector3D();
	for (auto q = begin(*p->neighbors); q != end(*p->neighbors); q++) {
		n += 1.0 / (*q)->rho * fluid->grad_spiky_kernel(p->position - (*q)->position, fp->h);
	}
	
	if (n.norm() > dp->l) {
		// p is a surface particle
		n.normalize();
		p->normal = -n;
	} else {
		p->normal = Vector3D();
	}
}

double ParentFluid::wc_potential(Particle *p) {
	// Calculate wave crest potential
	if (p->normal.norm() == 0) {
		return 0;
	}
	
	Vector3D n_i = p->normal;
	Vector3D x_i = p->position;
	double kappa_tilde = 0;
	for (auto q = begin(*p->neighbors); q != end(*p->neighbors); q++) {
		if ((*q)->normal.norm() > 0) {
			Vector3D x_ji = (*q)->position - x_i;
			x_ji.normalize();
			if (dot(x_ji, n_i) < 0) {
				kappa_tilde += (1 - dot(n_i, (*q)->normal)) * symm_kernel(x_i - (*q)->position, dp->h);
			}
		}
	}
	
	if (dot(p->velocity / p->velocity.norm(), p->normal) >= 0.6) {
		return clamp2(kappa_tilde, dp->t_wc_min, dp->t_wc_max);
	}
	
    return 0;
}

void ParentFluid::generate(Particle *p, int n) {
	// Get e_1 and e_2 as in Fig. 3 of Ihmsen et al. by performing Gram-Schmidt
    Vector3D e1;
    if (p->velocity.x != 0) {
        e1 = Vector3D(1,0,0) - p->velocity.x / p->velocity.norm2() * p->velocity;
    } else if (p->velocity.y != 0) {
        e1 = Vector3D(0,1,0) - p->velocity.y / p->velocity.norm2() * p->velocity;
    } else if (p->velocity.z != 0) {
        e1 = Vector3D(0,0,1) - p->velocity.z / p->velocity.norm2() * p->velocity;
    }
    e1.normalize();
    Vector3D e2 = cross(p->velocity, e1);
    e2.normalize();
	
	// Generate diffuse particles
    for (int _ = 0; _ < n; _++) {
        double xr = double(rand()) / RAND_MAX;
        double r = vol_radius * sqrt(xr);
        double xtheta = double(rand()) / RAND_MAX;
        double theta = xtheta * 2 * PI;
        double xh = double(rand()) / RAND_MAX;
        double h = xh * (dp->delta_t * p->velocity).norm();
		
		Vector3D xd = p->position + r * cos(theta) * e1 + r * sin(theta) * e2 + h * p->velocity / p->velocity.norm();
        Vector3D vd = r * cos(theta) * e1 + r * sin(theta) * e2 + p->velocity;
        DiffuseParticle *d = new DiffuseParticle(xd, vd, BUBBLE);
		
		// Classify d based on number of fluid particle neighbors
        int neighbors = p->neighbors->size();
        if (neighbors < 6) {
            d->type = SPRAY;
        } else if (neighbors <= 20) {
            d->type = FOAM;
            d->ttl = dp->lifetime * n;
        }
        diffuse_particles->emplace_back(d);
    }
}

// Return the total number of particles in the list.
// Always required!
size_t ParentFluid::size() const {
	return diffuse_particles->size();
}

// Get the world-space position of the nth particle.
// Required by rasterizeSpheres().
void ParentFluid::getPos(size_t n, Vec3R& xyz) const {
	Vector3D pos = diffuse_particles->at(n)->position;
	xyz.x() = pos.x;
	xyz.y() = pos.y;
	xyz.z() = pos.z;
}

// Get the world-space position and radius of the nth particle.
// Required by rasterizeSpheres().
void ParentFluid::getPosRad(size_t n, Vec3R& xyz, Real& radius) const {
	getPos(n, xyz);
	//TODO: Possibly change radius
	radius = vol_radius;
}

// Get the world-space position, radius and velocity of the nth particle.
// Required by rasterizeTrails().
void ParentFluid::getPosRadVel(size_t n, Vec3R& xyz, Real& radius, Vec3R& velocity) const {
	getPosRad(n, xyz, radius);
	Vector3D v = diffuse_particles->at(n)->velocity;
	velocity.x() = v.x;
	velocity.y() = v.y;
	velocity.z() = v.z;
}
