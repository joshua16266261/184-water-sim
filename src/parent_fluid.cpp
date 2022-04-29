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
    this->fluid = new Fluid(particle_density * length / 10, particle_density * width / 10, particle_density * height / 10, particle_density * length, particle_density * width, particle_density * height);
    this->fp = new FluidParameters(6000, total_time, 60, 5);
    this->diffuse_particles = new vector<DiffuseParticle*>();
    this->dp = dp;
    this->vol_radius = 1.0/particle_density;
	
    this->diffuse_particles->clear();
}

ParentFluid::~ParentFluid() {
    this->diffuse_particles->clear();
}
/*
1. Simulate 1 step of fluid, calulate total acceleration
2. For all particles in fluid:
    Find neighbors
3. For all diffuse particles:
    Calculate dissolution
4. For all diffuse particles:
    Calculate advection
    Collision detection
    Update velocity and position
5. For all particles in fluid:
    Calculate trapped air, wave crest, kinetic energy potential,
    generate diffuse particles by sampling in cyl.,
    set foam lifetime proportional to generation potential
    
 */


void ParentFluid::simulate_step(vector<Vector3D> external_accelerations, vector<CollisionObject *> *collision_objects) {
    //TODO: Off by one error, advection always 1 time step behind.
    //Step 1
    fluid->simulate(fp, external_accelerations, collision_objects);
    Vector3D net_accel = Vector3D();
    for (auto p = begin(external_accelerations); p != end(external_accelerations); p++) {
        net_accel += *p;
    }
    
    //Step 2
    #pragma omp parallel for
    for (auto p = begin(fluid->particles); p != end(fluid->particles); p++) {
        fluid->set_neighbors(&*p, fp->h);
    }
    
    //Step 3
    #pragma omp parallel for
    for (auto p = begin(*diffuse_particles); p != end(*diffuse_particles);) {
        if ((*p)->type == FOAM) {
            bool dissolved = dissolve(*p);
            if (dissolved) {
                p = diffuse_particles->erase(p);
            } else {
                p++;
            }
        }
    }
    
    //Step 4
    #pragma omp parallel for
    for (auto p = begin(*diffuse_particles); p != end(*diffuse_particles);) {
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
    
    //Step 5.1
    #pragma omp parallel for
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
    Particle fake = Particle(p->last_position);
    fluid->set_neighbors(&fake, fp->h);
    if (fake.neighbors->size() == 0) {
        return true;
    }
    Vector3D num = Vector3D();
    double den = 0;
    for (auto q = begin(*(fake.neighbors)); q != end(*(fake.neighbors)); q++) {
		num += ((*q)->position - (*q)->last_position)/dp->delta_t * cubic_spline_kernel((fake.position - (*q)->last_position), dp->h);
        den += cubic_spline_kernel((fake.position - (*q)->last_position), dp->h);
	}
    Vector3D avg_velocity = num / den;
    p->position = p->last_position + dp->delta_t * avg_velocity;
    return false;
}

bool ParentFluid::advect_spray(DiffuseParticle *p, Vector3D external_accelerations) {
    p->velocity = p->last_velocity + dp->delta_t * external_accelerations;
    p->position = p->last_position + dp->delta_t * p->velocity;
    return false;
}

bool ParentFluid::advect_bubbles(DiffuseParticle *p, Vector3D external_accelerations) {
    Particle fake = Particle(p->last_position);
    fluid->set_neighbors(&fake, fp->h);
    if (fake.neighbors->size() == 0) {
        return true;
    }
    Vector3D num = Vector3D();
    double den = 0;
    for (auto q = begin(*(fake.neighbors)); q != end(*(fake.neighbors)); q++) {
		num += ((*q)->position - (*q)->last_position)/dp->delta_t * cubic_spline_kernel((fake.position - (*q)->last_position), dp->h);
        den += cubic_spline_kernel((fake.position - (*q)->last_position), dp->h);
	}
    Vector3D avg_velocity = num / den;
    p->velocity = p->last_velocity + dp->delta_t * (-dp->k_b * Vector3D(0, 0, -9.81) + dp->k_d * (avg_velocity - p->last_velocity) / dp->delta_t);
    p->position = p->last_position + dp->delta_t * p->velocity;
    return false;
}

bool ParentFluid::dissolve(DiffuseParticle *p) {
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
    //Assuming mass is 1
    double kin = 0.5 * p->velocity.norm2();
    return clamp2(kin, dp->t_k_min, dp->t_k_max);
}

double ParentFluid::ta_potential(Particle *p) {
    double v_diff = 0;
    for (auto q = begin(*p->neighbors); q != end(*p->neighbors); q++) {
        Vector3D vij = p->velocity - (*q)->velocity;
        Vector3D vij_hat = vij / vij.norm();
        Vector3D xij = p->position - (*q)->position;
        Vector3D xij_hat = xij / xij.norm();
        v_diff += vij.norm() * (1 - dot(vij_hat, xij_hat)) * symm_kernel(xij, dp->h);
    }
    return clamp2(v_diff, dp->t_ta_min, dp->t_ta_max);
}

double ParentFluid::wc_potential(Particle *p) {
    //TODO: Do this
    return 0;
}

void ParentFluid::generate(Particle *p, int n) {
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
    for (int _ = 0; _ < n; _++) {
        double xr = double(rand()) / RAND_MAX;
        double r = vol_radius * sqrt(xr);
        double xtheta = double(rand()) / RAND_MAX;
        double theta = xtheta * 2 * PI;
        double xh = double(rand()) / RAND_MAX;
        double h = xh * (dp->delta_t * p->velocity).norm();
        Vector3D xd = p->last_position + r * cos(theta) * e1 + r * sin(theta) * e2 + h * p->velocity / p->velocity.norm();
        Vector3D vd = r * cos(theta) * e1 + r * sin(theta) * e2 + p->velocity;
        DiffuseParticle *d = new DiffuseParticle(xd, vd, BUBBLE);
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

