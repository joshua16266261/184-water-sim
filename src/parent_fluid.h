#ifndef PARENT_FLUID_H
#define PARENT_FLUID_H

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

using namespace CGL;
using namespace std;


struct DiffuseParameters {
    DiffuseParameters() {}
    /*
    t_k_min, t_k_max: Parameters used in the clamping function for potential to generate diffuse particles due to kinetic energy. (see: page 4 of Ihmsen et al.)
    t_ta_min, t_ta_max: Same as above for trapped air
    t_wc_min, t_wc_max: Same as above for wave crests
    k_b: bouyancy of bubbles (see: page 5 of Imsen et al.)
    k_d: drag of bubbles. If k_d is 1, air bubbles are immediately dragged into the flow dir. of the fluid (see: page 5 of Imsen et al.)
    k_wc: max number of particles for wave crests per second
    k_ta: max number of particles for trapped air per second
    For default values, the paper runs their sim at 50 fps
     */

    DiffuseParameters(double h, double k_b, double k_d, double k_wc, double k_ta, double t_wc_min = 2, double t_wc_max = 8, double t_ta_min = 5, double t_ta_max = 20, double t_k_min = 5, double t_k_max = 50, double delta_t = 1.0/60, double lifetime = 1.0)
        : h(h), delta_t(delta_t), k_b(k_b), k_d(k_d), k_wc(k_wc), k_ta(k_ta), t_wc_min(t_wc_min), t_wc_max(t_wc_max), t_ta_min(t_ta_min), t_ta_max(t_ta_min), t_k_min(t_k_min), t_k_max(t_k_max), lifetime(lifetime) {}
    ~DiffuseParameters() {}

    double h;
    double delta_t;
    double k_b;
    double k_d;
    double k_wc;
    double k_ta;
    double t_wc_min;
    double t_wc_max;
    double t_ta_min;
    double t_ta_max;
    double t_k_min;
    double t_k_max;
    double lifetime;

};


struct ParentFluid {
  ParentFluid() {}
  ParentFluid(double length, double width, double height, int particle_density, int fps, double total_time, DiffuseParameters *dp);
  ~ParentFluid();

//Sim parameters
    int fps;
    double total_time;

//Fluid Parameters
    double length;
    double width;
    double height;
    int particle_density;
    Fluid *fluid;
    FluidParameters *fp;
    vector<DiffuseParticle*> *diffuse_particles;
    DiffuseParameters *dp;
    double vol_radius;


    void simulate_step(vector<Vector3D> external_accelerations,
                vector<CollisionObject *> *collision_objects);

    bool advect(DiffuseParticle *p, Vector3D external_accelerations);
    bool advect_foam(DiffuseParticle *p, Vector3D external_accelerations);
    bool advect_spray(DiffuseParticle *p, Vector3D external_accelerations);
    bool advect_bubbles(DiffuseParticle *p, Vector3D external_accelerations);
    bool dissolve(DiffuseParticle *p);
    void collide(DiffuseParticle *p, vector<CollisionObject *> *collision_objects);
    double k_potential(Particle *p);
    double ta_potential(Particle *p);
    double wc_potential(Particle *p);
    void generate(Particle *p, int n);
	
	double symm_kernel(Vector3D x, double h) {
		if (x.norm() <= h) {
			return 1 - x.norm() / h;
		}
		return 0;
	}

	double cubic_spline_kernel(Vector3D x, double h) {
		double q = x.norm()/h;
		double constant = 1.0 / (pow(h, 3) * PI);
		if (q <= 1) {
			return (1 - 1.5 * pow(q,2) + 0.75 * pow(q, 3)) * constant;
		} else if (q <= 2) {
			return 1.0/4 * pow(2-q, 3) * constant;
		} else {
			return 0;
		}
	}

	double clamp2(double I, double tao_min, double tao_max) {
		return (min(I, tao_max) - min(I, tao_min)) / (tao_max - tao_min);
	}

};





#endif /* PARENT_FLUID_H */
