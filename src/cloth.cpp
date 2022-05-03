//#include <iostream>
//#include <math.h>
//#include <random>
//#include <vector>
//
//#include "cloth.h"
//#include "collision/plane.h"
//#include "collision/sphere.h"
//
//using namespace std;
//
//Cloth::Cloth(double width, double height, int num_width_points,
//             int num_height_points, double thickness) {
//  this->width = width;
//  this->height = height;
//  this->num_width_points = num_width_points;
//  this->num_height_points = num_height_points;
//  this->thickness = thickness;
//
//  buildGrid();
//  buildClothMesh();
//}
//
//Cloth::~Cloth() {
//  point_masses.clear();
//  springs.clear();
//
//  if (clothMesh) {
//    delete clothMesh;
//  }
//}
//
//void Cloth::buildGrid() {
//  // TODO (Part 1): Build a grid of masses and springs.
//	double x, y, z;
//
//	for (int row = 0; row < num_height_points; row++) {
//		for (int col = 0; col < num_width_points; col++) {
//			x = width * col / (num_width_points - 1);
//			if (orientation == HORIZONTAL) {
//				y = 1;
//				z = height * row / (num_height_points - 1);
//			} else {
//				y = height * row / (num_height_points - 1);
//				z = double(rand()) / RAND_MAX * 2.0 / 1000 - 1.0 / 1000;
//			}
//
//			bool found = false;
//			for (vector<int> v : pinned) {
//				if (v[0] == col && v[1] == row) {
//					found = true;
//					break;
//				}
//			}
//			point_masses.emplace_back(PointMass(Vector3D(x, y, z), found));
//
//		}
//	}
//
//	// Structural
//	// Connect above
//
//	for (int row = 1; row < num_height_points; row++) {
//		for (int col = 0; col < num_width_points; col++) {
//			PointMass *me = &point_masses[row * num_width_points + col];
//			PointMass *above = &point_masses[(row - 1) * num_width_points + col];
//			springs.emplace_back(Spring(me, above, STRUCTURAL));
//		}
//	}
//	// Connect left
//
//	for (int row = 0; row < num_height_points; row++) {
//		for (int col = 1; col < num_width_points; col++) {
//			PointMass *me = &point_masses[row * num_width_points + col];
//			PointMass *left = &point_masses[row * num_width_points + col - 1];
//			springs.emplace_back(Spring(me, left, STRUCTURAL));
//		}
//	}
//
//	// Shearing
//	// Diagonal upper left
//
//	for (int row = 1; row < num_height_points; row++) {
//		for (int col = 1; col < num_width_points; col++) {
//			PointMass *me = &point_masses[row * num_width_points + col];
//			PointMass *up_left = &point_masses[(row - 1) * num_width_points + col - 1];
//			springs.emplace_back(Spring(me, up_left, SHEARING));
//		}
//	}
//	// Diagonal upper right
//
//	for (int row = 1; row < num_height_points; row++) {
//		for (int col = 0; col < num_width_points - 1; col++) {
//			PointMass *me = &point_masses[row * num_width_points + col];
//			PointMass *up_right = &point_masses[(row - 1) * num_width_points + col + 1];
//			springs.emplace_back(Spring(me, up_right, SHEARING));
//		}
//	}
//
//	// Bending
//	// Connect 2 above
//
//	for (int row = 2; row < num_height_points; row++) {
//		for (int col = 0; col < num_width_points; col++) {
//			PointMass *me = &point_masses[row * num_width_points + col];
//			PointMass *above = &point_masses[(row - 2) * num_width_points + col];
//			springs.emplace_back(Spring(me, above, BENDING));
//		}
//	}
//	// Connect 2 left
//
//	for (int row = 0; row < num_height_points; row++) {
//		for (int col = 2; col < num_width_points; col++) {
//			PointMass *me = &point_masses[row * num_width_points + col];
//			PointMass *left = &point_masses[row * num_width_points + col - 2];
//			springs.emplace_back(Spring(me, left, BENDING));
//		}
//	}
//
//}
//
//void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
//                     vector<Vector3D> external_accelerations,
//                     vector<CollisionObject *> *collision_objects) {
//  double mass = width * height * cp->density / num_width_points / num_height_points;
//  double delta_t = 1.0f / frames_per_sec / simulation_steps;
//
//  // TODO (Part 2): Compute total force acting on each point mass.
//	Vector3D total_external_force = Vector3D();
//	for (Vector3D a : external_accelerations) {
//		total_external_force += mass * a;
//	}
//
//	for (auto pm = begin(point_masses); pm != end(point_masses); pm++) {
//		pm->forces = total_external_force;
//	}
//
//	for (auto s = begin(springs); s != end(springs); s++) {
//		if (s->spring_type == BENDING && cp->enable_bending_constraints) {
//			Vector3D f = s->pm_a->position - s->pm_b->position;
//			double fs = 0.2 * cp->ks * (f.norm() - s->rest_length);
//			f.normalize();
//			f *= fs;
//			s->pm_a->forces -= f;
//			s->pm_b->forces += f;
//		} else if (s->spring_type == SHEARING && cp->enable_shearing_constraints) {
//			Vector3D f = s->pm_a->position - s->pm_b->position;
//			double fs = cp->ks * (f.norm() - s->rest_length);
//			f.normalize();
//			f *= fs;
//			s->pm_a->forces -= f;
//			s->pm_b->forces += f;
//		} else if (s->spring_type == STRUCTURAL && cp->enable_structural_constraints) {
//			Vector3D f = s->pm_a->position - s->pm_b->position;
//			double fs = cp->ks * (f.norm() - s->rest_length);
//			f.normalize();
//			f *= fs;
//			s->pm_a->forces -= f;
//			s->pm_b->forces += f;
//		}
//	}
//
//  // TODO (Part 2): Use Verlet integration to compute new point mass positions
//	for (auto pm = begin(point_masses); pm != end(point_masses); pm++) {
//		if (!pm->pinned) {
//			Vector3D xt = pm->position;
//			double d = cp->damping / 100;
//			Vector3D x_prev = pm->last_position;
//			Vector3D at = pm->forces / mass;
//			pm->last_position = pm->position;
//			pm->position = xt + (1 - d) * (xt - x_prev) + at * delta_t * delta_t;
//		}
//	}
//
//  // TODO (Part 4): Handle self-collisions.
//	build_spatial_map();
//	for (auto pm = begin(point_masses); pm != end(point_masses); pm++) {
//		self_collide(*pm, simulation_steps);
//	}
//
//
//  // TODO (Part 3): Handle collisions with other primitives.
//	for (auto c = begin(*collision_objects); c != end(*collision_objects); c++) {
//		for (auto pm = begin(point_masses); pm != end(point_masses); pm++) {
//			(*c)->collide(*pm);
//		}
//	}
//
//
//  // TODO (Part 2): Constrain the changes to be such that the spring does not change
//  // in length more than 10% per timestep [Provot 1995].
//	for (auto s = begin(springs); s != end(springs); s++) {
//		Vector3D a_to_b = s->pm_b->position - s->pm_a->position;
//		double diff = a_to_b.norm() - 1.1 * s->rest_length;
//		if (diff > 0) {
//			a_to_b.normalize();
//			if (!s->pm_a->pinned && !s->pm_b->pinned) {
//				a_to_b *= diff / 2;
//				s->pm_a->position += a_to_b;
//				s->pm_b->position -= a_to_b;
//			} else if (!s->pm_a->pinned) {
//				a_to_b *= diff;
//				s->pm_a->position += a_to_b;
//			} else if (!s->pm_b->pinned) {
//				a_to_b *= diff;
//				s->pm_b->position -= a_to_b;
//			}
//		}
//
//	}
//
//}
//
//void Cloth::build_spatial_map() {
//  for (const auto &entry : map) {
//    delete(entry.second);
//  }
//  map.clear();
//
//  // TODO (Part 4): Build a spatial map out of all of the point masses.
//	for (auto pm = begin(point_masses); pm != end(point_masses); pm++) {
//		double key = hash_position(pm->position);
//		if (map.count(key)) {
//			map[key]->emplace_back(&*pm);
//		} else {
//			vector<PointMass*> *v = new vector<PointMass*>;
//			map[key] = v;
//			v->emplace_back(&*pm);
//		}
//	}
//}
//
//void Cloth::self_collide(PointMass &pm, double simulation_steps) {
//  // TODO (Part 4): Handle self-collision for a given point mass.
//	double key = hash_position(pm.position);
//	Vector3D final_corr = Vector3D();
//	int count = 0;
//	for (auto p = begin(*(map.at(key))); p != end(*(map.at(key))); p++) {
//		if (*p != &pm) {
//			Vector3D diff = pm.position - (*p)->position;
//			double length = diff.norm();
//			if (length < 2 * thickness) {
//				diff.normalize();
//				diff *= (2 * thickness - length);
//				final_corr += diff;
//				count++;
//			}
//		}
//	}
//	if (count) {
//		final_corr /= count;
//		final_corr /= simulation_steps;
//		pm.position += final_corr;
//	}
//
//}
//
//double Cloth::hash_position(Vector3D pos) {
//  // TODO (Part 4): Hash a 3D position into a unique double identifier that represents membership in some 3D box volume.
//	double w = 3 * width / num_width_points;
//	double h = 3 * height / num_height_points;
//	double t = max(w, h);
//
//	int x_box = floor(pos.x / w);
//	int y_box = floor(pos.y / h);
//	int z_box = floor(pos.z / t);
//
//	return x_box + y_box * double(num_width_points) / 3 + z_box * double(num_width_points) / 3 * double(num_height_points) / 3;
//}
//
/////////////////////////////////////////////////////////
///// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
/////////////////////////////////////////////////////////
//
//void Cloth::reset() {
//  PointMass *pm = &point_masses[0];
//  for (int i = 0; i < point_masses.size(); i++) {
//    pm->position = pm->start_position;
//    pm->last_position = pm->start_position;
//    pm++;
//  }
//}
//
//void Cloth::buildClothMesh() {
//  if (point_masses.size() == 0) return;
//
//  ClothMesh *clothMesh = new ClothMesh();
//  vector<Triangle *> triangles;
//
//  // Create vector of triangles
//  for (int y = 0; y < num_height_points - 1; y++) {
//    for (int x = 0; x < num_width_points - 1; x++) {
//      PointMass *pm = &point_masses[y * num_width_points + x];
//      // Get neighboring point masses:
//      /*                      *
//       * pm_A -------- pm_B   *
//       *             /        *
//       *  |         /   |     *
//       *  |        /    |     *
//       *  |       /     |     *
//       *  |      /      |     *
//       *  |     /       |     *
//       *  |    /        |     *
//       *      /               *
//       * pm_C -------- pm_D   *
//       *                      *
//       */
//
//      double u_min = x;
//      u_min /= num_width_points - 1;
//      double u_max = x + 1;
//      u_max /= num_width_points - 1;
//      double v_min = y;
//      v_min /= num_height_points - 1;
//      double v_max = y + 1;
//      v_max /= num_height_points - 1;
//
//      PointMass *pm_A = pm                       ;
//      PointMass *pm_B = pm                    + 1;
//      PointMass *pm_C = pm + num_width_points    ;
//      PointMass *pm_D = pm + num_width_points + 1;
//
//      Vector3D uv_A = Vector3D(u_min, v_min, 0);
//      Vector3D uv_B = Vector3D(u_max, v_min, 0);
//      Vector3D uv_C = Vector3D(u_min, v_max, 0);
//      Vector3D uv_D = Vector3D(u_max, v_max, 0);
//
//
//      // Both triangles defined by vertices in counter-clockwise orientation
//      triangles.push_back(new Triangle(pm_A, pm_C, pm_B,
//                                       uv_A, uv_C, uv_B));
//      triangles.push_back(new Triangle(pm_B, pm_C, pm_D,
//                                       uv_B, uv_C, uv_D));
//    }
//  }
//
//  // For each triangle in row-order, create 3 edges and 3 internal halfedges
//  for (int i = 0; i < triangles.size(); i++) {
//    Triangle *t = triangles[i];
//
//    // Allocate new halfedges on heap
//    Halfedge *h1 = new Halfedge();
//    Halfedge *h2 = new Halfedge();
//    Halfedge *h3 = new Halfedge();
//
//    // Allocate new edges on heap
//    Edge *e1 = new Edge();
//    Edge *e2 = new Edge();
//    Edge *e3 = new Edge();
//
//    // Assign a halfedge pointer to the triangle
//    t->halfedge = h1;
//
//    // Assign halfedge pointers to point masses
//    t->pm1->halfedge = h1;
//    t->pm2->halfedge = h2;
//    t->pm3->halfedge = h3;
//
//    // Update all halfedge pointers
//    h1->edge = e1;
//    h1->next = h2;
//    h1->pm = t->pm1;
//    h1->triangle = t;
//
//    h2->edge = e2;
//    h2->next = h3;
//    h2->pm = t->pm2;
//    h2->triangle = t;
//
//    h3->edge = e3;
//    h3->next = h1;
//    h3->pm = t->pm3;
//    h3->triangle = t;
//  }
//
//  // Go back through the cloth mesh and link triangles together using halfedge
//  // twin pointers
//
//  // Convenient variables for math
//  int num_height_tris = (num_height_points - 1) * 2;
//  int num_width_tris = (num_width_points - 1) * 2;
//
//  bool topLeft = true;
//  for (int i = 0; i < triangles.size(); i++) {
//    Triangle *t = triangles[i];
//
//    if (topLeft) {
//      // Get left triangle, if it exists
//      if (i % num_width_tris != 0) { // Not a left-most triangle
//        Triangle *temp = triangles[i - 1];
//        t->pm1->halfedge->twin = temp->pm3->halfedge;
//      } else {
//        t->pm1->halfedge->twin = nullptr;
//      }
//
//      // Get triangle above, if it exists
//      if (i >= num_width_tris) { // Not a top-most triangle
//        Triangle *temp = triangles[i - num_width_tris + 1];
//        t->pm3->halfedge->twin = temp->pm2->halfedge;
//      } else {
//        t->pm3->halfedge->twin = nullptr;
//      }
//
//      // Get triangle to bottom right; guaranteed to exist
//      Triangle *temp = triangles[i + 1];
//      t->pm2->halfedge->twin = temp->pm1->halfedge;
//    } else {
//      // Get right triangle, if it exists
//      if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
//        Triangle *temp = triangles[i + 1];
//        t->pm3->halfedge->twin = temp->pm1->halfedge;
//      } else {
//        t->pm3->halfedge->twin = nullptr;
//      }
//
//      // Get triangle below, if it exists
//      if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
//        Triangle *temp = triangles[i + num_width_tris - 1];
//        t->pm2->halfedge->twin = temp->pm3->halfedge;
//      } else {
//        t->pm2->halfedge->twin = nullptr;
//      }
//
//      // Get triangle to top left; guaranteed to exist
//      Triangle *temp = triangles[i - 1];
//      t->pm1->halfedge->twin = temp->pm2->halfedge;
//    }
//
//    topLeft = !topLeft;
//  }
//
//  clothMesh->triangles = triangles;
//  this->clothMesh = clothMesh;
//}
