#ifndef __CONSTRAINT_H__
#define __CONSTRAINT_H__

#include "CelesteObject.h"
#include "PBC.h"
#include <cmath>
using namespace std;

class Constraint : public CelesteObject {
 private:
  
 protected:
  // shk
  static const int Pairs_idx[6][2];
  int    max_loops;
  real   tolerance;
  int    max_n_pair;
  int    n_pair;
  int**  pair_atomids;
  real*  pair_dist;
  int    max_n_trio;
  int    n_trio;
  int**  trio_atomids;
  real** trio_dist;
  int    max_n_quad;
  int    n_quad;
  int**  quad_atomids;
  real** quad_dist;

 public:
  Constraint();
  ~Constraint();
  int alloc_constraint();
  int free_constraint();
  int set_parameters(int in_max_loops, real in_tolerance);
  int set_max_n_constraints(int in_n_pair, int in_n_trio, int in_n_quad);
  int add_pair(int atom1, int atom2, real dist1);
  int add_trio(int atom1, int atom2, int atom3,
	       real dist1, real dist2, real dist3);
  int add_quad(int atom1, int atom2, int atom3, int atom4,
	       real dist1, real dist2, real dist3,
	       real dist4, real dist5, real dist6);
  int apply_constraint(real* in_crd, real* in_crd_prev, real* mass,
		       PBC* pbc);
  int calc_linear_eq(real** a, real* x, real* b, int size);
  int set_subset_constraint(Constraint& super,
			    int* atomids_rev);
  
  int get_n_pair() const { return n_pair; }
  int get_n_trio() const { return n_trio; }
  int get_n_quad() const { return n_quad; }
  const int** get_pair_atomids() const { return (const int**)pair_atomids; };
  const int** get_trio_atomids() const { return (const int**)trio_atomids; };
  const int** get_quad_atomids() const { return (const int**)quad_atomids; };
  const real* get_pair_dist() const { return (const real*)pair_dist; };
  const real** get_trio_dist() const { return (const real**)trio_dist; };
  const real** get_quad_dist() const { return (const real**)quad_dist; };
};

#endif
