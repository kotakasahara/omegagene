#ifndef __CONSTRAINT_H__
#define __CONSTRAINT_H__

#include "CelesteObject.h"
#include "PBC.h"
#include <cmath>

class ConstraintObject : public CelesteObject {
  private:
  protected:
    // shk
    static const int Pairs_idx[6][2];
    int              max_loops;
    real_cst         tolerance;
    int              max_n_pair;
    int              n_pair;
    int **           pair_atomids;
    real_cst *       pair_dist;
    int              max_n_trio;
    int              n_trio;
    int **           trio_atomids;
    real_cst **      trio_dist;
    int              max_n_quad;
    int              n_quad;
    int **           quad_atomids;
    real_cst **      quad_dist;

  public:
    ConstraintObject();
    ~ConstraintObject();
    int alloc_constraint();
    int free_constraint();
    int set_parameters(int in_max_loops, real_cst in_tolerance);
    int set_max_n_constraints(int in_n_pair, int in_n_trio, int in_n_quad);
    int add_pair(int atom1, int atom2, real_cst dist1);
    int add_trio(int atom1, int atom2, int atom3, real_cst dist1, real_cst dist2, real_cst dist3);
    int add_quad(int      atom1,
                 int      atom2,
                 int      atom3,
                 int      atom4,
                 real_cst dist1,
                 real_cst dist2,
                 real_cst dist3,
                 real_cst dist4,
                 real_cst dist5,
                 real_cst dist6);
    virtual int apply_constraint(real *in_crd, real *in_crd_prev, real_pw *mass_inv, PBC *pbc);
    int calc_linear_eq(real_cst a[6][6], real_cst x[6], real_cst b[6], int size);
    int set_subset_constraint(ConstraintObject &super, int *atomids_rev);

    int              get_n_pair() const { return n_pair; }
    int              get_n_trio() const { return n_trio; }
    int              get_n_quad() const { return n_quad; }
    const int **     get_pair_atomids() const { return (const int **)pair_atomids; };
    const int **     get_trio_atomids() const { return (const int **)trio_atomids; };
    const int **     get_quad_atomids() const { return (const int **)quad_atomids; };
    const real_cst * get_pair_dist() const { return (const real_cst *)pair_dist; };
    const real_cst **get_trio_dist() const { return (const real_cst **)trio_dist; };
    const real_cst **get_quad_dist() const { return (const real_cst **)quad_dist; };
};

#endif
