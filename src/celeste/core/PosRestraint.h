#ifndef __POS_RESTRAINT_H__
#define __POS_RESTRAINT_H__

#include "CelesteObject.h"
#include "PBC.h"

class PRUnit : public CelesteObject {
  private:
  protected:
    int  atomid;
    real crd[3];
    real dist_margin;
    real coef;
    int rest_type;
    //  POSRESTUNIT_NORMAL = 0, POSRESTUNIT_Z = 1
  public:
    PRUnit();
    ~PRUnit();
    int set_parameters(int in_atomid, real crd_x, real crd_y, real crd_z, real in_dist_margin, real in_coef,
		       int rest_type);
    int   get_atomid() { return atomid; };
    real *get_crd() { return crd; };
    real  get_dist_margin() { return dist_margin; };
    real  get_coef() { return coef; };
    real  get_rest_type() { return rest_type; };
};

class PosRestraintObject : public CelesteObject {
  private:
  protected:
    real    weight;
    int     max_n_prunits;
    int     n_prunits;
    PRUnit *prunits;

  public:
    PosRestraintObject();
    ~PosRestraintObject();
    void set_weight(real in_w) { weight = in_w; };
    real                 get_weight() { return weight; };
    int                  get_n_prunits() { return n_prunits; };
    int alloc_prunits(const int in_n);
    int free_prunits();
    int add_prunit(int in_aid, real in_x, real in_y, real in_z, real in_margin, real in_coef, int in_type);
    virtual real_fc apply_restraint(int n_atoms, real **crd, PBC &pbc, real **force);
};

////////////////////////////////////////////////////////////

class PosRestraintHarmonic : public PosRestraintObject {
  private:
  protected:
    // const real boltz;
  public:
    PosRestraintHarmonic();
    ~PosRestraintHarmonic();
    virtual real_fc apply_restraint(int n_atoms, real **crd, PBC &pbc, real **force);
};
////////////////////////////////////////////////////////////
#endif
