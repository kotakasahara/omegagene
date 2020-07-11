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
    int n_params;
    real params[MAX_N_POSRES_PARAMS];
    //  POSRESTUNIT_NORMAL = 0, POSRESTUNIT_Z = 1
    //  POSRESTUNIT_MULTIWELL01 = 2
  public:
    PRUnit();
    ~PRUnit();
    int set_parameters(int in_atomid, real crd_x, real crd_y, real crd_z, real in_dist_margin, real in_coef,
		       int rest_type, int in_n_params, real* in_params);
    int   get_atomid() { return atomid; };
    real *get_crd() { return crd; };
    real  get_dist_margin() { return dist_margin; };
    real  get_coef() { return coef; };
    real  get_rest_type() { return rest_type; };
    int get_n_params(){ return n_params; };
    real get_params(int i){ return params[i]; };
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
    int add_prunit(int in_aid, real in_x, real in_y, real in_z, real in_margin, real in_coef, int in_type, int in_n_params, real *in_params);
    virtual real_fc apply_restraint(int n_atoms, real *crd, PBC &pbc, real **force);
};

////////////////////////////////////////////////////////////

class PosRestraintHarmonic : public PosRestraintObject {
  private:
  protected:
    // const real boltz;
  public:
    PosRestraintHarmonic();
    ~PosRestraintHarmonic();
    virtual real_fc apply_restraint(int n_atoms, real *crd, PBC &pbc, real **force);
};
////////////////////////////////////////////////////////////
#endif
