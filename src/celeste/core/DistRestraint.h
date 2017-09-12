#ifndef __DIST_RESTRAINT_H__
#define __DIST_RESTRAINT_H__

#include "CelesteObject.h"
#include "PBC.h"

class DRUnit : public CelesteObject {
  private:
  protected:
    int  atomid1, atomid2;
    real coef_low, coef_high;
    real dist_low, dist_high;

  public:
    DRUnit();
    ~DRUnit();
    int set_parameters(int  in_aid1,
                       int  in_aid2,
                       real in_coef_low,
                       real in_coef_high,
                       real in_dist_low,
                       real in_dist_high);
    int  get_atomid1() { return atomid1; };
    int  get_atomid2() { return atomid2; };
    real get_coef_low() { return coef_low; };
    real get_coef_high() { return coef_high; };
    real get_dist_low() { return dist_low; };
    real get_dist_high() { return dist_high; };
};

class DistRestraintObject : public CelesteObject {
  private:
  protected:
    real    weight;
    int     max_n_drunits;
    int     n_drunits;
    DRUnit *drunits;

  public:
    DistRestraintObject();
    virtual ~DistRestraintObject();
    void set_weight(real in_w) { weight = in_w; };
    real                 get_weight() { return weight; };
    int                  get_n_drunits() { return n_drunits; };
    int alloc_drunits(const int in_n);
    int free_drunits();
    int add_drunit(int in_aid1, int in_aid2, real in_coef_low, real in_coef_high, real in_dist_low, real in_dist_high);
    virtual real_fc apply_restraint(int n_atoms, real **crd, PBC &pbc, real **force);
};

////////////////////////////////////////////////////////////

class DistRestraintHarmonic : public DistRestraintObject {
  private:
  protected:
    // const real boltz;
  public:
    DistRestraintHarmonic();
    ~DistRestraintHarmonic();
    virtual real_fc apply_restraint(int n_atoms, real **crd, PBC &pbc, real **force);
};
#endif
