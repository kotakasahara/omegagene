#ifndef __FORCE_FIELD_OBJECT_H__
#define __FORCE_FIELD_OBJECT_H__

#include "CelesteObject.h"
#include "Config.h"
#include "PBC.h"

class ForceFieldObject : public CelesteObject {
  private:
  protected:
    real       cutoff;
    const PBC *pbc;
    // MmSystem* mmsys;

  public:
    ForceFieldObject();
    ~ForceFieldObject();
    virtual int set_config_parameters(const Config *cfg);
    virtual int initial_preprocess(const PBC *in_pbc);
    virtual int cal_self_energy(const int &     n_atoms,
				const int &  n_excess,
				const int **&excess_pairs,
                                //const int &     n_bonds,
                                //const int **&   bond_atomid_pairs,
                                //const int &     n_angles,
                                //const int **&   angle_atomid_triads,
                                //const int &     n_torsions,
                                //const int **&   torsion_atomid_quads,
                                //const int *&    torsion_nb14,
                                real_pw *&charge,
                                real *&         energy_self,
                                real &          energy_self_sum);
};

#endif
