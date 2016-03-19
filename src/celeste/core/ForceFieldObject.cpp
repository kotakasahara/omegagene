#include "ForceFieldObject.h"

using namespace std;

ForceFieldObject::ForceFieldObject() : CelesteObject() {
    if (DBG >= 1) cout << "DBG1: ForceFieldObject::ForceFieldObject()" << endl;
    // mmsys = in_mmsys;
    // mmsys->pbc.print_pbc();
}

ForceFieldObject::~ForceFieldObject() {}

int ForceFieldObject::set_config_parameters(const Config *cfg) {
    if (DBG >= 1) cout << "DBG1: ForceFieldObject::set_config_parameters()" << endl;
    cutoff = cfg->cutoff;

    // pbc = in_pbc;
    return 0;
}

int ForceFieldObject::initial_preprocess(const PBC *in_pbc) {
    pbc = in_pbc;

    return 0;
}

int ForceFieldObject::cal_self_energy(const int &n_atoms, const int &n_bonds, const int **&bond_atomid_pairs, const int &n_angles,
                                      const int **&angle_atomid_triads, const int &n_torsions, const int **&torsion_atomid_quads,
                                      const int *&torsion_nb14, const real_pw *&charge, real *&energyself, real &energy_self_sum) {
    return 0;
}
