#include "COMMotion.h"

COMMotion::COMMotion() {}

COMMotion::~COMMotion() {}

int COMMotion::set_groups(const int in_n_mod_groups,
                          int *     in_mod_group_ids,
                          int *     in_n_atoms_in_groups,
                          int **    in_groups,
                          real_pw * in_mass_inv_groups,
                          real_pw * in_mass) {
    n_mod_groups      = in_n_mod_groups;
    mod_group_ids     = in_mod_group_ids;
    n_atoms_in_groups = in_n_atoms_in_groups;
    groups            = in_groups;
    mass_inv_groups   = in_mass_inv_groups;
    mass              = in_mass;

    // cout << "Canceling the center-of-mass motion" << endl;
    // for (int grp=0; grp < n_groups; grp++){
    // cout << "Group : " << group_ids[grp]  << " "
    //<< n_atoms_in_groups[group_ids[grp]] << " atoms, "
    //<< "mass: " << 1.0/mass_inv_groups[group_ids[grp]] << endl;
    //}

    return 0;
}

int COMMotion::cancel_translation(int *atomids_rev, real *vel_next) {
    for (int i_grp = 0; i_grp < n_mod_groups; i_grp++) {
        int grp_id = mod_group_ids[i_grp];
        // real center[3] = {0.0, 0.0, 0.0};
        real moment[3] = {0.0, 0.0, 0.0};
        for (int i_atom = 0; i_atom < n_atoms_in_groups[grp_id]; i_atom++) {
            int idx   = atomids_rev[groups[grp_id][i_atom]];
            int idx_3 = idx * 3;
            for (int d = 0; d < 3; d++) {
                // center[d] += crd[idx_3+d] * mass[idx];
                moment[d] += vel_next[idx_3 + d] * mass[idx];
            }
        }
        for (int d = 0; d < 3; d++) {
            // center[d] *= mass_inv_groups[grp_id];
            moment[d] *= mass_inv_groups[grp_id];
        }
        for (int i_atom = 0; i_atom < n_atoms_in_groups[grp_id]; i_atom++) {
            int idx   = atomids_rev[groups[grp_id][i_atom]];
            int idx_3 = idx * 3;
            for (int d = 0; d < 3; d++) { vel_next[idx_3 + d] -= moment[d]; }
        }
    }

    return 0;
}
