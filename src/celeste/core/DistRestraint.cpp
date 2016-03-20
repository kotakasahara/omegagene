#include "DistRestraint.h"

DRUnit::DRUnit() {}
DRUnit::~DRUnit() {}
int DRUnit::set_parameters(int  in_aid1,
                           int  in_aid2,
                           real in_coef_low,
                           real in_coef_high,
                           real in_dist_low,
                           real in_dist_high) {
    atomid1   = in_aid1;
    atomid2   = in_aid2;
    coef_low  = in_coef_low;
    coef_high = in_coef_high;
    dist_low  = in_dist_low;
    dist_high = in_dist_high;
    return 0;
}

////////////////////////////////////////////////////////

DistRestraintObject::DistRestraintObject() {
    n_drunits     = 0;
    max_n_drunits = 0;
    // boltz = -GAS_CONST*FORCE_VEL;
}
DistRestraintObject::~DistRestraintObject() {
    free_drunits();
}
int DistRestraintObject::alloc_drunits(const int in_n) {
    max_n_drunits = in_n;
    drunits       = new DRUnit[max_n_drunits];
    return 0;
}
int DistRestraintObject::free_drunits() {
    delete[] drunits;
    return 0;
}
int DistRestraintObject::add_drunit(int  in_aid1,
                                    int  in_aid2,
                                    real in_coef_low,
                                    real in_coef_high,
                                    real in_dist_low,
                                    real in_dist_high) {
    drunits[n_drunits].set_parameters(in_aid1, in_aid2, in_coef_low, in_coef_high, in_dist_low, in_dist_high);
    n_drunits++;
    return n_drunits;
}
real_fc DistRestraintObject::apply_restraint(int n_atoms, real **crd, PBC &pbc, real **force) {
    return 0;
}

///////////////////////////////////////////////////////////////

DistRestraintHarmonic::DistRestraintHarmonic() : DistRestraintObject() {}
DistRestraintHarmonic::~DistRestraintHarmonic() {
    free_drunits();
}

real_fc DistRestraintHarmonic::apply_restraint(int n_atoms, real **crd, PBC &pbc, real **force) {

    for (int i = 0; i < n_atoms; i++) {
        for (int d = 0; d < 3; d++) { force[i][d] = 0.0; }
    }
    real_fc ene = 0.0;
    for (int i = 0; i < n_drunits; i++) {
        real diff[3];
        pbc.diff_crd_minim_image(diff, crd[drunits[i].get_atomid1()], crd[drunits[i].get_atomid2()]);
        real r_sq = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
        real r    = sqrt(r_sq);

        real coef;
        real ref_dist;
        if (r < drunits[i].get_dist_low()) {
            coef     = drunits[i].get_coef_low();
            ref_dist = drunits[i].get_dist_low();
        } else if (r > drunits[i].get_dist_high()) {
            coef     = drunits[i].get_coef_high();
            ref_dist = drunits[i].get_dist_high();
        } else
            continue;

        real k = weight * coef;

        real dist_diff = r - ref_dist;
        ene += k * dist_diff * dist_diff;

        // force
        real    k_g = 2.0 * k * dist_diff;
        real_fc frc[3];
        for (int d = 0; d < 3; d++) {
            // force[drunits[i].atomid1] += diff[d] / r;
            force[drunits[i].get_atomid1()][d] -= k_g * diff[d] / r;
            force[drunits[i].get_atomid2()][d] += k_g * diff[d] / r;
            // frc[d] += diff[d] / r;
        }
    }
    return ene;
}
