#include "PBC.h"

using namespace std;

PBC::PBC() : CelesteObject() {
    for (int i = 0; i < 3; i++) {
        L[i]           = 0.0;
        angle[i]       = 0.0;
        L_half[i]      = 0.0;
        L_inv[i]       = 0.0;
        lower_bound[i] = 0.0;
        upper_bound[i] = 0.0;
    }
}

PBC::~PBC() {}

int PBC::set_pbc(real val[]) {
    if (DBG >= 1) {
        // cout << "PBC::set_pbc()"<<endl;
        for (int i = 0; i < 12; i++) {
            // cout << val[i] << " ";
        }
    }
    L[0]     = val[0];
    L[1]     = val[4];
    L[2]     = val[8];
    angle[0] = 90.0;
    angle[1] = 90.0;
    angle[2] = 90.0;

    for (int i = 0; i < 3; i++) {
        L_half[i]      = L[i] * 0.5;
        L_inv[i]       = 1.0 / L[i];
        L_half_inv[i]  = 1.0 / L_half[i];
	lower_bound[i] = val[9 + i] - L_half[i];
        upper_bound[i] = lower_bound[i] + L[i];
        cout << "L[i]: " << L[i] << endl;
    }
    return 0;
}
// template <typename TYPE> int PBC::diff_crd_minim_image(TYPE d[], const real crd1[], const real crd2[]) const{

int PBC::diff_crd_minim_image(float d[], const float crd1[], const float crd2[]) const {
    for (int i = 0; i < 3; i++) {
        d[i] = crd1[i] - crd2[i];
        d[i] = d[i] - L[i] * (int)(d[i] * L_half_inv[i]);
        d[i] = d[i] - L[i] * (int)(d[i] * L_half_inv[i]);
    }
    return 0;
}
int PBC::diff_crd_minim_image(double d[], const float crd1[], const float crd2[]) const {
    for (int i = 0; i < 3; i++) {
        d[i] = crd1[i] - crd2[i];
        d[i] = d[i] - L[i] * (int)(d[i] * L_half_inv[i]);
        d[i] = d[i] - L[i] * (int)(d[i] * L_half_inv[i]);
    }
    return 0;
}
int PBC::diff_crd_minim_image(float d[], const double crd1[], const double crd2[]) const {
    for (int i = 0; i < 3; i++) {
        d[i] = crd1[i] - crd2[i];
        d[i] = d[i] - L[i] * (int)(d[i] * L_half_inv[i]);
        d[i] = d[i] - L[i] * (int)(d[i] * L_half_inv[i]);
    }
    return 0;
}

int PBC::diff_crd_minim_image(double d[], const double crd1[], const double crd2[]) const {
    for (int i = 0; i < 3; i++) {
        d[i] = crd1[i] - crd2[i];
        // if (d[i] > L_half[i]) d[i] -= L[i];
        // else if (d[i] < -L_half[i]) d[i] += L[i];

        // d[i] *= L_inv[i] + m;
        // int k = int(d[i] + 0.5);
        // d[i] -= L[i] * k;
        d[i] = d[i] - L[i] * (int)(d[i] * L_half_inv[i]);
        d[i] = d[i] - L[i] * (int)(d[i] * L_half_inv[i]);
    }
    return 0;
}
int PBC::mid_crd_minim_image(float d[], const real crd1[], const real crd2[]) const {
    float diff[3];
    diff_crd_minim_image(diff, crd1, crd2);
    for (int i = 0; i < 3; i++) {
        d[i] = crd1[i] + diff[i] * 0.5;
        if (d[i] < lower_bound[i])
            d[i] += L[i];
        else if (d[i] >= upper_bound[i])
            d[i] -= L[i];
    }
    return 0;
}

int PBC::print_pbc() {
    cout << "PBC: " << L[0] << " , " << L[1] << " , " << L[2] << endl;
    cout << "PBC-lower_bound: " << lower_bound[0] << " , " << lower_bound[1] << " , " << lower_bound[2] << endl;
    return 0;
}

real PBC::cal_volume() {
    return L[0] * L[1] * L[2];
}

int PBC::fix_pbc_image(float *crd, const int image) {
    if (image == 0) return 1;
    if ((image & 1) == 1)
        crd[0] -= L[0];
    else if ((image & 2) == 2)
        crd[0] += L[0];
    if ((image & 4) == 4)
        crd[1] -= L[1];
    else if ((image & 8) == 8)
        crd[1] += L[1];
    if ((image & 16) == 16)
        crd[2] -= L[2];
    else if ((image & 32) == 32)
        crd[2] += L[2];
    return 0;
}
int PBC::fix_pbc_image(double *crd, const int image) {
    if (image == 0) return 1;
    if ((image & 1) == 1)
        crd[0] -= L[0];
    else if ((image & 2) == 2)
        crd[0] += L[0];
    if ((image & 4) == 4)
        crd[1] -= L[1];
    else if ((image & 8) == 8)
        crd[1] += L[1];
    if ((image & 16) == 16)
        crd[2] -= L[2];
    else if ((image & 32) == 32)
        crd[2] += L[2];
    return 0;
}
