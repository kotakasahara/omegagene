#ifndef __PBC_H__
#define __PBC_H__

#include "CelesteObject.h"
#include <cmath>
#include <iostream>

/*
References:

Minimum image convention:

Hloucha M, Deiters UK (1998) Fast coding of the minimum image convention. Molecular Simulation 20: 239–244.

 */

class PBC : public CelesteObject{
 private:
 protected:
 public:
  real_pw L[3];
  real_pw angle[3];
  real_pw L_half[3];
  real_pw L_inv[3];
  real_pw L_half_inv[3];
  real_pw lower_bound[3];
  real_pw upper_bound[3];

  PBC();
  ~PBC();
  int set_pbc(real val[]);
  int diff_crd_minim_image(float d[], const float crd1[], const float crd2[]) const;
  int diff_crd_minim_image(double d[], const float crd1[], const float crd2[]) const;
  int diff_crd_minim_image(float d[], const double crd1[], const double crd2[]) const;
  int diff_crd_minim_image(double d[], const double crd1[], const double crd2[]) const;
  int mid_crd_minim_image(float d[], const real crd1[], const real crd2[]) const;
  int print_pbc();
  real cal_volume();
  int fix_pbc_image(float* crd, const int image);
  int fix_pbc_image(double* crd, const int image);
};

#endif
