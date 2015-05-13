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
  real L[3];
  real angle[3];
  real L_half[3];
  real L_inv[3];
  real L_half_inv[3];
  real lower_bound[3];
  real upper_bound[3];

  PBC();
  ~PBC();
  int set_pbc(real val[]);
  int diff_crd_minim_image(float d[], const real crd1[], const real crd2[]) const;
  int diff_crd_minim_image(double d[], const real crd1[], const real crd2[]) const;
  int mid_crd_minim_image(float d[], const real crd1[], const real crd2[]) const;
  int print_pbc();
  real cal_volume();
  int fix_pbc_image(real* crd, const int image);
};

#endif
