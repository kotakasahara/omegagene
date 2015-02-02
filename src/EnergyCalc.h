#ifndef __ENERGY_CALC_H__
#define __ENERGY_CALC_H__

#include <set>
#include <algorithm>
#include <cmath>
using namespace std;

#include "MmSystem.h"
#include "EnergyCalcObject.h"
#include "ForceField.h"

class EnergyCalc : public EnergyCalcObject{
 private:
 protected:
  ForceField* ff;
 public:
  EnergyCalc(MmSystem* in_mmsys,
	     SubBox* in_subbox);
  int set_config_parameters(const Config* in_cfg);
  int set_dev_pointers(real_pw *in_d_crd,
		       int *in_d_grid_atom_index,
		       int *in_d_grid_n_atoms);
  int initial_preprocess();
  int calc_energy();
  int calc_energy_bonds();
  int calc_energy_angles();
  int calc_energy_torsions();
  int calc_energy_impros();
  int calc_energy_14nb();
  int calc_energy_pairwise_nogrid();
  int calc_energy_pairwise();
  int calc_energy_ele_excess();
  bool check_nb15off(const int& a1, const int& a2, const int* bitmask);
  int fix_pbc_image(real* crd, const int image);
};

#endif 
