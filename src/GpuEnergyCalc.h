#ifndef __GPU_ENERGY_CALC_H__
#define __GPU_ENERGY_CALC_H__

#include "EnergyCalc.h"
using namespace std;

class GpuEnergyCalc : public EnergyCalc{
 private:
 protected:
  int n_threads;
  real_pw *d_crd;
  int *d_grid_atom_index;
  int *d_grid_n_atoms;
 public:
  GpuEnergyCalc(MmSystem* in_mmsys);
  ~GpuEnergyCalc();
  int initial_preprocess();
  int calc_energy(bool grid_update);
  int calc_energy_pairwise();
  int finish_energy_pairwise();
  int set_dev_pointers(real_pw *in_d_crd,
		       int *in_d_grid_atom_index,
		       int *in_d_grid_n_atoms);
};

#endif
