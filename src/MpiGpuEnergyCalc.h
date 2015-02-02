#ifndef __MPI_GPU_ENERGY_CALC_H__
#define __MPI_GPU_ENERGY_CALC_H__

#include "GpuEnergyCalc.h"
#include <mpi.h>
using namespace std;

class MpiGpuEnergyCalc : public GpuEnergyCalc{
 private:
 protected:
  int n_threads;
  real_pw *d_crd;
  int *d_grid_atom_index;
  int *d_grid_n_atoms;

  // Information for each MPI process
  int offset_gridpairs;
  int n_cal_gridpairs;
  int offset_grids;
  int n_cal_grids;

  real* buf_work;

  int mpi_myrank;
  int mpi_n_ranks;

 public:
  
  MpiGpuEnergyCalc(MmSystem* in_mmsys);
  ~MpiGpuEnergyCalc();
  int initial_preprocess();
  int calc_energy();
  int calc_energy_pairwise();
  int finish_energy_pairwise();
  int set_dev_pointers(real_pw *in_d_crd,
		       int *in_d_grid_atom_index,
		       int *in_d_grid_n_atoms);
  // Information for each MPI process
  int set_range_in_mpi_process();
  
  int alloc_buffers();
  int free_buffers();
};

#endif
