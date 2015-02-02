#ifndef __MPI_GPU_DYNAMICS_MODE_H__
#define __MPI_GPU_DYNAMICS_MODE_H__

#include <mpi.h>
#include "GpuDynamicsMode.h"
#include "GpuEnergyCalc.h"

//#include <chrono>

class MpiGpuDynamicsMode : public GpuDynamicsMode {
 private:
 protected:
  //super: real time_step;
  //super: int thermostat;
  //super: real temperature;
  //super: int integrator;
  //super: int barostat;
  //super: real pressure;
  //super: real temperature_coeff;
  //super: real nsgrid_cutoff;
  //super: real nsgrid_min_width;
  //super: int nsgrid_max_n_atoms;
  //super2: int n_steps;
  //super2: int integrator;
  //super2: int electrostatic;
  //super2: int center_of_motion;
  //super2: int print_intvl_crd;
  //super2: int print_intvl_vel;
  //super2: int print_intvl_log;
  //super2: int print_intvl_energy;
  //super2: int print_intvl_energyflow;
  //super2: string fn_o_crd;
  //super2: string fn_o_log;
  //super2: string fn_o_energy;
  //super2: string fn_o_energyflow;
  
  // Cuda Interface
  // Information about sizes
  //int n_atoms;
  //int n_all_grids;
  //int n_grid_pairs;

  int n_processes;
  int myrank;
 public:
  MpiGpuDynamicsMode();
  ~MpiGpuDynamicsMode();  
  int initial_preprocess();
  int main_stream();
  int broadcast_each_step();
};

#endif
