#ifndef __GPU_DYNAMICS_MODE_H__
#define __GPU_DYNAMICS_MODE_H__

#include "DynamicsMode.h"
#include "GpuEnergyCalc.h"
//#include <chrono>

class GpuDynamicsMode : public DynamicsMode {
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
  // Pointers for device memory

  //int* d_nb15off;
  //int* d_n_nb15off;

  // Information about sizes
  int n_atoms;
  //int n_all_grids;
  //int n_grid_pairs;


  
 public:
  GpuDynamicsMode();
  ~GpuDynamicsMode();  
  int initial_preprocess();
  int update_device_cell_info();
  int terminal_process();
  int set_config_parameters(Config* in_cfg);
  int main_stream();
  int calc_in_each_step();
};

#endif
