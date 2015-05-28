#ifndef __DYNAMICS_MODE_H__
#define __DYNAMICS_MODE_H__

#include <sstream>
#include <cstdio>

#include "RunMode.h"
#include "MmSystem.h"
#include "SubBox.h"

class DynamicsMode : public RunMode {
 private:
  
 protected:
  real time_step;
  real time_step_half;
  real temperature;
  int integrator;
  int barostat;
  real pressure;
  real temperature_coeff;
  real nsgrid_cutoff;
  real nsgrid_min_width;
  int nsgrid_max_n_atoms;
  //super: int n_steps;
  //super: int integrator;
  //super: int electrostatic;
  //super: int center_of_motion;
  //super: int print_intvl_crd;
  //super: int print_intvl_vel;
  //super: int print_intvl_log;
  //super: int print_intvl_energy;
  //super: int print_intvl_energyflow;
  //super: string fn_o_crd;
  //super: string fn_o_log;
  //super: string fn_o_energy;
  //super: string fn_o_energyflow;
  SubBox subbox;

 public:
  // super: MmSystem mmsys;
  // super: EnergyCalcObject* enecal;
  DynamicsMode();
  ~DynamicsMode();
  int test(Config* in_cfg);
  int set_config_parameters(Config* in_cfg);
  int initial_preprocess();
  int terminal_process();
  int main_stream();
  virtual int calc_in_each_step();
  virtual int apply_constraint();
  int apply_dist_restraint();
  int sub_output();
  int sub_output_log();
  int cal_kinetic_energy(const real** vel);
  int subbox_setup();
  int subbox_set_bonding_potentials();
  int gather_energies();
  int output_restart();
};

class DynamicsModePresto : public DynamicsMode {
 private:
  
 protected:
 public:
  DynamicsModePresto();
  ~DynamicsModePresto();
  virtual int calc_in_each_step();
  virtual int apply_constraint();

};
class DynamicsModeZhang : public DynamicsMode {
 private:
  
 protected:
 public:
  DynamicsModeZhang();
  ~DynamicsModeZhang();
  virtual int calc_in_each_step();
  virtual int apply_constraint();

};

#endif
