#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <vector>
#include <set>
#include <cstdlib>
using namespace std;

#include "CelesteObject.h"

class Config : public CelesteObject {
 private:
 public:
  int mode;
  string fn_cfg;
  string fn_inp;
  int processor;
  int integrator;
  int constraint;
  int constraint_tolerance;
  int constraint_max_loops;
  real cutoff;
  real cutoff_buf;
  int n_steps;
  real time_step;
  int electrostatic;
  real ele_alpha;
  int thermostat;
  real temperature;
  int center_of_motion;
  int random_seed;
  int expanded_ensemble;

  int box_div[3];
  
  int print_intvl_crd;
  int print_intvl_vel;
  int print_intvl_log;
  int print_intvl_force;
  int print_intvl_energy;
  int print_intvl_energyflow;
  int print_intvl_expand_lambda;

  string fn_o_crd;
  string fn_o_log;
  string fn_o_energy;
  string fn_o_vmcmd_log;
  string fn_o_expand_lambda;
  string fn_o_energyflow;

  real nsgrid_cutoff;
  //real nsgrid_min_width;
  //real nsgrid_max_n_atoms;
  int nsgrid_update_intvl;

  Config();
  ~Config();
  int set_defaults();
  void setAll(int argn,char* argv[]);
  void setAll(vector<string> arg);
  void operator=(Config op);
};

#endif
