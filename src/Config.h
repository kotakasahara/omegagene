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
  enum { MAX_N_COM_GROUPS = 32 };
  int mode;
  string fn_cfg;
  string fn_inp;
  int processor;
  int integrator_type;
  int constraint_type;
  real constraint_tolerance;
  int constraint_max_loops;
  int thermostat_type;
  real temperature;
  real thermo_const_tolerance;
  int thermo_const_max_loops;

  real_pw cutoff;
  real_pw cutoff_buf;
  int n_steps;
  real time_step;
  int electrostatic;
  real ele_alpha;
  int center_of_motion;
  int n_com_cancel_groups;
  int com_cancel_groups[MAX_N_COM_GROUPS];
  int n_com_cancel_groups_name;
  string com_cancel_groups_name[MAX_N_COM_GROUPS];
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

  string fn_o_restart;
  string fn_o_crd;
  int format_o_crd;
  string fn_o_log;
  string fn_o_energy;
  string fn_o_vmcmd_log;
  string fn_o_expand_lambda;
  int format_o_expand_lambda;
  string fn_o_energyflow;

  int init_vel_just;
  // 0: initial velocity is 0-dt
  // 1: initial velocity is 0

  real_pw nsgrid_cutoff;
  //real nsgrid_min_width;
  //real nsgrid_max_n_atoms;
  int nsgrid_update_intvl;

  int dist_restraint_type;
  real dist_restraint_weight;

  Config();
  ~Config();
  int set_defaults();
  void setAll(int argn,char* argv[]);
  void setAll(vector<string> arg);
  void operator=(Config op);
};

#endif
