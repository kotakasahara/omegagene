#include "define.h"
#include "Config.h"
#include <iostream>
using namespace std;


Config::Config(){
  mode=M_TEST;
  fn_cfg="";
}

Config::~Config(){
}
int Config::set_defaults(){
  mode = M_TEST;
  fn_cfg = "md_i.cfg";
  fn_inp = "md_i.inp";
  gpu_device_id = -1;
  processor = PRCS_SINGLE;
  integrator_type = INTGRTR_LEAPFROG_PRESTO;
  constraint_type = CONST_NONE;
  cutoff = 12.0;
  n_steps = 1;
  time_step = 0.0005;
  electrostatic =  ELCTRST_ZERODIPOLE;
  ele_alpha = 0.0;
  thermostat_type = THMSTT_NONE;
  temperature = 300;

  com_motion = COM_NONE;
  n_com_cancel_groups = 0;
  n_com_cancel_groups_name = 0;
  
  box_div[0] = 1;
  box_div[1] = 1;
  box_div[2] = 1;

  print_intvl_crd = 10000;
  print_intvl_vel = 0;
  print_intvl_force = 0;
  fn_o_restart = "md_o.restart";
  fn_o_crd = "md_o.trr";
  format_o_crd = CRDOUT_GROMACS;
  fn_o_log = "md_o.log";
  fn_o_energy = "md_o.erg";
  fn_o_energyflow = "md_o.efl";
  nsgrid_cutoff = cutoff + 1.0;
  //nsgrid_min_width = cutoff * 0.5;
  //nsgrid_max_n_atoms = 100;
  nsgrid_update_intvl = 1;

  constraint_tolerance = 0.000001;
  constraint_max_loops = 1000;
  thermo_const_tolerance = 1000;
  thermo_const_max_loops = 0.000001;

  init_vel_just = 0;
  expanded_ensemble = EXPAND_NONE;
  format_o_expand_lambda = LAMBDAOUT_BIN;

  dist_restraint_type = DISTREST_NONE;
  dist_restraint_weight = 0.0;

  return 0;
}
void Config::setAll(int argn, char* argv[]){
  vector<string> arg;
  int i;
  for(i=1;i<argn;i++)
    arg.push_back(string(argv[i]));
  setAll(arg);
}
void Config::setAll(vector<string> arg){
  vector<string>::iterator itr;
  string type,val;
  for(itr=arg.begin(); itr!=arg.end(); itr++){
    if(*itr=="--mode"){
      itr++;
      if(*itr=="test")        { mode=M_TEST; }
      else if(*itr=="md")     { mode=M_DYNAMICS; }
      else{
	cerr<<"invalid mode ["<<(*itr)<<"]\n"; exit(1);
      }
    }

    else if(*itr=="--cfg"){ fn_cfg = *++itr; }
    else if(*itr=="--inp"){ fn_inp = *++itr; }
    else if(*itr=="--processor"){ 
      itr++;
      if(*itr=="single"){ processor = PRCS_SINGLE; }
      else if(*itr=="mpi"){ processor = PRCS_MPI; }
      else if(*itr=="cuda"){ processor = PRCS_CUDA; }
      else if(*itr=="mpi-cuda"){ processor = PRCS_MPI_CUDA; }
      else{ processor = PRCS_DUMMY; }
    }
    else if(*itr=="--integrator"){ 
      itr++;
      if(*itr=="zhang"){ integrator_type = INTGRTR_ZHANG; }
      else if(*itr=="leapfrog-presto"){ integrator_type = INTGRTR_LEAPFROG_PRESTO; }
      else { integrator_type = INTGRTR_DUMMY; }
    }
    else if(*itr=="--constraint"){ 
      itr++;
      if(*itr=="none"){ constraint_type = CONST_NONE; }
      else if(*itr=="shake"){ constraint_type = CONST_SHAKE; }
      else if(*itr=="shake-settle"){ constraint_type = CONST_SHAKE_SETTLE; }
      else { constraint_type = CONST_DUMMY; }
    }
    else if(*itr=="--const-max-loops"){ constraint_max_loops = atoi((*++itr).c_str()); }
    else if(*itr=="--const-tolerance"){ constraint_tolerance = atof((*++itr).c_str()); }
    else if(*itr=="--thermo-const-max-loops"){ thermo_const_max_loops = atoi((*++itr).c_str()); }
    else if(*itr=="--thermo-const-tolerance"){ thermo_const_tolerance = atof((*++itr).c_str()); }
    else if(*itr=="--cutoff"){ cutoff=atof((*++itr).c_str()); }
    else if(*itr=="--n-steps"){ n_steps=atoi((*++itr).c_str()); }
    else if(*itr=="--time-step"){ time_step=atof((*++itr).c_str()); }

    else if(*itr=="--electrostatic"){
      itr++;
      if(*itr=="zero-dipole"){ electrostatic = ELCTRST_ZERODIPOLE; }
      else{ electrostatic = ELCTRST_DUMMY; }
    }
    else if(*itr=="--ele-alpha"){ ele_alpha = atof((*++itr).c_str()); }

    else if(*itr=="--thermostat"){
      itr++;
      if(*itr=="none"){ thermostat_type = THMSTT_NONE; }
      else if(*itr=="scaling"){ thermostat_type = THMSTT_SCALING; }
      else if(*itr=="hoover-evans"){ thermostat_type = THMSTT_HOOVER_EVANS; }
      else { thermostat_type = THMSTT_DUMMY; }
    }
    else if(*itr=="--expanded-ensemble"){
      itr++;
      if(*itr == "none"){ expanded_ensemble = EXPAND_NONE; }
      else if(*itr == "v-mcmd"){ expanded_ensemble = EXPAND_VMCMD; }
      else{ expanded_ensemble = EXPAND_DUMMY; } 
    }
    else if(*itr=="--temperature"){ temperature = atof((*++itr).c_str()); }
    else if(*itr=="--com-motion"){
      itr++;
      if(*itr=="none"){ com_motion = COM_NONE; }
      else if(*itr=="cancel"){ com_motion = COM_CANCEL; }
    }
    else if(*itr=="--com-cancel-group-name"){
      com_cancel_groups_name[n_com_cancel_groups] = ((*++itr).c_str());
      n_com_cancel_groups_name++;
    }
    else if(*itr=="--com-cancel-group-id"){
      com_cancel_groups[n_com_cancel_groups] = atoi((*++itr).c_str());
      n_com_cancel_groups++;
    }
    else if(*itr=="--random-seed"){ random_seed = atoi((*++itr).c_str()); }
    else if(*itr=="--box-division"){ 
      box_div[0] = atoi((*++itr).c_str());
      box_div[1] = atoi((*++itr).c_str());
      box_div[2] = atoi((*++itr).c_str());
    }
    else if(*itr=="--nsgrid-cutoff"){ nsgrid_cutoff = atof((*++itr).c_str()); }
    //else if(*itr=="--nsgrid-min-width"){ nsgrid_min_width= atof((*++itr).c_str()); }
    //else if(*itr=="--nsgrid-max-n-atoms"){ nsgrid_max_n_atoms = atof((*++itr).c_str()); }
    else if(*itr=="--nsgrid-update-intvl"){ nsgrid_update_intvl = atoi((*++itr).c_str()); }

    else if(*itr=="--print-interval-coord"){ print_intvl_crd = atoi((*++itr).c_str()); }
    else if(*itr=="--print-interval-velo"){ print_intvl_vel = atoi((*++itr).c_str()); }
    else if(*itr=="--print-interval-force"){ print_intvl_force = atoi((*++itr).c_str()); }
    else if(*itr=="--print-interval-log"){ print_intvl_log = atoi((*++itr).c_str()); }
    else if(*itr=="--print-interval-energy"){ print_intvl_energy = atoi((*++itr).c_str()); }
    else if(*itr=="--print-interval-energyflow"){ print_intvl_energyflow = atoi((*++itr).c_str()); }
    else if(*itr=="--print-interval-expand-lambda"){ print_intvl_expand_lambda = atoi((*++itr).c_str()); }
    else if(*itr=="--fn-o-restart"){ fn_o_restart = *++itr; }
    else if(*itr=="--fn-o-coord"){ fn_o_crd = *++itr; }
    else if(*itr=="--format-o-coord"){
      itr++;
      if(*itr == "gromacs"){ format_o_crd = CRDOUT_GROMACS; }
      else if(*itr == "presto"){  format_o_crd = CRDOUT_PRESTO; }
      else{ format_o_crd = CRDOUT_DUMMY; }
    }
    else if(*itr=="--fn-o-log"){ fn_o_log = *++itr; }
    else if(*itr=="--fn-o-energy"){ fn_o_energy = *++itr; }
    else if(*itr=="--fn-o-vmcmd-log"){ fn_o_vmcmd_log = *++itr; }
    else if(*itr=="--fn-o-expand-lambda"){ fn_o_expand_lambda = *++itr; }
    else if(*itr=="--format-o-expand-lambda"){
      itr++;
      if(*itr == "binary"){ format_o_expand_lambda = LAMBDAOUT_BIN; }
      else if(*itr == "ascii"){  format_o_expand_lambda = LAMBDAOUT_ASC; }
      else{ format_o_expand_lambda = LAMBDAOUT_DUMMY; }
    }

    else if(*itr=="--fn-o-energyflow"){ fn_o_energyflow = *++itr; }
    else if(*itr=="--dist-restraint"){
      itr++;
      if(*itr == "none"){ dist_restraint_type = DISTREST_NONE; }
      else if(*itr == "harmonic"){ dist_restraint_type = DISTREST_HARMONIC; }
      else{ dist_restraint_type = DISTREST_DUMMY; }
    }
    else if(*itr=="--dist-restraint-weight"){ dist_restraint_weight = atof((*++itr).c_str()); }
   else if(*itr=="--gpu-device-id"){ gpu_device_id = atof((*++itr).c_str()); }
    else{
      cerr<<"unknown keyword <"<<(*itr)<<">"<<endl;
    }
  }
}

