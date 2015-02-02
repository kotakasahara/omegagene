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
  processor = PRCS_SINGLE;
  integrator = INTGRTR_LEAPFLOG;
  cutoff = 12.0;
  n_steps = 1;
  time_step = 0.0005;
  electrostatic =  ELCTRST_ZERODIPOLE;
  ele_alpha = 0.0;
  thermostat = THMSTT_NONE;
  temperature = 300;
  center_of_motion = COM_NONE;

  box_div[0] = 1;
  box_div[1] = 1;
  box_div[2] = 1;

  print_intvl_crd = 10000;
  print_intvl_vel = 0;
  print_intvl_force = 0;
  fn_o_crd = "md_o.crd";
  fn_o_log = "md_o.log";
  fn_o_energy = "md_o.erg";
  fn_o_energyflow = "md_o.efl";
  nsgrid_cutoff = cutoff + 1.0;
  //nsgrid_min_width = cutoff * 0.5;
  //nsgrid_max_n_atoms = 100;
  nsgrid_update_intvl = 1;

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
      if(*itr=="leapflog"){ integrator = INTGRTR_LEAPFLOG; }
      else { integrator = INTGRTR_DUMMY; }
    }
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
      if(*itr=="none"){ thermostat = THMSTT_NONE; }
      else if(*itr=="berendsen"){ thermostat = THMSTT_BERENDSEN; }
      else if(*itr=="hoover_evans"){ thermostat = THMSTT_HOOVER_EVANS; }
      else { thermostat = THMSTT_DUMMY; }
    }
    else if(*itr=="--temperature"){ temperature = atof((*++itr).c_str()); }
    else if(*itr=="--center-of-motion"){
      itr++;
      if(*itr=="none"){ center_of_motion = COM_NONE; }
      else if(*itr=="cancel"){ center_of_motion = COM_CANCEL; }
      //itr++;
      //if(*itr=="all"){ center_of_motion = ; }
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
    else if(*itr=="--fn-o-coord"){ fn_o_crd = *++itr; }
    else if(*itr=="--fn-o-log"){ fn_o_log = *++itr; }
    else if(*itr=="--fn-o-energy"){ fn_o_energy = *++itr; }
    else if(*itr=="--fn-o-energyflow"){ fn_o_energyflow = *++itr; }
    else{
      cerr<<"unknown keyword <"<<(*itr)<<">"<<endl;
    }
  }
}

