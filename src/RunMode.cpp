#include "RunMode.h"

RunMode::RunMode()
  : CelesteObject(){
}

RunMode::~RunMode(){
  if (DBG>=1)
    cout << "DBG1 RunMode::~RunMode()"<<endl;
}

int RunMode::initial_preprocess(){
  writer_trr.set_fn(cfg->fn_o_crd);
  writer_trr.open();
  return 0;
}
int RunMode::terminal_process(){
  writer_trr.close();
  return 0;
}

int RunMode::set_config_parameters(Config* in_cfg){
  cfg = in_cfg;
  if(DBG>=1)
    cout << "DBG1: RunMode::set_config_parameters()"<<endl;
  //#if defined(F_CUDA) && defined(F_MPI)
  //  enecal = new MpiGpuEnergyCalc(&mmsys);
  //#elif defined(F_CUDA)
  //  enecal = new GpuEnergyCalc(&mmsys);
  //#else
  //enecal = new EnergyCalc(&mmsys);
  //#endif
  /*
  switch(cfg->processor){
  case PRCS_SINGLE:
    enecal = new EnergyCalc(&mmsys);
    break;
  case PRCS_MPI:
  case PRCS_CUDA:
    enecal = new GpuEnergyCalc(&mmsys);    
  case PRCS_MPI_CUDA:
  default:
    enecal = new EnergyCalcObject(&mmsys);    
    break;
  }
  */
  enecal->set_config_parameters(cfg);

  n_steps = cfg->n_steps;
  
  print_intvl_crd = cfg->print_intvl_crd;
  print_intvl_vel = cfg->print_intvl_vel;
  print_intvl_log = cfg->print_intvl_log;
  print_intvl_energy = cfg->print_intvl_energy;
  print_intvl_energyflow = cfg->print_intvl_energyflow;

  fn_o_crd = cfg->fn_o_crd;
  fn_o_crd = cfg->fn_o_log;
  fn_o_energy = cfg->fn_o_energy;
  fn_o_energyflow = cfg->fn_o_energyflow;

  return 0;
}

