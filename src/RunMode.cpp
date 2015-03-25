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
  //enecal->set_config_parameters(cfg);
  
  return 0;
}

