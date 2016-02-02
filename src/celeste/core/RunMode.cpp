#include "RunMode.h"

RunMode::RunMode()
  : CelesteObject(){
}

RunMode::~RunMode(){
  if (DBG>=1)
    cout << "DBG1 RunMode::~RunMode()"<<endl;
}

int RunMode::initial_preprocess(){
  if(cfg->format_o_crd == CRDOUT_GROMACS){
    writer_trr = new WriteTrrGromacs();
  }else if (cfg->format_o_crd == CRDOUT_PRESTO){
    writer_trr = new WriteTrrPresto();
  }
  writer_trr->set_fn(cfg->fn_o_crd);
  writer_trr->open();

  return 0;
}
int RunMode::terminal_process(){
  writer_trr->close();
  delete writer_trr;
  delete mmsys.dist_restraint;
  return 0;
}

int RunMode::set_config_parameters(Config* in_cfg){
  cfg = in_cfg;

  if(cfg->dist_restraint_type==DISTREST_HARMONIC){
    mmsys.dist_restraint = new DistRestraintHarmonic();
  }else{
    mmsys.dist_restraint = new DistRestraintObject();
  }
  mmsys.dist_restraint->set_weight(cfg->dist_restraint_weight);

  if(cfg->pos_restraint_type==POSREST_HARMONIC){
    mmsys.pos_restraint = new PosRestraintHarmonic();
  }else{
    mmsys.pos_restraint = new PosRestraintObject();
  }
  mmsys.pos_restraint->set_weight(cfg->pos_restraint_weight);

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

