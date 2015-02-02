#include "EnergyCalc.h"

EnergyCalcObject::EnergyCalcObject(MmSystem* in_mmsys,
				   SubBox* in_subbox)
  : CelesteObject(){
  mmsys = in_mmsys;
  subbox = in_subbox;
}

int EnergyCalcObject::set_config_parameters(const Config* cfg){
  cout << "EnergyCalcObject::set_config_parameters" << endl;
  cutoff = cfg->cutoff;
  return 0;
}
int EnergyCalcObject::initial_preprocess(){
  if (DBG >= 1)
    cout << "DBG1: EnergyCalcObject::initial_preprocess"<<endl;
  // ff->initial_preprocess();
  // ele->initial_preprocess();
  return 0;
}

int EnergyCalcObject::calc_energy(){
  return 0;
}

int EnergyCalcObject::calc_energy(bool grid_update){
  return 0;
}
