#include "ElectrostaticObject.h"

ElectrostaticObject::ElectrostaticObject()
  : CelesteObject(){

  if (DBG>=1)
    cout << "DBG1: ElectrostaticObject::ElectrostaticObject()" << endl;
  //mmsys = in_mmsys;
}
ElectrostaticObject::~ElectrostaticObject(){
}

int ElectrostaticObject::set_config_parameters(const Config* in_cfg){
  if (DBG>=1)
    cout << "DBG1: ElectrostaticObject::set_config_parameters()" << endl;
  return 0;
}

int ElectrostaticObject::initial_preprocess(){
  return 0;
}
