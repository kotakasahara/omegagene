#ifndef __ENERGY_CALC_OBJECT_H__
#define __ENERGY_CALC_OBJECT_H__

#include "CelesteObject.h"
#include "Config.h"
#include "SubBox.h"
#include "MmSystem.h"

class EnergyCalcObject : public CelesteObject{
 private:
 protected:
  MmSystem* mmsys;
  SubBox* subbox;
  real cutoff;
  
 public:
  EnergyCalcObject(MmSystem* in_mmsys,
		   SubBox* in_subbox);
  
  //void set_mmsystem(MmSystem& in_mmsys){ mmsys = in_mmsys; };
  virtual int set_config_parameters(const Config* in_cfg);
  virtual int initial_preprocess();
  virtual int calc_energy();
  virtual int calc_energy(bool grid_update);

};

#endif 

