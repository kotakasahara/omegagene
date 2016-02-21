#ifndef __ELECTROSTATIC_OBJECT_H__
#define __ELECTROSTATIC_OBJECT_H__

#include "CelesteObject.h"
#include "Config.h"

class ElectrostaticObject : public CelesteObject {
 private:
 protected:
  //MmSystem* mmsys;
 public:
  ElectrostaticObject();
  virtual ~ElectrostaticObject();

  virtual int set_config_parameters(const Config* in_cfg);
  virtual int initial_preprocess();
};

#endif
