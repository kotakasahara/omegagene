#ifndef __RUN_MODE_H__
#define __RUN_MODE_H__

#include "CelesteObject.h"
#include "MmSystem.h"
#include "Config.h"

#include "WriteTrr.h"

class RunMode : public CelesteObject {
 private:
 protected:
  Config* cfg;
  int integrator;
  int electrostatic;

  int com_cancel;
  int print_intvl_crd;
  int print_intvl_vel;
  int print_intvl_log;
  int print_intvl_energy;
  int print_intvl_energyflow;

  string fn_o_crd;
  string fn_o_log;
  string fn_o_energy;
  string fn_o_energyflow;

  WriteTrr* writer_trr;
  WriteRestart writer_restart;

 public:
  MmSystem mmsys;

  RunMode();
  virtual ~RunMode();
  virtual int initial_preprocess();
  virtual int terminal_process();
  virtual int set_config_parameters(Config* in_cfg);
};


#endif
