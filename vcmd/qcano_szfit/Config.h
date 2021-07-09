#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <cstdlib>
#include <vector>
#include "define.h"

using namespace std;

class Config{
private:
public:
  int mode;
  string fname_i_cfg;
  string fname_i_params;
  string fname_i_qraw_is;
  string fname_o_qcano;
  string fname_o_qweight_opt;
  
  double target_error;

  double mc_temp;
  double mc_delta_x;
  double mc_delta_x_max;
  size_t mc_steps;
  int mc_log_interval;
  double mc_target_acc_ratio;
  int greedy_max_steps;
  size_t greedy_pivot;
  int mc_acc_duration;

  int qweight_write_mode;
  
  Config();
  ~Config();

  void setAll(int argn,char* argv[]);
  void setAll(vector<string> arg);
  void operator=(Config op);
};

#endif
