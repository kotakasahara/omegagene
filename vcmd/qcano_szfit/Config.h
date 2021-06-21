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
  string fname_o_qcano;
  Config();
  ~Config();

  void setAll(int argn,char* argv[]);
  void setAll(vector<string> arg);
  void operator=(Config op);

};

#endif
