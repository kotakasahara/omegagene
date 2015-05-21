#ifndef __EXPAND_H__
#define __EXPAND_H__

#include "CelesteObject.h"
using namespace std;

class Expand : public CelesteObject {
 private:

 protected:
  int write_lambda_interval;
 public:
  Expand();
  ~Expand();
  void set_lambda_interval(int in_lambda_interval);
};

#endif
