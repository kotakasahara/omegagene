#ifndef __WRITE_H__
#define __WRITE_H__

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <vector>
using namespace std;

#include "CelesteObject.h"



class Write : private CelesteObject {
 private:
  string filename;
  bool op;
 public:
  ofstream ofs;
  Write();
  Write(string inFn);
  void set_fn(string in_fn){filename = in_fn;};
  string getFn(){return filename;};
  bool is_open(){return op;};
  int open();
  int openApp();
  int close();
};

#endif
