#ifndef __GENERAL_H__
#define __GENERAL_H__

#include <random>
#include <iostream>
using namespace std;

class RandomNum{
 private:
 protected:
  mt19937 mt;
  //static uniform_float dist_uni_f();
 public:
  RandomNum();
  ~RandomNum();

  int set_seed(int in_seed);
  float get_float_01();
};

#endif
