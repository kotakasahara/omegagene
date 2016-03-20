#ifndef __GENERAL_H__
#define __GENERAL_H__

#include <iostream>
#include <random>

class RandomNum {
  private:
  protected:
    std::mt19937 mt;
    // static uniform_float dist_uni_f();
  public:
    RandomNum();
    ~RandomNum();

    int set_seed(int in_seed);
    float get_float_01();
};

#endif
