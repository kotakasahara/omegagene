#ifndef __CELESTE_H__
#define __CELESTE_H__

#include <iostream>
#include <set>
#include <map>
#include <sstream>
#include <iomanip>
#include <cstdlib>

using namespace std;

#include "CelesteObject.h"
#include "Config.h"
#include "Read.h"
#include "Write.h"
#include "MmSystem.h"
#include "RunMode.h"
#include "DynamicsMode.h"

#if defined(F_CUDA) && defined(F_MPI) 
  #include "MpiGpuDynamicsMode.h"
#endif
#if defined(F_CUDA)
  #include "GpuDynamicsMode.h"
#endif


class Celeste : private CelesteObject{
 private:
  Config cfg;
 public:
  Celeste();
  int setup(int argn, char* argv[]);
  int main_stream();
  int test_mode();
  int dynamics_mode();
};

#endif
