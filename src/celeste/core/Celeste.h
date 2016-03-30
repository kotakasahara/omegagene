#ifndef __CELESTE_H__
#define __CELESTE_H__

#include "CelesteObject.h"
#include "Config.h"
#include "DynamicsMode.h"
#include "MmSystem.h"
#include "Read.h"
#include "RunMode.h"
#include "Write.h"

class Celeste : private CelesteObject {
    Config config;

  public:
    Celeste() = default;
    int setup(int argn, char *argv[]);
    int main_stream();
    int test_mode();
    int dynamics_mode();
};

#endif
