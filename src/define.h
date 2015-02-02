#ifndef __DEFINE_H__
#define __DEFINE_H__

//#include <string>
//#include <ctime>
//#include <map>
//using namespace std;

//About this program.
//#define EXE "celeste"
//#define ABOUT_ME "celeste ver.0.01 25-nov-2013"
//#define PI 3.14159265
//#define MAGIC_NUMBER 66261

enum {
  M_TEST=0,
  M_DYNAMICS,
  M_DUMMY
};
enum {
  PRCS_SINGLE = 0,
  PRCS_MPI,
  PRCS_CUDA,
  PRCS_MPI_CUDA,
  PRCS_DUMMY
};
enum {
  INTGRTR_LEAPFLOG = 0,
  INTGRTR_DUMMY
};
enum {
  ELCTRST_ZERODIPOLE = 0,
  ELCTRST_ZEROMULTIPOLE,
  ELCTRST_DUMMY
};
enum {
  THMSTT_NONE = 0,
  THMSTT_BERENDSEN,
  THMSTT_HOOVER_EVANS,
  THMSTT_DUMMY
};

enum {
  COM_NONE = 0,
  COM_CANCEL
};


#endif
  
