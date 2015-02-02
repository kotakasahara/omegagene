#ifndef __WRITE_TRR_H__
#define __WRITE_TRR_H__

#include "Write.h"
#include "MmSystem.h"

class WriteTrr : public Write {
 private:
 public:
  int write_trr(MmSystem& mmsys,
		bool out_box,
		bool out_crd, bool out_vel, bool out_force);

};

#endif
