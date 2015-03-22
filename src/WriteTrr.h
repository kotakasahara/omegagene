#ifndef __WRITE_TRR_H__
#define __WRITE_TRR_H__

#include "Write.h"

class WriteTrr : public Write {
 private:
 public:
  WriteTrr();
  ~WriteTrr();
  int write_trr(int n_atoms,
		unsigned long cur_step, real cur_time,
		real lx, real ly, real lz, 
		real** crd, real** vel_just, real_fc** force,
		bool out_box,
		bool out_crd, bool out_vel, bool out_force);

};

class WriteTTPVMcMDLog : public Write {
 private:
 public:
  WriteTTPVMcMDLog();
  ~WriteTTPVMcMDLog();
  int write_ttpvMcMDLog(int step, int vstate);
};

class WriteTableLog : public Write {
 private:
 protected:
  int n_cols;
 public:
  WriteTableLog();
  ~WriteTableLog();  
  void set_ncolumns(int in_n_cols){ n_cols = in_n_cols; };
  int write_header();
  int write_row(real* values);
};

#endif
