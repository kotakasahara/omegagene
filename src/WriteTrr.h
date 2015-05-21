#ifndef __WRITE_TRR_H__
#define __WRITE_TRR_H__
#include <cstring>
#include "Write.h"

class WriteTrr : public Write {
 private:
 public:
  WriteTrr();
  ~WriteTrr();
  virtual int write_trr(int n_atoms,
			int cur_step, real cur_time,
			real lx, real ly, real lz, 
			real** crd, real** vel_just, real_fc** force,
			bool out_box,
			bool out_crd, bool out_vel, bool out_force);
};

class WriteTrrGromacs : public WriteTrr {
 private:
 public:
  WriteTrrGromacs();
  ~WriteTrrGromacs();
  virtual int write_trr(int n_atoms,
			int cur_step, real cur_time,
			real lx, real ly, real lz, 
			real** crd, real** vel_just, real_fc** force,
			bool out_box,
			bool out_crd, bool out_vel, bool out_force);
  
};
class WriteTrrPresto : public WriteTrr {
 private:
 public:
  WriteTrrPresto();
  ~WriteTrrPresto();
  virtual int write_trr(int n_atoms,
			int cur_step, real cur_time,
			real lx, real ly, real lz, 
			real** crd, real** vel_just, real_fc** force,
			bool out_box,
			bool out_crd, bool out_vel, bool out_force);
  
};

class WriteRestart : public Write{
 private:
 protected:
 public:
  WriteRestart();
  ~WriteRestart();
  int write_restart(int n_atoms, int n_steps,
		    double time, double e_potential, double e_kinetic,
		    real** crd, real** vel);
};

#endif
