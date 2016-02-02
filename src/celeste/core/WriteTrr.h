#ifndef __WRITE_TRR_H__
#define __WRITE_TRR_H__
#include <cstring>
#include "Write.h"
#include <vector>

class WriteTrr : public Write {
 private:
 public:
  WriteTrr();
  ~WriteTrr();
  virtual int write_trr(int n_atoms,
			int cur_step, real cur_time,
			real lx, real ly, real lz, 
			real** crd, real** vel_just, real_fc** force,
			real cpu_time, real total_e, real kinetic_e,
			real temperature, real potential_e,
			real vdw_e,
			bool out_box,
			bool out_crd, bool out_vel, bool out_force,
			int n_atoms_group, int* atom_group);
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
			real cpu_time, real total_e, real kinetic_e,
			real temperature, real potential_e,
			real vdw_e,
			bool out_box,
			bool out_crd, bool out_vel, bool out_force,
			int n_atoms_group, int* atom_group);
  
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
			real cpu_time, real total_e, real kinetic_e,
			real temperature, real potential_e,
			real vdw_e,
			bool out_box,
			bool out_crd, bool out_vel, bool out_force,
			int n_atoms_group, int* atom_group);
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

class WriteGroupCoord : public Write{
 private:
 protected:
 public:
  WriteGroupCoord();
  ~WriteGroupCoord();
  int write_aus_restart(const int aus_type,
			int n_enhance_groups,
			vector<int>  enhance_groups,
			int* n_atoms_in_groups,
			real*** crd_groups);
};

#endif
