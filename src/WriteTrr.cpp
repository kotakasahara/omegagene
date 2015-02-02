#include "WriteTrr.h"

int WriteTrr::write_trr(MmSystem& mmsys,
			bool out_box,
			bool out_crd, bool out_vel, bool out_force){


  int box_size = 0;
  if(out_box)
    box_size = 9 * sizeof(real);
  int x_size = 0;
  if(out_crd)
    x_size = mmsys.n_atoms * 3 * sizeof(real);
  int v_size = 0;
  if(out_vel)
    v_size = mmsys.n_atoms * 3 * sizeof(real);
  int f_size = 0;
  if(out_force)
    f_size = mmsys.n_atoms * 3 * sizeof(real);
  
  int magic = 1993;
  int nchar1 = 13;
  int nchar2 = 12;
  real dummy_real = 0.0;
  int dummy = 0;

  if (!out_crd && !out_vel && !out_force) {
    return 1;
  }
  ofs.write((const char*)&magic, sizeof magic);
  ofs.write((const char*)&nchar1, sizeof(int));
  ofs.write((const char*)&nchar2, sizeof(int));
  ofs.write("GMX_trn_file", 12);
  ofs.write((const char*)&dummy, sizeof(int));
  ofs.write((const char*)&dummy, sizeof(int));
  ofs.write((const char*)&box_size, sizeof box_size);
  ofs.write((const char*)&dummy, sizeof(int));
  ofs.write((const char*)&dummy, sizeof(int));
  ofs.write((const char*)&dummy, sizeof(int));
  ofs.write((const char*)&dummy, sizeof(int));
  ofs.write((const char*)&x_size, sizeof x_size);
  ofs.write((const char*)&v_size, sizeof v_size);
  ofs.write((const char*)&f_size, sizeof f_size);
  ofs.write((const char*)&mmsys.n_atoms, sizeof mmsys.n_atoms);
  ofs.write((const char*)&mmsys.cur_step, sizeof mmsys.cur_step);
  ofs.write((const char*)&dummy, sizeof(int));
  ofs.write((const char*)&mmsys.cur_time, sizeof mmsys.cur_time);
  ofs.write((const char*)&dummy_real, sizeof(real));

  if(out_box){
    real x = mmsys.pbc.L[0]*0.1;
    real y = mmsys.pbc.L[1]*0.1;
    real z = mmsys.pbc.L[2]*0.1;

    ofs.write((const char*)&x, sizeof x);
    ofs.write((const char*)&dummy_real, sizeof(real));    
    ofs.write((const char*)&dummy_real, sizeof(real));    
    ofs.write((const char*)&dummy_real, sizeof(real));    
    ofs.write((const char*)&y, sizeof y);
    ofs.write((const char*)&dummy_real, sizeof(real));    
    ofs.write((const char*)&dummy_real, sizeof(real));    
    ofs.write((const char*)&dummy_real, sizeof(real));    
    ofs.write((const char*)&z, sizeof z);
  }
  if(out_crd){
    cout << "WRITETRR CRD" <<endl;
    for(int atomid=0; atomid < mmsys.n_atoms; atomid++){
      real x = mmsys.crd[atomid][0]*0.1;
      real y = mmsys.crd[atomid][1]*0.1;
      real z = mmsys.crd[atomid][2]*0.1;
      ofs.write((const char*)&x, sizeof(real));
      ofs.write((const char*)&y, sizeof(real));
      ofs.write((const char*)&z, sizeof(real));
    }
  }
  if(out_vel){
    for(int atomid=0; atomid < mmsys.n_atoms; atomid++){
      real x = mmsys.vel_just[atomid][0]*0.1;
      real y = mmsys.vel_just[atomid][1]*0.1;
      real z = mmsys.vel_just[atomid][2]*0.1;
      ofs.write((const char*)&x, sizeof(real));
      ofs.write((const char*)&y, sizeof(real));
      ofs.write((const char*)&z, sizeof(real));
    }
  }
  if(out_force){
    for(int atomid=0; atomid < mmsys.n_atoms; atomid++){
      real x = mmsys.force[atomid][0]*0.1;
      real y = mmsys.force[atomid][1]*0.1;
      real z = mmsys.force[atomid][2]*0.1;
      ofs.write((const char*)&x, sizeof(real));
      ofs.write((const char*)&y, sizeof(real));
      ofs.write((const char*)&z, sizeof(real));

    }
  }
  return 0;
}


