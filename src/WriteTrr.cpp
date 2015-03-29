#include "WriteTrr.h"

WriteTrr::WriteTrr()
  : Write(){
}

WriteTrr::~WriteTrr(){
}

int WriteTrr::write_trr(int n_atoms,
			int cur_step, real cur_time,
			real lx, real ly, real lz, 
			real** crd, real** vel_just, real_fc** force,
			bool out_box,
			bool out_crd, bool out_vel, bool out_force){


  int box_size = 0;
  if(out_box)
    box_size = 9 * sizeof(real);
  int x_size = 0;
  if(out_crd)
    x_size = n_atoms * 3 * sizeof(real);
  int v_size = 0;
  if(out_vel)
    v_size = n_atoms * 3 * sizeof(real);
  int f_size = 0;
  if(out_force)
    f_size = n_atoms * 3 * sizeof(real);
  
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
  ofs.write((const char*)&n_atoms, sizeof(int));
  ofs.write((const char*)&cur_step, sizeof(int));
  ofs.write((const char*)&dummy, sizeof(int));
  ofs.write((const char*)&cur_time, sizeof(real));
  ofs.write((const char*)&dummy_real, sizeof(real));

  if(out_box){
    real x = lx*0.1;
    real y = ly*0.1;
    real z = lz*0.1;

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
    for(int atomid=0; atomid < n_atoms; atomid++){
      real x = crd[atomid][0]*0.1;
      real y = crd[atomid][1]*0.1;
      real z = crd[atomid][2]*0.1;
      ofs.write((const char*)&x, sizeof(real));
      ofs.write((const char*)&y, sizeof(real));
      ofs.write((const char*)&z, sizeof(real));
    }
  }
  if(out_vel){
    for(int atomid=0; atomid < n_atoms; atomid++){
      real x = vel_just[atomid][0]*0.1;
      real y = vel_just[atomid][1]*0.1;
      real z = vel_just[atomid][2]*0.1;
      ofs.write((const char*)&x, sizeof(real));
      ofs.write((const char*)&y, sizeof(real));
      ofs.write((const char*)&z, sizeof(real));
    }
  }
  if(out_force){
    for(int atomid=0; atomid < n_atoms; atomid++){
      real x = force[atomid][0]*0.1;
      real y = force[atomid][1]*0.1;
      real z = force[atomid][2]*0.1;
      ofs.write((const char*)&x, sizeof(real));
      ofs.write((const char*)&y, sizeof(real));
      ofs.write((const char*)&z, sizeof(real));
    }
  }
  return 0;
}

WriteTTPVMcMDLog::WriteTTPVMcMDLog()
  : Write() {
}
WriteTTPVMcMDLog::~WriteTTPVMcMDLog(){
}


int WriteTTPVMcMDLog::write_ttpvMcMDLog(int step, int vstate){
  ofs << step << "\t" << vstate << endl;
  return 0;
}


WriteTableLog::WriteTableLog()
  : Write(){
}
WriteTableLog::~WriteTableLog(){
}

int WriteTableLog::write_header(){
  ofs.write((const char*)&MAGIC_NUMBER, sizeof(int));
  ofs.write((const char*)&n_cols, sizeof(int));
  return 0;
}
int WriteTableLog::write_row(real* values){
  ofs.write((const char*)values, sizeof(real)*n_cols);
  return 0;
}

WriteRestart::WriteRestart()
  : Write(){
}
WriteRestart::~WriteRestart(){
}
int WriteRestart::write_restart(int n_atoms, int n_steps,
				double time, double e_potential, double e_kinetic,
				real** crd, real** vel){
  int buf_int;
  double buf_dbl;
  open();
  // title
  buf_int = 80;
  char title[80];
  int i;
  strcpy(title, ABOUT_ME.c_str());
  
  ofs.write((const char*)&buf_int, sizeof(int));
  ofs.write(title, 80);
  ofs.write((const char*)&buf_int, sizeof(int));
  // #atoms, #velocities
  buf_int = 8;
  ofs.write((const char*)&buf_int, sizeof(int));
  ofs.write((const char*)&n_atoms, sizeof(int));
  ofs.write((const char*)&n_atoms, sizeof(int));
  ofs.write((const char*)&buf_int, sizeof(int));
  //
  buf_int = 36;
  ofs.write((const char*)&buf_int, sizeof(int));
  ofs.write((const char*)&n_steps, sizeof(int));
  ofs.write((const char*)&time, sizeof(double));
  buf_dbl = e_kinetic + e_potential;
  ofs.write((const char*)&buf_dbl, sizeof(double));
  ofs.write((const char*)&e_kinetic, sizeof(double));
  ofs.write((const char*)&e_potential, sizeof(double));
  ofs.write((const char*)&buf_int, sizeof(int));  

  buf_int = n_atoms*3*8;
  ofs.write((const char*)&buf_int, sizeof(int));
  for(int i=0; i < n_atoms; i++){
    for(int d=0; d < 3; d++){
      buf_dbl = (double)crd[i][d];
      ofs.write((const char*)&buf_dbl, sizeof(double));
    }
  }
  ofs.write((const char*)&buf_int, sizeof(int));  

  ofs.write((const char*)&buf_int, sizeof(int));
  for(int i=0; i < n_atoms; i++){
    for(int d=0; d < 3; d++){
      buf_dbl = (double)vel[i][d];
      ofs.write((const char*)&buf_dbl, sizeof(double));
    }
  }
  ofs.write((const char*)&buf_int,sizeof(int));  
  close();
  return 0;
}
