#include "GpuDynamicsMode.h"

extern "C" int cuda_print_device_info(int myrank=0, bool verbose=false);
extern "C" int cuda_alloc_atom_info(int n_atoms,
				    int n_atom_array,
				    int max_n_cells,
				    int max_n_cell_pairs);
extern "C" int cuda_free_atom_info();
extern "C" int cuda_memcpy_htod_cell_pairs(CellPair*& h_cell_pairs,
					   int*& h_idx_head_cell_pairs,
					   int n_cell_pairs,
					   int n_cells);
extern "C" int cuda_memcpy_htod_atom_info(real_pw*& h_charge_orig,
					  int*& h_atomtype_orig,
					  int n_atoms);
extern "C" int cuda_set_cell_constant(int n_cells, int n_cell_pairs, int n_atom_array);
extern "C" int cuda_set_constant(int n_atoms, real cutoff, int n_atomtypes);
extern "C" int cuda_alloc_set_lj_params(real_pw* h_lj_6term,
					real_pw* h_lj_12term,
					int n_lj_types);
extern "C" int cuda_free_lj_params();

extern "C" int cuda_alloc_set_nb15off(int* h_nb15off1,
				      int* h_nb15off2,
				      int n_atoms);
extern "C" int cuda_memcpy_htod_atomids(int*& h_atomids,
					int n_atoms);
extern "C" int cuda_set_pbc(real_pw* l);
extern "C" int cuda_reset_work_ene(int n_atoms);
//extern "C" int cuda_lj_test(int x, int y);
//extern "C" int cuda_alloc_test(int*& test);
//extern "C" int cuda_free_test(int*& test);

GpuDynamicsMode::GpuDynamicsMode()
 : DynamicsMode(){
} 

GpuDynamicsMode::~GpuDynamicsMode(){
} 

int GpuDynamicsMode::update_device_cell_info(){
  cuda_set_cell_constant(mmsys.nsgrid.get_max_n_cells(),
			 mmsys.nsgrid.get_max_n_cell_pairs(),
			 mmsys.nsgrid.get_n_atom_array());
  cuda_memcpy_htod_cell_pairs(mmsys.nsgrid.get_cell_pairs(),
			      mmsys.nsgrid.get_idx_head_cell_pairs(),
			      mmsys.nsgrid.get_n_cell_pairs(),
			      mmsys.nsgrid.get_n_cells());
  cuda_memcpy_htod_atomids(mmsys.nsgrid.get_atomids(),
			   mmsys.nsgrid.get_max_n_atom_array());

  return 0;
}

int GpuDynamicsMode::initial_preprocess(){
  cout << "GpuDynamicsMode::initial_preprocess()"<<endl;

  DynamicsMode::initial_preprocess();

  //  mmsys.nsgrid.reorder_grid_pairs_random();
  cuda_print_device_info();
  cuda_alloc_atom_info(mmsys.n_atoms,
		       mmsys.nsgrid.get_max_n_atom_array(),
		       mmsys.nsgrid.get_max_n_cells(),
		       mmsys.nsgrid.get_max_n_cell_pairs());
  cuda_memcpy_htod_atom_info(mmsys.charge, mmsys.atom_type,
			     mmsys.n_atoms);
  //cuda_memcpy_htod_grid_pairs(mmsys.nsgrid.grid_pairs,
  //mmsys.nsgrid.n_grid_pairs);
  cuda_set_constant(mmsys.n_atoms,
		    (real_pw)cfg->cutoff, mmsys.n_lj_types);

  cuda_alloc_set_lj_params(mmsys.lj_6term, mmsys.lj_12term,
			   mmsys.n_lj_types);

  update_device_cell_info();  

  //real_pw pbc[3] = {(real_pw)mmsys.pbc.L[0],
  //(real_pw)mmsys.pbc.L[1],
  //(real_pw)mmsys.pbc.L[2]};
  cuda_set_pbc(mmsys.pbc.L);

  //cuda_alloc_set_nb15off(mmsys.nb15off1, mmsys.nb15off2,
  //mmsys.n_atoms);
  /* test
  cout << "param 0-0 " << mmsys.lj_6term[0] << " " << mmsys.lj_12term[0] << endl;
  cout << "param6 0-1 " << mmsys.lj_6term[1]<< " " << mmsys.lj_12term[1] << endl;
  cout << "param6 1-0 " << mmsys.lj_6term[mmsys.n_lj_types]<< " " << mmsys.lj_12term[mmsys.n_lj_types] << endl;
  cout << "param6 4-5 " << mmsys.lj_6term[4*mmsys.n_lj_types+5]<< " " << mmsys.lj_12term[4*mmsys.n_lj_types+5] << endl;
  cout << "param6 5-4 " << mmsys.lj_6term[5*mmsys.n_lj_types+4]<< " " << mmsys.lj_12term[5*mmsys.n_lj_types+4] << endl;
  cuda_lj_test(0,0);
  cuda_lj_test(0,1);
  cuda_lj_test(1,0);
  cuda_lj_test(4,5);
  cuda_lj_test(5,4);

  for (int i = 0 ; i < mmsys.n_lj_types; i++){
    cuda_lj_test(i,0);    
    cout << "param "<<i<<"-0 " << mmsys.lj_6term[i*mmsys.n_lj_types] << " " << mmsys.lj_12term[i*mmsys.n_lj_types] << endl;
  }
  */
  return 0;
}

int GpuDynamicsMode::terminal_process(){
  cout << "GpuDynamicsMode::terminal_process()"<<endl;
  DynamicsMode::terminal_process();
  cuda_free_atom_info();
  cuda_free_lj_params();
  return 0;
}

int GpuDynamicsMode::set_config_parameters(Config* in_cfg){
  DynamicsMode::set_config_parameters(in_cfg);
  return 0;
}

int GpuDynamicsMode::main_stream(){
  for(mmsys.cur_step = 0;
      mmsys.cur_step <= n_steps;
      mmsys.cur_step++){
    calc_in_each_step();
  }
  mmsys.output_ctimes();
  return 0;
}

int GpuDynamicsMode::calc_in_each_step(){
  const clock_t startTimeStep = clock();
  //cout << "(GPU) step " << mmsys.cur_step << endl;
  
  const clock_t startTimeReset = clock();
  cuda_reset_work_ene(mmsys.n_atoms);
  mmsys.cur_time = mmsys.cur_step * time_step;
  mmsys.reset_energy();
  const clock_t endTimeReset = clock();
  mmsys.ctime_cuda_reset_work_ene += endTimeReset - startTimeReset;
  
  const clock_t startTimeEne = clock();
  enecal->calc_energy(mmsys.cur_step%cfg->nsgrid_update_intvl==0);
  const clock_t endTimeEne = clock();
  mmsys.ctime_calc_energy += endTimeEne - startTimeEne;
  
  const clock_t startTimeVel = clock();
  update_velocities(1.0);
  const clock_t endTimeVel = clock();
  mmsys.ctime_update_velo += endTimeVel - startTimeVel;
  
  const clock_t startTimeKine = clock();
  cal_kinetic_energy(mmsys.vel_just);
  const clock_t endTimeKine = clock();
  mmsys.ctime_calc_kinetic += endTimeKine - startTimeKine;
  
  sub_output();
  sub_output_log();
  
  const clock_t startTimeCoord = clock();
  update_coordinates();
  const clock_t endTimeCoord = clock();
  mmsys.ctime_update_coord += endTimeCoord - startTimeCoord;
  revise_coordinates_pbc();

  if(mmsys.cur_step%cfg->nsgrid_update_intvl==0){
    const clock_t startTimeHtod = clock();
    //mmsys.nsgrid_update();
    mmsys.nsgrid_update_receiver();
    update_device_cell_info();
    const clock_t endTimeHtod = clock();
    mmsys.ctime_cuda_htod_atomids += endTimeHtod - startTimeHtod;
  }

  const clock_t endTimeStep = clock();
  mmsys.ctime_per_step += endTimeStep - startTimeStep;
  
  return 0;
}

