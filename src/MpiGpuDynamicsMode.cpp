#include "MpiGpuDynamicsMode.h"

extern "C" int cuda_alloc_atom_info(int n_atoms,
				    int n_grids,
				    int n_grid_pairs);
extern "C" int cuda_free_atom_info();
extern "C" int cuda_memcpy_htod_grid_pairs(GridPair*& h_grid_pairs,
					   int n_grid_pairs);
extern "C" int cuda_memcpy_htod_atom_info(real_pw*& h_charge_orig,
					  int*& h_atomtype_orig,
					  int n_atoms);
extern "C" int cuda_set_constant(int n_atoms,
				 int n_all_grids, int n_grid_pairs,
				 int max_n_atoms_grid,
				 real_pw cutoff, int n_atomtypes);

extern "C" int cuda_alloc_set_lj_params(real_pw* h_lj_6term,
					real_pw* h_lj_12term,
					int n_lj_types);
extern "C" int cuda_free_lj_params();

extern "C" int cuda_alloc_set_nb15off(int* h_nb15off1,
				      int* h_nb15off2,
				      int n_atoms);
extern "C" int cuda_free_nb15off();

extern "C" int cuda_memcpy_htod_atomids(int*& h_atomids,int n_atoms);
extern "C" int cuda_set_pbc(real_pw* l);
extern "C" int cuda_reset_work_ene(int n_atoms);
//extern "C" int cuda_lj_test(int x, int y);
//extern "C" int cuda_alloc_test(int*& test);
//extern "C" int cuda_free_test(int*& test);


MpiGpuDynamicsMode::MpiGpuDynamicsMode()
 : GpuDynamicsMode(){
} 

MpiGpuDynamicsMode::~MpiGpuDynamicsMode(){
} 

int MpiGpuDynamicsMode::initial_preprocess(){
  cout << "GpuDynamicsMode::initial_preprocess()"<<endl;
  MPI_Comm_size(MPI_COMM_WORLD, &n_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  GpuDynamicsMode::initial_preprocess();
  return 0;
}

int MpiGpuDynamicsMode::main_stream(){
  real_pw pbc[3] = {(real_pw)mmsys.pbc.L[0],
		    (real_pw)mmsys.pbc.L[1],
		    (real_pw)mmsys.pbc.L[2]};
  cuda_set_pbc(pbc);
  
  for(mmsys.cur_step = 0;
      mmsys.cur_step <= n_steps;
      mmsys.cur_step++){
    broadcast_each_step();
    calc_in_each_step();
  }
  
  mmsys.output_ctimes();

  return 0;
}

int MpiGpuDynamicsMode::broadcast_each_step(){
  MPI_Bcast(mmsys.nsgrid.atomids,
	    mmsys.n_atoms,
	    MPI_INT,
	    0,
	    MPI_COMM_WORLD);
  MPI_Bcast(mmsys.nsgrid.crd,
	    mmsys.n_atoms * 3,
	    mpi_real_pw,
	    0,
	    MPI_COMM_WORLD);
  MPI_Bcast(mmsys.nsgrid.grid_atom_index,
	    mmsys.nsgrid.n_all_grids+1,
	    MPI_INT,
	    0,
	    MPI_COMM_WORLD);
  return 0;
}
