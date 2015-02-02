#include "MpiGpuEnergyCalc.h"

extern "C" int cuda_memcpy_htod_crd(real_pw*& h_crd,
				    int*& h_grid_atom_index,
				    int n_atoms,
				    int n_grids);
extern "C" int cuda_set_atominfo(int n_atoms);
extern "C" int cuda_pairwise_ljzd(const int n_grid_pairs,     const int n_all_grids,
				  const int offset_gridpairs, const int n_cal_gridpairs,
				  const int offset_grids,     const int n_cal_grids);

extern "C" int cuda_memcpy_dtoh_work(real*& h_work, real_pw*& h_energy,
				     int n_atoms);

extern "C" int cuda_pair_sync();
extern "C" int cuda_thread_sync();
extern "C" int cuda_test(real*& h_work, real_pw*& h_energy,
			 int n_atoms);
extern "C" int cuda_zerodipole_constant(real_pw zcore,
					real_pw bcoeff,
					real_pw fcoeff);

MpiGpuEnergyCalc::MpiGpuEnergyCalc(MmSystem* in_mmsys)
  : GpuEnergyCalc(in_mmsys){
}
MpiGpuEnergyCalc::~MpiGpuEnergyCalc(){
  if (mpi_myrank == 0)  free_buffers();
}

int MpiGpuEnergyCalc::initial_preprocess(){
  GpuEnergyCalc::initial_preprocess();
  set_range_in_mpi_process();
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_n_ranks);  
  if (mpi_myrank == 0) alloc_buffers();
  return 0;
}
int MpiGpuEnergyCalc::alloc_buffers(){
  // Allocating member variables:
  //   real** buf_work
  // for buffering work calculated by
  // other mpi ranks.
  buf_work = new real[mmsys->n_atoms*3];
  for(int i=0; i < mmsys->n_atoms*3; i++) buf_work[i] = 0.0;
  return 0;
}
int MpiGpuEnergyCalc::free_buffers(){
  delete[] buf_work;
  return 0;
}

int MpiGpuEnergyCalc::set_range_in_mpi_process(){
  // set variables 
  // - offset_gridpairs;
  // -- ID of gridpairs for the first pair
  //    in this MPI process
  // - n_cal_gridpairs;
  // -- Number of pairs
  // - offset_grids;
  // - n_cal_grids;

  n_cal_gridpairs = (mmsys->nsgrid.n_grid_pairs + mpi_n_ranks -1) / mpi_n_ranks;
  offset_gridpairs = mpi_myrank * n_cal_gridpairs;
  if(offset_gridpairs + n_cal_gridpairs > mmsys->nsgrid.n_grid_pairs){
    n_cal_gridpairs = mmsys->nsgrid.n_grid_pairs - offset_gridpairs;
  }

  n_cal_grids = (mmsys->nsgrid.n_all_grids + mpi_n_ranks - 1) / mpi_n_ranks;
  offset_grids = mpi_myrank * n_cal_grids;
  if(offset_grids + n_cal_grids > mmsys->nsgrid.n_all_grids){
    n_cal_grids = mmsys->nsgrid.n_all_grids - offset_grids;
  }
  return 0;
}

int MpiGpuEnergyCalc::calc_energy(){
  int a =1;
  GpuEnergyCalc::calc_energy();
  /*
  //cout << "MpiGpuEnergyCalc::calc_energy" << endl;

  const clock_t start_time_pair = clock();
  calc_energy_pairwise();
  const clock_t end_time_pair = clock();
  mmsys->ctime_calc_energy_pair += end_time_pair - start_time_pair;
  const clock_t start_time_bonded = clock();

  //cout << "MpiGpuEnergyCalc::calc_energy_bonds" << endl;
  calc_energy_bonds();
  calc_energy_angles();
  calc_energy_torsions();
  calc_energy_impros();
  calc_energy_14nb();
  calc_energy_ele_excess();
  const clock_t end_time_bonded = clock();
  mmsys->ctime_calc_energy_bonded += end_time_bonded - start_time_bonded;
  //cout << "MpiGpuEnergyCalc::finish_energy_pairwise" << endl;

  finish_energy_pairwise();

  mmsys->pote_ele += mmsys->energy_self_sum;
  //cout << "end MpiGpuEnergyCalc::calc_energy" << endl;
  */
  return 0;
}


int MpiGpuEnergyCalc::calc_energy_pairwise(){
  cuda_memcpy_htod_crd(mmsys->nsgrid.crd,
		       mmsys->nsgrid.grid_atom_index,
		       mmsys->n_atoms,
		       mmsys->nsgrid.n_all_grids);
  cuda_set_atominfo(mmsys->n_atoms);
  // MPI
  // n_grids_cal_diag = (n_all_girds + ) / N_MPI_PROCESSES
  // ooset_grid_cal_diag = rank * n_grids_cal_diag

  cuda_pairwise_ljzd(mmsys->nsgrid.n_grid_pairs,
		     mmsys->nsgrid.n_all_grids,
		     offset_gridpairs,
		     n_cal_gridpairs,
		     offset_grids,
		     n_cal_grids);


  return 0;
}
int MpiGpuEnergyCalc::finish_energy_pairwise(){
  //cout << "cuda_thread_sync();"<<endl;
  cuda_pair_sync();
  
  //cout << "cuda_memcpy_dtoh_work();"<<endl;
  //cout << mmsys->nsgrid.work[1] << endl;
  //cout << mmsys->nsgrid.energy[1] << endl;
  //cout << mmsys->n_atoms << endl;
  //cout << "test"<<endl;
  //cuda_test(mmsys->nsgrid.work,
  //mmsys->nsgrid.energy,
  //mmsys->n_atoms);
  //cout << "test" << endl;
  cuda_memcpy_dtoh_work(mmsys->nsgrid.work,
			mmsys->nsgrid.energy,
			mmsys->n_atoms);

  //cout << "cuda_thread_sync();"<<endl;
  cuda_thread_sync();
  //cout << "copy forces;"<<endl;
  //cout << "copy energies;"<<endl;
  //EnergyCalc::calc_energy_pairwise();

  if (mpi_myrank == 0){
    for(int i=0; i < mmsys->n_atoms*3; i++) buf_work[i] = 0.0;    

    MPI_Reduce(&mmsys->nsgrid.work, buf_work,
	       mmsys->n_atoms*3,
	       MPI_DOUBLE,
	       MPI_SUM, 0, MPI_COMM_WORLD);
    real_pw buf_vdw = 0.0;
    MPI_Reduce(&mmsys->nsgrid.energy[0], &buf_vdw,
	       1,
	       MPI_FLOAT,
	       MPI_SUM, 0, MPI_COMM_WORLD);
    real_pw buf_ele = 0.0;
    MPI_Reduce(&mmsys->nsgrid.energy[1], &buf_ele,
	       1,
	       MPI_FLOAT,
	       MPI_SUM, 0, MPI_COMM_WORLD);
    
    mmsys->pote_vdw = buf_vdw;
    mmsys->pote_ele = buf_ele;
    for (int i = 0; i < mmsys->n_atoms; i++){
      mmsys->force[i][0] += buf_work[i*3];
      mmsys->force[i][1] += buf_work[i*3+1];
      mmsys->force[i][2] += buf_work[i*3+2];
    }
  }

  return 0;
}
