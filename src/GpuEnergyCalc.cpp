#include "GpuEnergyCalc.h"

extern "C" int cuda_memcpy_htod_crd(real*& h_crd,
				    int n_atom_array);
extern "C" int cuda_set_atominfo(int n_atom_array);

extern "C" int cuda_pairwise_ljzd(const int offset_cellpairs, const int n_cal_cellpairs,
				  const int offset_cells,     const int n_cal_cells);

extern "C" int cuda_memcpy_dtoh_work(real_fc*& h_work, real_fc*& h_energy,
				     int n_atoms, int n_atom_array);

extern "C" int cuda_pair_sync();
extern "C" int cuda_thread_sync();
extern "C" int cuda_test(real*& h_work, real_pw*& h_energy,
			 int n_atoms);
extern "C" int cuda_zerodipole_constant(real_pw zcore,
					real_pw bcoeff,
					real_pw fcoeff);

GpuEnergyCalc::GpuEnergyCalc(MmSystem* in_mmsys)
  : EnergyCalc(in_mmsys){
}
GpuEnergyCalc::~GpuEnergyCalc(){
}

int GpuEnergyCalc::initial_preprocess(){
  EnergyCalc::initial_preprocess();
  cuda_zerodipole_constant(ff->ele->get_zcore(),
			   ff->ele->get_bcoeff(),
			   ff->ele->get_fcoeff());
  return 0;
}

int GpuEnergyCalc::calc_energy(bool grid_update){

  cout << "GpuEnergyCalc::calc_energy" << endl;
  const clock_t start_time_pair = clock();
  calc_energy_pairwise();
  const clock_t end_time_pair = clock();
  mmsys->ctime_calc_energy_pair += end_time_pair - start_time_pair;
  const clock_t start_time_bonded = clock();
  //cout << "GpuEnergyCalc::calc_energy_bonds" << endl;
  calc_energy_bonds();
  calc_energy_angles();
  calc_energy_torsions();
  calc_energy_impros();
  calc_energy_14nb();
  calc_energy_ele_excess();
  const clock_t end_time_bonded = clock();
  mmsys->ctime_calc_energy_bonded += end_time_bonded - start_time_bonded;
  //cout << "GpuEnergyCalc::finish_energy_pairwise" << endl;
  finish_energy_pairwise();
  mmsys->pote_ele += mmsys->energy_self_sum;

  if(grid_update){
    mmsys->nsgrid_update();
  }

  //cout << "end GpuEnergyCalc::calc_energy" << endl;
  return 0;
}

int GpuEnergyCalc::calc_energy_pairwise(){
  mmsys->nsgrid.init_energy_work();
  cuda_memcpy_htod_crd(mmsys->nsgrid.get_crd(),
		       mmsys->nsgrid.get_n_atom_array());
  cuda_set_atominfo(mmsys->nsgrid.get_n_atom_array());

  cuda_pairwise_ljzd(0, // offset_paridpairs,
		     mmsys->nsgrid.get_n_cell_pairs(), // n_cal_gridpairs,
		     0, // offset_grids,
		     mmsys->nsgrid.get_n_cells() ); // n_cal_grids);

  return 0;
}
int GpuEnergyCalc::finish_energy_pairwise(){
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
  cuda_memcpy_dtoh_work(mmsys->nsgrid.get_work(),
			mmsys->nsgrid.get_energy(),
			mmsys->n_atoms,
			mmsys->nsgrid.get_n_atom_array());

  //cout << "cuda_thread_sync();"<<endl;
  cuda_thread_sync();
  //cout << "copy forces;"<<endl;
  mmsys->nsgrid.get_ene_forces(mmsys->pote_vdw, mmsys->pote_ele, mmsys->force);
  //cout << "//copy forces;"<<endl;
  return 0;
}
