#include "SubBox.h"

#ifdef F_CUDA
extern "C" int cuda_set_device(int device_id);

extern "C" int cuda_alloc_atom_info(int max_n_atoms_exbox,
				    //int n_atom_array,
				    int max_n_cells,
				    int max_n_cell_pairs,
				    int n_columns);

extern "C" int cuda_free_atom_info();
extern "C" int cuda_memcpy_htod_atom_info(real_pw*& h_charge_orig,
					  int*& h_atomtype_orig);

extern "C" int cuda_set_cell_constant(const int  in_n_cells,
				      const int in_n_atoms_box,
				      const int  in_n_atom_array,
				      const int*  in_n_cells_xyz,
				      const int  in_n_columns,
				      const real_pw* in_l_cell_xyz,
				      const int* in_n_neighbor_xyz);

extern "C" int cuda_set_constant(
				 real_pw cutoff, 
				 real_pw cutoff_pairlist, int n_atomtypes);
extern "C" int cuda_alloc_set_lj_params(real_pw* h_lj_6term,
					real_pw* h_lj_12term,
					int n_lj_types,
					int* h_nb15off,
					const int max_n_nb15off,
					const int max_n_atoms,
					const int max_n_atom_array);

extern "C" int cuda_free_lj_params();

extern "C" int cuda_memcpy_htod_atomids(int*& h_atomids,
					int*& h_idx_xy_head_cell);
extern "C" int cuda_set_pbc(real_pw* l, real_pw* lb);
extern "C" int cuda_memcpy_htod_crd(real_pw*& h_crd);
				    
extern "C" int cuda_set_atominfo();
extern "C" int cuda_set_crd();
extern "C" int cuda_pairwise_ljzd(const bool flg_mod_15mask);

extern "C" int cuda_memcpy_dtoh_work(real_fc*& h_work, real_fc*& h_energy,
				     int n_atoms, int n_atom_array);

extern "C" int cuda_pair_sync();
extern "C" int cuda_thread_sync();
//extern "C" int cuda_test(real_pw*& h_work, real_pw*& h_energy,
//			 int n_atoms);
extern "C" int cuda_zerodipole_constant(real_pw zcore,
					real_pw bcoeff,
					real_pw fcoeff);

extern "C" int cuda_hostfree_atom_type_charge(int* h_atom_type,
					      real_pw* h_charge);


extern "C" int cuda_hostalloc_atom_type_charge(int*& h_atom_type,
					       real_pw*& h_charge,
					       const int n_atoms_system);

extern "C" int cuda_init_cellinfo(const int n_cells);
extern "C" int cuda_enumerate_cell_pairs(int*& h_atomids,
					 const int n_cells,// const int n_uni,
					 const int n_neighbor_col,
					 const int* idx_atom_cell_xy);
  #ifdef F_ECP
  extern "C" int cuda_memcpy_htod_cell_pairs(CellPair*& h_cell_pairs,
					   int*& h_idx_head_cell_pairs,
					     int n_cell_pairs);


  #endif

#endif

SubBox::SubBox(){
  ctime_setgrid = 0;
  ctime_enumerate_cellpairs = 0;
  ctime_calc_energy_pair = 0;
  ctime_calc_energy_bonded = 0;
  flg_constraint = 0;
  flg_settle = 0;
}
SubBox::~SubBox(){
#if defined(F_CUDA)
  cuda_free_atom_info();
  cuda_free_lj_params();
#endif
  free_variables();
  delete constraint;
  delete settle;

  delete thermostat;
}

int SubBox::alloc_variables(){
  //cout << "SubBox::alloc_variables"<<endl;
  crd = new real[max_n_atoms_exbox*3];
  crd_prev = new real[max_n_atoms_exbox*3];
  vel = new real[max_n_atoms_exbox*3];
  vel_next = new real[max_n_atoms_exbox*3];
  vel_just = new real[max_n_atoms_exbox*3];
  work = new real_fc[max_n_atoms_exbox*3];

  //
  if(cfg->thermostat_type!=THMSTT_NONE){
    buf_crd = new real[max_n_atoms_exbox*3];
  }

#if defined(F_CUDA)
  cout << "cuda_hostalloc_atom_type_charge " << max_n_atoms_exbox << endl;
  cuda_hostalloc_atom_type_charge(atom_type, charge, n_atoms);

#else
  charge = new real_pw[max_n_atoms_exbox];
 
  atom_type = new int[max_n_atoms_exbox];
#endif
  mass = new real_pw[max_n_atoms_exbox];
  mass_inv = new real_pw[max_n_atoms_exbox];
  //crd_box = new real*[n_boxes];
  //for(int i = 0; i < n_boxes; i++){
  //crd_box[i] = new real[max_n_atoms_box*3];
  //}
  atomids = new int[max_n_atoms_exbox];
  atomids_rev = new int[n_atoms];
  //atomids_box = new int*[n_boxes];
  //for(int i = 0; i < n_boxes; i++){
  //atomids_box[i] = new int[max_n_atoms_box];
  //}
  //crd = new real[n_atoms*3];  
  //atomids = new int[n_atoms];
  /*
  region_atoms = new int*[27];
  for(int i=0; i < 27; i++){
    region_atoms[i] = new int[max_n_atoms_region[i]];
  }
  n_region_atoms = new int[27];
  */

  /* for bonded potentials */
  

  if(rank==0)
    rank0_alloc_variables();

  return 0;
}
int SubBox::alloc_variables_for_bonds(int in_n_bonds){
  max_n_bonds = in_n_bonds;

  bond_epsiron = new real[max_n_bonds];
  bond_r0 = new real[max_n_bonds];
  bond_atomid_pairs = new int*[max_n_bonds];
  for(int i=0; i < max_n_bonds; i++){
    bond_atomid_pairs[i] = new int[2];
    bond_atomid_pairs[i][0] = 0;
    bond_atomid_pairs[i][1] = 0;
    bond_epsiron[i] = 0.0;
    bond_r0[i] = 0.0;
  }
  return 0;
}
int SubBox::alloc_variables_for_angles(int in_n_angles){
  max_n_angles = in_n_angles;
  angle_epsiron = new real[max_n_angles];
  angle_theta0 = new real[max_n_angles];
  angle_atomid_triads = new int*[max_n_angles];
  for(int i=0; i < max_n_angles; i++){
    angle_atomid_triads[i] = new int[3];
    angle_atomid_triads[i][0] = 0;
    angle_atomid_triads[i][1] = 0;
    angle_atomid_triads[i][2] = 0;
    angle_epsiron[i] = 0.0;
    angle_theta0[i] = 0.0;
  }
  return 0;
}
int SubBox::alloc_variables_for_torsions(int in_n_torsions){
  max_n_torsions = in_n_torsions;
  torsion_energy = new real[max_n_torsions];
  torsion_overlaps = new int[max_n_torsions];
  torsion_symmetry = new int[max_n_torsions];
  torsion_phase = new real[max_n_torsions];
  torsion_nb14 = new int[max_n_torsions];
  torsion_atomid_quads = new int*[max_n_torsions];
  for(int i=0; i<max_n_torsions; i++){
    torsion_atomid_quads[i] = new int[4];
    torsion_atomid_quads[i][0] = 0;
    torsion_atomid_quads[i][1] = 0;
    torsion_atomid_quads[i][2] = 0;
    torsion_atomid_quads[i][3] = 0;
    torsion_energy[i] = 0.0;
    torsion_overlaps[i] = 0;
    torsion_symmetry[i] = 0;
    torsion_phase[i] = 0.0;
    torsion_nb14[i] = 0;
  }
  return 0;
}
int SubBox::alloc_variables_for_impros(int in_n_impros){
  max_n_impros = in_n_impros;
  impro_energy = new real[max_n_impros];
  impro_overlaps = new int[max_n_impros];
  impro_symmetry = new int[max_n_impros];
  impro_phase = new real[max_n_impros];
  impro_nb14 = new int[max_n_impros];
  impro_atomid_quads = new int*[max_n_impros];
  for(int i=0; i<max_n_impros; i++){
    impro_atomid_quads[i] = new int[4];
    impro_atomid_quads[i][0] = 0;
    impro_atomid_quads[i][1] = 0;
    impro_atomid_quads[i][2] = 0;
    impro_atomid_quads[i][3] = 0;
    impro_energy[i] = 0.0;
    impro_overlaps[i] = 0;
    impro_symmetry[i] = 0;
    impro_phase[i] = 0.0;
    impro_nb14[i] = 0;
  }
  return 0;
}
int SubBox::alloc_variables_for_nb14(int in_n_nb14){
  max_n_nb14 = in_n_nb14;
  nb14_coeff_vdw = new real[max_n_nb14];
  nb14_coeff_ele = new real[max_n_nb14];
  nb14_atomid_pairs = new int*[max_n_nb14];
  nb14_atomtype_pairs = new int*[max_n_nb14];
  for(int i=0; i<max_n_nb14; i++){
    nb14_atomid_pairs[i] = new int[2];
    nb14_atomid_pairs[i][0] = 0;
    nb14_atomid_pairs[i][1] = 0;
    nb14_atomtype_pairs[i] = new int[2];
    nb14_atomtype_pairs[i][0] = 0;
    nb14_atomtype_pairs[i][1] = 0;
    nb14_coeff_vdw[i] = 0.0;
    nb14_coeff_ele[i] = 0.0;
  }
  return 0;
}
int SubBox::alloc_variables_for_excess(int in_n_excess){
  max_n_excess = in_n_excess;
  excess_pairs = new int*[max_n_excess];
  for(int i=0; i<max_n_excess; i++){
    excess_pairs[i] = new int[2];
    excess_pairs[i][0] = -1;
    excess_pairs[i][1] = -1;
  }
  return 0;
}
int SubBox::alloc_variables_for_nb15off(int in_max_n_nb15off){
  max_n_nb15off = in_max_n_nb15off;
  nb15off = new int[max_n_atoms_exbox * max_n_nb15off];
  n_nb15off = 0;
  for(int i=0; i < max_n_atoms_exbox; i++){
    for(int j=0; j < max_n_nb15off; j++){
      nb15off[i*max_n_nb15off + j] = -1;
    }
  }
  return 0;
}

int SubBox::init_variables(){
  //for(int i=0; i < 27; i++){
  //n_region_atoms[i]  = 0;
  //for(int j=0; j < max_n_atoms_region[27]; j++){
  //region_atoms[i][j] = -1;
  //}
  //} 
  init_energy();
  init_work();
  return 0;
}
int SubBox::free_variables(){
  delete[] crd;
  delete[] crd_prev;
  delete[] vel;
  delete[] vel_next;
  delete[] vel_just;
  delete[] work;
  //delete[] frc;
  delete[] atomids;

  if(cfg->thermostat_type!=THMSTT_NONE){
    delete[] buf_crd;
  }

#if defined(F_CUDA)
  cuda_hostfree_atom_type_charge(atom_type, charge);
#else
  delete[] charge;
  delete[] atom_type;
#endif
  delete[] mass;
  delete[] mass_inv;
  free_variables_for_bonds();
  free_variables_for_angles();
  free_variables_for_torsions();
  free_variables_for_impros();
  free_variables_for_nb14();
  free_variables_for_excess();
  free_variables_for_nb15off();

  //for(int i=0; i < n_boxes; i++){
  //delete[] crd_box[i];
  //}
  //delete[] crd_box;

  //for(int i=0; i < n_boxes; i++){
  //delete[] atomids_box[i];
  //}
  //delete[] atomids_box;

  //delete[] crd;
  //delete[] atomids;

  //for(int i=0; i < 27; i++){
  //delete[] region_atoms[i];
  //}
  //delete[] region_atoms;
  //delete[] n_region_atoms;

  delete[] bp_bonds;
  delete[] bp_angles;
  delete[] bp_torsions;
  delete[] bp_impros;
  delete[] bp_nb14;

  if(rank==0)
    rank0_free_variables();

  return 0;
}
int SubBox::free_variables_for_bonds(){
  for(int i=0; i < max_n_bonds; i++){
    delete[] bond_atomid_pairs[i];
  }
  delete[] bond_atomid_pairs;
  delete[] bond_epsiron;
  delete[] bond_r0;
  return 0;
}
int SubBox::free_variables_for_angles(){
  for(int i=0; i < max_n_angles; i++){
    delete[] angle_atomid_triads[i];
  }
  delete[] angle_atomid_triads;
  delete[] angle_epsiron;
  delete[] angle_theta0;
  return 0;
}
int SubBox::free_variables_for_torsions(){
  for(int i=0; i < max_n_torsions; i++){
    delete[] torsion_atomid_quads[i];
  }
  delete[] torsion_atomid_quads;
  delete[] torsion_energy;
  delete[] torsion_overlaps;
  delete[] torsion_symmetry;
  delete[] torsion_phase;
  delete[] torsion_nb14;
  return 0;
}
int SubBox::free_variables_for_impros(){
  for(int i=0; i < max_n_impros; i++){
    delete[] impro_atomid_quads[i];
  }
  delete[] impro_atomid_quads;
  delete[] impro_energy;
  delete[] impro_overlaps;
  delete[] impro_symmetry;
  delete[] impro_phase;
  delete[] impro_nb14;
  return 0;
}
int SubBox::free_variables_for_nb14(){
  for(int i=0; i < max_n_nb14; i++){
    delete[] nb14_atomid_pairs[i];
    delete[] nb14_atomtype_pairs[i];
  }
  delete[] nb14_atomid_pairs;
  delete[] nb14_atomtype_pairs;
  delete[] nb14_coeff_vdw;
  delete[] nb14_coeff_ele;
  return 0;

}

int SubBox::free_variables_for_excess(){
  for(int i=0; i < max_n_excess; i++){
    delete[] excess_pairs[i];
  }
  delete[] excess_pairs;
  return 0;
}

int SubBox::free_variables_for_nb15off(){
  delete[] nb15off;
  return 0;
}

int SubBox::set_parameters(int in_n_atoms, PBC* in_pbc, 
			   Config* in_cfg,
			   real in_cutoff_pair,
			   int in_n_boxes_x, int in_n_boxes_y, int in_n_boxes_z){
  // set
  //   n_atoms, pbc, cutoff_pair
  //   n_boxes_xyz,  n_boxes
  //   max_n_atoms_box
  //   box_l, exbox_l
  //   box_crd <- set_box_crd()
  //   box_lower, box_upper, exbox_lower, exbox_upper
  //   

  //cout << "SubBox::set_parameters" <<endl;
  n_atoms = in_n_atoms;
  pbc = in_pbc;
  cfg = in_cfg;
  time_step = cfg->time_step;
  time_step_inv = 1.0 / time_step;
  time_step_inv_sq = time_step_inv * time_step_inv;
  //temperature_coeff = 1.0 / (GAS_CONST * (real)d_free) * JOULE_CAL * 1e3 * 2.0;
  //temperature_coef = (2.0 * JOULE_CAL * 1e+3) / (GAS_CONST * (real)d_free);
  //cout << "DBG coef1: " << temperature_coef << " " << d_free << endl;

  cutoff_pair = in_cutoff_pair;
  n_boxes_xyz[0] = in_n_boxes_x;
  n_boxes_xyz[1] = in_n_boxes_y;
  n_boxes_xyz[2] = in_n_boxes_z;
  n_boxes = n_boxes_xyz[0] * n_boxes_xyz[1] * n_boxes_xyz[2];

  max_n_atoms_box = n_atoms / n_boxes * COEF_MAX_N_ATOMS_BOX;

  // temporary:
  max_n_atoms_exbox = n_atoms / n_boxes * COEF_MAX_N_ATOMS_BOX;
  // 

  rank = 0;

#if defined(F_MPI)
  //rank = MPI
#endif
  //cout << "rank: " << rank << " max_n_atoms_box: " << max_n_atoms_box << endl;
  for(int d=0; d<3; d++){
    box_l[d] = pbc->L[d] / n_boxes_xyz[d];  
    exbox_l[d] = box_l[d] + cutoff_pair;
  }  
  set_box_crd();
  for(int d=0; d<3; d++){
    box_lower[d] = pbc->lower_bound[d] + box_crd[d] * box_l[d];
    box_upper[d] = box_lower[d] + box_l[d];
    exbox_lower[d] = box_lower[d] - cutoff_pair_half;
    exbox_upper[d] = box_upper[d] + cutoff_pair_half;
  }
  //cout << "DBG SubBox box_lower " << box_lower[0] << " " << box_lower[1] << " " << box_lower[2] << endl;
  //cout << "DBG SubBox box_upper " << box_upper[0] << " " << box_upper[1] << " " << box_upper[2] << endl;
  ff = ForceField();
  ff.set_config_parameters(cfg);
  ff.initial_preprocess((const PBC*)pbc);

  return 0;
}

int SubBox::set_nsgrid(){
  // #ifdef F_CUDA  
  //  cuda_print_device_info(0, true);
  //#endif

  //cout << "set_grid_parameters" << endl;
  nsgrid.set_grid_parameters(n_atoms, cfg->nsgrid_cutoff, 
			     pbc, max_n_nb15off, nb15off);
  //cout << "set_box_info" << endl;
  nsgrid.set_box_info(n_boxes_xyz, box_l);

  nsgrid.set_max_n_atoms_region();

  //nsgrid.setup_crd_into_grid(crd, charge, atom_type);
  nsgrid.set_grid_xy();
  nsgrid.alloc_variables();
  revise_coordinates_pbc();
  nsgrid.set_crds_to_homebox(get_crds(),
			     get_atomids(),
			     get_n_atoms_box());


  //nsgrid.setup_replica_regions();
  //nsgrid.alloc_variables_box();

  nsgrid.set_atoms_into_grid_xy();
  nsgrid.set_atomids_buf();

#ifdef F_CUDA
  if(cfg->gpu_device_id >= 0)
    cuda_set_device(cfg->gpu_device_id);
  cout << "gpu_device_setup()"<<endl;
  gpu_device_setup();

  cuda_memcpy_htod_atom_info(charge, atom_type);

  update_device_cell_info();  

  nsgrid_crd_to_gpu();
  

  #ifdef F_ECP
    nsgrid.enumerate_cell_pairs();
    cuda_memcpy_htod_cell_pairs(nsgrid.get_cell_pairs(),
				nsgrid.get_idx_head_cell_pairs(),
				nsgrid.get_n_cell_pairs());

  #else

    cuda_enumerate_cell_pairs(nsgrid.get_atomids(),
			      nsgrid.get_n_cells(),
			      //nsgrid.get_n_uni(),
			      nsgrid.get_n_neighbor_cols(),
			      nsgrid.get_idx_atom_cell_xy());

  #endif

#else
    nsgrid.enumerate_cell_pairs();
#endif

  return 0;
}
int SubBox::nsgrid_crd_to_gpu(){
#ifdef F_CUDA
  //nsgrid.init_energy_work();
  
  cuda_memcpy_htod_crd(nsgrid.get_crd());
  cuda_set_crd();

  //nsgrid.update_crd((const real**)crd);
#endif
  return 0;
}

int SubBox::nsgrid_update(){
  const clock_t startTimeSet = clock();
  //nsgrid.init_variables_box();
  nsgrid.set_crds_to_homebox(crd,
			     atomids,
			     n_atoms_box);

  nsgrid.set_atoms_into_grid_xy();
  const clock_t endTimeSet = clock();
  //nsgrid.enumerate_cell_pairs();

#if defined(F_CUDA)  
  update_device_cell_info();  
  nsgrid_crd_to_gpu();
  #ifdef F_ECP
    nsgrid.enumerate_cell_pairs();
    cuda_memcpy_htod_cell_pairs(nsgrid.get_cell_pairs(),
				nsgrid.get_idx_head_cell_pairs(),
				nsgrid.get_n_cell_pairs());
			
  #else
  cuda_enumerate_cell_pairs(nsgrid.get_atomids(),
			    nsgrid.get_n_cells(),
			    //nsgrid.get_n_uni(),
			    nsgrid.get_n_neighbor_cols(),
			    nsgrid.get_idx_atom_cell_xy());

  #endif

#else
  nsgrid.enumerate_cell_pairs();
#endif

  const clock_t endTimePair = clock();
  ctime_setgrid += endTimeSet - startTimeSet;
  ctime_enumerate_cellpairs += endTimePair - endTimeSet;

  flg_mod_15mask = true;
  return 0;
}

int SubBox::nsgrid_update_receiver(){
  nsgrid.init_energy_work();
  nsgrid.set_atomids_buf();
  return 0;
}

int SubBox::rank0_alloc_variables(){
  // alloc
  all_atomids = new int*[n_boxes];
  for(int i=0; i<n_boxes; i++){
    all_atomids[i] = new int[max_n_atoms_exbox];
  }
  all_n_atoms = new int[n_boxes];
  for(int i=0; i<n_boxes; i++){
    all_n_atoms[i] = 0;
  }
  return 0;
}
int SubBox::rank0_free_variables(){
  //free
  for(int i=0; i<n_boxes; i++){
    delete[] all_atomids[i];
  }
  delete[] all_atomids;

  delete[] all_n_atoms;
  return 0;
}

int SubBox::initial_division(real** in_crd,
			     real** in_vel,
			     real_pw* in_charge,
			     real_pw* in_mass,
			     int* in_atom_type){
  if(rank==0){
    rank0_div_box(in_crd, in_vel);
    rank0_send_init_data(in_crd, in_vel, in_charge, in_mass, in_atom_type);
  }else{
    recv_init_data();
  }
  return 0;
}

int SubBox::rank0_div_box(real** in_crd,
			  real** in_vel){
  // called only the begining of the simulation
  // divide
  // set
  //   box_assign
  //   all_atomids
  //   atomids_rev
  //   all_n_atoms
  //     
  
  for(int atomid = 0; atomid < n_atoms; atomid++){
    int box_assign[3];
    for(int d=0; d<3; d++)
      box_assign[d] = floor((in_crd[atomid][d] - pbc->lower_bound[d]) / box_l[d]);
    int box_id = get_box_id_from_crd(box_assign);
    all_atomids[box_id][all_n_atoms[box_id]] = atomid;
    atomids_rev[atomid] = all_n_atoms[box_id];
    all_n_atoms[box_id]++;
  }

  return 0;
}

int SubBox::rank0_send_init_data(real** in_crd,
				 real** in_vel,
				 real_pw* in_charge,
				 real_pw* in_mass,
				 int* in_atom_type){
  // set
  //   crd
  //   vel_next
  //   atomids_rev

  for(int i_box=1; i_box < n_boxes; i_box++){
    // MPI_SEND  => rank=i,
    // MPI_SEND  => rank=i, all_crd[i], all_n_atoms[i]*3, real
    // MPI_SEND  => rank=i, all_atomids[i], all_n_atoms[i]3, int

    for(int i_atom = 0; i_atom < n_atoms; i_atom++){
      atomids_rev[i_atom] = -1;
    }
    // velocities
    for(int i_atom = 0, i_atom3=0;
	i_atom < all_n_atoms[i_box]; 
	i_atom++, i_atom3+=3){
      crd[i_atom3] = in_crd[all_atomids[i_box][i_atom]][0];
      crd[i_atom3+1] = in_crd[all_atomids[i_box][i_atom]][1];
      crd[i_atom3+2] = in_crd[all_atomids[i_box][i_atom]][2];
      vel_next[i_atom3] = in_vel[all_atomids[i_box][i_atom]][0];
      vel_next[i_atom3+1] = in_vel[all_atomids[i_box][i_atom]][1];
      vel_next[i_atom3+2] = in_vel[all_atomids[i_box][i_atom]][2];
      charge[i_atom] = in_charge[all_atomids[i_box][i_atom]];
      atom_type[i_atom] = in_atom_type[all_atomids[i_box][i_atom]];
      atomids_rev[all_atomids[i_box][i_atom]] = i_atom;
      //Atomids 
      // write the MPI code here!!
    }
  }
  for(int i_atom = 0; i_atom < n_atoms; i_atom++){
    atomids_rev[i_atom] = -1;
  }
  for(int i_atom = 0, i_atom3=0;
      i_atom < all_n_atoms[0]; 
      i_atom++, i_atom3+=3){
    crd[i_atom3] = in_crd[all_atomids[0][i_atom]][0];
    crd[i_atom3+1] = in_crd[all_atomids[0][i_atom]][1];
    crd[i_atom3+2] = in_crd[all_atomids[0][i_atom]][2];
    vel_next[i_atom3] = in_vel[all_atomids[0][i_atom]][0];
    vel_next[i_atom3+1] = in_vel[all_atomids[0][i_atom]][1];
    vel_next[i_atom3+2] = in_vel[all_atomids[0][i_atom]][2];
    charge[i_atom] = in_charge[all_atomids[0][i_atom]];
    mass[i_atom] = in_mass[all_atomids[0][i_atom]];
    mass_inv[i_atom] = 1.0 / mass[i_atom];
    atom_type[i_atom] = in_atom_type[all_atomids[0][i_atom]];
    atomids_rev[all_atomids[0][i_atom]] = i_atom;
  }
  //memcpy(crd,     all_crd[0],     all_n_atoms[rank]*3*sizeof(real));
  memcpy(atomids, all_atomids[0], all_n_atoms[0]*sizeof(int));
  n_atoms_box = all_n_atoms[0];

  //temporary
  n_atoms_exbox = all_n_atoms[0];

  cpy_crd_prev();

  return 0;
}
int SubBox::recv_init_data(){
  return 0;
}
int SubBox::set_bond_potentials(const int** in_bond_atomid_pairs,
				const real* in_bond_epsiron,
				const real* in_bond_r0){
  n_bonds = 0;
  for(int i=0; i < max_n_bonds; i++){
    int aid1_b = atomids_rev[in_bond_atomid_pairs[i][0]];
    int aid2_b = atomids_rev[in_bond_atomid_pairs[i][1]];
    if(aid1_b != -1 && aid2_b != -1){
      real mid_crd[3];
      for(int d=0; d < 3; d++){
	mid_crd[d] = (crd[aid1_b*3+d]
		      + crd[aid2_b*3+d])*0.5;
      }

      if(is_in_box(mid_crd)){
	bond_atomid_pairs[n_bonds][0] = aid1_b;
	bond_atomid_pairs[n_bonds][1] = aid2_b;
	bond_epsiron[n_bonds] = in_bond_epsiron[i];
	bond_r0[n_bonds] = in_bond_r0[i];
	n_bonds++;
      }
    }
  }
  return 0;
}
int SubBox::set_angle_potentials(const int** in_angle_atomid_triads,
				 const real* in_angle_epsiron,
				 const real* in_angle_theta0){
  n_angles = 0;
  for(int i=0; i < max_n_angles; i++){
    int aid1_b = atomids_rev[in_angle_atomid_triads[i][0]];
    int aid2_b = atomids_rev[in_angle_atomid_triads[i][1]];
    int aid3_b = atomids_rev[in_angle_atomid_triads[i][2]];
    if( aid1_b != -1 && aid2_b != -1 && aid3_b != -1){
      real mid_crd[3];
      for(int d=0; d < 3; d++){
	mid_crd[d] = (crd[aid1_b*3+d]
		      + crd[aid2_b*3+d]
		      + crd[aid3_b*3+d]) /3;
      }
      if(is_in_box(mid_crd)){
	angle_atomid_triads[n_angles][0] = aid1_b;
	angle_atomid_triads[n_angles][1] = aid2_b;
	angle_atomid_triads[n_angles][2] = aid3_b;
	angle_epsiron[n_angles] = in_angle_epsiron[i];
	angle_theta0[n_angles] = in_angle_theta0[i];
	n_angles++;
      }
    }
  }
  return 0;
}
int SubBox::set_torsion_potentials(const int** in_torsion_atomid_quads,
				   const real* in_torsion_energy,
				   const int* in_torsion_overlaps,
				   const int* in_torsion_symmetry,
				   const real* in_torsion_phase,
				   const int* in_torsion_nb14){
  n_torsions = 0;
  for(int i=0; i < max_n_torsions; i++){
    int aid1_b = atomids_rev[in_torsion_atomid_quads[i][0]];
    int aid2_b = atomids_rev[in_torsion_atomid_quads[i][1]];
    int aid3_b = atomids_rev[in_torsion_atomid_quads[i][2]];
    int aid4_b = atomids_rev[in_torsion_atomid_quads[i][3]];
    if(aid1_b != -1 && aid2_b != -1 && aid3_b != -1 && aid4_b != -1){
      real mid_crd[3];
      for(int d=0; d < 3; d++){
	mid_crd[d] = (crd[aid1_b*3+d]
		      + crd[aid2_b*3+d]
		      + crd[aid3_b*3+d]
		      + crd[aid4_b*3+d]) *0.25;
      }
      if(is_in_box(mid_crd)){
	torsion_atomid_quads[n_torsions][0] = aid1_b;
	torsion_atomid_quads[n_torsions][1] = aid2_b;
	torsion_atomid_quads[n_torsions][2] = aid3_b;
	torsion_atomid_quads[n_torsions][3] = aid4_b;
	torsion_energy[n_torsions] = in_torsion_energy[i];
	torsion_overlaps[n_torsions] = in_torsion_overlaps[i];
	torsion_symmetry[n_torsions] = in_torsion_symmetry[i];
	torsion_phase[n_torsions] = in_torsion_phase[i];
	torsion_nb14[n_torsions] = in_torsion_nb14[i];
	n_torsions++;
      }
    }
  }

  return 0;
}

int SubBox::set_impro_potentials(const int** in_impro_atomid_quads,
				 const real* in_impro_energy,
				 const int* in_impro_overlaps,
				 const int* in_impro_symmetry,
				 const real* in_impro_phase,
				 const int* in_impro_nb14){
  n_impros=0;
  for(int i=0; i < max_n_impros; i++){
    int aid1_b = atomids_rev[in_impro_atomid_quads[i][0]];
    int aid2_b = atomids_rev[in_impro_atomid_quads[i][1]];
    int aid3_b = atomids_rev[in_impro_atomid_quads[i][2]];
    int aid4_b = atomids_rev[in_impro_atomid_quads[i][3]];
    if(aid1_b != -1 && aid2_b != -1 && aid3_b != -1 && aid4_b != -1){
      real mid_crd[3];
      for(int d=0; d < 3; d++){
	mid_crd[d] = (crd[aid1_b*3+d]
		      + crd[aid2_b*3+d]
		      + crd[aid3_b*3+d]
		      + crd[aid4_b*3+d]) *0.25;
      }
      if(is_in_box(mid_crd)){
	impro_atomid_quads[n_impros][0] = aid1_b;
	impro_atomid_quads[n_impros][1] = aid2_b;
	impro_atomid_quads[n_impros][2] = aid3_b;
	impro_atomid_quads[n_impros][3] = aid4_b;
	impro_energy[n_impros] = in_impro_energy[i];
	impro_overlaps[n_impros] = in_impro_overlaps[i];
	impro_symmetry[n_impros] = in_impro_symmetry[i];
	impro_phase[n_impros] = in_impro_phase[i];
	impro_nb14[n_impros] = in_impro_nb14[i];

	n_impros++;
      }
    }
  }
  return 0;
}

int SubBox::set_nb14_potentials(const int** in_nb14_atomid_pairs,
				const int** in_nb14_atomtype_pairs,
				const real* in_nb14_coeff_vdw,
				const real* in_nb14_coeff_ele){
  n_nb14 = 0;
  for(int i=0; i < max_n_nb14; i++){
    int aid1_b = atomids_rev[in_nb14_atomid_pairs[i][0]];
    int aid2_b = atomids_rev[in_nb14_atomid_pairs[i][1]];
    if(aid1_b != -1 && aid2_b != -1){
      real mid_crd[3];
      for(int d=0; d < 3; d++){
	mid_crd[d] = (crd[aid1_b*3+d]
		      + crd[aid2_b*3+d]) * 0.5;
      }
      if(is_in_box(mid_crd)){
	nb14_atomid_pairs[n_nb14][0] = aid1_b;
	nb14_atomid_pairs[n_nb14][1] = aid2_b;
	nb14_atomtype_pairs[n_nb14][0] = in_nb14_atomtype_pairs[i][0];
	nb14_atomtype_pairs[n_nb14][1] = in_nb14_atomtype_pairs[i][1];
	nb14_coeff_vdw[n_nb14] = in_nb14_coeff_vdw[i];
	nb14_coeff_ele[n_nb14] = in_nb14_coeff_ele[i];
	n_nb14++;
      }
    }
  }
  return 0;
}

int SubBox::set_ele_excess(const int** in_excess_pairs){
  n_excess = 0;
  for(int i=0; i < max_n_excess; i++){
    int aid1_b = atomids_rev[in_excess_pairs[i][0]];
    int aid2_b = atomids_rev[in_excess_pairs[i][1]];
    if(aid1_b != -1 && aid2_b != -1){
      real mid_crd[3];
      for(int d=0; d < 3; d++){
	mid_crd[d] = (crd[aid1_b*3+d]
		      + crd[aid2_b*3+d]) * 0.5;
      }
      if(is_in_box(mid_crd)){
	excess_pairs[n_excess][0] = aid1_b;
	excess_pairs[n_excess][1] = aid2_b;
	n_excess++;
      }
    }
  }
  return 0;
}

int SubBox::set_nb15off(const int* in_nb15off){
  int nball = 0;
  for(int i=0; i < n_atoms* max_n_nb15off; i++)
    if(in_nb15off[i] != -1) nball++;
  for(int i=0; i < max_n_atoms_exbox * max_n_nb15off; i++)
    nb15off[i] = -1;
  int nb=0;
  for(int atomid = 0; atomid < n_atoms; atomid++){
    if(atomids_rev[atomid] > -1){
      for(int i = 0; i < max_n_nb15off; i++){
	int dest_atomid = in_nb15off[atomid * max_n_nb15off + i];
	if(dest_atomid > -1 && atomids_rev[dest_atomid] > -1){
	  int j;
	  for(j=0; j<max_n_nb15off; j++)
	    if(nb15off[atomids_rev[atomid] * max_n_nb15off + j] == -1) break;
	  nb15off[atomids_rev[atomid] * max_n_nb15off + j] = atomids_rev[dest_atomid];
	  //for(j=0; j<max_n_nb15off; j++)
	    //if(nb15off[atomids_rev[dest_atomid] * max_n_nb15off + j] == -1) break;
	  //nb15off[atomids_rev[dest_atomid] * max_n_nb15off + j] = atomids_rev[atomid];
	  //cout << "nb15off " << atomids_rev[atomid] * max_n_nb15off + n_off << " - "
	  //<<atomids_rev[dest] << " / " << n_atoms_exbox * max_n_nb15off << endl;;
	  nb+=1;
	}
      }
    }
  }
  //cout << "set_nb15off " << nb <<  " / " << nball << endl;
  return 0;
}

int SubBox::set_lj_param(const int in_n_lj_types,
			 real_pw* in_lj_6term,
			 real_pw* in_lj_12term){
  n_lj_types = in_n_lj_types;
  lj_6term = in_lj_6term;
  lj_12term = in_lj_12term;
  return 0;
}

int SubBox::calc_energy(){
  init_energy();
  init_work();

  const clock_t start_time_pair = clock();
#if defined(F_WO_NS)
  calc_energy_pairwise_wo_neighborsearch();
#elif defined(F_CUDA)
  calc_energy_pairwise_cuda();
#else
  calc_energy_pairwise();
#endif
  const clock_t end_time_pair = clock();
  ctime_calc_energy_pair += end_time_pair - start_time_pair;

  const clock_t start_time_bonded = clock();
  calc_energy_bonds();
  calc_energy_angles();
  calc_energy_torsions();
  calc_energy_impros();
  calc_energy_14nb();
  calc_energy_ele_excess();
  const clock_t end_time_bonded = clock();
  ctime_calc_energy_bonded += end_time_bonded - start_time_bonded;

#if defined(F_CUDA)
  cuda_pair_sync();
  cuda_memcpy_dtoh_work(nsgrid.get_work(),
			nsgrid.get_energy(),
			n_atoms_exbox,
			nsgrid.get_n_atom_array());
  cuda_thread_sync();
#endif
#if !defined(F_WO_NS)
  add_work_from_minicell();
  pote_vdw += nsgrid.get_energy()[0];
  pote_ele += nsgrid.get_energy()[1];
#endif

  for ( int i = 0; i < n_atoms_box*3; i++){
    work[i] *= -FORCE_VEL;
  }


  //cout << "nsgrid.get_energy()[1]" << nsgrid.get_energy()[1] << endl;


  //  cout << "SubBox::celc_energy excess:" << pote_ele << endl;
  //cout << "SubBox::calc_energy " << pote_vdw << " " << pote_ele << endl;
  return 0;
}
int SubBox::calc_energy_pairwise(){
  // for debug
  /*cout << " E : " << pote_vdw << ", " << pote_ele << endl;
    double sum_dist = 0.0;
    double sum_dist_incut = 0.0;
    int atomid1sum = 0;
    int atomid2sum = 0;
    int atomid12mult = 0;
    double lj6mult = 0.0;
    double lj12mult = 0.0;
    double chgmult = 0.0;
    int n_pairs=0;
    int n_pairs_incutoff=0;
    int n_pairs_15off = 0;
    int n_pairs_nonzero = 0;
    double p_vdw = 0.0;
    double p_ele = 0.0;
  */
  nsgrid.init_energy_work();
  CellPair* cellpairs = nsgrid.get_cell_pairs();
  for(int cp=0; cp < nsgrid.get_n_cell_pairs(); cp++){
    //CellPair cellpair = nsgrid.get_cell_pair(cp);
    
    int c1 = cellpairs[cp].cell_id1;
    int c2 = cellpairs[cp].cell_id2;
    int n_atoms_c1 = nsgrid.get_n_atoms_in_cell(c1);
    int n_atoms_c2 = nsgrid.get_n_atoms_in_cell(c2);
    int atoms_index_c1 = nsgrid.get_idx_cell_head_atom(c1);
    int atoms_index_c2 = nsgrid.get_idx_cell_head_atom(c2);
    int a2 = 0;
    for (a2=0; a2 < N_ATOM_CELL; a2++){
      int atomid_grid2 = atoms_index_c2 + a2;
      int atomid2 = nsgrid.get_atomid_from_gridorder(atomid_grid2);
      real_pw crd2[3];
      nsgrid.get_crd(atomid_grid2, crd2[0], crd2[1], crd2[2]);
      pbc->fix_pbc_image(crd2, cellpairs[cp].image);
      for (int a1=0; a1 < N_ATOM_CELL; a1++){
	int atomid_grid1 = atoms_index_c1 + a1;
	int atomid1 = nsgrid.get_atomid_from_gridorder(atomid_grid1);
	//n_pairs ++;
	int mask_id;
	int interact_bit;
	if (check_nb15off(a1, a2, cellpairs[cp].pair_mask, mask_id, interact_bit) ){ 
	  //n_pairs_15off++;
	  continue; }
	real_pw crd1[3];
	nsgrid.get_crd(atomid_grid1, crd1[0], crd1[1], crd1[2]);
	
	real_pw tmp_ene_vdw = 0.0;
	real_pw tmp_ene_ele = 0.0;
	real_fc tmp_work[3] = {0.0, 0.0, 0.0};
	real_pw param_6term  = lj_6term[atom_type[atomid1]  * n_lj_types + atom_type[atomid2]];
	real_pw param_12term = lj_12term[atom_type[atomid1] * n_lj_types + atom_type[atomid2]];
	
	//real_pw r12 = sqrt(pow(crd2[0]-crd1[0],2)+pow(crd2[1]-crd1[1],2)+pow(crd2[2]-crd1[2],2));
	//sum_dist += r12;
	//if(sum_dist > 100000) sum_dist -= 100000;
	
	real_pw r12 = ff.calc_pairwise(tmp_ene_vdw, tmp_ene_ele, tmp_work,
				       crd1, crd2,
				       param_6term, param_12term,
				       charge[atomid1],
				       charge[atomid2]);
	if(r12 > cfg->nsgrid_cutoff){
	  cellpairs[cp].pair_mask[mask_id] &= ~interact_bit;
	}
	
	nsgrid.add_energy(tmp_ene_vdw, tmp_ene_ele);
	/*
	  if(isnan(tmp_ene_vdw)){
	  cout << "Error! nonbond " << atomid1 << " "  << atomid2 <<" "
	  << " g:" << atomid_grid1 << " "  << atomid_grid2 <<" ";
	  cout << tmp_ene_vdw << " " <<  tmp_ene_ele << " " << r12 << " "
	  << crd1[0] << "," << crd1[1] << "," << crd1[2] <<" " 
	  << crd2[0] << "," << crd2[1] << "," << crd2[2] <<" " 
	  <<endl;
	  }
	*/
	nsgrid.add_work(atomid_grid1, tmp_work[0], tmp_work[1], tmp_work[2]);
	nsgrid.add_work(atomid_grid2, -tmp_work[0], -tmp_work[1], -tmp_work[2]);
	
	//n_pairs_incutoff++;
      /*
	p_vdw += tmp_ene_vdw;
	p_ele += tmp_ene_ele;
	if (tmp_ene_vdw != 0.0 || tmp_ene_ele != 0.0){
	  n_pairs_nonzero++;
	  sum_dist_incut += r12;
	  if(sum_dist_incut > 100000) sum_dist_incut -= 100000;	
	  if(atomid1 > atomid2){
	atomid1sum+=atomid2;
	atomid2sum+=atomid1;
	}else{
	atomid1sum+=atomid1;
	    atomid2sum+=atomid2;
	    }
	  lj6mult += param_6term;
	  while(lj6mult > 100000) lj6mult -= 100000;
	  lj12mult += param_12term;
	  while(lj12mult > 100000) lj12mult -= 100000;
	  chgmult += charge[atomid1] * charge[atomid2];
	  if(chgmult > 100000) chgmult -= 100000;
	  atomid12mult+=atomid1*atomid2;
	  atomid1sum = atomid1sum%100000;
	  atomid2sum = atomid2sum%100000;
	  atomid12mult = atomid12mult%100000;
	}
	*/
      }
    }
  }
  /*
  cout << "nb15off pairs " << n_pairs_15off << endl;
  cout << "15 pairs: " << n_pairs_nonzero << " / " << n_pairs_incutoff << " / " << n_pairs << endl;
  cout << " E : " << pote_vdw << ", " << pote_ele << endl;
  cout << " E2 : " << p_vdw << ", " << p_ele << endl;
  cout << " atomidsum : " << atomid1sum << " " << atomid2sum << " " << atomid1sum + atomid2sum << " " << atomid12mult << endl;
  cout << " sum_dist: " <<  sum_dist << " - " << sum_dist_incut << endl;
  cout << " lj6: " << lj6mult << " lj12: "<<lj12mult <<endl;
  cout << " chg: " << chgmult << endl;
  */
  return 0;
}
int SubBox::calc_energy_pairwise_wo_neighborsearch(){
  /*
    cout << " E : " << pote_vdw << ", " << pote_ele << endl;
    double sum_dist = 0.0;
    double sum_dist_incut = 0.0;
    int atomid1sum = 0;
    int atomid2sum = 0;
    int atomid12mult = 0;
    double lj6mult = 0.0;
    double lj12mult = 0.0;
    double chgmult = 0.0;
    int n_pairs=0;
    int n_pairs_incutoff=0;
    int n_pairs_15off = 0;
    double p_vdw = 0.0;
    double p_ele = 0.0;
  */
  for(int atomid1 = 0, atomid1_3=0; atomid1 < n_atoms_exbox; atomid1++, atomid1_3+=3){
    real_pw crd1[3] = {(real_pw)crd[atomid1_3], (real_pw)crd[atomid1_3+1], (real_pw)crd[atomid1_3+2]};
    for(int atomid2 = 0, atomid2_3=0; atomid2 < atomid1; atomid2++, atomid2_3+=3){
      //n_pairs++;
      real_pw crd2[3] = {(real_pw)crd[atomid2_3], (real_pw)crd[atomid2_3+1], (real_pw)crd[atomid2_3+2]};
      
      bool flg=true;
      for(int i = atomid1 * max_n_nb15off;
	  i < atomid1 * max_n_nb15off + max_n_nb15off;
	  i++){
	if(nb15off[i] == atomid2){
	  //n_pairs_15off++;
	  flg = false; 
	}
      }
      if(!flg) continue;

      real_pw tmp_ene_vdw = 0.0;
      real_pw tmp_ene_ele = 0.0;
      real_fc tmp_work[3] = {0.0, 0.0, 0.0};
      real_pw param_6term  = lj_6term[atom_type[atomid1]  * n_lj_types + atom_type[atomid2]];
      real_pw param_12term = lj_12term[atom_type[atomid1] * n_lj_types + atom_type[atomid2]];
      
      for(int d=0; d < 3; d++){
	if(crd2[d]-crd1[d] >= pbc->L_half[d])
	  crd2[d] -= pbc->L[d];
	else if(crd2[d]-crd1[d] <= -pbc->L_half[d])
	  crd2[d] += pbc->L[d];
      }
      
      //real_pw r12 = sqrt(pow(crd2[0]-crd1[0],2)+pow(crd2[1]-crd1[1],2)+pow(crd2[2]-crd1[2],2));
      //sum_dist += r12;
      //if(sum_dist > 100000) sum_dist -= 100000;
      
      if(ff.calc_pairwise(tmp_ene_vdw, tmp_ene_ele, tmp_work,
			  crd1, crd2,
			  param_6term, param_12term,
			  charge[atomid1],
			  charge[atomid2])==0){
	pote_vdw += tmp_ene_vdw;
	pote_ele += tmp_ene_ele;
	work[atomid1_3] += tmp_work[0];
	work[atomid1_3+1] += tmp_work[1];
	work[atomid1_3+2] += tmp_work[2];
	work[atomid2_3] -= tmp_work[0];
	work[atomid2_3+1] -= tmp_work[1];
	work[atomid2_3+2] -= tmp_work[2];
	//nsgrid.add_energy(tmp_ene_vdw, tmp_ene_ele);
	/*if(isnan(tmp_ene_vdw)){
	  cout << "Error! nonbond " << atomid1 << " "  << atomid2 <<" "
	       << " g:" << atomid_grid1 << " "  << atomid_grid2 <<" ";
	  cout << tmp_ene_vdw << " " <<  tmp_ene_ele << " " << r12 << " "
	       << crd1[0] << "," << crd1[1] << "," << crd1[2] <<" " 
	       << crd2[0] << "," << crd2[1] << "," << crd2[2] <<" " 
	       <<endl;
	}
	nsgrid.add_work(atomid_grid1, tmp_work[0], tmp_work[1], tmp_work[2]);
	  nsgrid.add_work(atomid_grid2, -tmp_work[0], -tmp_work[1], -tmp_work[2]);
	*/
	//n_pairs_incutoff++;
      }
      /*
	p_vdw += tmp_ene_vdw;
	p_ele += tmp_ene_ele;
	if (tmp_ene_vdw != 0.0 || tmp_ene_ele != 0.0){
	sum_dist_incut += r12;
	if(sum_dist_incut > 100000) sum_dist_incut -= 100000;	
	if(atomid1 > atomid2){
	atomid1sum+=atomid2;
	    atomid2sum+=atomid1;
	    }else{
	    atomid1sum+=atomid1;
	    atomid2sum+=atomid2;
	    }
	    lj6mult += param_6term;
	  while(lj6mult > 100000) lj6mult -= 100000;
	  lj12mult += param_12term;
	  while(lj12mult > 100000) lj12mult -= 100000;
	  chgmult += charge[atomid1] * charge[atomid2];
	  if(chgmult > 100000) chgmult -= 100000;
	  atomid12mult+=atomid1*atomid2;
	  atomid1sum = atomid1sum%100000;
	  atomid2sum = atomid2sum%100000;
	  atomid12mult = atomid12mult%100000;
	  }
      */
      
      
    }
  }
  /*
  cout << "nb15off pairs " << n_pairs_15off << endl;
  cout << "15 pairs: " << n_pairs_incutoff << " / " << n_pairs << endl;
  cout << " E : " << pote_vdw << ", " << pote_ele << endl;
  cout << " E2 : " << p_vdw << ", " << p_ele << endl;
  cout << " atomidsum : " << atomid1sum << " " << atomid2sum << " " << atomid1sum + atomid2sum << " " << atomid12mult << endl;
  cout << " sum_dist: " <<  sum_dist << " - " << sum_dist_incut << endl;
  cout << " lj6: " << lj6mult << " lj12: "<<lj12mult <<endl;
  cout << " chg: " << chgmult << endl;
  */
  return 0;
}

bool SubBox::check_nb15off(const int& a1, const int& a2, const int* bitmask,
			   int& mask_id, int& interact_bit){
  int bit_pos = a2 * N_ATOM_CELL + a1;
  mask_id =  bit_pos / 32;
  interact_bit = 1 << bit_pos % 32;
  return (bitmask[mask_id] & interact_bit) == interact_bit;
}

int SubBox::calc_energy_bonds(){
  for(int i=0; i < n_bonds; i++){
    //int i = bp_bonds[i];
    real_bp tmp_ene;
    real_bp tmp_work[3];
    int atomidx1 = bond_atomid_pairs[i][0]*3;
    int atomidx2 = bond_atomid_pairs[i][1]*3;
    real c1[3] = {crd[atomidx1], crd[atomidx1+1], crd[atomidx1+2]};
    real c2[3] = {crd[atomidx2], crd[atomidx2+1], crd[atomidx2+2]};
    ff.calc_bond(tmp_ene, tmp_work,
		 c1, c2,
		 bond_epsiron[i],
		 bond_r0[i]);
    pote_bond += tmp_ene;
    for(int d=0; d<3; d++){
      work[atomidx1+d] += tmp_work[d];
      work[atomidx2+d] -= tmp_work[d];
    }
    //cout << "BOND " << atomidx1 << "(" << atomids[atomidx1/3] << ") "
      //    << atomidx2 << "(" << atomids[atomidx2/3] << ") "
      //    << pote_bond << endl;
  }
  return 0;
}
int SubBox::calc_energy_angles(){
  for(int i=0; i < n_angles; i++){
    //  int angle_idx = bp_angles[i];
    real tmp_ene;
    real_fc tmp_work1[3], tmp_work2[3];

    int atomidx1 = angle_atomid_triads[i][0]*3;
    int atomidx2 = angle_atomid_triads[i][1]*3;
    int atomidx3 = atomids_rev[angle_atomid_triads[i][2]]*3;
    real c1[3] = {crd[atomidx1], crd[atomidx1+1], crd[atomidx1+2]};
    real c2[3] = {crd[atomidx2], crd[atomidx2+1], crd[atomidx2+2]};
    real c3[3] = {crd[atomidx3], crd[atomidx3+1], crd[atomidx3+2]};
    ff.calc_angle(tmp_ene, tmp_work1, tmp_work2,
		  c1, c2, c3,
		  angle_epsiron[i],
		  angle_theta0[i]);
    pote_angle += tmp_ene;
    for(int d=0; d<3; d++){
      work[atomidx1+d] += tmp_work1[d];
      work[atomidx2+d] -= tmp_work1[d] + tmp_work2[d];
      work[atomidx3+d] += tmp_work2[d];
    }
  }
  return 0;
}

int SubBox::calc_energy_torsions(){
  for(int i=0; i < n_torsions; i++){
    //int torsion_idx = bp_torsions[i];
    real tmp_ene;
    real_fc tmp_work1[3], tmp_work2[3], tmp_work3[3];
    int atomidx1 = torsion_atomid_quads[i][0]*3;
    int atomidx2 = torsion_atomid_quads[i][1]*3;
    int atomidx3 = torsion_atomid_quads[i][2]*3;
    int atomidx4 = torsion_atomid_quads[i][3]*3;
    real c1[3] = {crd[atomidx1], crd[atomidx1+1], crd[atomidx1+2]};
    real c2[3] = {crd[atomidx2], crd[atomidx2+1], crd[atomidx2+2]};
    real c3[3] = {crd[atomidx3], crd[atomidx3+1], crd[atomidx3+2]};
    real c4[3] = {crd[atomidx4], crd[atomidx4+1], crd[atomidx4+2]};
    ff.calc_torsion(tmp_ene, tmp_work1, tmp_work2, tmp_work3,
		    c1, c2, c3, c4,
		    torsion_energy[i],
		    torsion_overlaps[i],
		    torsion_symmetry[i],
		    torsion_phase[i]);
    pote_torsion += tmp_ene;
    for(int d=0; d<3; d++){
      work[atomidx1+d] += tmp_work1[d];
      work[atomidx2+d] += tmp_work2[d];
      work[atomidx3+d] -= tmp_work1[d] + tmp_work2[d] + tmp_work3[d];
      work[atomidx4+d] += tmp_work3[d];
    }
  }

  return 0;
}
int SubBox::calc_energy_impros(){
  for(int i=0; i < n_impros; i++){
    //int impro_idx = bp_impros[i];
    real tmp_ene;
    real_fc tmp_work1[3], tmp_work2[3], tmp_work3[3];
    int atomidx1 = impro_atomid_quads[i][0]*3;
    int atomidx2 = impro_atomid_quads[i][1]*3;
    int atomidx3 = impro_atomid_quads[i][2]*3;
    int atomidx4 = impro_atomid_quads[i][3]*3;
    real c1[3] = {crd[atomidx1], crd[atomidx1+1], crd[atomidx1+2]};
    real c2[3] = {crd[atomidx2], crd[atomidx2+1], crd[atomidx2+2]};
    real c3[3] = {crd[atomidx3], crd[atomidx3+1], crd[atomidx3+2]};
    real c4[3] = {crd[atomidx4], crd[atomidx4+1], crd[atomidx4+2]};
    ff.calc_torsion(tmp_ene, tmp_work1, tmp_work2, tmp_work3,
		    c1, c2, c3, c4,
		    impro_energy[i],
		    impro_overlaps[i],
		    impro_symmetry[i],
		    impro_phase[i]);
    pote_impro += tmp_ene;
    for(int d=0; d<3; d++){
      work[atomidx1+d] += tmp_work1[d];
      work[atomidx2+d] += tmp_work2[d];
      work[atomidx3+d] -= tmp_work1[d] + tmp_work2[d] + tmp_work3[d];
      work[atomidx4+d] += tmp_work3[d];
    }
  }
  return 0;
}
int SubBox::calc_energy_14nb(){
  for(int i=0; i < n_nb14; i++){
    real tmp_ene_vdw;
    real tmp_ene_ele;
    real_fc tmp_work[3];
    int atomidx1 = nb14_atomid_pairs[i][0]*3;
    int atomidx2 = nb14_atomid_pairs[i][1]*3;
    int atomtype1 = nb14_atomtype_pairs[i][0];
    int atomtype2 = nb14_atomtype_pairs[i][1];
    real c1[3] = {crd[atomidx1], crd[atomidx1+1], crd[atomidx1+2]};
    real c2[3] = {crd[atomidx2], crd[atomidx2+1], crd[atomidx2+2]};
    ff.calc_14pair(tmp_ene_vdw, tmp_ene_ele,
		   tmp_work,
		   c1, c2,
		   lj_6term[atomtype1 * n_lj_types + atomtype2],
		   lj_12term[atomtype1 * n_lj_types + atomtype2],
		   charge[nb14_atomid_pairs[i][0]],
		   charge[nb14_atomid_pairs[i][1]],
		   nb14_coeff_vdw[i],
		   nb14_coeff_ele[i]
		   );
    pote_14vdw += tmp_ene_vdw;
    pote_14ele += tmp_ene_ele;
    for(int d=0; d<3; d++){
      work[atomidx1+d] += tmp_work[d];
      work[atomidx2+d] -= tmp_work[d];
    }
  }
  return 0;
}
int SubBox::calc_energy_ele_excess(){
  real tmp = pote_ele;
  for (int i=0; i < n_excess; i++){
    real tmp_ene;
    real_fc tmp_work[3];
    int atomidx1 = excess_pairs[i][0]*3;
    int atomidx2 = excess_pairs[i][1]*3;
    real c1[3] = {crd[atomidx1], crd[atomidx1+1], crd[atomidx1+2]};
    real c2[3] = {crd[atomidx2], crd[atomidx2+1], crd[atomidx2+2]};
    ff.calc_zms_excess(tmp_ene, tmp_work,
		       c1, c2,
		       charge[excess_pairs[i][0]],
		       charge[excess_pairs[i][1]]);
    pote_ele += tmp_ene;
    for(int d=0; d<3; d++){
      work[atomidx1+d] += tmp_work[d];
      work[atomidx2+d] -= tmp_work[d];
    }
    //cout << "Excess : " << i << " "<<  pote_ele - tmp <<  " " << pote_ele << " " << tmp << endl;
  }
  //cout << "Excess : "  << pote_ele - tmp <<  " " << pote_ele << " " << tmp << endl;
  return 0;
}

int SubBox::get_box_id_from_crd(const int box_crd[]){
  return n_boxes_xyz[0]*n_boxes_xyz[1]*box_crd[2] +
    n_boxes_xyz[0]*box_crd[1] + box_crd[0];
}

int SubBox::get_box_crd_from_id(const int box_id, 
				int *box_crd){
  box_crd[0] = box_id%(n_boxes_xyz[0]);
  box_crd[1] = (box_id%(n_boxes_xyz[0]*n_boxes_xyz[1])) / n_boxes_xyz[0];
  box_crd[2] = (box_id/(n_boxes_xyz[0]*n_boxes_xyz[1]));
  return 0;
}

int SubBox::init_energy(){
  pote_bond = 0.0;
  pote_angle = 0.0;
  pote_torsion = 0.0;
  pote_impro = 0.0;
  pote_14vdw = 0.0;
  pote_14ele = 0.0;
  pote_vdw = 0.0;
  pote_ele = 0.0;
  return 0;
}

int SubBox::init_work(){
  int n_atoms3 = all_n_atoms[rank]*3;
  for(int atomid_b3=0;
      atomid_b3 < n_atoms3;
      atomid_b3++){
    work[atomid_b3] = 0.0;
  }
  return 0;
}

int SubBox::add_work_from_minicell(){
  //real_fc*& in_work,
  //int*& in_atomids_rev){
  for(int atomid_b=0, atomid_b3=0;
      atomid_b < all_n_atoms[rank];
      atomid_b++, atomid_b3+=3){
    int atomid_g = nsgrid.get_atomids_rev()[atomids[atomid_b]];
    int atomid_g3 = atomid_g * 3;
    for(int d=0; d < 3; d++){
      work[atomid_b3+d] += nsgrid.get_work()[atomid_g3+d];
    }
  }
  return 0;
}

int SubBox::cpy_crd_prev(){
  memcpy(crd_prev, crd, sizeof(real) * max_n_atoms_exbox*3);
  return 0;
}
int SubBox::cpy_crd_from_prev(){
  //memcpy(crd, crd_prev, sizeof(real) * max_n_atoms_exbox*3);
  for( int i = 0; i < n_atoms_box *3 ; i ++){
    crd[i] = crd_prev[i];
  }
  return 0;
}
int SubBox::cpy_vel_prev(){
  memcpy(vel, vel_next, sizeof(real) * max_n_atoms_exbox*3);
  return 0;
}
int SubBox::cpy_vel_buf_from_prev(){
  memcpy(buf_crd, vel, sizeof(real) * max_n_atoms_exbox*3);
  return 0;
}
int SubBox::cpy_vel_prev_from_buf(){
  memcpy(vel, buf_crd, sizeof(real) * max_n_atoms_exbox*3);
  return 0;
}
int SubBox::swap_velocity_buffer(){
  real* tmp = vel;
  vel = vel_next;
  vel_next = tmp;
  return 0;
}

int SubBox::update_velocities(const real time_step){
  //swap_velocity_buffer();
  for(int atomid_b=0, atomid_b3=0;
      atomid_b < all_n_atoms[rank];
      atomid_b++, atomid_b3+=3){
    for(int d=0; d<3; d++){    
      vel_next[atomid_b3+d] = vel[atomid_b3+d]  +
	( time_step * 
	 //	 FORCE_VEL * 
	 work[atomid_b3+d] 
	  * mass_inv[atomid_b]);
    }
  }
  return 0;
}

int SubBox::velocity_average(){
  for(int atomid_b=0, atomid_b3=0;
      atomid_b < all_n_atoms[rank];
      atomid_b++, atomid_b3+=3){
    for(int d=0; d<3; d++){    
      vel_just[atomid_b3+d] = (vel[atomid_b3+d] + vel_next[atomid_b3+d])*0.5;
    }
  }
  return 0;
}
int SubBox::set_velocity_from_crd(){
  
  for(int atomid_b = 0, atomid_b3 = 0;
      atomid_b < all_n_atoms[rank];
      atomid_b++, atomid_b3+=3){
    real norm1 = 0.0;
    real norm2 = 0.0;
    real d_crd[3];
    //real crd1[3] = {crd[atomid_b3], crd[atomid_b3+1], crd[atomid_b3+2]};
   //real crd2[3] = {crd_prev[atomid_b3], crd_prev[atomid_b3+1], crd_prev[atomid_b3+2]};
    //pbc->diff_crd_minim_image(d_crd,crd1,crd2);

    for(int d=0; d<3; d++){
      norm1 += vel_next[atomid_b3+d] * vel_next[atomid_b3+d];
      vel_next[atomid_b3+d] = (crd[atomid_b3+d] - crd_prev[atomid_b3+d]) * time_step_inv;
      //vel_next[atomid_b3+d] = d_crd[d] * ts_inv;
      norm2 += vel_next[atomid_b3+d] * vel_next[atomid_b3+d];
    }
    real diff = fabs(norm1 - norm2)*mass[atomid_b];
    if(diff > 0.01)
      cout << "diff " << atomid_b << " " << diff
	   << "(" << crd[atomid_b3] << ", " << crd[atomid_b3+1] << ", " << crd[atomid_b3+2] << ") "
	   << "(" << crd_prev[atomid_b3] << ", " << crd_prev[atomid_b3+1] << ", " << crd_prev[atomid_b3+2] << ") "
	   << endl;
  }
  return 0;
}
int SubBox::revise_coordinates_pbc(){
  for(int atomid_b=0, atomid_b3=0;
      atomid_b < all_n_atoms[rank];
      atomid_b++, atomid_b3+=3){
    for(int d=0; d<3; d++){    
      if(crd[atomid_b3+d] >= pbc->upper_bound[d]){
	crd[atomid_b3+d] -= pbc->L[d];
      }else if(crd[atomid_b3+d] < pbc->lower_bound[d]){
	crd[atomid_b3+d] += pbc->L[d];
      }
    }
  }
  return 0;
}

int SubBox::copy_vel_just(real** p_vel){
  for(int atomid_b=0, atomid_b3=0;
      atomid_b < all_n_atoms[rank];
      atomid_b++, atomid_b3+=3){
    for(int d=0; d<3; d++){
      p_vel[atomids[atomid_b]][d] = (vel[atomid_b3+d] + vel_next[atomid_b3+d])*0.5;
    }
  }
  return 0;
  /// copy velocities from this object to the mmsys
  //  this.vel_just => mmsys.p_vel
  //  it called from DynamicsMoce.cpp
  //     subbox.set_vel_just(mmsys.vel)
  //return get_vel(vel_just, p_vel);
}
int SubBox::set_force_from_velocity(const real time_step){
  time_step_inv = 1.0 / time_step;
  for(int atomid_b = 0, atomid_b3 = 0;
      atomid_b < all_n_atoms[rank];
      atomid_b++, atomid_b3+=3){
    for(int d=0; d<3; d++){
      work[atomid_b3+d] = (vel_next[atomid_b3+d] - vel[atomid_b3+d])*mass[atomid_b] * time_step_inv;
    }
  }  
  return 0;
}
int SubBox::copy_crd(real** p_crd){
  return copy_crdvel_to_mmsys(crd, p_crd);
}
int SubBox::copy_vel(real** p_vel){
  return copy_crdvel_to_mmsys(vel, p_vel);
}
int SubBox::copy_vel_next(real** p_vel){
  return copy_crdvel_to_mmsys(vel_next, p_vel);
}
int SubBox::copy_crdvel_to_mmsys(real* src, real** dst){
  for(int atomid_b=0, atomid_b3=0;
      atomid_b < all_n_atoms[rank];
      atomid_b++, atomid_b3+=3){
    for(int d=0; d<3; d++){
      dst[atomids[atomid_b]][d] = src[atomid_b3+d];
    }
  }
  return 0;
}
int SubBox::add_force_from_mmsys(real_fc** in_force){
  add_crdvel_from_mmsys(work, in_force);
  return 0;
}
template <typename TYPE> int SubBox::add_crdvel_from_mmsys(TYPE* dst, TYPE** src){
  for(int atomid_b=0, atomid_b3=0;
      atomid_b < all_n_atoms[rank];
      atomid_b++, atomid_b3+=3){
    for(int d=0; d<3; d++){
      dst[atomid_b3+d] += src[atomids[atomid_b]][d];
    }
  }
  return 0;
}
int SubBox::update_coordinates_prev(const real time_step){
  update_coordinates(time_step, crd_prev, vel);
  return 0;
}
int SubBox::update_coordinates_cur(const real time_step){
  update_coordinates(time_step, crd, vel_next);
  return 0;
}

int SubBox::update_coordinates(const real time_step, real *p_crd, real* p_vel){
  //subbox.cpy_crd_prev();
  for(int atomid_b=0, atomid_b3=0;
      atomid_b < all_n_atoms[rank];
      atomid_b++, atomid_b3+=3){
    for(int d=0; d<3; d++){
      real diff = p_vel[atomid_b3+d] * time_step;
      crd[atomid_b3+d] = p_crd[atomid_b3+d] + diff;
      //#ifndef F_WO_NS      
      //      nsgrid.move_atom(atomid_b, d, diff);
      //#endif
    }
  }
  return 0;
}

int SubBox::update_coordinates_nsgrid(){
  // this method should be called 
  //   after update_coordinates(), before revise_coordinates_pbc()

  for(int atomid_b=0, atomid_b3=0;
      atomid_b < all_n_atoms[rank];
      atomid_b++, atomid_b3+=3){
    //real d_crd[3];
    //real crd1[3] = {crd[atomid_b3], crd[atomid_b3+1], crd[atomid_b3+2]};
    //real crd2[3] = {crd_prev[atomid_b3], crd_prev[atomid_b3+1], crd_prev[atomid_b3+2]};
    //pbc->diff_crd_minim_image(d_crd,crd1,crd2);
    for(int d=0; d<3; d++){
      nsgrid.move_atom(atomid_b, d, crd[atomid_b3+d]-crd_prev[atomid_b3+d]);      
      //nsgrid.move_atom(atomid_b, d, d_crd[d]);
    }
  }

  return 0;
}

bool SubBox::is_in_box(real* in_crd){
  return in_crd[0]>=box_lower[0] && in_crd[0]<box_upper[0] &&
    in_crd[1]>=box_lower[1] && in_crd[1]<box_upper[1] &&
    in_crd[2]>=box_lower[2] && in_crd[2]<box_upper[2];
}
bool SubBox::is_in_exbox(real* in_crd){
  for(int d=0; d<3; d++){
    if(in_crd[d] < exbox_lower[d]) return false;
    if(in_crd[d] >= exbox_upper[d]) return false;
  }
  return true;
}
int SubBox::set_box_crd(){
  box_crd[0] = rank % n_boxes_xyz[0];
  box_crd[2] = rank / (n_boxes_xyz[0] * n_boxes_xyz[1]);
  box_crd[1] = (rank % (n_boxes_xyz[0] * n_boxes_xyz[1])) / n_boxes_xyz[0];
  return 0;
}

int SubBox::init_constraint(int in_constraint,
			    int in_max_loops, real in_tolerance,
			    int max_n_pair,
			    int max_n_trio,
			    int max_n_quad,
			    int max_n_settle){
  switch (in_constraint){
  case CONST_SHAKE:
    flg_constraint = 1;
    break;
  case CONST_SHAKE_SETTLE:    
    flg_constraint = 1;
    flg_settle = 1;
    break;
  default:
    return in_constraint;
  }
  if (flg_constraint != 0){
    constraint = new ConstraintShake();
    constraint->set_parameters(in_max_loops, in_tolerance);
    constraint->set_max_n_constraints(max_n_pair, max_n_trio, max_n_quad);
    //cout << "max_n_const: " << max_n_pair << " " << max_n_trio << " " << max_n_quad  << endl;
    constraint->alloc_constraint();
  }else{
    constraint = new ConstraintObject();
  }
  if (flg_settle != 0){
    settle = new ConstraintSettle();
    settle->set_max_n_constraints(0, max_n_settle, 0);
    settle->alloc_constraint();
  }else{
    settle = new ConstraintObject();
  }

  return in_constraint;
}

int SubBox::set_subset_constraint(ConstraintObject& in_cst,
				  ConstraintObject& in_settle){
  if(flg_constraint != 0)
    //cout << "set_subset_constraint" << endl;
    constraint->set_subset_constraint(in_cst, atomids_rev);

  if(flg_settle != 0)
    settle->set_subset_constraint(in_settle, atomids_rev);

  //temperature_coef = (2.0 * JOULE_CAL * 1e+3) / (GAS_CONST * d_free);
  return 0;
}
int SubBox::init_thermostat(const int in_thermostat_type,
			    const real in_temperature_init,
			    const int d_free){
  if(in_thermostat_type == THMSTT_NONE){
    thermostat = new ThermostatObject();
  }else if(in_thermostat_type == THMSTT_SCALING){
    thermostat = new ThermostatScaling();
  }else if(in_thermostat_type == THMSTT_HOOVER_EVANS){
    thermostat = new ThermostatHooverEvans();
  }
  thermostat->set_temperature(in_temperature_init);
  thermostat->set_temperature_coeff(d_free);

  thermostat->set_time_step(cfg->time_step);
  thermostat->set_constant(n_atoms_box, mass_inv, vel_next, work);
  return 0;
}

#ifdef F_CUDA
int SubBox::gpu_device_setup(){
  
  //cuda_print_device_info();
  //cuda_memcpy_htod_grid_pairs(mmsys.nsgrid.grid_pairs,
  //mmsys.nsgrid.n_grid_pairs);
  cuda_alloc_atom_info(max_n_atoms_exbox,
		       //nsgrid.get_max_n_atom_array(),
		       nsgrid.get_max_n_cells(),
		       nsgrid.get_max_n_cell_pairs(),
		       nsgrid.get_n_columns()+1);
		       //     nsgrid.get_n_uni());

  cuda_alloc_set_lj_params(lj_6term,
			   lj_12term,
			   n_lj_types,
			   nb15off,
			   max_n_nb15off,
			   max_n_atoms_exbox,
			   nsgrid.get_max_n_atom_array());

  
  real_pw tmp_l[3] = {(real_pw)pbc->L[0], (real_pw)pbc->L[1], (real_pw)pbc->L[2]};
  real_pw tmp_lb[3] = {(real_pw)pbc->lower_bound[0], (real_pw)pbc->lower_bound[1], (real_pw)pbc->lower_bound[2]};
  //cuda_set_pbc((const real_pw*)pbc->L);
  cuda_set_pbc(tmp_l, tmp_lb);
  cuda_set_constant(
		    (real_pw)cfg->cutoff,
		    (real_pw)cfg->nsgrid_cutoff,
		    n_lj_types);
  cuda_zerodipole_constant(ff.ele->get_zcore(),
			   ff.ele->get_bcoeff(),
			   ff.ele->get_fcoeff());

  return 0;
}

int SubBox::update_device_cell_info(){
  cuda_set_cell_constant(nsgrid.get_n_cells(),
			 n_atoms_exbox,
			 nsgrid.get_n_atom_array(),
			 nsgrid.get_n_cells_xyz(),
			 nsgrid.get_n_columns(),
			 nsgrid.get_L_cell_xyz(),
			 nsgrid.get_n_neighbors_xyz());

  cuda_memcpy_htod_atomids(nsgrid.get_atomids(),
			   nsgrid.get_idx_xy_head_cell());


  //cuda_init_cellinfo(nsgrid.get_n_cells());
  cuda_set_atominfo();

  return 0;
}

int SubBox::calc_energy_pairwise_cuda(){
  nsgrid.init_energy_work();
  cuda_pairwise_ljzd(flg_mod_15mask);
  flg_mod_15mask = false;
  return 0;
}

#endif
int SubBox::apply_constraint(){
  constraint->apply_constraint(crd, crd_prev, mass_inv, pbc);
  set_velocity_from_crd();
  return 0;
}

int SubBox::update_thermostat(const int cur_step){
  //cout << "update thermostat " << cur_step << " / " << cfg->heating_steps<< endl;
  if(cur_step < cfg->heating_steps){
    real new_temp = cfg->temperature_init + 
      ((cfg->temperature - cfg->temperature_init) /
       (real)cfg->heating_steps) * cur_step;
    thermostat->set_temperature(new_temp);
    //cout << "set temperature : " << new_temp << " " << cfg->temperature_init << " - " << cfg->temperature << endl;
  }else{
    thermostat->set_temperature(cfg->temperature);
  }
  return 0;
}

int SubBox::apply_thermostat(){

  thermostat->apply_thermostat(n_atoms_box,
			       work, vel, vel_next,
			       mass, mass_inv);
  return 0;
}

int SubBox::apply_thermostat_with_shake(const int max_loops,
					const real tolerance){
  thermostat->apply_thermostat_with_shake(n_atoms_box,
					  work, crd, crd_prev,
					  vel, vel_next,
					  mass, mass_inv,
					  constraint,
					  pbc, buf_crd,
					  max_loops, tolerance,
					  &commotion, atomids_rev);

  return 0;
}
int SubBox::expand_apply_bias(unsigned long cur_step,  real in_lambda){
  expand->apply_bias(cur_step, in_lambda, work, n_atoms_box);
  return 0;
}
int SubBox::expand_apply_bias_struct_param(unsigned long cur_step){
  real param = expand->cal_struct_parameters(crd, pbc);
  expand->apply_bias(cur_step, param, work, n_atoms_box);
  return 0;
}
void SubBox::expand_enable_vs_transition(){
  expand->enable_vs_transition();
}

int SubBox::cancel_com_motion(){
  return commotion.cancel_translation(atomids_rev, vel_next);
}

int SubBox::set_com_motion(int n_groups, int* group_ids,
			   int*  n_atoms_in_groups, 
			   int** groups,
			   real* mass_inv_groups){
  return commotion.set_groups(n_groups,
			      group_ids, n_atoms_in_groups,
			      groups, mass_inv_groups,
			      mass);
  
}
/*
int SubBox::set_bonding_info
( int** in_bond_atomid_pairs,
  real* in_bond_epsiron,
  real* in_bond_r0,
  int** in_angle_atomid_triads,
  real* in_angle_epsiron,
  real* in_angle_theta0,
  int** in_torsion_atomid_quads,
  real* in_torsion_energy,
  int* in_torsion_overlaps,
  int* in_torsion_symmetry,
  real* in_torsion_phase,
  int* in_torsion_nb14,
  int** in_impro_atomid_quads,
  real* in_impro_energy,
  int* in_impro_overlaps,
  int* in_impro_symmetry,
  real* in_impro_phase,
  int* in_impro_nb14,
  int** in_nb14_atomid_pairs,
  int** in_nb14_atomtype_pairs,
  real* in_nb14_coeff_vdw,
  real* in_nb14_coeff_ele
){
  bond_atomid_pairs = in_bond_atomid_pairs;
  bond_epsiron = in_bond_epsiron;
  bond_r0 = in_bond_r0;
  angle_atomid_triads = in_angle_atomid_triads;
  angle_epsiron = in_angle_epsiron;
  angle_theta0 = in_angle_theta0;
  torsion_atomid_quads = in_torsion_atomid_quads;
  torsion_energy = in_torsion_energy;
  torsion_overlaps = in_torsion_overlaps;
  torsion_symmetry = in_torsion_symmetry;
  torsion_phase = in_torsion_phase;
  torsion_nb14 = in_torsion_nb14;
  impro_atomid_quads = in_impro_atomid_quads;
  impro_energy = in_impro_energy;
  impro_overlaps = in_impro_overlaps;
  impro_symmetry = in_impro_symmetry;
  impro_phase = in_impro_phase;
  impro_nb14 = in_impro_nb14;
  nb14_atomid_pairs = in_nb14_atomid_pairs;
  nb14_atomtype_pairs = in_nb14_atomtype_pairs;
  nb14_coeff_vdw = in_nb14_coeff_vdw;
  nb14_coeff_ele = in_nb14_coeff_ele;
  return 0;
}
*/
/*
qint SubBox::assign_regions(){
  // set 
  //   region_atoms
  //   n_region_atoms
  // require
  //   box, box_l

  for(int atom_id_box = 0; atom_id_box < n_atoms_box; atom_id_box++){
    int region[3];
    for(int d = 0; d < 3; d++){
      real tmp_crd = (pbc->lower_bound[d] + crd[atom_id_box*3+d]);
	box[d] = floor(tmp_crd / box_l[d]);
      if( tmp_crd < pbc->lower_bound[d] + box[d]*box_l[d] + cutoff_pair_half )
	region[d] = -1;
      else if( tmp_crd >= pbc->lower_bound[d] + (box[d]+1)*box_l[d] - cutoff_pair_half )
	region[d] = 1;
      else region[d] = 0;
    }
    //int boxid = get_box_id_from_crd(box[0],box[1],box[2]);
    int regionid = get_region_id_from_crd(3, region[0],region[1],region[2]);
    //int g_regionid = get_global_region_id_from_box_region(boxid, regionid);
    //global_region_atoms[g_regionid][n_global_region_atoms[g_regionid]] = atom_id_box;
    n_region_atoms[regionid]++;
  }
  return 0;
}

int SubBox::set_max_n_atoms_region(){
  // set
  //   max_n_atoms_region
  //   max_n_atoms_box
  //   max_n_atoms_exbox 
  max_n_atoms_box = 0;
  int tmp_max_n_atoms_box = (n_atoms + n_boxes-1)/n_boxes * COEF_MAX_N_ATOMS;
  real vol_box = box_l[0] * box_l[1] * box_l[2];
  for(int i=0; i<27; i++){
    int rx[3];
    get_region_crd_from_id(3, i, rx[0], rx[1], rx[2]);

    int nx[3];
    for(int j=0; j<3; j++){
      if(rx[j]==0)
	nx[j] = box_l[j] - cutoff_pair;
      else
	nx[j] = cutoff_pair*0.5;
    }
    real vol = nx[0] * nx[1] * nx[2];
    max_n_atoms_region[i] = vol/vol_box * tmp_max_n_atoms_box;
    max_n_atoms_box += max_n_atoms_region[i];
    cout << "region_id : " << i << " (" <<get_region_id_from_crd(3, rx[0],rx[1],rx[2])<<") "
	 << " [" <<rx[0]<<","<<rx[1]<<","<<rx[2]<<"]" <<endl;
      //<< " max_n_cell_reg:"<<max_n_cells_region[i]<< endl;
  }
  max_n_atoms_exbox = max_n_atoms_box;
  // -1,-1,-1 => 56 regions
  max_n_atoms_exbox += max_n_atoms_region[get_region_id_from_crd(3, -1, -1, -1)] * 56;
  // -1,-1, 0 => 12 regions
  max_n_atoms_exbox += max_n_atoms_region[get_region_id_from_crd(3, -1, -1, 0)] * 12;
  // -1, 0,-1 => 12 regions
  max_n_atoms_exbox += max_n_atoms_region[get_region_id_from_crd(3, -1, 0, -1)] * 12;
  //  0,-1,-1 => 12 regions
  max_n_atoms_exbox += max_n_atoms_region[get_region_id_from_crd(3, 0, -1, -1)] * 12;
  // -1, 0, 0 => 2
  max_n_atoms_exbox += max_n_atoms_region[get_region_id_from_crd(3, -1,  0,  0)] * 2;
  //  0,-1, 0 => 2
  max_n_atoms_exbox += max_n_atoms_region[get_region_id_from_crd(3, 0, -1,  0)] * 2;
  //  0, 0,-1 => 2
  max_n_atoms_exbox += max_n_atoms_region[get_region_id_from_crd(3,  0,  0, -1)] * 2;
  return 0;
}
int SubBox::get_region_id_from_crd(int width, int rx, int ry, int rz){
  // for 27  regions: width = 3;
  // for 125 regions: width = 5;
  return (rx+width/2) * width*width + (ry+width/2)*width + rz+width/2;
}
int SubBox::get_region_crd_from_id(int width, int regid,
				     int& rx, int& ry, int& rz){
  int width2 = width*width;
  int width_d2= width / 2;
  rx = regid/(width2) - width_d2;
  ry = (regid%width2)/width - width_d2;
  rz = regid%width - width_d2;
  return 0;
}

int SubBox::get_global_region_id_from_box_region(int boxid, int regionid){
  return boxid*27+regionid;
}

*/
