#ifndef __SUB_BOX_H__
#define __SUB_BOX_H__

#include "CelesteObject.h"
#include "MiniCell.h"
#include "PBC.h"
#include "Config.h"
#include "ForceField.h"
#include <ctime>
using namespace std;

#define COEF_MAX_N_ATOMS_BOX 1.2

class SubBox : public CelesteObject {
 private:
  // only for rank=0
  int n_atoms;
  
  int n_boxes;
  int n_boxes_xyz[3];
  real box_l[3];
  real exbox_l[3];
  real box_lower[3];
  real box_upper[3];
  real exbox_lower[3];
  real exbox_upper[3];
  real cutoff_pair;
  real cutoff_pair_half;
  PBC* pbc;
  Config* cfg;
  int box_crd[3];

  ForceField ff;

  // the maximum number of atoms in 27 region (-1 ~ +1)
  int max_n_atoms_box;
  // the maximum number of atoms in 125 regions (-2 ~ +2)
  int max_n_atoms_exbox;

  //int **region_atoms;
  //int *n_region_atoms;
  //int max_n_atoms_region[27];

  // INFO for each box

  // crd, force, atomids in each box
  // crd[max_n_atoms_box*3]
  int rank;
  real *crd;
  real *vel;
  real *vel_next;
  real *vel_just;
  real *work;
  real *charge;
  int  *atom_type;

  int  *atomids;
  // atomids_rev[atomid] = -1: it is not in the box
  int  *atomids_rev;
  int n_atoms_box;
  int n_atoms_exbox;

  // variables for rank0
  // all_atomids[box_id][atomid in the box] = original atomid
  //    set in the function SubBox::rnak0_div_box()
  int** all_atomids;
  // all_n_atoms
  int* all_n_atoms;

  // for bonded potential
  int* bp_bonds;
  int* bp_angles;
  int* bp_torsions;
  int* bp_impros;
  int* bp_nb14;
  int max_n_bonds;

  int max_n_angles;
  int max_n_torsions;
  int max_n_impros;
  int max_n_nb14;
  int max_n_excess;
  int n_bonds;
  int n_angles;
  int n_torsions;
  int n_impros;
  int n_nb14;
  int n_excess;

  int** bond_atomid_pairs;
  real* bond_epsiron;
  real* bond_r0;
  int** angle_atomid_triads;
  real* angle_epsiron;
  real* angle_theta0;
  int** torsion_atomid_quads;
  real* torsion_energy;
  int* torsion_overlaps;
  int* torsion_symmetry;
  real* torsion_phase;
  int* torsion_nb14;
  int** impro_atomid_quads;
  real* impro_energy;
  int* impro_overlaps;
  int* impro_symmetry;
  real* impro_phase;
  int* impro_nb14;
  int** nb14_atomid_pairs;
  int** nb14_atomtype_pairs;
  real* nb14_coeff_vdw;
  real* nb14_coeff_ele;
  int** excess_pairs;
  
  int max_n_nb15off;
  int n_nb15off;
  int* nb15off;

  int n_lj_types;
  real_pw* lj_6term;
  real_pw* lj_12term;
  
  real_fc pote_vdw;
  real_fc pote_ele;
  real_fc pote_bond;
  real_fc pote_angle;
  real_fc pote_torsion;
  real_fc pote_impro;
  real_fc pote_14ele;
  real_fc pote_14vdw;

  MiniCell nsgrid;

  clock_t ctime_setgrid;
  clock_t ctime_enumerate_cellpairs;
  clock_t ctime_calc_energy_pair;
  clock_t ctime_calc_energy_bonded;

 public:
  SubBox();
  ~SubBox();
  int alloc_variables();
  int init_variables();
  int alloc_variables_for_bonds(int in_n_bonds);
  int alloc_variables_for_angles(int in_n_angles);
  int alloc_variables_for_torsions(int in_n_torsions);
  int alloc_variables_for_impros(int in_n_impros);
  int alloc_variables_for_nb14(int in_n_nb14);
  int alloc_variables_for_excess(int in_n_excess);
  int alloc_variables_for_nb15off(int in_max_n_nb15off);
  int free_variables();
  int free_variables_for_bonds();
  int free_variables_for_angles();
  int free_variables_for_torsions();
  int free_variables_for_impros();
  int free_variables_for_nb14();
  int free_variables_for_excess();
  int free_variables_for_nb15off();
  
  int set_parameters(int in_n_atomds, PBC* in_pbc, Config* in_cfg,
		     real in_cutoff_pair,
		     int in_n_boxes_x, int in_n_boxes_y, int in_n_boxes_z);
  int set_nsgrid();
  int nsgrid_crd_update();
  int nsgrid_update();
  int nsgrid_update_receiver();
  int rank0_alloc_variables();
  int rank0_free_variables();
  int initial_division(const real** in_crd,
		       const real** in_vel,
		       const real* in_charge,
		       const int* in_atom_type);
  int rank0_div_box(const real** in_crd,
		    const real** in_vel);
		    
  int rank0_send_init_data(const real** in_crd,
			   const real** in_vel,
			   const real* in_charge,
			   const int* in_atom_type);
  int recv_init_data();

  int set_bond_potentials(const int** in_bond_atomid_pairs,
			  const real* in_bond_epsiron,
			  const real* in_bond_r0);
  int set_angle_potentials(const int** in_angle_atomid_triads,
			   const real* in_angle_epsiron,
			   const real* in_angle_theta0);
  int set_torsion_potentials(const int** in_torsion_atomid_quads,
			     const real* in_torsion_energy,
			     const int* in_torsion_overlaps,
			     const int* in_torsion_symmetry,
			     const real* in_torsion_phase,
			     const int* in_torsion_nb14);
  int set_impro_potentials(const int** in_impro_atomid_quads,
			   const real* in_impro_energy,
			   const int* in_impro_overlaps,
			   const int* in_impro_symmetry,
			   const real* in_impro_phase,
			   const int* in_impro_nb14);
  int set_nb14_potentials(const int** in_nb14_atomid_pairs,
			  const int** in_nb14_atomtype_pairs,
			  const real* in_nb14_coeff_vdw,
			  const real* in_nb14_coeff_ele);
  int set_ele_excess(const int** in_excess_pairs);
  int set_nb15off(const int* in_nb15off);
  int set_lj_param(const int in_n_lj_types,
		   real_pw* in_lj_6term,
		   real_pw* in_lj_12term);
  int calc_energy();
  int calc_energy_pairwise();
  int calc_energy_pairwise_wo_neighborsearch();
  bool check_nb15off(const int& a1, const int& a2, const int* bitmask);
  int calc_energy_bonds();
  int calc_energy_angles();
  int calc_energy_torsions();
  int calc_energy_impros();
  int calc_energy_14nb();
  int calc_energy_ele_excess();
  int get_box_id_from_crd(const int box_crd[]);
  int get_box_crd_from_id(const int box_id, 
			  int *box_crd);
  int init_energy();
  int init_work();
  int add_work_from_minicell();
  
  real* get_box_l(){return box_l;};
  int* get_n_boxes_xyz(){return n_boxes_xyz;};
  const real* get_crds(){return (const real*)crd;};
  const int* get_atomids(){return (const int*)atomids;};
  int get_n_atoms_box(){return (const int)n_atoms_box;};

  int swap_velocity_buffer();
  int update_velocities(const real firstcoeff,
			const real time_step,
			const real* mass);
  int velocity_average();
  int revise_coordinates_pbc();
  int set_vel_just(real** vel);
  real cal_kinetic_energy();
  int  update_coordinates(const real time_step);
  bool is_in_box(real* in_crd);
  bool is_in_exbox(real* in_crd);
  int set_box_crd();

  real_fc get_pote_vdw(){return pote_vdw;};
  real_fc get_pote_ele(){return pote_ele;};
  real_fc get_pote_bond(){return pote_bond;};
  real_fc get_pote_angle(){return pote_angle;};
  real_fc get_pote_torsion(){return pote_torsion;};
  real_fc get_pote_impro(){return pote_impro;};
  real_fc get_pote_14vdw(){return pote_14vdw;};
  real_fc get_pote_14ele(){return pote_14ele;};
  
#ifdef F_CUDA
  int gpu_device_setup();
  int update_device_cell_info();
  int calc_energy_pairwise_cuda();
#endif

  //int set_box_region_info(const real** in_crd);  
  //int set_max_n_atoms_region();
  //int get_region_id_from_crd(int width, int rx, int ry, int rz);
  //int get_region_crd_from_id(int width, int regid,
  //int& rx, int& ry, int& rz);
  //int get_box_id_from_crd(int bx, int by, int bz);
  //int get_box_crd_from_id(int boxid, int& bx, int& by, int&bz);
  //int get_global_region_id_from_box_region(int boxid, int regionid);
};

#endif
