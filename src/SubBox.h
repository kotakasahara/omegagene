#ifndef __SUB_BOX_H__
#define __SUB_BOX_H__

#include "CelesteObject.h"
#include "MiniCell.h"
#include "PBC.h"
#include "Config.h"
#include "ForceField.h"
#include "ConstraintShake.h"
#include "Extend.h"
#include "Thermostat.h"
#include "COMMotion.h"
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

  real time_step;
  real time_step_inv;
  real time_step_inv_sq;
  real temperature_coef;

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
  real *crd_prev;
  real *vel;
  real *vel_next;
  real *vel_just;
  real_fc *work;
  real_pw *charge;
  real_pw *mass;
  real_pw *mass_inv;
  int  *atom_type;

  // buffer for thermostat with shake 
  real *buf_crd;
  //

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
  //int* bp_bonds;
  //int* bp_angles;
  //int* bp_torsions;
  //int* bp_impros;
  //int* bp_nb14;
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

  //int constraint_type;
  int  flg_constraint;
  int  flg_settle;
  ConstraintObject* constraint;
  ConstraintObject* settle;
  ExtendedVMcMD* extended;
  
  int flg_thermostat;
  ThermostatObject* thermostat;

  bool flg_mod_15mask;
  COMMotion commotion;
  
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
  int nsgrid_crd_to_gpu();
  int nsgrid_update();
  int nsgrid_update_receiver();
  int rank0_alloc_variables();
  int rank0_free_variables();
  int initial_division(real** in_crd,
		       real** in_vel,
		       real_pw* in_charge,
		       real_pw* in_mass,
		       int* in_atom_type);
  int rank0_div_box(real** in_crd,
		    real** in_vel);
		    
  int rank0_send_init_data(real** in_crd,
			   real** in_vel,
			   real_pw* in_charge,
			   real_pw* in_mass,
			   int* in_atom_type);
  int recv_init_data();

  int set_bond_potentials(int** in_bond_atomid_pairs,
			  real* in_bond_epsiron,
			  real* in_bond_r0);
  int set_angle_potentials(int** in_angle_atomid_triads,
			   real* in_angle_epsiron,
			   real* in_angle_theta0);
  int set_torsion_potentials(int** in_torsion_atomid_quads,
			     real* in_torsion_energy,
			     int* in_torsion_overlaps,
			     int* in_torsion_symmetry,
			     real* in_torsion_phase,
			     int* in_torsion_nb14);
  int set_impro_potentials(int** in_impro_atomid_quads,
			   real* in_impro_energy,
			   int* in_impro_overlaps,
			   int* in_impro_symmetry,
			   real* in_impro_phase,
			   int* in_impro_nb14);
  int set_nb14_potentials(int** in_nb14_atomid_pairs,
			  int** in_nb14_atomtype_pairs,
			  real* in_nb14_coeff_vdw,
			  real* in_nb14_coeff_ele);
  int set_ele_excess(int** in_excess_pairs);
  int set_nb15off(int* in_nb15off);
  int set_lj_param(const int in_n_lj_types,
		   real_pw* in_lj_6term,
		   real_pw* in_lj_12term);
  int calc_energy();
  int calc_energy_pairwise();
  int calc_energy_pairwise_wo_neighborsearch();
  bool check_nb15off(const int& a1, const int& a2, const int* bitmask,
		     int& mask_id, int& interact_bit);
  int calc_energy_bonds();
  int calc_energy_angles();
  int calc_energy_torsions();
  int calc_energy_impros();
  int calc_energy_14nb();
  int calc_energy_ele_excess();
  real_pw* get_mass(){return mass;};
  int get_box_id_from_crd(const int box_crd[]);
  int get_box_crd_from_id(const int box_id, 
			  int *box_crd);
  int init_energy();
  int init_work();
  int add_work_from_minicell();
  
  real* get_box_l(){return box_l;};
  int* get_n_boxes_xyz(){return n_boxes_xyz;};
  real* get_crds(){return crd;};
  int* get_atomids(){return atomids;};
  int get_n_atoms_box(){return n_atoms_box;};

  int cpy_crd_prev();
  int cpy_crd_from_prev();
  int cpy_vel_buf_from_prev();
  int cpy_vel_prev_from_buf();
  int cpy_vel_prev();
  int swap_velocity_buffer();
  int update_velocities(const real time_step);
  int velocity_average();
  int set_velocity_from_crd();
  //int set_force_from_velocity(const real time_step);
  int revise_coordinates_pbc();
  int copy_crd(real** p_crd);
  int copy_vel_just(real** p_vel);
  int copy_vel(real** p_vel);
  int copy_vel_next(real** p_vel);
  int copy_crdvel_to_mmsys(real* src, real** dst);
  int add_force_from_mmsys(real_fc** in_force);
  template <typename TYPE> int add_crdvel_from_mmsys(TYPE* dst, TYPE** src);
  real cal_kinetic_energy();
  int update_coordinates_cur(const real time_step);
  int update_coordinates_prev(const real time_step);
  int  update_coordinates(const real time_step, real *p_crd, real* p_vel);
  int update_coordinates_nsgrid();
  bool is_in_box(real* in_crd);
  bool is_in_exbox(real* in_crd);
  int set_box_crd();
  
  int init_constraint(int in_constraint,
		      int in_max_loops,
		      real in_tolerance,
		      int max_n_pair,
		      int max_n_trio,
		      int max_n_quad,
		      int max_n_settle);
  int set_subset_constraint(ConstraintObject& in_cst,
			    ConstraintObject& in_settle);
  int init_thermostat(const int in_thermostat_type,
		      const real in_temperature,
		      const int d_free);
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
  
  int apply_constraint();
  int update_thermostat(const int cur_step);
  int apply_thermostat();
  int apply_thermostat_with_shake(const int max_loop,
				  const real tolerance);
  void set_extended(ExtendedVMcMD* in_exp){ extended =in_exp; };
  int extended_apply_bias(unsigned long cur_step, real in_lambda);
  int extended_apply_bias_struct_param(unsigned long cur_step);
  void extended_enable_vs_transition();
  int cancel_com_motion();

  int set_com_motion(int n_groups, int* group_ids,
			   int*  n_atoms_in_groups, 
			   int** groups,
			   real* mass_inv_groups);

  int print_work(int atom_id);

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
