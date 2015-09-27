#ifndef __MM_SYSTEM_H__
#define __MM_SYSTEM_H__

#include "define.h"
#include <iostream>
#include <set>
#include <ctime>

//#include <random>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
using namespace std;

#include "CelesteObject.h"
#include "Config.h"
#include "PBC.h"
#include "ForceField.h"
#include "Constraint.h"
#include "Extend.h"
#include "DistRestraint.h"

class MmSystem : public CelesteObject{
 private:
 public:
  // topology information
  // nonbonds
  int n_lj_types;
  int n_lj_type_pairs;
  real_pw* lj_6term;
  real_pw* lj_12term;
  // bonds
  int n_bonds;
  int** bond_atomid_pairs;
  real* bond_epsiron;
  real* bond_r0;
  // angle
  int n_angles;
  int** angle_atomid_triads;
  real* angle_epsiron;
  real* angle_theta0;
  // torsions
  int n_torsions;
  int** torsion_atomid_quads;
  real* torsion_energy;
  int*  torsion_overlaps;
  int*  torsion_symmetry;
  real* torsion_phase;
  int*  torsion_nb14;
  // impropers
  int n_impros;
  int** impro_atomid_quads;
  real* impro_energy;
  int*  impro_overlaps;
  int*  impro_symmetry;
  real* impro_phase;
  int*  impro_nb14;
  // nonbonds 1-4
  int n_nb14;
  int** nb14_atomid_pairs;
  int** nb14_atomtype_pairs;
  real* nb14_coeff_vdw;
  real* nb14_coeff_ele;

  // 1-5 pairs excluded from 1-5 interactions
  // nb15off_atomid_pairs[atom_id] = [atom_id1, atom_id2, ...]
  // n_nb15off_atom[atom_id] = number of atoms without 1-5 interactions
  int* nb15off1;  
  int* nb15off2;  
  int n_nb15off;
  int* nb15off;
  int max_n_nb15off;

  // excess pairs
  int max_n_excess;
  int n_excess;
  int** excess_pairs;

  string launchset_version;

  real pbc_val[9];
  PBC pbc;

  real leapfrog_coef;
  int n_atoms;
  int d_free;
  real** crd;
  real_fc** force;
  //real** vel;
  real** vel_just;
  //real** vel_next;
  real_pw* charge;
  real_pw* mass;
  int* atom_type;
  // self energy for ZD
  real* energy_self;  
  real energy_self_sum;

  real_fc potential_e;
  real_fc pote_bond;
  real_fc pote_angle;
  real_fc pote_torsion;
  real_fc pote_impro;
  real_fc pote_14vdw;
  real_fc pote_14ele;
  real_fc pote_vdw;
  real_fc pote_ele;

  real_fc pote_dist_rest;

  real kinetic_e;
  real temperature;

  unsigned long cur_step;
  real cur_time;

  ForceField ff;

  ConstraintObject constraint;
  ConstraintObject settle;
  ExtendedVMcMD* vmcmd;

  int n_groups;
  int* n_atoms_in_groups;
  int** atom_groups;
  real_pw* mass_groups;
  real_pw* mass_inv_groups;
  vector<string> atom_group_names;

  int n_com_cancel_groups;
  int com_cancel_groups[MAX_N_COM_GROUPS];
  int n_enhance_groups;
  int enhance_groups[MAX_N_COM_GROUPS];
  int out_group;
  
  DistRestraintObject*  dist_restraint;

  //mt19937 random_mt;

  //time

  clock_t ctime_per_step;
  clock_t ctime_cuda_htod_atomids;
  clock_t ctime_cuda_reset_work_ene;
  clock_t ctime_calc_energy;
  clock_t ctime_calc_energy_pair;
  clock_t ctime_calc_energy_bonded;
  clock_t ctime_update_velo;
  clock_t ctime_calc_kinetic;
  clock_t ctime_update_coord;
  clock_t ctime_setgrid;
  clock_t ctime_enumerate_cellpairs;

  ///// METHODS //////

  MmSystem();
  ~MmSystem();
  int alloc_atom_vars();
  int alloc_lj_params();
  int alloc_bonds();
  int alloc_angles();
  int alloc_torsions();
  int alloc_impros();
  int alloc_nb14();
  int alloc_nb15off();
  int alloc_excess_pairs();
  int alloc_atom_groups(int in_n_groups,
			int* in_n_atoms_in_groups);

  //int alloc_pcluster_vars();

  int free_all();
  int free_atom_vars();
  int free_lj_params();
  int free_bonds();
  int free_angles();
  int free_torsions();
  int free_impros();
  int free_nb14();
  int free_nb15off();
  int free_excess_pairs();
  int free_atom_groups();

  //int free_pcluster_vars();
  // parameter setter
  int set_lj_pair_param(int type1, int type2, real_pw param6, real_pw param12);
  int set_bond_param(int bond_id, 
		     int atomid1, int atomid2,
		     real eps, real r0);
  int set_angle_param(int angle_id, 
		      int atomid1, int atomid2, int atomid3,
		      real eps, real theta0);
  int set_torsion_param(int torsion_id, 
			int atomid1, int atomid2, int atomid3, int atomid4,
			real ene, int overlaps,
			int symmetry, real phase,
			int flag_14nb);
  int set_impro_param(int torsion_id, 
			int atomid1, int atomid2, int atomid3, int atomid4,
			real ene, int overlaps,
			int symmetry, real phase,
			int flag_14nb);
  int set_nb14_param(int nb14_id, 
		     int atomid1, int atomid2,
		     int atomtype1, int atomtype2,
		     real coeff_vdw, real coeff_ele);
  bool search_nb15off(int atomid1, int atomid2);
  int set_nb15off(int atomid1, int atomid2);

  int add_excess_pairs(int atomid1, int atomid2);
  int set_excess_pairs();
  int set_atom_group_info(Config* cfg);
  int get_atom_group_id_from_name(const string name);
  int print_com_cancel_groups();
  int set_com_cancel_groups(Config* cfg);
  int set_enhance_groups(Config* cfg);
  int print_enhance_groups();
  int set_out_group(Config* cfg);
  int print_out_group();
  int get_out_group(){return out_group;};

  // calc
  
  int reset_energy();
  //int velocity_swap();
  //int velocity_average();
  //int update_velocities();
  int write_data();
  int nsgrid_setup(real nsgrid_cutoff); // , int box_div[]);		

  int nsgrid_update();
  int nsgrid_update_receiver();
  int nsgrid_crd_update();
  int revise_crd_inbox();
  int set_random(int seed);
  int output_ctimes();
  int set_vel_from_box();

  int ff_setup(const Config* cfg);

  real_fc set_potential_e();
};

#endif

