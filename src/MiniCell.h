#ifndef __MINI_CELL_H__
#define __MINI_CELL_H__

#include "CelesteObject.h"
#include "PBC.h"
#include <algorithm>
#include <cstdio>
#include <cstring>
using namespace std;

#define N_ATOM_CELL 8
#define N_BITMASK (N_ATOM_CELL * N_ATOM_CELL + 31) / 32
#define COEF_MAX_N_ATOMS_CELL 1.2
#define MAX_N_CELL_UNI 10
#define MAX_N_PAIR_Z 100
#define COEF_MAX_N_CELL_PAIRS 1.3

typedef struct cell_pair_info{
  int cell_id1;
  int cell_id2;
  int image;
  int pair_mask[N_BITMASK];
  //real cx;
  //real cy;
  //real cz;

} CellPair;

class MiniCell : public CelesteObject{
 private:
  // n_atoms
  //  the number of atoms in the whole system
  int n_atoms;

  int n_cells;
  int n_columns;
  int n_cells_xyz[3];
  int *n_cells_z;

  int max_n_cells;
  // L_cell_x, L_cell_y
  real_pw L_cell_xyz[3];

  // n_neighbors_x, _y
  //   the number of neighboring cells in X and Y axes
  //   "neighboring" means it can be included in 
  //   the cutoff sphere
  int n_neighbors_xyz[3];
  
  // cutoff_pair
  //   cutoff length for pair interactions
  real_pw cutoff_pair_half;
  real_pw cutoff_pair;
  real_pw cutoff_pair2;


  // cell_crd
  //  cell_crd[cell_id] = cell_id_x, cell_id_y, cell_id_z
  int **cell_crd;

  //  <->
  // idx_crd_cell
  int ***idx_crd_cell;

  // idx_atom_cell
  //  idx_atom_cell[atom_id_grid] = cell_id  
  int *idx_atom_cell;
  //  idx_cell_head_atom[cell_id] = atom_id  _grid
  //int *idx_cell_head_atom;
  int *idx_cell_n_atoms;
  
  // flg_dummy[atom_id_grid] = 1 ... dummy, 0 ... not dummy
  int *flg_dummy;

  // idx_atom_cell_xy
  //  idx_atom_cell_xy[atom_id (original)] = cell_id_xy
  //    id for xy cloumn
  int *idx_atom_cell_xy;

  // idx_head_cell_xy
  //  idx_xy_head_cell[column_id] = cell_id
  //    the first cell in the column_id
  int *idx_xy_head_cell;
  int *idx_xy_head_atom;

  // number of atoms in each xy column
  //   n_atoms_xy[cell xy id] = N
  int *n_atoms_xy;


  // the maximum number of atos for each column
  int max_n_atoms_column;
  
  PBC* pbc;  

  ///// About Box division /////

  //  the number of atoms in the box
  int n_atoms_box;
  //  the number of atoms in the box including replica
  int n_atoms_exbox;
  //  the maximum number of atoms for allocation
  int max_n_atoms_box;
  int max_n_atoms_exbox;
  int max_n_atoms_region[125];
  
  int  box_id;
  int  box_crd[3];
  int  n_boxes;
  int  n_boxes_xyz[3];
  real_pw box_lower_bound[3];
  real_pw box_upper_bound[3];
  real_pw box_l[3];
  real_pw exbox_lower_bound[3];
  real_pw exbox_upper_bound[3];
  real_pw exbox_l[3];

  int n_atoms_region[125];
  int** region_atoms; //atomid_g

  int* nb15off;
  int max_n_nb15off;
  
  // for mpi
  int mpi_n_atoms;
  real_pw* mpi_sendbuf_crd;
  real_pw* mpi_recvbuf_crd;
  real_fc* mpi_sendbuf_frc;
  real_fc* mpi_recvbuf_frc;
  int* mpi_sendbuf_atomids;
  int* mpi_recvbuf_atomids;
  
  // atominfo
  real_pw **crd_in_cell;

  int *atomids_rev;
  int *atomids;

  // atomids_buf
  //   2. Temporary buffer for atomids 
  //      (atomids for space decompositions)
  //        set_crds_to_homebox
  int *atomids_buf;

  real_pw *crd;
  real_fc *work;
  real_fc *energy;

  int max_n_cell_pairs;
  int n_cell_pairs;
  CellPair* cell_pairs;  

  // idx_head_cell_pairs[cell_id] = 
  //   index in the array cell_pairs
  //   for the first element of 
  //   cell_pairs[ ].cell1_id = cell_id
  int *idx_head_cell_pairs;

  // Uniform z-grid
  //   index for enumerating pairs of 8-atom grid cells
  //   z crd divided into PBC.L[2]/(cutoff/2)
  //  uni_id = 
  //   uni_id(z) = uni_id / ( n_uni(x) * n_uni(y) )
  //   uni_id(y) = uni_id % ( n_uni(x) * n_uni(y) )
  //   uni_id(x) = uni_id % ( n_uni(x) * n_uni(y) )
  //  uni2cell[uni_id] = [cell1, cell2, ...]
  //  cell2uni[cell_id] = [uni_z1, uni_z2, ...]
  int n_uni;
  real_pw L_z_uni;
  int n_uni_z;
  int* n_uni2cell_z;
  int** uni2cell_z;
  int** cell2uni_z;
  int get_uni_z(int uni_id);
  int get_uni_id_from_crd(int x, int y, int z);

 public:
  MiniCell();
  ~MiniCell();
  int alloc_variables();
  int free_variables();
  int init_variables();
  int init_energy_work();
  // set_grid()
  //   reset grid mini-cells
  int set_grid_xy();
  int set_atoms_into_grid_xy();
  int reorder_atominfo_for_columns();
  int set_dummy_atoms();
  int reorder_atominfo_in_columns();
  int reset_cell_assignment();
  int swap_atomid_crd(const int i, const int j);
  int quick_sort_in_columns(const int l, const int r);
  int debug_set_atoms_into_grid();
  //int update_cell_assign(const real** in_crd);
  int update_crd(real** in_crd);
  int set_idx_xy_head_atom_from_n_atoms_xy();
  int set_atomids_rev();

  int set_atomids_buf();
  int set_atomids_from_buf();

  // pair
  real_pw get_cell_z_min(int cell_id);
  real_pw get_cell_z_max(int cell_id);
  const int* get_n_neighbors_xyz(){ return (const int*)n_neighbors_xyz;};
  int get_n_neighbor_cols(){ return (n_neighbors_xyz[0]*2+1) * (n_neighbors_xyz[1]*2+1);};
  int set_uniform_grid();
  int enumerate_cell_pairs();
  bool check_valid_pair(const int cell1_id, const int cell2_id,
			const bool cell1_odd, const bool cell2_odd);
  int add_cell_pair(const int cell_id1, const int cell_id2, const int image[3]);
  int set_cell_pair_bitmask(const int cell_id1, const int cell_id2, int* bitmask);

  // getter
  const real_pw* get_L_cell_xyz() {return L_cell_xyz;};
  int get_n_cell_pairs(){return n_cell_pairs; };
  const CellPair get_cell_pair(const int cpid){return cell_pairs[cpid]; };
  CellPair*& get_cell_pairs(){return cell_pairs; };
  int*& get_idx_head_cell_pairs(){return idx_head_cell_pairs;};
  real_pw*& get_crd(){return crd;};
  void get_crd(int atomid_grid, real_pw& x, real_pw& y, real_pw& z);
  real_fc*& get_work(){return work;};
  real_fc*& get_energy(){return energy;};
  int* get_idx_atom_cell_xy(){return idx_atom_cell_xy;};  
  int get_idx_cell_head_atom(const int cid){return cid*N_ATOM_CELL; }
  int get_n_atoms_in_cell(const int cid){
    int col = get_column_id_from_crd(cell_crd[cid][0], cell_crd[cid][1]);
    if(cell_crd[cid][2] < n_cells_z[col]) return N_ATOM_CELL;
    else return n_atoms_xy[col]%N_ATOM_CELL;
  }
  int*& get_atomids(){return atomids;}
  int*& get_atomids_rev(){return atomids_rev;}
  int get_atomid_from_gridorder(const int aid){return atomids[aid];}
  //  void get_crd_from_gridorder(const int aid, real_pw* ret_crd){
  //memcpy(ret_crd, crd+aid*3, sizeof(real_pw)*3);};
  int get_column_id_from_crd(int x, int y);
  //int get_column_crd_from_id(int x, int y);
  int get_cell_id_from_crd(int x, int y, int z);
  
  int is_dummy(int atom_id_grid){ return atomids[atom_id_grid]==-1; };

  int get_n_atom_array(){ return n_cells*N_ATOM_CELL; }
  int get_max_n_atom_array(){ return max_n_atoms_exbox+n_columns*N_ATOM_CELL; };
  int get_n_cells(){ return n_cells; };
  int get_max_n_cells(){ return max_n_cells; };
  int get_max_n_cell_pairs(){ return max_n_cell_pairs; };

  const int* get_n_cells_xyz(){ return (const int*)n_cells_xyz; };
  int get_n_cells_x(){ return n_cells_xyz[0]; };
  int get_n_cells_y(){ return n_cells_xyz[1]; };
  int get_n_columns(){ return n_columns; };

  //int get_ene_forces(real_fc& in_ene_vdw, real_fc& in_ene_ele, real_fc**& in_force);

  /////// setter
  // set_grid_parameters()
  //  executed only at once, at before initial step
  int set_grid_parameters(const int in_n_atoms,
			  const real in_cutoff_pair,
			  const PBC* in_pbc,
			  const int in_max_n_nb15off,
			  int* in_nb15off);
  
  //real move_crd_in_cell(const int atomid, const int dim, const real val);
  void add_energy(const real_fc in_vdw, const real_fc in_ele){ energy[0] += in_vdw; energy[1] += in_ele; };
  int add_work(const int atomid_grid,
	       const real_fc in_w1, const real_fc in_w2, const real_fc in_w3);
  int move_atom(const int& atomid, const int& d, const real& diff);
  //int update_coordinates(const real** vel_next, const real& time_step);

  int set_box_info(int* in_n_boxes_xyz, real* in_box_l);
  int print_box_info();
  int set_crds_to_homebox(real* in_crd,
			  int* in_atomids,
			  int in_n_atoms_box);
  int set_max_n_atoms_region();
  int get_region_id_from_crd(int width, int rx, int ry, int rz);
  int get_region_crd_from_id(int width, int regid,
			 int& rx, int& ry, int& rz);
  int get_box_id_from_crd(const int box_crd[]);
  int get_box_crd_from_id(const int& box_id, int* box_crd);
  int assign_regions(int from_atomid_g, int to_atomid_g);
  int add_atom_to_region(int atomid_g, int region_id);
  int setup_replica_regions();
  int mpi_send_crd_replica_regions(int dir_x, int dir_y, int dir_z);

  int get_n_uni(){ return n_uni; };
  int get_n_uni_z(){ return n_uni_z; };
  real_pw get_l_uni_z(){ return L_z_uni; };
  int*& get_idx_xy_head_cell(){return idx_xy_head_cell;};
};

#endif 

