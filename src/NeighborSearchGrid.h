#ifndef __NEIGHBOR_SEARCH_GRID_H__
#define __NEIGHBOR_SEARCH_GRID_H__

#include <algorithm>
//#include <random>
#include <cstdio>
#include "CelesteObject.h"
#include "PBC.h"

using namespace std;

typedef struct grid_pair_info{
  int grid_id1;
  int grid_id2;
  real cx;
  real cy;
  real cz;
  int n_pairs;
  //bool flg_all;
} GridPair;

//typedef struct atom_info{
//  real_pw x;
//  real_pw y;
//  real_pw z;
//  real_pw charge;
//  int atom_type;
//} AtomInfo;

class NeighborSearchGrid : public CelesteObject {
 private:
 public:

  // the maximum number of atoms in each grid cell
  int max_n_atoms;
  
  // n_ng: number of neighbor grid for each axis.
  // when n_ng=2 ... grids in gx-2, gx-1, gx+1, gx+2
  // took into account
  int min_grid_width;

  // previous n grids in each axis
  // were searched from each grid
  int n_neighbors_axis[3];

  // grid_pairs:
  int max_n_grid_pairs;
  int n_grid_pairs;
  // grid_pairs[i] = [grid_id1, grid_id2, cx, cy, cz, n_atom_pairs]
  //  cx, cy, dz ... -1, 0, +1 
  //    that means straddling of PBC boundaries of grid 2
  GridPair* grid_pairs;
  //  GridPair* grid_pairs_diag;

  int n_atoms;

  // Cutoff length to make neighbor paris .
  // It needs offset to include all of particle clusters
  //   in cutoff.
  // (the offset has to be set as diameter of the large cluster)
  real cutoff_pair;

  // Atomic coordinates for each grid;
  // crd[3 * atom_id in grid array] = 
  // atom_id in grid array ... grid_atom_id[grid_id] ~ grid_n_atoms[grid_id]
  // AtomInfo* atominfo;
  real_pw *crd;
  // atom ids in order of *crd
  int *atomids;

  // page locked
  //  for return from cuda_pairwise_ljzd()
  real_fc *work;
  real_pw *energy;

  // internal coordinates in the neighbor search grids
  //   crd_grid[atom_id] = [x, y, z]
  real** crd_grid;

  // assignment to ns grid cells;
  //   atom_id => grid cell
  //   atom_grid[atom_id] = [gird_x, grid_y, grid_z]
  int** atom_grid;

  //   grid cell => atom
  //   grid_atom[grid_id] = [atom_id, ... ]
  int** grid_atom;

  // grid_atom_id[grid_id]
  //  index of the first atom of each grid in 
  //  the atominfo array
  int* grid_atom_index;

  //   grid cell => number of clusters
  //   grid_n_atom[grid_id] = number of atoms
  int* grid_n_atoms;

  //  n_grid = [number of cells_x, _y, _z]
  int n_grid[3];
  //  n_all_grids = n_grid[0] * n_grid[1] * n_grid[2];
  int n_all_grids;

  //   ns_length_grid = [cell_lenght_x, _y, _z]
  real L_grid[3];

  const PBC* pbc;

  // METHODS
  NeighborSearchGrid();
  ~NeighborSearchGrid();
  inline int get_grid_id(int x, int y, int z){ return z*n_grid[0]*n_grid[1] + y*n_grid[0] + x; }
  inline void get_grid_xyz(int grid_id, int& x, int& y, int& z){
    int xy = grid_id % (n_grid[0] * n_grid[1]);
    x = xy % n_grid[0];
    y = xy / n_grid[0];
    z = grid_id / (n_grid[0] * n_grid[1]);
  }

  int set_grid(int in_n_atoms, 
	       real in_cutoff_pair, 
	       real in_min_grid_width,
	       int in_max_n_atoms,
	       const PBC* in_pbc);
  int alloc_variables();
  int free_variables();

  // setup_crd_into_grid called by MmSystem::setup_nsgrid
  //  only at once at begining of simulation
  int setup_crd_into_grid(real** in_crd, real_pw* in_charge, int* in_atom_type);
  int print_max_n_atoms_in_grid();
  // functions called every steps
  int set_crd_grid(real** in_crd, real_pw* in_charge, int* in_atom_type);
  int set_grid_atom_index_from_grid_n_atoms();
  int reflesh_atom_assign(real** in_crd, real_pw* in_charge, int* in_atom_type);
  int enumerate_grid_pairs();
  int get_grid_positions(real positions[8][3],int gx, int gy, int gz);
  int add_grid_pairs(int grid_id1, int grid_id2, int image[3]);
  int check_n_atoms_in_grid();
  bool check_grid_pairs_exist(int grid_id1, int grid_id2);
  //int reorder_grid_pairs_random(mt19937 random_mt);
  int reorder_grid_pairs_random();
};

#endif
