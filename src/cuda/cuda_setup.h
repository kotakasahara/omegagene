#ifndef __CUDA_SETUP_H__
#define __CUDA_SETUP_H__

#include <stdio.h>
#include <stdlib.h>
#include "cuda_common.h"
#include "../CelesteObject.h"
#include "../MiniCell.h"

#define real4 float4
#define real_fc3 double3

#ifdef  THREADS128
  #define PW_THREADS 128
#elif   THREADS256
  #define PW_THREADS 256
#elif   THREADS512
  #define PW_THREADS 512
#else
  #define PW_THREADS 128
#endif

#define PW_THREADS_DIAG PW_THREADS
#define REORDER_THREADS 512

#define N_MULTI_WORK 8
#define MAX_INT 99999999

cudaStream_t stream_pair_home;

#define D_N_ATOM_CELL  8
#define D_N_ATOM_CELL_3  24

#define WARPSIZE 32

__constant__ real_pw  D_CHARGE_COEFF;
__constant__ real_pw  PBC_L[3];
__constant__ real_pw  PBC_LOWER_BOUND[3];
__constant__ int      D_N_CELL_PAIRS;
__constant__ int      D_N_CELLS;
__constant__ int      D_N_CELLS_X;
__constant__ int      D_N_CELLS_Y;
__constant__ int      D_N_COLUMNS;
__constant__ int      D_N_NEIGHBOR_COL_X;
__constant__ int      D_N_NEIGHBOR_COL_Y;
__constant__ int      D_N_ATOMS;
__constant__ int      D_N_ATOM_ARRAY;
__constant__ int      D_N_ATOMTYPES;
__constant__ int      D_N_UNI;
__constant__ int      D_N_UNI_Z;
__constant__ real_pw  D_L_CELL_X;
__constant__ real_pw  D_L_CELL_Y;
__constant__ real_pw  D_L_UNI_Z;
//__constant__ int      D_MAX_N_NB15OFF;
//__constant__ int      D_MAX_N_ATOMS_GRID;

__constant__ real_pw  D_CUTOFF;
__constant__ real_pw  D_CUTOFF_PAIRLIST;

__constant__ real_pw  D_ZCORE;
__constant__ real_pw  D_BCOEFF;
__constant__ real_pw  D_FCOEFF;

// x,y,z: Cartesian coordinate,
// w: charge
real4* d_crd_chg;

real_pw* d_crd;

CellPair* d_cell_pairs;
int* d_idx_head_cell_pairs;
int* d_cell_pair_removed;

// atomid ordered by atomid_grid order
//  integer values in this array are original atomid
//  index of this array is atomid_grid index
//   d_atomids[d_atomid] = atomid
//   d_atomids_rev[atomid] = d_atomid
int* d_atomids;
int* d_atomids_rev;

// info reordered by atomid_grid
// x: atomid in original order
// y: atomtype
int* d_atomtype;

// info in original atomid order
real_pw* d_charge_orig;
int* d_atomtype_orig;

//int* d_nb15off1_orig;
//int* d_nb15off2_orig;

real_pw* d_lj_6term;
real_pw* d_lj_12term;

real_fc* d_energy;
real_fc* d_work;

int*  d_idx_xy_head_cell;
int2* d_uni2cell_z;

//texture<real, 2> tex_lj_6term;
//texture<real, 2> tex_lj_12term;
//texture<int, 2> tex_nb15off;
//texture<int> tex_n_nb15off;

//float h_test_buf[2];
//float *d_test_buf;

//test

#endif
