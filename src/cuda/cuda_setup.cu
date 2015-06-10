#include "cuda_setup.h"

__device__   __inline__    double  shfl_xor ( double value,  int const lane, int const warpsize ){
  return  __hiloint2double( __shfl_xor(__double2hiint(value),lane, warpsize),
			    __shfl_xor(__double2loint(value),lane, warpsize)); 
}

__device__ double atomicAdd(double* address, double val){
  unsigned long long int* address_as_ull =
    (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
		    __double_as_longlong(val +
					 __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}


extern "C" int cuda_alloc_atom_info(int n_atoms,
				    int n_atom_array,
				    int max_n_cells,
				    int max_n_cell_pairs,
				    int n_columns,
				    int n_uni){
  printf("cuda_alloc_atom_info\n");
  // max_n_atoms ... maximum number of atoms for each grid cell
  HANDLE_ERROR( cudaMalloc((void**)&d_crd_chg,
			   n_atom_array * sizeof(real4)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_crd,
			   n_atom_array * 3 * sizeof(real_pw)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_charge_orig,
			   n_atoms * sizeof(real_pw) ));
  HANDLE_ERROR( cudaMalloc((void**)&d_atomtype,
			   n_atom_array * sizeof(int) ));
  HANDLE_ERROR( cudaMalloc((void**)&d_atomids,
			   n_atom_array * sizeof(int)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_atomids_rev,
			   n_atoms * sizeof(int)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_atomtype_orig,
			   n_atoms * sizeof(int) ));
  HANDLE_ERROR( cudaMalloc((void**)&d_cell_pairs,
			   max_n_cell_pairs * sizeof(CellPair)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_cell_pairs_buf,
			   max_n_cell_pairs * sizeof(CellPair)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_idx_head_cell_pairs,
			   (max_n_cells+1) * sizeof(int)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_cell_pair_removed,
			   (max_n_cells+1) * sizeof(int)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_n_cell_pairs,
			    (max_n_cells) * sizeof(int)) );
  //HANDLE_ERROR( cudaMalloc((void**)&d_grid_atom_index,
  //(max_n_ + 1) * sizeof(int)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_energy,
			   N_MULTI_WORK * 2 * sizeof(real_fc) ) );
  HANDLE_ERROR( cudaMalloc((void**)&d_work,
			   N_MULTI_WORK * n_atom_array * 3 * sizeof(real_fc) ) );
  HANDLE_ERROR( cudaMalloc((void**)&d_idx_xy_head_cell,
			   n_columns * sizeof(int) ) );
  HANDLE_ERROR( cudaMalloc((void**)&d_uni2cell_z,
			   n_uni * sizeof(int2) ) );

  HANDLE_ERROR( cudaMemcpyToSymbol(D_MAX_N_CELL_PAIRS,
				   &max_n_cell_pairs,
				   sizeof(int) ) );  

  //HANDLE_ERROR( cudaMalloc((void**)&d_work_orig,
  //n_atom_array * 3 * sizeof(real_fc) ) );
  return 0;
}

extern "C" int cuda_free_atom_info(){
  //printf("cuda_free_device_atom_info\n");
  HANDLE_ERROR( cudaFree(d_crd_chg) );
  HANDLE_ERROR( cudaFree(d_crd) );
  HANDLE_ERROR( cudaFree(d_atomids) );
  HANDLE_ERROR( cudaFree(d_atomids_rev) );
  HANDLE_ERROR( cudaFree(d_charge_orig) );
  HANDLE_ERROR( cudaFree(d_atomtype) );
  HANDLE_ERROR( cudaFree(d_atomtype_orig) );
  HANDLE_ERROR( cudaFree(d_cell_pairs) );
  HANDLE_ERROR( cudaFree(d_cell_pairs_buf) );
  HANDLE_ERROR( cudaFree(d_idx_head_cell_pairs) );
  HANDLE_ERROR( cudaFree(d_cell_pair_removed) );
  HANDLE_ERROR( cudaFree(d_n_cell_pairs) );  
  HANDLE_ERROR( cudaFree(d_energy) );
  HANDLE_ERROR( cudaFree(d_work) );
  HANDLE_ERROR( cudaFree(d_idx_xy_head_cell));
  HANDLE_ERROR( cudaFree(d_uni2cell_z));
  //HANDLE_ERROR( cudaFree(d_work_orig) );
  return 0;
}

//extern "C" int cuda_memcpy_htod_cell_pairs(CellPair*& h_cell_pairs,
////int*& h_idx_head_cell_pairs,
//int n_cell_pairs,
//int n_cells){
  //printf("cuda_memcpy_htod_cell_pairs\n");
  //HANDLE_ERROR(cudaMemcpy(d_cell_pairs, h_cell_pairs,
  //n_cell_pairs * sizeof(CellPair),
  //cudaMemcpyHostToDevice));
  //HANDLE_ERROR(cudaMemcpy(d_idx_head_cell_pairs,
  //h_idx_head_cell_pairs,
  //(n_cells+1) * sizeof(int),
  //cudaMemcpyHostToDevice));
//  return 0;
//}
extern "C" int cuda_memcpy_htod_atomids(int*& h_atomids,
					int*& h_idx_xy_head_cell,
					int n_atom_array,
					int n_columns){
  HANDLE_ERROR(cudaMemcpy(d_atomids, h_atomids,
			  n_atom_array * sizeof(int),
			  cudaMemcpyHostToDevice));
  HANDLE_ERROR(cudaMemcpy(d_idx_xy_head_cell, h_idx_xy_head_cell,
  			  n_columns * sizeof(int),
  			  cudaMemcpyHostToDevice));
  return 0;
}

// cuda_memcpy_htod_atom_info
//   Arrays of charges and atomtypes of all atoms in the process are sent to
//   the device.
extern "C" int cuda_memcpy_htod_atom_info(real_pw*& h_charge_orig,
					  int*& h_atomtype_orig,
					  int n_atoms){
  HANDLE_ERROR(cudaMemcpy(d_charge_orig, h_charge_orig,
			  n_atoms * sizeof(real_pw),
			  cudaMemcpyHostToDevice));
  HANDLE_ERROR(cudaMemcpy(d_atomtype_orig, h_atomtype_orig,
			  n_atoms * sizeof(int),
			  cudaMemcpyHostToDevice));
  return 0;
}

// cuda_memcpy_htod_crd
//  Sending nsgrid.crd to device
extern "C" int cuda_memcpy_htod_crd(real_pw*& h_crd,
				    int n_atom_array){
  HANDLE_ERROR(cudaMemcpy(d_crd, h_crd,
			  n_atom_array * 3 * sizeof(real_pw),
			  cudaMemcpyHostToDevice));
  //HANDLE_ERROR(cudaMemcpy(d_grid_atom_index, h_grid_atom_index,
  //(n_grids+1) * sizeof(int),
  //cudaMemcpyHostToDevice));
  //printf("cuda_memcpy_htod_atom_info(AtomInfo *d_atominfo\n");
  //n_grids * sizeof(int),
  //cudaMemcpyHostToDevice));
  HANDLE_ERROR( cudaMemset(d_energy, 0.0, sizeof(real_fc)*2*N_MULTI_WORK ));
  HANDLE_ERROR( cudaMemset(d_work, 0.0, sizeof(real_fc)*n_atom_array*3*N_MULTI_WORK ));
  
  return 0;
}

extern "C" int cuda_set_pbc(real_pw* l, real_pw* lb){
  HANDLE_ERROR( cudaMemcpyToSymbol(PBC_L,
				   l, sizeof(real_pw) * 3) );
  HANDLE_ERROR( cudaMemcpyToSymbol(PBC_LOWER_BOUND,
				   lb, sizeof(real_pw) * 3) );
  //printf("dbg cuda l %f %f %f , %f %f %f \n",PBC_L[0], PBC_L[1], PBC_L[2],
  //l[0], l[1], l[2]);
  return 0;
}

extern "C" int cuda_zerodipole_constant(real_pw zcore,
					real_pw bcoeff,
					real_pw fcoeff){
  HANDLE_ERROR( cudaMemcpyToSymbol(D_ZCORE,
				   &zcore,
				   sizeof(real_pw) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_BCOEFF,
				   &bcoeff,
				   sizeof(real_pw) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_FCOEFF,
				   &fcoeff,
				   sizeof(real_pw) ) );
  return 0;
}

// cuda_set_cell_constant
//  These constants are updated when the cell grid is updated
extern "C" int cuda_set_cell_constant(const int  n_cells,
				      const int max_n_cell_pairs,
				      const int  n_cell_pairs,
				      const int  n_atom_array,
				      const int  n_cells_x,
				      const int  n_cells_y,
				      const int  n_columns,
				      const int  n_uni,
				      const int  n_uni_z,
				      const real_pw l_cell_x,
				      const real_pw l_cell_y,
				      const real_pw L_uni_z,
				      const int n_neighbor_col_x,
				      const int n_neighbor_col_y){
  HANDLE_ERROR( cudaMemcpyToSymbol(D_N_CELL_PAIRS,
				   &n_cell_pairs,
				   sizeof(int) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_N_CELLS,
				   &n_cells,
				   sizeof(int) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_N_ATOM_ARRAY,
				   &n_atom_array,
				   sizeof(int) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_N_CELLS_X,
				   &n_cells_x,
				   sizeof(int) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_N_CELLS_Y,
				   &n_cells_y,
				   sizeof(int) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_N_COLUMNS,
				   &n_columns,
				   sizeof(int) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_N_UNI,
				   &n_uni,
				   sizeof(int) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_N_UNI_Z,
				   &n_uni_z,
				   sizeof(int) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_L_CELL_X,
				   &l_cell_x,
				   sizeof(real_pw) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_L_CELL_Y,
				   &l_cell_y,
				   sizeof(real_pw) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_L_UNI_Z,
				   &L_uni_z,
				   sizeof(real_pw) ));
  HANDLE_ERROR( cudaMemcpyToSymbol(D_N_NEIGHBOR_COL_X,
				   &n_neighbor_col_x,
				   sizeof(int) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_N_NEIGHBOR_COL_Y,
				   &n_neighbor_col_y,
				   sizeof(int) ) );
  const int n_neighbor_col = (n_neighbor_col_x*2+1) * (n_neighbor_col_y*2+1);
  HANDLE_ERROR( cudaMemcpyToSymbol(D_N_NEIGHBOR_COL,
				   &n_neighbor_col,
				   sizeof(int) ) );
  const int max_n_cell_pairs_per_column = max_n_cell_pairs / (n_cells * n_neighbor_col);
  HANDLE_ERROR( cudaMemcpyToSymbol(D_MAX_N_CELL_PAIRS_PER_COLUMN,
				   &max_n_cell_pairs_per_column,
				   sizeof(int) ) );  
  const int max_n_cell_pairs_per_cell = max_n_cell_pairs_per_column * n_neighbor_col;
  HANDLE_ERROR( cudaMemcpyToSymbol(D_MAX_N_CELL_PAIRS_PER_CELL,
				   &max_n_cell_pairs_per_cell,
				   sizeof(int) ) );  
  return 0;
}

// cuda_set_constant
//   called only onece at the beginning of simulation
extern "C" int cuda_set_constant(int n_atoms, real_pw cutoff,
				 real_pw cutoff_pairlist, int n_atomtypes){
  real_pw tmp_charge_coeff = (real_pw)CelesteObject::CHARGE_COEFF;
  HANDLE_ERROR( cudaMemcpyToSymbol(D_CHARGE_COEFF,
				   &tmp_charge_coeff,
				   sizeof(real_pw) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_N_ATOMS,
				   &n_atoms,
				   sizeof(int) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_N_ATOMTYPES,
				   &n_atomtypes,
				   sizeof(int) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_CUTOFF,
				   &cutoff,
				   sizeof(real_pw) ));
  HANDLE_ERROR( cudaMemcpyToSymbol(D_CUTOFF_PAIRLIST,
				   &cutoff_pairlist,
				   sizeof(real_pw) ));
  return 0;
}


extern "C" int cuda_alloc_set_lj_params(real_pw* h_lj_6term,
					real_pw* h_lj_12term,
					int n_lj_types,
					int* h_nb15off,
					const int max_n_nb15off,
					const int max_n_atoms,
					const int max_n_atom_array){
  //printf("threads : %d\n", PW_THREADS);
  printf("cuda_alloc_set_lj_params\n");
  const unsigned int size_lj_matrix = sizeof(real_pw) * n_lj_types * n_lj_types;
  // cudaMalloc
  HANDLE_ERROR( cudaMalloc( (void**)&d_lj_6term,
			    size_lj_matrix ) );
  HANDLE_ERROR( cudaMalloc( (void**)&d_lj_12term,
			    size_lj_matrix ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_MAX_N_NB15OFF,
				   &max_n_nb15off,
				   sizeof(int) ) );  
  const unsigned int size_nb15off_orig = sizeof(int)*max_n_nb15off*max_n_atoms;
  HANDLE_ERROR( cudaMalloc( (void**)&d_nb15off_orig,
			    size_nb15off_orig));
  const unsigned int size_nb15off = sizeof(int)*max_n_nb15off*max_n_atom_array;
  HANDLE_ERROR( cudaMalloc( (void**)&d_nb15off,
			    size_nb15off));
  // cudaBindTexture2D
  /*
    cudaChannelFormatDesc desc = cudaCreateChannelDesc<real>();
    HANDLE_ERROR( cudaBindTexture2D( NULL,
    tex_lj_6term,
    d_lj_6term,
    desc, n_lj_types, n_lj_types,
    sizeof(real)*n_lj_types) );
    HANDLE_ERROR( cudaBindTexture2D( NULL,
				   tex_lj_12term,
				   d_lj_12term,
				   desc, n_lj_types, n_lj_types,
				   sizeof(real)*n_lj_types) );
  */
  // cudaMempcpy
  HANDLE_ERROR( cudaMemcpy( d_lj_6term, h_lj_6term,
			    size_lj_matrix,
			    cudaMemcpyHostToDevice) );
  HANDLE_ERROR( cudaMemcpy( d_lj_12term, h_lj_12term,
			    size_lj_matrix,
			    cudaMemcpyHostToDevice) );
  HANDLE_ERROR( cudaMemcpy( d_nb15off_orig, h_nb15off,
			    size_nb15off_orig,
			    cudaMemcpyHostToDevice) );
  
  return 0;
}

extern "C" int cuda_free_lj_params(){
  //printf("cuda_free_lj_param\n");
  //cudaUnbindTexture(tex_lj_6term);
  //cudaUnbindTexture(tex_lj_12term);
  HANDLE_ERROR( cudaFree(d_lj_6term) );
  HANDLE_ERROR( cudaFree(d_lj_12term) );
  HANDLE_ERROR( cudaFree(d_nb15off_orig) );
  HANDLE_ERROR( cudaFree(d_nb15off) );
  return 0;
}

// cuda_hostalloc_atom_type_charge
extern "C" int cuda_hostalloc_atom_type_charge(int*& h_atom_type,
					       real_pw*& h_charge,
					       int n_atoms){
  printf("hostalloc atom_type_charge cu\n");
  HANDLE_ERROR( cudaHostAlloc( (void**)&h_atom_type,
			       n_atoms * sizeof(int),
			       cudaHostAllocDefault));
  HANDLE_ERROR( cudaHostAlloc( (void**)&h_charge,
			       n_atoms * sizeof(real_pw),
			       cudaHostAllocDefault));
  return 0;
}

// cuda_hostalloc_atom_info
//   Allocation for MiniCell members
extern "C" int cuda_hostalloc_atom_info(real_pw*& h_crd, int*& h_atomids,
					real_fc*& h_work, real_fc*& h_energy,
					int n_atom_array){
  printf("hostalloc_atom_info : \n");
  printf("hostalloc_atom_info %d\n", n_atom_array);
  HANDLE_ERROR( cudaHostAlloc( (void**)&h_crd,
			       n_atom_array * 3 * sizeof(real_pw),
			       cudaHostAllocDefault));
  HANDLE_ERROR( cudaHostAlloc( (void**)&h_atomids,
			       n_atom_array * sizeof(int),
			       cudaHostAllocDefault));
  HANDLE_ERROR( cudaHostAlloc( (void**)&h_work,
			       n_atom_array * 3 * sizeof(real_fc),
			       cudaHostAllocDefault));
  HANDLE_ERROR( cudaHostAlloc( (void**)&h_energy,
			       2 * sizeof(real_fc),
			       cudaHostAllocDefault));
 return 0;
}

extern "C" int cuda_hostalloc_cell_info(//CellPair*& h_cell_pairs, 
					//int*& h_idx_head_cell_pairs,
					int*& h_idx_xy_head_cell,
					//int max_n_cell_pairs,
					//int max_n_cells,
					int n_columns){
  printf("cuda_hostalloc_cell_info cu\n");
  //HANDLE_ERROR( cudaHostAlloc( (void**)&h_cell_pairs,
  //max_n_cell_pairs * sizeof(CellPair),
  //cudaHostAllocDefault));
  //HANDLE_ERROR( cudaHostAlloc( (void**)&h_idx_head_cell_pairs,
  //(max_n_cells) * sizeof(int),
  //cudaHostAllocDefault));
  HANDLE_ERROR( cudaHostAlloc( (void**)&h_idx_xy_head_cell,
			       (n_columns) * sizeof(int),
			       cudaHostAllocDefault));

  return 0;
}



extern "C" int cuda_hostfree_atom_type_charge(int* h_atom_type, real_pw* h_charge){
  HANDLE_ERROR( cudaFreeHost(h_atom_type));
  HANDLE_ERROR( cudaFreeHost(h_charge));
  return 0;
}

extern "C" int cuda_hostfree_atom_info(real_pw* h_crd, int* h_atomids,
				       real_fc*& h_work, real_fc*& h_energy){
  HANDLE_ERROR( cudaFreeHost(h_crd));
  HANDLE_ERROR( cudaFreeHost(h_atomids));
  HANDLE_ERROR( cudaFreeHost(h_work));
  HANDLE_ERROR( cudaFreeHost(h_energy));
  return 0;
}
extern "C" int cuda_hostfree_cell_info(//CellPair* h_cell_pairs,
				       //int* h_idx_head_cell_pairs,
				       int* h_idx_xy_head_cell){
  //HANDLE_ERROR( cudaFreeHost(h_cell_pairs) );
  //HANDLE_ERROR( cudaFreeHost(h_idx_head_cell_pairs) );
  //HANDLE_ERROR( cudaFreeHost(h_idx_xy_head_cell) );
  return 0;
}
//__global__ void kernel_set_cellinfo(int* d_cell_pair_removed,
//				    const int n_cells){
//  const int cell_id = blockDim.x * blockIdx.x + threadIdx.x;
//  if(cell_id >= n_cells) return;
//  d_cell_pair_removed[cell_id] = 0;
//}
__global__ void kernel_set_nb15off(const int* d_atomids,
				   const int* d_atomids_rev,
				   const int* d_nb15off_orig,
				   int* d_nb15off){
  const int g_thread_idx = threadIdx.x + blockDim.x * blockIdx.x;
  const int atomid = g_thread_idx/D_MAX_N_NB15OFF;
  const int idx = g_thread_idx%D_MAX_N_NB15OFF;
  if(atomid >= D_N_ATOM_ARRAY) return;
  if(d_atomids[atomid] < 0){
    d_nb15off[g_thread_idx] = atomid;
  }else{
    const int orig = d_nb15off_orig[d_atomids[atomid]*D_MAX_N_NB15OFF + idx];
    if(orig == -1){
      d_nb15off[g_thread_idx] = -1;
    }else{
      d_nb15off[g_thread_idx] = d_atomids_rev[orig];
    }
  }
  //if(atomid == 0){
  //printf("nb15off[%d] atom: %d (%d, rev:%d) nb[%d]: %d (%d) \n",
  //g_thread_idx,
  //atomid, d_atomids[atomid], d_atomids_rev[d_atomids[atomid]],
  //idx,
  //d_nb15off[g_thread_idx], d_nb15off_orig[d_atomids[atomid]]);
  //}
}
__global__ void kernel_set_atominfo(const int* d_atomids,
				    const int* d_atomtype_orig,
				    const real_pw* d_charge_orig,
				    int* d_atomtype,
				    real4* d_crd_chg,
				    int* d_atomids_rev){
  int atomid = threadIdx.x + blockIdx.x * blockDim.x;
  if(atomid < D_N_ATOM_ARRAY && d_atomids[atomid] >= 0){
    d_atomtype[atomid] = d_atomtype_orig[d_atomids[atomid]];
    d_crd_chg[atomid].w = d_charge_orig[d_atomids[atomid]];
    d_atomids_rev[d_atomids[atomid]] = atomid;
  }
}
__global__ void kernel_set_crd(const int* d_atomids,
			       const real_pw* d_crd,
			       real4* d_crd_chg){
  int atomid = threadIdx.x + blockIdx.x * blockDim.x;
  if(atomid < D_N_ATOM_ARRAY){
    d_crd_chg[atomid].x = d_crd[atomid*3+0];
    d_crd_chg[atomid].y = d_crd[atomid*3+1];
    d_crd_chg[atomid].z = d_crd[atomid*3+2];
    /*if(d_crd_chg[atomid].x <= PBC_LOWER_BOUND[0] ||
       d_crd_chg[atomid].x > PBC_L[0]+  PBC_LOWER_BOUND[0] ||
       d_crd_chg[atomid].y <= PBC_LOWER_BOUND[1] ||
       d_crd_chg[atomid].y > PBC_L[1]+  PBC_LOWER_BOUND[2] ||
       d_crd_chg[atomid].z <= PBC_LOWER_BOUND[2] ||
       d_crd_chg[atomid].z > PBC_L[2]+  PBC_LOWER_BOUND[2]){
      printf("crd !? %d %f %f %f\n",atomid,d_crd_chg[atomid].x,
	     d_crd_chg[atomid].y, d_crd_chg[atomid].z);
	     }*/
  }
}
/*
__Global__ void kernel_set_atomids_rev(const int* d_atomids, int* d_atomids_rev){
  int atomid = threadIdx.x + blockIdx.x * blockDim.x;
  if(atomid < D_N_ATOMS){
    d_atomids_rev[d_atomids[atomid]] = atomid;
  }
}
*/
extern "C" int cuda_set_atominfo(const int n_atom_array,
				 const int max_n_nb15off, const int max_n_cells){
  HANDLE_ERROR( cudaMemset(d_atomids_rev, -1, sizeof(int)*n_atom_array ));
  HANDLE_ERROR( cudaMemset(d_nb15off, -1, sizeof(int)*n_atom_array*max_n_nb15off ));
  HANDLE_ERROR( cudaMemset(d_n_cell_pairs, -1, sizeof(int)*max_n_cells ));
  HANDLE_ERROR( cudaMemset(d_cell_pair_removed, 0, sizeof(int)*max_n_cells ));

  int blocks1 = (n_atom_array + REORDER_THREADS-1) / REORDER_THREADS;
  kernel_set_atominfo<<<blocks1, REORDER_THREADS>>>(d_atomids,
						   d_atomtype_orig,
						   d_charge_orig,
						   d_atomtype,
						    d_crd_chg,
						    d_atomids_rev);
  int blocks2 = (n_atom_array*max_n_nb15off + REORDER_THREADS-1) / REORDER_THREADS;
  kernel_set_nb15off<<<blocks2, REORDER_THREADS>>>(d_atomids,
						   d_atomids_rev,
						   d_nb15off_orig,
						   d_nb15off);
  return 0;
}
extern "C" int cuda_set_crd(int n_atom_array){
  int blocks = (n_atom_array + REORDER_THREADS-1) / REORDER_THREADS;
  kernel_set_crd<<<blocks, REORDER_THREADS>>>(d_atomids,
					      d_crd,
					      d_crd_chg);
  return 0;
}

__global__ void kernel_set_work_orig(real_fc* d_work, //real_fc* d_work_orig,
				     const int* d_atomids){
  int atomid = threadIdx.x + blockIdx.x * blockDim.x;
  if(atomid < D_N_ATOM_ARRAY){// && d_atomids[atomid] >= 0){
    //printf("atomid %d %d %f %f %f\n", atomid, d_atomids[atomid], d_work[atomid*3], d_work[atomid*3+1], d_work[atomid*3+2]);
    int index_orig = atomid*3;
    //d_work_orig[index_orig+0] = 0.0;
    //d_work_orig[index_orig+1] = 0.0;
    //d_work_orig[index_orig+2] = 0.0;
    for(int n=1; n < N_MULTI_WORK; n++){
      int index = (atomid + D_N_ATOM_ARRAY*n) * 3;
      d_work[index_orig+0] += d_work[index+0];
      d_work[index_orig+1] += d_work[index+1];
      d_work[index_orig+2] += d_work[index+2];
    }
  }
}

__global__ void kernel_reduction_energy(real_fc* d_energy){
  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  for(int i=1; i< N_MULTI_WORK; i++){
    if(tid==1){
      d_energy[0] += d_energy[i*2];
      d_energy[1] += d_energy[i*2+1];
      //printf("ene %f %f\n",d_energy[0],d_energy[1]);
    }
  }
  /*
  for(unsigned int s = N_MULTI_WORK;
      s > 1;  s >>= 1){
    if(tid < s){
      d_energy[tid] += d_energy[tid + s];
      //printf("reduction, %d: %f\n", tid, );
      
      //printf("d_energy[%d] (%f) += d_energy[%d] (%f) \n", tid, d_energy[tid], tid+s, d_energy[tid+s]);
    }
    //__syncthreads();
  }
  */
}

__device__ bool check_15off64(const int atom_idx1, const int atom_idx2,
			      const int* bitmask, int& mask_id, int& interact_bit){
  int bit_pos = atom_idx2 * N_ATOM_CELL + atom_idx1;
  mask_id =  bit_pos / 32;
  //int interact_pos =  bit_pos % 32;
  interact_bit = 1 << ( bit_pos % 32 ) ;
  return (bitmask[mask_id] & interact_bit) == interact_bit;
}
//__device__ int disable_15pair(const int interact_bit,
//			      int& bitmask_i){
//  bitmask_i = bitmask_i & ~interact_bit;
//  return 0;
//}

__device__ real_pw check_15off(const int atomid1, const int atomid2,
			       const int tmask_a1,
			       const int tmask_a2){
  int aid_diff = atomid2 - atomid1;
  int target = tmask_a1;
  if(aid_diff < 0){aid_diff = -aid_diff; target=tmask_a2; }
  
  int mask = 0;
  if(aid_diff <= 32)  mask = 1 << (aid_diff-1);
  real_pw valid_pair = 1.0;
  if(mask != 0 && (mask & target) == mask)
    valid_pair = 0.0;
  return valid_pair;
}
__device__ real_pw cal_pair(real_pw& w1,
			real_pw& w2,
			real_pw& w3,
			real_pw& ene_vdw,
			real_pw& ene_ele,
			const real4& crd_chg1, 
			const real4& crd_chg2, 
			const int& atomtype1,
			const int& atomtype2,
			const real_pw* __restrict__ d_lj_6term,
			const real_pw* __restrict__ d_lj_12term
			){

  const real_pw d12[3] = {
    crd_chg1.x - crd_chg2.x,
    crd_chg1.y - crd_chg2.y,
    crd_chg1.z - crd_chg2.z
  };
  const real_pw r12_2 = d12[0]*d12[0] + d12[1]*d12[1] + d12[2]*d12[2];
  const real_pw r12 = sqrt(r12_2);
  const real_pw r12_inv = 1.0 / r12;
  const real_pw r12_2_inv = r12_inv * r12_inv;
  const real_pw r12_3_inv = r12_inv * r12_2_inv;
  const real_pw r12_6_inv = r12_3_inv * r12_3_inv;
  const real_pw r12_12_inv = r12_6_inv * r12_6_inv;
  const real_pw term6 = d_lj_6term[atomtype1*D_N_ATOMTYPES + atomtype2] * r12_6_inv;
  const real_pw term12 = d_lj_12term[atomtype1*D_N_ATOMTYPES + atomtype2] * r12_12_inv;
  real_pw work_coef = r12_2_inv * (-12.0 * term12 + 6.0 * term6);
  const real_pw cc = crd_chg1.w * crd_chg2.w * D_CHARGE_COEFF;
  work_coef -= cc * (r12_3_inv - D_FCOEFF);
  //  printf("dbgpair %10e %d %d %d %d %10e %10e %10e\n", r12, atomid1, atomid2,  a1, a2, D_CHARGE_COEFF, D_ZCORE, D_BCOEFF);

  if(r12 >= D_CUTOFF) {return r12;}
  //pairs[0]++;
  
  //work_coef *= valid_pair;
  w1 = (work_coef) * d12[0];
  w2 = (work_coef) * d12[1];
  w3 = (work_coef) * d12[2];
  //real_pw tmp1 = D_BCOEFF*r12_2;
  //real_pw tmp2 = tmp1 - D_ZCORE;
  //real_pw tmp3 = tmp2 + r12_inv;
  //ene_ele = cc * tmp3;
  ene_ele = cc * (r12_inv - D_ZCORE + D_BCOEFF * r12_2);// * valid_pair;

  ene_vdw = (-term6 + term12);// * valid_pair;
  //if(valid_pair > 0.5 && (ene_vdw!=0.0 || ene_ele!=0.0)) pairs[1]++;
  //printf("dbgpair %d-%d (%15e, %15e, %15e) (%15e, %15e %15e) %15e %15e %15e\n", a1, a2,
  //crd_chg1.x,crd_chg1.y,crd_chg1.z,
  //	 crd_chg2.x,crd_chg2.y,crd_chg2.z,
  //r12, cc, D_CHARGE_COEFF);
  /*
  printf("dbgcrd %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\ndbgpair %d-%d %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
	 crd_chg1.x,crd_chg1.y,crd_chg1.z,
	 crd_chg2.x,crd_chg2.y,crd_chg2.z,
	 a1,a2,
	 ene_ele,
	 D_BCOEFF, r12_2, tmp1, tmp2, tmp3,
	 (r12_inv + (D_BCOEFF * r12_2))-D_ZCORE,
	 r12_inv + ((D_BCOEFF * r12_2)-D_ZCORE),
	 r12_inv+(-D_ZCORE+( D_BCOEFF*r12_2)),
	 (r12_inv-D_ZCORE)+( D_BCOEFF*r12_2),
	 (( D_BCOEFF*r12_2)+r12_inv)-D_ZCORE,
	 ( D_BCOEFF*r12_2)+(r12_inv-D_ZCORE)

	 );
  */
  return r12;
}

__global__ void kernel_pairwise_ljzd(const real4* d_crd_chg,
				     CellPair* d_cell_pairs,
				     const int* d_idx_head_cell_pairs,
				     const int* d_atomtype,
				     const real_pw* __restrict__ d_lj_6term,
				     const real_pw* __restrict__ d_lj_12term,
				     real_fc* d_energy, real_fc* d_work,
				     //int* d_pairs, real_pw* d_realbuf,
				     const int offset_cells,
				     const bool flg_mod_15mask,
				     const int* d_n_cell_pairs){
  real_fc ene_vdw = 0.0;
  real_fc ene_ele = 0.0;

  const int global_threadIdx = blockDim.x * blockIdx.x + threadIdx.x;
  
  const int c1 = global_threadIdx >> 5;
  const int warpIdx = threadIdx.x >> 5;  
  if(c1 >= D_N_CELLS){ return; }
  const int laneIdx = global_threadIdx & 31;
  //;const int n_loops = (d_idx_head_cell_pairs[c1+1] - d_idx_head_cell_pairs[c1] - d_cell_pair_removed[c1])*2;
  const int n_loops = d_n_cell_pairs[c1]*2;//(D_MAX_N_CELL_PAIRS_PER_CELL - d_cell_pair_removed[c1])*2;

  const int ene_index_offset = global_threadIdx%N_MULTI_WORK;  
  
  real_fc work_c1[3] = {0.0, 0.0, 0.0};

  const int atom_idx1 = (laneIdx & 7);  // laneIdx%8
  const int a1 = c1 * N_ATOM_CELL + atom_idx1;
  
  __shared__ real4 crd_chg1[N_ATOM_CELL * (PW_THREADS >> 5)];
  __shared__ int   atomtype1[N_ATOM_CELL * (PW_THREADS >> 5)];

  //__shared__ real4 crd_chg2[N_ATOM_CELL * PW_THREADS / 32];
  const int sharedmem_idx  = N_ATOM_CELL * warpIdx + atom_idx1;
  if(laneIdx<N_ATOM_CELL){
    crd_chg1[sharedmem_idx] = d_crd_chg[c1*N_ATOM_CELL + laneIdx];
    atomtype1[sharedmem_idx] = d_atomtype[c1*N_ATOM_CELL + laneIdx];
  }
  __syncthreads();

  CellPair cellpair;
  int cp;
  for(int loopIdx=0; loopIdx < n_loops; loopIdx++){
    if(loopIdx%2 == 0){
      cp = d_idx_head_cell_pairs[c1] + (loopIdx >> 1);
      if(cp >= D_MAX_N_CELL_PAIRS) break;
      cellpair = d_cell_pairs[cp];
    }
    //atomicAdd(&d_pairs[3],1);
    const int c2 = cellpair.cell_id2;

    // atom_idx ... index in cell, 0-7
    const int atom_idx2 = (laneIdx / 8)  + 4 * (loopIdx % 2);  // laneIdx/8 + 4*(warpIdx%2)

    // remove 1-2, 1-3, 1-4 pairs
    const int a2 = c2 * N_ATOM_CELL + atom_idx2;
    //real4 crd_chg2;
    //int2 atominfo2;
    //if(atom_idx1 == 0){
    real4 crd_chg2 = d_crd_chg[a2];
    const int atomtype2 = d_atomtype[a2];
    //}
    //int atomid2_top = laneIdx - laneIdx%8;
    //crd_chg2.x = __shfl(crd_chg2.x, laneIdx - atom_idx1);
    //crd_chg2.y = __shfl(crd_chg2.y, laneIdx - atom_idx1);
    //crd_chg2.z = __shfl(crd_chg2.z, laneIdx - atom_idx1);
    //crd_chg2.w = __shfl(crd_chg2.w, laneIdx - atom_idx1);
    //atominfo2.x = __shfl(atominfo2.x, laneIdx - atom_idx1);
    //atominfo2.y = __shfl(atominfo2.y, laneIdx - atom_idx1);

    if      ( (cellpair.image & 1) == 1 )   crd_chg2.x -= PBC_L[0];
    else if ( (cellpair.image & 2) == 2 )   crd_chg2.x += PBC_L[0];
    if      ( (cellpair.image & 4) == 4 )   crd_chg2.y -= PBC_L[1];
    else if ( (cellpair.image & 8) == 8 )   crd_chg2.y += PBC_L[1];
    if      ( (cellpair.image & 16) == 16 ) crd_chg2.z -= PBC_L[2];
    else if ( (cellpair.image & 32) == 32 ) crd_chg2.z += PBC_L[2];

    real_pw w1=0.0, w2=0.0, w3=0.0;
    real_pw cur_ene_ele=0.0;
    real_pw cur_ene_vdw=0.0;
    int mask_id;
    int interact_bit;
    //if (threadIdx.x + blockIdx.x * blockDim.x == 0){
    //printf("cp: %d (%d-%d) at: %d-%d mask: %d %d \n", cp, c1, c2, atom_idx1, atom_idx2,
    //d_cell_pairs[cp].pair_mask[0], d_cell_pairs[cp].pair_mask[1]);
    //}      

    if(!check_15off64(atom_idx1, atom_idx2, d_cell_pairs[cp].pair_mask,
		      mask_id, interact_bit)){
      //if(d_atominfo[a1].x == -1 || d_atominfo[a2].x == -1) continue;
      real_pw r12 = cal_pair(w1, w2, w3, cur_ene_vdw, cur_ene_ele,
			     crd_chg1[sharedmem_idx],
			     crd_chg2,
			     atomtype1[sharedmem_idx],
			     atomtype2,
			     d_lj_6term, d_lj_12term);
      //if (threadIdx.x + blockIdx.x * blockDim.x == 0){
      //printf("r %f ene %f %f\n", r12, cur_ene_vdw, cur_ene_ele);
      //}      


      if(flg_mod_15mask && r12 < D_CUTOFF_PAIRLIST) interact_bit = 0;
      //atomicAdd(&d_pairs[1],1);
      
      //debug
      /*
      const real_pw d12[3] = {
	d_crd_chg[a1].x - crd_chg2.x,
	d_crd_chg[a1].y - crd_chg2.y,
	d_crd_chg[a1].z - crd_chg2.z
      };
      const real_pw r12_2 = d12[0]*d12[0] + d12[1]*d12[1] + d12[2]*d12[2];
      const real_pw r12 = sqrt(r12_2);
      if(r12 < D_CUTOFF){
	atomicAdd(&d_pairs[4],1);
      }
      //atomicAdd(&d_realbuf[0], r12);
      //if(global_threadIdx == 1){
      //while(d_realbuf[0] > 100000){
      //atomicAdd(&d_realbuf[0], -100000);
      //}
      //}
      atomicAdd(&d_realbuf[1], d_crd_chg[a1].w*d_crd_chg[a2].w);
      if(global_threadIdx == 1){
	while(d_realbuf[1] > 100000){
	  atomicAdd(&d_realbuf[1], -100000);
	}
      }
      */
      //
  
      //if(cur_ene_vdw!=0.0 || cur_ene_ele!=0.0) atomicAdd(&d_pairs[2],1);
      //printf("dbg02ene %d a(%d-%d) c(%d-%d): %10e %10e %10e %10e\n",global_threadIdx, a1, a2,c1,c2,ene_vdw, ene_ele,cur_ene_vdw, cur_ene_ele);
      //return;

      ene_vdw += cur_ene_vdw;
      ene_ele += cur_ene_ele;
      //printf("dbg2 %d %d \n",atom_idx1, atom_idx2);
      
      work_c1[0] += w1;
      work_c1[1] += w2;
      work_c1[2] += w3;
      //work_c1[atom_idx1 * 3 + 0] += w1;
      //work_c1[atom_idx1 * 3 + 1] += w2;
    //work_c1[atom_idx1 * 3 + 2] += w3;
    
    //work_c2[0] -= w1;
    //work_c2[1] -= w2;
    //work_c2[2] -= w3;
    //printf("dbg3\n");
    //    if(w1 != 0.0) printf("cp %d: %d-%d\n",cp, c1, c2);
    //printf("SHFL wapridx:%d lane:%d a2:%d w:%12.8e %12.8e %12.8e\n", warpIdx, laneIdx, a2, w1, w2, w3);
    }    
    if(flg_mod_15mask){
      for(int i = 32; i >= 1; i/=2){
	interact_bit |= __shfl_xor(interact_bit, i);
      }
      if(laneIdx == 0)
	d_cell_pairs[cp].pair_mask[mask_id]  |= interact_bit;
    }
    for(int i = 4; i >= 1; i/=2){
      w1 += shfl_xor(w1, i, 8);
      w2 += shfl_xor(w2, i, 8);
      w3 += shfl_xor(w3, i, 8);
      //printf("SHFL[%d] wapridx:%d lane:%d a2:%d w:%12.8e %12.8e %12.8e\n",i, warpIdx, laneIdx,  a2, w1, w2, w3);
    }
    if(laneIdx % 8 == 0 && w1 != 0.0 && w2 != 0.0 && w3 != 0.0){
      const int tmp_index = (((global_threadIdx/32)%N_MULTI_WORK)*D_N_ATOM_ARRAY + a2) * 3;
      //const int tmp_index = ((ene_index_offset*D_N_ATOM_ARRAY + a2) * 3);
      atomicAdd(&(d_work[tmp_index+0]), -w1);
      atomicAdd(&(d_work[tmp_index+1]), -w2);
      atomicAdd(&(d_work[tmp_index+2]), -w3);
    }
      //printf("dbg4\n");
  }
  //printf("dbg5\n");
  //int a1 = c1 * N_ATOM_CELL;
  //for(int atom_idx1 = 0; atom_idx1 < N_ATOM_CELL; atom_idx1++, a1++){
  //for(int i = 16; i <= 8; i/=2){
  //work_c1[0] += shfl_xor(work_c1, i, 32);    
  //}  
  //const int tmp_index =  (((global_threadIdx%N_MULTI_WORK)*D_N_ATOM_ARRAY) + a1)*3;
  for(int i = 16; i >= 8; i/=2){
    work_c1[0] += shfl_xor(work_c1[0], i, 32);
    work_c1[1] += shfl_xor(work_c1[1], i, 32);
    work_c1[2] += shfl_xor(work_c1[2], i, 32);
  }
  if(laneIdx < 8){
    const int tmp_index =  ((ene_index_offset*D_N_ATOM_ARRAY) + a1)*3;
    atomicAdd(&(d_work[tmp_index+0]), work_c1[0]);
    atomicAdd(&(d_work[tmp_index+1]), work_c1[1]);
    atomicAdd(&(d_work[tmp_index+2]), work_c1[2]);
  }
  //}
  //printf("dbg8\n");
  //printf("dbg01b_ene %d %d: %10e %10e %10e %10e\n",ene_index_offset, global_threadIdx, d_energy[ene_index_offset*2], d_energy[ene_index_offset*2+1], ene_vdw, ene_ele);
  for(int i = 16; i >= 1; i/=2){
    ene_vdw += shfl_xor(ene_vdw, i, 32);
    ene_ele += shfl_xor(ene_ele, i, 32);
  }
  if(laneIdx == 0){
    const int tmp_index = ((global_threadIdx/32)%N_MULTI_WORK)*2;
    atomicAdd(&(d_energy[tmp_index+0]), ene_vdw);
    atomicAdd(&(d_energy[tmp_index+1]), ene_ele);
  }
  //printf("dbg10\n");  
  //printf("dbg01ene %d: %f %f\n",global_threadIdx, d_energy[ene_index_offset*2], d_energy[ene_index_offset*2+1]);
}
//extern "C" int cuda_init_cellinfo(const int n_cells){
//  int blocks = (n_cells+PW_THREADS/32-1) / (PW_THREADS/32);
//  kernel_set_cellinfo<<<blocks, PW_THREADS>>>(d_cell_pair_removed, n_cells);
//  return 0;
//}

__global__ void set_idx_head_cell_pairs(const int* d_n_cell_pairs,
					int* d_idx_head_cell_pairs){
  
  const int g_thread_id = blockDim.x * blockIdx.x + threadIdx.x;    
  if(g_thread_id==0){
    int idx_cp = 0;
    for(int cell_id = 0; cell_id < D_N_CELLS; cell_id++){
      d_idx_head_cell_pairs[cell_id] = idx_cp;
      idx_cp += d_n_cell_pairs[cell_id];
      idx_cp = (idx_cp + CP_PER_THREAD -1)/CP_PER_THREAD;
    }
  }
}
__global__ void pack_cellpairs_array(CellPair* d_cell_pairs,
				     CellPair* d_cell_pairs_buf,
				     int* d_n_cell_pairs,
				     int* d_idx_head_cell_pairs){
  
  const int cp = blockDim.x * blockIdx.x + threadIdx.x;
  const CellPair cellpair = d_cell_pairs_buf[cp];
  const int cp_in_cell1 = cp - cellpair.cell_id1*D_MAX_N_CELL_PAIRS_PER_CELL;
  if(cp_in_cell1 >= d_n_cell_pairs[cellpair.cell_id1]) return;
  const int dest = d_idx_head_cell_pairs[cellpair.cell_id1] + cp_in_cell1;
  d_cell_pairs[dest] = cellpair;
}

__global__ void kernel_reset_cellpairs(CellPair* d_cell_pairs,
				       //int* d_idx_head_cell_pairs,
				       int* d_n_cell_pairs,
				       const int n_cells){
  const int cell1_id = blockDim.x * blockIdx.x + threadIdx.x;
  if(cell1_id >= n_cells){ return; }
  int n_cp1 = d_n_cell_pairs[cell1_id];
  for(int cell2=0; cell2<n_cp1; cell2++){
    bool flg = true;
    int n_mask_int = (N_ATOM_CELL * N_ATOM_CELL + 31) / 32;
    const int cp = D_MAX_N_CELL_PAIRS_PER_CELL*cell1_id + cell2;
    for(int i=0; i < n_mask_int; i++)
      flg &= (d_cell_pairs[cp].pair_mask[i] == ~0);
    if(flg){
      d_n_cell_pairs[cell1_id]--;
      int cp_src = D_MAX_N_CELL_PAIRS_PER_CELL*cell1_id + --n_cp1;
      d_cell_pairs[cp] = d_cell_pairs[cp_src];
    }
  }
  d_n_cell_pairs[cell1_id] = n_cp1;  
}

extern "C" int cuda_pairwise_ljzd(const int offset_cellpairs, const int n_cal_cellpairs,
				  const int offset_cells,     const int n_cal_cells,
				  const bool flg_mod_15mask){
  // test
  //real_pw w1=0.0, w2=0.0, w3=0.0;
  //real_pw cur_ene_ele=0.0;
  //real_pw cur_ene_vdw=0.0;
  //printf("test pair\n");
  
  //cal_pair(w1,w2,w3, cur_ene_vdw, cur_ene_ele, a1, a2,
  //d_crd_chg, d_atominfo, crd2,
  //d_lj_6term, d_lj_12term);

  cudaStreamCreate(&stream_pair_home);

  //int *d_pairs;
  //HANDLE_ERROR(cudaMalloc((void**)&d_pairs, sizeof(int)*5));
  //HANDLE_ERROR(cudaMemset(d_pairs, 0, sizeof(int)*5));
  //int h_pairs[5] = {0,0,0,0,0};
  //real_pw *d_realbuf;
  //HANDLE_ERROR(cudaMalloc((void**)&d_realbuf, sizeof(real_pw)*4));
  //HANDLE_ERROR(cudaMemset(d_realbuf, 0.0, sizeof(real_pw)*4));
  //real_pw h_realbuf[4] = {0,0,0,0};
  //HANDLE_ERROR(cudaMemcpy(d_pairs, h_pairs, sizeof(int)*5, cudaMemcpyHostToDevice));
  //HANDLE_ERROR(cudaMemcpy(d_realbuf, h_realbuf, sizeof(real)*4, cudaMemcpyHostToDevice));
  //

  const int blocks = (n_cal_cells+PW_THREADS/32-1) / (PW_THREADS/32);
  //printf("kernel_pairwize_ljzd %d\n", n_cal_cells);
  kernel_pairwise_ljzd<<<blocks, PW_THREADS,
    0, stream_pair_home>>>(d_crd_chg,
			   d_cell_pairs,
			   d_idx_head_cell_pairs,
			   d_atomtype,
			   d_lj_6term, d_lj_12term,
			   d_energy, d_work,// d_pairs, d_realbuf,
			   offset_cells, flg_mod_15mask,
			   d_n_cell_pairs);

  
  if(flg_mod_15mask){
    const int blocks2 = (n_cal_cells+PW_THREADS-1) / PW_THREADS;
    //kernel_reset_cellpairs<<<blocks2, PW_THREADS, 0, stream_pair_home>>>
    //(d_cell_pairs, d_n_cell_pairs, n_cal_cells);
  }

  //HANDLE_ERROR(cudaMemcpy(h_pairs, d_pairs, sizeof(int)*5, cudaMemcpyDeviceToHost));
  //HANDLE_ERROR(cudaMemcpy(h_realbuf, d_realbuf, sizeof(real)*4, cudaMemcpyDeviceToHost));
  //printf("n 15 pairs: nb:%d / interact:%d / all:%d / cellpairs:%d / incut:%d\n", h_pairs[0],h_pairs[1],h_pairs[2], h_pairs[3], h_pairs[4]); 
  //printf("sum_dist: %f  charge: %f\n",h_realbuf[0], h_realbuf[1]);
  //HANDLE_ERROR(cudaFree(d_pairs));
  //HANDLE_ERROR(cudaFree(d_realbuf));

  return 0;
}
extern "C" int cuda_thread_sync(){
  cudaThreadSynchronize();
  return 0;
}
extern "C" int cuda_pair_sync(){
  cudaStreamSynchronize(stream_pair_home);
  cudaStreamDestroy(stream_pair_home);
  return 0;
}

extern "C" int cuda_memcpy_dtoh_work(real_fc*& h_work, real_fc*& h_energy,
				     int n_atoms, int n_atom_array){
  //printf("! cuda_memcpy_dtoh_work\n");
  int blocks = (n_atom_array + REORDER_THREADS-1) / REORDER_THREADS;
  //printf("kernel_set_work_orig\n");

  cudaStream_t stream_reduction1;
  cudaStream_t stream_reduction2;
  cudaStreamCreate(&stream_reduction1);
  cudaStreamCreate(&stream_reduction2);
  kernel_set_work_orig<<<blocks, REORDER_THREADS, 0, stream_reduction1>>>(d_work, d_atomids);
  kernel_reduction_energy<<<1, REORDER_THREADS, 0, stream_reduction2>>>(d_energy);
  cudaStreamSynchronize(stream_reduction1);
  cudaStreamSynchronize(stream_reduction2);
  cudaStreamDestroy(stream_reduction1);
  cudaStreamDestroy(stream_reduction2);

  //kernel_set_work_orig<<<blocks, REORDER_THREADS>>>(d_work, d_work_orig, d_atomids);
  //printf("threadsync\n");
  //cudaThreadSynchronize();  
  //printf("memcpy\n");

  HANDLE_ERROR(cudaMemcpy(h_work, d_work,
			  sizeof(real_fc) * n_atom_array * 3,
			  cudaMemcpyDeviceToHost));
  HANDLE_ERROR(cudaMemcpy(h_energy, d_energy,
			  sizeof(real_fc) * 2,
			  cudaMemcpyDeviceToHost));

  //printf("memcpy ene %f %f\n", h_energy[0], h_energy[1]);
  return 0;
}
extern "C" int cuda_reset_work_ene(int n_atoms){
  HANDLE_ERROR(cudaMemset(d_work, 0.0, sizeof(real_fc)*n_atoms*3*N_MULTI_WORK));
  HANDLE_ERROR(cudaMemset(d_energy, 0.0, sizeof(real_pw)*2*N_MULTI_WORK));
  return 0;
}
__device__ int get_column_id_from_crd(const int x, const int y){
  return y*D_N_CELLS_X + x;
}
__device__ int get_uni_id_from_crd(const int x, const int y, const int z){
  return x * D_N_UNI_Z + y * (D_N_UNI_Z * D_N_CELLS_X) + z;  
}
__global__ void kernel_set_uniform_grid(const real4* d_crd_chg,
					int2* d_uni2cell_z
					){
  
  const int cell_id = threadIdx.x + blockIdx.x * blockDim.x;  
  if ( cell_id >= D_N_CELLS ) return;
  const int laneIdx = threadIdx.x%WARPSIZE;
  const int warpIdx = threadIdx.x/WARPSIZE;

  const real4 crd_chg = d_crd_chg[cell_id * N_ATOM_CELL];
  const int col_x = floor((crd_chg.x-PBC_LOWER_BOUND[0]) / D_L_CELL_X);
  const int col_y = floor((crd_chg.y-PBC_LOWER_BOUND[1]) / D_L_CELL_Y);

  const int uni_z_min = (int)((d_crd_chg[cell_id*N_ATOM_CELL].z-PBC_LOWER_BOUND[2]) / D_L_UNI_Z);
  //const int uni_id_min = get_uni_id_from_crd(col_x, col_y, uni_z_min);
  const int uni_z_max = (int)((d_crd_chg[cell_id*N_ATOM_CELL+N_ATOM_CELL-1].z - PBC_LOWER_BOUND[2]) / D_L_UNI_Z);
  
  for(int z = uni_z_min; z <= uni_z_max; z++){
    const int uni_id = get_uni_id_from_crd(col_x, col_y, z);
    //if (uni_id_min >= D_N_UNI || uni_id_min < 0){// || uni_id_max >= D_N_UNI || uni_id_max < 0){
    if(uni_id >= D_N_UNI || uni_id < 0){
      printf("DBG!! cell:%d/%d uni_z:%d-%d %d/%d x:%d/%d y:%d/%d z:%d/%d x:%f/%f y:%f/%f z: %f-%f %f\n",
	     cell_id, D_N_CELLS,
	     uni_z_min, uni_z_max,
	     uni_id, D_N_UNI,
	     col_x, D_N_CELLS_X,
	     col_y, D_N_CELLS_Y,
	     z, D_N_UNI_Z,
	     crd_chg.x-PBC_LOWER_BOUND[0],PBC_L[0],
	     crd_chg.y-PBC_LOWER_BOUND[1],PBC_L[1],
	     d_crd_chg[cell_id*N_ATOM_CELL].z-PBC_LOWER_BOUND[2],
	     d_crd_chg[cell_id*N_ATOM_CELL+N_ATOM_CELL-1].z - PBC_LOWER_BOUND[2],
	     D_L_UNI_Z);
      
    }
    atomicMin(&d_uni2cell_z[uni_id].x, cell_id);
    atomicMax(&d_uni2cell_z[uni_id].y, cell_id);
  }
  
}

__device__ bool check_valid_pair(const int cell1_id, const int cell2_id){
  const bool cell1_odd = cell1_id%2!=0;
  const bool cell2_odd = cell2_id%2!=0;
  if ( cell1_odd){
    if ((cell2_id < cell1_id && !cell2_odd ) || 
	(cell2_id > cell1_id && cell2_odd )) return false;
  }else{
    if ((cell2_id < cell1_id && cell2_odd ) || 
	(cell2_id > cell1_id && !cell2_odd )) return false;
  }
  return true;
}
__device__ int set_cell_pair_bitmask(const int cell_id1, const int cell_id2,
				     const int cell1_id_in_block,
				     const int* d_atomids,
				     const int* d_nb15off_orig,
				     const int* d_nb15off,
				     const int* s_nb15off,
				     const int  n_atom_cell2,
				     int* pair_mask){
  for(int i = 0; i < N_BITMASK; i++) pair_mask[i] = 0;
  int a1 = N_ATOM_CELL*cell_id1;
  for(int a1_cell = 0; a1_cell < N_ATOM_CELL; a1++, a1_cell++){
    bool flg1 = false;
    //    if(d_atomids[a1] == -1) flg1 = true;
    //if(d_nb15off[a1*D_MAX_N_NB15OFF] == a1) flg1 = true;
    if(s_nb15off[(cell1_id_in_block*N_ATOM_CELL+a1_cell)*D_MAX_N_NB15OFF] == a1) flg1 = true;
    int a2 = N_ATOM_CELL*cell_id2;
    for(int a2_cell = 0; 
	a2_cell < N_ATOM_CELL; a2++, a2_cell++){
	//a2_cell < n_atom_cell2; a2++, a2_cell++){
      const int bit_pos = a2_cell * N_ATOM_CELL + a1_cell;
      const int mask_id =  bit_pos / 32;
      const int mask_pos =  bit_pos % 32;
      const int add_bit = 1 << mask_pos;
      bool flg12 = false;
      if(flg1) flg12 = true;
      if (a2_cell >= n_atom_cell2) flg12=true;
      //else if (d_atomids[a2] == -1) flg12=true;
      //else if(d_nb15off[a2*D_MAX_N_NB15OFF] == a2) flg12=true;
      else if((cell_id1 == cell_id2 && a1 >= a2)) flg12=true;
      else{
	/*
	  const int tail = (d_atomids[a1]+1) * D_MAX_N_NB15OFF;
	  for(int i = d_atomids[a1] * D_MAX_N_NB15OFF;
	  i < tail && d_nb15off_orig[i] != -1; i++){
	  if(d_nb15off_orig[i] == d_atomids[a2]){ 
	*/   
	//const int tail = (a1+1) * D_MAX_N_NB15OFF;
	const int tail = (cell1_id_in_block*N_ATOM_CELL+a1_cell+1) * D_MAX_N_NB15OFF;
	//for(int i = a1 * D_MAX_N_NB15OFF;
	for(int i = tail-D_MAX_N_NB15OFF;
	    i < tail && s_nb15off[i] != -1; i++){
	  if(s_nb15off[i] == a2){
	    flg12=true;
	    break;}
	}
      }
      //if(atomids[a1]==1 || atomids[a2]==1) flg=true;
      if(flg12){
	//int tmp = pair_mask[mask_id];
	pair_mask[mask_id] |= add_bit;
      }
    }
    
  }
  return 0;
}
__device__ CellPair get_new_cell_pair(const int cell1_id, const int cell2_id,
				      const int cell1_id_in_block,
				      const int image[3],
				      const int* d_atomids,
				      const int* d_nb15off_orig,
				      const int* d_nb15off,
				      const int* s_nb15off,
				      const int n_atom_cell2){
  CellPair new_cp;
  new_cp.cell_id1 = cell1_id;
  new_cp.cell_id2 = cell2_id;
  int bit_image = 0;
  if(image[0] == -1)      bit_image = bit_image | 1;
  else if(image[0] == 1)  bit_image = bit_image | 2;
  if(image[1] == -1)      bit_image = bit_image | 4;
  else if(image[1] == 1)  bit_image = bit_image | 8;
  if(image[2] == -1)      bit_image = bit_image | 16;
  else if(image[2] == 1)  bit_image = bit_image | 32;
  new_cp.image = bit_image;
  //new_cp.pair_mask[0] = 0;
  //new_cp.pair_mask[1] = 0;
  set_cell_pair_bitmask(cell1_id, cell2_id,
			cell1_id_in_block,
			d_atomids,
			d_nb15off_orig, d_nb15off, s_nb15off,
			n_atom_cell2,
			new_cp.pair_mask);
  return new_cp;
}


__global__ void kernel_enumerate_cell_pair(const int2* d_uni2cell_z,
					   const real4* d_crd_chg,
					   const int* d_idx_xy_head_cell,
					   const int* d_atomids,
					   const int* d_nb15off_orig,
					   const int* d_nb15off,
					   int * d_n_cell_pairs,
					   CellPair* d_cell_pairs
					   ){
  // 1 warp calculates pairs with a cell
  const int g_thread_id = (threadIdx.x + blockIdx.x * blockDim.x);
  const int cell1_id = g_thread_id/D_N_NEIGHBOR_COL; 
  if(cell1_id >= D_N_CELLS) return;
  const int neighbor_col_id = g_thread_id%D_N_NEIGHBOR_COL; 
  const int cell1_id_in_block = cell1_id - blockIdx.x*blockDim.x/D_N_NEIGHBOR_COL;

  if(cell1_id_in_block >= MAX_N_CELL_BLOCK) {
    printf("The number of cells in each block exceeds the constant MAX_N_CELL_BLOCK: %d / %d",
	   cell1_id_in_block, MAX_N_CELL_BLOCK);
  }
  __shared__ int s_nb15off[MAX_N_CELL_BLOCK*N_ATOM_CELL*MAX_N_NB15OFF];
  //__shared__ int s_nb15off[(REORDER_THREADS/D_N_NEIGHBOR_COL +2) * D_MAX_N_NB15OFF];

  if(threadIdx.x == 0 || neighbor_col_id == 0){
    for(int i=0; i<N_ATOM_CELL*MAX_N_NB15OFF; i++){
      s_nb15off[cell1_id_in_block*N_ATOM_CELL*MAX_N_NB15OFF + i] = d_nb15off[cell1_id*N_ATOM_CELL*MAX_N_NB15OFF+i];
    }
  }
  __syncthreads();
  if(neighbor_col_id >= D_N_NEIGHBOR_COL) return;
  //const int laneIdx = threadIdx.x%WARPSIZE;
  //const int warpIdx = threadIdx.x/WARPSIZE;
  const real4 crd_chg11 = d_crd_chg[cell1_id*N_ATOM_CELL];
  const int cell1_x = floor((crd_chg11.x - PBC_LOWER_BOUND[0]) / D_L_CELL_X);
  const int cell1_y = floor((crd_chg11.y - PBC_LOWER_BOUND[1]) / D_L_CELL_Y);
  int first_uni_z[3] = {0, 0, 0};
  int last_uni_z[3] = {0, 0, 0};
  int tmp_first = floor((crd_chg11.z-PBC_LOWER_BOUND[2]) / D_L_UNI_Z) - 2;
  if(tmp_first < 0){
    first_uni_z[0] = D_N_UNI_Z + tmp_first;
    last_uni_z[0] = D_N_UNI_Z -1;
    tmp_first = 0;
  }else{
    first_uni_z[0] = -1;
    last_uni_z[0] = -1;
  }
  first_uni_z[1] = tmp_first;
  int tmp_last = floor((d_crd_chg[cell1_id*N_ATOM_CELL+N_ATOM_CELL-1].z - PBC_LOWER_BOUND[2])/D_L_UNI_Z) + 2;
  if(tmp_last >= D_N_UNI_Z){
    first_uni_z[2] = 0;
    last_uni_z[2] = tmp_last - D_N_UNI_Z;
    last_uni_z[1] = D_N_UNI_Z -1;
  }else{
    first_uni_z[2] = -1;
    last_uni_z[2] = -1;
    last_uni_z[1] = tmp_last;
  }
  int image[3] = {0,0,0};
  //if(cell1_id==0 && laneIdx==0)
  //printf("DBG n_col_thread:%d %d-%d %d-%d %d-%d\n", n_col_thread,
  //first_uni_z[0], last_uni_z[0],
  //first_uni_z[1], last_uni_z[1],
  //first_uni_z[2], last_uni_z[2]);

  int added_col = 0;
  //const int idx_cell_pair_head = D_MAX_N_CELL_PAIRS_PER_CELL * cell1_id +
  //D_MAX_N_CELL_PAIRS_PER_COLUMN * neighbor_col_id;
  const int idx_cell_pair_head = D_MAX_N_CELL_PAIRS_PER_CELL * cell1_id;
  
  //if(cell1_id == 0){
  //  printf("neighbor col: %d / %d \n",
  //neighbor_col_id, D_N_NEIGHBOR_COL);
    //}

  const int dx = neighbor_col_id % (D_N_NEIGHBOR_COL_X*2+1) - D_N_NEIGHBOR_COL_X;
  const int dy = neighbor_col_id / (D_N_NEIGHBOR_COL_X*2+1) - D_N_NEIGHBOR_COL_X;
  // for x
  image[0] = 0;
  const int rel_x = cell1_x + dx;
  int cell2_x = rel_x;
  if(rel_x < 0){
    image[0] = -1;
    cell2_x = D_N_CELLS_X + rel_x;
  }else if(rel_x >= D_N_CELLS_X){
    image[0] = 1;
    cell2_x = rel_x - D_N_CELLS_X;
  }
  // for y
  image[1] = 0;
  const int rel_y = cell1_y + dy;
  int cell2_y = rel_y;
  if(rel_y < 0){
    image[1] = -1;
    cell2_y = D_N_CELLS_Y + rel_y;
  }else if(rel_y >= D_N_CELLS_Y){
    image[1] = 1;
    cell2_y = rel_y - D_N_CELLS_Y;
  }
  //const int head_cell_id = d_idx_xy_head_cell[get_column_id_from_crd(cell2_x, cell2_y)];
  
  for(int i_img = 0; i_img < 3; i_img++){
    image[2] = i_img-1;
    if(first_uni_z[i_img] < 0) continue;
    
    const int first_uni_id = get_uni_id_from_crd(cell2_x, cell2_y,
						 first_uni_z[i_img]);
    const int last_uni_id = get_uni_id_from_crd(cell2_x, cell2_y,
						last_uni_z[i_img]);
    
    const int first_cell = d_uni2cell_z[first_uni_id].x;
    const int last_cell = d_uni2cell_z[last_uni_id].y;
    //if(first_cell > last_cell || first_cell < 0|| last_cell >= D_N_CELLS || first_cell < 0 || last_cell >= D_N_CELLS){
    //if(cell1_id == 0){
    //printf("DGB!! cell:%d/%d img:%d uni:%d-%d cell:%d-%d\n", 
    //cell1_id, D_N_CELLS, 
    //i_img,
    //first_uni_id, last_uni_id,
    //first_cell, last_cell);
    //}
    //}
    
    for(int cell2_id = first_cell;
	cell2_id <= last_cell; cell2_id++){      
      int n_atom_cell2 = 0;
      const int tail = (cell2_id+1) * N_ATOM_CELL;
      for(int at2 = tail-N_ATOM_CELL; at2 < tail; at2++){
	if(d_atomids[at2] >= 0) n_atom_cell2++;
	else break;
      }
      //if (cell1_id == 0){
      //printf("cell: %d-%d\n", cell1_id, cell2_id);
      //}
      if(added_col >=  D_MAX_N_CELL_PAIRS_PER_COLUMN){
	printf("Index exceeds the maximum value. %d / %d",added_col, D_MAX_N_CELL_PAIRS_PER_COLUMN);
      }
      if(check_valid_pair(cell1_id, cell2_id)){
	//d_cell_pairs[idx_cell_pair_head + added_col] = 
	const int cp_idx_cell = atomicAdd(&d_n_cell_pairs[cell1_id], 1);
	d_cell_pairs[idx_cell_pair_head + cp_idx_cell] = 
	  get_new_cell_pair(cell1_id, cell2_id,
			    cell1_id_in_block,
			    image,
			    d_atomids, d_nb15off_orig,
			    d_nb15off,
			    s_nb15off, n_atom_cell2);
	added_col++;
	//if (cell1_id==0){
	//printf("nb_col:%d i_col:%d i_img:%d cp_idx:%d mask: %d %d\n",
	//neighbor_col_id,
	//i_col, i_img,
	//idx_cell_pair_head + added_col-1,
	//d_cell_pairs[idx_cell_pair_head + added_col-1].pair_mask[0],
	//d_cell_pairs[idx_cell_pair_head + added_col-1].pair_mask[1]);
	//}
	
      }//else{
      //	  if (threadIdx.x + blockIdx.x * blockDim.x == 0){
      //printf("pair invalid\n");
      //}
      //}
    }
  }
}


__global__ void kernel_init_cell_pairs(CellPair* d_cell_pairs){
  const int cp_id = threadIdx.x + blockDim.x * blockIdx.x;  
  if (cp_id >= D_MAX_N_CELL_PAIRS) return;
  d_cell_pairs[cp_id].cell_id1 = -1;
  d_cell_pairs[cp_id].cell_id2 = -1;
  d_cell_pairs[cp_id].image = 0;
  for(int i = 0; i < N_BITMASK; i++){
    d_cell_pairs[cp_id].pair_mask[i] = ~0;
  }
}
__global__ void kernel_init_uni2cell(const int n_uni, int2* d_uni2cell_z){
  const int uni_id = threadIdx.x + blockDim.x * blockIdx.x;
  if(uni_id >= n_uni) return;
  d_uni2cell_z[uni_id].x = MAX_INT;
  d_uni2cell_z[uni_id].y = -1;
}

extern "C" int cuda_enumerate_cell_pairs(const int n_cells, const int n_uni,
					 const int max_n_cell_pairs,
					 const int n_neighbor_col){

  cudaStream_t stream1;
  cudaStream_t stream2;
  cudaStreamCreate(&stream1);
  cudaStreamCreate(&stream2);

  const int blocks1 = (n_uni + REORDER_THREADS-1) / REORDER_THREADS;  
  kernel_init_uni2cell<<<blocks1, REORDER_THREADS, 0, stream1>>>(n_uni, d_uni2cell_z);
  const int blocks2 = (n_cells + REORDER_THREADS-1) / REORDER_THREADS;  
  kernel_set_uniform_grid<<<blocks2, REORDER_THREADS, 0, stream1>>>(d_crd_chg, d_uni2cell_z);
  const int blocks3 = (max_n_cell_pairs + REORDER_THREADS-1) / REORDER_THREADS;  
  kernel_init_cell_pairs<<<blocks3, REORDER_THREADS, 0, stream2>>>(d_cell_pairs);

  cudaStreamSynchronize(stream1);
  cudaStreamSynchronize(stream2);
  cudaStreamDestroy(stream1);

  const int blocks4 = (n_neighbor_col * n_cells + REORDER_THREADS-1) / REORDER_THREADS;  
  kernel_enumerate_cell_pair<<<blocks4, REORDER_THREADS,0,stream2>>>(d_uni2cell_z, d_crd_chg,
								     d_idx_xy_head_cell,
								     d_atomids, d_nb15off_orig,
								     d_nb15off,
								     d_n_cell_pairs,
								     d_cell_pairs_buf);
  
  
  set_idx_head_cell_pairs<<<1, 128, 0, stream2>>>
    (d_n_cell_pairs, d_idx_head_cell_pairs);

  pack_cellpairs_array<<<blocks3, REORDER_THREADS, 0, stream2>>>
    (d_cell_pairs, d_cell_pairs_buf, d_idx_xy_head_cell, d_n_cell_pairs);

  cudaStreamDestroy(stream2);

  
  return 0;
}
		