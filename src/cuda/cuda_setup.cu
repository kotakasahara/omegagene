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
				    int max_n_cell_pairs){
  //printf("cuda_alloc_atom_info\n");
  // max_n_atoms ... maximum number of atoms for each grid cell
  HANDLE_ERROR( cudaMalloc((void**)&d_crd_chg,
			   n_atom_array * sizeof(real4)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_crd,
			   n_atom_array * 3 * sizeof(real_pw)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_charge_orig,
			   n_atoms * sizeof(real_pw) ));
  HANDLE_ERROR( cudaMalloc((void**)&d_atominfo,
			   n_atom_array * sizeof(int2) ));
  HANDLE_ERROR( cudaMalloc((void**)&d_atomids,
			   n_atom_array * sizeof(int)) );
  //  HANDLE_ERROR( cudaMalloc((void**)&d_atomids_rev,
  //			   n_atoms * sizeof(int)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_atomtype_orig,
			   n_atoms * sizeof(int) ));
  HANDLE_ERROR( cudaMalloc((void**)&d_cell_pairs,
			   max_n_cell_pairs * sizeof(CellPair)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_idx_head_cell_pairs,
			   (max_n_cells+1) * sizeof(int)) );
  //HANDLE_ERROR( cudaMalloc((void**)&d_grid_atom_index,
  //(max_n_ + 1) * sizeof(int)) );
  HANDLE_ERROR( cudaMalloc((void**)&d_energy,
			   N_MULTI_WORK * 2 * sizeof(real_fc) ) );
  HANDLE_ERROR( cudaMalloc((void**)&d_work,
			   N_MULTI_WORK * n_atom_array * 3 * sizeof(real_fc) ) );
  //HANDLE_ERROR( cudaMalloc((void**)&d_work_orig,
  //n_atom_array * 3 * sizeof(real_fc) ) );
  return 0;
}

extern "C" int cuda_free_atom_info(){
  //printf("cuda_free_device_atom_info\n");
  HANDLE_ERROR( cudaFree(d_crd_chg) );
  HANDLE_ERROR( cudaFree(d_crd) );
  HANDLE_ERROR( cudaFree(d_atomids) );
  //  HANDLE_ERROR( cudaFree(d_atomids_rev) );
  HANDLE_ERROR( cudaFree(d_charge_orig) );
  HANDLE_ERROR( cudaFree(d_atominfo) );
  HANDLE_ERROR( cudaFree(d_atomtype_orig) );
  HANDLE_ERROR( cudaFree(d_cell_pairs) );
  HANDLE_ERROR( cudaFree(d_idx_head_cell_pairs) );
  HANDLE_ERROR( cudaFree(d_energy) );
  HANDLE_ERROR( cudaFree(d_work) );
  //HANDLE_ERROR( cudaFree(d_work_orig) );
  return 0;
}

extern "C" int cuda_memcpy_htod_cell_pairs(CellPair*& h_cell_pairs,
					   int*& h_idx_head_cell_pairs,
					   int n_cell_pairs,
					   int n_cells){
  //printf("cuda_memcpy_htod_cell_pairs\n");
  HANDLE_ERROR(cudaMemcpy(d_cell_pairs, h_cell_pairs,
			  n_cell_pairs * sizeof(CellPair),
			  cudaMemcpyHostToDevice));
  HANDLE_ERROR(cudaMemcpy(d_idx_head_cell_pairs,
			  h_idx_head_cell_pairs,
			  (n_cells+1) * sizeof(int),
			  cudaMemcpyHostToDevice));
  return 0;
}
extern "C" int cuda_memcpy_htod_atomids(int*& h_atomids,
					int n_atom_array){
  HANDLE_ERROR(cudaMemcpy(d_atomids, h_atomids,
			  n_atom_array * sizeof(int),
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

extern "C" int cuda_set_pbc(real_pw* l){
  HANDLE_ERROR( cudaMemcpyToSymbol(PBC_L,
				   l, sizeof(real_pw) * 3) );
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
extern "C" int cuda_set_cell_constant(int n_cells, int n_cell_pairs, int n_atom_array){
  HANDLE_ERROR( cudaMemcpyToSymbol(D_N_CELL_PAIRS,
				   &n_cell_pairs,
				   sizeof(int) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_N_CELLS,
				   &n_cells,
				   sizeof(int) ) );
  HANDLE_ERROR( cudaMemcpyToSymbol(D_N_ATOM_ARRAY,
				   &n_atom_array,
				   sizeof(int) ) );
  return 0;
}

// cuda_set_constant
//   called only onece at the beginning of simulation
extern "C" int cuda_set_constant(int n_atoms, real_pw cutoff, int n_atomtypes){
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
  return 0;
}

//extern "C" int cuda_alloc_set_nb15off(int* h_nb15off1,
//			      int* h_nb15off2,
//				      int n_atoms){
  //cudaMalloc
//  HANDLE_ERROR( cudaMalloc((void**)&d_nb15off1_orig,
//			   n_atoms *  sizeof(int)) );
//  HANDLE_ERROR( cudaMalloc((void**)&d_nb15off2_orig,
//			   n_atoms *  sizeof(int)) );

  // cudaMemcpy
//  HANDLE_ERROR( cudaMemcpy( d_nb15off1_orig, h_nb15off1,
//			    sizeof(int) * n_atoms,
//			    cudaMemcpyHostToDevice) );
//  HANDLE_ERROR( cudaMemcpy( d_nb15off2_orig, h_nb15off2,
//			    sizeof(int) * n_atoms,
//			    cudaMemcpyHostToDevice) );
//  return 0;
//}
//extern "C" int cuda_free_nb15off(){
//  HANDLE_ERROR( cudaFree(d_nb15off1_orig) );
//  HANDLE_ERROR( cudaFree(d_nb15off2_orig) );
//  return 0;
//}

extern "C" int cuda_alloc_set_lj_params(real_pw* h_lj_6term,
					real_pw* h_lj_12term,
					int n_lj_types){
  //printf("threads : %d\n", PW_THREADS);
  unsigned int size_lj_matrix = sizeof(real_pw) * n_lj_types * n_lj_types;
  // cudaMalloc
  HANDLE_ERROR( cudaMalloc( (void**)&d_lj_6term,
			    size_lj_matrix ) );
  HANDLE_ERROR( cudaMalloc( (void**)&d_lj_12term,
			    size_lj_matrix ) );
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
  
  return 0;
}

extern "C" int cuda_free_lj_params(){
  //printf("cuda_free_lj_param\n");
  //cudaUnbindTexture(tex_lj_6term);
  //cudaUnbindTexture(tex_lj_12term);
  HANDLE_ERROR( cudaFree(d_lj_6term) );
  HANDLE_ERROR( cudaFree(d_lj_12term) );
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
  printf("test\n");
 return 0;
}

extern "C" int cuda_hostalloc_cell_info(CellPair*& h_cell_pairs, 
					int*& h_idx_head_cell_pairs,
					int max_n_cell_pairs,
					int max_n_cells){
  printf("cuda_hostalloc_cell_info cu\n");
  HANDLE_ERROR( cudaHostAlloc( (void**)&h_cell_pairs,
			       max_n_cell_pairs * sizeof(CellPair),
			       cudaHostAllocDefault));
  HANDLE_ERROR( cudaHostAlloc( (void**)&h_idx_head_cell_pairs,
			       (max_n_cells) * sizeof(int),
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
extern "C" int cuda_hostfree_cell_info(CellPair* h_cell_pairs, int* h_idx_head_cell_pairs){
  HANDLE_ERROR( cudaFreeHost(h_cell_pairs) );
  HANDLE_ERROR( cudaFreeHost(h_idx_head_cell_pairs) );
  return 0;
}
__global__ void kernel_set_atominfo(const int* d_atomids,
				    const int* d_atomtype_orig,
				    const real_pw* d_crd,
				    const real_pw* d_charge_orig,
				    int2* d_atominfo,
				    real4* d_crd_chg){
  int atomid = threadIdx.x + blockIdx.x * blockDim.x;
  if(atomid < D_N_ATOM_ARRAY){
    d_atominfo[atomid].x = d_atomids[atomid];
    if(d_atomids[atomid] >= 0){
      d_atominfo[atomid].y = d_atomtype_orig[d_atomids[atomid]];
      d_crd_chg[atomid].w = d_charge_orig[d_atomids[atomid]];
    }else{
      d_atominfo[atomid].y = 0;
      d_crd_chg[atomid].w = 0.0;
    }
    //d_atominfo[atomid].z = d_nb15off1_orig[d_atomids[atomid]];
    //d_atominfo[atomid].w = d_nb15off2_orig[d_atomids[atomid]];
    d_crd_chg[atomid].x = d_crd[atomid*3+0];
    d_crd_chg[atomid].y = d_crd[atomid*3+1];
    d_crd_chg[atomid].z = d_crd[atomid*3+2];
  }
  
}
/*
__global__ void kernel_set_atomids_rev(const int* d_atomids, int* d_atomids_rev){
  int atomid = threadIdx.x + blockIdx.x * blockDim.x;
  if(atomid < D_N_ATOMS){
    d_atomids_rev[d_atomids[atomid]] = atomid;
  }
}
*/
extern "C" int cuda_set_atominfo(int n_atom_array){
  int blocks = (n_atom_array + REORDER_THREADS-1) / REORDER_THREADS;
  kernel_set_atominfo<<<blocks, REORDER_THREADS>>>(d_atomids,
						   d_atomtype_orig,
						   d_crd,
						   d_charge_orig,
						   d_atominfo,
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
			      const int*bitmask){
  int bit_pos = atom_idx2 * N_ATOM_CELL + atom_idx1;
  int mask_id =  bit_pos / 32;
  int interact_pos =  bit_pos % 32;
  int interact = 1 << interact_pos;
  return (bitmask[mask_id] & interact) == interact;
}

__device__ real_pw check_15off(const int atomid1, const int atomid2,
			       const int tmask_a1,
			       const int tmask_a2
				 ){
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
__device__ int cal_pair(real_pw& w1,
			real_pw& w2,
			real_pw& w3,
			real_pw& ene_vdw,
			real_pw& ene_ele,
			const real4& crd_chg1, 
			const real4& crd_chg2, 
			const int& atomtype1,
			const int& atomtype2,
			const real_pw* __restrict__ d_lj_6term,
			const real_pw* __restrict__ d_lj_12term,
			const int a1, const int a2
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

  if(r12 >= D_CUTOFF) {return 1;}
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
  return 0;
}

__global__ void kernel_pairwise_ljzd(const real4* d_crd_chg,
				     const CellPair* d_cell_pairs,
				     const int* d_idx_head_cell_pairs,
				     const int2* d_atominfo,
				     const real_pw* __restrict__ d_lj_6term,
				     const real_pw* __restrict__ d_lj_12term,
				     real_fc* d_energy, real_fc* d_work,
				     int* d_pairs, real_pw* d_realbuf,
				     const int offset_cells,
				     const int n_cells){
  real_fc ene_vdw = 0.0;
  real_fc ene_ele = 0.0;

  const int global_threadIdx = blockDim.x * blockIdx.x + threadIdx.x;

  const int c1 = global_threadIdx >> 5;
  const int warpIdx = threadIdx.x >> 5;  
  if(c1 >= n_cells){ return; }
  const int laneIdx = global_threadIdx & 31;
  const int n_loops = (d_idx_head_cell_pairs[c1+1] - d_idx_head_cell_pairs[c1])*2;
  const int ene_index_offset = global_threadIdx%N_MULTI_WORK;  
  
  real_fc work_c1[3] = {0.0, 0.0, 0.0};

  const int atom_idx1 = (laneIdx & 7);  // laneIdx%8
  const int a1 = c1 * D_N_ATOM_CELL + atom_idx1;
  
  __shared__ real4 crd_chg1[D_N_ATOM_CELL * (PW_THREADS >> 5)];
  __shared__ int2  atominfo1[D_N_ATOM_CELL * (PW_THREADS >> 5)];

  //__shared__ real4 crd_chg2[D_N_ATOM_CELL * PW_THREADS / 32];
  const int sharedmem_idx  = D_N_ATOM_CELL * warpIdx + atom_idx1;
  if(laneIdx<D_N_ATOM_CELL){
    crd_chg1[sharedmem_idx] = d_crd_chg[c1*D_N_ATOM_CELL + laneIdx];
    atominfo1[sharedmem_idx] = d_atominfo[c1*D_N_ATOM_CELL + laneIdx];
  }
  __syncthreads();

  CellPair cellpair;
  int cp;
  for(int loopIdx=0; loopIdx < n_loops; loopIdx++){
    if(loopIdx%2 == 0){
      cp = d_idx_head_cell_pairs[c1] + (loopIdx >> 1);
      if(cp >= D_N_CELL_PAIRS) break;
      cellpair = d_cell_pairs[cp];
    }
    
    //CellPair cellpair;
    //if(laneIdx==0){
    //cellpair = cellpair;
    //}
    //for(int i=1; i < 32; i*=2){
    //cellpair = __shfl_up(cellpair, i);
    //}

    //atomicAdd(&d_pairs[3],1);
    const int c2 = cellpair.cell_id2;

    // atom_idx ... index in cell, 0-7
    const int atom_idx2 = (laneIdx / 8)  + 4 * (loopIdx % 2);  // laneIdx/8 + 4*(warpIdx%2)

    // remove 1-2, 1-3, 1-4 pairs
    const int a2 = c2 * D_N_ATOM_CELL + atom_idx2;
    //real4 crd_chg2;
    //int2 atominfo2;
    //if(atom_idx1 == 0){
    real4 crd_chg2 = d_crd_chg[a2];
    const int2 atominfo2 = d_atominfo[a2];
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

    if(!check_15off64(atom_idx1, atom_idx2, cellpair.pair_mask)){
      //if(d_atominfo[a1].x == -1 || d_atominfo[a2].x == -1) continue;
      cal_pair(w1, w2, w3, cur_ene_vdw, cur_ene_ele,
	       crd_chg1[sharedmem_idx],
	       crd_chg2,
	       atominfo1[sharedmem_idx].y,
	       atominfo2.y,
	       d_lj_6term, d_lj_12term,
	       atominfo1[sharedmem_idx].x,
	       atominfo2.x);
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
  //int a1 = c1 * D_N_ATOM_CELL;
  //for(int atom_idx1 = 0; atom_idx1 < D_N_ATOM_CELL; atom_idx1++, a1++){
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

extern "C" int cuda_pairwise_ljzd(const int offset_cellpairs, const int n_cal_cellpairs,
				  const int offset_cells,     const int n_cal_cells){
  // test
  //real_pw w1=0.0, w2=0.0, w3=0.0;
  //real_pw cur_ene_ele=0.0;
  //real_pw cur_ene_vdw=0.0;
  //printf("test pair\n");
  
  //cal_pair(w1,w2,w3, cur_ene_vdw, cur_ene_ele, a1, a2,
  //d_crd_chg, d_atominfo, crd2,
  //d_lj_6term, d_lj_12term);

  cudaStreamCreate(&stream_pair_home);

  //test
  int *d_pairs;
  HANDLE_ERROR(cudaMalloc((void**)&d_pairs, sizeof(int)*5));
  HANDLE_ERROR(cudaMemset(d_pairs, 0, sizeof(int)*5));
  int h_pairs[5] = {0,0,0,0,0};
  real_pw *d_realbuf;
  HANDLE_ERROR(cudaMalloc((void**)&d_realbuf, sizeof(real_pw)*4));
  HANDLE_ERROR(cudaMemset(d_realbuf, 0.0, sizeof(real_pw)*4));
  real_pw h_realbuf[4] = {0,0,0,0};
  //HANDLE_ERROR(cudaMemcpy(d_pairs, h_pairs, sizeof(int)*5, cudaMemcpyHostToDevice));
  //HANDLE_ERROR(cudaMemcpy(d_realbuf, h_realbuf, sizeof(real)*4, cudaMemcpyHostToDevice));
  //

  int blocks = (n_cal_cells+PW_THREADS/32-1) / (PW_THREADS/32);
  kernel_pairwise_ljzd<<<blocks, PW_THREADS,
    0, stream_pair_home>>>(d_crd_chg,
			   d_cell_pairs,
			   d_idx_head_cell_pairs,
			   d_atominfo,
			   d_lj_6term, d_lj_12term,
			   d_energy, d_work, d_pairs, d_realbuf,
			   offset_cells, n_cal_cells);
  //HANDLE_ERROR(cudaMemcpy(h_pairs, d_pairs, sizeof(int)*5, cudaMemcpyDeviceToHost));
  //HANDLE_ERROR(cudaMemcpy(h_realbuf, d_realbuf, sizeof(real)*4, cudaMemcpyDeviceToHost));
  //printf("n 15 pairs: nb:%d / interact:%d / all:%d / cellpairs:%d / incut:%d\n", h_pairs[0],h_pairs[1],h_pairs[2], h_pairs[3], h_pairs[4]); 
  //printf("sum_dist: %f  charge: %f\n",h_realbuf[0], h_realbuf[1]);
  HANDLE_ERROR(cudaFree(d_pairs));
  HANDLE_ERROR(cudaFree(d_realbuf));

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

