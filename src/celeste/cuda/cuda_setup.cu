#include "cuda_setup.h"

__device__ __inline__ double shfl_xor(double value, int const lane, int const warpsize) {
    return __hiloint2double(__shfl_xor(__double2hiint(value), lane, warpsize),
                            __shfl_xor(__double2loint(value), lane, warpsize));
}

__device__ double atomicAdd2(double *address, double val) {
  unsigned long long int *address_as_ull = (unsigned long long int *)address;
  unsigned long long int  old            = *address_as_ull, assumed;
  do {
    assumed = old;
    old     = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}

extern "C" int cuda_alloc_atom_info(int in_max_n_atoms_exbox,
                                    // int in_max_n_atom_array,
                                    int in_max_n_cells,
                                    int in_max_n_cell_pairs,
                                    int in_n_columns) {
  printf("cuda_alloc_atom_info\n");
  max_n_atoms_exbox = in_max_n_atoms_exbox;
  // max_n_atom_array = in_max_n_atom_array;
  max_n_cell_pairs = in_max_n_cell_pairs;
  max_n_cells      = in_max_n_cells;
  HANDLE_ERROR(cudaMalloc((void **)&d_crd_chg, max_n_atom_array * sizeof(real4)));
  HANDLE_ERROR(cudaMalloc((void **)&d_cell_z, max_n_cells * sizeof(real2)));
  HANDLE_ERROR(cudaMalloc((void **)&d_crd, max_n_atom_array * 3 * sizeof(real_pw)));
  HANDLE_ERROR(cudaMalloc((void **)&d_charge_orig, max_n_atoms_exbox * sizeof(real_pw)));
  HANDLE_ERROR(cudaMalloc((void **)&d_atomtype, max_n_atom_array * sizeof(int)));
  HANDLE_ERROR(cudaMalloc((void **)&d_atomids, max_n_atom_array * sizeof(int)));
  HANDLE_ERROR(cudaMalloc((void **)&d_atomids_rev, max_n_atoms_exbox * sizeof(int)));
  HANDLE_ERROR(cudaMalloc((void **)&d_atomtype_orig, max_n_atoms_exbox * sizeof(int)));
  HANDLE_ERROR(cudaMalloc((void **)&d_cell_pairs, max_n_cell_pairs * sizeof(CellPair)));
  HANDLE_ERROR(cudaMalloc((void **)&d_cell_pairs_buf, max_n_cell_pairs * sizeof(CellPair)));
  HANDLE_ERROR(cudaMalloc((void **)&d_idx_head_cell_pairs, (in_max_n_cells + 1) * sizeof(int)));
  HANDLE_ERROR(cudaMalloc((void **)&d_idx_cell_column, (in_max_n_cells) * sizeof(int)));
  HANDLE_ERROR(cudaHostAlloc((void **)&h_idx_cell_column, in_max_n_cells * sizeof(int), cudaHostAllocDefault));
  // HANDLE_ERROR( cudaMalloc((void**)&d_cell_pair_removed,
  //(in_max_n_cells+1) * sizeof(int)) );
  HANDLE_ERROR(cudaMalloc((void **)&d_n_cell_pairs, (max_n_cells) * sizeof(int)));
  // HANDLE_ERROR( cudaMalloc((void**)&d_grid_atom_index,
  //(max_n_ + 1) * sizeof(int)) );
  HANDLE_ERROR(cudaMalloc((void **)&d_energy, N_MULTI_WORK * 2 * sizeof(real_fc)));
  HANDLE_ERROR(cudaMalloc((void **)&d_work, N_MULTI_WORK * max_n_atom_array * 3 * sizeof(real_fc)));
  HANDLE_ERROR(cudaMalloc((void **)&d_idx_xy_head_cell, in_n_columns * sizeof(int)));
  HANDLE_ERROR(cudaMemcpyToSymbol(D_MAX_N_CELL_PAIRS, &in_max_n_cell_pairs, sizeof(int)));
  
  return 0;
}

extern "C" int cuda_free_atom_info() {
  // printf("cuda_free_device_atom_info\n");
  HANDLE_ERROR(cudaFree(d_crd_chg));
  HANDLE_ERROR(cudaFree(d_cell_z));
  HANDLE_ERROR(cudaFree(d_crd));
  HANDLE_ERROR(cudaFree(d_atomids));
  HANDLE_ERROR(cudaFree(d_atomids_rev));
  HANDLE_ERROR(cudaFree(d_charge_orig));
  HANDLE_ERROR(cudaFree(d_atomtype));
  HANDLE_ERROR(cudaFree(d_atomtype_orig));
  HANDLE_ERROR(cudaFree(d_cell_pairs));
  HANDLE_ERROR(cudaFree(d_cell_pairs_buf));
  HANDLE_ERROR(cudaFree(d_idx_head_cell_pairs));
  HANDLE_ERROR(cudaFree(d_idx_cell_column));
  HANDLE_ERROR(cudaFreeHost(h_idx_cell_column));
  // HANDLE_ERROR( cudaFree(d_cell_pair_removed) );
  HANDLE_ERROR(cudaFree(d_n_cell_pairs));
  HANDLE_ERROR(cudaFree(d_energy));
  HANDLE_ERROR(cudaFree(d_work));
  HANDLE_ERROR(cudaFree(d_idx_xy_head_cell));
  // HANDLE_ERROR( cudaFree(d_uni2cell_z));
  // HANDLE_ERROR( cudaFree(d_work_orig) );
  return 0;
}

extern "C" int cuda_memcpy_htod_atomids(int *&h_atomids, int *&h_idx_xy_head_cell) {
    HANDLE_ERROR(cudaMemset(d_atomids, -1, sizeof(int) * max_n_atom_array));
    HANDLE_ERROR(cudaMemcpy(d_atomids, h_atomids, n_atom_array * sizeof(int), cudaMemcpyHostToDevice));
    HANDLE_ERROR(
        cudaMemcpy(d_idx_xy_head_cell, h_idx_xy_head_cell, (n_columns + 1) * sizeof(int), cudaMemcpyHostToDevice));
    return 0;
}

// cuda_memcpy_htod_atom_info
//   Arrays of charges and atomtypes of all atoms in the process are sent to
//   the device.
extern "C" int cuda_memcpy_htod_atom_info(real_pw *&h_charge_orig, int *&h_atomtype_orig) {

    HANDLE_ERROR(cudaMemcpy(d_charge_orig, h_charge_orig, n_atoms_system * sizeof(real_pw), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(d_atomtype_orig, h_atomtype_orig, n_atoms_system * sizeof(int), cudaMemcpyHostToDevice));
    return 0;
}

// cuda_memcpy_htod_crd
//  Sending nsgrid.crd to device
extern "C" int cuda_memcpy_htod_crd(real_pw *&h_crd) {
    HANDLE_ERROR(cudaMemcpy(d_crd, h_crd, n_atom_array * 3 * sizeof(real_pw), cudaMemcpyHostToDevice));
    return 0;
}

extern "C" int cuda_set_pbc(real_pw *l, real_pw *lb) {
    HANDLE_ERROR(cudaMemcpyToSymbol(PBC_L, l, sizeof(real_pw) * 3));
    HANDLE_ERROR(cudaMemcpyToSymbol(PBC_LOWER_BOUND, lb, sizeof(real_pw) * 3));
    return 0;
}

extern "C" int cuda_zerodipole_constant(real_pw zcore, real_pw bcoeff, real_pw fcoeff) {
    HANDLE_ERROR(cudaMemcpyToSymbol(D_ZCORE, &zcore, sizeof(real_pw)));
    HANDLE_ERROR(cudaMemcpyToSymbol(D_BCOEFF, &bcoeff, sizeof(real_pw)));
    HANDLE_ERROR(cudaMemcpyToSymbol(D_FCOEFF, &fcoeff, sizeof(real_pw)));
    return 0;
}

// cuda_set_cell_constant
//  These constants are updated when the cell grid is updated
extern "C" int cuda_set_cell_constant(const int      in_n_cells,
                                      const int      in_n_atoms_exbox,
                                      const int      in_n_atom_array,
                                      const int *    in_n_cells_xyz,
                                      const int      in_n_columns,
                                      const real_pw *in_l_cell_xyz,
                                      const int *    in_n_neighbor_xyz) {
    n_atoms_exbox = in_n_atoms_exbox;
    n_cells       = in_n_cells;
    n_atom_array  = in_n_atom_array;
    n_columns     = in_n_columns;
    HANDLE_ERROR(cudaMemcpyToSymbol(D_N_CELLS, &n_cells, sizeof(int)));
    HANDLE_ERROR(cudaMemcpyToSymbol(D_N_ATOM_ARRAY, &n_atom_array, sizeof(int)));
    HANDLE_ERROR(cudaMemcpyToSymbol(D_N_CELLS_XYZ, in_n_cells_xyz, sizeof(int) * 3));
    HANDLE_ERROR(cudaMemcpyToSymbol(D_N_COLUMNS, &n_columns, sizeof(int)));
    HANDLE_ERROR(cudaMemcpyToSymbol(D_L_CELL_XYZ, in_l_cell_xyz, sizeof(real_pw) * 3));
    HANDLE_ERROR(cudaMemcpyToSymbol(D_N_NEIGHBOR_XYZ, in_n_neighbor_xyz, sizeof(int) * 3));
    const int n_neighbor_col = (in_n_neighbor_xyz[0] * 2 + 1) * (in_n_neighbor_xyz[1] * 2 + 1);
    HANDLE_ERROR(cudaMemcpyToSymbol(D_N_NEIGHBOR_COL, &n_neighbor_col, sizeof(int)));
    const int max_n_cell_pairs_per_cell = max_n_cell_pairs / n_cells;
    HANDLE_ERROR(cudaMemcpyToSymbol(D_MAX_N_CELL_PAIRS_PER_CELL, &max_n_cell_pairs_per_cell, sizeof(int)));
    return 0;
}

// cuda_set_constant
//   called only onece at the beginning of simulation
extern "C" int cuda_set_constant(real_pw cutoff, real_pw cutoff_pairlist, int n_atomtypes) {
    real_pw tmp_charge_coeff = (real_pw)332.06378; // CelesteObject::CHARGE_COEFF;
    HANDLE_ERROR(cudaMemcpyToSymbol(D_CHARGE_COEFF, &tmp_charge_coeff, sizeof(real_pw)));

    HANDLE_ERROR(cudaMemcpyToSymbol(D_N_ATOMTYPES, &n_atomtypes, sizeof(int)));
    HANDLE_ERROR(cudaMemcpyToSymbol(D_CUTOFF, &cutoff, sizeof(real_pw)));
    HANDLE_ERROR(cudaMemcpyToSymbol(D_CUTOFF_PAIRLIST, &cutoff_pairlist, sizeof(real_pw)));
    const real_pw cutoff_pairlist_2 = cutoff_pairlist * cutoff_pairlist;
    HANDLE_ERROR(cudaMemcpyToSymbol(D_CUTOFF_PAIRLIST_2, &cutoff_pairlist_2, sizeof(real_pw)));
    return 0;
}

extern "C" int cuda_alloc_set_lj_params(real_pw * h_lj_6term,
                                        real_pw * h_lj_12term,
                                        int       n_lj_types,
                                        int *     h_nb15off,
                                        const int in_max_n_nb15off) {
    // printf("threads : %d\n", PW_THREADS);
    printf("cuda_alloc_set_lj_params\n");
    const unsigned int size_lj_matrix = sizeof(real_pw) * n_lj_types * n_lj_types;
    // cudaMalloc
    HANDLE_ERROR(cudaMalloc((void **)&d_lj_6term, size_lj_matrix));
    HANDLE_ERROR(cudaMalloc((void **)&d_lj_12term, size_lj_matrix));
    max_n_nb15off = in_max_n_nb15off;
    HANDLE_ERROR(cudaMemcpyToSymbol(D_MAX_N_NB15OFF, &max_n_nb15off, sizeof(int)));

    const unsigned int size_nb15off_orig = sizeof(int) * max_n_nb15off * max_n_atoms_exbox;
    HANDLE_ERROR(cudaMalloc((void **)&d_nb15off_orig, size_nb15off_orig));
    size_nb15off = max_n_nb15off * max_n_atom_array;
    HANDLE_ERROR(cudaMalloc((void **)&d_nb15off, sizeof(int) * size_nb15off));
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
    HANDLE_ERROR(cudaMemcpy(d_lj_6term, h_lj_6term, size_lj_matrix, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(d_lj_12term, h_lj_12term, size_lj_matrix, cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpy(d_nb15off_orig, h_nb15off, size_nb15off_orig, cudaMemcpyHostToDevice));

    return 0;
}

extern "C" int cuda_free_lj_params() {
    // printf("cuda_free_lj_param\n");
    // cudaUnbindTexture(tex_lj_6term);
    // cudaUnbindTexture(tex_lj_12term);
    HANDLE_ERROR(cudaFree(d_lj_6term));
    HANDLE_ERROR(cudaFree(d_lj_12term));
    HANDLE_ERROR(cudaFree(d_nb15off_orig));
    HANDLE_ERROR(cudaFree(d_nb15off));
    return 0;
}
// cuda_hostalloc_atom_type_charge
extern "C" int cuda_hostalloc_atom_type_charge(int *&h_atom_type, real_pw *&h_charge, const int in_n_atoms_system) {
    n_atoms_system = in_n_atoms_system;
    HANDLE_ERROR(cudaMemcpyToSymbol(D_N_ATOMS_SYSTEM, &in_n_atoms_system, sizeof(int)));
    printf("hostalloc atom_type_charge cu %d\n", in_n_atoms_system);
    HANDLE_ERROR(cudaHostAlloc((void **)&h_atom_type, n_atoms_system * sizeof(int), cudaHostAllocDefault));
    HANDLE_ERROR(cudaHostAlloc((void **)&h_charge, n_atoms_system * sizeof(real_pw), cudaHostAllocDefault));
    return 0;
}

// cuda_hostalloc_atom_info
//   Allocation for MiniCell members
extern "C" int cuda_hostalloc_atom_info(real_pw *&h_crd,
                                        int *&    h_atomids,
                                        real_fc *&h_work,
                                        real_fc *&h_energy,
                                        int       in_max_n_atom_array) {
    max_n_atom_array = in_max_n_atom_array;
    printf("hostalloc_atom_info %d\n", max_n_atom_array);
    HANDLE_ERROR(cudaHostAlloc((void **)&h_crd, max_n_atom_array * 3 * sizeof(real_pw), cudaHostAllocDefault));
    HANDLE_ERROR(cudaHostAlloc((void **)&h_atomids, max_n_atom_array * sizeof(int), cudaHostAllocDefault));
    HANDLE_ERROR(cudaHostAlloc((void **)&h_work, max_n_atom_array * 3 * sizeof(real_fc), cudaHostAllocDefault));
    HANDLE_ERROR(cudaHostAlloc((void **)&h_energy, 2 * sizeof(real_fc), cudaHostAllocDefault));
    return 0;
}

extern "C" int cuda_hostalloc_cell_info(int *&h_idx_xy_head_cell, int n_columns) {
    printf("cuda_hostalloc_cell_info cu\n");
    HANDLE_ERROR(cudaHostAlloc((void **)&h_idx_xy_head_cell, (n_columns) * sizeof(int), cudaHostAllocDefault));

    return 0;
}

extern "C" int cuda_hostalloc_cellpair_info(CellPair *&h_cell_pairs,
                                            int *&     h_idx_head_cell_pairs,
                                            int *&     h_n_cells_z,
                                            int        max_n_cell_pairs,
                                            int        max_n_cells,
                                            int        n_columns) {
    printf("cuda_hostalloc_cellpair_info cu\n");
    HANDLE_ERROR(cudaHostAlloc((void **)&h_cell_pairs, max_n_cell_pairs * sizeof(CellPair), cudaHostAllocDefault));
    HANDLE_ERROR(cudaHostAlloc((void **)&h_idx_head_cell_pairs, (max_n_cells) * sizeof(int), cudaHostAllocDefault));
    HANDLE_ERROR(cudaHostAlloc((void **)&h_n_cells_z, (n_columns) * sizeof(int), cudaHostAllocDefault));
    return 0;
}
extern "C" int cuda_hostfree_cellpair_info(CellPair *h_cell_pairs, int *h_idx_head_cell_pairs, int *&h_n_cells_z) {
    HANDLE_ERROR(cudaFreeHost(h_cell_pairs));
    HANDLE_ERROR(cudaFreeHost(h_idx_head_cell_pairs));
    HANDLE_ERROR(cudaFreeHost(h_n_cells_z));
}

extern "C" int cuda_hostfree_atom_type_charge(int *h_atom_type, real_pw *h_charge) {
    HANDLE_ERROR(cudaFreeHost(h_atom_type));
    HANDLE_ERROR(cudaFreeHost(h_charge));
    return 0;
}

extern "C" int cuda_hostfree_atom_info(real_pw *h_crd, int *h_atomids, real_fc *&h_work, real_fc *&h_energy) {
    HANDLE_ERROR(cudaFreeHost(h_crd));
    HANDLE_ERROR(cudaFreeHost(h_atomids));
    HANDLE_ERROR(cudaFreeHost(h_work));
    HANDLE_ERROR(cudaFreeHost(h_energy));

    return 0;
}

extern "C" int cuda_hostfree_cell_info(int *h_idx_xy_head_cell) {
    HANDLE_ERROR(cudaFreeHost(h_idx_xy_head_cell));
    return 0;
}

__global__ void kernel_set_nb15off(const int *d_atomids,
                                   const int *d_atomids_rev,
                                   const int *d_nb15off_orig,
                                   int *      d_nb15off) {
    const int g_thread_idx = threadIdx.x + blockDim.x * blockIdx.x;
    const int atomid       = g_thread_idx / D_MAX_N_NB15OFF;
    const int idx          = g_thread_idx % D_MAX_N_NB15OFF;
    if (atomid >= D_N_ATOM_ARRAY) return;
    if (d_atomids[atomid] < 0) {
        d_nb15off[g_thread_idx] = atomid;
    } else {
        const int orig = d_nb15off_orig[d_atomids[atomid] * D_MAX_N_NB15OFF + idx];
        if (orig == -1) {
            d_nb15off[g_thread_idx] = -1;
        } else {
            d_nb15off[g_thread_idx] = d_atomids_rev[orig];
        }
    }
}

__global__ void kernel_set_atominfo(const int *    d_atomids,
                                    const int *    d_atomtype_orig,
                                    const real_pw *d_charge_orig,
                                    int *          d_atomtype,
                                    real4 *        d_crd_chg,
                                    int *          d_atomids_rev) {
    int atomid = threadIdx.x + blockIdx.x * blockDim.x;
    if (atomid < D_N_ATOM_ARRAY && d_atomids[atomid] >= 0) {
        d_atomtype[atomid]               = d_atomtype_orig[d_atomids[atomid]];
        d_crd_chg[atomid].w              = d_charge_orig[d_atomids[atomid]];
        d_atomids_rev[d_atomids[atomid]] = atomid;
    }
}
__global__ void kernel_set_crd(const int *d_atomids, const real_pw *d_crd, real4 *d_crd_chg, real2 *d_cell_z) {
    int atomid = threadIdx.x + blockIdx.x * blockDim.x;
    if (atomid < D_N_ATOM_ARRAY) {
        int at_idx          = atomid * 3;
        d_crd_chg[atomid].x = d_crd[at_idx];
        d_crd_chg[atomid].y = d_crd[at_idx + 1];
        d_crd_chg[atomid].z = d_crd[at_idx + 2];
        /*if(atomid % N_ATOM_CELL == 0){
          d_cell_z[atomid/N_ATOM_CELL].x = d_crd[at_idx+2];
        }else if(atomid % N_ATOM_CELL == N_ATOM_CELL-1){
          d_cell_z[atomid/N_ATOM_CELL].y = d_crd[at_idx+2];
          }*/
    }
}

extern "C" int cuda_set_atominfo() {

    HANDLE_ERROR(cudaMemset(d_atomtype, -1, sizeof(int) * max_n_atom_array));
    HANDLE_ERROR(cudaMemset(d_atomids_rev, -1, sizeof(int) * max_n_atoms_exbox));
    HANDLE_ERROR(cudaMemset(d_nb15off, -1, sizeof(int) * size_nb15off));
    HANDLE_ERROR(cudaMemset(d_n_cell_pairs, 0, sizeof(int) * max_n_cells));
    // HANDLE_ERROR( cudaMemset(d_cell_pair_removed, 0, sizeof(int)*max_n_cells ));

    int blocks1 = (n_atom_array + REORDER_THREADS - 1) / REORDER_THREADS;
    kernel_set_atominfo<<<blocks1, REORDER_THREADS>>>(d_atomids, d_atomtype_orig, d_charge_orig, d_atomtype, d_crd_chg,
                                                      d_atomids_rev);
    int blocks2 = (n_atom_array * max_n_nb15off + REORDER_THREADS - 1) / REORDER_THREADS;
    kernel_set_nb15off<<<blocks2, REORDER_THREADS>>>(d_atomids, d_atomids_rev, d_nb15off_orig, d_nb15off);
    return 0;
}

extern "C" int cuda_set_crd() {
    int blocks = (n_atom_array + REORDER_THREADS - 1) / REORDER_THREADS;
    kernel_set_crd<<<blocks, REORDER_THREADS>>>(d_atomids, d_crd, d_crd_chg, d_cell_z);

    return 0;
}

__global__ void kernel_set_work_orig(real_fc *d_work, const int *d_atomids) {
    int atomid = threadIdx.x + blockIdx.x * blockDim.x;
    if (atomid < D_N_ATOM_ARRAY) { // && d_atomids[atomid] >= 0){
        int index_orig = atomid * 3;
        for (int n = 1; n < N_MULTI_WORK; n++) {
            int index = (atomid + D_N_ATOM_ARRAY * n) * 3;
            d_work[index_orig + 0] += d_work[index + 0];
            d_work[index_orig + 1] += d_work[index + 1];
            d_work[index_orig + 2] += d_work[index + 2];
        }
    }
}

__global__ void kernel_reduction_energy(real_fc *d_energy) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    for (int i = 1; i < N_MULTI_WORK; i++) {
        if (tid == 1) {
            d_energy[0] += d_energy[i * 2];
            d_energy[1] += d_energy[i * 2 + 1];
            // printf("ene %f %f\n",d_energy[0],d_energy[1]);
        }
    }
}

__device__ bool check_15off64(const int  atom_idx1,
                              const int  atom_idx2,
                              const int *bitmask,
                              int &      mask_id,
                              int &      interact_bit) {
    int bit_pos  = atom_idx2 * N_ATOM_CELL + atom_idx1;
    mask_id      = bit_pos / 32;
    interact_bit = 1 << (bit_pos % 32);
    return (bitmask[mask_id] & interact_bit) == interact_bit;
}

__device__ real_pw check_15off(const int atomid1, const int atomid2, const int tmask_a1, const int tmask_a2) {
    int aid_diff = atomid2 - atomid1;
    int target   = tmask_a1;
    if (aid_diff < 0) {
        aid_diff = -aid_diff;
        target   = tmask_a2;
    }

    int mask                                             = 0;
    if (aid_diff <= 32) mask                             = 1 << (aid_diff - 1);
    real_pw valid_pair                                   = 1.0;
    if (mask != 0 && (mask & target) == mask) valid_pair = 0.0;
    return valid_pair;
}

__device__ real_pw cal_pair(real_pw &    w1,
                            real_pw &    w2,
                            real_pw &    w3,
                            real_pw &    ene_vdw,
                            real_pw &    ene_ele,
                            const real4 &crd_chg1,
                            const real4 &crd_chg2,
                            const int &  atomtype1,
                            const int &  atomtype2,
                            const real_pw *__restrict__ d_lj_6term,
                            const real_pw *__restrict__ d_lj_12term) {

    const real_pw d12[3]     = {crd_chg1.x - crd_chg2.x, crd_chg1.y - crd_chg2.y, crd_chg1.z - crd_chg2.z};
    const real_pw r12_2      = d12[0] * d12[0] + d12[1] * d12[1] + d12[2] * d12[2];
    const real_pw r12        = sqrt(r12_2);

    if (r12 >= D_CUTOFF) { return r12; }

    const real_pw r12_inv    = 1.0 / r12;
    const real_pw r12_2_inv  = r12_inv * r12_inv;
    const real_pw r12_3_inv  = r12_inv * r12_2_inv;
    const real_pw r12_6_inv  = r12_3_inv * r12_3_inv;
    const real_pw r12_12_inv = r12_6_inv * r12_6_inv;
    const real_pw term6      = d_lj_6term[atomtype1 * D_N_ATOMTYPES + atomtype2] * r12_6_inv;
    const real_pw term12     = d_lj_12term[atomtype1 * D_N_ATOMTYPES + atomtype2] * r12_12_inv;
    real_pw       work_coef  = r12_2_inv * (-12.0 * term12 + 6.0 * term6);
    const real_pw cc         = crd_chg1.w * crd_chg2.w * D_CHARGE_COEFF;
    work_coef -= cc * (r12_3_inv - D_FCOEFF);

    w1 = (work_coef)*d12[0];
    w2 = (work_coef)*d12[1];
    w3 = (work_coef)*d12[2];

    ene_ele = cc * (r12_inv - D_ZCORE + D_BCOEFF * r12_2);
    ene_vdw = (-term6 + term12);

    return r12;
}

__global__ void kernel_pairwise_ljzd(const real4 *d_crd_chg,
                                     CellPair *   d_cell_pairs,
                                     const int *  d_idx_head_cell_pairs,
                                     const int *  d_atomtype,
                                     const real_pw *__restrict__ d_lj_6term,
                                     const real_pw *__restrict__ d_lj_12term,
                                     real_fc *d_energy,
                                     real_fc *d_work) {
    // const bool flg_mod_15mask){
    real_fc ene_vdw = 0.0;
    real_fc ene_ele = 0.0;

    const int global_threadIdx = blockDim.x * blockIdx.x + threadIdx.x;
    const int c1               = global_threadIdx >> 5;
    const int warpIdx          = threadIdx.x >> 5;
    if (c1 >= D_N_CELLS) { return; }
    const int laneIdx = global_threadIdx & 31;
    const int n_loops = (d_idx_head_cell_pairs[c1 + 1] - d_idx_head_cell_pairs[c1]) * 2;

    const int ene_index_offset = global_threadIdx % N_MULTI_WORK;

    real_fc work_c1[3] = {0.0, 0.0, 0.0};

    const int atom_idx1 = (laneIdx & 7); // laneIdx%8
    const int a1        = c1 * N_ATOM_CELL + atom_idx1;

    __shared__ real4 crd_chg1[N_ATOM_CELL * (PW_THREADS >> 5)];
    __shared__ int   atomtype1[N_ATOM_CELL * (PW_THREADS >> 5)];

    const int sharedmem_idx = N_ATOM_CELL * warpIdx + atom_idx1;
    if (laneIdx < N_ATOM_CELL) {
        crd_chg1[sharedmem_idx]  = d_crd_chg[c1 * N_ATOM_CELL + laneIdx];
        atomtype1[sharedmem_idx] = d_atomtype[c1 * N_ATOM_CELL + laneIdx];
    }
    __syncthreads();

    CellPair cellpair;
    int      cp;

    for (int loopIdx = 0; loopIdx < n_loops; loopIdx++) {
        if (loopIdx % 2 == 0) {
            if (laneIdx == 0) {
                cp = d_idx_head_cell_pairs[c1] + (loopIdx >> 1);
                if (cp >= D_N_CELL_PAIRS) break;
                cellpair = d_cell_pairs[cp];
            }
            cp                    = __shfl(cp, 0);
            cellpair.cell_id1     = __shfl(cellpair.cell_id1, 0);
            cellpair.cell_id2     = __shfl(cellpair.cell_id2, 0);
            cellpair.image        = __shfl(cellpair.image, 0);
            cellpair.pair_mask[0] = __shfl(cellpair.pair_mask[0], 0);
            cellpair.pair_mask[1] = __shfl(cellpair.pair_mask[1], 0);
        }
        if (cellpair.cell_id1 != c1) break;
        const int c2 = cellpair.cell_id2;

        // atom_idx ... index in cell, 0-7
        const int atom_idx2 = (laneIdx >> 3) + 4 * (loopIdx % 2); // laneIdx/8 + 4*(warpIdx%2)

        // remove 1-2, 1-3, 1-4 pairs
        const int a2 = c2 * N_ATOM_CELL + atom_idx2;
        real4     crd_chg2;
        int       atomtype2;
        if (atom_idx1 == 0) {
            crd_chg2  = d_crd_chg[a2];
            atomtype2 = d_atomtype[a2];
            if ((cellpair.image & 1) == 1)
                crd_chg2.x -= PBC_L[0];
            else if ((cellpair.image & 2) == 2)
                crd_chg2.x += PBC_L[0];
            if ((cellpair.image & 4) == 4)
                crd_chg2.y -= PBC_L[1];
            else if ((cellpair.image & 8) == 8)
                crd_chg2.y += PBC_L[1];
            if ((cellpair.image & 16) == 16)
                crd_chg2.z -= PBC_L[2];
            else if ((cellpair.image & 32) == 32)
                crd_chg2.z += PBC_L[2];
        }
        int atomid2_top = laneIdx - laneIdx % 8;
        crd_chg2.x      = __shfl(crd_chg2.x, laneIdx - atom_idx1);
        crd_chg2.y      = __shfl(crd_chg2.y, laneIdx - atom_idx1);
        crd_chg2.z      = __shfl(crd_chg2.z, laneIdx - atom_idx1);
        crd_chg2.w      = __shfl(crd_chg2.w, laneIdx - atom_idx1);
        atomtype2       = __shfl(atomtype2, laneIdx - atom_idx1);

        real_pw w1 = 0.0, w2 = 0.0, w3 = 0.0;
        real_pw cur_ene_ele = 0.0;
        real_pw cur_ene_vdw = 0.0;
        int     mask_id;
        int     interact_bit;

        if (!check_15off64(atom_idx1, atom_idx2, cellpair.pair_mask, mask_id, interact_bit)) {

            real_pw r12 = cal_pair(w1, w2, w3, cur_ene_vdw, cur_ene_ele,
                                   // d_crd_chg[a1],
                                   crd_chg1[sharedmem_idx], crd_chg2,
                                   // d_atomtype[a1],
                                   atomtype1[sharedmem_idx], atomtype2, d_lj_6term, d_lj_12term);

            // if(flg_mod_15mask && r12 < D_CUTOFF_PAIRLIST) interact_bit = 0;

            ene_vdw += cur_ene_vdw;
            ene_ele += cur_ene_ele;

            work_c1[0] += w1;
            work_c1[1] += w2;
            work_c1[2] += w3;
        }
        /*if(flg_mod_15mask){
          for(int i = 32; i >= 1; i/=2){
          interact_bit |= __shfl_xor(interact_bit, i);
          }
          if(laneIdx == 0)
          d_cell_pairs[cp].pair_mask[mask_id]  |= interact_bit;
          }*/
        for (int i = 4; i >= 1; i /= 2) {
            w1 += shfl_xor(w1, i, 8);
            w2 += shfl_xor(w2, i, 8);
            w3 += shfl_xor(w3, i, 8);
        }

        if (laneIdx % 8 == 0) { // && (w1 != 0.0 || w2 != 0.0 || w3 != 0.0)){
            const int tmp_index = (((global_threadIdx / WARPSIZE) % N_MULTI_WORK) * D_N_ATOM_ARRAY + a2) * 3;
            atomicAdd2(&(d_work[tmp_index + 0]), -w1);
            atomicAdd2(&(d_work[tmp_index + 1]), -w2);
            atomicAdd2(&(d_work[tmp_index + 2]), -w3);
        }
    }
    for (int i = 16; i >= 8; i /= 2) {
        work_c1[0] += shfl_xor(work_c1[0], i, 32);
        work_c1[1] += shfl_xor(work_c1[1], i, 32);
        work_c1[2] += shfl_xor(work_c1[2], i, 32);
    }
    if (laneIdx < 8) {
        const int tmp_index = ((ene_index_offset * D_N_ATOM_ARRAY) + a1) * 3;
        atomicAdd2(&(d_work[tmp_index + 0]), work_c1[0]);
        atomicAdd2(&(d_work[tmp_index + 1]), work_c1[1]);
        atomicAdd2(&(d_work[tmp_index + 2]), work_c1[2]);
    }
    for (int i = 16; i >= 1; i /= 2) {
        ene_vdw += shfl_xor(ene_vdw, i, 32);
        ene_ele += shfl_xor(ene_ele, i, 32);
    }
    if (laneIdx == 0) {
        const int tmp_index = ((global_threadIdx / 32) % N_MULTI_WORK) * 2;
        atomicAdd2(&(d_energy[tmp_index + 0]), ene_vdw);
        atomicAdd2(&(d_energy[tmp_index + 1]), ene_ele);
    }
}

__global__ void set_idx_head_cell_pairs(const int *d_n_cell_pairs, int *d_idx_head_cell_pairs) {

    int n_rep = (D_N_CELLS + REORDER_THREADS - 1) / REORDER_THREADS;
    for (int i = 0; i < n_rep; i++) {
        const int idx_head  = REORDER_THREADS * i;
        const int idx_write = idx_head + threadIdx.x;
        if (idx_write < D_N_CELLS) {
            if (idx_write > 0) {
                const int idx = ((d_n_cell_pairs[idx_write - 1] + CP_PER_THREAD - 1) / CP_PER_THREAD) * CP_PER_THREAD;
                d_idx_head_cell_pairs[idx_write] = idx;
            } else {
                d_idx_head_cell_pairs[idx_write] = 0;
            }
        }
        for (int j = 1; j < REORDER_THREADS; j *= 2) {
            const int idx = (threadIdx.x / j);
            if (idx_write < D_N_CELLS && idx % 2 == 1) {
                d_idx_head_cell_pairs[idx_write] += d_idx_head_cell_pairs[idx_head + idx * j - 1];
            }
            __syncthreads();
        }
        if (i > 0) { d_idx_head_cell_pairs[idx_write] += d_idx_head_cell_pairs[idx_head - 1]; }

        __syncthreads();
    }

    if (threadIdx.x == 0) {
        const int idx = ((d_n_cell_pairs[D_N_CELLS - 1] + CP_PER_THREAD - 1) / CP_PER_THREAD) * CP_PER_THREAD;
        d_idx_head_cell_pairs[D_N_CELLS] = d_idx_head_cell_pairs[D_N_CELLS - 1] + idx;
        // d_idx_head_cell_pairs[D_N_CELLS] = d_idx_head_cell_pairs[D_N_CELLS-1] + d_n_cell_pairs[D_N_CELLS-1];
        // for(int i=0; i <= D_N_CELLS; i++){
        // printf("n:%d head:%d\n",
        // d_n_cell_pairs[i] , d_idx_head_cell_pairs[i]);
        //}
    }
    // printf("max cp: %d\n",idx_cp);
}

__global__ void pack_cellpairs_array(CellPair *      d_cell_pairs,
                                     const CellPair *d_cell_pairs_buf,
                                     const int *     d_n_cell_pairs,
                                     const int *     d_idx_head_cell_pairs) {

    const int cp = blockDim.x * blockIdx.x + threadIdx.x;

    if (cp >= D_MAX_N_CELL_PAIRS) return;
    const CellPair cellpair = d_cell_pairs_buf[cp];

    if (cellpair.cell_id1 < 0 || cellpair.cell_id2 < 0 || cellpair.cell_id1 >= D_N_CELLS
        || cellpair.cell_id2 >= D_N_CELLS) {
        return;
    }
    const int cp_in_cell1 = cp - cellpair.cell_id1 * D_MAX_N_CELL_PAIRS_PER_CELL;
    if (cp_in_cell1 >= d_n_cell_pairs[cellpair.cell_id1]) {
        printf("Error: cp_in_cell1:%d d_n_cell_pairs:%d cp:%d c1:%d head:%d\n", cp_in_cell1,
               d_n_cell_pairs[cellpair.cell_id1], cp, cellpair.cell_id1, d_idx_head_cell_pairs[cellpair.cell_id1]);
        return;
    }
    const int dest = d_idx_head_cell_pairs[cellpair.cell_id1] + cp_in_cell1;
    if (dest < cellpair.cell_id1) {
        printf("!!?? dest: %d (%d-%d) cp_in_cell:%d head:%d\n", dest, cellpair.cell_id1, cellpair.cell_id2, cp_in_cell1,
               d_idx_head_cell_pairs[cellpair.cell_id1]);
    }
    d_cell_pairs[dest] = cellpair;
}

__global__ void kernel_reset_cellpairs(CellPair *d_cell_pairs, int *d_n_cell_pairs, const int n_cells) {
    const int cell1_id = blockDim.x * blockIdx.x + threadIdx.x;
    if (cell1_id >= n_cells) { return; }
    int n_cp1 = d_n_cell_pairs[cell1_id];
    for (int cell2 = 0; cell2 < n_cp1; cell2++) {
        bool      flg        = true;
        int       n_mask_int = (N_ATOM_CELL * N_ATOM_CELL + 31) / 32;
        const int cp         = D_MAX_N_CELL_PAIRS_PER_CELL * cell1_id + cell2;
        for (int i = 0; i < n_mask_int; i++) flg &= (d_cell_pairs[cp].pair_mask[i] == ~0);
        if (flg) {
            d_n_cell_pairs[cell1_id]--;
            int cp_src       = D_MAX_N_CELL_PAIRS_PER_CELL * cell1_id + --n_cp1;
            d_cell_pairs[cp] = d_cell_pairs[cp_src];
        }
    }
    d_n_cell_pairs[cell1_id] = n_cp1;
}

extern "C" int cuda_pairwise_ljzd(const bool flg_mod_15mask) {
    HANDLE_ERROR(cudaMemset(d_energy, 0.0, sizeof(real_fc) * 2 * N_MULTI_WORK));
    HANDLE_ERROR(cudaMemset(d_work, 0.0, sizeof(real_fc) * max_n_atom_array * 3 * N_MULTI_WORK));

    cudaStreamCreate(&stream_pair_home);
    const int blocks = (n_cells + PW_THREADS / 32 - 1) / (PW_THREADS / 32);
    kernel_pairwise_ljzd<<<blocks, PW_THREADS, 0, stream_pair_home>>>(d_crd_chg, d_cell_pairs, d_idx_head_cell_pairs,
                                                                      d_atomtype, d_lj_6term, d_lj_12term, d_energy,
                                                                      d_work);
    // if(flg_mod_15mask){
    // const int blocks2 = (n_cal_cells+PW_THREADS-1) / PW_THREADS;
    // kernel_reset_cellpairs<<<blocks2, PW_THREADS, 0, stream_pair_home>>>
    //(d_cell_pairs, d_n_cell_pairs, n_cal_cells);
    //}
    return 0;
}
extern "C" int cuda_thread_sync() {
    cudaThreadSynchronize();
    return 0;
}
extern "C" int cuda_pair_sync() {
    cudaStreamSynchronize(stream_pair_home);
    cudaStreamDestroy(stream_pair_home);
    return 0;
}

extern "C" int cuda_memcpy_dtoh_work(real_fc *&h_work, real_fc *&h_energy, int n_atoms, int n_atom_array) {
    // printf("! cuda_memcpy_dtoh_work\n");
    int blocks = (n_atom_array + REORDER_THREADS - 1) / REORDER_THREADS;
    // printf("kernel_set_work_orig\n");
    if(N_MULTI_WORK > 1) {
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
    }
    HANDLE_ERROR(cudaMemcpy(h_work, d_work, sizeof(real_fc) * n_atom_array * 3, cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaMemcpy(h_energy, d_energy, sizeof(real_fc) * 2, cudaMemcpyDeviceToHost));
    //printf("cuda ene %f %f\n",h_energy[0], h_energy[1]);
    return 0;
}
int cuda_reset_work_ene() {
    HANDLE_ERROR(cudaMemset(d_work, 0.0, sizeof(real_fc) * max_n_atom_array * 3 * N_MULTI_WORK));
    HANDLE_ERROR(cudaMemset(d_energy, 0.0, sizeof(real_fc) * 2 * N_MULTI_WORK));
    return 0;
}
__device__ int get_column_id_from_crd(const int x, const int y) {
    return y * D_N_CELLS_XYZ[0] + x;
}

__device__ bool check_valid_pair(const int cell1_id, const int cell2_id) {
    const bool cell1_odd = cell1_id % 2 != 0;
    const bool cell2_odd = cell2_id % 2 != 0;
    if (cell1_odd) {
        if ((cell2_id < cell1_id && !cell2_odd) || (cell2_id > cell1_id && cell2_odd)) return false;
    } else {
        if ((cell2_id < cell1_id && cell2_odd) || (cell2_id > cell1_id && !cell2_odd)) return false;
    }
    return true;
}

__device__ int set_cell_pair_bitmask(const int  cell_id1,
                                     const int  cell_id2,
                                     const int  cell1_id_in_block,
                                     const int *s_nb15off,
                                     const int  n_atom_cell2,
                                     int *      pair_mask) {
    for (int i = 0; i < N_BITMASK; i++) pair_mask[i] = 0;
    int      a1                                      = N_ATOM_CELL * cell_id1;
    for (int a1_cell = 0; a1_cell < N_ATOM_CELL; a1++, a1_cell++) {
        bool flg1 = false;

        if (s_nb15off[(cell1_id_in_block * N_ATOM_CELL + a1_cell) * D_MAX_N_NB15OFF] == a1) flg1 = true;

        int a2 = N_ATOM_CELL * cell_id2;
        for (int a2_cell = 0; a2_cell < N_ATOM_CELL; a2++, a2_cell++) {
            const int bit_pos  = a2_cell * N_ATOM_CELL + a1_cell;
            const int mask_id  = bit_pos / 32;
            const int mask_pos = bit_pos % 32;
            const int add_bit  = 1 << mask_pos;
            bool      flg12    = false;
            if (flg1) flg12    = true;
            if (a2_cell >= n_atom_cell2)
                flg12 = true;

            else if ((cell_id1 == cell_id2 && a1 >= a2))
                flg12 = true;
            else {
                const int tail = (cell1_id_in_block * N_ATOM_CELL + a1_cell + 1) * D_MAX_N_NB15OFF;
                for (int i = tail - D_MAX_N_NB15OFF; i < tail && s_nb15off[i] != -1; i++) {
                    if (s_nb15off[i] == a2) {
                        flg12 = true;
                        break;
                    }
                }
            }
            if (flg12) { pair_mask[mask_id] |= add_bit; }
        }
    }
    return 0;
}

__device__ CellPair get_new_cell_pair(const int cell1_id, const int cell2_id, const int image[3]) {
    CellPair new_cp;
    new_cp.cell_id1 = cell1_id;
    new_cp.cell_id2 = cell2_id;
    int bit_image   = 0;
    if (image[0] == -1)
        bit_image = bit_image | 1;
    else if (image[0] == 1)
        bit_image = bit_image | 2;
    if (image[1] == -1)
        bit_image = bit_image | 4;
    else if (image[1] == 1)
        bit_image = bit_image | 8;
    if (image[2] == -1)
        bit_image = bit_image | 16;
    else if (image[2] == 1)
        bit_image       = bit_image | 32;
    new_cp.image        = bit_image;
    new_cp.pair_mask[0] = ~0;
    new_cp.pair_mask[1] = ~0;

    return new_cp;
}

__global__ void kernel_enumerate_cell_pair(const real4 *d_crd_chg,
                                           const real2 *d_cell_z,
                                           const int *  d_idx_xy_head_cell,
                                           const int *  d_idx_cell_column,
                                           const int *  d_atomids,
                                           const int *  d_nb15off_orig,
                                           const int *  d_nb15off,
                                           int *        d_n_cell_pairs,
                                           CellPair *   d_cell_pairs) {
    // 1 warp calculates pairs with a cell
    const int g_thread_id = (threadIdx.x + blockIdx.x * blockDim.x);
    const int cell1_id    = g_thread_id / D_N_NEIGHBOR_COL;
    if (cell1_id >= D_N_CELLS) return;
    const int neighbor_col_id = g_thread_id % D_N_NEIGHBOR_COL;

    int       cell1_crd[3];
    const int col1_id = d_idx_cell_column[cell1_id];
    cell1_crd[0]      = col1_id % D_N_CELLS_XYZ[0];
    cell1_crd[1]      = col1_id / D_N_CELLS_XYZ[0];
    cell1_crd[2]      = cell1_id - d_idx_xy_head_cell[col1_id];

    const real4   crd_chg11      = d_crd_chg[cell1_id * N_ATOM_CELL];
    const real_pw cell1_z_bottom = crd_chg11.z;
    // const real_pw cell1_z_top = d_cell_z[cell1_id].y;
    const real_pw cell1_z_top = d_crd_chg[(cell1_id + 1) * N_ATOM_CELL - 1].z;
    int           image[3]    = {0, 0, 0};

    const int idx_cell_pair_head = D_MAX_N_CELL_PAIRS_PER_CELL * cell1_id;
    int       d_cell[3];
    d_cell[0] = neighbor_col_id % (D_N_NEIGHBOR_XYZ[0] * 2 + 1) - D_N_NEIGHBOR_XYZ[0];
    d_cell[1] = neighbor_col_id / (D_N_NEIGHBOR_XYZ[0] * 2 + 1) - D_N_NEIGHBOR_XYZ[1];

    int     rel_x[3];
    int     cell2_crd[3];
    real_pw dx[3] = {0.0, 0.0, 0.0};
    for (int d = 0; d < 2; d++) {
        image[d]     = 0;
        rel_x[d]     = cell1_crd[d] + d_cell[d];
        cell2_crd[d] = rel_x[d];
        if (rel_x[d] < 0) {
            image[d]     = -1;
            cell2_crd[d] = D_N_CELLS_XYZ[d] + rel_x[d];
        } else if (rel_x[d] >= D_N_CELLS_XYZ[d]) {
            image[d]     = 1;
            cell2_crd[d] = rel_x[d] - D_N_CELLS_XYZ[d];
        }
        if (d_cell[d] > 0)
            dx[d] = (d_cell[d] - 1) * D_L_CELL_XYZ[d];
        else if (d_cell[d] < 0)
            dx[d] = (d_cell[d] + 1) * D_L_CELL_XYZ[d];
        dx[d]     = dx[d] * dx[d];
    }

    /*  20150612-- */
    const int col2_id         = cell2_crd[0] + cell2_crd[1] * D_N_CELLS_XYZ[0];
    const int n_cells_in_col2 = d_idx_xy_head_cell[col2_id + 1] - d_idx_xy_head_cell[col2_id];
    int       inc             = 1;
    bool      flg_up          = true;
    bool      flg_down        = true;

    for (d_cell[2] = 0; flg_up || flg_down; d_cell[2] += inc) {
        rel_x[2]     = cell1_crd[2] + d_cell[2];
        image[2]     = 0;
        cell2_crd[2] = rel_x[2];
        if (rel_x[2] < 0) {
            image[2] = -1;
            cell2_crd[2] += n_cells_in_col2;
        } else if (rel_x[2] >= n_cells_in_col2) {
            image[2] = 1;
            cell2_crd[2] -= n_cells_in_col2;
        }

        const int cell2_id = d_idx_xy_head_cell[col2_id] + cell2_crd[2];

        // const real_pw cell2_z_bottom = d_cell_z[cell2_id].x + image[2] * PBC_L[2];
        // const real_pw cell2_z_top = d_cell_z[cell2_id].y + image[2] * PBC_L[2];
        const real_pw cell2_z_bottom = d_crd_chg[cell2_id * N_ATOM_CELL].z + image[2] * PBC_L[2];
        const real_pw cell2_z_top    = d_crd_chg[(cell2_id + 1) * N_ATOM_CELL - 1].z + image[2] * PBC_L[2];

        if (cell2_z_top < cell1_z_bottom) {
            dx[2] = cell1_z_bottom - cell2_z_top;
            dx[2] = dx[2] * dx[2];
            if (inc == -1 && dx[0] + dx[1] + dx[2] > D_CUTOFF_PAIRLIST_2) { flg_down = false; }
        } else if (cell2_z_bottom > cell1_z_top) {
            dx[2] = cell2_z_bottom - cell1_z_top;
            dx[2] = dx[2] * dx[2];
            if (inc == 1 && dx[0] + dx[1] + dx[2] > D_CUTOFF_PAIRLIST_2) {
                d_cell[2] = 0;
                inc       = -1;
                flg_up    = false;
            }
        } else {
            dx[2] = 0.0;
        }

        if (dx[0] + dx[1] + dx[2] < D_CUTOFF_PAIRLIST_2) {

            if (check_valid_pair(cell1_id, cell2_id)) {
                const int cp_idx_cell = atomicAdd(&d_n_cell_pairs[cell1_id], 1);
                if (cp_idx_cell >= D_MAX_N_CELL_PAIRS_PER_CELL) {
                    printf("Index exceeds the maximum value. %d / %d\n", cp_idx_cell, D_MAX_N_CELL_PAIRS_PER_CELL);
                }

                d_cell_pairs[idx_cell_pair_head + cp_idx_cell] = get_new_cell_pair(cell1_id, cell2_id, image);
            }
        }
    }
}

__global__ void kernel_init_cell_pairs(CellPair *d_cell_pairs, CellPair *d_cell_pairs_buf) {
    const int cp_id = threadIdx.x + blockDim.x * blockIdx.x;
    if (cp_id >= D_MAX_N_CELL_PAIRS) return;
    CellPair new_cp;
    new_cp.cell_id1 = -1;
    new_cp.cell_id2 = -1;
    new_cp.image    = 0;
    for (int i = 0; i < N_BITMASK; i++) { new_cp.pair_mask[i] = ~0; }
    d_cell_pairs[cp_id]                                       = new_cp;
    d_cell_pairs_buf[cp_id]                                   = new_cp;
}
__global__ void set_cell_pairs_nb15off(const int *d_idx_head_cell_pairs,
                                       const int *d_nb15off,
                                       const int *d_atomids,
                                       const int  n_cell_pairs,
                                       CellPair * d_cell_pairs) {
    const int cp_block_first = blockIdx.x * blockDim.x;
    const int cp             = threadIdx.x + cp_block_first;
    if (cp >= n_cell_pairs) return;

    const int cell1_id = d_cell_pairs[cp].cell_id1;
    const int cell2_id = d_cell_pairs[cp].cell_id2;
    if (cell1_id < 0 || cell2_id < 0) return;
    const int cell1_id_in_block = cell1_id - d_cell_pairs[cp_block_first].cell_id1;
    if (cell1_id_in_block >= MAX_N_CELL_BLOCK) {
        printf("The number of cells in each block exceeds the constant MAX_N_CELL_BLOCK: %d / %d\ncell: %d %d - %d "
               "cp_block_first:%d "
               "c1_in_first:%d\n",
               cell1_id_in_block, MAX_N_CELL_BLOCK, cp, cell1_id, cell2_id, cp_block_first,
               d_cell_pairs[cp_block_first].cell_id1);
    }

    __shared__ int s_nb15off[MAX_N_CELL_BLOCK * N_ATOM_CELL * MAX_N_NB15OFF];

    if (threadIdx.x == 0 || d_cell_pairs[cp - 1].cell_id1 != cell1_id) {
        for (int i = 0; i < N_ATOM_CELL * MAX_N_NB15OFF; i++) {
            s_nb15off[cell1_id_in_block * N_ATOM_CELL * MAX_N_NB15OFF + i] =
                d_nb15off[cell1_id * N_ATOM_CELL * MAX_N_NB15OFF + i];
        }
    }
    __syncthreads();

    int       n_atom_cell2 = 0;
    const int tail         = (cell2_id + 1) * N_ATOM_CELL;
    for (int at2 = tail - N_ATOM_CELL; at2 < tail; at2++) {
        if (d_atomids[at2] >= 0)
            n_atom_cell2++;
        else
            break;
    }
    set_cell_pair_bitmask(cell1_id, cell2_id, cell1_id_in_block, s_nb15off, n_atom_cell2, d_cell_pairs[cp].pair_mask);
}

int set_idx_cell_column(const int *idx_atom_cell_xy, const int n_cells, const int *h_atomids) {
    for (int cell_id = 0; cell_id < n_cells; cell_id++) {
        h_idx_cell_column[cell_id] = idx_atom_cell_xy[cell_id * N_ATOM_CELL];
        // printf("cell: %d col: %d atid: %d \n",
        // cell_id, h_idx_cell_column[cell_id],
        // h_atomids[cell_id*N_ATOM_CELL]);
    }
    return 0;
}

extern "C" int cuda_enumerate_cell_pairs(int *&     h_atomids,
                                         const int  n_cells, // const int n_uni,
                                         const int  n_neighbor_col,
                                         const int *idx_atom_cell_xy) {

    set_idx_cell_column(idx_atom_cell_xy, n_cells, h_atomids);
    HANDLE_ERROR(cudaMemcpy(d_idx_cell_column, h_idx_cell_column, n_cells * sizeof(int), cudaMemcpyHostToDevice));
    cudaStream_t stream1;
    // cudaStream_t stream2;
    cudaStreamCreate(&stream1);
    // cudaStreamCreate(&stream2);

    // HANDLE_ERROR( cudaMemset(d_cell_pairs, -1, sizeof(int)*5*D_MAX_N_CELL_PAIRS));
    // HANDLE_ERROR( cudaMemset(d_cell_pairs_buf, -1, sizeof(int)*5*D_MAX_N_CELL_PAIRS));
    const int blocks3 = (max_n_cell_pairs + REORDER_THREADS - 1) / REORDER_THREADS;
    kernel_init_cell_pairs<<<blocks3, REORDER_THREADS, 0, stream1>>>(d_cell_pairs, d_cell_pairs_buf);

    const int blocks4 =
        (n_neighbor_col * n_cells + REORDER_THREADS - 1) / REORDER_THREADS; //  printf("bbb %d\n", max_n_cell_pairs);
    kernel_enumerate_cell_pair<<<blocks4, REORDER_THREADS, 0, stream1>>>(   // d_uni2cell_z,
        d_crd_chg, d_cell_z, d_idx_xy_head_cell, d_idx_cell_column, d_atomids, d_nb15off_orig, d_nb15off,
        d_n_cell_pairs, d_cell_pairs_buf);
    // d_cell_pairs);

    set_idx_head_cell_pairs<<<1, REORDER_THREADS, 0, stream1>>>(d_n_cell_pairs, d_idx_head_cell_pairs);
    // n_cell_pairs = d_idx_head_cell_pairs[n_cells+1];
    HANDLE_ERROR(cudaMemcpy(&n_cell_pairs, &d_idx_head_cell_pairs[n_cells], sizeof(int), cudaMemcpyDeviceToHost));
    HANDLE_ERROR(cudaMemcpyToSymbol(D_N_CELL_PAIRS, &n_cell_pairs, sizeof(int)));

    pack_cellpairs_array<<<blocks3, REORDER_THREADS, 0, stream1>>>(d_cell_pairs, d_cell_pairs_buf, d_n_cell_pairs,
                                                                   d_idx_head_cell_pairs);

    const int blocks5 = (n_cell_pairs + REORDER_THREADS - 1) / REORDER_THREADS;
    set_cell_pairs_nb15off<<<blocks5, REORDER_THREADS, 0, stream1>>>(d_idx_head_cell_pairs, d_nb15off, d_atomids,
                                                                     n_cell_pairs, d_cell_pairs);

    cudaStreamDestroy(stream1);

    return 0;
}

extern "C" int cuda_memcpy_htod_cell_pairs(CellPair *&h_cell_pairs, int *&h_idx_head_cell_pairs, int n_cell_pairs) {
    // printf("cuda_memcpy_htod_cell_pairs\n");
    HANDLE_ERROR(cudaMemcpy(d_cell_pairs, h_cell_pairs, n_cell_pairs * sizeof(CellPair), cudaMemcpyHostToDevice));
    // HANDLE_ERROR( cudaMemset(d_idx_head_cell_pairs, -1, sizeof(int)*(max_n_cells+1));
    HANDLE_ERROR(
        cudaMemcpy(d_idx_head_cell_pairs, h_idx_head_cell_pairs, (n_cells + 1) * sizeof(int), cudaMemcpyHostToDevice));
    HANDLE_ERROR(cudaMemcpyToSymbol(D_N_CELL_PAIRS, &n_cell_pairs, sizeof(int)));

    return 0;
}


extern "C" int cuda_alloc_set_hps_params(real_pw* h_hps_cutoff,
					 real_pw* h_hps_lambda,
					 int       n_lj_types){
    // printf("threads : %d\n", PW_THREADS);
  printf("cuda_alloc_set_hps_params\n");
  const unsigned int size_lj_matrix = sizeof(real_pw) * n_lj_types * n_lj_types;
  // cudaMalloc
  HANDLE_ERROR(cudaMalloc((void **)&d_hps_cutoff, size_lj_matrix));
  HANDLE_ERROR(cudaMalloc((void **)&d_hps_lambda, size_lj_matrix));
  
  HANDLE_ERROR(cudaMemcpy(d_hps_cutoff, h_hps_cutoff, size_lj_matrix, cudaMemcpyHostToDevice));
  HANDLE_ERROR(cudaMemcpy(d_hps_lambda, h_hps_lambda, size_lj_matrix, cudaMemcpyHostToDevice));
  //printf("cudadbg %d\n", d_hps_lambda[0]);
  return 0;
}
extern "C" int cuda_free_hps_params() {
  // printf("cuda_free_lj_param\n");
  // cudaUnbindTexture(tex_lj_6term);
  // cudaUnbindTexture(tex_lj_12term);
  HANDLE_ERROR(cudaFree(d_hps_cutoff));
  HANDLE_ERROR(cudaFree(d_hps_lambda));
  return 0;
}

extern "C" int cuda_hps_constant(real_pw hps_eps){
  printf("set_cuda_hps_constant\n");
  HANDLE_ERROR(cudaMemcpyToSymbol(D_HPS_EPS, &hps_eps, sizeof(real_pw)));      
  return 0;
}
extern "C" int cuda_debye_huckel_constant(real_pw in_dielect, real_pw in_temperature,
					  real_pw in_ionic_strength){
  printf("cuda_debye_huckel_constant (cuda)\n");
  debye_length_inv = (1.0/(sqrt(PERMITTIVITY*in_dielect*BOLTZMAN*in_temperature/(2*AVOGADRO*ELEM_CHARGE*ELEM_CHARGE*in_ionic_strength))*1e10));
  dielect_inv = 1.0 / in_dielect;
  HANDLE_ERROR(cudaMemcpyToSymbol(D_DEBYE_LEN_INV, &debye_length_inv, sizeof(real_pw)));      
  printf("cuda debye_len_inv : %f\n", debye_length_inv);
  HANDLE_ERROR(cudaMemcpyToSymbol(D_DIELECT_INV, &dielect_inv, sizeof(real_pw)));      
  printf("cuda dielect_inv : %f\n", dielect_inv);

  return 0;
}


__device__ real_pw cal_pair_hps_dh(real_pw &    w1,
				   real_pw &    w2,
				   real_pw &    w3,
				   real_pw &    ene_vdw,
				   real_pw &    ene_ele,
				   const real4 &crd_chg1,
				   const real4 &crd_chg2,
				   const int &  atomtype1,
				   const int &  atomtype2,
				   const real_pw *__restrict__ d_lj_6term,
				   const real_pw *__restrict__ d_lj_12term,
				   const real_pw *__restrict__ d_hps_cutoff,
				   const real_pw *__restrict__ d_hps_lambda){
    const real_pw d12[3]     = {crd_chg1.x - crd_chg2.x, crd_chg1.y - crd_chg2.y, crd_chg1.z - crd_chg2.z};
    const real_pw r12_2      = d12[0] * d12[0] + d12[1] * d12[1] + d12[2] * d12[2];
    const real_pw r12        = sqrt(r12_2);

    if (r12 >= D_CUTOFF) { return r12; }

    const real_pw r12_inv    = 1.0 / r12;
    const real_pw r12_2_inv  = r12_inv * r12_inv;
    const real_pw r12_3_inv  = r12_inv * r12_2_inv;
    const real_pw r12_6_inv  = r12_3_inv * r12_3_inv;
    const real_pw r12_12_inv = r12_6_inv * r12_6_inv;
    const int pairtype = atomtype1 * D_N_ATOMTYPES + atomtype2;
    const real_pw hps_cutoff = d_hps_cutoff[pairtype];
    const real_pw hps_lambda = d_hps_lambda[pairtype];
    const real_pw term6      = d_lj_6term[pairtype] * r12_6_inv;
    const real_pw term12     = d_lj_12term[pairtype] * r12_12_inv;
    real_pw       work_coef  = r12_2_inv * (-12.0 * term12 + 6.0 * term6);

    ene_vdw = (-term6 + term12);

    if(r12 >  hps_cutoff){
      ene_vdw *= hps_lambda;
      work_coef *= hps_lambda;
    }else{
      ene_vdw += (1-hps_lambda) * D_HPS_EPS;
    }

    const real_pw cc         = crd_chg1.w * crd_chg2.w * D_CHARGE_COEFF;
    const real_pw r12_ld_exp = exp(-r12 * D_DEBYE_LEN_INV);

    ene_ele = cc * r12_inv * D_DIELECT_INV * r12_ld_exp;
    work_coef -= ene_ele * (r12_inv + D_DEBYE_LEN_INV);

    w1 = (work_coef)*d12[0];
    w2 = (work_coef)*d12[1];
    w3 = (work_coef)*d12[2];

    return r12;
}

__global__ void kernel_pairwise_hps_dh(const real4 *d_crd_chg,
				       CellPair *   d_cell_pairs,
				       const int *  d_idx_head_cell_pairs,
				       const int *  d_atomtype,
				       const real_pw *__restrict__ d_lj_6term,
				       const real_pw *__restrict__ d_lj_12term,
				       const real_pw *__restrict__ d_hps_cutoff,
				       const real_pw *__restrict__ d_hps_lambda,
				       real_fc *d_energy,
				       real_fc *d_work) {
    // const bool flg_mod_15mask){
    real_fc ene_vdw = 0.0;
    real_fc ene_ele = 0.0;

    const int global_threadIdx = blockDim.x * blockIdx.x + threadIdx.x;
    const int c1               = global_threadIdx >> 5;
    const int warpIdx          = threadIdx.x >> 5;
    if (c1 >= D_N_CELLS) { return; }
    const int laneIdx = global_threadIdx & 31;
    const int n_loops = (d_idx_head_cell_pairs[c1 + 1] - d_idx_head_cell_pairs[c1]) * 2;

    const int ene_index_offset = global_threadIdx % N_MULTI_WORK;

    real_fc work_c1[3] = {0.0, 0.0, 0.0};

    const int atom_idx1 = (laneIdx & 7); // laneIdx%8
    const int a1        = c1 * N_ATOM_CELL + atom_idx1;

    __shared__ real4 crd_chg1[N_ATOM_CELL * (PW_THREADS >> 5)];
    __shared__ int   atomtype1[N_ATOM_CELL * (PW_THREADS >> 5)];

    const int sharedmem_idx = N_ATOM_CELL * warpIdx + atom_idx1;
    if (laneIdx < N_ATOM_CELL) {
        crd_chg1[sharedmem_idx]  = d_crd_chg[c1 * N_ATOM_CELL + laneIdx];
        atomtype1[sharedmem_idx] = d_atomtype[c1 * N_ATOM_CELL + laneIdx];
    }
    __syncthreads();

    CellPair cellpair;
    int      cp;

    for (int loopIdx = 0; loopIdx < n_loops; loopIdx++) {
        if (loopIdx % 2 == 0) {
            if (laneIdx == 0) {
                cp = d_idx_head_cell_pairs[c1] + (loopIdx >> 1);
                if (cp >= D_N_CELL_PAIRS) break;
                cellpair = d_cell_pairs[cp];
            }
            cp                    = __shfl(cp, 0);
            cellpair.cell_id1     = __shfl(cellpair.cell_id1, 0);
            cellpair.cell_id2     = __shfl(cellpair.cell_id2, 0);
            cellpair.image        = __shfl(cellpair.image, 0);
            cellpair.pair_mask[0] = __shfl(cellpair.pair_mask[0], 0);
            cellpair.pair_mask[1] = __shfl(cellpair.pair_mask[1], 0);
        }
        if (cellpair.cell_id1 != c1) break;
        const int c2 = cellpair.cell_id2;

        // atom_idx ... index in cell, 0-7
        const int atom_idx2 = (laneIdx >> 3) + 4 * (loopIdx % 2); // laneIdx/8 + 4*(warpIdx%2)

        // remove 1-2, 1-3, 1-4 pairs
        const int a2 = c2 * N_ATOM_CELL + atom_idx2;
        real4     crd_chg2;
        int       atomtype2;
        if (atom_idx1 == 0) {
            crd_chg2  = d_crd_chg[a2];
            atomtype2 = d_atomtype[a2];
            if ((cellpair.image & 1) == 1)
                crd_chg2.x -= PBC_L[0];
            else if ((cellpair.image & 2) == 2)
                crd_chg2.x += PBC_L[0];
            if ((cellpair.image & 4) == 4)
                crd_chg2.y -= PBC_L[1];
            else if ((cellpair.image & 8) == 8)
                crd_chg2.y += PBC_L[1];
            if ((cellpair.image & 16) == 16)
                crd_chg2.z -= PBC_L[2];
            else if ((cellpair.image & 32) == 32)
                crd_chg2.z += PBC_L[2];
        }
        int atomid2_top = laneIdx - laneIdx % 8;
        crd_chg2.x      = __shfl(crd_chg2.x, laneIdx - atom_idx1);
        crd_chg2.y      = __shfl(crd_chg2.y, laneIdx - atom_idx1);
        crd_chg2.z      = __shfl(crd_chg2.z, laneIdx - atom_idx1);
        crd_chg2.w      = __shfl(crd_chg2.w, laneIdx - atom_idx1);
        atomtype2       = __shfl(atomtype2, laneIdx - atom_idx1);

        real_pw w1 = 0.0, w2 = 0.0, w3 = 0.0;
        real_pw cur_ene_ele = 0.0;
        real_pw cur_ene_vdw = 0.0;
        int     mask_id;
        int     interact_bit;

        if (!check_15off64(atom_idx1, atom_idx2, cellpair.pair_mask, mask_id, interact_bit)) {

            real_pw r12 = cal_pair_hps_dh(w1, w2, w3, cur_ene_vdw, cur_ene_ele,
					  // d_crd_chg[a1],
					  crd_chg1[sharedmem_idx], crd_chg2,
					  // d_atomtype[a1],
					  atomtype1[sharedmem_idx], atomtype2, d_lj_6term, d_lj_12term,
					  d_hps_cutoff, d_hps_lambda);
            // if(flg_mod_15mask && r12 < D_CUTOFF_PAIRLIST) interact_bit = 0;
            ene_vdw += cur_ene_vdw;
            ene_ele += cur_ene_ele;

            work_c1[0] += w1;
            work_c1[1] += w2;
            work_c1[2] += w3;
        }
        /*if(flg_mod_15mask){
          for(int i = 32; i >= 1; i/=2){
          interact_bit |= __shfl_xor(interact_bit, i);
          }
          if(laneIdx == 0)
          d_cell_pairs[cp].pair_mask[mask_id]  |= interact_bit;
          }*/
        for (int i = 4; i >= 1; i /= 2) {
            w1 += shfl_xor(w1, i, 8);
            w2 += shfl_xor(w2, i, 8);
            w3 += shfl_xor(w3, i, 8);
        }

        if (laneIdx % 8 == 0) { // && (w1 != 0.0 || w2 != 0.0 || w3 != 0.0)){
            const int tmp_index = (((global_threadIdx / WARPSIZE) % N_MULTI_WORK) * D_N_ATOM_ARRAY + a2) * 3;
            atomicAdd2(&(d_work[tmp_index + 0]), -w1);
            atomicAdd2(&(d_work[tmp_index + 1]), -w2);
            atomicAdd2(&(d_work[tmp_index + 2]), -w3);
        }
    }
    for (int i = 16; i >= 8; i /= 2) {
        work_c1[0] += shfl_xor(work_c1[0], i, 32);
        work_c1[1] += shfl_xor(work_c1[1], i, 32);
        work_c1[2] += shfl_xor(work_c1[2], i, 32);
    }
    if (laneIdx < 8) {
        const int tmp_index = ((ene_index_offset * D_N_ATOM_ARRAY) + a1) * 3;
        atomicAdd2(&(d_work[tmp_index + 0]), work_c1[0]);
        atomicAdd2(&(d_work[tmp_index + 1]), work_c1[1]);
        atomicAdd2(&(d_work[tmp_index + 2]), work_c1[2]);
    }
    for (int i = 16; i >= 1; i /= 2) {
        ene_vdw += shfl_xor(ene_vdw, i, 32);
        ene_ele += shfl_xor(ene_ele, i, 32);
    }
    if (laneIdx == 0) {
        const int tmp_index = ((global_threadIdx / 32) % N_MULTI_WORK) * 2;
        atomicAdd2(&(d_energy[tmp_index + 0]), ene_vdw);
        atomicAdd2(&(d_energy[tmp_index + 1]), ene_ele);
    }
}
extern "C" int cuda_pairwise_hps_dh(const bool flg_mod_15mask) {
    HANDLE_ERROR(cudaMemset(d_energy, 0.0, sizeof(real_fc) * 2 * N_MULTI_WORK));
    HANDLE_ERROR(cudaMemset(d_work, 0.0, sizeof(real_fc) * max_n_atom_array * 3 * N_MULTI_WORK));
    cudaStreamCreate(&stream_pair_home);
    const int blocks = (n_cells + PW_THREADS / 32 - 1) / (PW_THREADS / 32);
    kernel_pairwise_hps_dh<<<blocks, PW_THREADS, 0, stream_pair_home>>>(d_crd_chg, d_cell_pairs, d_idx_head_cell_pairs,
									d_atomtype, d_lj_6term, d_lj_12term,
									d_hps_cutoff, d_hps_lambda,
									d_energy,
									d_work);
    

    return 0;
}

