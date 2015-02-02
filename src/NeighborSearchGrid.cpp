#include "NeighborSearchGrid.h"

#ifdef F_CUDA
extern "C" int cuda_hostalloc_atom_info(real_pw*& h_crd, int*& h_atomids,
					real_fc*& h_work, real_pw*& h_energy,
					int n_atoms);
extern "C" int cuda_hostfree_atom_info(real_pw* h_crd, int* h_atomids);
#endif

NeighborSearchGrid::NeighborSearchGrid(){
}

NeighborSearchGrid::~NeighborSearchGrid(){
  free_variables();
}

int NeighborSearchGrid::set_grid(int in_n_atoms, 
				 real in_cutoff_pair,
				 real in_min_grid_width,
				 int in_max_n_atoms,
				 const PBC* in_pbc){
  n_atoms = in_n_atoms;
  cutoff_pair = in_cutoff_pair;
  pbc = in_pbc;
  min_grid_width = in_min_grid_width;
  cout << "NeighborSearchGrid::set_grid" << endl;
  for (int d=0; d < 3; d++){
    n_grid[d] = floor(pbc->L[d] / min_grid_width);
    L_grid[d] = pbc->L[d] / (real)n_grid[d];
  }
  n_all_grids = n_grid[0] * n_grid[1] * n_grid[2];
  
  cout << "Grid ... " << n_grid[0] << " * " << n_grid[1];
  cout << " * " << n_grid[2] << " = " << n_all_grids << endl;
  cout << "Grid size ... " << L_grid[0] << ", ";
  cout << L_grid[1] << ", " << L_grid[2] << endl;
  cout << "Cutoff for grid search ... " << cutoff_pair << endl;

  max_n_atoms = ceil(ceil((float)n_atoms / (float)n_all_grids) * 1.2);
  if (max_n_atoms < in_max_n_atoms)  max_n_atoms = in_max_n_atoms;
  int tmp_max_n_atoms = ceil(L_grid[0]*L_grid[1]*L_grid[2] * 0.1);
  if (max_n_atoms < tmp_max_n_atoms)  max_n_atoms = tmp_max_n_atoms;
  
  cout << "Maximum number of atoms for each grid = " << max_n_atoms << endl;

  for(int d=0; d < 3; d++){
    n_neighbors_axis[d] = ceil(cutoff_pair / L_grid[d]);
  }
  max_n_grid_pairs = ((n_neighbors_axis[0]*2+1) *
		      (n_neighbors_axis[1]*2+1) * (n_neighbors_axis[2]*2+1) /2) * n_all_grids;


  n_grid_pairs = 0;
  return 0;
}

int NeighborSearchGrid::alloc_variables(){

  crd_grid = new real*[n_atoms];
  for(int i=0; i < n_atoms; i++){
    crd_grid[i] = new real[3];
    for(int j=0; j < 3; j++) crd_grid[i][j] = 0.0;
  }

  atom_grid = new int*[n_atoms];
  for(int i=0; i < n_atoms; i++){
    atom_grid[i] = new int[3];
  }

  grid_atom = new int*[n_all_grids];
  for(int i=0; i < n_all_grids; i++){
    grid_atom[i] = new int[max_n_atoms];
  }

  grid_atom_index = new int[n_all_grids + 1];
  grid_n_atoms = new int[n_all_grids];
  grid_pairs = new GridPair[max_n_grid_pairs];

#ifdef F_CUDA
  cout << "cuda_hostalloc_atom_info"<<endl;
  cuda_hostalloc_atom_info(crd, atomids, work, energy, n_atoms);
#else
  crd = new real_pw[n_atoms*3];
  atomids = new int[n_atoms];
  work = new real_fc[n_atoms*3];
  energy = new real_pw[2];
#endif
  return 0;
}

int NeighborSearchGrid::free_variables(){
  for(int i=0; i < n_atoms; i++){
    delete[] crd_grid[i];
  }
  delete[] crd_grid;

  for(int i=0; i < n_atoms; i++){
    delete[] atom_grid[i];
  }
  delete[] atom_grid;
  
  for(int i=0; i < n_all_grids; i++){
    delete[] grid_atom[i];
  }
  delete[] grid_atom;
  
  delete[] grid_atom_index;

  delete[] grid_n_atoms;
  
  delete[] grid_pairs;

#ifdef F_CUDA
  cuda_hostfree_atom_info(crd, atomids);
#else
  delete[] crd;
  delete[] atomids;
  delete[] work;
  delete[] energy;
#endif

  return 0;
}

int NeighborSearchGrid::set_crd_grid(real** in_crd, real_pw* in_charge,
				     int* in_atom_type){
  // set atominfo, grid_n_atoms, grid_atom
  //cout << "set_crd_grid" << endl;
  for(int gid=0; gid < n_all_grids; gid++){
    grid_n_atoms[gid] = 0;
  }
  for (int atom_id=0; atom_id < n_atoms; atom_id++){
    int gx = atom_grid[atom_id][0];
    int gy = atom_grid[atom_id][1];
    int gz = atom_grid[atom_id][2];
    int grid_id = get_grid_id(gx, gy, gz);
    // int i = grid_n_atoms[grid_id];
    //  cout << "atom_id:" << atom_id << endl;
    int i = (grid_atom_index[grid_id] + grid_n_atoms[grid_id]);
    crd[i*3]   = (real_pw)in_crd[atom_id][0];
    crd[i*3+1] = (real_pw)in_crd[atom_id][1];
    crd[i*3+2] = (real_pw)in_crd[atom_id][2];
    atomids[i] = atom_id;
    //atominfo[i].charge = (real_pw)in_charge[atom_id];
    //atominfo[i].atom_type = in_atom_type[atom_id];
    //cout << "grid_atom" <<endl;
    grid_atom[grid_id][grid_n_atoms[grid_id]] = atom_id;
    //cout << "grid_n_atoms" << endl;
    grid_n_atoms[grid_id] += 1;
  }  
  return 0;
}

int NeighborSearchGrid::setup_crd_into_grid(real** in_crd, real_pw* in_charge,
					    int* in_atom_type){
  for (int atom_id=0; atom_id < n_atoms; atom_id++){
    int grid_crd[3];
    //cout << "set_crd_into_grid atomid="<<atom_id<<endl;
    for (int d=0; d < 3; d++){
      real x = in_crd[atom_id][d] - pbc->lower_bound[d];
      while(x < 0.0) x += pbc->L[d];
      while(x >= pbc->L[d]) x -= pbc->L[d];
      //cout << "dim : " << d << endl;
      //cout << "in_crd[atom_id][d]: " << in_crd[atom_id][d] << endl;
      //cout << "in_crd[atom_id][d]-origin: " << in_crd[atom_id][d] - pbc->origin[d] << endl;
      grid_crd[d] = floor(x / L_grid[d]);
      //cout << "grid_crd[d]: " << grid_crd[d] << endl;

      //cout << "in_crd[atom_id][d]: " << in_crd[atom_id][d] << endl;
      //cout << "crd_grid[atom_id][d]: " << crd_grid[atom_id][d] << endl;
      //cout << "tmp1 " << d <<endl;
      atom_grid[atom_id][d] = grid_crd[d];
      crd_grid[atom_id][d] = (x - L_grid[d] * atom_grid[atom_id][d]) / L_grid[d];
      if(crd_grid[atom_id][d] < 0 || crd_grid[atom_id][d] >= L_grid[d]){
	cout << "!! Grid error: atom:" << atom_id << " : " << crd_grid[atom_id][0] <<endl;
      }
    }

    //cout << " : " << crd_grid[atom_id][0] << "," 
    //<< crd_grid[atom_id][1] <<", "
    //<< crd_grid[atom_id][2] << endl;
    int gx = grid_crd[0];
    int gy = grid_crd[1];
    int gz = grid_crd[2];
    int grid_id = get_grid_id(gx, gy, gz);
    grid_n_atoms[grid_id] += 1;
  }
  set_grid_atom_index_from_grid_n_atoms();
  set_crd_grid(in_crd, in_charge, in_atom_type);
  return 0;
}
int NeighborSearchGrid::print_max_n_atoms_in_grid(){
  int max_grid_id = -1;
  int max_num = 0;
  for(int grid_id=0; grid_id < n_all_grids; grid_id++){
    if(grid_n_atoms[grid_id] >= max_num){
      max_num = grid_n_atoms[grid_id];
      max_grid_id = grid_id;
    }
  }
  cout << "Max num of atoms in grid: " << max_num << " atoms in grid " << max_grid_id << endl;
  return max_num;
}
int NeighborSearchGrid::set_grid_atom_index_from_grid_n_atoms(){
  int atom_index = 0;
  for(int grid_id=0; grid_id < n_all_grids; grid_id++){
    grid_atom_index[grid_id] = atom_index;
    atom_index += grid_n_atoms[grid_id];
  }
  grid_atom_index[n_all_grids] = atom_index;
  return 0;
}
int NeighborSearchGrid::reflesh_atom_assign(real** in_crd, real_pw* in_charge, int* in_atom_type){
  // check_n_atoms_in_grid();
  // crd_grid
  //cout << "reflesh_atom_assign" << endl;
  for (int atom_id=0; atom_id < n_atoms; atom_id++){
    //cout << "  atom_id : " << atom_id << endl;
    // move ... -1, 0, or 1
    int mg[3] = {0,0,0};
    for (int d=0; d < 3; d++){
      if (crd_grid[atom_id][d] < 0.0){
	//cout << d << " : " << crd_grid[atom_id][d] << endl;
	crd_grid[atom_id][d] += 1.0;
	mg[d] -= 1;
      } else if(crd_grid[atom_id][d] > 1.0){
	//cout << d << " : " << crd_grid[atom_id][d] << endl;
	crd_grid[atom_id][d] -= 1.0;
	mg[d] += 1;
      }
    }

    // atom_grid
    //cout << "  update atom_grid : " << endl;
    int fg[3]; // from
    int tg[3]; // to 
    for(int d=0; d < 3; d++){
      fg[d] = atom_grid[atom_id][d];
      tg[d] = fg[d] + mg[d];
      if(tg[d] == -1) tg[d] = n_grid[d] - 1;
      else if(tg[d] == n_grid[d]) tg[d] = 0;
      atom_grid[atom_id][d] = tg[d];
    }
    
    // remove from the previous grid
    //cout << "  remove from the previous grid : " << endl;
    int grid_id_from = get_grid_id(fg[0], fg[1], fg[2]);
    int grid_id_to = get_grid_id(tg[0], tg[1], tg[2]);
    grid_n_atoms[grid_id_from] -= 1; 
    grid_n_atoms[grid_id_to] += 1; 

    //int from_g_n = grid_n_atoms[grid_id];
    //for (int i=0; i < from_g_n; i++){
    //if(grid_atom[fg[0]][fg[1]][fg[2]][i] == atom_id){
    //grid_atom[fg[0]][fg[1]][fg[2]][i] = 
    //grid_atom[fg[0]][fg[1]][fg[2]][from_g_n-1];
    //grid_atom[fg[0]][fg[1]][fg[2]][from_g_n-1] = -1;
    //}
    //}
    // grid_n_atoms[fg[0]][fg[1]][fg[2]] -= 1;
    
    // add to the next grid
    //cout << "  add to the grid : " <<tg[0]<<","<<tg[1]<<","<<tg[2]<< endl;
    // int to_g_n = grid_n_atoms[tg[0]][tg[1]][tg[2]];
    //cout << "    to_g_n: " << to_g_n << endl; 
    // grid_atom[tg[0]][tg[1]][tg[2]][to_g_n] = atom_id;
    //cout << "    grid_atom: " << grid_atom[tg[0]][tg[1]][tg[2]][to_g_n] << endl; 
    // grid_n_atoms[tg[0]][tg[1]][tg[2]] += 1;
  }
  set_grid_atom_index_from_grid_n_atoms();
  set_crd_grid(in_crd, in_charge, in_atom_type);
  //check_n_atoms_in_grid();
  return 0;
}

int NeighborSearchGrid::enumerate_grid_pairs(){
  // Enumurating pairs of grids
  // Each grid sees only left, bottom, and front grids
  int n_pairs = 0;
  int grid1[3];
  for(grid1[0]=0; grid1[0] < n_grid[0]; grid1[0]++){
    for(grid1[1]=0; grid1[1] < n_grid[1]; grid1[1]++){
      for(grid1[2]=0; grid1[2] < n_grid[2]; grid1[2]++){
	int tmp_image[3] = {0, 0, 0};
	//add_grid_pairs(grid1, grid1, tmp_image);
	real pos1[8][3];
	get_grid_positions(pos1, grid1[0], grid1[1], grid1[2]);

	int g_id1 = get_grid_id(grid1[0], grid1[1], grid1[2]);
	int d_grid[3];
	int grid2rel[3];

	for(d_grid[0] = -n_neighbors_axis[0];
	    d_grid[0] <= n_neighbors_axis[0]; d_grid[0]++){
	  grid2rel[0] = grid1[0] + d_grid[0];

	  //cout << "d_grid[0]:"<<d_grid[0] <<endl;
	  for(d_grid[1] = -n_neighbors_axis[1];
	      d_grid[1] <= n_neighbors_axis[1]; d_grid[1]++){
	    //cout << "d_grid[1]:"<<d_grid[1] <<endl;
	    grid2rel[1] = grid1[1] + d_grid[1];
	    
	    for(d_grid[2] = -n_neighbors_axis[2];
		d_grid[2] <= n_neighbors_axis[2]; d_grid[2]++){
	      grid2rel[2] = grid1[2] + d_grid[2];
	      int image[3];
	      int grid2[3];
	      for (int d=0; d < 3; d++){
		if (grid2rel[d] < 0){
		  image[d] = -1; grid2[d] = n_grid[d] + grid2rel[d];
		}else if (grid2rel[d] >= n_grid[d]){
		  image[d] = 1;  grid2[d] = grid2rel[d] - n_grid[d];
		}else{
		  image[d] = 0;  grid2[d] = grid2rel[d];
		}
	      }

	      int g_id2 = get_grid_id(grid2[0], grid2[1], grid2[2]);
	      if (g_id1 > g_id2) continue;

	      real pos2[8][3];
	      get_grid_positions(pos2, grid2rel[0], grid2rel[1], grid2rel[2]);
	      bool flg_add = false;
	      bool flg_all = true;
	      for(int gp1=0; gp1 < 8 && !flg_add; gp1++){
		for(int gp2=0; gp2 < 8 && !flg_add; gp2++){
		  real diff_x = pos1[gp1][0]-pos2[gp2][0];
		  real diff_y = pos1[gp1][1]-pos2[gp2][1];
		  real diff_z = pos1[gp1][2]-pos2[gp2][2];
		  real dist = sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z);
		  if(dist <= cutoff_pair)
		    flg_add = true;
		  else
		    flg_all = false;
		}
	      }
	      
	      //if (flg_add and !check_grid_pairs_exist(g_id1, g_id2)){
	      if(flg_add)
		add_grid_pairs(g_id1, g_id2, image);
	    }
	  }
	}
      }
    }
  }
  cout << "Number of grid pairs ... " << n_grid_pairs << endl;
  return n_pairs;
}
int NeighborSearchGrid::get_grid_positions(real positions[8][3],
					   int gx, int gy, int gz){
  int i = 0;
  
  for (int dx = 0; dx < 2; dx++){
    int ngx = gx+dx;
    //if(ngx < 0)  ngx = n_grid[0] + ngx + 1;
    for (int dy = 0; dy < 2; dy++){
      int ngy = gy+dy;
      //if(ngy < 0)  ngy = n_grid[1] + ngy + 1;
      for (int dz = 0; dz < 2; dz++){
	int ngz = gz+dz;
	//if(ngz < 0)  ngz = n_grid[2] + ngz + 1;
	positions[i][0] = ngx * L_grid[0];
	positions[i][1] = ngy * L_grid[1];
	positions[i][2] = ngz * L_grid[2];
	i++;
      }
    }
  }
  return 0;
}

bool NeighborSearchGrid::check_grid_pairs_exist(int grid_id1, int grid_id2){
  for(int i=0; i < n_grid_pairs; i++){
    bool flg1 = true;
    bool flg2 = true;
    if(grid_pairs[i].grid_id1 == grid_id1 && 
       grid_pairs[i].grid_id2 == grid_id2)
      return true;
  }
  return false;
}

int NeighborSearchGrid::add_grid_pairs(int grid_id1, int grid_id2, int image[3]){
  if (grid_id1 == grid_id2) return 1;
  grid_pairs[n_grid_pairs].grid_id1 = grid_id1;
  grid_pairs[n_grid_pairs].grid_id2 = grid_id2;
  grid_pairs[n_grid_pairs].cx = image[0];
  grid_pairs[n_grid_pairs].cy = image[1];
  grid_pairs[n_grid_pairs].cz = image[2];

  int natoms1 = grid_n_atoms[grid_id1];
  int natoms2 = grid_n_atoms[grid_id2];
  grid_pairs[n_grid_pairs].n_pairs = natoms1* natoms2;
  //cout <<  "add_grid_pairs " << n_grid_pairs;
  //for(int d=0; d < 9; d++){
  //cout << " " <<     grid_pairs[n_grid_pairs][d];
  //}
  //cout << endl;
  n_grid_pairs++;
  return 0;
}

int NeighborSearchGrid::check_n_atoms_in_grid(){
  int n_atoms_in_grid = 0;
  for(int i=0; i < n_all_grids; i++){
    n_atoms_in_grid += grid_n_atoms[i];
  }
  cout << "check_n_atoms_in_grid() " << n_atoms_in_grid << endl;
  return n_atoms_in_grid;
}

int NeighborSearchGrid::reorder_grid_pairs_random(){
  GridPair tmp;
  for(int i=0;i<n_grid_pairs-1;i++){
    // uniform_int_distribution was not recognized by icpc (?)
    // uniform_int_distribution<int> dist(i, n_grid_pairs-1);
    // int j = dist(random_mt);
    int j = rand()%(i-n_grid_pairs)+i;
    tmp = grid_pairs[i];
    grid_pairs[i] = grid_pairs[j];
    grid_pairs[j] = tmp;
  }
  return 0;
}
