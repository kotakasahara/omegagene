#include "MiniCell.h"
#include <bitset>

#ifdef F_CUDA
extern "C" int cuda_hostalloc_atom_info(real_pw*& h_crd, int*& h_atomids,
					real_fc*& h_work, real_fc*& h_energy,
					int n_atom_array);
extern "C" int cuda_hostalloc_cell_info(//CellPair*& h_cell_pairs, 
					//int*& h_idx_head_cell_pairs,
					int*& h_idx_xy_head_cell,
					//int max_n_cell_pairs,
					//int max_n_cells,
					int n_columns);
extern "C" int cuda_hostfree_atom_info(real_pw* h_crd, int* h_atomids,
				       real_fc*& h_work, real_fc*& h_energy);
extern "C" int cuda_hostfree_cell_info(//CellPair* h_cell_pairs,
				       //int* h_idx_head_cell_pairs,
				       int* h_idx_xy_head_cell);
#endif

MiniCell::MiniCell(){
}

MiniCell::~MiniCell(){
  free_variables();
}

int MiniCell::alloc_variables(){
  cell_crd = new int*[max_n_cells];
  for(int i=0; i < max_n_cells; i++){
    cell_crd[i] = new int[3];
  }
  idx_crd_cell = new int**[n_cells_xyz[0]];
  for(int i=0; i < n_cells_xyz[0]; i++){
    idx_crd_cell[i] = new int*[n_cells_xyz[1]];    
    for(int j=0; j < n_cells_xyz[1]; j++){
      idx_crd_cell[i][j] = new int[n_cells_xyz[2]];
      //for(int k=0; k < n_cells_xyz[2]; k++)
      //idx_crd_cell[i][j][k] = 0;
    }
  }
  //memset(idx_crd_cell, 0, sizeof(int) * max_n_cells);
  idx_atom_cell = new int[max_n_atoms_exbox + n_columns * N_ATOM_CELL];
  n_atoms_xy = new int[n_columns];
  n_cells_z = new int[n_columns];
  //crd_in_cell = new real*[max_n_atoms_exbox + n_columns * N_ATOM_CELL];
  //for (int i=0; i < max_n_atoms_exbox + n_columns * N_ATOM_CELL; i++){
  //crd_in_cell[i] = new real[2];
  //}

  idx_atom_cell_xy = new int[max_n_atoms_exbox];
  idx_xy_head_atom = new int[n_columns+1];

  //idx_cell_head_atom = new int[max_n_cells+1];
  idx_cell_n_atoms = new int[max_n_cells];
  //dummy_particles = new int*[n_columns];
  //for(int i=0; i < n_columns; i++){
  //dummy_particles[i] = new int[N_ATOM_CELL];
  
  atomids_rev = new int[get_max_n_atom_array()];
  atomids_buf = new int[get_max_n_atom_array()];
#ifdef F_CUDA
  cout << "cuda_hostalloc_atom_info"<<endl;
  cuda_hostalloc_atom_info(crd, atomids, work, energy, get_max_n_atom_array());
  cout << "cuda_hostalloc_cell_info"<<endl;
  cuda_hostalloc_cell_info(//cell_pairs, idx_head_cell_pairs,
			   idx_xy_head_cell,
			   //max_n_cell_pairs, max_n_cells+1,
			   n_columns+1);
#else
  crd = new real_pw[get_max_n_atom_array()*3];
  atomids = new int[get_max_n_atom_array()];
  work = new real_fc[get_max_n_atom_array()*3];
  energy = new real_fc[2];
  cell_pairs = new CellPair[max_n_cell_pairs];
  idx_head_cell_pairs = new int[max_n_cells+1];
  idx_xy_head_cell = new int[n_columns+1];
#endif
  // region
  region_atoms = new int*[125];
  for(int i=0; i < 125; i++){
    region_atoms[i] = new int[max_n_atoms_region[i]];
  }

  uni2cell_z = new int*[n_uni];
  for(int i=0; i < n_uni; i++){
    uni2cell_z[i] = new int[2];
  } 
  cell2uni_z = new int*[max_n_cells];
  for(int i=0; i < max_n_cells; i++){
    cell2uni_z[i] = new int[2];
  }

  return 0;
}
int MiniCell::init_variables(){
  for(int i=0; i < max_n_cells; i++)
    for(int j=0; j < 3; j++)
      cell_crd[i][j] = 0;
  for(int i=0; i < n_cells_xyz[0]; i++)
    for(int j=0; j < n_cells_xyz[1]; j++)
      for(int k=0; k < n_cells_xyz[2]; k++)
	idx_crd_cell[i][j][k] = 0;
  for(int i=0; i < get_max_n_atom_array(); i++)
    idx_atom_cell[i] = 0;
  for(int i=0; i < n_columns; i++)
    n_atoms_xy[i]=0;
  for(int i=0; i < n_columns; i++)
    n_cells_z[i] = 0;
  
  //for(int i = 0; i < (get_max_n_atom_array()); i++){
  //crd_in_cell[i][0]=0.0;
  //crd_in_cell[i][1]=0.0;
  //}
  for(int i = 0; i < get_max_n_atom_array(); i++){
    atomids[i] = -1;
  }
  for(int i = 0; i < get_max_n_atom_array(); i++){
    atomids_buf[i] = -1;
  }
  
  for(int i=0; i < max_n_atoms_exbox; i++)
    idx_atom_cell_xy[i] = -1;
  for(int i=0; i < n_columns+1; i++)
    idx_xy_head_cell[i] = 0;
  for(int i=0; i < n_columns+1; i++)
    idx_xy_head_atom[i] = 0;
  for(int i=0; i < max_n_cells; i++)
    idx_cell_n_atoms[i] = 0;
  
  //#ifdef F_CUDA
  //  cout << "cuda_hostalloc_atom_info"<<endl;
  //  cuda_hostalloc_atom_info(crd, atomids, work, energy, n_atoms);
  //#else
  //  for(int i=0; i < (n_atoms + n_columns * N_ATOM_CELL) * 3; i++)
  //crd[i] = 0.0;
  //for(int i=0; i < (n_atoms + n_columns * N_ATOM_CELL); i++)
  //atomids[i] = 0;
  //init_energy_work();
  //#endif
  return 0;
}
int MiniCell::init_energy_work(){
  for (int i=0; i<get_n_atom_array()*3; i++){
    work[i] = 0.0;
  }
  energy[0] = 0.0;
  energy[1] = 0.0;
  return 0;
}

int MiniCell::free_variables(){
  for(int i=0; i < max_n_cells; i++){
    delete[] cell_crd[i];
  }
  delete[] cell_crd;
  for(int i=0; i < n_cells_xyz[0]; i++){
    for(int j=0; j < n_cells_xyz[1]; j++){
      delete[] idx_crd_cell[i][j];
    }
    delete[] idx_crd_cell[i];
  }
  delete[] idx_crd_cell;
  delete[] idx_atom_cell;
  
  delete[] n_atoms_xy;
  delete[] n_cells_z;

  //for(int i=0; i < n_atoms; i++){
  //delete[] crd_in_cell[i];
  //}
  //delete[] crd_in_cell;
  delete[] idx_atom_cell_xy;
  delete[] idx_xy_head_atom;
  //delete[] idx_cell_head_atom;
  delete[] idx_cell_n_atoms;

  //for(int i=0; i < n_columns; i++){
    //delete[] dummy_particles[i];
  //}
  //delete[] dummy_particles;

  delete[] atomids_rev;
  delete[] atomids_buf;
#ifdef F_CUDA
  cuda_hostfree_atom_info(crd, atomids, work, energy);
  cuda_hostfree_cell_info(//cell_pairs, idx_head_cell_pairs,
			  idx_xy_head_cell);
#else
  delete[] crd;
  delete[] work;
  delete[] atomids;
  delete[] energy;
  delete[] idx_xy_head_cell;
  delete[] cell_pairs;
  delete[] idx_head_cell_pairs;
#endif

  for(int i=0; i < 125; i++){
    delete[] region_atoms[i];
  }
  delete[] region_atoms;

  for(int i=0; i<n_uni; i++){
    delete[] uni2cell_z[i];
  }
  delete[] uni2cell_z;

  for(int i=0; i<max_n_cells; i++){
    delete[] cell2uni_z[i];
  }
  delete[] cell2uni_z;
  return 0;
}

int MiniCell::set_grid_xy(){
  // set 
  //   n_cells_xyz
  //   L_cell_xy
  //   exbox_l

  //cout << "set_grid_xy()" << endl;
  real_pw volume = exbox_l[0] * exbox_l[1] * exbox_l[2];

  real_pw min_cell_width_xy = pow((real_pw)((real_pw)N_ATOM_CELL / ((real_pw)n_atoms/volume )), (real_pw)(1.0/3.0));
  //cout << "min_cell_width_xy : " << min_cell_width_xy << endl;;
  //cout << " ... " << volume << ", "  << (real)n_atoms/volume << ", " << (real)N_ATOM_CELL / ((real)(n_atoms)/pbc->cal_volume()) << endl;
  n_cells_xyz[0] = floor(exbox_l[0] / min_cell_width_xy);
  n_cells_xyz[1] = floor(exbox_l[1] / min_cell_width_xy);
  L_cell_xy[0] = exbox_l[0] / (real_pw)n_cells_xyz[0];
  L_cell_xy[1] = exbox_l[1] / (real_pw)n_cells_xyz[1];
  L_cell_xy[2] = 0.0;
  n_columns = (n_cells_xyz[0] * n_cells_xyz[1]);

  // the maximum number of atoms in each column is estimated as
  //  1.3 times larger than the expected number. 
  //  there is no data the factor 1.3 is good or not.
  //  this point should be tested.
  max_n_atoms_column = ceil((real_pw)max_n_atoms_exbox/(real_pw)n_columns * 1.3);

  //cout << "max_n_atoms_column " << max_n_atoms_column << " max_n_atoms_exbox: " << max_n_atoms_exbox << endl;
  n_cells_xyz[2] = (max_n_atoms_column-N_ATOM_CELL+1)/N_ATOM_CELL;

  n_neighbors_xy[0] = ceil(cutoff_pair/L_cell_xy[0]);
  n_neighbors_xy[1] = ceil(cutoff_pair/L_cell_xy[1]);

  //cout << " N_columns: " << n_columns << "  max_n_atoms_column: " << max_n_atoms_column <<endl;
  //cout << " XY grid (n_cell): " << n_cells_xyz[0] << " - " << n_cells_xyz[1]<< " - " << n_cells_xyz[2] << endl;
  //cout << " XY grid (L_cell): " << L_cell_xy[0] << " - " << L_cell_xy[1] << endl;
  //cout << " neighbors: " << n_neighbors_xy[0] << " - " << n_neighbors_xy[1] << endl;
  // 
  max_n_cells = n_cells_xyz[0] * n_cells_xyz[1] * n_cells_xyz[2];
  max_n_cell_pairs = ((n_neighbors_xy[0]*2+1) *
		      (n_neighbors_xy[1]*2+1) * (n_neighbors_xy[0]+n_neighbors_xy[1]+2) )
    * 0.5 * max_n_cells;
  //max_n_cell_pairs = max_n_cells * max_n_cells;
  //cout << "max_n_cell_pairs : " << max_n_cell_pairs << endl;

  n_cell_pairs = 0;

  n_uni = n_columns * n_uni_z;

  return 0;
}

int MiniCell::set_atoms_into_grid_xy(){
  // set variables
  //   cell_crd
  //   idx_atom_cell_xy
  //   n_atoms_xy
  //   idx_xy_head_atom
  //   L_cell_xy[3]
  // reorder atomids

  // idx_cell_atom
  //   idx_cell_crd_xy_atoms[cell xy id] = (atom1, atom2, ...)
  //    crd x = ID%n_cells_xyz[0]
  //    crd y = ID/n_cells_xyz[0]

  //int **idx_cell_xy_atom;
  //idx_cell_xy_atom = new int[n_columns];
  //for(int i=0; i < n_columns; i++){
  //    idx_cell_xy_atom[i] = new int[max_n_atoms_column];
  //  }

  //memset(idx_cell_xy_atom, 0, sizeof(int) * n_columns * max_n_atoms_column);
  //memset(n_atoms_xy, 0, sizeof(int) * n_columns);

  for(int i=0; i<n_columns; i++) n_atoms_xy[i] = 0;
  //memset(n_atoms_xy, 0, sizeof(int) * n_columns);

  // Assigning each atom to each column (column means XY grid cell)
  for (int atom_id_g = 0; atom_id_g < n_atoms_exbox; atom_id_g++){
    int cell_crd_xy[2];

    for (int d=0; d < 2; d++){
      
      if(crd[atom_id_g*3+d] > exbox_upper_bound[d] ||
	 crd[atom_id_g*3+d] < exbox_lower_bound[d]){
	cout << "CRDERROR! EXCEEDED PB: " << atom_id_g << "-" << d << " "
	     << crd[atom_id_g*3+d] << " , " << box_upper_bound[d] << " , " 
	     << box_lower_bound[d] << endl;
      }

      real_pw x = crd[atom_id_g*3+d] - exbox_lower_bound[d];
      while(x < 0.0) x += exbox_l[d];
      while(x >= exbox_l[d]) x -= exbox_l[d];

      cell_crd_xy[d] = floor(x / L_cell_xy[d]);
      //crd_in_cell[atom_id_g][d] = (x - L_cell_xy[d] * cell_crd_xy[d]) / L_cell_xy[d];
    }

    //int column_id = cell_crd_xy[0] + cell_crd_xy[1] * n_cells_xyz[0];
    int column_id = get_column_id_from_crd(cell_crd_xy[0], cell_crd_xy[1]);

    idx_atom_cell_xy[atom_id_g] = column_id;
    //idx_cell_xy_atom[column_id][n_atoms_xy[column_id]] = atom_id;
    n_atoms_xy[column_id]++;
  }

  //cout << "set_idx_xy_head_atom_from_n_atoms_xy()" << endl;
  set_idx_xy_head_atom_from_n_atoms_xy();  
  
  set_atomids_from_buf();

  //cout << "reorder_atominfo_for_columns()" << endl;
  reorder_atominfo_for_columns();

  //cout << "set_dummy_atoms()" << endl;
  set_dummy_atoms();

  //cout << "reorder_atominfo_in_columns() " <<endl;
  reorder_atominfo_in_columns();
  set_n_neighbors_z_cells();
  set_atomids_rev();
  //cout << " // set_atomids_rev"<<endl;
  //debug_set_atoms_into_grid();

   // TEST
  /*
 int cid = 1202;
   cout << "test " <<cid<< endl;
   cout << cell_crd[cid][0] <<"-"<<cell_crd[cid][1] <<"-"<<cell_crd[cid][2] <<endl;
   for(int i=0; i < N_ATOM_CELL; i++){
     cout << (cid*N_ATOM_CELL+i) << "- " 
 	 << crd[(cid*N_ATOM_CELL+i)*3] << ", " 
 	 << crd[(cid*N_ATOM_CELL+i)*3+1]<< ", " 
 	 << crd[(cid*N_ATOM_CELL+i)*3+2]<<endl;
   }
   cid = 1976;
   cout << "test " <<cid<< endl;
   cout << cell_crd[cid][0] <<"-"<<cell_crd[cid][1] <<"-"<<cell_crd[cid][2] <<endl;
   for(int i=0; i < N_ATOM_CELL; i++){
   cout << (cid*N_ATOM_CELL+i) << "- " 
   << crd[(cid*N_ATOM_CELL+i)*3] << ", " 
   << crd[(cid*N_ATOM_CELL+i)*3+1]<< ", " 
   << crd[(cid*N_ATOM_CELL+i)*3+2]<<endl;
   }
  */
  return 0;
}

int MiniCell::reorder_atominfo_for_columns(){
  // required
  //   idx_atom_cell_xy
  //   idx_xy_head_cell
  // re set
  //   n_atoms_xy
  for (int i_col = 0; i_col < n_columns; i_col++){
    n_atoms_xy[i_col] = 0;
  }
  real_pw tmp1_x[3];
  real_pw tmp2_x[3];
  int tmp1_atomid;
  int tmp2_atomid;
  int tmp1_column;
  int tmp2_column;
  int atomid_g = 0;
  tmp2_column = idx_atom_cell_xy[atomid_g];
  int to_atomid_g = idx_xy_head_atom[tmp2_column] + n_atoms_xy[tmp2_column];
  tmp2_atomid = atomids[atomid_g];
  tmp2_x[0] = crd[atomid_g*3];
  tmp2_x[1] = crd[atomid_g*3+1];
  tmp2_x[2] = crd[atomid_g*3+2];
  idx_atom_cell_xy[atomid_g] = -1;
  int idx_to_atomid_g = to_atomid_g*3;
  for(int n_reordered = 1; n_reordered < n_atoms_exbox; n_reordered++){
    tmp1_x[0] = crd[idx_to_atomid_g+0];
    tmp1_x[1] = crd[idx_to_atomid_g+1];
    tmp1_x[2] = crd[idx_to_atomid_g+2];
    crd[idx_to_atomid_g+0] = tmp2_x[0];
    crd[idx_to_atomid_g+1] = tmp2_x[1];
    crd[idx_to_atomid_g+2] = tmp2_x[2];
    tmp1_atomid = atomids[to_atomid_g];
    atomids[to_atomid_g] = tmp2_atomid;
    tmp1_column = idx_atom_cell_xy[to_atomid_g];
    idx_atom_cell_xy[to_atomid_g] = tmp2_column;
    n_atoms_xy[tmp2_column]++;
    atomid_g = to_atomid_g;
    // seek atom which has not sorted yet
    if(atomid_g >= n_atoms_exbox || tmp1_column == -1){
      bool flg_seek = true;
      for(int tmp_col = 0; flg_seek && tmp_col < n_columns; tmp_col++){
	int tmp_atomid_g;
	for(tmp_atomid_g = idx_xy_head_atom[tmp_col] + n_atoms_xy[tmp_col];
	    flg_seek && tmp_atomid_g < idx_xy_head_atom[tmp_col+1];
	    tmp_atomid_g++){
	  if(idx_atom_cell_xy[tmp_atomid_g] >= 0){
	    atomid_g = tmp_atomid_g;
	    flg_seek = false;
	  }
	}
      }
      tmp1_column = idx_atom_cell_xy[atomid_g];
      idx_atom_cell_xy[atomid_g] = -1;
      tmp1_x[0] = crd[atomid_g*3];
      tmp1_x[1] = crd[atomid_g*3+1];
      tmp1_x[2] = crd[atomid_g*3+2];
      tmp1_atomid = atomids[atomid_g];
    }
    to_atomid_g = idx_xy_head_atom[tmp1_column] + n_atoms_xy[tmp1_column];
    idx_to_atomid_g = to_atomid_g*3;
    // set next info
    tmp2_column = tmp1_column;
    tmp2_x[0] = tmp1_x[0];
    tmp2_x[1] = tmp1_x[1];
    tmp2_x[2] = tmp1_x[2];
    tmp2_atomid = tmp1_atomid;
  }

  crd[idx_to_atomid_g+0] = tmp2_x[0];
  crd[idx_to_atomid_g+1] = tmp2_x[1];
  crd[idx_to_atomid_g+2] = tmp2_x[2];
  idx_atom_cell_xy[to_atomid_g] = tmp2_column;
  n_atoms_xy[tmp2_column]++;
  atomids[to_atomid_g] = tmp2_atomid;

  //for (int i_col = 0; i_col < n_columns; i_col++){
  //n_atoms_xy[i_col] = 0;
  //}

  return 0;
}

int MiniCell::set_dummy_atoms(){
  int valid_atomid_g = -1;
  for(int i_col = 0; i_col < n_columns; i_col++){
    for(int i_atom = 0; i_atom < n_atoms_xy[i_col]; i_atom++){
      int atomid_g = idx_xy_head_atom[i_col] + i_atom;
      if(idx_atom_cell_xy[atomid_g] >= 0){
	valid_atomid_g = atomid_g;
      }else{
	idx_atom_cell_xy[atomid_g] = i_col;
	//crd[atomid_g*3]   = crd[valid_atomid_g*3];
	//crd[atomid_g*3+1] = crd[valid_atomid_g*3+1];
	//crd[atomid_g*3+2] = crd[valid_atomid_g*3+2];
	atomids[atomid_g] = -1;
      }
    }
  }
  return 0;
}
int MiniCell::reorder_atominfo_in_columns(){
  //set
  //   atomids
  //require
  //  idx_atom_cell_xy
  //  idx_xy_head_atom
  //  idx_xy_head_cell

  //for(int i=0; i<n_columns; i++) n_atoms_xy[i]=0;
  //memset(n_atoms_xy, 0, sizeof(int)*n_columns);
  
  // set atomids[atom_id_grid] = atomids_buf[atom_id_g]
  //for(int atom_id_g=0; atom_id_g < n_atoms_exbox; atom_id_g++){
  //int column_id = idx_atom_cell_xy[atom_id_g];
  //int atom_id_grid = idx_xy_head_atom[column_id] + n_atoms_xy[column_id];
  //atomids[atom_id_grid] = atomids_buf[atom_id_g];
  //if(atom_id_grid >= n_atoms_exbox+n_columns*N_ATOM_CELL ||
  //       atom_id_g >= n_atoms_exbox ||
  //column_id >= n_columns ||
  //n_atoms_xy[column_id] >= max_n_atoms_column){
  //cout << "atomids["<<atom_id_grid<<"/" << n_atoms_exbox+n_columns*N_ATOM_CELL
  //	   <<"]="<<atomids_buf[atom_id_g]<< "/" << n_atoms 
  //<<" column:"<<column_id<<"/"<<n_columns
  //<<" n_atoms_xy:"<<n_atoms_xy[column_id]<<"/"<<max_n_atoms_column<<endl;
  //exit(1);
  //}
  //n_atoms_xy[column_id]++;
  //}
  //cout << "reorder_atominfo 2:" <<endl; 
   //for(int atom_id_grid=0; atom_id_grid < n_atoms; atom_id_grid++){
   //cout << atom_id_grid << " " <<  atomids[atom_id_grid] << endl;
   //}
   //cout << endl;
  for(int column_id=0; column_id < n_columns; column_id++){  
    // sort atomids in column by order of z coordinates
    quick_sort_in_columns(idx_xy_head_atom[column_id],
 			  idx_xy_head_atom[column_id] + n_atoms_xy[column_id] - 1);
    //int last_ = atomids[idx_xy_head_atom[column_id] + n_atoms_xy[column_id] - 1];
    //for(int i=0;i < n_atoms_xy[column_id]; i++){
    //int aid = idx_xy_head_atom[column_id] + i;
    //exit(1);
    int last_id = idx_xy_head_atom[column_id] + n_atoms_xy[column_id] -1;
    
    // set crd 
    for(int i = idx_xy_head_atom[column_id] + n_atoms_xy[column_id];
	i <  idx_xy_head_atom[column_id] + n_cells_z[column_id] * N_ATOM_CELL; i++){
      crd[i*3+0] = crd[last_id*3+0];
      crd[i*3+1] = crd[last_id*3+1];
      crd[i*3+2] = crd[last_id*3+2];
      atomids[i] = -1;
    }
      //int last = 0;
      //for(int i=0; i < n_cells_z[column_id] * N_ATOM_CELL; i++){// n_atoms_xy[column_id]; i++){
      //int atom_id_grid = idx_xy_head_atom[column_id] + i;	
      //if(i < n_atoms_xy[column_id]) {last = atom_id_grid;}
      //else cout << "last : "<< last << endl;;
      //crd[atom_id_grid*3] = in_crd[atomids[last]][0];
      //crd[atom_id_grid*3+1] = in_crd[atomids[last]][1];
      //crd[atom_id_grid*3+2] = in_crd[atomids[last]][2];
      /*
      //       int cid = idx_xy_head_cell[column_id] + i/N_ATOM_CELL;
      //       cout << "set: cell:" << cid << "["
      // 	   << cell_crd[cid][0] <<"-"
      // 	   << cell_crd[cid][1] <<"-"
      // 	   << cell_crd[cid][2] <<"] "
      // 	   << " atom:" << atom_id_grid << " "
      // 	   << " (" << N_ATOM_CELL*cid << ") "
      // 	   <<crd[atom_id_grid*3] << ", "
      // 	   <<crd[atom_id_grid*3+1] << ", "
      // 	   <<crd[atom_id_grid*3+2]
      // 	   << endl;
      //       */
      //crd[atom_id_grid*3] = (real)in_crd[atomids[atom_id_grid]][0];
      //crd[atom_id_grid*3+1] = (real)in_crd[atomids[atom_id_grid]][1];
      //crd[atom_id_grid*3+2] = (real)in_crd[atomids[atom_id_grid]][2];
    for(int i=0; i < n_cells_z[column_id] * N_ATOM_CELL; i++){
      int atom_id_grid = idx_xy_head_atom[column_id] + i;
      idx_atom_cell[atom_id_grid] = idx_xy_head_cell[column_id] + i / N_ATOM_CELL;
    }
  }
  n_neighbors_xy[2] = ceil(cutoff_pair/L_cell_xy[2]);
  //cout << "L_cell_xy[2]:" << L_cell_xy[2] << " neighbors[2]:"<<n_neighbors_xy[2]<< endl;

  return 0;
}
int MiniCell::set_n_neighbors_z_cells(){
  real_pw lower_z = exbox_lower_bound[2];    
  L_cell_xy[2] = 1e10;
  real_pw ave_n_cells = 0.0;
  for(int column_id=0; column_id < n_columns; column_id++){  
    ave_n_cells += n_cells_z[column_id];
  }
  ave_n_cells /= n_columns;
  real_pw ave_l_cell_z = exbox_l[2] / ave_n_cells;
  L_cell_xy[2] = ave_l_cell_z;
  n_neighbors_xy[2] = ceil(cutoff_pair/ave_l_cell_z * COEF_MAX_N_NEIGHBORS);
  return n_neighbors_xy[2];
}
int MiniCell::swap_atomid_crd(const int i, const int j){
  //  cout << "swap " << i << "-" << j << endl;
  int t;
  real_pw x0, x1, x2;
  t = atomids[i];
  atomids[i] = atomids[j];
  atomids[j] = t;
  x0 = crd[i*3+0];
  x1 = crd[i*3+1];
  x2 = crd[i*3+2];
  crd[i*3+0] = crd[j*3+0];
  crd[i*3+1] = crd[j*3+1];
  crd[i*3+2] = crd[j*3+2];
  crd[j*3+0] = x0;
  crd[j*3+1] = x1;
  crd[j*3+2] = x2;
  //cout << " crd : " << i << " " << crd[i*3+2] << " " << j << " " << crd[j*3+2] << endl;
  return 0;
}
int MiniCell::quick_sort_in_columns(const int l, const int r){
  if (r > l){
    real_pw v = crd[r*3+2];
    int i = l-1;
    int j = r+1;
    int t;
    //cout << "qs " << l << " - " << r << " v: " << v <<endl;
    for (;;){
      while(crd[(++i)*3+2] < v);// cout << "qs1.1 " << i  << " - " << crd[i*3+2] << endl;
      while(crd[(--j)*3+2] > v && j > i);// cout << "qs1.2 " << j  << " - " << crd[j*3+2] << endl;
      if(i>=j) break;
      swap_atomid_crd(i,j);
    }
    //swap_atomid_crd(i,j);
    //cout << "qs2 " << l << " - " << i-1 << " " << i << " - " << r << endl;
    quick_sort_in_columns(l, i-1);
    quick_sort_in_columns(i, r);
  }

  /*
  cout << "qs " << l << " - " << r << endl;
  if (r > l){
    
    int pivot = l;
    cout << "pivot : " << pivot << "\t" << atomids[pivot] << endl;
    cout << "crd " << crd[pivot*3+2] << endl;    
    real v = crd[pivot*3+2];
    
    int i = l-1;
    int j = r+1;

    for (;;){
      while(i<r and crd[(++i)*3+2] < v) cout << "i-j : " << i << " " << j << endl;;
      while(j>l and crd[(--j)*3+2] >= v) cout << "j-i : " << j << " "  << i << endl;;
      cout << " qs2 : " << i << " - " << j << endl;
      if(i>=j) break;
      swap_atomid_crd(i, j);
      cout << " swp " << endl;
    }
    cout << " qs3 " << l << " - " << i-1 << " " << i << " - " << r << endl;
    //swap_atomid_crd(i, r);    
    //int t = atomids[i];
    //atomids[i] = atomids[r];
    //atomids[r] = t;
    if (l < i-1){
      quick_sort_in_columns(l, i-1);
      cout << " qs4 " << endl;
      if (i < r) quick_sort_in_columns(i, r);
      cout << " qs5 " << endl;
    }
  }
  */
  return 0;
}

int MiniCell::debug_set_atoms_into_grid(){
  int *rev_atomids = new int[n_atoms];
  for(int i = 0; i< n_atoms; i++){ rev_atomids[i] = -1; }   
  cout << "!!!!!!!!!!! DEBUG" << endl;
  for(int atom_id_grid=0; atom_id_grid < N_ATOM_CELL*n_cells; atom_id_grid++){
    int cell_id = idx_atom_cell[atom_id_grid];
    if(atomids[atom_id_grid] != -1){
      if(rev_atomids[atomids[atom_id_grid]] != -1){
 	cout << "!!! " << atomids[atom_id_grid] << " " << atom_id_grid << " "
 	     << rev_atomids[atomids[atom_id_grid]]  << endl;
      }
      rev_atomids[atomids[atom_id_grid]] = atom_id_grid;
    }else{
    }
    cout << atom_id_grid << " " << atomids[atom_id_grid] << " ("<<crd[atom_id_grid*3]<<","
 	 << crd[atom_id_grid*3+1] << "," << crd[atom_id_grid*3+2]
 	 << ") cell:" << cell_id << ":"
 	 << cell_crd[cell_id][0] << "," << cell_crd[cell_id][1] << ","
 	 << cell_crd[cell_id][2] << " " 
 	 <<endl;
  }
  //cout << "atomids rev : ";
  //for(int i = 0; i< n_atoms; i++){
  //cout << " " << rev_atomids[i] ;
  //}
  //cout << endl;
  delete[] rev_atomids;
  return 0;
}

int MiniCell::update_crd(real** in_crd){
  /*

  for(int column_id=0; column_id < n_columns; column_id++){  
    int last = 0;
    for(int i=0; i < n_cells_z[column_id] * N_ATOM_CELL; i++){
      int atom_id_grid = idx_xy_head_atom[column_id] + i;	
      if(i < n_atoms_xy[column_id]) {last = atom_id_grid;}
      crd[atom_id_grid*3] = in_crd[atomids[last]][0];
      crd[atom_id_grid*3+1] = in_crd[atomids[last]][1];
      crd[atom_id_grid*3+2] = in_crd[atomids[last]][2];
    }
    }*/

  int atomid = 0;
  int i_grid = 0;
  for(int atomid_grid=0;
      atomid_grid < get_n_atom_array();
      atomid_grid++, i_grid+=3){
    if(atomids[atomid_grid] >= 0)  atomid = atomids[atomid_grid];
    crd[i_grid] = in_crd[atomid][0];
    crd[i_grid+1] = in_crd[atomid][1];
    crd[i_grid+2] = in_crd[atomid][2];
  }

  return 0;
}
/*
int MiniCell::update_cell_assign(const real** in_crd){
  // update variables
  //   idx_atom_cell_xy
  //   n_atoms_xy
  //   idx_xy_head_atom
  
  for (int atom_id=0; atom_id < n_atoms; atom_id++){
    int mg[2] = {0,0};
    for (int d=0; d < 2; d++){
      if (crd_in_cell[atom_id][d] < 0.0){
 	crd_in_cell[atom_id][d] += 1.0;
 	mg[d] -= 1;
      }else if(crd_in_cell[atom_id][d] > 1.0){
 	crd_in_cell[atom_id][d] -= 1.0;
 	mg[d] += 1;
      }
    }
    int fg[2]; //from
    int tg[2]; //to
    int from  = idx_atom_cell_xy[atom_id];
    fg[0] = from%n_cells_xyz[0];
    fg[1] = from/n_cells_xyz[0];
    for(int d=0; d < 2; d++){
      tg[d] = fg[d] + mg[d];
      if(tg[d] == -1) tg[d] = n_cells_xyz[d] - 1;
      else if(tg[d] == n_cells_xyz[d]) tg[d] = 0;
    }
    int to = tg[1]*n_cells_xyz[0] + tg[0];
    idx_atom_cell_xy[atom_id] = to;
    n_atoms_xy[from]--;
    n_atoms_xy[to]++;
    
  }
  set_idx_xy_head_atom_from_n_atoms_xy();
  reorder_atominfo(in_crd);  
  //reset_cell_assignment();
  return 0;
}
*/
int MiniCell::set_idx_xy_head_atom_from_n_atoms_xy(){
  // set 
  //   idx_xy_head_atom
  //   n_cells_z
  //   n_cells
  //   cell_crd
  //   idx_crd_cell
  //   L_cell_xy[2]
  // reorder crd, atomids
  // require
  //   n_atoms_xy

  n_cells = 0;
  int atom_index=0;
  int cur_cell_id = 0;
  int n_dummies_all = 0;
  //cout << "set_idx_xy_head_atom_from_n_atoms_xy" <<endl;
  for(int column_id = 0; column_id < n_columns; column_id++){
    //    cout << "column_id " << column_id << endl;
    int column_x = column_id%n_cells_xyz[0];
    int column_y = column_id/n_cells_xyz[0];
    //real_pw lower_z = exbox_lower_bound[2];
    //cout << "lower_z col:" << column_id << " " << lower_z <<endl;
    idx_xy_head_cell[column_id] = cur_cell_id;
    idx_xy_head_atom[column_id] = atom_index;
    n_cells_z[column_id] = (n_atoms_xy[column_id] + N_ATOM_CELL - 1)/N_ATOM_CELL;
    n_cells += n_cells_z[column_id];
    for (int cur_cell_z=0; cur_cell_z < n_cells_z[column_id]; cur_cell_z++, cur_cell_id++){
      // cout << "set " << cur_cell_id << "("<<column_x<<","<<column_y<<","<<cur_cell_z<<")"<<endl;
      idx_crd_cell[column_x][column_y][cur_cell_z] = cur_cell_id;
      cell_crd[cur_cell_id][0] = column_x;
      cell_crd[cur_cell_id][1] = column_y;
      cell_crd[cur_cell_id][2] = cur_cell_z;
    }
    int n_dummies = n_cells_z[column_id]*N_ATOM_CELL - n_atoms_xy[column_id];
    n_dummies_all += n_dummies;
    int i;
    
    //for(i = 0; i < n_dummies; i++){
      //cout << "dummy " << column_id << "-" << n_atoms_xy[column_id] << " " << i <<" " << idx_xy_head_atom[column_id] << endl;
    //int cur_atomid_grid = idx_xy_head_atom[column_id]+n_atoms_xy[column_id]+i;
    //atomids[cur_atomid_grid] = -1;
    //}
    //dummy_particles[column_id][i] = atom_index + n_atoms_xy[column_id] + i;
    //atom_index += n_atoms_xy[column_id];
    atom_index += n_cells_z[column_id]*N_ATOM_CELL;
    //cout << "idx_xy_head_atom["<<column_id<<"] = " << atom_index << " , " << n_atoms_xy[column_id]<<endl;
  }
  idx_xy_head_cell[n_columns] = n_cells;
  idx_xy_head_atom[n_columns] = atom_index;
  //for(int cell_id=0; cell_id < n_cells; cell_id++)
  //idx_cell_head_atom[cell_id] = N_ATOM_CELL*cell_id;
  //cout << "n_dummies_all: " << n_dummies_all << endl;

  return 0;
}

int MiniCell::set_atomids_rev(){
  for(int atom_id_grid = 0; 
      atom_id_grid < get_n_atom_array(); 
      atom_id_grid++){
    if(atomids[atom_id_grid] >= 0){
      atomids_rev[atomids[atom_id_grid]] = atom_id_grid;
    }
    //cout << "set_atomids_rev " << atomids[atom_id_grid] << " " << atom_id_grid << endl;
  }
  return 0;
}

int MiniCell::set_atomids_buf(){
  for (int i=0; i < get_n_atom_array(); i++){
    atomids_buf[i] = atomids[i];
  }
  return 0;
}
int MiniCell::set_atomids_from_buf(){
  for (int i=0; i < get_max_n_atom_array(); i++){
    atomids[i] = atomids_buf[i];
  }
  return 0;
}

real_pw MiniCell::get_cell_z_min(int cell_id){  
  //return crd[(idx_cell_head_atom[cell_id])*3+2];
  return crd[cell_id * N_ATOM_CELL * 3 + 2];
}
real_pw MiniCell::get_cell_z_max(int cell_id){
  //cout << "get_cell_z_max " << cell_id <<" "<< cell_id*N_ATOM_CELL+N_ATOM_CELL-1
  //<<" " << (cell_id*N_ATOM_CELL+N_ATOM_CELL-1)*3+2<<" "<<n_atoms*3+n_columns*N_ATOM_CELL*3 << endl;
  //cout << "  "<< crd[(cell_id*N_ATOM_CELL+N_ATOM_CELL-1)*3+2] << endl;
  return crd[(cell_id * N_ATOM_CELL + N_ATOM_CELL-1) * 3 + 2];
  //return crd[(idx_cell_head_atom[cell_id+1]-1)*3+2];
}

int MiniCell::set_uniform_grid(){
  for(int i=0; i<n_uni; i++){
    uni2cell_z[i][0] = -1;
    uni2cell_z[i][1] = -1;
  }
  for(int i=0; i<max_n_cells; i++){
    cell2uni_z[i][0] = -1;
    cell2uni_z[i][1] = -1;
  }

  int cell[3];
  for(int cell_id=0; cell_id < n_cells; cell_id++){
    cell[0] = cell_crd[cell_id][0];
    cell[1] = cell_crd[cell_id][1];
    cell[2] = cell_crd[cell_id][2];
    real_pw cell_z_min = get_cell_z_min(cell_id) - pbc->lower_bound[2];
    real_pw cell_z_max = get_cell_z_max(cell_id) - pbc->lower_bound[2];
    int uni_z_min = cell_z_min/L_z_uni;
    int uni_z_max = cell_z_max/L_z_uni;    

    cell2uni_z[cell_id][0] = uni_z_min;
    cell2uni_z[cell_id][1] = uni_z_max;
    for(int z = uni_z_min; z <= uni_z_max; z++){
      int uni_id = get_uni_id_from_crd(cell[0], cell[1], z);
      if(uni2cell_z[uni_id][0] < 0 || cell[2] < uni2cell_z[uni_id][0])
	uni2cell_z[uni_id][0] = cell[2];
      if(cell[2] > uni2cell_z[uni_id][1])
	uni2cell_z[uni_id][1] = cell[2];
      
    }
  }

  return 0;
}

int MiniCell::enumerate_cell_pairs(){
  // set
  //   cel_pairs
  //   n_cell_pairs
  // require
  //   cell_crd
  //   idx_crd_cell

  //cout << "enumerate_cell_pairs()"<<endl;
  set_uniform_grid();

  n_cell_pairs = 0;
  int cell1[3];
  int tmp_cell1 = 0;
  int tmp_column_pairs = 0;
  int tmp_n_cell_pair = 0;
  for(int cell1_id=0; cell1_id < n_cells; cell1_id++){
    cell1[0] = cell_crd[cell1_id][0];
    cell1[1] = cell_crd[cell1_id][1];
    cell1[2] = cell_crd[cell1_id][2];
    int column1_id = get_column_id_from_crd(cell1[0], cell1[1]);
    //real_pw cell1_z_max = get_cell_z_max(cell1_id);
    //real_pw cell1_z_min = get_cell_z_min(cell1_id);

    idx_head_cell_pairs[cell1_id] = n_cell_pairs;
    //cout << "dbg idx_head["<< cell1_id<< "] ("<<cell1[0]<<","<<cell1[1]
    //<<","<<cell1[2]<<") = " << n_cell_pairs << endl;
    //real_pw cell1_z_min = get_cell_z_min(cell1_id);
    //real_pw cell1_z_max = get_cell_z_max(cell1_id);
    bool cell1_odd = cell1_id%2!=0;
    
    int d_cell[3];
    int cell_rel[3];
    int image[3];
    int cell2[3];
    //cout << "dbg 4" << endl;
    for(d_cell[0] = -n_neighbors_xy[0];
	d_cell[0] <= n_neighbors_xy[0]; d_cell[0]++){
      cell_rel[0] = cell1[0] + d_cell[0];
      //cout << "dbg 5 " << d_cell[0] << endl;
      ////  val
      //real_pw dx = 0.0;
      //if (d_cell[0] > 1)        dx = (d_cell[0] - 1) * L_cell_xy[0];
      //else if (d_cell[0] < -1)  dx = (d_cell[0] + 1) * L_cell_xy[0];
      //real_pw dx2 = dx * dx;
      image[0] = 0; cell2[0] = cell_rel[0];
      if(cell_rel[0] < 0){
	image[0] = -1; cell2[0] = n_cells_xyz[0] + cell_rel[0];
      }else if (cell_rel[0] >= n_cells_xyz[0]){
	image[0] = 1;  cell2[0] = cell_rel[0] - n_cells_xyz[0];
      }


      ////  --val
      for(d_cell[1] = -n_neighbors_xy[1];
	  d_cell[1] <= n_neighbors_xy[1]; d_cell[1]++){
	cell_rel[1] = cell1[1] + d_cell[1];
	//cout << "dbg 30 " << d_cell[1] << endl;
	////  val
	//real_pw dy = 0.0;
	//if (d_cell[1] > 1)        dy = (d_cell[1] - 1) * L_cell_xy[1];
	//else if (d_cell[1] < -1)  dy = (d_cell[1] + 1) * L_cell_xy[1];
	//real_pw dy2 = dy * dy;
	image[1] = 0; cell2[1] = cell_rel[1];
	if(cell_rel[1] < 0){
	  image[1] = -1; cell2[1] = n_cells_xyz[1] + cell_rel[1];
	}else if (cell_rel[1] >= n_cells_xyz[1]){
	  image[1] = 1;  cell2[1] = cell_rel[1] - n_cells_xyz[1];
	}
	////  --val
	int column2_id = cell2[0] + cell2[1] * n_cells_xyz[0];
	tmp_column_pairs++;
	
	
	
	// 1. i1 ... image = -1
	// 2. i2 ... image = 0
	// 3. i3 ... image = +1
	int first_uni_z[3] = {0, 0, 0};
	int last_uni_z[3] = {0, 0, 0};
	int tmp_first = cell2uni_z[cell1_id][0] - 2;
	if(tmp_first < 0){
	  first_uni_z[0] = n_uni_z + tmp_first;
	  last_uni_z[0] = n_uni_z - 1;
	  tmp_first = 0;
	}else{
	  first_uni_z[0] = -1;
	  last_uni_z[0] = -1;
	}
	first_uni_z[1] = tmp_first;
	int tmp_last = cell2uni_z[cell1_id][1] + 2;
	if(tmp_last >= n_uni_z){
	  first_uni_z[2] = 0;
	  last_uni_z[2] = tmp_last - n_uni_z;
	  last_uni_z[1] = n_uni_z - 1;
	}else{
	  first_uni_z[2] = -1;
	  last_uni_z[2] = -1;	  
	  last_uni_z[1] = tmp_last;
	}
	
	int tmp_img = -2;
	for(int i_img = 0; i_img < 3; i_img++){
	  image[2] = i_img-1;	  
	  if (first_uni_z[i_img] < 0) continue;
	  int first_uni_id = get_uni_id_from_crd(cell2[0], cell2[1], 
						 first_uni_z[i_img]);
	  int last_uni_id = get_uni_id_from_crd(cell2[0], cell2[1], 
						last_uni_z[i_img]);
	  int first_cell = get_cell_id_from_crd(cell2[0], cell2[1], uni2cell_z[first_uni_id][0]);
	  int last_cell = get_cell_id_from_crd(cell2[0], cell2[1], uni2cell_z[last_uni_id][1]);
	  
	  for(int cell2_id = first_cell;
	      cell2_id <= last_cell; cell2_id++){
	    //real_pw cell2_z_max = get_cell_z_max(cell2_id);
	    //real_pw cell2_z_min = get_cell_z_min(cell2_id);	    
	    //real_pw dz = cell2_z_max - cell1_z_min;
	    //if((cell2_z_max <= cell1_z_max && cell2_z_max >= cell1_z_min) ||
	    //(cell2_z_min <= cell1_z_max && cell2_z_min >= cell1_z_min))
	    //dz = 0;
	    //real_pw dz2 = dz*dz;
	    //real_pw dz_bt_up = cell2_z_min - cell1_z_max;
	    //real_pw dz_bt_up2 = dz_bt_up*dz_bt_up;
	    //bool flg_rev = true;
	    //if(dz2 > dz_bt_up2) {flg_rev=false; dz2 = dz_bt_up2; }
	             
	    //if(dx2+dy2+dz2 > cutoff_pair2){
	      //if(!flg_rev) break;
	      //continue;
	      //}

	    bool cell2_odd = cell2_id%2!=0;
	    if(check_valid_pair(cell1_id, cell2_id, cell1_odd, cell2_odd))
	      add_cell_pair(cell1_id, cell2_id, image);	    
	  }
	}
      }
    }
  }
  idx_head_cell_pairs[n_cells] = n_cell_pairs;
  //cout << "checked cells: " << tmp_cell1 << " / " << n_cells << endl ;
  //cout << "checked column_pairs: " << tmp_column_pairs << " / " 
  //<< n_columns * (n_neighbors_xy[0]*2+1)*(n_neighbors_xy[1]*2+1)   <<" "
  //<< n_neighbors_xy[0] <<":"<<n_neighbors_xy[1] <<":"<<n_columns<<endl;
  //cout << "n_cell_pairs : " << n_cell_pairs << " / " << tmp_n_cell_pair<<endl;
  
  return 0;
}
bool MiniCell::check_valid_pair(const int cell1_id, const int cell2_id,
				const bool cell1_odd, const bool cell2_odd){
  //avoid taking both i-j, j-i pairs
  //bool cell2_odd = cell2_id%2!=0;
  if (cell1_odd){
    if ((cell2_id < cell1_id && !cell2_odd ) || 
	(cell2_id > cell1_id && cell2_odd )) return false;
  }else{
    if ((cell2_id < cell1_id && cell2_odd ) || 
	(cell2_id > cell1_id && !cell2_odd )) return false;
  }
  return true;
  //return cell1_id <= cell2_id;
}
int MiniCell::add_cell_pair(const int cell_id1, const int cell_id2, const int image[3]){
  //cout << "addpair " << cell_id1 <<"-"<< cell_id2 <<endl;
  // for debug
  //for(int i=0; i< n_cell_pairs; i++){
  //if(cell_pairs[i].cell_id1 == cell_id1 && cell_pairs[i].cell_id2 == cell_id2)
  //cout <<"DUPLICATION OF CELL PAIR!!!: " << cell_id1 << " " << cell_id2 << endl;
  //}
  if(n_cell_pairs >= max_n_cell_pairs){
    cout << "n_cell_pairs >= max_n_cell_pairs " << n_cell_pairs <<" /  " << max_n_cell_pairs << endl;
    exit(0);
  }
  cell_pairs[n_cell_pairs].cell_id1 = cell_id1;
  cell_pairs[n_cell_pairs].cell_id2 = cell_id2;
  int bit_image = 0;
  if(image[0] == -1)      bit_image = bit_image | 1;
  else if(image[0] == 1)  bit_image = bit_image | 2;
  if(image[1] == -1)      bit_image = bit_image | 4;
  else if(image[1] == 1)  bit_image = bit_image | 8;
  if(image[2] == -1)      bit_image = bit_image | 16;
  else if(image[2] == 1)  bit_image = bit_image | 32;
  cell_pairs[n_cell_pairs].image = bit_image;
  //cout << "cell pair [" << n_cell_pairs << "].image= "<< bit_image << endl;
  //cell_pairs[n_cell_pairs].cx = image[0];
  //cell_pairs[n_cell_pairs].cy = image[1];
  //cell_pairs[n_cell_pairs].cz = image[2];

  set_cell_pair_bitmask(cell_id1, cell_id2, cell_pairs[n_cell_pairs].pair_mask);
  n_cell_pairs++;
  return 0;
}

int MiniCell::set_cell_pair_bitmask(const int cell_id1, const int cell_id2, int* bitmask){
  //initialize bitmask
  for(int i=0; i<N_BITMASK; i++) bitmask[i] = 0;

  int a1 = N_ATOM_CELL*cell_id1;
  for(int a1_cell = 0;
      a1_cell < N_ATOM_CELL;   a1++, a1_cell++){
    bool flg1 = false;
    if(atomids[a1] == -1) flg1 = true;
    int a2 = N_ATOM_CELL*cell_id2;
    //int a2_cell = 0;
    for(int a2_cell = 0; 
	a2_cell < N_ATOM_CELL; a2++, a2_cell++){
      int bit_pos = a2_cell * N_ATOM_CELL + a1_cell;
      int mask_id =  bit_pos / 32;
      int mask_pos =  bit_pos % 32;
      int add_bit = 1 << mask_pos;
      bool flg12 = false;
      if(flg1) flg12 = true;
      else if (atomids[a2] == -1) flg12=true;
      else if (cell_id1 == cell_id2 && a1 >= a2) flg12=true;
      else{
	int tail = atomids[a1] * max_n_nb15off + max_n_nb15off;
	for(int i = atomids[a1] * max_n_nb15off;
	    i < tail && nb15off[i] != -1; i++){
	  if(nb15off[i] == atomids[a2]){ 
	    flg12=true;
	    break;}
	}
      }
      //if(atomids[a1]==1 || atomids[a2]==1) flg=true;
      if(flg12){
	int tmp = bitmask[mask_id];
	bitmask[mask_id] = tmp | add_bit;
      }
    }
  }
  return 0;
}

void MiniCell::get_crd(int atomid_grid, real_pw& x, real_pw& y, real_pw& z){
  x = crd[atomid_grid*3];
  y = crd[atomid_grid*3+1];
  z = crd[atomid_grid*3+2];
}

int MiniCell::get_column_id_from_crd(int x, int y){
  while(x < 0) { x += n_cells_xyz[0];}
  while(x >= n_cells_xyz[0]) { x -= n_cells_xyz[0];}  
  while(y < 0) { y += n_cells_xyz[1];}
  while(y >= n_cells_xyz[1]) { y -= n_cells_xyz[1];}    
  int col = x + y * n_cells_xyz[0];
  if (col >= n_columns) {
    cout << "ERROR!! get_column_id_from_crd " << x << "-" << y << ":" << col<<endl;
  }
  return x + y * n_cells_xyz[0];
}
//int MiniCell::get_column_crd_from_id(int* xy, int col_id){
//  xy[0] = col_id % n_cells_xyz[0];
//  xy[1] = col_id / n_cells_xyz[0];
//  return 0;
//}
 int MiniCell::get_cell_id_from_crd(int x, int y, int z){
   //while(x < 0) { x += n_cells_xyz[0];}
  //while(x >= n_cells_xyz[0]) { x -= n_cells_xyz[0];}  
  //while(y < 0) { y += n_cells_xyz[1];}
  //while(y >= n_cells_xyz[1]) { y -= n_cells_xyz[1];}    
  
  //int n_z = n_cells_z[get_column_id_from_crd(x,y)];
  //while(z < 0) { z += n_z; }
  //while(z >= n_z) { z -= n_z; }
  //cout << "dbg get_cell_id_from_crd " << x <<"-"<<y<<"-"<<z<<"="<<idx_crd_cell[x][y][z] << endl;
  return idx_crd_cell[x][y][z];
}

int MiniCell::set_grid_parameters(const int in_n_atoms,
				  const real in_cutoff_pair,
				  const PBC* in_pbc,
				  const int in_max_n_nb15off,
				  int* in_nb15off){
  // set member variables
  //  n_atoms
  //  cutoff_pair
  //  pbc
  n_atoms = in_n_atoms;
  cutoff_pair = in_cutoff_pair;
  cutoff_pair_half = cutoff_pair * 0.5;
  cutoff_pair2 = cutoff_pair * cutoff_pair;
  pbc = (PBC*)in_pbc;
  nb15off = in_nb15off;

  max_n_nb15off = in_max_n_nb15off;

  n_uni_z = pbc->L[2] / cutoff_pair_half;
  L_z_uni = pbc->L[2] / (float)n_uni_z;
  cout << "uniform grid ... n_uni_z: " << n_uni_z << " L_z_uni: " << L_z_uni << endl;
  return 0;
}

/*
real MiniCell::move_crd_in_cell(const int atomid, const int dim, const real val){
  // atomid : atomid (original)
  // dim : dimension (0-1 : x-y)
  // val : distance to be moved in Angstrome

  crd_in_cell[atomid][dim] += val / L_cell_xy[dim];
  if(val /  L_cell_xy[dim] >= 1.0 ){
    cout << "update_coordinates error: atom:" << atomid << " dim:" << dim << endl;
    cout << val << " - " << val/L_cell_xy[dim] << " - " << L_cell_xy[dim] << endl;
  }
  return crd_in_cell[atomid][dim];
}
*/
/*
int MiniCell::get_ene_forces(real_fc& in_ene_vdw,
			     real_fc& in_ene_ele,
			     real_fc**& in_force){
  int i_grid=0;
  for(int i=0; i < get_n_atom_array(); i++, i_grid+=3){
    if(atomids_buf[i] != -1){
      in_force[atomids_buf[i]][0] += work[i_grid];
      in_force[atomids_buf[i]][1] += work[i_grid+1];
      in_force[atomids_buf[i]][2] += work[i_grid+2];
    }
  }
  in_ene_vdw += energy[0];
  in_ene_ele += energy[1];
  return 0;
}
*/
int MiniCell::add_work(const int atomid_grid,
			const real_fc in_w1, const real_fc in_w2, const real_fc in_w3){
  work[atomid_grid*3] += in_w1;
  work[atomid_grid*3+1] += in_w2;
  work[atomid_grid*3+2] += in_w3;
  return 0;
}


int MiniCell::move_atom(const int& atomid, const int& d, const real& diff){
  ///
  /// Called by SubBox::update_coordinates()
  ///
  crd[atomids_rev[atomid]*3+d] += diff;
  return 0;
}

int MiniCell::set_box_info(int in_n_boxes_xyz[], real in_box_l[]){
  // set 
  //   box_id,      box_crd
  //   n_boxes, n_boxes_xyz
  //   box_l,   box_upper_bound,   box_lower_bound
  //   exbox_l, exbox_upper_bound, exbox_lower_bound
  //
  cout << "MiniCell::set_box_info"<<endl;

  box_id = 0;
#if defined(F_MPI)
  //home_box_id = rank;
#endif
  for(int d=0; d < 3 ; d++){
    n_boxes_xyz[d] = in_n_boxes_xyz[d];
  }
  n_boxes = n_boxes_xyz[0] * n_boxes_xyz[1] * n_boxes_xyz[2];

  get_box_crd_from_id((const int&)box_id, box_crd);
  for(int d=0; d < 3 ; d++){
    box_l[d] = in_box_l[d];

    box_lower_bound[d] = pbc->lower_bound[d] + box_l[d]*box_crd[d];
    box_upper_bound[d] = box_lower_bound[d] + box_l[d];
    if(n_boxes_xyz[d] > 1){
      exbox_lower_bound[d] = box_lower_bound[d] - cutoff_pair_half;
      exbox_upper_bound[d] = box_upper_bound[d] + cutoff_pair_half;
      exbox_l[d] = box_l[d] + cutoff_pair;
    }else{
      exbox_lower_bound[d] = box_lower_bound[d];
      exbox_upper_bound[d] = box_upper_bound[d];
      exbox_l[d] = box_l[d];
    }
  }
  print_box_info();
  return 0;
}
int MiniCell::print_box_info(){
  cout << "print_box_info(): box_id=" << box_id << endl;
  cout << "box_crd : " << box_crd[0] << ", "
       << box_crd[1] << ", " << box_crd[2] << " / " 
       << n_boxes_xyz[0] << ", " << n_boxes_xyz[1]
       << ", " << n_boxes_xyz[2] << endl;
  cout << "box_l : " << box_l[0] << ", " << box_l[1]
       << ", " << box_l[2] << endl;
  cout << "exbox_l : " << exbox_l[0] << ", " << exbox_l[1]
       << ", " << exbox_l[2] << endl;
  cout << "box_upper_bound : " << box_upper_bound[0] << ", " << box_upper_bound[1]
       << ", " << box_upper_bound[2] << endl;
  cout << "exbox_upper_bound : " << exbox_upper_bound[0] << ", " << exbox_upper_bound[1]
       << ", " << exbox_upper_bound[2] << endl;
  cout << "box_lower_bound : " << box_lower_bound[0] << ", " << box_lower_bound[1]
       << ", " << box_lower_bound[2] << endl;
  cout << "exbox_lower_bound : " << exbox_lower_bound[0] << ", " << exbox_lower_bound[1]
       << ", " << exbox_lower_bound[2] << endl;
  return 0;
}

int MiniCell::set_crds_to_homebox(real* in_crd,
				  int* in_atomids,
				  int in_n_atoms_box){
  // set 
  //   n_atoms_box, n_atoms_exbox
  //   a part of crd, atomids_buf
  //   only home atoms
  n_atoms_exbox = in_n_atoms_box;
  n_atoms_box = in_n_atoms_box;
  int idx = 0;
  for(int i = 0; i < n_atoms_exbox; i++, idx+=3){
    crd[idx+0] = in_crd[idx+0];
    crd[idx+1] = in_crd[idx+1];
    crd[idx+2] = in_crd[idx+2];
    atomids_buf[i] = in_atomids[i];
    if(crd[idx+0] <= pbc->lower_bound[0] || crd[idx+0] >= pbc->upper_bound[0] || 
       crd[idx+1] <= pbc->lower_bound[1] || crd[idx+1] >= pbc->upper_bound[1] || 
       crd[idx+2] <= pbc->lower_bound[2] || crd[idx+2] >= pbc->upper_bound[2]){
      cout << "CRD!!! " << crd[idx] << " " << crd[idx+1] << " " << crd[idx+2] <<endl; 
    }
  }
  return 0;
}

int MiniCell::set_max_n_atoms_region(){
  // set
  //   max_n_atoms_region
  //   max_n_atoms_box
  //   max_n_atoms_exbox 
  // requires
  //   n_atoms
  //   n_boxes
  //   box_l
  cout << "MiniCell::set_max_n_atoms_region()"<<endl;
  max_n_atoms_box = 0;
  int tmp_max_n_atoms_box = (n_atoms + n_boxes-1)/n_boxes * COEF_MAX_N_ATOMS_CELL;
  real_pw vol_box = box_l[0] * box_l[1] * box_l[2];
  for(int i=0; i<125; i++){
    //cout << "i " << i<<endl;
    int rx[3];
    get_region_crd_from_id(5, i, rx[0], rx[1], rx[2]);
    int nx[3];
    for(int j=0; j<3; j++){
      if(rx[j]==0)
	nx[j] = box_l[j] - cutoff_pair;
      else
	nx[j] = cutoff_pair*0.5;
    }
    real_pw vol = nx[0] * nx[1] * nx[2];
    max_n_atoms_region[i] = vol/vol_box * tmp_max_n_atoms_box;
    max_n_atoms_exbox += max_n_atoms_region[i];
    if(rx[0] != 2 and rx[1] != 2 and rx[2] != 2 and 
       rx[0] != -2 and rx[1] != -2 and rx[2] != -2)
      max_n_atoms_box += max_n_atoms_region[i];
    //cout << "region_id : " << i << " (" <<get_region_id_from_crd(3, rx[0],rx[1],rx[2])<<") "
    //<< " [" <<rx[0]<<","<<rx[1]<<","<<rx[2]<<"]" <<endl;
      //<< " max_n_cell_reg:"<<max_n_cells_region[i]<< endl;
    
  }
  // -1,-1,-1 => 56 regions
  // -1,-1, 0 => 12 regions
  // -1, 0,-1 => 12 regions
  //  0,-1,-1 => 12 regions
  // -1, 0, 0 => 2
  //  0,-1, 0 => 2
  //  0, 0,-1 => 2
  cout << "max_n_atoms_exbox : " << max_n_atoms_exbox << endl;
  cout << "max_n_atoms_box : " << max_n_atoms_box << endl;
  return 0;
}
int MiniCell::get_region_id_from_crd(int width, int rx, int ry, int rz){
  // for 27  regions: width = 3;
  // for 125 regions: width = 5;
  return (rx+width/2) * width*width + (ry+width/2)*width + rz+width/2;
}
int MiniCell::get_region_crd_from_id(int width, int regid,
				     int& rx, int& ry, int& rz){
  int width2 = width*width;
  int width_d2= width / 2;
  rx = regid/(width2) - width_d2;
  ry = (regid%width2)/width - width_d2;
  rz = regid%width - width_d2;
  return 0;
}
int MiniCell::get_box_id_from_crd(const int box_crd[]){
  return n_boxes_xyz[0]*n_boxes_xyz[1]*box_crd[2] +
    n_boxes_xyz[0]*box_crd[1] + box_crd[0];
}

int MiniCell::get_box_crd_from_id(const int& box_id, 
				  int* box_crd){
  box_crd[0] = box_id%(n_boxes_xyz[0]);
  box_crd[1] = (box_id%(n_boxes_xyz[0]*n_boxes_xyz[1])) / n_boxes_xyz[0];
  box_crd[2] = (box_id/(n_boxes_xyz[0]*n_boxes_xyz[1]));
  return 0;
}

int MiniCell::assign_regions(int from_atomid_g=0, int to_atomid_g=-1){
  /*
  if(to_atomid_g<0) to_atomid_g = n_atoms_exbox;
  for(int atomid_g = from_atomid_g; atomid_g < to_atomid_g; atomid_g++){
    int region[3]; 
    for(int d = 0; d < 3; d++){
      if(crd[atomid_g*3+d] < box_lower_bound[d])
	region[d] = -2;
      else if(crd[atomid_g*3+d] < box_lower_bound[d]+cutoff_pair_half)
	region[d] = -1;
      else if(crd[atomid_g*3+d] < box_upper_bound[d]-cutoff_pair_half)
	region[d] = 0;
      else if(crd[atomid_g*3+d] < box_upper_bound[d])
	region[d] = 1;
      else if(crd[atomid_g*3+d] < box_upper_bound[d]+cutoff_pair_half)
	region[d] = 2;
      else{
	cerr << "ASSIGN REGION ERROR: ATOM " << atomid_g << " " << atomids[atomid_g] << " " << d << endl;
	cerr << crd[atomid_g*3] << ", " << crd[atomid_g*3+1] << ", " << crd[atomid_g*3+2] << endl;
	exit(1);
      }
    }
    int region_id = get_region_id_from_crd(5, region[0], region[1], region[2]);
    add_atom_to_region(atomid_g, region_id);
  }
  */
  return 0;
}

int MiniCell::add_atom_to_region(int atomid_g, int region_id){
  region_atoms[region_id][n_atoms_region[region_id]] = atomid_g;
  n_atoms_region[region_id]++;
  return 0;
}
/*
int MiniCell::copy_pbc_replica(int dir_x, int dir_y, int dir_z){
  // This function is used when there are no space decomposition 
  //   in the axis (n_boxes_xyz[d] == 1)
  // Atoms in the regions with +1(-1) position to the regions
  //   with -2(+2) postion.
  // Only one of dir_x, dir_y, dir_z was -1 or +1.
  for(int region_id=0; region_id < 125; region_id++){
    int rx,ry,rz;
    get_region_crd_from_id(5, region_id,
			   rx, ry, rz);
    if(dir_x != 0 && dir_x != rx) continue;			   
    if(dir_y != 0 && dir_y != ry) continue;
    if(dir_z != 0 && dir_z != rz) continue;
    for(int i=0; i < n_atoms_region[region_id]; i++){
      for(int d=0; d < 3; d++){

      }
    } 
  }
  return 0;
}
*/

int MiniCell::setup_replica_regions(){
  /*
  assign_regions();
  if(n_boxes_xyz[0] == 1){
  }else{
    mpi_send_crd_replica_regions(1,0,0);
    mpi_send_crd_replica_regions(-1,0,0);
  }
  if(n_boxes_xyz[1] == 1){
  }else{
    mpi_send_crd_replica_regions(0,1,0);
    mpi_send_crd_replica_regions(0,-1,0);
  }
  if(n_boxes_xyz[2] == 1){  
  }else{
    mpi_send_crd_replica_regions(0,0,1);
    mpi_send_crd_replica_regions(0,0,-1);
  }
  assign_regions(n_atoms_box, n_atoms_exbox);
  */
  return 0;
}

int MiniCell::mpi_send_crd_replica_regions(int dir_x, int dir_y, int dir_z){
  // args
  //  dir_x, _y, _z : -1, 0, or +1
  //    1, 0, 0 ... send to +1x, receive from -1x
  
  mpi_n_atoms = 0;
  for(int region_id=0; region_id < 125; region_id++){
    int rx,ry,rz;
    get_region_crd_from_id(5, region_id,
			   rx, ry, rz);
    if(dir_x != 0 && dir_x != rx) continue;			   
    if(dir_y != 0 && dir_y != ry) continue;
    if(dir_z != 0 && dir_z != rz) continue;
    for(int i=0; i < n_atoms_region[region_id]; i++){
      for(int d=0; d < 3; d++)
	mpi_sendbuf_crd[mpi_n_atoms*3+d] = crd[region_atoms[region_id][i]*3+d];
      mpi_sendbuf_atomids[mpi_n_atoms] = atomids[region_atoms[region_id][i]];
      mpi_n_atoms++;
    }
  }
#if defined(F_MPI)
  // MPI
  
  
#endif  
  return 0;
}

/*
int MiniCell:update_coordinates(const real** vel_next, const real& time_step){
  for(int atomid_grid=0, i_grid=0;
      atom_id_grid < get_n_atom_array();
      atomid_grid++, i_grid += 3){
    if(atomids[atom_id_grid] < 0) continue;
    for(int d=0; d<3; d++){
      real diff = vel_next[atomids[atom_id_grid]][0] * time_step;
      
    }
  }
  return 0;
  }*/


int MiniCell::get_uni_z(int uni_id){
  return (uni_id / n_columns);
}

int MiniCell::get_uni_id_from_crd(int x, int y, int z){
  return x * n_uni_z + y * (n_uni_z * n_cells_xyz[0]) + z;
}
