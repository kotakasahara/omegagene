#include "Constraint.h"

const int Constraint::Pairs_idx[6][2] = {{0,1},{1,2},{2,0},{3,0},{2,3},{1,3}};

Constraint::Constraint() 
  : CelesteObject() {
  n_pair = 0;
  n_trio = 0;
  n_quad = 0;
  max_n_pair = 0;
  max_n_trio = 0;
  max_n_quad = 0;
  max_loops = 1000;
  tolerance = 1e-5;
}

Constraint::~Constraint(){
  free_constraint();
}

int Constraint::alloc_constraint(){
  if(max_n_pair > 0){
    pair_atomids = new int*[max_n_pair];
    pair_dist    = new real[max_n_pair];
    for(int i=0; i < max_n_pair; i++){
      pair_atomids[i] = new int[2];
      pair_atomids[i][0] = -1;
      pair_atomids[i][1] = -1;
    }
  }
  if(max_n_trio > 0){
    trio_atomids = new int*[max_n_trio];
    trio_dist    = new real*[max_n_trio];
    for(int i=0; i < max_n_trio; i++){
      trio_atomids[i] = new int[3];
      trio_dist[i]    = new real[3];
      trio_atomids[i][0] = -1;
      trio_atomids[i][1] = -1;
      trio_atomids[i][2] = -1;
    }
  }
  if(max_n_quad > 0){
    quad_atomids = new int*[max_n_quad];
    quad_dist    = new real*[max_n_quad];
    for(int i=0; i < max_n_quad; i++){
      quad_atomids[i] = new int[4];
      quad_dist[i]    = new real[6];
      quad_atomids[i][0] = -1;
      quad_atomids[i][1] = -1;
      quad_atomids[i][2] = -1;
      quad_atomids[i][3] = -1;
    }
  }
  return 0;
}

int Constraint::free_constraint(){
  if(max_n_pair > 0){
    for(int i=0; i < max_n_pair; i++){
      delete[] pair_atomids[i];
    }
    delete[] pair_atomids;
    delete[] pair_dist;
  }
  if(max_n_trio > 0){
    for(int i=0; i < max_n_trio; i++){
      delete[] trio_atomids[i];
      delete[] trio_dist[i];
    }
    delete[] trio_atomids;
    delete[] trio_dist;
  }
  if(max_n_quad > 0){
    for(int i=0; i < max_n_quad; i++){
      delete[] quad_atomids[i];
      delete[] quad_dist[i];
    }
    delete[] quad_atomids;
    delete[] quad_dist;
  }
  return 0;
}
int Constraint::set_parameters(int in_max_loops, real in_tolerance){
  max_loops = in_max_loops;
  tolerance = in_tolerance;
  return 0;
}
int Constraint::set_max_n_constraints(int in_n_pair, int in_n_trio, int in_n_quad){
  max_n_pair = in_n_pair;
  max_n_trio = in_n_trio;
  max_n_quad = in_n_quad;
  return 0;
}

int Constraint::add_pair(int atom1, int atom2, real dist1){
  pair_atomids[n_pair][0] = atom1;
  pair_atomids[n_pair][1] = atom2;
  pair_dist[n_pair] = dist1;
  n_pair++;
  return 0;
}

int Constraint::add_trio(int atom1, int atom2, int atom3,
			 real dist1, real dist2, real dist3){
  trio_atomids[n_trio][0] = atom1;
  trio_atomids[n_trio][1] = atom2;
  trio_atomids[n_trio][2] = atom3;
  trio_dist[n_trio][0] = dist1;
  trio_dist[n_trio][1] = dist2;
  trio_dist[n_trio][2] = dist3;
  n_trio++;
  return 0;
}

int Constraint::add_quad(int atom1, int atom2, int atom3, int atom4,
			 real dist1, real dist2, real dist3,
			 real dist4, real dist5, real dist6){
  quad_atomids[n_quad][0] = atom1;
  quad_atomids[n_quad][1] = atom2;
  quad_atomids[n_quad][2] = atom3;
  quad_atomids[n_quad][3] = atom4;
  quad_dist[n_quad][0] = dist1;
  quad_dist[n_quad][1] = dist2;
  quad_dist[n_quad][2] = dist3;
  quad_dist[n_quad][3] = dist4;
  quad_dist[n_quad][4] = dist5;
  quad_dist[n_quad][5] = dist6;
  n_quad++;
  return 0;
}

int Constraint::apply_constraint(real* in_crd, real* in_crd_prev, real* mass,
				 PBC* pbc){
  return 0;
}
int Constraint::calc_linear_eq(real a[6][6],
			       real x[6],
			       real b[6],
			       int size){
  int* pivot = new int[size];
  for( int i=0; i < size; i++){
    pivot[i] = i;
  }
  int index;
  for( int i=0; i < size; i++){
    index = i;
    real max_a = fabs(a[i][pivot[i]]);
    for( int j=0; j < size; j++){
      real tmp_a = fabs(a[i][pivot[j]]);
	if(  tmp_a > max_a ){
	  index = j;
	  max_a = tmp_a;
	}
    }

    if (index != i){
      int tmp = pivot[i];
      pivot[i] = pivot[index];
      pivot[index] = tmp;
    }
    
    if (fabs(a[i][pivot[i]]) < EPS){
      cerr << "matrix is singular" << endl;
      return 1;
    }
    a[i][pivot[i]] = 1.0 / a[i][pivot[i]];
    for(int j = i+1; j < size; j++){
      a[i][pivot[j]] *= a[i][pivot[i]];
      for(int k = i+1; k < size; k++){
	a[k][pivot[j]] -= a[i][pivot[j]] * a[k][pivot[i]];
      }
    }
  }
  // 3 forward substitution
  for( int i = 0; i < size; i++){
    real new_x = b[pivot[i]];
    for( int j = 0; j < i; j++){
      new_x -= a[j][pivot[i]] * x[j];
    } 
    x[i] = new_x;
  }
  // 4 backward substitution
  for( int i = size-1; i >= 0; i--){
    real new_x = x[i];
    for( int j = i+1; j < size; j++){
      new_x -= a[j][pivot[i]] * x[j];
    } 
    x[i] = new_x * a[i][pivot[i]];
  }

  delete[] pivot;
  return 0;
}

int Constraint::set_subset_constraint(Constraint& super,
				      int* atomids_rev){
  n_pair = 0;
  n_trio = 0;
  n_quad = 0;
  //cout << "subset pair" << endl;
  for(int i = 0; i < super.get_n_pair(); i++){
    int a1 = atomids_rev[super.get_pair_atomids()[i][0]];
    int a2 = atomids_rev[super.get_pair_atomids()[i][1]];
    if(a1 != -1){
      add_pair(a1, a2, super.get_pair_dist()[i]);
    }
  }
  //cout << "subset trio" << endl;
  for(int i = 0; i < super.get_n_trio(); i++){
    int a1 = atomids_rev[super.get_trio_atomids()[i][0]];
    int a2 = atomids_rev[super.get_trio_atomids()[i][1]];
    int a3 = atomids_rev[super.get_trio_atomids()[i][2]];
    if(a1 != -1){
      add_trio(a1, a2, a3,
	       super.get_trio_dist()[i][0],
	       super.get_trio_dist()[i][1],
	       super.get_trio_dist()[i][2]);
    }
  }
  //cout << "subset quad" << endl;
  for(int i = 0; i < super.get_n_quad(); i++){
    int a1 = atomids_rev[super.get_quad_atomids()[i][0]];
    int a2 = atomids_rev[super.get_quad_atomids()[i][1]];
    int a3 = atomids_rev[super.get_quad_atomids()[i][2]];
    int a4 = atomids_rev[super.get_quad_atomids()[i][3]];
    if(a1 != -1){
      add_quad(a1, a2, a3, a4,
	       super.get_quad_dist()[i][0],
	       super.get_quad_dist()[i][1],
	       super.get_quad_dist()[i][2],
	       super.get_quad_dist()[i][3],
	       super.get_quad_dist()[i][4],
	       super.get_quad_dist()[i][5]);
    }
  }
  //cout << "subset end"<< endl;
  return 0;
}
