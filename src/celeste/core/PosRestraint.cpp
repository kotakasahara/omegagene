#include "PosRestraint.h"

using namespace std;

PRUnit::PRUnit(){
}
PRUnit::~PRUnit(){
}
int PRUnit::set_parameters(int in_atomid,
			   real crd_x, real crd_y, real crd_z,
			   real in_dist_margin, real in_coef){
  atomid = in_atomid;
  crd[0] = crd_x;
  crd[1] = crd_y;
  crd[2] = crd_z;
  dist_margin = in_dist_margin;
  coef = in_coef;
  return 0;
}

////////////////////////////////////////////////////////

PosRestraintObject::PosRestraintObject(){
  n_prunits = 0;
  max_n_prunits = 0;
  //boltz = -GAS_CONST*FORCE_VEL;
}
PosRestraintObject::~PosRestraintObject(){
  free_prunits();
}
int PosRestraintObject::alloc_prunits(const int in_n){
  max_n_prunits = in_n;
  cout << "alloc " << max_n_prunits << endl;
  prunits = new PRUnit[max_n_prunits];
  return 0;
}
int PosRestraintObject::free_prunits(){
  delete[] prunits;
  return 0;
}
int PosRestraintObject::add_prunit(int in_aid, real in_x, real in_y, real in_z,
				   real in_margin, real in_coef){
  prunits[n_prunits].set_parameters(in_aid, in_x, in_y, in_z,
				    in_margin, in_coef);
  n_prunits++;
  return n_prunits;
}
real_fc PosRestraintObject::apply_restraint(int n_atoms, real** crd, PBC& pbc,
					    real** force){

  return 0;
}

///////////////////////////////////////////////////////////////

PosRestraintHarmonic::PosRestraintHarmonic()
  : PosRestraintObject(){
}
PosRestraintHarmonic::~PosRestraintHarmonic(){
  free_prunits();
}

real_fc PosRestraintHarmonic::apply_restraint(int n_atoms,
					      real** crd, PBC& pbc,
					      real** force){


  for ( int i = 0; i < n_atoms; i++){
    for (int d=0; d<3; d++){
      force[i][d] = 0.0;
    }
  }
  real_fc ene = 0.0;
  for( int i = 0; i < n_prunits; i++){
    real diff[3];
    pbc.diff_crd_minim_image(diff, crd[prunits[i].get_atomid()],
			     prunits[i].get_crd());
    real r_sq = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
    real r = sqrt(r_sq);

    real coef;
    real ref_dist;
    if(r > prunits[i].get_dist_margin()){
      coef = prunits[i].get_coef();
      ref_dist = prunits[i].get_dist_margin();
    } else continue;
    real k = weight * coef;
    real dist_diff = r - ref_dist;
    ene += k * dist_diff * dist_diff;

    // force
    real k_g = 2.0 * k * dist_diff;
    real_fc frc[3];
    for (int d = 0; d < 3; d++){
      //force[drunits[i].atomid1] += diff[d] / r;
      force[prunits[i].get_atomid()][d] -= k_g * diff[d]/r;
      //frc[d] += diff[d] / r;
    }
  }
  return ene;
}
