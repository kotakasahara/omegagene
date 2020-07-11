#include "PosRestraint.h"

using namespace std;

PRUnit::PRUnit() {}

PRUnit::~PRUnit() {}

int PRUnit::set_parameters(int in_atomid, real crd_x, real crd_y, real crd_z, real in_dist_margin, real in_coef, int in_rest_type, int in_n_params, real* in_params) {
    atomid      = in_atomid;
    crd[0]      = crd_x;
    crd[1]      = crd_y;
    crd[2]      = crd_z;
    dist_margin = in_dist_margin;
    coef        = in_coef;
    rest_type   = in_rest_type;
    n_params    = in_n_params;
    for(int i=0; i<n_params; i++){
      params[i] = in_params[i];
    }
    //params      = in_params;
    //cout << "dbg0708b "  << params[0]<< " " << params[1] << endl;
    return 0;
}

////////////////////////////////////////////////////////

PosRestraintObject::PosRestraintObject() {
    n_prunits     = 0;
    max_n_prunits = 0;
    // boltz = -GAS_CONST*FORCE_VEL;
}

PosRestraintObject::~PosRestraintObject() {
    free_prunits();
}

int PosRestraintObject::alloc_prunits(const int in_n) {
    max_n_prunits = in_n;
    cout << "alloc " << max_n_prunits << endl;
    prunits = new PRUnit[max_n_prunits];
    return 0;
}

int PosRestraintObject::free_prunits() {
    delete[] prunits;
    return 0;
}

int PosRestraintObject::add_prunit(int in_aid, real in_x, real in_y, real in_z, real in_margin, real in_coef, int in_type, int in_n_params, real* in_params) {
  prunits[n_prunits].set_parameters(in_aid, in_x, in_y, in_z, in_margin, in_coef, in_type,
				    in_n_params, in_params);
  n_prunits++;
  return n_prunits;
}

real_fc PosRestraintObject::apply_restraint(int n_atoms, real *crd, PBC &pbc, real **force) {

    return 0;
}

///////////////////////////////////////////////////////////////

PosRestraintHarmonic::PosRestraintHarmonic() : PosRestraintObject() {}

PosRestraintHarmonic::~PosRestraintHarmonic() {
    free_prunits();
}

real_fc PosRestraintHarmonic::apply_restraint(int n_atoms, real *crd, PBC &pbc, real **force) {
  
    for (int i = 0; i < n_atoms; i++) {
      for (int d = 0; d < 3; d++) { force[i][d] = 0.0; }
    }
    real_fc ene = 0.0;
    for (int i = 0; i < n_prunits; i++) {
      real diff[3];
      real t_crd[3] = {crd[prunits[i].get_atomid()*3],
		       crd[prunits[i].get_atomid()*3+1],
		       crd[prunits[i].get_atomid()*3+2]};
      pbc.diff_crd_minim_image(diff, t_crd, prunits[i].get_crd());
      real r_sq = diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
      if(prunits[i].get_rest_type() == POSRESTUNIT_Z) { r_sq = diff[2] * diff[2]; } 
      real r    = sqrt(r_sq);
      real coef;
      real ref_dist;
      if (r > prunits[i].get_dist_margin()) {
	coef     = prunits[i].get_coef();
	ref_dist = prunits[i].get_dist_margin();
      } else
	continue;
      real k         = weight * coef;
      real dist_diff = r - ref_dist;
      if (dist_diff < EPS) continue;
      real dist_diff_sq = dist_diff * dist_diff;
      real c_ene  = k * dist_diff_sq;
      
      if(prunits[i].get_rest_type() == POSRESTUNIT_MULTIWELL01){
	for (int i_param=0; i_param < prunits[i].get_n_params(); i_param++){
	  c_ene *= dist_diff_sq - prunits[i].get_params(i_param);
	  //if(i_param == 1) break;
	}
      }else if(prunits[i].get_rest_type() == POSRESTUNIT_INV_HARMONIC){
	c_ene = k * 1.0/((dist_diff+1.0)*(dist_diff+1.0));
      }
      ene += c_ene;
      // force
      if(prunits[i].get_rest_type() == POSRESTUNIT_NORMAL || prunits[i].get_rest_type() == POSRESTUNIT_Z){
	real    k_g = 2.0 * k * dist_diff;
	real_fc frc[3];
	for (int d = 0; d < 3; d++) {
	  // force[drunits[i].atomid1] += diff[d] / r;
	  real f_d =  k_g * diff[d] / dist_diff;
	  if(prunits[i].get_rest_type() == POSRESTUNIT_Z && d != 2) f_d = 0.0;
	  force[prunits[i].get_atomid()][d] += f_d;
	}
      }else if(prunits[i].get_rest_type() == POSRESTUNIT_MULTIWELL01){
	real dist_diff_3 = dist_diff_sq * dist_diff;
	real dist_diff_5 = dist_diff_sq * dist_diff_3;
	real k_g = 6.0 * k * dist_diff_5 +
	  4.0 * k * dist_diff_3 * (prunits[i].get_params(0) + prunits[i].get_params(1)) +
	  2.0 * k * dist_diff * prunits[i].get_params(0) * prunits[i].get_params(1);
	
	for (int d = 0; d < 3; d++) {
	  real f_d = (k_g * 1.0 / dist_diff * diff[d]);
	  force[prunits[i].get_atomid()][d] += f_d;
	  //cout << "dbg0708a " << dist_diff << " " << k << " " 
	  //<< prunits[i].get_params(0) << " " 
	  //<< prunits[i].get_params(1) << " "  << endl;
	}
      }else if(prunits[i].get_rest_type() == POSRESTUNIT_INV_HARMONIC){
      }
      // Frc[d] += diff[d] / r;
    }

    return ene;
}
///////////////////////////////////////////////////////////////
