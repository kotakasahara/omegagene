#ifndef __ZERO_MULTIPOLE_SUM_H__
#define __ZERO_MULTIPOLE_SUM_H__

#include <iostream>
#include <cmath>
#include "ElectrostaticObject.h"

using namespace std;

class ZeroMultipoleSum : public ElectrostaticObject {
 private:
 protected:

  real_pw cutoff;
  real_pw ewald_alpha;
  real_pw bcoeff;
  real_pw fcoeff;
  real_pw scoeff;
  real_pw zcore;
  real_pw piewald;
  real_pw cutoff_2;
  real_pw cutoff_3;

  real_pw d_self;
  real_pw d_self_mon;

 public:

  ZeroMultipoleSum();
  int (ZeroMultipoleSum::*func_calc_zms)(real_pw&, real_pw&, const real_pw&, const real_pw&, const real_pw&, const real_pw&, const real_pw&, const real_pw&);
  int (ZeroMultipoleSum::*func_calc_zms_excess)(real&, real&, const real&, const real&, const real&, const real&);

  virtual int set_config_parameters(const Config* in_cfg);
  virtual int initial_preprocess();
  int set_zms_params();
  int cal_self_energy(const int& n_atoms,
		      const int& n_excess,
		      const int**& excess_pairs,
		      /*
		      const int& n_bonds,
		      const int**& bond_atomid_pairs,
		      const int& n_angles,
		      const int**& angle_atomid_triads,
		      const int& n_torsions,
		      const int**& torsion_atomid_quads,
		      const int*& torsion_nb14,*/
		      real_pw*& charge,
		      real*& energy_self,
		      real& energy_self_sum);
  inline real_pw get_bcoeff(){ return bcoeff; }
  inline real_pw get_fcoeff(){ return fcoeff; }
  inline real_pw get_scoeff(){ return scoeff; }
  inline real_pw get_zcore(){ return zcore; }
  int calc_zms_alpha0(real_pw& ene_ele, real_pw& grad_coeff,
		      const real_pw& r12,      const real_pw& r12sq,
		      const real_pw& r12inv,   const real_pw& r12inv_2,
		      const real_pw& r12inv_3, const real_pw& cc);
  int calc_zms_excess_alpha0(real& ene_ele, real& grad_coeff,
			     const real& r12,      const real& r12_2,
			     const real& r12_3_inv, const real& cc);
  int calc_zms(real_pw& ene_ele, real_pw& grad_coeff,
	       const real_pw& r12,      const real_pw& r12_2,
	       const real_pw& r12_inv,   const real_pw& r12_2_inv,
	       const real_pw& r12_3_inv, const real_pw& cc);
  int calc_zms_excess(real& ene_ele, real& grad_coeff,
		      const real& r12,      const real& r12_2,
		      const real& r12_3_inv, const real& cc);
  
  
};

#endif
