#include "ZeroMultipoleSum.h"
#include <iomanip>

using namespace std;

ZeroMultipoleSum::ZeroMultipoleSum(const int in_zms_mode, const real in_alpha,
				   const real in_cutoff,
				   const ZmsCalc in_func, const ZmsCalcExcess in_func_excess)
  //int (ZeroMultipoleSum::*func_excess)(real&, real&, const real&, const real&, const real&, const real&, const real&, const real&))
  : zms_mode(in_zms_mode), ewald_alpha(in_alpha), cutoff(in_cutoff),
    func_calc_zms(in_func), func_calc_zms_excess(in_func_excess),
    ElectrostaticObject(){
  if (DBG>=1)
    cout << "DBG1: ZeroMultipoleSum::ZeroMultipoleSum()" << endl;
  cout << "cutoff : " << cutoff << endl;
  //cutoff = in_cutoff;
  /*
  if (ewald_alpha < EPS){
    switch(zms_mode){
    case ELCTRST_ZERODIPOLE:
      func_calc_zms = &ZeroMultipoleSum::calc_zero02pole_alpha0;
      func_calc_zms_excess = &ZeroMultipoleSum::calc_zero02pole_excess_alpha0;
      break;
    case ELCTRST_ZEROQUADRUPOLE:
      func_calc_zms = &ZeroMultipoleSum::calc_zero04pole_alpha0;
      func_calc_zms_excess = &ZeroMultipoleSum::calc_zero04pole_excess_alpha0;
      break;
    case ELCTRST_ZEROOCTUPOLE:
      func_calc_zms = &ZeroMultipoleSum::calc_zero08pole_alpha0;
      func_calc_zms_excess = &ZeroMultipoleSum::calc_zero08pole_excess_alpha0;
      break;
    case ELCTRST_ZEROHEXADECAPOLE:
      func_calc_zms = &ZeroMultipoleSum::calc_zero16pole_alpha0;
      func_calc_zms_excess = &ZeroMultipoleSum::calc_zero16pole_excess_alpha0;
      break;
    }
  }else{
    switch(zms_mode){
    case ELCTRST_ZERODIPOLE:
      func_calc_zms = &ZeroMultipoleSum::calc_zero02pole;
      func_calc_zms_excess = &ZeroMultipoleSum::calc_zero02pole_excess;
      break;
    case ELCTRST_ZEROQUADRUPOLE:
      func_calc_zms = &ZeroMultipoleSum::calc_zero04pole;
      func_calc_zms_excess = &ZeroMultipoleSum::calc_zero04pole_excess;
      break;
    case ELCTRST_ZEROOCTUPOLE:
      func_calc_zms = &ZeroMultipoleSum::calc_zero08pole;
      func_calc_zms_excess = &ZeroMultipoleSum::calc_zero08pole_excess;
      break;
    case ELCTRST_ZEROHEXADECAPOLE:
      func_calc_zms = &ZeroMultipoleSum::calc_zero16pole;
      func_calc_zms_excess = &ZeroMultipoleSum::calc_zero16pole_excess;
      break;
    }
  }
  */
}

int ZeroMultipoleSum::set_config_parameters(const Config* cfg){
  if (DBG>=1)
    cout << "DBG1: ZeroMultipoleSum::set_config_parameters()" << endl;
  //ElectrostaticObject::set_config_parameters(cfg);
  return 0;
}
int ZeroMultipoleSum::initial_preprocess(){
  // self energy
  //zms_mode = 1;
  set_zms_params();
  return 0;
}

int ZeroMultipoleSum::cal_self_energy(const int& n_atoms,
				      const int& n_excess,
				      const int**& excess_pairs,
				      real_pw*& charge,
				      real*& energy_self,
				      real& energy_self_sum){

				      /*int ZeroMultipoleSum::cal_self_energy(const int& n_atoms,
				      const int& n_bonds,
				      const int**& bond_atomid_pairs,
				      const int& n_angles,
				      const int**& angle_atomid_triads,
				      const int& n_torsions,
				      const int**& torsion_atomid_quads,
				      const int*& torsion_nb14,
				      const real_pw*& charge,
				      real_pw*& energy_self,
				      real_pw& energy_self_sum){*/

  if (DBG >= 1)
    cout << "DBG1: ZeroMultipoleSum::cal_self_energy()"<<endl;
  real* dwork;
  dwork = new real[n_atoms];
  for (int i=0; i < n_atoms; i++)
    dwork[i] = 0.0;
  for (int i=0; i < n_excess; i++){
    dwork[excess_pairs[i][0]] += charge[excess_pairs[i][1]];
    dwork[excess_pairs[i][1]] += charge[excess_pairs[i][0]];
  }
  real_pw chgsum = 0.0;
  real_pw chgsmm = 0.0;
  d_self_mon = 0.0;
  real zmcore = 0.0;
  switch(zms_mode){
  case ELCTRST_ZERODIPOLE:       zmcore = zcore;  break;
  case ELCTRST_ZEROQUADRUPOLE:   zmcore = zqcore; break;
  case ELCTRST_ZEROOCTUPOLE:     zmcore = zocore; break;
  case ELCTRST_ZEROHEXADECAPOLE: zmcore = zhcore; break;
  default: break;
  }
  for (int i=0; i < n_atoms; i++){
    energy_self[i] = 0.0;
    energy_self[i] =
      ( -piewald * charge[i] * charge[i]
	- zmcore * charge[i] * (charge[i]+dwork[i]) )
      * 0.5 * CHARGE_COEFF;
    chgsum += charge[i] * charge[i];
    chgsmm += charge[i] * (charge[i] + dwork[i]);
    d_self_mon += energy_self[i];
  }
  d_self = (-zmcore * chgsmm - piewald * chgsum) * 0.5 * CHARGE_COEFF;

  energy_self_sum = d_self_mon;
  if(DBG>=1){
    cout << setprecision(15);
    cout << "fcoeff : "  << fcoeff << endl;
    cout << "bcoeff : "  << bcoeff << endl;
    cout << "scoeff : "  << scoeff << endl;
    cout << "zmcore : "  << zmcore << endl;
    cout << "chgsum : "  << chgsum << endl;
    cout << "chgsmm: "  << chgsmm << endl;
    cout << "selfEnergy_dcore : "  << zmcore*chgsmm*0.5*CHARGE_COEFF << endl;
    cout << "d_self: "  << d_self << endl;
    cout << "sum of energy self: "  << d_self_mon << endl;
  }
  return 0;
}

int ZeroMultipoleSum::set_zms_params(){
  if (DBG >= 1)
    cout << "DBG1: ZeroMultipoleSum::set_zms_params()"<<endl;

  bcoeff = 0.0;
  fcoeff = 0.0;
  zcore = 0.0;
  piewald = 0.0;
  scoeff = 0.0;
  cutoff_2 = cutoff * cutoff;
  cutoff_3 = cutoff_2 * cutoff;
  cutoff_4 = cutoff_2 * cutoff_2;
  cutoff_5 = cutoff_3 * cutoff_2;
  cutoff_6 = cutoff_3 * cutoff_3;
  cutoff_7 = cutoff_3 * cutoff_4;
  real_pw tmp = cutoff * ewald_alpha;
  real_pw errorfc = erfc(tmp);
  real_pw expterm = exp(-tmp*tmp);

  if (ewald_alpha < EPS){
    fcoeff = 1.0 / cutoff_3;
    bcoeff = 0.5 * fcoeff;
    scoeff = 1.5 / cutoff_3;
    hzqrc = 2.0 / cutoff_3;
    coeffd1 = 1.0 / cutoff_2;
    coeffd2 = 2.0 / cutoff_3;
    coeffd3 = 6.0 / cutoff_4;
    coeffd4 = 24.0 / cutoff_5;
    //zcore = 1.5 / cutoff;
  }else{
    real_pw root_pi = sqrt(PI);

    piewald = 2.0 / root_pi * ewald_alpha;

    fcoeff = errorfc / cutoff_3 + piewald * expterm / cutoff_2;
    bcoeff = 0.5 * fcoeff;
    scoeff = bcoeff + errorfc / cutoff_3;
    real_pw tmp_3 = tmp * tmp * tmp;
    real_pw tmp_5 = tmp_3 * tmp * tmp;
    real_pw tmp_7 = tmp_5 * tmp * tmp;

    hzqrc = ( errorfc + expterm * (tmp + tmp_3) * 2.0 / root_pi ) * 2.0 / cutoff_3;
    coeffd1 = (errorfc + expterm * tmp * 2.0 / root_pi ) / cutoff_2;
    coeffd2 = (errorfc + expterm * (tmp + tmp_3) * 2.0 / root_pi) * 2.0 / cutoff_3;
    coeffd3 = (errorfc + expterm * (tmp + 2.0/3.0 * tmp_3 + 2.0/3.0 * tmp_5 )
	       * 2.0/root_pi) * 6.0/cutoff_4;
    coeffd4 = (errorfc + expterm * (tmp + 2.0/3.0 * tmp_3 + 1.0/6.0 * tmp_5
	       + 1.0/3.0 + tmp_7 ) * 2.0/root_pi) * 24.0/cutoff_5;

  }

  zqcoeff2 = -3.0 / 4.0 * fcoeff - hzqrc / 4.0;
  zqcoeff4 = fcoeff / 8.0 / cutoff_2 + hzqrc / 8.0 / cutoff_2;;
  zqcoeff22 = 2.0 * zqcoeff2;
  zqcoeff44 = 4.0 * zqcoeff4;

  zocoeff2 = - 15.0/16.0 * coeffd1 / cutoff
    - 7.0/16.0 * coeffd2 - 1.0/16.0 * coeffd3 * cutoff;
  zocoeff4 = 5.0/16.0 * coeffd1 / cutoff_3
    + 5.0/16.0 * coeffd2 / cutoff_2
    + 1.0/16.0 * coeffd3 / cutoff;
  zocoeff6 = -1.0/16.0 * coeffd1 / cutoff_5
    - 1.0/16.0 * coeffd2 / cutoff_4
    - 1.0/48.0 * coeffd3 / cutoff_3;
  zocoeff22 = 2.0 * zocoeff2;
  zocoeff44 = 4.0 * zocoeff4;
  zocoeff66 = 6.0 * zocoeff6;

  zhcoeff2 = -35.0/32.0 * coeffd1 / cutoff
    - 19.0/32.0 * coeffd2
    - 1.0/96.0 * coeffd4 * cutoff_2
    - 1.0/8.0 * coeffd3 * cutoff;
  zhcoeff4 = 1.0/64.0 * coeffd4
    + 35.0/64.0 * coeffd1 / cutoff_3
    + 35.0/64.0 * coeffd2 / cutoff_2
    + 5.0/32.0 * coeffd3 / cutoff;
  zhcoeff6 = -7.0/32.0 * coeffd1 / cutoff_5
    - 7.0/32.0 * coeffd2 / cutoff_4
    - 1.0/12.0 * coeffd3 / cutoff_3
    - 1.0/96.0 * coeffd4 / cutoff_2;
  zhcoeff8 = 5.0/128.0 * coeffd2 / cutoff_6
    + 5.0/128.0 * coeffd1 / cutoff_7
    + 1.0/64.0 * coeffd3 / cutoff_5
    + 1.0/384.0 * coeffd4 / cutoff_4;
  zhcoeff22 = 2.0 * zhcoeff2;
  zhcoeff44 = 4.0 * zhcoeff4;
  zhcoeff66 = 6.0 * zhcoeff6;
  zhcoeff88 = 8.0 * zhcoeff8;

  zcore = errorfc / cutoff + bcoeff * cutoff_2;
  zqcore = errorfc / cutoff - zqcoeff2 * cutoff_2
    - zqcoeff4 * cutoff_4;
  zocore = errorfc / cutoff - zocoeff2 * cutoff_2
    - zocoeff4 * cutoff_4
    - zocoeff6 * cutoff_6;
  zhcore = errorfc / cutoff - zhcoeff2 * cutoff_2
    - zhcoeff4 * cutoff_4
    - zhcoeff6 * cutoff_6
    - zhcoeff8 * cutoff_6 * cutoff_2;

  return 0;
}
int ZeroMultipoleSum::calc_zero02pole_excess_alpha0(real& ene_ele, real& grad_coeff,
						    const real& r12,      const real& r12_2,
						    const real& r12_inv,
						    const real& r12_2_inv, const real& r12_3_inv,
						    const real& cc){
  ene_ele = cc * bcoeff * r12_2;
  grad_coeff = cc * fcoeff;
  return 0;
}
int ZeroMultipoleSum::calc_zero02pole_excess(real& ene_ele, real& grad_coeff,
					     const real& r12,       const real& r12_2,
					     const real& r12_inv,
					     const real& r12_2_inv, const real& r12_3_inv,
					     const real& cc){
  real_pw tmp = ewald_alpha * r12;
  real_pw errorfc = erfc(tmp);
  real_pw errorfn = 1.0 - errorfc;
  ene_ele = cc * (-errorfn * r12 + bcoeff * r12_2);
  grad_coeff = -cc * ( errorfc * r12_3_inv + piewald * exp(-tmp*tmp)
		       * r12_2_inv - fcoeff - r12_3_inv);
  return 0.0;
}
int ZeroMultipoleSum::calc_zero04pole_excess_alpha0(real& ene_ele, real& grad_coeff,
						    const real& r12,      const real& r12_2,
						    const real& r12_inv,
						    const real& r12_2_inv, const real& r12_3_inv,
						    const real& cc){
  ene_ele = cc * (-zqcoeff2 * r12_2 - zqcoeff4 * r12_2 * r12_2);
  grad_coeff = -cc * (zqcoeff22 + zqcoeff44 * r12_2);
  return 0;
}
int ZeroMultipoleSum::calc_zero04pole_excess(real& ene_ele, real& grad_coeff,
					     const real& r12,      const real& r12_2,
					     const real& r12_inv,
					     const real& r12_2_inv, const real& r12_3_inv,
					     const real& cc){
  real_pw tmp = ewald_alpha * r12;
  real_pw errorfc = erfc(tmp);
  real_pw errorfn = 1.0 - errorfc;
  ene_ele = cc * (-errorfn * r12_inv - zqcoeff2 * r12_2 - zqcoeff4 * r12_2 * r12_2);
  grad_coeff = -cc * ( errorfc * r12_3_inv + piewald * exp(-tmp*tmp) * r12_2_inv
		       + zqcoeff22 + zqcoeff44 * r12_2 - r12_3_inv);
  return 0.0;
}
int ZeroMultipoleSum::calc_zero08pole_excess_alpha0(real& ene_ele, real& grad_coeff,
						    const real& r12,      const real& r12_2,
						    const real& r12_inv,
						    const real& r12_2_inv, const real& r12_3_inv,
						    const real& cc){
  const real r12_4 = r12_2 * r12_2;
  ene_ele = cc * (-zocoeff2 * r12_2 - zocoeff4 * r12_4 - zocoeff6 * r12_4 * r12);
  grad_coeff = -cc * (zocoeff22 + zocoeff44 * r12_2 + zocoeff66 * r12_4);
  return 0;
}
int ZeroMultipoleSum::calc_zero08pole_excess(real& ene_ele, real& grad_coeff,
					     const real& r12,      const real& r12_2,
					     const real& r12_inv,
					     const real& r12_2_inv, const real& r12_3_inv,
					     const real& cc){
  const real_pw tmp = ewald_alpha * r12;
  const real_pw errorfc = erfc(tmp);
  const real_pw errorfn = 1.0 - errorfc;
  const real_pw r12_4 = r12_2 * r12_2;
  ene_ele = cc * (-errorfn * r12_inv - zocoeff2 * r12_2 - zocoeff4 * r12_4
		  - zocoeff6 * r12_4 * r12);
  grad_coeff = -cc * ( errorfc * r12_3_inv + piewald * exp(-tmp*tmp) * r12_2_inv
		       + zocoeff22 + zocoeff44 * r12_2
		       + zocoeff66 * r12_4 - r12_3_inv);
  return 0.0;
}
int ZeroMultipoleSum::calc_zero16pole_excess_alpha0(real& ene_ele, real& grad_coeff,
						    const real& r12,      const real& r12_2,
						    const real& r12_inv,
						    const real& r12_2_inv, const real& r12_3_inv,
						    const real& cc){
  const real_pw r12_4 = r12_2 * r12_2;
  const real_pw r12_6 = r12_4 * r12_2;
  const real_pw r12_8 = r12_4 * r12_4;
  ene_ele = cc * ( -zhcoeff2 * r12_2 - zhcoeff4 * r12_4
		   - zhcoeff6 * r12_6 - zhcoeff8 * r12_8);
  grad_coeff = -cc * (zhcoeff22 + zhcoeff44 * r12_2
		      + zhcoeff66 * r12_2 + zhcoeff88 * r12_6 - r12_3_inv);
  return 0;
}
int ZeroMultipoleSum::calc_zero16pole_excess(real& ene_ele, real& grad_coeff,
					     const real& r12,      const real& r12_2,
					     const real& r12_inv,
					     const real& r12_2_inv, const real& r12_3_inv,
					     const real& cc){
  const real_pw r12_4 = r12_2 * r12_2;
  const real_pw r12_6 = r12_4 * r12_2;
  const real_pw tmp = ewald_alpha * r12;
  const real_pw errorfc = erfc(tmp);
  const real_pw errorfn = 1.0 - erfc(tmp);
  ene_ele = cc * (-errorfn * r12_inv - zhcoeff2 * r12_2 - zhcoeff4 * r12_4);
  grad_coeff = -cc * ( errorfc * r12_3_inv + piewald * exp(-tmp*tmp) * r12_2_inv
		       + zhcoeff22 + zhcoeff44 * r12_2
		       + zhcoeff66 * r12_4 + zhcoeff88 * r12_6 - r12_3_inv);
  return 0.0;
}
int ZeroMultipoleSum::calc_zero02pole_alpha0(real_pw& ene_ele, real_pw& grad_coeff,
					     const real_pw& r12,      const real_pw& r12_2,
					     const real_pw& r12_inv,   const real_pw& r12_2_inv,
					     const real_pw& r12_3_inv, const real_pw& cc){
  ene_ele = cc * (r12_inv - zcore + bcoeff * r12_2);
  grad_coeff = -cc * (r12_3_inv - fcoeff);
  /*
  printf("dbgpair  %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
	 //cc, r12_inv, zcore, bcoeff, r12_2,
	 ene_ele,
	 bcoeff, r12_2, tmp1, tmp2, tmp3,
	 (r12_inv + (bcoeff * r12_2))-zcore,
	 r12_inv + ((bcoeff * r12_2)-zcore),
	 r12_inv+(-zcore+( bcoeff*r12_2)),
	 (r12_inv-zcore)+( bcoeff*r12_2),
	 (( bcoeff*r12_2)+r12_inv)-zcore,
	 ( bcoeff*r12_2)+(r12_inv-zcore));
  */
  return 0;
}
int ZeroMultipoleSum::calc_zero02pole(real_pw& ene_ele, real_pw& grad_coeff,
				      const real_pw& r12,      const real_pw& r12_2,
				      const real_pw& r12_inv,   const real_pw& r12_2_inv,
				      const real_pw& r12_3_inv, const real_pw& cc){
  real_pw tmp = r12 * ewald_alpha;
  real_pw errorfc = erfc(tmp);
  ene_ele = cc * (r12_inv * errorfc - zcore + bcoeff * r12_2);
  grad_coeff = -cc * (errorfc * r12_3_inv + piewald * exp(-tmp*tmp) * r12_2_inv - fcoeff);
  return 0;
}
int ZeroMultipoleSum::calc_zero04pole_alpha0(real_pw& ene_ele, real_pw& grad_coeff,
					     const real_pw& r12,      const real_pw& r12_2,
					     const real_pw& r12_inv,   const real_pw& r12_2_inv,
					     const real_pw& r12_3_inv, const real_pw& cc){
  const real_pw r12_4 = r12_2 * r12_2;
  ene_ele    = cc * (r12_inv - zqcore - zqcoeff2 * r12_2 - zqcoeff4 * r12_4);
  grad_coeff = -cc * (r12_3_inv + zqcoeff22 + zqcoeff44 * r12_2);
  //cout << "kkdbg20151209 " << ene_ele << " " << cc << " "
  //<< r12_2 << " " << zqcore << " " << zqcoeff2 << " " << zqcoeff4
  //<< " " << zqcoeff22 << " " << zqcoeff44 << endl;
  return 0;
}
int ZeroMultipoleSum::calc_zero04pole(real_pw& ene_ele, real_pw& grad_coeff,
				      const real_pw& r12,      const real_pw& r12_2,
				      const real_pw& r12_inv,   const real_pw& r12_2_inv,
				      const real_pw& r12_3_inv, const real_pw& cc){
  const real_pw r12_4 = r12_2 * r12_2;
  real_pw tmp = r12 * ewald_alpha;
  real_pw errorfc = erfc(tmp);
  ene_ele = cc * (r12_inv * errorfc - zqcore - zqcoeff2 * r12_2 - zqcoeff4 * r12_4);
  grad_coeff = -cc * (r12_3_inv * errorfc + piewald * exp(-tmp*tmp) * r12_2_inv
		      + zqcoeff22 + zqcoeff44 * r12_2);
  return 0;
}
int ZeroMultipoleSum::calc_zero08pole_alpha0(real_pw& ene_ele, real_pw& grad_coeff,
					     const real_pw& r12,      const real_pw& r12_2,
					     const real_pw& r12_inv,   const real_pw& r12_2_inv,
					     const real_pw& r12_3_inv, const real_pw& cc){
  const real_pw r12_4 = r12_2 * r12_2;
  const real_pw r12_6 = r12_4 * r12_2;
  ene_ele    = cc * (r12_inv - zocore - zocoeff2 * r12_2 - zocoeff4 * r12_4
		     - zocoeff6 * r12_6);
  grad_coeff = -cc * (r12_3_inv + zocoeff22 + zocoeff44 * r12_2
		      + zocoeff66 * r12_4);
  return 0;
}
int ZeroMultipoleSum::calc_zero08pole(real_pw& ene_ele, real_pw& grad_coeff,
				      const real_pw& r12,      const real_pw& r12_2,
				      const real_pw& r12_inv,   const real_pw& r12_2_inv,
				      const real_pw& r12_3_inv, const real_pw& cc){
  const real_pw r12_4 = r12_2 * r12_2;
  const real_pw r12_6 = r12_4 * r12_2;
  real_pw tmp = r12 * ewald_alpha;
  real_pw errorfc = erfc(tmp);
  ene_ele = cc * (r12_inv * errorfc - zocore - zocoeff2 * r12_2
		  - zocoeff4 * r12_4 - zocoeff6 * r12_6);
  grad_coeff = -cc * (r12_3_inv * errorfc + piewald * exp(-tmp*tmp) * r12_2_inv
		      + zhcoeff22 + zhcoeff44 * r12_2
		      + zhcoeff66 * r12_4);
  return 0;
}
int ZeroMultipoleSum::calc_zero16pole_alpha0(real_pw& ene_ele, real_pw& grad_coeff,
					     const real_pw& r12,      const real_pw& r12_2,
					     const real_pw& r12_inv,   const real_pw& r12_2_inv,
					     const real_pw& r12_3_inv, const real_pw& cc){
  const real_pw r12_4 = r12_2 * r12_2;
  const real_pw r12_6 = r12_4 * r12_2;
  const real_pw r12_8 = r12_4 * r12_4;
  ene_ele    = cc * (r12_inv - zhcore - zhcoeff2 * r12_2 - zhcoeff4 * r12_4
		     - zhcoeff6 * r12_6 - zhcoeff8 * r12_8);
  grad_coeff = -cc * (r12_3_inv + zhcoeff22 + zhcoeff44 * r12_2
		      + zhcoeff66 * r12_4 + zhcoeff88 * r12_6);
  return 0;
}
int ZeroMultipoleSum::calc_zero16pole(real_pw& ene_ele, real_pw& grad_coeff,
				      const real_pw& r12,      const real_pw& r12_2,
				      const real_pw& r12_inv,   const real_pw& r12_2_inv,
				      const real_pw& r12_3_inv, const real_pw& cc){
  const real_pw r12_4 = r12_2 * r12_2;
  const real_pw r12_6 = r12_4 * r12_2;
  const real_pw r12_8 = r12_4 * r12_4;
  const real_pw tmp = r12 * ewald_alpha;
  const real_pw errorfc = erfc(tmp);
  ene_ele = cc * (r12_inv * errorfc - zhcore - zhcoeff2 * r12_2
		  - zhcoeff4 * r12_4 - zhcoeff6 * r12_6 - zhcoeff8 * r12_8);
  grad_coeff = -cc * (r12_3_inv * errorfc + piewald * exp(-tmp*tmp) * r12_2_inv
		      + zhcoeff22 + zhcoeff44 * r12_2
		      + zhcoeff66 * r12_4 + zhcoeff88 * r12_6);
  return 0;
}
