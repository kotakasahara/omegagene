#include "ZeroMultipoleSum.h"

ZeroMultipoleSum::ZeroMultipoleSum()
  : ElectrostaticObject(){
  if (DBG>=1)
    cout << "DBG1: ZeroMultipoleSum::ZeroMultipoleSum()" << endl;
}

int ZeroMultipoleSum::set_config_parameters(const Config* cfg){
  if (DBG>=1)
    cout << "DBG1: ZeroMultipoleSum::set_config_parameters()" << endl;
  ElectrostaticObject::set_config_parameters(cfg);
  cutoff = cfg->cutoff;
  ewald_alpha = cfg->ele_alpha;
  return 0;
}
int ZeroMultipoleSum::initial_preprocess(){
  // self energy
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
  real_pw* dwork;
  dwork = new real_pw[n_atoms];
  for (int i=0; i < n_atoms; i++)
    dwork[i] = 0.0;
  for (int i=0; i < n_excess; i++){
    dwork[excess_pairs[i][0]] += charge[excess_pairs[i][1]];
    dwork[excess_pairs[i][1]] += charge[excess_pairs[i][0]];
  }
  /*
  for (int i=0; i < n_bonds; i++){
    const int* pair = bond_atomid_pairs[i];
    dwork[pair[0]] += charge[pair[1]];
    dwork[pair[1]] += charge[pair[0]];
    //cout << "bond : " << i << " , " << pair[0] << " - " << pair[1] << endl;
    //cout << dwork[pair[0]] << " , " << dwork[pair[1]] << endl;
  }
  for (int i=0; i < n_angles; i++){  
    const int* trio = angle_atomid_triads[i];
    dwork[trio[0]] += charge[trio[2]];
    dwork[trio[2]] += charge[trio[0]];
    //cout << "angle : " << i << " , " << trio[0] << " - " << trio[2] << endl;
    //cout << dwork[trio[0]] << " , " << dwork[trio[2]] << endl;
  }
  for (int i=0; i < n_torsions; i++){  	 
    if(torsion_nb14[i] == 1){
      const int* quad = torsion_atomid_quads[i];
      dwork[quad[0]] += charge[quad[3]];
      dwork[quad[3]] += charge[quad[0]];
      //cout << "trosion : " << i << " , " << quad[0] << " - " << quad[3] << endl;
      //cout << dwork[quad[0]] << " , " << dwork[quad[3]] << endl;
    }
  }
  */
  /*
    for (int i=0; i < mmsys->n_nb14; i++){
    int* pair = mmsys->nb14_atomid_pairs[i];
    dwork[pair[0]] += mmsys->charge[pair[1]];
    dwork[pair[1]] += mmsys->charge[pair[0]];
    }
  */
  real_pw chgsum = 0.0;  
  real_pw chgsmm = 0.0;
  //  cout << "DWORK" << endl;
  d_self_mon = 0.0;
  for (int i=0; i < n_atoms; i++){  	 
    //cout << "loop " << i <<endl;
    //cout << "loop " << mmsys->charge[i] <<endl;
    energy_self[i] = 0.0;
    //cout << "dwork " << dwork[i] <<endl;
    energy_self[i] = 
      ( -piewald * charge[i] * charge[i]
	- zcore * charge[i] * (charge[i]+dwork[i]) ) 
      * 0.5 * CHARGE_COEFF;
    //cout << "test1 " <<endl;
    chgsum += charge[i] * charge[i];
    //cout << "test2 " <<endl;
    chgsmm += charge[i] * (charge[i] + dwork[i]);
    //cout << "test3 " <<endl;
    d_self_mon += energy_self[i];
    //cout << "i: " << i << endl;
  }
  d_self = (-zcore * chgsmm - piewald * chgsum) * 0.5 * CHARGE_COEFF;
  energy_self_sum = d_self_mon;
  if(DBG>=1){
    cout << setprecision(15);
    cout << "fcoeff : "  << fcoeff << endl;
    cout << "bcoeff : "  << bcoeff << endl;
    cout << "scoeff : "  << scoeff << endl;
    cout << "zcore : "  << zcore << endl;
    cout << "chgsum : "  << chgsum << endl;
    cout << "chgsmm: "  << chgsmm << endl;
    cout << "selfEnergy_dcore : "  << zcore*chgsmm*0.5*CHARGE_COEFF << endl;
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
  real_pw tmp = cutoff * ewald_alpha;
  real_pw errorfc = erfc(tmp);
  real_pw expterm = exp(-tmp*tmp);

  if (ewald_alpha < EPS){
    func_calc_zms = &ZeroMultipoleSum::calc_zms_alpha0;
    func_calc_zms_excess = &ZeroMultipoleSum::calc_zms_excess_alpha0;
    fcoeff = 1.0 / cutoff_3;
    bcoeff = 0.5 * fcoeff;
    scoeff = 1.5 / cutoff_3;
    zcore = 1.5 / cutoff;
  }else{
    func_calc_zms = &ZeroMultipoleSum::calc_zms;
    func_calc_zms_excess = &ZeroMultipoleSum::calc_zms_excess;
    fcoeff = errorfc / cutoff_3 + piewald * expterm / cutoff_2;
    bcoeff = 0.5 * fcoeff;
    scoeff = bcoeff + errorfc / cutoff_3;
    zcore = errorfc / cutoff + bcoeff * cutoff_2;
  }


  return 0;
}
int ZeroMultipoleSum::calc_zms_alpha0(real_pw& ene_ele, real_pw& grad_coeff,
				      const real_pw& r12,      const real_pw& r12_2,
				      const real_pw& r12_inv,   const real_pw& r12_2_inv,
				      const real_pw& r12_3_inv, const real_pw& cc){
  //#pragma optimize("", off)
  //real tmp1 = bcoeff * r12_2;
  //real tmp2 = tmp1 - zcore;
  //real tmp3 = tmp2 + r12_inv;
  //ene_ele = cc * tmp3;
  ene_ele = cc * (r12_inv - zcore + bcoeff * r12_2);
  //if(ene_ele > 1){
  //cout << "dbg 01 1 :"  << cc << " "<<  r12_inv << " " << zcore <<" " <<bcoeff << " " <<r12_2<<endl;
  //  }

  //cout << "DBG0203a: " << cc << " " << r12_inv << " "  << r12_2 << endl;
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
  grad_coeff = -cc * (r12_3_inv - fcoeff);
  //#pragma optimize("", on)
  return 0;
}
int ZeroMultipoleSum::calc_zms_excess_alpha0(real_pw& ene_ele, real_pw& grad_coeff,
					     const real_pw& r12,      const real_pw& r12_2,
					     const real_pw& r12_3_inv, const real_pw& cc){
  ene_ele = cc * bcoeff * r12_2;
  grad_coeff = cc * fcoeff;
  //cout << "DBG0203excess: " << cc << " " << bcoeff << " "  << r12_2 << endl;  
  return 0;
}
int ZeroMultipoleSum::calc_zms(real_pw& ene_ele, real_pw& grad_coeff,
			       const real_pw& r12,      const real_pw& r12_2,
			       const real_pw& r12_inv,   const real_pw& r12_2_inv,
			       const real_pw& r12_3_inv, const real_pw& cc){
  real_pw tmp = r12 * ewald_alpha;
  real_pw errorfc = erfc(tmp);
  ene_ele = cc * (r12_inv * errorfc - zcore + bcoeff * r12_2);
  grad_coeff = -cc * (errorfc * r12_3_inv + piewald * exp(tmp*tmp) * r12_2_inv - fcoeff);
  return 0;
}
int ZeroMultipoleSum::calc_zms_excess(real_pw& ene_ele, real_pw& grad_coeff,
				       const real_pw& r12,      const real_pw& r12_2,
				       const real_pw& r12_3_inv, const real_pw& cc){
  
  real_pw tmp = ewald_alpha * r12;
  real_pw errorfc = 1.0 - erfc(tmp);
  ene_ele = cc * (-errorfc * r12 + bcoeff * r12_2);
  grad_coeff = -cc * ( errorfc * r12_3_inv + piewald * exp(-tmp*tmp)
		       * r12 - fcoeff - r12_3_inv);
  return 0.0;
}
