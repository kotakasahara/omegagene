#include "ForceField.h"

ForceField::ForceField()
  : ForceFieldObject(){
  if(DBG >= 1)
    cout << "DBG1: ForceField::ForceField()"<<endl;
  
  gmul[0] = 0.0;  gmul[1] = 0.0;   gmul[2] = 2.0;
  gmul[3] = 0.0;  gmul[4] = 4.0;   gmul[5] = 0.0;
  gmul[6] = 6.0;  gmul[7] = 0.0;   gmul[8] = 8.0;
  gmul[9] = 0.0;  gmul[10] = 10.0;
  
  //if(DBG >= 1)
  //cout << "DBG1: ForceField::ForceField() owari"<<endl;    
}

ForceField::~ForceField(){
}

int ForceField::set_config_parameters(const Config* cfg){
  ForceFieldObject::set_config_parameters(cfg);
  switch(cfg->electrostatic){
  case ELCTRST_ZEROMULTIPOLE:
  case ELCTRST_ZERODIPOLE:
    ele = new ZeroMultipoleSum();
    break;
    //  case ELCTRST_DUMMY:
    //    ele = new ElectrostaticObject(mmsys);
    //    break;
  default:
    // ERROR
    break;
  }
  ele->set_config_parameters(cfg);
  
  return 0;
}

int ForceField::initial_preprocess(const PBC* in_pbc){
  if (DBG >= 1)
    cout << "DBG1: ForceField::initial_preprocess"<<endl;    
  ForceFieldObject::initial_preprocess(in_pbc);
  ele->initial_preprocess();
  return 0;
}

int ForceField::calc_bond(real_pw& ene, real_pw work[],
			  const real* crd1, const real* crd2,
			  const real& param_e, const real& param_r0){
  
  //mmsys->pbc.print_pbc();
  real d12[3];
  
  pbc->diff_crd_minim_image(d12, crd1, crd2);

  real_bp r12 = sqrt(d12[0]*d12[0] + d12[1]*d12[1] + d12[2]*d12[2]);
  //cout << "r12: " << r12 << endl;
  real_bp dif = r12 - param_r0;
  real_bp coef = param_e * dif;
  ene = coef * dif;
  real_bp grad_coef = 2.0 * coef / r12;
  for (int d=0; d < 3; d++)  work[d] = d12[d] * grad_coef;
  /*
  cout << setprecision(15);
  cout << "ForceField::calc_bond() r12:" << r12 << endl;
  cout << " param_e: " << param_e << "  param_r0: " << param_r0 << endl;
  cout << " coef: " << coef << " grad_coef:" << grad_coef << endl;
  cout << " r12: " << r12 << " diff: " << dif << " ene: " << ene << endl;
  */
  return 0;
}

int ForceField::calc_angle(real_pw& ene, real_pw work1[], real_pw work2[],
			   const real* crd1, const real* crd2, const real* crd3,
			   const real& param_e, const real& param_theta0){
  
  real d12[3];
  real d32[3];
  pbc->diff_crd_minim_image(d12, crd1, crd2);
  pbc->diff_crd_minim_image(d32, crd3, crd2);

  real_bp r12_2 = d12[0]*d12[0] + d12[1]*d12[1] + d12[2]*d12[2];
  real_bp r32_2 = d32[0]*d32[0] + d32[1]*d32[1] + d32[2]*d32[2];
  real_bp r12r32 = sqrt(r12_2 * r32_2);
  real_bp r12r32_inv = 1.0 / r12r32;
  real_bp cost = (d12[0]*d32[0] + d12[1]*d32[1] + d12[2]*d32[2]) * r12r32_inv;
  cost = max(min(cost, (real_bp)1.0), (real_bp)-1.0);
  real_bp theta = acos(cost);
  real_bp dif = theta - param_theta0;
  real_bp coef = param_e * dif;
  ene = coef * dif;  // Energy

  real_bp r12inv_2 = 1.0 / r12_2;
  real_bp r32inv_2 = 1.0 / r32_2;
  real_bp sin_theta = sin(theta);
  sin_theta = max((real_bp)EPS, sin_theta);

  real_bp grad_coef = - (2.0 * coef) / sin_theta;

  for (int i=0; i < 3; i++){
    work1[i] = grad_coef * (d32[i] * r12r32_inv - d12[i] * cost * r12inv_2);
    work2[i] = grad_coef * (d12[i] * r12r32_inv - d32[i] * cost * r32inv_2);
  }
  //cout << "angle work1 " << work1[0] << " " << work1[1] << " " << work1[2] << endl;
  //cout << "angle work2 " << work2[0] << " " << work2[1] << " " << work2[2] << endl;
  return 0;
}

int ForceField::calc_torsion(real_pw& ene, real_pw work1[], real_pw work2[], real_pw work3[],
			     const real* crd1, const real* crd2,
			     const real* crd3, const real* crd4,
			     const real& param_ene,      const real& param_overlaps,
			     const real& param_symmetry, const real& param_phase){
  ene = 0.0;
  for (int d=0; d<3; d++){
    work1[d] = 0.0; work2[d] = 0.0; work3[d] = 0.0;
  }
  if (param_ene <= EPS) return 0;

  real_bp d21[3];
  pbc->diff_crd_minim_image(d21, crd2, crd1);
  real_bp d32[3];
  pbc->diff_crd_minim_image(d32, crd3, crd2);
  real_bp d43[3];
  pbc->diff_crd_minim_image(d43, crd4, crd3);  
  
  real_bp p12[3];
  cross(d21, d32, p12);
  real_bp p23[3];
  cross(d32, d43, p23);
  real_bp p12_2 = p12[0]*p12[0] + p12[1]*p12[1] + p12[2]*p12[2];
  real_bp p23_2 = p23[0]*p23[0] + p23[1]*p23[1] + p23[2]*p23[2];
  real_bp p12p23 = sqrt(p12_2*p23_2);
  if (p12p23 <= 0.0) return 0;

  //cout << "dbg0204:p12 " << p12[0] << " " << p12[1] << " " << p12[2] << endl;
  real_bp p12p23inv = 1.0 / p12p23;
  real_bp s1223 = p12[0]*p23[0] + p12[1]*p23[1] + p12[2]*p23[2];
  real_bp cosp = s1223 * p12p23inv;
  cosp = max((real_bp)-1.0, min((real_bp)1.0, cosp));
  real_bp p123[3];
  cross(p12, p23, p123);
  //cout << "dbg0204:cosp " << cosp << " " << s1223 << " " <<  p12p23inv << endl;
  real_bp phi = acos(cosp);
  real_bp s32123 = d32[0]*p123[0] + d32[1]*p123[1] + d32[2]*p123[2];
  if (s32123 < 0.0) phi = 2.0 * PI - phi;
  real_bp phimod = param_symmetry * phi - param_phase;
  //cout << "dbg0204:ps, phi, phase " << param_symmetry << " ";
  //cout << phi << " " << param_phase << endl;
  
  real_bp div = 1.0 / real_bp(param_overlaps);
  ene = param_ene * (1.0 + cos(phimod) ) * div;
  
  //cout << "dbg0204:torsion " << ene << " " << param_ene << " ";
  //cout << phimod << " " << div << endl;
  real_bp sinp = sin(phi);
  int nrot = (int)param_symmetry;
  real_bp sinp_sq = sinp * sinp;
  real_bp cosp_sq = cosp * cosp;

  real_bp sin_fphs = sin(param_phase);
  real_bp cos_fphs = cos(param_phase);
  real_bp sin_rot = sin(param_symmetry * phi);
  real_bp cos_rot = cos(param_symmetry * phi);
  
  real_bp dums = sinp;
  real_bp dflim = cos_fphs * (param_symmetry - gmul[nrot] + gmul[nrot] * cosp);
  real_bp df0 = cos_fphs * sin_rot - sin_fphs * cos_rot;
  real_bp df1 = dflim;
  if (EPS <= abs(dums))
    df1 = df0 / dums;

  real_bp work_coeff = param_symmetry * param_ene * div * df1 * p12p23inv;
  real_bp r12 = 1.0 / p12_2 * s1223;
  real_bp r23 = 1.0 / p23_2 * s1223;
  real_bp op_p23d32[3];
  cross(p23, d32, op_p23d32);
  real_bp op_p12d32[3];
  cross(p12, d32, op_p12d32);

  for(int i=0; i<3; i++)
    work1[i] = work_coeff * (op_p23d32[i] - r12 * op_p12d32[i]);

  real_bp d2132[3];
  d2132[0] = d21[0] + d32[0];
  d2132[1] = d21[1] + d32[1];
  d2132[2] = d21[2] + d32[2];
  real_bp op_d2132_p23[3];
  cross(d2132, p23, op_d2132_p23);
  real_bp op_d2132_p12[3];
  cross(d2132, p12, op_d2132_p12);
  real_bp op_p12_d43[3];
  cross(p12, d43, op_p12_d43);
  real_bp op_p23_d43[3];
  cross(p23, d43, op_p23_d43);

  for(int i=0; i<3; i++)
    work2[i] = work_coeff * (op_d2132_p23[i] + op_p12_d43[i]
			     - r12 * op_d2132_p12[i]
			     - r23 * op_p23_d43[i]);
  
  real_bp op_p12_d32[3];
  cross(p12, d32, op_p12_d32);
  real_bp op_p23_d32[3];
  cross(p23, d32, op_p23_d32);

  for(int i=0; i<3; i++)
    work3[i] = work_coeff * (op_p12_d32[i] - r23 * op_p23_d32[i]);
  
  //cout << "dbg0204:torsion:work: "<<endl;
  //cout << work1[0] << " " << work1[1] << " " << work1[2] << endl;
  //cout << work2[0] << " " << work2[1] << " " << work2[2] << endl;
  //cout << work3[0] << " " << work3[1] << " " << work3[2] << endl;

  return 0;
}

int ForceField::calc_14pair(real_pw& ene_vdw,
			    real_pw& ene_ele,
			    real_fc work[],
			    const real* crd1, const real* crd4,
			    const real& lj_6term,
			    const real& lj_12term,
			    const real& charge1,
			    const real& charge4,
			    const real& param_coeff_vdw,
			    const real& param_coeff_ele){
  
  real d14[3];
  pbc->diff_crd_minim_image(d14, crd1, crd4);
  real r14_2 = d14[0] * d14[0] + d14[1]*d14[1] + d14[2]*d14[2];
  real r14_2_inv = 1.0 / r14_2;
  real r14_inv = sqrt(r14_2_inv);

  // VDW
  real r14_inv_6 = r14_2_inv * r14_2_inv * r14_2_inv;
  real r14_inv_12 = r14_inv_6 * r14_inv_6;
  real term6 = lj_6term * r14_inv_6;
  real term12 = lj_12term * r14_inv_12;
  ene_vdw = (-term6 + term12) * param_coeff_vdw;
  //cout << "dbg0204: t6,t23,coef: " << term6 << " " << term12 << " " << param_coeff_vdw << " " << sqrt(r14_2)<<endl;
  real coef_work = r14_2_inv * param_coeff_vdw * ( -12.0 * term12 + 6.0 * term6);
  real_bp work_vdw[3];
  for (int d=0; d<3; d++)
    work_vdw[d] = coef_work*d14[d];

  // ELE
  real cc = charge1 * charge4 * CHARGE_COEFF * param_coeff_ele;
  ene_ele = cc * r14_inv;
  real coef_ele = -cc * r14_2_inv * r14_inv;
  real_bp work_ele[3];
  for (int d=0; d<3; d++)
    work_ele[d] = coef_ele*d14[d];
  
  for (int d=0; d<3; d++)
    work[d] = work_vdw[d] + work_ele[d];
  
  return 0;
}

int ForceField::calc_pairwise(real_pw& ene_vdw, real_pw& ene_ele,
			      real_fc work[],
			      const real* crd1, const real* crd2,
			      const real& param_6term,
			      const real& param_12term,
			      const real& charge1,
			      const real& charge2){
  
  real d12[3] = {0.0, 0.0, 0.0};

  //pbc->diff_crd_minim_image(d12, crd1, crd2);
  
  for(int d=0; d<3; d++){
    d12[d] = crd1[d] - crd2[d];
  }

  real r12_2 = d12[0]*d12[0] + d12[1]*d12[1] + d12[2]*d12[2];
  
  real r12 = sqrt(r12_2);
  //  cout << "r12 : " << r12 << endl;
  ene_vdw = 0.0;
  ene_ele = 0.0;
  real work_vdw[3] = {0.0, 0.0, 0.0};
  real work_ele[3] = {0.0, 0.0, 0.0};
  //if( r12 < 0.1) 
  //cout << "  r12 " << r12 << endl;
  if (r12 >= cutoff) return 1;
  real r12_inv = 1.0 / r12;
  real r12_2_inv = 1.0 / r12_2;
  real r12_3_inv = r12_inv * r12_2_inv;
  real r12_6_inv = r12_3_inv * r12_3_inv;
  real r12_12_inv = r12_6_inv * r12_6_inv;
  real term6 = param_6term * r12_6_inv;
  real term12 = param_12term * r12_12_inv;
  ene_vdw = -term6 + term12;
  real_pw work_coef = r12_2_inv * (-12.0 * term12 + 6.0 * term6);
  for (int d=0; d<3; d++)
    work_vdw[d] = work_coef * d12[d];

  real cc = charge1 * charge2 * CHARGE_COEFF;

  real work_coef_ele;
  //cout << "dbg0204: cc " << cc << " " <<  charge1 << " " << charge2 << " " << r12_2 << " " << r12_inv << endl;
  //printf("dbgcrd %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n",
  //crd1[0],crd1[1],crd1[2],crd2[0],crd2[1],crd2[2]);
  (ele->*(ele->func_calc_zms))(ene_ele, work_coef_ele,
			       r12, r12_2, r12_inv,
			       r12_2_inv, r12_3_inv, cc);
  //printf("dbgpair %10e %10e %15e\n", r12, cc, CHARGE_COEFF);
  for(int d=0; d<3; d++)
    work_ele[d] = work_coef_ele * d12[d];

  for(int d=0; d<3; d++)
    work[d] = work_vdw[d] + work_ele[d];

  return 0;
}
int ForceField::calc_zms_excess(real_pw& ene, real_pw work[],
				const real* crd1,
				const real* crd2,
				const real& charge1,
				const real& charge2){
  
  real d12[3];
  pbc->diff_crd_minim_image(d12, crd1, crd2);
  real r12_2 = d12[0]*d12[0] + d12[1]*d12[1] + d12[2]*d12[2];
  real r12 = sqrt(r12_2);
  real r12_inv = 1.0 / r12;
  real r12_3_inv = r12_inv * r12_inv * r12_inv;
  real cc = charge1 * charge2 * CHARGE_COEFF;
  real work_coef;
  //cout << "dbg0204: cc " << cc << " " <<  charge1 << " " << charge2 << " " << r12_2 << " " << r12_inv << endl;
  (ele->*(ele->func_calc_zms_excess))(ene, work_coef,
				      r12, r12_2, r12_3_inv, cc);
  for (int d=0; d<3; d++)
    work[d] = work_coef * d12[d];
  return 0;
}

int ForceField::cal_self_energy(const int& n_atoms,
				const int& n_excess,
				const int**& excess_pairs,
				/*const int& n_atoms,
				const int& n_bonds,
				const int**& bond_atomid_pairs,
				const int& n_angles,
				const int**& angle_atomid_triads,
				const int& n_torsions,
				const int**& torsion_atomid_quads,
				const int*& torsion_nb14,
				*/
				const real_pw*& charge,
				real*& energy_self,
				real& energy_self_sum){
  ele->cal_self_energy(n_atoms, n_excess, excess_pairs,
		       /*n_atoms, 
		       n_bonds,    bond_atomid_pairs,
		       n_angles,   angle_atomid_triads,
		       n_torsions, torsion_atomid_quads,
		       torsion_nb14,*/
		       charge,
		       energy_self, energy_self_sum);
  return 0;
}
