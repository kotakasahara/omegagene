#include "MmSystem.h"

MmSystem::MmSystem()
  : CelesteObject(){
  max_n_nb15off = 64;
  ctime_per_step = 0;
  ctime_cuda_htod_atomids = 0;
  ctime_cuda_reset_work_ene = 0;
  ctime_calc_energy = 0;
  ctime_calc_energy_pair = 0;
  ctime_calc_energy_bonded = 0;
  ctime_update_velo = 0;
  ctime_calc_kinetic = 0;
  ctime_update_coord = 0;
}
MmSystem::~MmSystem(){
  free_all();
}

int MmSystem::free_all(){
  if(DBG==1) cout << "free_atom_vars()" << endl;
  if (n_atoms > 0) free_atom_vars();
  if(DBG==1) cout << "free_lj_params()" << endl;
  if (n_lj_type_pairs > 0) free_lj_params();
  if(DBG==1) cout << "free_bonds()" << endl;
  if (n_bonds > 0)    free_bonds();
  if(DBG==1) cout << "free_angles()" << endl;
  if (n_angles > 0)   free_angles();
  if(DBG==1) cout << "free_torsions()" << endl;
  if (n_torsions > 0) free_torsions();
  if(DBG==1) cout << "free_impros()" << endl;
  if (n_impros > 0)   free_impros();
  if(DBG==1) cout << "free_nb14()" << endl;
  if (n_nb14 > 0)     free_nb14();
  if(DBG==1) cout << "free_excess_pairs()" << endl;
  if (n_excess > 0)     free_excess_pairs();
  //if(DBG==1) cout << "free_pcluster_vars()" << endl;
  //if (n_nb14 > 0)     free_pcluster_vars();

  if(DBG==1) cout << "free_nb15off()" << endl;
  if (n_nb15off > 0)  free_nb15off();
  return 0;
}
int MmSystem::alloc_atom_vars(){
  crd = new real*[n_atoms];
  for(int i=0; i < n_atoms; i++){
    crd[i] = new real[3];
    for(int d=0; d < 3; d++) crd[i][d] = 0.0;
  }
  
  force = new real_fc*[n_atoms];
  for(int i=0; i < n_atoms; i++){
    force[i] = new real_fc[3];
    for(int d=0; d < 3; d++) force[i][d] = 0.0;
  }
  //vel = new real*[n_atoms];
  //for(int i=0; i < n_atoms; i++){
  //vel[i] = new real[3];
  //for(int d=0; d < 3; d++) vel[i][d] = 0.0;
  //}
  vel_just = new real*[n_atoms];
  for(int i=0; i < n_atoms; i++){
    vel_just[i] = new real[3];
    for(int d=0; d < 3; d++) vel_just[i][d] = 0.0;
  }
  //vel_next = new real*[n_atoms];
  //for(int i=0; i < n_atoms; i++){
  //vel_next[i] = new real[3];
  //for(int d=0; d < 3; d++) vel_next[i][d] = 0.0;
  //}

  charge = new real_pw[n_atoms];
  mass = new real[n_atoms];
  atom_type = new int[n_atoms];
  energy_self = new real[n_atoms];
  for(int i=0; i < n_atoms; i++){
    charge[i] = 0.0;
    mass[i] = 0.0;
    atom_type[i] = 0;
    energy_self[i] = 0.0;
  }
  return 0;
}

int MmSystem::alloc_lj_params(){
  lj_6term = new real_pw[n_lj_types * n_lj_types];
  lj_12term = new real_pw[n_lj_types * n_lj_types];
  return 0;
}
int MmSystem::alloc_bonds(){
  bond_epsiron = new real[n_bonds];
  bond_r0 = new real[n_bonds];
  bond_atomid_pairs = new int*[n_bonds];
  for(int i=0; i < n_bonds; i++){
    bond_atomid_pairs[i] = new int[2];
    bond_atomid_pairs[i][0] = 0;
    bond_atomid_pairs[i][1] = 0;
    bond_epsiron[i] = 0.0;
    bond_r0[i] = 0.0;
  }
  return 0;
}
int MmSystem::alloc_angles(){
  angle_epsiron = new real[n_angles];
  angle_theta0 = new real[n_angles];
  angle_atomid_triads = new int*[n_angles];
  for(int i=0; i < n_angles; i++){
    angle_atomid_triads[i] = new int[3];
    angle_atomid_triads[i][0] = 0;
    angle_atomid_triads[i][1] = 0;
    angle_atomid_triads[i][2] = 0;
    angle_epsiron[i] = 0.0;
    angle_theta0[i] = 0.0;
  }
  return 0;
}
int MmSystem::alloc_torsions(){
  torsion_energy = new real[n_torsions];
  torsion_overlaps = new int[n_torsions];
  torsion_symmetry = new int[n_torsions];
  torsion_phase = new real[n_torsions];
  torsion_nb14 = new int[n_torsions];
  torsion_atomid_quads = new int*[n_torsions];
  for(int i=0; i<n_torsions; i++){
    torsion_atomid_quads[i] = new int[4];
    torsion_atomid_quads[i][0] = 0;
    torsion_atomid_quads[i][1] = 0;
    torsion_atomid_quads[i][2] = 0;
    torsion_atomid_quads[i][3] = 0;
    torsion_energy[i] = 0.0;
    torsion_overlaps[i] = 0;
    torsion_symmetry[i] = 0;
    torsion_phase[i] = 0.0;
    torsion_nb14[i] = 0;
  }
  return 0;
}
int MmSystem::alloc_impros(){
  impro_energy = new real[n_impros];
  impro_overlaps = new int[n_impros];
  impro_symmetry = new int[n_impros];
  impro_phase = new real[n_impros];
  impro_nb14 = new int[n_impros];
  impro_atomid_quads = new int*[n_impros];
  for(int i=0; i<n_impros; i++){
    impro_atomid_quads[i] = new int[4];
    impro_atomid_quads[i][0] = 0;
    impro_atomid_quads[i][1] = 0;
    impro_atomid_quads[i][2] = 0;
    impro_atomid_quads[i][3] = 0;
    impro_energy[i] = 0.0;
    impro_overlaps[i] = 0;
    impro_symmetry[i] = 0;
    impro_phase[i] = 0.0;
    impro_nb14[i] = 0;
  }
  return 0;
}
int MmSystem::alloc_nb14(){
  nb14_coeff_vdw = new real[n_nb14];
  nb14_coeff_ele = new real[n_nb14];
  nb14_atomid_pairs = new int*[n_nb14];
  nb14_atomtype_pairs = new int*[n_nb14];
  for(int i=0; i<n_nb14; i++){
    nb14_atomid_pairs[i] = new int[2];
    nb14_atomid_pairs[i][0] = 0;
    nb14_atomid_pairs[i][1] = 0;
    nb14_atomtype_pairs[i] = new int[2];
    nb14_atomtype_pairs[i][0] = 0;
    nb14_atomtype_pairs[i][1] = 0;
    nb14_coeff_vdw[i] = 0.0;
    nb14_coeff_ele[i] = 0.0;
  }
  return 0;
}

int MmSystem::alloc_nb15off(){
  nb15off = new int[n_atoms * max_n_nb15off];
  nb15off1 = new int[n_atoms];
  nb15off2 = new int[n_atoms];
  n_nb15off = 0;
  for(int i=0; i < n_atoms; i++){
    //debug
    for(int j=0; j < max_n_nb15off; j++){
      nb15off[i*max_n_nb15off+j] = -1;
    }
   nb15off1[i] = 0;
    nb15off2[i] = 0;
  }
  return 0;
}

/*
int MmSystem::alloc_pcluster_vars(){
  atom_pclusters = new int[n_atoms];
  pcluster_atoms = new int[n_atoms];
  n_atoms_pclusters = new int[n_pclusters];
  pclusters_index = new int[n_pclusters];
  return 0;
}
*/

int MmSystem::free_atom_vars(){
  cout << "free atom vars" << endl;
  for(int i=0; i < n_atoms; i++){
    delete[] crd[i];
    delete[] force[i];
    //delete[] vel[i];
    delete[] vel_just[i];
    //delete[] vel_next[i];
  }
  delete[] crd;
  delete[] force;
  //  delete[] vel;
  delete[] vel_just;
  //delete[] vel_next;
  delete[] charge;
  delete[] mass;
  delete[] atom_type;
  delete[] energy_self;
  return 0;
}
int MmSystem::free_lj_params(){
  delete[] lj_6term;
  delete[] lj_12term;
  return 0;
}
int MmSystem::free_bonds(){
  for(int i=0; i < n_bonds; i++){
    delete[] bond_atomid_pairs[i];
  }
  delete[] bond_atomid_pairs;
  delete[] bond_epsiron;
  delete[] bond_r0;
  return 0;
}
int MmSystem::free_angles(){
  for(int i=0; i < n_angles; i++){
    delete[] angle_atomid_triads[i];
  }
  delete[] angle_atomid_triads;
  delete[] angle_epsiron;
  delete[] angle_theta0;
  return 0;
}
int MmSystem::free_torsions(){
  for(int i=0; i < n_torsions; i++){
    delete[] torsion_atomid_quads[i];
  }
  delete[] torsion_atomid_quads;
  delete[] torsion_energy;
  delete[] torsion_overlaps;
  delete[] torsion_symmetry;
  delete[] torsion_phase;
  delete[] torsion_nb14;
  return 0;
}
int MmSystem::free_impros(){
  for(int i=0; i < n_impros; i++){
    delete[] impro_atomid_quads[i];
  }
  delete[] impro_atomid_quads;
  delete[] impro_energy;
  delete[] impro_overlaps;
  delete[] impro_symmetry;
  delete[] impro_phase;
  delete[] impro_nb14;
  return 0;

}
int MmSystem::free_nb14(){
  for(int i=0; i < n_nb14; i++){
    delete[] nb14_atomid_pairs[i];
    delete[] nb14_atomtype_pairs[i];
  }
  delete[] nb14_atomid_pairs;
  delete[] nb14_atomtype_pairs;
  delete[] nb14_coeff_vdw;
  delete[] nb14_coeff_ele;
  return 0;
}

int MmSystem::free_nb15off(){
  delete[] nb15off;
  delete[] nb15off1;
  delete[] nb15off2;

  return 0;
}

// Parameter Setter
int MmSystem::set_lj_pair_param(int type1, int type2, real_pw param6, real_pw param12){
  lj_6term[type1 * n_lj_types + type2] = param6;
  lj_6term[type2 * n_lj_types + type1] = param6;
  lj_12term[type1 * n_lj_types + type2] = param12;
  lj_12term[type2 * n_lj_types + type1] = param12;
  return 0;
}
int MmSystem::set_bond_param(int bond_id, 
			     int atomid1, int atomid2,
			     real eps, real r0){
  bond_atomid_pairs[bond_id][0] = atomid1;
  bond_atomid_pairs[bond_id][1] = atomid2;
  bond_epsiron[bond_id] = eps;
  bond_r0[bond_id] = r0;

  //  if(DBG>=1)
  //cout << "setbond_param " << bond_id << " : "
  // << bond_atomid_pairs[bond_id][0] << " - "
  // << bond_atomid_pairs[bond_id][1];

  return 0;
}
int MmSystem::set_angle_param(int angle_id, 
			     int atomid1, int atomid2, int atomid3,
			     real eps, real theta0){
  angle_atomid_triads[angle_id][0] = atomid1;
  angle_atomid_triads[angle_id][1] = atomid2;
  angle_atomid_triads[angle_id][2] = atomid3;
  angle_epsiron[angle_id] = eps;
  angle_theta0[angle_id] = theta0;
  return 0;
}
int MmSystem::set_torsion_param(int torsion_id, 
				int atomid1, int atomid2, int atomid3, int atomid4,
				real ene, int overlaps,
				int symmetry, real phase,
				int flag_14nb){
  torsion_atomid_quads[torsion_id][0] = atomid1;
  torsion_atomid_quads[torsion_id][1] = atomid2;
  torsion_atomid_quads[torsion_id][2] = atomid3;
  torsion_atomid_quads[torsion_id][3] = atomid4;
  torsion_energy[torsion_id] = ene;
  torsion_overlaps[torsion_id] = overlaps;
  torsion_symmetry[torsion_id] = symmetry;
  torsion_phase[torsion_id] = phase;
  torsion_nb14[torsion_id] = flag_14nb;
  return 0;
}
int MmSystem::set_impro_param(int impro_id, 
			      int atomid1, int atomid2, int atomid3, int atomid4,
			      real ene, int overlaps,
			      int symmetry, real phase,
			      int flag_14nb){
  impro_atomid_quads[impro_id][0] = atomid1;
  impro_atomid_quads[impro_id][1] = atomid2;
  impro_atomid_quads[impro_id][2] = atomid3;
  impro_atomid_quads[impro_id][3] = atomid4;
  impro_energy[impro_id] = ene;
  impro_overlaps[impro_id] = overlaps;
  impro_symmetry[impro_id] = symmetry;
  impro_phase[impro_id] = phase;
  impro_nb14[impro_id] = flag_14nb;
  return 0;
}

int MmSystem::set_nb14_param(int nb14_id, 
			     int atomid1, int atomid2,
			     int atomtype1, int atomtype2,
			     real coeff_vdw, real coeff_ele){
  nb14_atomid_pairs[nb14_id][0] = atomid1;
  nb14_atomid_pairs[nb14_id][1] = atomid2;
  nb14_atomtype_pairs[nb14_id][0] = atomtype1;
  nb14_atomtype_pairs[nb14_id][1] = atomtype2;
  nb14_coeff_vdw[nb14_id] = coeff_vdw;
  nb14_coeff_ele[nb14_id] = coeff_ele;
  return 0;
}
void show_int(int x){
  int i;
  printf("%4d  %08x  ", x, (unsigned char)x);
  for (i = sizeof(int)*8-1; i >= 0; i--) {
    printf("%d", (x>>i) & 1);
    if(i%8==0) printf(" ");
  }
  printf("\n");
}

bool MmSystem::search_nb15off(int atomid1, int atomid2){
  // if the atompair (atomid1, atomid2) is 
  // in the list of nb15off, the function returns false
  bool flg1 = false;
  int aid_diff = atomid2 - atomid1;
  int bit_idx = atomid1;
  if(aid_diff < 0){
    aid_diff = -aid_diff;
    bit_idx = atomid2;
  }
  int mask1 = 0;
  int mask2 = 0;
  if(aid_diff <= 32)
    mask1 = 1 << (aid_diff-1);
  else if(aid_diff > 32 && aid_diff <= 64)
    mask2 = 1 << (aid_diff-33);
  if(aid_diff > 0 && mask1 != 0 &&
     (mask1 & nb15off1[bit_idx]) == mask1) {flg1=true;}
  if(aid_diff > 32 && mask2 != 0 &&
     (mask2 & nb15off2[bit_idx]) == mask2) {flg1=true;}
  if(aid_diff > 64) flg1=false;

  bool flg2=false;
  int tail = max_n_nb15off * atomid1 + max_n_nb15off;
  for(int i = max_n_nb15off * atomid1; i < tail; i++){
    if (nb15off[i] == atomid2)
      flg2=true;
  }
  /*
  if (flg1 != flg2){
    cout << "DBG: search_nb15off " << atomid1 << "-" << atomid2 << endl;
    if(flg1) cout << "new algo TRUE" << endl;
    if(flg2) cout << "prev algo TRUE" << endl;
    cout << "bit_idx: " << bit_idx << " aid_diff:" << aid_diff << endl;
    cout << "mask1 ";
    show_int(mask1);
    cout << "mask2 ";
    show_int(mask2);
    cout << "nb15off1 ";
    show_int(nb15off1[bit_idx]);
    cout << "nb15off2 ";
    show_int(nb15off2[bit_idx]);
    cout << "and 1 ";
    show_int(mask1 & nb15off1[bit_idx]);
    cout << "and 2 ";
    show_int(mask2 & nb15off2[bit_idx]);
    if((mask1 & nb15off1[bit_idx]) == mask1) cout << "mask1 true" <<endl;
    if((mask2 & nb15off2[bit_idx]) == mask2) cout << "mask2 true" <<endl;
    cout << "oldpair "; 
    for(int i = max_n_nb15off * atomid1; i < tail; i++){
      cout <<  nb15off[i] <<" ";
    }
    cout << endl;
  }
*/
  return flg1;
  //return flg2;
}

int MmSystem::set_nb15off(int atomid1, int atomid2){
  //cout << "set_nb15off"<< endl;
  //debug
  int i = 0;
  for(i=0; i<max_n_nb15off; i++)
    if(nb15off[atomid1*max_n_nb15off+i]==-1 ||
       nb15off[atomid1*max_n_nb15off+i]==atomid2 ) break;
  nb15off[atomid1 * max_n_nb15off + i] = atomid2;
  for(i=0; i<max_n_nb15off; i++)
    if(nb15off[atomid2*max_n_nb15off+i]==-1 ||
       nb15off[atomid2*max_n_nb15off+i]==atomid1 ) break;
  nb15off[atomid2 * max_n_nb15off + i] = atomid1;
  n_nb15off+=2;

  if(atomid1 > atomid2){
    int tmp = atomid1;
    atomid1 = atomid2;
    atomid2 = tmp;
  }else if(atomid1 == atomid2) return 1;
  int id_diff = atomid2 - atomid1;
  int add_bit1 = 0;
  int add_bit2 = 0;
  if(id_diff <= 32){
    add_bit1 = 1 << (id_diff-1);
  }else if(id_diff <= 64){
    add_bit2 = 1 << (id_diff-33);
  }else{
    cerr << "The atom pair [" << atomid2 << "-" << atomid1 << "] " ;
    cerr << "was within four covalent bonds and in exclusion list of non bonded pair potential. ";
    cerr << "However this software requires that differences of ";
    cerr << "atom ids of pairs in the exclusion list ";
    cerr << "must be equal or less than 64. ";
    cerr << "Reordering of atoms in the topology file and strucure file are needed.";
    cerr << endl;
    exit(1);
  }
  int tmp1 = nb15off1[atomid1];
  int tmp2 = nb15off2[atomid1];
  nb15off1[atomid1] = tmp1 | add_bit1;
  nb15off2[atomid1] = tmp2 | add_bit2;

  //if(nb15off1[atomid1] != tmp1 || nb15off2[atomid1] != tmp2) n_nb15off++;
  //  cout << "//set_nb15off"<< endl;
  return 0;
}

int MmSystem::write_data(){
  cout << "launchset_version: " << launchset_version << endl;
  cout << "pbc:";
  for (int i = 0; i < 9; i++){
    cout << " " << pbc_val[i];
  }
  cout << endl;
  cout << "n_lj_types: " << n_lj_types << endl;
  cout << "n_lj_type_pairs: " << n_lj_type_pairs << endl;
  cout << "n_bonds: " << n_bonds << endl;
  cout << "n_angles: " << n_angles << endl;
  cout << "n_torsions: " << n_torsions << endl;
  cout << "n_impros: " << n_impros << endl;
  cout << "n_nb14: " << n_nb14 << endl;
  cout << "n_nb15off: " << n_nb15off << endl;
  return 0;
}

int MmSystem::alloc_excess_pairs(){
  max_n_excess = n_bonds + n_angles + n_torsions;// + n_impros;

  if(DBG>=1) cout << "alloc_excess_pairs " << max_n_excess << endl;
  excess_pairs = new int*[max_n_excess];
  for(int i=0; i < max_n_excess; i++){
    excess_pairs[i] = new int[2];
    excess_pairs[i][0] = -1;
    excess_pairs[i][1] = -1;
  }
  return 0;
}
int MmSystem::free_excess_pairs(){
  for(int i=0; i < max_n_excess; i++){
    delete[] excess_pairs[i];
  }
  delete[] excess_pairs;
  return 0;
}
int MmSystem::add_excess_pairs(int atomid1, int atomid2){
  bool flg=true;
  for(int j=0; j < n_excess; j++){
    if(excess_pairs[j][0] == atomid1 && 
       excess_pairs[j][1] == atomid2){
      flg=false;
      return 1;
    }
  }
  if(flg){
    excess_pairs[n_excess][0] = atomid1;
    excess_pairs[n_excess][1] = atomid2;
    n_excess++;
  }
  return 0;
}
int MmSystem::set_excess_pairs(){
  n_excess = 0;
  for(int i=0; i < n_bonds; i++){
    int atomid1 = bond_atomid_pairs[i][0];
    int atomid2 = bond_atomid_pairs[i][1];
    add_excess_pairs(atomid1, atomid2);
  }
  for(int i=0; i < n_angles; i++){
    int atomid1 = angle_atomid_triads[i][0];
    int atomid2 = angle_atomid_triads[i][2];
    add_excess_pairs(atomid1, atomid2);
  }
  for(int i=0; i < n_torsions; i++){
    int atomid1 = torsion_atomid_quads[i][0];
    int atomid2 = torsion_atomid_quads[i][3];
    add_excess_pairs(atomid1, atomid2);
  }
  //for(int i=0; i < n_impros; i++,j++){
  //int atomid1 = impro_atomid_quads[i][0];
  //int atomid2 = impro_atomid_quads[i][3];
  //excess_pairs[j][0] = atomid1;
  //excess_pairs[j][1] = atomid2;
  //  }
  //cout << "excess set " << n_excess << " / "  << n_max_excess << endl;
  return 0;  
}

int MmSystem::reset_energy(){
  pote_bond = 0.0;
  pote_angle = 0.0;
  pote_torsion = 0.0;
  pote_impro = 0.0;
  pote_14vdw = 0.0;
  pote_14ele = 0.0;
  pote_vdw = 0.0;
  pote_ele = 0.0;
  for(int atomid=0; atomid < n_atoms; atomid++)
    for(int d=0; d<3; d++) force[atomid][d] = 0.0;

  return 0;
}
/*
int MmSystem::velocity_swap(){
  real** tmp = vel;
  vel = vel_next;
  vel_next = tmp;
  return 0;
}

int MmSystem::velocity_average(){
  for(int i=0; i < n_atoms; i++){
    for(int d=0; d < 3; d++)
      vel_just[i][d] = (vel[i][d] + vel_next[i][d])*0.5;
  }
  return 0;
}
*/
int MmSystem::revise_crd_inbox(){
  for(int i=0; i < n_atoms; i++){
    for(int d=0; d < 3; d++){
      while(crd[i][d] < pbc.lower_bound[d]) crd[i][d] += pbc.L[d];
      while(crd[i][d] >= pbc.upper_bound[d]) crd[i][d] -= pbc.L[d];
    }
  }
  return 0;
}

int MmSystem::set_random(int seed){
  //random_mt.seed(seed);
  //srand(seed);
  return 0;
}


int MmSystem::output_ctimes(){
  cout << "[Ctime] step: " << ctime_per_step / (double)CLOCKS_PER_SEC <<endl;
  cout << "[Ctime] cuda_htod_atomids: " << ctime_cuda_htod_atomids / (double)CLOCKS_PER_SEC <<endl;
  cout << "[Ctime] cuda_reset_work_ene: " << ctime_cuda_reset_work_ene / (double)CLOCKS_PER_SEC <<endl;
  cout << "[Ctime] calc_energy: " << ctime_calc_energy / (double)CLOCKS_PER_SEC <<endl;
  cout << "[Ctime] calc_energy_pair: " << ctime_calc_energy_pair / (double)CLOCKS_PER_SEC <<endl;
  cout << "[Ctime] calc_energy_bonded: " << ctime_calc_energy_bonded / (double)CLOCKS_PER_SEC <<endl;
  cout << "[Ctime] update_velo: " << ctime_update_velo / (double)CLOCKS_PER_SEC <<endl;
  cout << "[Ctime] calc_kinetic: " << ctime_calc_kinetic / (double)CLOCKS_PER_SEC <<endl;
  cout << "[Ctime] update_coord: " << ctime_update_coord / (double)CLOCKS_PER_SEC <<endl;
  cout << "[Ctime] set_grid: " << ctime_setgrid / (double)CLOCKS_PER_SEC <<endl;
  cout << "[Ctime] enumerate_cellpairs: " << ctime_enumerate_cellpairs / (double)CLOCKS_PER_SEC <<endl;

  return 0;
}

int MmSystem::ff_setup(const Config* cfg){
  ff = ForceField();
  ff.set_config_parameters(cfg);
  ff.initial_preprocess((const PBC*)&pbc);
  ff.cal_self_energy((const int&)n_atoms,
		     (const int&)n_bonds,
		     (const int**&)bond_atomid_pairs,
		     (const int&)n_angles,
		     (const int**&)angle_atomid_triads,
		     (const int&)n_torsions,
		     (const int**&)torsion_atomid_quads,
		     (const int*&)torsion_nb14,
		     (const real_pw*&)charge,
		     energy_self,
		     energy_self_sum);
  return 0;
}
