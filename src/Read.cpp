#include "Read.h"

const int Read::MAX_LEN_NAME = 256;

Read::Read(string inFn)
  : CelesteObject(){
  op=false;
  filename = inFn;
}

int Read::open(){
  ifs.open(filename.c_str());
  if(!ifs){
    cerr <<"Cannot open "<< filename << "." <<endl;
    return 1;
  }
  op=true;
  set_conv_endian_false();
  return 0;
}
int Read::close(){
  ifs.close();
  op=false;
  return 0;
}

vector<string> Read::load_config(){
  vector<string> vconf;
  open();
  string buf;
  while(ifs && getline(ifs, buf)){
    stringstream ss(buf);
    while(ss>>buf){
      if(buf[0] == '#' || buf[0] == ';') break;
      vconf.push_back(buf);
    }
  }
  close();
  return vconf;
}

vector<int> Read::load_integers(){
  vector<int> intvec;
  open();
  int tmp;
  string buf;
  while(getline(ifs,buf)){
    intvec.push_back(atoi(buf.c_str()));
  }
  close();
  return intvec;
}
vector<string> Read::load_strings(){
  vector<string> strvec;
  open();
  string buf;
  while(getline(ifs,buf)){
    stringstream ss(buf);
    string str;
    ss>>str;
    strvec.push_back(str);
  }
  close();
  return strvec;
}

int Read::load_launch_set(MmSystem& mmsys){
  open();
  cout << "-- Load Celeste input package." << endl;
  load_ls_header(mmsys);
  if(size_box>0){
    cout << "--- Load PBC definition   : " << size_box << " bytes."<< endl;
    load_ls_box(mmsys);
  }
  if(size_crd>0){
    cout << "--- Load atom coordinates : " << size_crd << " bytes." << endl;
    load_ls_crd(mmsys);
  }
  if (size_vel > 0){
    cout << "--- Load atom velocities  : " << size_vel << " bytes." << endl;
    load_ls_vel(mmsys);
  }
  if(size_topol > 0){
    cout << "--- Load topology data    : " << size_topol << " bytes." << endl;
    load_ls_tpl(mmsys);
  }
  if(size_constraint > 0){
    cout << "--- Load constraint definition : " << size_constraint << " bytes." << endl;
    load_ls_constraint(&mmsys.constraint);
  }
  if(size_expand > 0){
    cout << "--- Load expand ensemble definition : " << size_expand << " bytes." << endl;
    load_ls_vmcmd(&mmsys.vmcmd);
  }
  if(size_groups > 0){
    cout << "--- Load atom group definition : " << size_groups << " bytes." << endl;
    load_ls_atom_groups(mmsys);
  }
  //cout << "load_ls_pcluster()" << endl;
  //load_ls_pcluster(mmsys);
  close();
  return 0;
}

int Read::load_ls_header(MmSystem& mmsys){
  cout << "--- Load file header." << endl;

  int magic;
  ifs.read((char*)&magic, sizeof(int));
  if(magic != MAGIC_NUMBER){
    set_conv_endian_true();
    cerr << magic << endl;
    magic = reverse_endian(magic);
    if(magic != MAGIC_NUMBER){
      cerr << magic << endl;
      cerr << "ERROR: " << filename << " : the first 4 bytes were not " << MAGIC_NUMBER << endl;
      exit(1);
    }
  }

  read_bin_values(&mmsys.launchset_version, 1);
  cout << "---- Input file format : version " << mmsys.launchset_version << endl;

  read_bin_values(&size_box, 1);
  read_bin_values(&size_crd, 1);
  read_bin_values(&size_vel, 1);
  read_bin_values(&size_topol, 1);
  read_bin_values(&size_constraint, 1);
  read_bin_values(&size_expand, 1);
  read_bin_values(&size_groups, 1);
  //int size_pcluster;
  //read_bin_values(&size_pcluster, 1);

  if(DBG==1){
    cout << "size_box:        " << size_box << endl;
    cout << "size_crd:        " << size_crd << endl;
    cout << "size_vel:        " << size_vel << endl;
    cout << "size_topol:      " << size_topol << endl;
    cout << "size_constraint: " << size_constraint << endl;
    cout << "size_expand:     " << size_expand << endl;
    cout << "size_groups:     " << size_groups << endl;
    //cout << "size_pcluster: " << size_pcluster << endl;
  }

  return 0;
};

int Read::load_ls_box(MmSystem& mmsys){
  // BOX
  int size_box;
  //read_bin_values(&size_box, 1);
  double pbc_val[12];
  cout << "load_ls_box : ";
  for(int i=0; i<12; i++){
    read_bin_values(&pbc_val[i], 1);
    mmsys.pbc_val[i] = (real)pbc_val[i];
    cout << pbc_val[i] << " " ;
  }
  mmsys.pbc.set_pbc(mmsys.pbc_val);
  return 0;
}

int Read::load_ls_crd(MmSystem& mmsys){
  // COORDINATES
  //int size_crd;
  //  read_bin_values(&size_crd, 1);
  //if(size_crd <= 0){
    // ERROR: size of coordinates data is zero
  //return 1;
  //}
  read_bin_values(&mmsys.n_atoms, 1);
  if(DBG==1){
    cout << "n_atoms: " << mmsys.n_atoms << endl;
  }
  if(mmsys.n_atoms <= 0){
    // ERROR: the number of atoms is zero
    return 1;
  }  
  mmsys.n_free = mmsys.n_atoms * 3;

  mmsys.alloc_atom_vars();

  for(int i=0; i < mmsys.n_atoms; i++){
    double x, y, z;
    read_bin_values(&x, 1);
    read_bin_values(&y, 1);
    read_bin_values(&z, 1);
    mmsys.crd[i][0] = (real)x;
    mmsys.crd[i][1] = (real)y;
    mmsys.crd[i][2] = (real)z;
    
    //cout << "load_ls_crd " << i << " : ";
    //cout << mmsys.crd[i][0] << " ";
    //cout << mmsys.crd[i][1] << " ";
    //cout << mmsys.crd[i][2] << endl;
  }
  return 0;
}
int Read::load_ls_vel(MmSystem& mmsys){
  // VELOCITIES
  //int size_crd;
  //read_bin_values(&size_crd, 1);
  //if(size_crd <= 0){
    // ERROR: size of coordinates data is zero
  //return 1;
  //}
  int n_atoms;
  read_bin_values(&n_atoms, 1);

  if(DBG==1){
    cout << "read_ls_vel ... n_atoms: " << n_atoms << endl;
  }
  if(mmsys.n_atoms != n_atoms){
    cerr << "ERROR" << endl;
    // ERROR: the number of atoms is inconsistent to the coordinates field
    return 2;
  }

  for(int i=0; i < mmsys.n_atoms; i++){
    double x, y, z;
    read_bin_values(&x, 1);
    read_bin_values(&y, 1);
    read_bin_values(&z, 1);
    mmsys.vel_just[i][0] = (real)x;
    mmsys.vel_just[i][1] = (real)y;
    mmsys.vel_just[i][2] = (real)z;
    /*
    cout << "load_ls_vel " << i << " : ";
    cout << mmsys.vel_next[i][0] << " ";
    cout << mmsys.vel_next[i][1] << " ";
    cout << mmsys.vel_next[i][2] << endl;
    */
  }
  return 0;
}

int Read::load_ls_tpl(MmSystem& mmsys){
  // TOPOLOGIES
  //int size_tpl;
  //read_bin_values(&size_tpl, 1);
  //if(size_tpl <= 0){
    // ERROR: size of topology data is zero
  //return 1;
  //}  
  int n_atoms;
  read_bin_values(&n_atoms, 1);
  if(mmsys.n_atoms != n_atoms){
    // ERROR: the number of atoms is inconsistent to the coordinates field
    return 2;
  }  
  // charge
  for(int i=0; i < mmsys.n_atoms; i++){
    double charge;
    read_bin_values(&charge, 1);
    mmsys.charge[i] = (real)charge;
    //if (DBG==1){ cout << "charge " << i << " : " << mmsys.charge[i] << endl; }
  }
    // mass
  for(int i=0; i < mmsys.n_atoms; i++){  
    double mass;
    read_bin_values(&mass, 1);
    mmsys.mass[i] = (real)mass;
    //if (DBG==1){ cout << "mass " << i << " : " << mmsys.mass[i] << endl; }
  }
	// atom_type
  for(int i=0; i < mmsys.n_atoms; i++){
    int atom_type;
    read_bin_values(&atom_type, 1);
    mmsys.atom_type[i] = atom_type-1;
    //    if (DBG>=1){ cout << "atom_type " << i << " : " << mmsys.atom_type[i] << endl; }
  }  

  // nbpair
  int size_lj;
  read_bin_values(&size_lj, 1);
  read_bin_values(&mmsys.n_lj_types, 1);
  read_bin_values(&mmsys.n_lj_type_pairs, 1);
  //if(DBG==1)
  //cout << "ljpair " << mmsys.n_lj_types << " " << mmsys.n_lj_type_pairs << endl;

  mmsys.alloc_lj_params();

  for(int i=0; i < mmsys.n_lj_type_pairs; i++){
    int type1, type2;
    double lj6, lj12;
    read_bin_values(&type1, 1);
    read_bin_values(&type2, 1);
    read_bin_values(&lj6, 1);
    read_bin_values(&lj12, 1);
    //if(DBG==1)
    //cout << "nbpair " << type1 << " " << type2 << " " << lj6 << " " << lj12 << endl;
    mmsys.set_lj_pair_param(type1-1, type2-1, (real_pw)lj6, (real_pw)lj12);
  }

  // bond
  int size_bond;
  read_bin_values(&size_bond, 1);
  read_bin_values(&mmsys.n_bonds, 1);
  if(DBG==1)
    cout << "n_bonds " << mmsys.n_bonds << endl;

  mmsys.alloc_bonds();
  for(int i=0; i < mmsys.n_bonds; i++){
    int atomid1, atomid2;
    double eps, r0;
    read_bin_values(&atomid1, 1);
    read_bin_values(&atomid2, 1);
    read_bin_values(&eps, 1);
    read_bin_values(&r0, 1);
    mmsys.set_bond_param(i, atomid1, atomid2, (real)eps, (real)r0);
    //if (DBG>=1){ cout << "bond " << i << " : " << atomid1 << " - " << atomid2 << endl;}
  }  
  
  // angle
  int size_angle;
  read_bin_values(&size_angle, 1);
  read_bin_values(&mmsys.n_angles, 1);
  if(DBG==1)
    cout << "n_angles " << mmsys.n_angles << endl;

  mmsys.alloc_angles();
  for(int i=0; i < mmsys.n_angles; i++){
    int atomid1, atomid2, atomid3;
    double eps, theta0;
    read_bin_values(&atomid1, 1);
    read_bin_values(&atomid2, 1);
    read_bin_values(&atomid3, 1);
    read_bin_values(&eps, 1);
    read_bin_values(&theta0, 1);
    mmsys.set_angle_param(i, atomid1, atomid2, atomid3, (real)eps, (real)theta0);
  }

  // torsion
  int size_torsion;
  read_bin_values(&size_torsion, 1);
  read_bin_values(&mmsys.n_torsions, 1);
  if(DBG==1)
    cout << "n_torsions " << mmsys.n_torsions << endl;

  mmsys.alloc_torsions();
  for(int i=0; i < mmsys.n_torsions; i++){
    int atomid1, atomid2, atomid3, atomid4;
    double ene, phase;
    int overlaps, symmetry, flag_14nb;
    read_bin_values(&atomid1, 1);
    read_bin_values(&atomid2, 1);
    read_bin_values(&atomid3, 1);
    read_bin_values(&atomid4, 1);
    read_bin_values(&ene, 1);
    read_bin_values(&overlaps, 1);
    read_bin_values(&symmetry, 1);
    read_bin_values(&phase, 1);
    read_bin_values(&flag_14nb, 1);
    mmsys.set_torsion_param(i, atomid1, atomid2, atomid3, atomid4, 
			    (real)ene, overlaps, symmetry,
			    (real)phase, flag_14nb);
    //cout << "torsion: " << atomid1 << "-" << atomid4 << " " << flag_14nb << endl;
  }

  // improper
  int size_impro;
  read_bin_values(&size_impro, 1);
  read_bin_values(&mmsys.n_impros, 1);
  if(DBG==1)
    cout << "n_impros " << mmsys.n_impros << endl;

  mmsys.alloc_impros();
  for(int i=0; i < mmsys.n_impros; i++){
    int atomid1, atomid2, atomid3, atomid4;
    double ene, phase;
    int overlaps, symmetry, flag_14nb;
    read_bin_values(&atomid1, 1);
    read_bin_values(&atomid2, 1);
    read_bin_values(&atomid3, 1);
    read_bin_values(&atomid4, 1);
    read_bin_values(&ene, 1);
    read_bin_values(&overlaps, 1);
    read_bin_values(&symmetry, 1);
    read_bin_values(&phase, 1);
    read_bin_values(&flag_14nb, 1);
    mmsys.set_impro_param(i, atomid1, atomid2, atomid3, atomid4, 
			    (real)ene, overlaps, symmetry,
			    (real)phase, flag_14nb);
  }

  // 14 nonbond
  int size_nb14;
  read_bin_values(&size_nb14, 1);
  read_bin_values(&mmsys.n_nb14, 1);
  if(DBG==1)
    cout << "n_nb14 " << mmsys.n_nb14 << endl;

  mmsys.alloc_nb14();
  for(int i=0; i < mmsys.n_nb14; i++){
    int atomid1, atomid2;
    int atomtype1, atomtype2;
    double coeff_vdw, coeff_ele;
    read_bin_values(&atomid1, 1);
    read_bin_values(&atomid2, 1);
    read_bin_values(&atomtype1, 1);
    read_bin_values(&atomtype2, 1);
    read_bin_values(&coeff_vdw, 1);
    read_bin_values(&coeff_ele, 1);
    mmsys.set_nb14_param(i, atomid1, atomid2,
			 atomtype1-1, atomtype2-1, 
			 (real)coeff_ele, (real)coeff_vdw);
    //cout << "nb14:" << atomid1 <<"-"<< atomid2 << " " << coeff_vdw << " " << coeff_ele << endl;
  }

  // without 15
  int size_nb15off;
  int n_nb15off;
  read_bin_values(&size_nb15off, 1);
  read_bin_values(&n_nb15off, 1);
  if(DBG==1)
    cout << "n_nb15off " << n_nb15off << endl;
  mmsys.alloc_nb15off();
  for(int i=0; i < n_nb15off; i++){
    int atomid1, atomid2;
    read_bin_values(&atomid1, 1);
    read_bin_values(&atomid2, 1);
    mmsys.set_nb15off(atomid1, atomid2);
    //mmsys.set_nb15off(atomid2, atomid1);
  }  
  mmsys.alloc_excess_pairs();
  mmsys.set_excess_pairs();
  return 0;
}

int Read::load_ls_constraint(Constraint* cst){
  int n_const_2;
  int n_const_3;
  int n_const_4;

  int atom[4];
  float dist[6];

  read_bin_values(&n_const_2, 1);
  read_bin_values(&n_const_3, 1);
  read_bin_values(&n_const_4, 1);
  
  cst->set_max_n_constraints(n_const_2, n_const_3, n_const_4);
  cst->alloc_constraint();

  int atomid1, atomid2, atomid3, atomid4;
  double dist1, dist2, dist3, dist4, dist5, dist6;
  // 2 atoms
  for(int i=0; i < n_const_2; i++){
    read_bin_values(&atomid1, 1);    
    read_bin_values(&atomid2, 1);    
    read_bin_values(&dist1, 1);    
    cst->add_pair(atomid1, atomid2, (real)dist1);
  }
  // 3 atoms
  for(int i=0; i < n_const_3; i++){
    read_bin_values(&atomid1, 1);    
    read_bin_values(&atomid2, 1);    
    read_bin_values(&atomid3, 1);    
    read_bin_values(&dist1, 1);    
    read_bin_values(&dist2, 1);    
    read_bin_values(&dist3, 1);    
    cst->add_trio(atomid1, atomid2, atomid3,
		  (real)dist1, (real)dist2, (real)dist3);
  }
  // 4 atoms
  for(int i=0; i < n_const_4; i++){
    read_bin_values(&atomid1, 1);    
    read_bin_values(&atomid2, 1);    
    read_bin_values(&atomid3, 1);    
    read_bin_values(&atomid4, 1);    
    read_bin_values(&dist1, 1);    
    read_bin_values(&dist2, 1);    
    read_bin_values(&dist3, 1);    
    read_bin_values(&dist4, 1);    
    read_bin_values(&dist5, 1);    
    read_bin_values(&dist6, 1);    
    cst->add_quad(atomid1, atomid2, atomid3, atomid4,
		  (real)dist1, (real)dist2, (real)dist3,
		 (real)dist4, (real)dist5, (real)dist6);
  }

  return 0;
}
int Read::load_ls_vmcmd(ExpandVMcMD* vmcmd){
  int n_vs;
  read_bin_values(&n_vs, 1);
  int interval;
  read_bin_values(&interval, 1);
  double temperature;
  read_bin_values(&temperature, 1);

  vmcmd->set_n_vstates(n_vs);
  vmcmd->set_trans_interval(interval);
  vmcmd->set_temperature((real)temperature);

  for(int i = 0; i < n_vs; i++){
    int ord;
    read_bin_values(&ord, 1);
    vmcmd->set_vs_order(i, ord);
    double lambda_low, lambda_high;
    double prob_low, prob_high;
    read_bin_values(&lambda_low, 1);
    read_bin_values(&lambda_high, 1);
    read_bin_values(&prob_low, 1);
    read_bin_values(&prob_high, 1);
    for(int j = 0; j < ord+1; j++){
      double buf;
      read_bin_values(&buf, 1);
      vmcmd->set_vs_poly_param(i, j, (real)buf);
    }
    double alpha_low, alpha_high;
    read_bin_values(&alpha_low, 1);
    read_bin_values(&alpha_high, 1);
    vmcmd->set_vs_params(i,
			 (real)lambda_low, (real)lambda_high,
			 (real)prob_low, (real)prob_high,
			 (real)alpha_low, (real)alpha_high);
  }
  int init, seed;
  read_bin_values(&init, 1);
  read_bin_values(&seed, 1);
  vmcmd->set_init_vs(init);
  vmcmd->set_random_seed(seed);
  return 0;
}
int Read::load_ls_atom_groups(MmSystem& mmsys){
  int n_groups;
  int* n_atoms_in_group;
  read_bin_values(&n_groups, 1);
  n_atoms_in_group = new int[n_groups];

  for(int i=0; i < n_groups; i++){
    int len_name;
    char name[MAX_LEN_NAME];
    read_bin_values(&len_name, 1);
    ifs.read(name, len_name);
    int n_atoms;
    read_bin_values(&n_atoms_in_group[i], 1);
  }
  mmsys.alloc_atom_groups(n_groups, n_atoms_in_group);
  
  for(int i=0; i < n_groups; i++){  
    for(int j = 0; j < n_atoms_in_group[i]; j++){
      read_bin_values(&mmsys.atom_groups[i][j], 1);      
    }
  }  
  delete[] n_atoms_in_group;
  return 0;
}

/*
int Read::load_ls_pcluster(MmSystem& mmsys){
  int n_clusters1;
  int n_atoms;
  int n_distances;
  int tmp_int;
  double tmp_dbl;
  read_bin_values(&n_clusters1, 1);
  mmsys.n_pclusters = n_clusters1 - 1;
  read_bin_values(&n_atoms, 1);
  mmsys.alloc_pcluster_vars();

  for (int i=0; i < n_clusters1; i++){
    read_bin_values(&mmsys.pcluster_index[i], 1);
    mmsys.n_atoms_pclusters = 0;
  }
  for (int i=0; i < n_atoms; i++){
    read_bin_values(&mmsys.pcluster_atoms[i], 1);
  }
  int pclst_id = 0;
  for(int idx=0; idx < n_atoms; idx++){
    if (mmsys.pcluster_index[pclst_id+1] == idx){
      pclst_id += 1;
      mmsys.pcluster_head_atom[pclst_id] = idx;
    }
    mmsys.atom_pclusters[idx] = pclst_id;
    mmsys.n_atoms_pclusters[pclst_id] += 1;
  }
  return 0;
}
*/
template <typename TYPE> int Read::read_bin_values(TYPE *recept, int len){
  ifs.read((char*)recept, sizeof(TYPE)*len);
  if(is_conv_endian()){
    for(int i=0;i<len;i++){
      recept[i] = reverse_endian(recept[i]);
    }
  }
  return 0;
}

