#include "Extend.h"

Extended::Extended()
  : CelesteObject(){
  write_lambda_interval = 0;
}

Extended::~Extended(){
}

void Extended::set_lambda_interval(int in_lambda_interval){
  write_lambda_interval = in_lambda_interval;
}


VirtualState::VirtualState()
  : CelesteObject(){
  poly_order = 0;
}

VirtualState::~VirtualState(){
  if(poly_order > 0)
    delete[] poly_params;
}
int VirtualState::set_order(int in_order){
  poly_order = in_order;
  poly_params = new real[poly_order+1];
  return 0;
}

int VirtualState::set_poly_param(int ord, real param){
  poly_params[ord] = param;
  return ord;
}
int VirtualState::set_params(real in_lambda_low,
			     real in_lambda_high,
			     real in_prob_low,
			     real in_prob_high){
  lambda_range[0] = in_lambda_low;
  lambda_range[1] = in_lambda_high;
  trans_prob[0] = in_prob_low;
  trans_prob[1] = in_prob_high;
  return 0;
}
int VirtualState::set_alpha(real in_alpha_low, real in_alpha_high){
  alpha[0] = in_alpha_low;
  alpha[1] = in_alpha_high;
  return 0;
}
bool VirtualState::is_in_range(real lambda){
  return (lambda >= lambda_range[0] and lambda <= lambda_range[1]);
}


ExtendedVMcMD::ExtendedVMcMD()
  : Extended(){
  n_vstates = 0;
  n_enhance_groups=0;
  //flg_vs_transition = false;
  flg_vs_transition = true;
}

ExtendedVMcMD::~ExtendedVMcMD(){
  if(n_vstates > 0)
    delete[] vstates;
  delete writer_lambda;
  for(int i=0; i < n_enhance_group_pairs; i++){
    delete[] enhance_group_pairs[i];
  }
  delete[] enhance_group_pairs;
}

int ExtendedVMcMD::set_n_vstates(int in_n_vstates){
  n_vstates = in_n_vstates;
  vstates = new VirtualState[n_vstates];
  return 0;
}

void ExtendedVMcMD::set_trans_interval(int in_trans_interval){
  trans_interval = in_trans_interval;
}

void ExtendedVMcMD::set_temperature(real in_tmp){
  temperature = in_tmp;
  const_k = (GAS_CONST / JOULE_CAL) * 1e-3 * temperature;
}

int ExtendedVMcMD::get_trans_interval(){
  return trans_interval;
}
int ExtendedVMcMD::get_temperature(){
  return temperature;
}

int ExtendedVMcMD::apply_bias(unsigned long cur_step,
			    real in_lambda,
			    real_fc* work,
			    int n_atoms_box){
  if(cur_step%trans_interval == 0){
    if(flg_vs_transition) set_current_vstate(in_lambda);
    write_vslog(cur_step);
  }
  scale_force(in_lambda, work, n_atoms_box);
  if(cur_step%write_lambda_interval == 0){
    write_lambda(in_lambda);
  }
  return 0;
}
int ExtendedVMcMD::set_current_vstate(real lambda){
  int dest_vs;
  dest_vs = trial_transition(cur_vs, 1, lambda);
  if( dest_vs != cur_vs ){
    cur_vs = dest_vs;
    return 0;
  }
  dest_vs = trial_transition(cur_vs, -1, lambda);
  if( dest_vs != cur_vs ){
    cur_vs = dest_vs;
    return 0;
  }
  return 0;
}
int ExtendedVMcMD::trial_transition(int source, int rel_dest,
				  real lambda){
  // source ... vs_id of current state
  // rel_dest ... -1 or 1, down or up
  // lambda

  // return ...
  if (source == 0 and rel_dest==-1)          return source;
  if (source == n_vstates-1 and rel_dest==1) return source;
  int up_down = rel_dest;
  if(rel_dest == -1) up_down = 0;
  if(vstates[source+rel_dest].is_in_range(lambda)){
    real dice = 1.0; // random value
    if(dice > (1.0 - vstates[source].get_trans_prob(up_down))){
      return source + rel_dest;
    }
  }
  return source;
}
int ExtendedVMcMD::scale_force(real lambda, real_fc* work, int n_atoms){
  
  // case 1 : under the lower limit
  real param = lambda;
  if (lambda <= vstates[cur_vs].get_lambda_low()){
    param = vstates[cur_vs].get_lambda_low();
  }else if (lambda >= vstates[cur_vs].get_lambda_high()){
    param = vstates[cur_vs].get_lambda_high();
  }
  real d_ln_p = vstates[cur_vs].get_poly_param(0);
  //cout << "dbg0522 1 " << param << " " << d_ln_p << endl;
  real tmp = 1.0;
  for(int i=1; i < vstates[cur_vs].get_order()+1; i++){
    tmp *= param;
    d_ln_p += vstates[cur_vs].get_poly_param(i) * tmp;
    //cout << "dbg0522 2 "<<i << " " << vstates[cur_vs].get_poly_param(i) << " " << d_ln_p << endl;
  }
  
  //real k = (GAS_CONST / JOULE_CAL) * 1e-3;
  real dew = const_k * d_ln_p;
  //cout << "dbg0522 "<<dew << endl;
  int n_atoms_3 = n_atoms * 3;
  for(int i = 0; i < n_atoms_3; i++){
    work[i] *= dew;
  }
  
  return 0;
}


int ExtendedVMcMD::set_files(string fn_vslog, string fn_lambda, int format_lambda){
  writer_vslog.set_fn(fn_vslog);
  writer_vslog.open();
  if (format_lambda == LAMBDAOUT_BIN){
    writer_lambda = new WriteTableLogBinary();
  }else if(format_lambda == LAMBDAOUT_ASC){
    writer_lambda = new WriteTableLogAscii();
  }else{
    writer_lambda = new WriteTableLog();
  }
  writer_lambda->set_fn(fn_lambda);
  writer_lambda->open();
  writer_lambda->set_ncolumns(1);
  writer_lambda->write_header();
  write_vslog(0);
  return 0;
}
int ExtendedVMcMD::close_files(){
  writer_vslog.close();
  writer_lambda->close();
  return 0;
}
int ExtendedVMcMD::write_vslog(int cur_steps){
  writer_vslog.write_ttpvMcMDLog(cur_steps, cur_vs);
  return 0;
}
int ExtendedVMcMD::write_lambda(real lambda){
  writer_lambda->write_row(&lambda);
  return 0;
}
int ExtendedVMcMD::set_vs_order(int vs_id, int ord){
  return vstates[vs_id].set_order(ord);
}

int ExtendedVMcMD::set_vs_params(int vs_id,
			       real lambda_low, real lambda_high,
			       real prob_low, real prob_high,
			       real alpha_low, real alpha_high){
  vstates[vs_id].set_params(lambda_low, lambda_high,
			    prob_low, prob_high);
  vstates[vs_id].set_alpha(alpha_low, alpha_high);
  return 0;
}
int ExtendedVMcMD::set_vs_poly_param(int vs_id, int ord, real param){
  return vstates[vs_id].set_poly_param(ord, param);
}

int ExtendedVMcMD::print_info(){
  
  cout << "V-McMD parameters" << endl;
  for(int i=0; i < n_vstates; i++){
    cout << "  Virtual state: " << i+1 << " ... ";
    cout << vstates[i].get_lambda_low() << " ~ ";
    cout << vstates[i].get_lambda_high() << endl;
    for(int j=0; j < vstates[i].get_order()+1; j++){
      cout << "    " << j << ": " << vstates[i].get_poly_param(j) << endl;
    }
  }
  return 0;
}

real ExtendedVMcMD::cal_struct_parameters(real* crd, PBC* pbc){
  return 0.0;
}

int ExtendedVMcMD::set_enhance_groups(int* in_n_atoms_in_groups,
				    int** in_atom_groups,
				    int in_n_enhance_groups,
				    int* in_enhance_groups){
  n_atoms_in_groups = in_n_atoms_in_groups;
  atom_groups = in_atom_groups;
  n_enhance_groups = in_n_enhance_groups;
  enhance_groups = in_enhance_groups;

  n_enhance_group_pairs = (n_enhance_groups * (n_enhance_groups-1)) / 2;
  enhance_group_pairs = new int*[n_enhance_group_pairs];
  int i_pair =  0;
  for(int i=0; i < n_enhance_groups; i++){
    for(int j=i+1; i < n_enhance_groups; i++){
      enhance_group_pairs[i_pair] = new int[2];
      enhance_group_pairs[i_pair][0] = i;
      enhance_group_pairs[i_pair][1] = j;
      i_pair++;
    }
  }
  return 0;
}

int ExtendedVMcMD::set_mass(real_pw* in_mass, real_pw* in_mass_groups, real_pw* in_mass_groups_inv){
  mass = in_mass;

  mass_sum = 0.0;
  for ( int i_grp=0; i_grp < n_enhance_groups; i_grp++){
    int grp_id = enhance_groups[i_grp];
    for (int i_atm=1; i_atm < n_atoms_in_groups[grp_id]; i_atm++){
      int aid = atom_groups[grp_id][i_atm];
      mass_sum += mass[aid];
    }
  }
  mass_groups = in_mass_groups;
  mass_groups_inv = in_mass_groups_inv;

  return 0;
}
int ExtendedVMcMD::set_params(real in_sigma){
  sigma = in_sigma;
  sigma_half = sigma * 0.5;  
  sigma_sq_inv = 1.0 / (sigma * sigma);
  return 0;
}

///////////////// ExtendedVAUS //////////////////

ExtendedVAUS::ExtendedVAUS(){
  n_enhance_groups=0;
}
ExtendedVAUS::~ExtendedVAUS(){
  free_crd_centers();
}
int ExtendedVAUS::alloc_crd_centers(){
  crd_centers = new real*[n_enhance_groups];
  for(int i=0; i < n_enhance_groups; i++){
    crd_centers[i] = new real[3];
  }
  unit_vec = new real*[n_enhance_group_pairs];
  for(int i=0; i < n_enhance_group_pairs; i++){
    unit_vec[i] = new real[3];
  }
  return 0;
}
int ExtendedVAUS::free_crd_centers(){
  for(int i=0; i < n_enhance_groups; i++){
    delete[] crd_centers[i];
  }
  delete[] crd_centers;
  for(int i=0; i < n_enhance_group_pairs; i++){ 
    delete[] unit_vec[i];
  }
  delete[] unit_vec;
  return 0;
}
real ExtendedVAUS::set_crd_centers(real* crd, PBC* pbc){
  for ( int i_grp=0; i_grp < n_enhance_groups; i_grp++){
    int grp_id = enhance_groups[i_grp];
    int aid0 = atom_groups[grp_id][0];
    int aid0_3 = aid0*3;
    for ( int d=0; d<3; d++)
      crd_centers[i_grp][d] = crd[aid0_3+d]*mass[aid0];
    for (int i_atm=1; i_atm < n_atoms_in_groups[grp_id]; i_atm++){
      int aid1 = atom_groups[grp_id][i_atm];
      //int aid1 = enhance_groups[i_grp][i_atm];
      int aid1_3 = aid1*3;
      for(int d=0; d<3; d++){
	real diff = crd[aid0_3+d] - crd[aid1_3+d];
	real mod_crd = crd[aid1_3+d];
	if(diff > pbc->L_half[d])       mod_crd += pbc->L[d];
	else if(-diff > pbc->L_half[d]) mod_crd -= pbc->L[d];
	crd_centers[i_grp][d] += mod_crd * mass[aid1];
      }
    }
    for(int d=0; d<3; d++) crd_centers[i_grp][d] *= mass_groups_inv[grp_id];
  }
  return 0;
}

real ExtendedVAUS::cal_struct_parameters(real* crd, PBC* pbc){
  // center of mass for each groups
  set_crd_centers(crd, pbc);
  //  real dist = 0.0;
  int i_pair=0;
  real lambda = 0.0;
  for (int i_grp = 0; i_grp < n_enhance_groups; i_grp++){
    for (int j_grp = i_grp+1; j_grp < n_enhance_groups; j_grp++){    
      real diff[3];
      pbc->diff_crd_minim_image(diff, crd_centers[i_grp], crd_centers[j_grp]);
      real dist = sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
      for(int d=0; d<3; d++) unit_vec[i_pair][d] = diff[d]/dist;
      lambda += dist;
      i_pair++;
    }
  }
  return lambda;
}

int ExtendedVAUS::scale_force(real lambda, real_fc* work, int n_atoms){
  
  // case 1 : under the lower limit
  real param = lambda;
  real recovery = 0.0;
  
  if (param <= vstates[cur_vs].get_lambda_low()){
    recovery = ( param - vstates[cur_vs].get_lambda_low() 
		       - sigma);
    param = vstates[cur_vs].get_lambda_low();
      //reparam = vstates[cur_vs].get_lambda_low();
  }else if(param >= vstates[cur_vs].get_lambda_high()){
    recovery = ( param - vstates[cur_vs].get_lambda_high() 
		 + sigma);
    param = vstates[cur_vs].get_lambda_high();
  }
  
  //cout << " param " << param << endl;
  //cout << "dbg0522 1 " << param << " recov: " << recovery << endl;
  real tmp_lambda = 1.0;
  real d_ln_p = vstates[cur_vs].get_poly_param(0);
  for(int i=1; i < vstates[cur_vs].get_order()+1; i++){
    tmp_lambda *= param;
    d_ln_p += vstates[cur_vs].get_poly_param(i) * tmp_lambda;
  }
  
  //real k = (GAS_CONST / JOULE_CAL) * 1e-3;
  real dew = 2.0 * const_k * ( d_ln_p + recovery );
  //cout << "dbg0522 "<<dew << endl;
  int n_atoms_3 = n_atoms * 3;
  for(int i_pair = 0; i_pair < n_enhance_group_pairs; i_pair++){
    real direction = -1.0;
    for(int pair_ab=0; pair_ab < 2; pair_ab++){
      int i_grp = enhance_group_pairs[i_pair][pair_ab];
      int grp_id = enhance_groups[i_grp];
      //cout << "pair : " << i_pair
      //<< " " << enhance_group_pairs[i_pair][0]
      //<< "-" << enhance_group_pairs[i_pair][1]
      //<< " grp_id: " << grp_id 
      //<< " dew: " << dew
      //<< " unit: " << unit_vec[i_pair][0]
      //<< " " << unit_vec[i_pair][1]
      //<< " " << unit_vec[i_pair][2] << endl;
      real bias[3];
      for(int d=0; d<3; d++)
	bias[d] = dew * (real)mass_groups_inv[grp_id] 
	  * unit_vec[i_pair][d];
      
      for(int i_at = 0; i_at < n_atoms_in_groups[grp_id]; i_at++){
	int atom_id3 = atom_groups[grp_id][i_at] * 3;
	for(int d=0; d<3; d++){
	  work[atom_id3+d] += direction * bias[d] * (real)mass[atom_groups[grp_id][i_at]];
	  
	}
	//cout << " bias : " << atom_groups[grp_id][i_at]  << " " 
	//<< direction << " " << grp_id << " "
	//<< bias[0] <<" " << bias[1] << " " << bias[2] << endl;

      }
      direction *= -1.0;
    }
    /*
    i_grp = enhance_group_pairs[i_pair][1];
    grp_id = enhance_groups[i_grp];
    for(int d=0; d<3; d++)
      bias[d] = dew / (real)n_atoms_in_groups[grp_id] * -unit_vec[i_pair][d];

    for(int i_at = 0; i_at < n_atoms_in_groups[grp_id]; i_at++){
      int atom_id3 = atom_groups[grp_id][i_at] * 3;
      for(int d=0; d<3; d++)
	work[atom_id3+d] += bias[d];
      
	}
    */
  }
  return 0;
}
int ExtendedVAUS::set_enhance_groups(int* in_n_atoms_in_groups,
				    int** in_atom_groups,
				    int in_n_enhance_groups,
				    int* in_enhance_groups){
  ExtendedVMcMD::set_enhance_groups(in_n_atoms_in_groups,
				  in_atom_groups,
				  in_n_enhance_groups,
				  in_enhance_groups);
  alloc_crd_centers();
  return 0;
}

