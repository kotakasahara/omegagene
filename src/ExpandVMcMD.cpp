#include "ExpandVMcMD.h"

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
  poly_params = new real[poly_order];
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


ExpandVMcMD::ExpandVMcMD()
  : Expand(){
  n_vstates = 0;
  flg_vs_transition = false;
}

ExpandVMcMD::~ExpandVMcMD(){
  if(n_vstates > 0)
    delete[] vstates;
}

int ExpandVMcMD::set_n_vstates(int in_n_vstates){
  n_vstates = in_n_vstates;
  vstates = new VirtualState[n_vstates];
  return 0;
}

void ExpandVMcMD::set_trans_interval(int in_trans_interval){
  trans_interval = in_trans_interval;
}

void ExpandVMcMD::set_temperature(real in_tmp){
  temperature = in_tmp;
  const_k = (GAS_CONST / JOULE_CAL) * 1e-3 * temperature;
}

int ExpandVMcMD::get_trans_interval(){
  return trans_interval;
}
int ExpandVMcMD::get_temperature(){
  return temperature;
}

int ExpandVMcMD::apply_bias(unsigned long cur_step,
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
int ExpandVMcMD::set_current_vstate(real lambda){
  int dest_vs = -1;
  int dest_vs1 = -1;
  int dest_vs2 = -1;
  if(cur_vs == 0){
    dest_vs = trial_transition(cur_vs, 1, lambda);
  }else if(cur_vs == n_vstates-1){
    dest_vs = trial_transition(cur_vs, -1, lambda);
  }else{
    dest_vs1 = trial_transition(cur_vs, 1, lambda);
    dest_vs2 = trial_transition(cur_vs, -1, lambda);
  }
  if(dest_vs1 != dest_vs2){
    cerr << "Error: set_current_vstate" << endl;
    return 1;
  }
  return 0;
}
int ExpandVMcMD::trial_transition(int source, int rel_dest,
				  real lambda){
  // source ... vs_id of current state
  // rel_dest ... -1 or 1, down or up
  // lambda

  // return ...

  int up_down = rel_dest;
  if(rel_dest == -1) up_down = 0;

  if(vstates[source+rel_dest].is_in_range(lambda)){
    real dice = 1.0; // random value
    if(dice > (1.0 - vstates[source+rel_dest].get_trans_prob(up_down))){
      return source + rel_dest;
    }
  }
  return source;
}

int ExpandVMcMD::scale_force(real lambda, real_fc* work, int n_atoms){
  
  // case 1 : under the lower limit
  real param = lambda;
  if (lambda <= vstates[cur_vs].get_lambda_low()){
    param = vstates[cur_vs].get_lambda_low();
  }else if (lambda >= vstates[cur_vs].get_lambda_high()){
    param = vstates[cur_vs].get_lambda_high();
  }
  real d_ln_p = vstates[cur_vs].get_poly_param(0);
  for(int i=1; i < vstates[cur_vs].get_order()+1; i++){
    real tmp = pow(param, i);
    d_ln_p += vstates[cur_vs].get_poly_param(i) * tmp;
  }
  
  //real k = (GAS_CONST / JOULE_CAL) * 1e-3;
  real dew = const_k * d_ln_p;

  int n_atoms_3 = n_atoms * 3;
  for(int i = 0; i < n_atoms_3; i++){
    work[i] *= dew;
  }
  
  return 0;
}

int ExpandVMcMD::set_files(string fn_vslog, string fn_lambda, int format_lambda){
  writer_vslog.set_fn(fn_vslog);
  writer_vslog.open();
  if (format_lambda == LAMBDAOUT_BIN){
    writer_lambda = new WriteTableLogBinary();
  }else if(format_lambda == LAMBDAOUT_ASC){
    writer_lambda = new WriteTableLogAscii();
  }
  writer_lambda->set_fn(fn_lambda);
  writer_lambda->open();
  writer_lambda->set_ncolumns(1);
  writer_lambda->write_header();
  write_vslog(0);
  return 0;
}
int ExpandVMcMD::close_files(){
  writer_vslog.close();
  writer_lambda->close();
  return 0;
}
int ExpandVMcMD::write_vslog(int cur_steps){
  writer_vslog.write_ttpvMcMDLog(cur_steps, cur_vs);
  return 0;
}
int ExpandVMcMD::write_lambda(real lambda){
  writer_lambda->write_row(&lambda);
  return 0;
}
int ExpandVMcMD::set_vs_order(int vs_id, int ord){
  return vstates[vs_id].set_order(ord);
}

int ExpandVMcMD::set_vs_params(int vs_id,
			       real lambda_low, real lambda_high,
			       real prob_low, real prob_high,
			       real alpha_low, real alpha_high){
  vstates[vs_id].set_params(lambda_low, lambda_high,
			    prob_low, prob_high);
  vstates[vs_id].set_alpha(alpha_low, alpha_high);
  return 0;
}
int ExpandVMcMD::set_vs_poly_param(int vs_id, int ord, real param){
  return vstates[vs_id].set_poly_param(ord, param);
}
