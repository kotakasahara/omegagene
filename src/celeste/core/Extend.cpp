#include "Extend.h"
#include <ciso646>

using namespace std;
using namespace celeste;

Extended::Extended() : CelesteObject() {
    write_lambda_interval = 0;
}

Extended::~Extended() {}

void Extended::set_lambda_interval(int in_lambda_interval) {
    write_lambda_interval = in_lambda_interval;
}
void Extended::set_com_interval(int in_com_interval) {
    write_com_interval = in_com_interval;
}

VirtualState::VirtualState() : CelesteObject() {
    poly_order = 0;
}

VirtualState::~VirtualState() {
    if (poly_order > 0) delete[] poly_params;
}

int VirtualState::set_order(int in_order) {
    poly_order  = in_order;
    poly_params = new real[poly_order + 1];
    return 0;
}

int VirtualState::set_poly_param(int ord, real param) {
    poly_params[ord] = param;
    return ord;
}
int VirtualState::set_params(real in_lambda_low, real in_lambda_high, real in_prob_low, real in_prob_high) {
  lambda_range[0] = in_lambda_low;
  lambda_range[1] = in_lambda_high;
  trans_prob[0]   = in_prob_low;
  trans_prob[1]   = in_prob_high;
  // cout << "dbg1225 vs set_param "<< trans_prob[0] << " " << trans_prob[1] << endl;
  return 0;
}
int VirtualState::set_alpha(real in_alpha_low, real in_alpha_high) {
    alpha[0] = in_alpha_low;
    alpha[1] = in_alpha_high;
    return 0;
}
bool VirtualState::is_in_range(real lambda) {
    return (lambda >= lambda_range[0] and lambda <= lambda_range[1]);
}

ExtendedVMcMD::ExtendedVMcMD() : Extended() {
    n_vstates        = 0;
    n_enhance_groups = 0;
    // flg_vs_transition = false;
    flg_vs_transition = true;
}

ExtendedVMcMD::~ExtendedVMcMD() {
    if (n_vstates > 0) delete[] vstates;
    delete writer_lambda;
    for (int i = 0; i < n_enhance_group_pairs; i++) { delete[] enhance_group_pairs[i]; }
    delete[] enhance_group_pairs;
    free_crd_centers();
}

int ExtendedVMcMD::alloc_crd_centers() {
    crd_groups = new real **[n_enhance_groups];
    for (int i_grp = 0; i_grp < n_enhance_groups; i_grp++) {
        int grp_id        = enhance_groups[i_grp];
        crd_groups[i_grp] = new real *[n_atoms_in_groups[grp_id]];
        for (int j = 0; j < n_atoms_in_groups[grp_id]; j++) { crd_groups[i_grp][j] = new real[3]; }
    }

    crd_centers = new real *[n_enhance_groups];
    for (int i = 0; i < n_enhance_groups; i++) { crd_centers[i] = new real[3]; }
    unit_vec = new real *[n_enhance_group_pairs];
    for (int i = 0; i < n_enhance_group_pairs; i++) { unit_vec[i] = new real[3]; }
    return 0;
}
int ExtendedVMcMD::free_crd_centers() {
    for (int i = 0; i < n_enhance_groups; i++) {
        int grp_id = enhance_groups[i];
        for (int j = 0; j < n_atoms_in_groups[grp_id]; j++) { delete[] crd_groups[i][j]; }
        delete[] crd_groups[i];
    }
    delete[] crd_groups;

    for (int i = 0; i < n_enhance_groups; i++) { delete[] crd_centers[i]; }
    delete[] crd_centers;
    for (int i = 0; i < n_enhance_group_pairs; i++) { delete[] unit_vec[i]; }
    delete[] unit_vec;
    return 0;
}

int ExtendedVMcMD::set_n_vstates(int in_n_vstates) {
    n_vstates = in_n_vstates;
    vstates   = new VirtualState[n_vstates];
    return 0;
}


void ExtendedVMcMD::set_trans_interval(int in_trans_interval) {
    trans_interval = in_trans_interval;
}

void ExtendedVMcMD::set_temperature(real in_tmp) {
    temperature = in_tmp;
    const_k     = (GAS_CONST / JOULE_CAL) * 1e-3 * temperature;
}

int ExtendedVMcMD::get_trans_interval() {
    return trans_interval;
}
int ExtendedVMcMD::get_temperature() {
    return temperature;
}

int ExtendedVMcMD::apply_bias(unsigned long cur_step, real in_lambda, real_fc *work, int n_atoms_box) {
  //cout << "test in_lambda "<< in_lambda<<endl;
    if (cur_step > 0 && cur_step % trans_interval == 0) {
        if (flg_vs_transition) set_current_vstate(in_lambda);
        if (cur_step <= n_steps) { write_vslog(cur_step); }
    }
    //cout << " test 1 2 " << endl;
    scale_force(in_lambda, work, n_atoms_box);
    //cout << " test 1 3 " << endl;
    if (cur_step > 0 && cur_step % write_lambda_interval == 0 && cur_step <= n_steps) { write_lambda(in_lambda); }
    //cout << "dbg write_com_interval " << write_com_interval << endl;
    if (write_com_interval > 0 && cur_step > 0 &&
      cur_step % write_com_interval == 0 &&
      cur_step <= n_steps) { write_com(); }
  
    return 0;
}
int ExtendedVMcMD::set_current_vstate(real lambda) {
    int dest_vs;
    dest_vs = trial_transition(cur_vs, 1, lambda);
    if (dest_vs != cur_vs) {
        cur_vs = dest_vs;
        return 0;
    }
    dest_vs = trial_transition(cur_vs, -1, lambda);
    if (dest_vs != cur_vs) {
        cur_vs = dest_vs;
        return 0;
    }
    return 0;
}
int ExtendedVMcMD::trial_transition(int source, int rel_dest, real lambda) {
    // source ... vs_id of current state
    // rel_dest ... -1 or 1, down or up
    // lambda

    // return ...

    // std::uniform_real_distribution<> random_gen(0.0, 1.0);
    if (!vstates[source].is_in_range(lambda)) return source;
    if (source == 0 and rel_dest == -1) return source;
    if (source == n_vstates - 1 and rel_dest == 1) return source;
    int up_down                 = rel_dest;
    if (rel_dest == -1) up_down = 0;
    if (vstates[source + rel_dest].is_in_range(lambda)) {
        if ((*random_mt)() > (1.0 - vstates[source].get_trans_prob(up_down))) { return source + rel_dest; }
    }
    return source;
}
int ExtendedVMcMD::scale_force(real lambda, real_fc *work, int n_atoms) {

    // case 1 : under the lower limit
    real param = lambda;
    if (lambda <= vstates[cur_vs].get_lambda_low()) {
      param = vstates[cur_vs].get_lambda_low();
    } else if (lambda >= vstates[cur_vs].get_lambda_high()) {
      param = vstates[cur_vs].get_lambda_high();
    }
    //cout << "dbg0522" << endl;
    real d_ln_p = vstates[cur_vs].get_poly_param(0);
    //cout << "dbg0522 1 " << param << " " << d_ln_p << endl;
    real tmp = 1.0;
    for (int i = 1; i < vstates[cur_vs].get_order() + 1; i++) {
        tmp *= param;
        d_ln_p += vstates[cur_vs].get_poly_param(i) * tmp;
	//cout << "dbg0522 2 "<<i << " " << vstates[cur_vs].get_poly_param(i) << " " << d_ln_p << endl;
    }

    // real k = (GAS_CONST / JOULE_CAL) * 1e-3;
    real dew = const_k * d_ln_p;
    //cout << "dbg0522 "<<dew << endl;
    int n_atoms_3 = n_atoms * 3;
    for (int i = 0; i < n_atoms_3; i++) { work[i] *= dew; }

    return 0;
}

int ExtendedVMcMD::set_files(string fn_vslog, string fn_lambda, int format_lambda,
			     string fn_com) {
  //cout << "set_files" << endl;
    writer_vslog.set_fn(fn_vslog);
    writer_vslog.open();
    if (format_lambda == LAMBDAOUT_BIN) {
        writer_lambda = new WriteTableLogBinary();
    } else if (format_lambda == LAMBDAOUT_ASC) {
        writer_lambda = new WriteTableLogAscii();
    } else {
        writer_lambda = new WriteTableLog();
    }
    writer_lambda->set_fn(fn_lambda);
    writer_lambda->open();
    writer_lambda->set_ncolumns(1);
    writer_lambda->write_header();

    if (format_lambda == LAMBDAOUT_BIN) {
      writer_com = new WriteTableLogBinary();
    } else if (format_lambda == LAMBDAOUT_ASC) {
      writer_com = new WriteTableLogAscii();
    } else {
      writer_com = new WriteTableLog();
    }
    //cout << "test1" << endl;
    if(write_com_interval>0){
      //cout << "test2" << endl;
      writer_com->set_fn(fn_com);
      writer_com->open();
      writer_com->set_ncolumns(6);
      writer_com->write_header();
    }
    //cout << "test3" << endl;
    // write_vslog(0);
    return 0;
}
int ExtendedVMcMD::close_files() {
    writer_vslog.close();
    writer_lambda->close();
    if(writer_com->is_open())
      writer_com->close();
    return 0;
}
int ExtendedVMcMD::write_vslog(int cur_steps) {
    writer_vslog.write_ttpvMcMDLog(cur_steps, cur_vs);
    return 0;
}
int ExtendedVMcMD::write_lambda(real lambda) {
    writer_lambda->write_row(&lambda);
    return 0;
}
int ExtendedVMcMD::write_com() {
  real val[6] = {crd_centers[0][0],
		 crd_centers[0][1],
		 crd_centers[0][2],
		 crd_centers[1][0],
		 crd_centers[1][1],
		 crd_centers[1][2]};
  writer_com->write_row(val);
  return 0;
}
int ExtendedVMcMD::set_vs_order(int vs_id, int ord) {
    return vstates[vs_id].set_order(ord);
}

int ExtendedVMcMD::set_vs_params(int  vs_id,
                                 real lambda_low,
                                 real lambda_high,
                                 real prob_low,
                                 real prob_high,
                                 real alpha_low,
                                 real alpha_high) {
    vstates[vs_id].set_params(lambda_low, lambda_high, prob_low, prob_high);
    vstates[vs_id].set_alpha(alpha_low, alpha_high);
    return 0;
}
int ExtendedVMcMD::set_vs_poly_param(int vs_id, int ord, real param) {
    return vstates[vs_id].set_poly_param(ord, param);
}

int ExtendedVMcMD::print_info() {

    cout << "V-McMD parameters" << endl;
    for (int i = 0; i < n_vstates; i++) {
        cout << "  Virtual state: " << i + 1 << " ... ";
        cout << vstates[i].get_lambda_low() << " ~ ";
        cout << vstates[i].get_lambda_high() << " , ";
        cout << vstates[i].get_trans_prob(0) << " - ";
        cout << vstates[i].get_trans_prob(1) << endl;

        for (int j = 0; j < vstates[i].get_order() + 1; j++) {
            cout << "    " << j << ": " << vstates[i].get_poly_param(j) << endl;
        }
    }
    return 0;
}

real ExtendedVMcMD::cal_struct_parameters(real *crd, PBC *pbc) {
    return 0.0;
}

int ExtendedVMcMD::set_enhance_groups(int *       in_n_atoms_in_groups,
                                      int **      in_atom_groups,
                                      int         in_n_enhance_groups,
                                      std::vector<int> in_enhance_groups) {
    n_atoms_in_groups = in_n_atoms_in_groups;
    atom_groups       = in_atom_groups;
    n_enhance_groups  = in_n_enhance_groups;
    enhance_groups    = in_enhance_groups;

    n_enhance_group_pairs = (n_enhance_groups * (n_enhance_groups - 1)) / 2;
    enhance_group_pairs   = new int *[n_enhance_group_pairs];
    int i_pair            = 0;
    for (int i = 0; i < n_enhance_groups; i++) {
        for (int j = i + 1; i < n_enhance_groups; i++) {
            enhance_group_pairs[i_pair]    = new int[2];
            enhance_group_pairs[i_pair][0] = i;
            enhance_group_pairs[i_pair][1] = j;
            i_pair++;
        }
    }
    alloc_crd_centers();
    return 0;
}

int ExtendedVMcMD::set_mass(real_pw *in_mass, real_pw *in_mass_groups, real_pw *in_mass_groups_inv) {
    mass = in_mass;

    mass_groups     = in_mass_groups;
    mass_groups_inv = in_mass_groups_inv;
    mass_sum        = 0.0;
    for (int i_grp = 0; i_grp < n_enhance_groups; i_grp++) {
        int grp_id = enhance_groups[i_grp];
        for (int i_atm = 1; i_atm < n_atoms_in_groups[grp_id]; i_atm++) {
            int aid = atom_groups[grp_id][i_atm];
            mass_sum += mass[aid];
        }
        // cout << "dbg1130 mass " << i_grp << "-" << grp_id << " " << mass_groups_inv[grp_id] << endl;;
    }

    return 0;
}
int ExtendedVMcMD::set_params(random::Random *in_mt, real in_sigma, real in_recov_coef, int in_n_steps) {
    random_mt = in_mt;
    sigma     = in_sigma;
    // sigma_half = sigma * 0.5;
    // sigma_sq_inv = 1.0 / (sigma * sigma);
    recov_coef = in_recov_coef;
    n_steps    = in_n_steps;
    // aus_type = in_aus_type;
    return 0;
}

int ExtendedVMcMD::write_aus_restart(string fn_out) {
    WriteGroupCoord writer;
    writer.set_fn(fn_out);
    writer.open();
    writer.write_aus_restart(aus_type, n_enhance_groups, enhance_groups, n_atoms_in_groups, crd_groups);
    writer.close();
    return 0;
}
///////////////// ExtendedVAUS //////////////////

ExtendedVAUS::ExtendedVAUS() {
    n_enhance_groups = 0;
}
ExtendedVAUS::~ExtendedVAUS() {
    free_crd_centers();
}

real ExtendedVAUS::set_crd_centers(real *crd, PBC *pbc) {
    for (int i_grp = 0; i_grp < n_enhance_groups; i_grp++) {
      int grp_id = enhance_groups[i_grp];
      //int aid0   = atom_groups[grp_id][0];
      //int aid0_3 = aid0 * 3;
        // cout << "dbg1130 grp " << i_grp << " " << grp_id << " " << n_atoms_in_groups[grp_id] << endl;
        // real crd0[3] = {crd[aid0_3], crd[aid0_3 + 1], crd[aid0_3 + 2]};
      for (int d = 0; d < 3; d++) { crd_centers[i_grp][d] = 0.0; }
      for (int i_atm = 0; i_atm < n_atoms_in_groups[grp_id]; i_atm++) {
	int aid   = atom_groups[grp_id][i_atm];
	int aid_3 = aid * 3;
	for (int d = 0; d < 3; d++) {
	  real tmp_crd = crd[aid_3 + d];
	  real diff    = tmp_crd - crd_groups[i_grp][i_atm][d];
	  while (diff > pbc->L_half[d]) {
	    diff -= pbc->L[d];
	    tmp_crd -= pbc->L[d];
	  }
	  while (-diff > pbc->L_half[d]) {
	    diff += pbc->L[d];
	    tmp_crd += pbc->L[d];
	  }
	  crd_groups[i_grp][i_atm][d] = tmp_crd;
	  crd_centers[i_grp][d] += tmp_crd * mass[aid];
	}
	// cout << "dbg1130 crd " << crd_groups[i_grp][i_atm][0] << " "
	//<< crd_groups[i_grp][i_atm][1] << " "
	//<< crd_groups[i_grp][i_atm][2] << " "
	//<< endl;
      }
      for (int d = 0; d < 3; d++) { crd_centers[i_grp][d] *= mass_groups_inv[grp_id]; }
    }
    //  cout << "dbg1126 center " << crd_centers[0][0] << " "
    //<< crd_centers[0][1] << " "
    //<< crd_centers[0][2] << " "
    //<< crd_centers[1][0] << " "
    //<< crd_centers[1][1] << " "<< crd_centers[1][2] << endl;
    return 0;
}

int ExtendedVAUS::set_init_crd_groups(real *crd) {
    for (int i_grp = 0; i_grp < n_enhance_groups; i_grp++) {
        int grp_id = enhance_groups[i_grp];
        for (int i_atm_grp = 0; i_atm_grp < n_atoms_in_groups[grp_id]; i_atm_grp++) {
            int aid0   = atom_groups[grp_id][0];
            int aid0_3 = aid0 * 3;
            for (int d = 0; d < 3; d++) { crd_groups[i_grp][i_atm_grp][d] = crd[aid0_3 + d]; }
        }
    }
    return 0;
}

real ExtendedVAUS::cal_struct_parameters(real *crd, PBC *pbc) {
    // center of mass for each groups
    set_crd_centers(crd, pbc);
    //  real dist = 0.0;
    int  i_pair = 0;
    real lambda = 0.0;
    for (int i_grp = 0; i_grp < n_enhance_groups; i_grp++) {
        for (int j_grp = i_grp + 1; j_grp < n_enhance_groups; j_grp++) {
	  real diff[3];
	  pbc->diff_crd_minim_image(diff, crd_centers[i_grp], crd_centers[j_grp]);
	  real dist = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
	  for (int d = 0; d < 3; d++) unit_vec[i_pair][d] = diff[d] / dist;
	  lambda += dist;
	  i_pair++;
        }
    }
    return lambda;
}
int ExtendedVAUS::scale_force(real lambda, real_fc *work, int n_atoms) {
  
  // case 1 : under the lower limit
  real param    = lambda;
   real recovery = 0.0;

    if (param <= vstates[cur_vs].get_lambda_low()) {
        if (param <= vstates[cur_vs].get_lambda_low() - sigma)
            recovery = recov_coef * (param - (vstates[cur_vs].get_lambda_low() - sigma));
	//cout << "recov- " << param << endl;
        param        = vstates[cur_vs].get_lambda_low();
    } else if (param >= vstates[cur_vs].get_lambda_high()) {
        if (param >= vstates[cur_vs].get_lambda_high() + sigma)
            recovery = recov_coef * (param - (vstates[cur_vs].get_lambda_high() + sigma));
	//cout << "recov+ " << param << endl;
        param        = vstates[cur_vs].get_lambda_high();
    }
    // cout << " param " << param << endl;
    // cout << "dbg0522 1 " << param << " recov: " << recovery << endl;
    real tmp_lambda = 1.0;
    real d_ln_p     = vstates[cur_vs].get_poly_param(0);
    for (int i = 1; i <= vstates[cur_vs].get_order(); i++) {
        tmp_lambda *= param;
        d_ln_p += vstates[cur_vs].get_poly_param(i) * tmp_lambda;
    }

    // real k = (GAS_CONST / JOULE_CAL) * 1e-3;
    real dew = const_k * (d_ln_p + recovery);
    // cout << "dbg0522 "<<dew << endl;
    // int n_atoms_3 = n_atoms * 3;
    for (int i_pair = 0; i_pair < n_enhance_group_pairs; i_pair++) {
        real direction = 1.0;
        for (int pair_ab = 0; pair_ab < 2; pair_ab++) {
            int i_grp  = enhance_group_pairs[i_pair][pair_ab];
            int grp_id = enhance_groups[i_grp];
            /*
            cout << "lambda: " << lambda
          << " dlnp: " << d_ln_p
          << " recov: " << recovery
          <<" const_k: " << const_k
          << endl;
          cout   << " pair : " << i_pair
          << " " << enhance_group_pairs[i_pair][0]
          << "-" << enhance_group_pairs[i_pair][1]
          << " grp_id: " << grp_id
          << " dew: " << dew
          << " unit: " << unit_vec[i_pair][0]
          << " " << unit_vec[i_pair][1]
          << " " << unit_vec[i_pair][2] << endl;
            */
            real bias[3];
            for (int d = 0; d < 3; d++) bias[d] = dew * unit_vec[i_pair][d] * (real)mass_groups_inv[grp_id];
            /// n_atoms_in_groups[grp_id];

            for (int i_at = 0; i_at < n_atoms_in_groups[grp_id]; i_at++) {
                int atom_id3 = atom_groups[grp_id][i_at] * 3;
                /*
                  cout 	  << " bias : " << atom_groups[grp_id][i_at] << " at:" << i_at
                  << " ("<< atom_groups[grp_id][i_at] * 3  << ")"
                  << " mass:"<< (real)mass[atom_groups[grp_id][i_at]]
                  << direction << " " << grp_id << " "
                  <<endl;
                  cout << "prev:   " << work[atom_id3+0] << " " << work[atom_id3+1] << " "
                  << work[atom_id3+2] << endl;
                */
                for (int d = 0; d < 3; d++) {
                    work[atom_id3 + d] += direction * bias[d] * (real)mass[atom_groups[grp_id][i_at]];
                }

		//cout << "biased: " << work[atom_id3+0] << " " << work[atom_id3+1] << " "
		//<< work[atom_id3+2] << endl;
		//cout << "bias:   " << direction * bias[0] * (real)mass[atom_groups[grp_id][i_at]] << " "
		//<< direction * bias[1] * (real)mass[atom_groups[grp_id][i_at]] << " "
		//<< direction * bias[2] * (real)mass[atom_groups[grp_id][i_at]] <<endl;

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

///////////// VcMD //////////////

ExtendedVcMD::ExtendedVcMD() : Extended() {
    flg_vs_transition = true;
    q_raw_sum = 0;
}

ExtendedVcMD::~ExtendedVcMD() {
    delete writer_lambda;
    free_crd_centers();
}
int ExtendedVcMD::set_default_q_cano(){
  std::vector<int> key;
  for(int i=0; i < n_dim; i++){
    key.push_back(-1);
  }
  if(q_cano.find(key) == q_cano.end()){
    real min_q_cano = 1e10;
    for(auto i_q : q_cano){
      if (i_q.second > 0 && i_q.second < min_q_cano)
	min_q_cano = i_q.second;
    }
    if (min_q_cano == 1e10) min_q_cano=1;
    default_q_cano = min_q_cano;
  }else{
    default_q_cano = q_cano[key];
  }
  cout << "dbg 039j default_q_cano = "  << default_q_cano << endl;
  return 0;
}

void ExtendedVcMD::set_trans_interval(int in_trans_interval) {
    trans_interval = in_trans_interval;
}
void ExtendedVcMD::set_temperature(real in_tmp) {
    temperature = in_tmp;
    const_k     = (GAS_CONST / JOULE_CAL) * 1e-3 * temperature;
}
int ExtendedVcMD::set_params(random::Random *in_mt, real in_sigma,
			     real in_recov_coef, int in_n_steps,
			     int in_begin_count_q_raw,
			     int in_drift) {
  lambda.resize(n_dim);
  random_mt = in_mt;
  sigma     = in_sigma;
  // sigma_half = sigma * 0.5;
  // sigma_sq_inv = 1.0 / (sigma * sigma);
  recov_coef = in_recov_coef;
  n_steps    = in_n_steps;
  // aus_type = in_aus_type;
  begin_count_q_raw = in_begin_count_q_raw;
  drift = in_drift;
  return 0;
}
void ExtendedVcMD::set_n_dim(int in_n_dim){
  n_dim = in_n_dim;
  lambda.resize(n_dim);
  vs_next.resize(pow(2, n_dim));
  for ( auto vs1 : vs_next ) {
    vs1.resize(n_dim);
  }
}

int ExtendedVcMD::close_files() {
  write_q();
  writer_vslog.close();
  writer_lambda->close();
  return 0;
}
int ExtendedVcMD::push_vs_range(std::vector<real> new_min,
				std::vector<real> new_max){
  vc_range_min.push_back(new_min);
  vc_range_max.push_back(new_max);
  return 0;
}
void ExtendedVcMD::set_q_cano(std::map< std::vector<int>, real > in_q){
  q_cano = in_q; 
  for ( const auto itr : q_cano ) {
    q_raw[itr.first] = 0;
  }
}
int ExtendedVcMD::set_struct_parameters(real *crd, PBC *pbc) {
  //cout << "set_struct_parameters " << reactcrd_type << endl;

  if(reactcrd_type == AUSTYPE_MASSCENTER){
    set_struct_parameters_mass_center(crd, pbc);
  }else if(reactcrd_type == AUSTYPE_MIN){
    set_struct_parameters_mass_center(crd, pbc);
    set_struct_parameters_min(crd, pbc);    
  }
  return 0;
}
int ExtendedVcMD::set_struct_parameters_mass_center(real *crd, PBC *pbc) {
      //cout << "dbg 0303 set_struct_parameters" << endl;
      // center of mass for each groups
  set_crd_centers(crd, pbc);
  ///cout << "dbg 0303 set_crd_centers (finished)" << endl;
  //  real dist = 0.0;
  int  i_pair = 0;
  for(const auto itr_dim : enhance_group_pairs){
    real diff[3];
    int i_grp = itr_dim[0];
    int j_grp = itr_dim[1];
    //cout << "test i_grp: "<<i_grp<<" j_grp: " << j_grp << endl;
    pbc->diff_crd_minim_image(diff, crd_centers[i_grp], crd_centers[j_grp]);
    real dist = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
    for (int d = 0; d < 3; d++) unit_vec[i_pair][d] = diff[d] / dist;
    //cout << "dbg uv " << i_grp << " " << j_grp << " " << unit_vec[i_pair][0] << " "
    //<< unit_vec[i_pair][1] << " "
    //<< unit_vec[i_pair][2] <<endl;
    lambda[i_pair]=dist;
    //cout << lambda[i_pair] << endl;
    i_pair++;
  }
  //cout << "dbg 0303 set_struct_parameters (finished)" << endl;  
  return 0;
}
int ExtendedVcMD::set_struct_parameters_min(real *crd, PBC *pbc) {

  int  i_pair = 0;
  for(const auto itr_dim : enhance_group_pairs){
    int i_grp1 = itr_dim[0];
    int i_grp2 = itr_dim[1];
    int grp_id1 = enhance_groups[i_grp1];
    int grp_id2 = enhance_groups[i_grp2];
    real min_dist = 1e10;

    for (int i_atm1 = 0; i_atm1 < n_atoms_in_groups[grp_id1]; i_atm1++) {
      //int aid1 = atom_groups[grp_id1][i_atm1];
      int aid1_3 = atom_groups[grp_id1][i_atm1]*3;
      real crd1[3] = {crd[aid1_3], crd[aid1_3+1], crd[aid1_3+2]};
      for (int i_atm2 = 0; i_atm2 < n_atoms_in_groups[grp_id2]; i_atm2++) {
	//int aid2 = atom_groups[grp_id2][i_atm2];
	int aid2_3 = atom_groups[grp_id2][i_atm2]*3;
	real crd2[3] = {crd[aid2_3], crd[aid2_3+1], crd[aid2_3+2]};
	real d_crd[3];
	pbc->diff_crd_minim_image(d_crd, crd1, crd2);
	real dist = sqrt(d_crd[0]*d_crd[0]+d_crd[1]*d_crd[1]+d_crd[2]*d_crd[2]);
	if(dist <= min_dist) min_dist=dist;
      }
    }

    lambda[i_pair] = min_dist;
    cout << "lambda " << i_pair << " " << lambda[i_pair] << endl;
    i_pair++;
    //// unit_vec has been set in the set_struct_parameters_mass_center ////
  }

  return 0;
}

real ExtendedVcMD::set_crd_centers(real *crd, PBC *pbc) {
  for (int i_grp = 0; i_grp < n_enhance_groups; i_grp++) {
    // crd_centers[ grp id in the enhance groups ]

    int grp_id = enhance_groups[i_grp];
    //int aid0   = atom_groups[grp_id][0];
    // int aid0_3 = aid0 * 3;
    // cout << "dbg1130 grp " << i_grp << " " << grp_id << " " << n_atoms_in_groups[grp_id] << endl;
    // real crd0[3] = {crd[aid0_3], crd[aid0_3 + 1], crd[aid0_3 + 2]};
    for (int d = 0; d < 3; d++) { crd_centers[i_grp][d] = 0.0; }
    for (int i_atm = 0; i_atm < n_atoms_in_groups[grp_id]; i_atm++) {
      int aid   = atom_groups[grp_id][i_atm];
      int aid_3 = aid * 3;
      for (int d = 0; d < 3; d++) {
	real tmp_crd = crd[aid_3 + d];
	real diff    = tmp_crd - crd_groups[i_grp][i_atm][d];
	while (diff > pbc->L_half[d]) {
	  diff -= pbc->L[d];
	  tmp_crd -= pbc->L[d];
	}
	while (-diff > pbc->L_half[d]) {
	  diff += pbc->L[d];
	  tmp_crd += pbc->L[d];
	}
	crd_groups[i_grp][i_atm][d] = tmp_crd;
	crd_centers[i_grp][d] += tmp_crd * mass[aid];
      }
      // cout << "dbg1130 crd " << crd_groups[i_grp][i_atm][0] << " "
      //<< crd_groups[i_grp][i_atm][1] << " "
      //<< crd_groups[i_grp][i_atm][2] << " "
      //<< endl;
    }
    for (int d = 0; d < 3; d++) { crd_centers[i_grp][d] *= mass_groups_inv[grp_id]; }
  }
  //cout << "dbg1126 center " << crd_centers[0][0] << " "
  //<< crd_centers[0][1] << " "
  //<< crd_centers[0][2] << " "
  //<< crd_centers[1][0] << " "
  //<< crd_centers[1][1] << " "<< crd_centers[1][2] << endl;
  return 0;
}

bool ExtendedVcMD::is_in_range(){
  bool flg=true;
  for(int d=0; d<n_dim; d++){
    if(lambda[d] >= vc_range_max[d][cur_vs[d]])     { return false; }
    else if(lambda[d] < vc_range_min[d][cur_vs[d]]) { return false; }
  }
  return flg;
}

int ExtendedVcMD::apply_bias(unsigned long cur_step,
			     real_fc *work, int n_atoms_box) {
  
  if(is_in_range() || drift > 0){
    q_raw[cur_vs] += 1;
    q_raw_sum += 1;
  }
  if (cur_step > 0 && cur_step % trans_interval == 0) {
    //if ((cur_step+1) % trans_interval == 0) {
    write_vslog(cur_step);

    if (flg_vs_transition) trial_transition();
  }
  scale_force(work, n_atoms_box);
  if (cur_step > 0 && cur_step % write_lambda_interval == 0 &&
      cur_step <= n_steps) { write_lambda(); }
  
  return 0;
}
int ExtendedVcMD::set_vs_next(){
  std::vector< std::vector<int> > vs_next_crd(n_dim);
  // for(int i=0; i<n_dim; i++) vs_next_crd[i]=-1;
  // vs_next_crd[dimension] = ID of the overlapping virtual states
  for(int d=0; d<n_dim; d++){
    vs_next_crd[d].push_back(cur_vs[d]);
    int vc_size = vc_range_min[d].size();
    if ( cur_vs[d] > 0 ){

      if (lambda[d] < vc_range_max[d][cur_vs[d]-1]){
	vs_next_crd[d].push_back(cur_vs[d]-1);
      }else if(cur_vs[d] < vc_size-1 && lambda[d] >= vc_range_min[d][cur_vs[d]+1] ){
	vs_next_crd[d].push_back(cur_vs[d]+1);
      }
    }else if(vc_range_min[d].size() > 1 && lambda[d] >= vc_range_min[d][cur_vs[d]+1]){
      vs_next_crd[d].push_back(cur_vs[d]+1);
    }
    //cout << "dbg 0719 d:" <<d << " " << cur_vs[d] << " - " << *(--vs_next_crd[d].end()) <<endl;
  }
  vs_next.clear();
  std::vector<int> tmp_vs_next;
  set_vs_next_sub(tmp_vs_next, vs_next_crd);
  return 0;
}
int ExtendedVcMD::set_vs_next_sub(const std::vector<int> tmp_vs_next,
				  const std::vector< std::vector<int> > vs_next_crd){  
  int cur_dim = tmp_vs_next.size();
  for(const auto next_crd : vs_next_crd[cur_dim]){
    std::vector<int> tmp_vs_next_sub;
    copy(tmp_vs_next.begin(), tmp_vs_next.end(), back_inserter(tmp_vs_next_sub));
    tmp_vs_next_sub.push_back(next_crd);
    if(cur_dim == n_dim-1){
      vs_next.push_back(tmp_vs_next_sub);
      //cout << "dbg 0719 push-next ";
      //for ( const auto x: tmp_vs_next_sub )
      //std::cout << x << "-" ;
      //std::cout << endl;
    }else{
      set_vs_next_sub(tmp_vs_next_sub, vs_next_crd);
    }
  }

  return 0;
}

int ExtendedVcMD::trial_transition(){  // source ... vs_id of current state
  // rel_dest ... -1 or 1, down or up
  // lambda
  
  // return ...
  //cout << "dbg 0304 trial [1]" << endl;
  //  bool flg = true;
  
  //for(int d = 0; d < n_dim; d++){
    //cout << "dbg 0304 trial d[" << d << "] " << lambda[d]
    //<< " cur_vs: " <<cur_vs[d] 
    //<<" "<<vc_range_min[d][cur_vs[d]]<<" ~ "
    //<<vc_range_max[d][cur_vs[d]]<< endl;
  //if     (lambda[d] >= vc_range_max[d][cur_vs[d]] and drift == 0) { flg = false; break; }
  //else if(lambda[d] <  vc_range_min[d][cur_vs[d]] and drift == 0) { flg = false; break; }
  //}
  if (!is_in_range() and drift==0) { return 0; }
  //cout << "dbg 0304 trial [2]" << endl;
  /*
  std::vector< std::vector<int> > vs_next;
  // vs_next[0,1,2,3] = [x1,y1],[x1,y2],[x2,y1],[x2,y2]
  // the candidate states for the next step
  if (n_dim == 1){
    std::vector<int> tmp_vs(1);
    tmp_vs[0] = cur_vs[0];
    vs_next.push_back(tmp_vs);
    if(vs_next_crd[0] >= 0){
      tmp_vs[0] = vs_next_crd[0];
      vs_next.push_back(tmp_vs);
    }
    //std::cout << "dbg 0304 vs_next : "  << endl;
    //for(const auto x: vs_next)
    //std::cout << x[0] <<  endl;
    
  }else if (n_dim == 2){
    std::vector<int> tmp_vs(2);
    tmp_vs[0] = cur_vs[0];      tmp_vs[1] = cur_vs[1];
    vs_next.push_back(tmp_vs);
    if ( vs_next_crd[0] >= 0 ){
      tmp_vs[0] = vs_next_crd[0]; tmp_vs[1] = cur_vs[1];
      vs_next.push_back(tmp_vs);
      if ( vs_next_crd[1] >= 0 ){
	tmp_vs[0] = vs_next_crd[0]; tmp_vs[1] = vs_next_crd[1];
	vs_next.push_back(tmp_vs);
      }
    }
    if ( vs_next_crd[1] >= 0 ){
      tmp_vs[0] = cur_vs[0];      tmp_vs[1] = vs_next_crd[1];
      vs_next.push_back(tmp_vs);
    }
    //std::cout << "dbg 0304 vs_next : "  << endl;
    //for(const auto x: vs_next)
    //std::cout << x[0] << "-" << x[1] << endl;
    
  }else{
    error_exit(string("In this version, the VcMD allows upto 2 dimension."), "1A00006");
  }
  */

  set_vs_next();
  
  std::vector<real> i_val;
  i_val.resize(vs_next.size());
  std::vector<real> i_val_d;
  i_val_d.resize(vs_next.size());
  int idx_vs1=0;
  real sum_i_val = 0.0;

  int idx_vs2=0;
  real sum_q2 = 0.0;
  real sum_q2d = 0.0;
  //    cout << "dbg 0606a q1 " << q1 << endl ;
  for ( const auto vs : vs_next ) {
    real q2 = q_cano[vs];
    if(q2 == 0) q2 = default_q_cano;
    real q2d = 1.0;
    if(drift >= 1) q2d = ((float)q_raw[vs]+1) / ((float)q_raw_sum+1);
    //i_val[idx_vs1] += q_cano[vs1] / q_cano[vs2];
    //0606
    //i_val[idx_vs1] += q1/q2;
    //i_val_d[idx_vs1] += q1d/q2d;
    sum_q2 += 1/q2;
    sum_q2d += 1/q2d;
    idx_vs2++;
  }
  for ( const auto vs : vs_next ) {
    real q1 = q_cano[vs];
    if(q1==0) q1 = default_q_cano;
    i_val[idx_vs1] = -(log(q1) + log(sum_q2));

    real q1d = 1.0;
    if(drift >= 1){
      q1d = ((float)q_raw[vs]+1) / ((float)q_raw_sum+1);
      i_val[idx_vs1] += -(log(q1d) + log(sum_q2d)); 
    }

    i_val[idx_vs1] = exp(i_val[idx_vs1]);
    sum_i_val += i_val[idx_vs1];

    if(i_val[idx_vs1] < 0){
      std::cout << "!! (Dbg0605) Transition probability is negative." << endl;
      std::cout << "q1 : " << q1 << endl;
      for ( const auto vs2 : vs_next ) {
	std::cout << "q2 : " << q_cano[vs2] << endl;
      }
      exit(1);
    }
    idx_vs1++;
  }
  
  real rnd = (*random_mt)() * sum_i_val;
  real i_acc = 0.0;
  int idx = 0;
  for ( const auto val : i_val ){
    i_acc += val;
    if ( rnd < i_acc ) break;
    idx++;
  }
  cur_vs = vs_next[idx];

  //cout << "dbg 0909 cur_vs: ";
  //for (int i = 0; i < n_dim; i++)
  //cout << cur_vs[i] << " " ;
  //cout <<endl;

  return idx;
}

int ExtendedVcMD::scale_force(real_fc *work, int n_atoms) {
  for ( int d = 0; d < n_dim; d++){
    real param    = lambda[d];
    real recovery = 0.0;
    //cout << "dbg 0304 scale d:"<<d<< " lambda:"<<lambda[d]
    //<< " cur_vs:" << cur_vs[d] 
    //<< " " << vc_range_min[d][cur_vs[d]]<< "~" 
    //<< " " << vc_range_max[d][cur_vs[d]]<< endl;
    
    // case 1 : under the lower limit
    if (lambda[d] < vc_range_min[d][cur_vs[d]] - sigma){
      //cout << "dbg 0304 scale d:"<<d<< " min lambda:"<<lambda[d] << " < " << vc_range_min[d][cur_vs[d]] << " vs:" << cur_vs[d] <<endl;
      recovery = recov_coef * (lambda[d] - (vc_range_min[d][cur_vs[d]] - sigma));
      param = vc_range_min[d][cur_vs[d]];
    }else if (lambda[d] >= vc_range_max[d][cur_vs[d]] + sigma){
      //cout << "dbg 0304 scale d:"<<d<< " max lambda:"<<lambda[d] << " > " << vc_range_max[d][cur_vs[d]] << " vs:" << cur_vs[d]<< endl;
      recovery = recov_coef * (lambda[d] - (vc_range_max[d][cur_vs[d]] + sigma));
      param = vc_range_max[d][cur_vs[d]];
    }
    int i_pair = 0;
    real dew = const_k * recovery;
    //    cout << "dbg 0304 scale d:"<<d<<" recov: " << recovery << " dew:" << dew << endl;
    real direction = 1.0;
    for (int pair_ab = 0; pair_ab < 2; pair_ab++) {
      int i_grp = enhance_group_pairs[d][pair_ab];
      int grp_id = enhance_groups[i_grp];
      real bias[3];
      //cout << "dbg 0304 scale pair " << pair_ab << " i_grp: " << i_grp << " grp_id: " << grp_id << " mass: " << mass_groups[grp_id] << " mass_inv: " << mass_groups_inv[grp_id] << endl;
      for (int xd = 0; xd < 3; xd++)
	bias[xd] = dew * unit_vec[d][xd] * (real)mass_groups_inv[grp_id];
      //cout << " bias " << d << " " << pair_ab << " " << bias[0]*direction << " " << bias[1]*direction << " " << bias[2]*direction << endl; 
      for (int i_at = 0; i_at < n_atoms_in_groups[grp_id]; i_at++) {
	//cout << " " << atom_groups[grp_id][i_at]; 
	int atom_id3 = atom_groups[grp_id][i_at] * 3;
	for (int xd = 0; xd < 3; xd++) {
	  work[atom_id3 + xd] += direction * bias[xd] * (real)mass[atom_groups[grp_id][i_at]];
	}
      }
      //cout << endl;
      direction *= -1.0;
    }
    // cout << "dbg 0304 scale " << endl;
    i_pair++;
  }
  return 0;
}

int ExtendedVcMD::set_files(string fn_vslog, string fn_lambda, int format_lambda,
			    string fn_qraw, string fn_start) {
    writer_vslog.set_fn(fn_vslog);
    writer_vslog.open();
    if (format_lambda == LAMBDAOUT_BIN) {
        writer_lambda = new WriteTableLogBinary();
    } else if (format_lambda == LAMBDAOUT_ASC) {
        writer_lambda = new WriteTableLogAscii();
    } else {
        writer_lambda = new WriteTableLog();
    }
    writer_lambda->set_fn(fn_lambda);
    writer_lambda->open();
    writer_lambda->set_ncolumns(n_dim);
    writer_lambda->write_header();
    // write_vslog(0);
    writer_qraw.set_fn(fn_qraw);
    writer_start.set_fn(fn_start);
    return 0;
}
int ExtendedVcMD::write_vslog(int cur_steps) {
  writer_vslog.write_VcMDLog(cur_steps, cur_vs);
  return 0;
}
int ExtendedVcMD::write_lambda(){
  writer_lambda->write_row(lambda);
    return 0;
}
int ExtendedVcMD::write_q(){
  writer_qraw.open();
  writer_qraw.write(trans_interval,
		    grp_names,
		    vc_range_min, vc_range_max,
		    q_raw);
  writer_qraw.close();

  writer_start.open();
  writer_start.set_ncolumns(1);
  int rnd = (int)((*random_mt)() * 1000000);
  writer_start.write_row(&n_dim);
  for(const auto itr : cur_vs){
    int vs = itr + 1;
    writer_start.write_row(&vs);
  }
  writer_start.write_row(&rnd);
  writer_start.close();
  return 0;
}
int ExtendedVcMD::print_info() {
  std::cout << "V-McMD parameters" << endl;
  for (int d = 0; d < n_dim; d++){
    std::cout << "  Dimension : " << d+1 << endl;
    for (unsigned int i = 0; i < vc_range_min[d].size(); i++) {
      std::cout << "    Virtual state: " << i+1 << " ... ";
      std::cout << vc_range_min[d][i] << " - ";
      std::cout << vc_range_max[d][i] << endl;
    }
  }
  cout << "initial virtual states :" << endl;
  for ( const auto vs : init_vs ) {
    cout << " " << vs+1;
  }
  cout << endl;
  cout << "seed : " << random_seed << endl;
  cout << "param size: " << q_cano.size() << endl;
  for ( const auto pair : q_cano ) {
    for ( const auto vs : pair.first ) { 
      std::cout << " " << vs;
    }
    std::cout << " : " << pair.second << endl;
  }
  cout << "end print" << endl;
  return 0;
}

//real ExtendedVcMD::cal_struct_parameters(real *crd, PBC *pbc) {
//    return 0.0;
//}

int ExtendedVcMD::set_enhance_groups(int *       in_n_atoms_in_groups,
				     int **      in_atom_groups,
				     int         in_n_enhance_groups,
				     std::vector<int> in_enhance_groups) {
  
  n_atoms_in_groups = in_n_atoms_in_groups;
  atom_groups       = in_atom_groups;
  n_enhance_groups  = in_n_enhance_groups;
  enhance_groups    = in_enhance_groups;
  
  //n_enhance_group_pairs = (n_enhance_groups * (n_enhance_groups - 1)) / 2;
  int i_pair            = 0;
  for ( const auto itr_pair : grp_ids ){
    std::vector<int> grp_pair;
    for ( const auto itr_grpid : itr_pair ){    
      int grpid_in_enhance = 0;
      for ( const auto itr_grp_vc  : enhance_groups ){
	if ( itr_grp_vc == itr_grpid ) break;	  
	grpid_in_enhance++;
      }
      grp_pair.push_back(grpid_in_enhance);
    }
    enhance_group_pairs.push_back(grp_pair);
    i_pair++;
  }
  //cout << "dbg 0304[1a] group_pairs size : " << enhance_group_pairs.size() << endl;
  int tmp = 0;
  for( const auto itr : enhance_group_pairs ){
    //cout << "dbg 0304 1b enhance_gropu_pairs: " << tmp << " :";
    for (const auto itr2 : itr){
      cout << " " << itr2;
    }
    cout << endl;
    tmp++;
  }
  alloc_crd_centers();
  return 0;
}

int ExtendedVcMD::alloc_crd_centers() {
  cout << "dbg 0304 1a[2]" << endl;
  crd_groups = new real **[n_enhance_groups];
  for (int i_grp = 0; i_grp < n_enhance_groups; i_grp++) {
    int grp_id        = enhance_groups[i_grp];
    crd_groups[i_grp] = new real *[n_atoms_in_groups[grp_id]];
    for (int j = 0; j < n_atoms_in_groups[grp_id]; j++) { crd_groups[i_grp][j] = new real[3]; }
  }
  cout << "dbg 0304 1a[3] " << n_enhance_groups << endl;  
  crd_centers = new real *[n_enhance_groups];
  for (int i = 0; i < n_enhance_groups; i++) { crd_centers[i] = new real[3]; }
  cout << "dbg 0304 1a[4] " << n_dim << endl;  
  unit_vec = new real *[n_dim];
  for (int i = 0; i < n_dim; i++) { unit_vec[i] = new real[3]; }
  cout << "dbg 0304 1a[5]" << endl;
  return 0;
}
int ExtendedVcMD::free_crd_centers() {
  for (int i = 0; i < n_enhance_groups; i++) {
    int grp_id = enhance_groups[i];
    for (int j = 0; j < n_atoms_in_groups[grp_id]; j++) { delete[] crd_groups[i][j]; }
    delete[] crd_groups[i];
  }
  delete[] crd_groups;
  
  for (int i = 0; i < n_enhance_groups; i++) { delete[] crd_centers[i]; }
  delete[] crd_centers;
  for (int i = 0; i < n_dim; i++) { delete[] unit_vec[i]; }
  delete[] unit_vec;
  return 0;
}

int ExtendedVcMD::set_mass(real_pw *in_mass, real_pw *in_mass_groups, real_pw *in_mass_groups_inv) {
  mass = in_mass;
  mass_groups     = in_mass_groups;
  mass_groups_inv = in_mass_groups_inv;
  mass_sum        = 0.0;
  //for (int i_grp = 0; i_grp < n_enhance_groups; i_grp++) {
  //int grp_id = enhance_groups[i_grp];
  //for (int i_atm = 0; i_atm < n_atoms_in_groups[grp_id]; i_atm++) {
  //int aid = atom_groups[grp_id][i_atm];
  //mass_sum += mass[aid];
  //}
  //}
  
  return 0;
}
int ExtendedVcMD::write_aus_restart(std::string fn_out) {
    WriteGroupCoord writer;
    writer.set_fn(fn_out);
    writer.open();
    writer.write_aus_restart(reactcrd_type, n_enhance_groups,
			     enhance_groups,
			     n_atoms_in_groups, crd_groups);
    writer.close();
    return 0;
}
