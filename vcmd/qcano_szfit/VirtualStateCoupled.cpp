#include "VirtualStateCoupled.h"
  
VirtualStateCoupling::VirtualStateCoupling()
{
  default_weight = 1.0;
}  
  
VirtualStateCoupling::~VirtualStateCoupling(){

}

std::string VirtualStateCoupling::get_str_state_definition(int param_mode){
  std::stringstream ss;
  ss << exchange_interval;
  if ( param_mode == 1 )
    ss << " LOG";
  ss << std::endl;
  ss << n_dim << std::endl;

  for ( int d = 0; d < n_dim; d ++){
    ss << lowers_vaxis[d].size();
    for ( const auto itr : aus_groups[d] ){
      ss << " " << itr;
    }
    ss  << std::endl;
    for ( int l = 0; l < lowers_vaxis[d].size(); l++){
      ss << lowers_vaxis[d][l] << " " << uppers_vaxis[d][l] << std::endl;
    }
  }
  return ss.str();
}

void VirtualStateCoupling::set_cur_vsis(std::vector<int> v_crd, std::vector<double> arg){
  //cout << "dbg state : " << current_state << endl;
  // Updating current virtual state with intersection area
  //cur_vsis.clear();    
  int i=0;
  for( const auto crd : v_crd ) { cur_vsis[i++] = crd; }
  for( int d = 0; d < n_dim; d++){
    if( arg[d] < lowers_vaxis[d][v_crd[d]] ){
      cur_vsis[i++] = -1;
    }else if( arg[d] >= uppers_vaxis[d][v_crd[d]] ){
      cur_vsis[i++] = 2;
    }else{
      if( arg[d] < uppers_vaxis[d][v_crd[d]-1] ){
	cur_vsis[i++] = 0;
      }else{
	cur_vsis[i++] = 1;
      }
    }
  }
}  
void VirtualStateCoupling::parse_vstate(const string &fname){
  //(void) fname;
  ifstream ifs(fname.c_str());
  //ifstream ifs("state.dat");
  //ifs.open(fname.c_str());
  if(!ifs){
    printf("Failed to open the file %s\n",fname.c_str());
    exit(1);
  }

  std::vector<string> args;
  string         buf;
  string         cur;

  //printf(" dbg 1\n");

  int i_line=0;
  while (ifs && getline(ifs, buf)) {
    //printf("%s\n", buf.c_str());
    unsigned int pos1 = buf.find_first_of("#;");
    if (pos1 != string::npos) { buf = buf.substr(0, pos1); }
    if (buf.size() == 0) continue;
    stringstream ss(buf);
    while (ss >> buf) {
      args.push_back(buf);
      //printf("%s\n", buf.c_str());
    }
  }
  ifs.close();

  std::vector<int> vs_crd;
  for(size_t d = 0; d < n_dim; ++d) {
    // (-1) convert to 0-origin
    int cur_state=atoi(args[d+1].c_str()) - 1;
    vs_crd.push_back(cur_state);
    if(cur_state >= nstates_vaxis[d]){
      std::cerr << "Error: Invalid virtual system coordinate is specified in " << fname << std::endl;
      exit(1);
    }
  }
  int vs_id = conv_vstate_crd2id(vs_crd);
  //printf( " dbg 2\n" );
  current_state = vs_id;
  //printf( " dbg 3\n" );
  //current_state = 0;
}

void VirtualStateCoupling::parse_params(const string &fname)
{
  ifstream ifs;
  ifs.open(fname.c_str());
  //dimension

  nstates = 1;
  int param_mode = parse_params_state_definition(ifs);
  parse_params_qweight(ifs, param_mode);
  ifs.close();
}

void VirtualStateCoupling::parse_qraw_is(const string &fname)
{
  ifstream ifs;
  ifs.open(fname.c_str());
  //dimension
  nstates = 1;

  parse_params_state_definition(ifs);
  parse_params_qraw_is(ifs);
  ifs.close();
}


int VirtualStateCoupling::parse_params_state_definition(ifstream &ifs){
  vector<string> args;
  string         buf, buf2;
  string         cur, cur1, cur2;
  getline(ifs, buf);
  stringstream buf_ss(buf);
  int param_mode = 0;
  buf_ss >> buf2;
  //cout << "dbg : " <<  buf << endl;;
  exchange_interval = atoi(buf2.c_str());
  if( buf_ss >> buf2 ){
    if (buf2 == "LOG"){
      param_mode = 1;
    }
  }
  printf("exchange_interval: %d \n", exchange_interval);
  getline(ifs, buf);
  n_dim = atoi(buf.c_str());
  printf("n_dim: %d \n", n_dim);
  // read lower and upper bounds for each state in each axis
  for (size_t c_dim=0; c_dim < n_dim; c_dim++){
    int cur_nstates; 
    getline(ifs, buf);
    stringstream ss(buf);
    ss >> cur;  cur_nstates = atoi(cur.c_str());
    vector<string> groups;
    while(ss >> cur) groups.push_back(cur);
    aus_groups.push_back(groups);
    //printf("%d-th dim : %d states\n", c_dim, cur_nstates);
    vector<double> buf_lowers;
    vector<double> buf_uppers;
    for(int c_state = 0; c_state < cur_nstates; c_state++){
      getline(ifs, buf);
      stringstream ss2(buf);
      //ss2 << buf;
      ss2 >> cur1 >> cur2;
      //printf("dbg param read: %s %s\n", cur1.c_str(), cur2.c_str());
      buf_lowers.push_back(atof(cur1.c_str()));
      buf_uppers.push_back(atof(cur2.c_str()));
    }
    lowers_vaxis.push_back(buf_lowers);
    uppers_vaxis.push_back(buf_uppers);
    nstates_vaxis.push_back(cur_nstates);
    nstates *= cur_nstates;
  }

  // the information of lower/upper bounds for each axis 
  //  are converted into vector for each state-id
  lowers = vector< vector<double> >(nstates);
  uppers = vector< vector<double> >(nstates);
  //kappa = vector<double>(nstates);
  state_weights = vector<double>(nstates);
  state_adj_qw = vector<double>(nstates);
  state_adj_qw_opt = vector<double>(nstates);
  state_flg = vector<int>(nstates);

  for(size_t v_id=0; v_id < nstates; v_id++){
    state_adj_qw[v_id] = log(1.0/nstates);
    state_adj_qw_opt[v_id] = log(1.0/nstates);
  }


  state_qraw = vector<double>(nstates);
  transition_candidates = vector< vector<size_t> >(nstates);
  
  for(size_t v_id=0; v_id < nstates; v_id++){
    vector<int> v_crd = conv_vstate_id2crd(v_id);
    vector<double> tmp_upp(n_dim);
    vector<double> tmp_low(n_dim);
    for(size_t d=0; d<n_dim; d++){
      tmp_upp[d] = uppers_vaxis[d][v_crd[d]];
      tmp_low[d] = lowers_vaxis[d][v_crd[d]];
    }
    uppers[v_id] = tmp_upp;
    lowers[v_id] = tmp_low;
    //printf("dbg read param: v_id %ld %d ... [%lf ... %lf]\n",v_id, v_crd[0], lowers[v_id][0], uppers[v_id][0]);
    state_weights[v_id] = 1.0;
    state_qraw[v_id] = 0.0;
    //printf("v_id %d (", v_id);
    //for(int d=0; d<n_dim; d++){
    //printf("%d ", v_crd[d]);
    //}
    //printf(")\n");
  }

  // read information about parameters
  return param_mode;

}
void VirtualStateCoupling::parse_params_qweight(ifstream &ifs, int param_mode){
  vector<string> args;
  string         buf;
  string         cur, cur1, cur2;
  
  while(ifs && getline(ifs, buf)){
    size_t pos1 = buf.find_first_of("#;");
    if (pos1 != string::npos) { buf = buf.substr(0, pos1); }
    if (buf.size() == 0) continue;

    if(buf=="END") break;
    vector<int> v_crd(n_dim);

    stringstream ss(buf);

    bool flg_default = false;
    for (size_t c_dim=0; c_dim < n_dim; c_dim++){
      ss >> cur;
      v_crd[c_dim] = atoi(cur.c_str())-1;
      if(atoi(cur.c_str())-1 < 0) flg_default=true;
    }
    ss >> cur;
    double cur_param = atof(cur.c_str()) ;
    if(param_mode == 0) cur_param = log(cur_param);
    if(flg_default){
      default_weight = cur_param;
      continue;
    }
    
    size_t v_id = conv_vstate_crd2id(v_crd);
    state_weights[v_id] = cur_param;
  }
  //

  for(size_t v_id=0; v_id < nstates; v_id++){
    if(state_weights[v_id] >= 1.0){
      state_weights[v_id] = default_weight;      
    }
  }  
  
  //for(size_t v_id=0; v_id < nstates; v_id++){
  //vector<int> v_crd = conv_vstate_id2crd(v_id);
  //for(int d=0; d<n_dim; d++){
  //    }
  //}
  
  cur_vsis.reserve(n_dim*2);
  cur_vsis.assign(n_dim*2, 0);
}

void VirtualStateCoupling::parse_params_qraw_is(ifstream &ifs){
  vector<string> args;
  string         buf;
  string         cur, cur1, cur2;
  state_qraw_sum = 0;


  while(ifs && getline(ifs, buf)){
    size_t pos1 = buf.find_first_of("#;");
    if (pos1 != string::npos) { buf = buf.substr(0, pos1); }
    if (buf.size() == 0) continue;
    if(buf=="END") break;

    vector<int> v_crd(n_dim);

    stringstream ss(buf);

    for (size_t c_dim=0; c_dim < n_dim; c_dim++){
      ss >> cur;
      v_crd[c_dim] = atoi(cur.c_str())-1;
    }
    size_t v_id = conv_vstate_crd2id(v_crd);
    vector<size_t> is_crd(n_dim+1);
    is_crd[0] = v_id;
    bool flg_out = false;    
    for (size_t c_dim=0; c_dim < n_dim; c_dim++){
      ss >> cur;
      int is_crd_in_state = atoi(cur.c_str())-1;
      if ( is_crd_in_state != 0 && is_crd_in_state != 1 )
	flg_out = true;
      is_crd[c_dim+1] = is_crd_in_state;
    }
    
    ss >> cur;
    double cur_param = atof(cur.c_str()) ;
    if (!flg_out) {
      state_qraw_is[is_crd] = cur_param;
      state_qraw[v_id] += cur_param;
      state_qraw_sum += cur_param;
    }
  }


  // set default value
  if(process_unsampled_zone == UNSAMPLE_MIN){ 
    default_weight = 1e10;
    for(size_t v_id=0; v_id < nstates; v_id++){
      if(state_qraw[v_id] < default_weight && state_qraw[v_id] > 0)
	default_weight = state_qraw[v_id];
    }  

    for(size_t v_id=0; v_id < nstates; v_id++){
      if(state_qraw[v_id] == 0.0){
	cout << "dbg set default state_qraw " << v_id << " "<< default_weight << endl;
	state_qraw[v_id] = default_weight;      
      }
    }  
  }else if(process_unsampled_zone == UNSAMPLE_OMIT){ 
    for(size_t v_id=0; v_id < nstates; v_id++){
      if(state_qraw[v_id] == 0.0){
	state_adj_qw[v_id] = 1.0;
	state_adj_qw_opt[v_id] = 1.0;
      }
    }
  }

  for (const auto& [key, value] : state_qraw_is){
    if (state_qraw[key[0]] == 0){
      if(process_unsampled_zone == UNSAMPLE_MIN){       
	state_qraw_is[key] = log(1.0/pow(2, n_dim));
      }
    } else if(value == 0){
      state_qraw_is[key] = 1.0;
    }else{
      state_qraw_is[key] = log(value/state_qraw[key[0]]);
    }
    //cout << "dbg parse : ";
    //for ( const auto& a : key ){
    //cout << a << " ";
    //}
    //cout <<  state_qraw_is[key] << endl;
  }
  for (int st=0; st < nstates; st++){
    if (state_qraw[st] > 0)
      state_qraw[st] = log(state_qraw[st]/state_qraw_sum);
    else
      state_qraw[st] = 1.0;
    //    cout << "dbg " << st << " "  << state_qraw[st] << endl;
  }
  

}


void VirtualStateCoupling::init_transition_table()
{
  // FIXME TODO: this loop requires O(nstate^2), may be a problem for nstate >= 10^6
  // One workaround is to input transition information as well from the external script.
  transition_candidates.clear();
  for(size_t i = 0; i < nstates; ++i) {
    //printf("dbg: transition from state %ld : ", i);
    transition_candidates.push_back(vector<size_t>());
    for(size_t j = 0; j < nstates; ++j) {
      bool transitionable = true;
      for(size_t d = 0; d < n_dim; ++d) {
        // FIXME TODO: Do this condition enough?
        if((uppers[i][d] <= lowers[j][d]) ||
           (uppers[j][d] <= lowers[i][d])){ 
          transitionable = false;
          break;
        }
      }
      if(transitionable) {
        transition_candidates[i].push_back(j);
	//printf("%ld ", j);
      }

    }
    //printf("\n");
    //for(size_t d = 0; d < n_dim; ++d) {
      //printf("[%lf ... %lf] \n", lowers[i][d], uppers[i][d]);
    //}
  }

  size_t tottransition = 0;
  for(size_t i = 0; i < nstates; ++i) {
    tottransition += transition_candidates[i].size();
  }
  printf("  Total No. of possible transitions: %ld\n", tottransition);

}

void VirtualStateCoupling::init_data(){
  opt_error = 1e10;
  state_qraw_sum = 0;
  
}

size_t VirtualStateCoupling::conv_vstate_crd2id(vector<int> vcrd){
  size_t state_id = 0;
  size_t offset = 1;
  for(size_t i=0; i < n_dim; i++){
    state_id += offset * vcrd[i];
    offset *= nstates_vaxis[i];
  }
  return state_id;
}
  
vector<int> VirtualStateCoupling::conv_vstate_id2crd(size_t v_id){
  vector<int> v_crd(n_dim);
  
  size_t offset = 1;
  //std::cout  << "conv vstate id2crd : " << v_id << " " << nstates << endl;
  for(int i=0; i < n_dim; i++){
    //std::cout  << " " << nstates_vaxis[i] << " " << offset  << endl;
    v_crd[i] = (v_id / offset) % nstates_vaxis[i];
    offset *= nstates_vaxis[i];
  }

  /*size_t offset = nstates;
  std::cout  << "conv vstate id2crd : " << v_id << " " << nstates << endl;
  for(int i=n_dim-1; i >= 0; i--){
    size_t buf_v_id = v_id % offset;
    offset /= nstates_vaxis[i];
    std::cout  << buf_v_id << " " << nstates_vaxis[i] << " " << offset  << endl;
    v_crd[i] = buf_v_id / offset;
  }
  */
  return v_crd;
}


bool VirtualStateCoupling::is_in_range(size_t state, std::vector<double> arg){
  //std::cout << "dbg iir : " << state << " " << std::endl;
  for( int d = 0; d < n_dim; d++){
    //std::cout << d << " " << lowers[state][d] << " " << uppers[state][d] << " " << arg[d] << std::endl;
    if( lowers[state][d] > arg[d] ||
	uppers[state][d] <=  arg[d] ){
      //std::cout << "false" << std::endl;
      return false;
    }
  }
  return true;
}

void VirtualStateCoupling::write_qweight(std::string fname, std::vector<double> in_qw, bool def_val, int param_mode){
  //cout << "write_qweight " << param_mode << endl;

  std::string state_defs = get_str_state_definition(qweight_write_mode);
  ofstream ofs(fname);
  ofs << state_defs;
  
  if (def_val){
    double min_qw = 1e10;
    for ( int l = 0; l < nstates; l++){
      //cout << "min : " << l << " "  <<  in_qw[l] << " " << min_qw << endl;
      if ( in_qw[l] < min_qw and in_qw[l] < 0){
	min_qw = in_qw[l];
	//cout << "min update: " << l << " "  <<  in_qw[l] << " " << min_qw <<  endl;
      }
    }
    for ( int d = 0; d < n_dim; d++){
      ofs << "0 ";
    }
    double put_qw = min_qw;
    if(param_mode == 0) put_qw = exp(min_qw);
    ofs << scientific << put_qw << std::endl;
  }
  for ( int l = 0; l < nstates; l++){
    if( in_qw[l] > 0.0 ){ continue;  }
    std::vector<int> vs_crd = conv_vstate_id2crd(l);
    for ( const auto vsc : vs_crd )
      // (+1) convert to 1-origin integer
      ofs << vsc+1 << " ";
    double put_qw = in_qw[l];
    if(param_mode == 0) put_qw = exp(in_qw[l]);
    ofs << scientific << put_qw << std::endl;
  }
  ofs << "END" << std::endl;
  ofs.close();
  return;
}
/*
void VirtualStateCoupling::write_qrawstat_is(){
  std::string state_defs = get_str_state_definition();
  ofstream ofs_qrawstat_is(fname_o_qrawstat_is);
  ofs_qrawstat_is << state_defs;
  for ( const auto q : state_qraw_is ){
    stringstream ss;
    bool flg = true;
    int cur_dim = 0;
    for ( const auto vsc : q.first ){
      // (+1) convert to 1-origin integer
      ss << vsc + 1 << " ";
      //ofs_qrawstat_is << vsc+1 << " ";
      if (cur_dim >= n_dim && (vsc == -1 || vsc == 2)){
	flg = false;
      }
      cur_dim++;
    }
    if(verbose > 0 || flg)
      ofs_qrawstat_is << ss.str() << q.second << std::endl;
  }
  ofs_qrawstat_is << "END" << std::endl;
  ofs_qrawstat_is.close();
  return;
}
void VirtualStateCoupling::write_transit_count(){
  //std::string state_defs = get_str_state_definition();
  ofstream ofs(fname_o_transit_count);
  for ( const auto q : state_transit_count ){
    for ( const auto vsc : q.first )
      // (+1) convert to 1-origin integer
      ofs << vsc+1 << " ";
    ofs << q.second << std::endl;
  }
  ofs.close();
  return;
}
*/
int VirtualStateCoupling::setup(Config cfg){
  cout << "dbg setup()" << endl;
  //cout << " dbg setup : mc_steps: " << mc_steps <<  " " << cfg.mc_steps << endl;
  fname_i_params = cfg.fname_i_params;
  fname_o_qcano = cfg.fname_o_qcano;
  fname_o_qweight_opt = cfg.fname_o_qweight_opt;
  fname_i_qraw_is = cfg.fname_i_qraw_is;
  qweight_write_mode = cfg.qweight_write_mode;
  process_unsampled_zone = cfg.process_unsampled_zone;

  target_error = cfg.target_error;
  // for mc
  mc_temp = cfg.mc_temp;
  mc_delta_x = cfg.mc_delta_x;
  mc_delta_x_max = cfg.mc_delta_x_max;
  mc_steps = cfg.mc_steps;
  mc_log_interval = cfg.mc_log_interval;
  mc_target_acc_ratio = cfg.mc_target_acc_ratio;
  mc_acc_duration = cfg.mc_acc_duration;

  // for greedy search
  greedy_max_steps = cfg.greedy_max_steps;
  greedy_pivot = cfg.greedy_pivot;

  mc_n_window_trend = cfg.mc_n_window_trend;
  mc_error_ave_window_size = cfg.mc_error_ave_window_size;
  mc_max_temp = cfg.mc_max_temp;
  mc_min_temp = cfg.mc_min_temp;
  mc_delta_temp = cfg.mc_delta_temp;
  return 0;
}

int VirtualStateCoupling::enum_overlapping_subzones(){
  // set
  ////  overlapping subzones

  cout << "dbg: enum_overlapping_subzones() " << n_dim << endl;
  overlapping_subzones = vector< std::map< size_t, std::vector< std::vector<int> > > >(nstates);


  for ( int st_i = 0; st_i < nstates; st_i++ ){
    for (const auto& st_j : transition_candidates[st_i] ) {
      if (st_j == st_i) continue;
      std::vector<int> overlap_positions(n_dim);
      for ( int c_dim = 0 ; c_dim < n_dim ; c_dim++ ){
	if(uppers[st_j][c_dim] < uppers[st_i][c_dim])
	  overlap_positions[c_dim] = -1;
	else if(uppers[st_i][c_dim] < uppers[st_j][c_dim])
	  overlap_positions[c_dim] = 1;
	else
	  overlap_positions[c_dim] = 0;
      }

      std::vector<int> tmp_vec;
      std::vector< std::vector<int> > subzones;
      enum_overlapping_subzones_sub(overlap_positions,
				tmp_vec, subzones);
      overlapping_subzones[st_i][st_j] = subzones;
    }
  }
  
  return 0;
}
int VirtualStateCoupling::enum_overlapping_subzones_sub
(std::vector<int>& overlap_positions, std::vector<int>& cur_vec,
 std::vector< std::vector<int> >& subzones){
  int cur_depth = n_dim - overlap_positions.size();
  if(overlap_positions.empty()){
    subzones.push_back(cur_vec);
    return 1;
  }
  int pos = overlap_positions[0];
  overlap_positions.erase(overlap_positions.begin());
  std::vector<int> olp1(overlap_positions.size());
  std::vector<int> olp2(overlap_positions.size());
  for ( int i=0; i < overlap_positions.size(); i++) {
    olp1[i] = overlap_positions[i];
    olp2[i] = overlap_positions[i];
  }
  
  std::vector<int> v1(n_dim*2);
  std::vector<int> v2(n_dim*2);
  for ( int i=0; i < cur_vec.size(); i++){
    v1[i] = cur_vec[i];
    v2[i] = cur_vec[i];
  }

  if ( pos == -1 ){
    v1[cur_depth]   = 0;
    v1[n_dim+cur_depth] = 1;
    enum_overlapping_subzones_sub(olp1, v1, subzones);
  }else if( pos == 1 ){
    v1[cur_depth]   = 1;
    v1[n_dim+cur_depth] = 0;
    enum_overlapping_subzones_sub(olp1, v1, subzones);
  }else{
    v1[cur_depth] = 0;
    v1[n_dim+cur_depth] = 0;    
    enum_overlapping_subzones_sub(olp1, v1, subzones);
    v2[cur_depth] = 1;
    v2[n_dim+cur_depth] = 1;    
    enum_overlapping_subzones_sub(olp2, v2, subzones);
  }
  return 0;
}

int VirtualStateCoupling::print_overlapping_subzones(){
  for ( int st_i = 0 ; st_i < nstates; st_i++) {
    for (const auto& [st_j, subzones] : overlapping_subzones[st_i]){
      for (const auto& sz : subzones ) {
	cout << "ovlp " << st_i << "-" << st_j << " [" ;
	for (const auto& crd : sz ) {
	  cout << crd << " ";
	}
	cout << "]" <<endl;
      }
    }
  }
  return 0;
}
double VirtualStateCoupling::calc_qraw_error_all(){
  // set
  //  state_adj_qw_error
  //  total_err

  state_adj_qw_error = vector<double>(nstates);
  total_error = 0.0;

  for (int i = 0 ; i < nstates; i++){
    if(state_qraw[i] == 0) continue;
    state_adj_qw_error[i] = calc_qraw_error(i, state_adj_qw[i], false);
    total_error += state_adj_qw_error[i];
  }
  total_error *= 0.5;
  //cout << "dbg0626 sqer_sum: " << total_error << endl;
  return total_error;
}

double VirtualStateCoupling::calc_qraw_error(size_t st_i, double qw_i, bool skip_lower_id=false){
  //cout << "abs"<< endl;
  double sqer_sum=0.0;
  for (const auto& [st_j, subzones] : overlapping_subzones[st_i]){
    if(skip_lower_id && st_j < st_i) continue;
    for (const auto& sz : subzones ) {
      std::vector<size_t> sz_crd_i(n_dim+1), sz_crd_j(n_dim+1);
      sz_crd_i[0]=st_i;
      sz_crd_j[0]=st_j;
      for ( int i = 0; i < n_dim; i++){
	sz_crd_i[i+1] = sz[i];
	sz_crd_j[i+1] = sz[n_dim+i];
      }
      //cout << "dbg err sz_crd_i : ";
      //for ( int i = 0; i < n_dim+1; i++){
      //cout << sz_crd_i[i] << " ";
      //}
      //cout << endl;
      //cout << "dbg err x1 "<< qw_i ;
      //cout << " " << state_qraw[st_i] << " " << state_qraw_is[sz_crd_i];
      //cout << " " << state_adj_qw[st_j] << " " << state_qraw[st_j] << " " << state_qraw_is[sz_crd_j] << endl;;
      
      double qcano_i = qw_i + state_qraw[st_i] + state_qraw_is[sz_crd_i];
      double qcano_j = state_adj_qw[st_j] + state_qraw[st_j] + state_qraw_is[sz_crd_j];
      //cout << "dbg err qcano " << qcano_i << " " << qcano_j << " " << endl;
      //if(qcano_i > 0 || qcano_j > 0) continue;
      if(qw_i > 0.0 ||  state_qraw[st_i] > 0.0 || state_qraw_is[sz_crd_i] > 0.0 ||
	 state_adj_qw[st_j] > 0.0 ||  state_qraw[st_j] > 0.0 ||  state_qraw_is[sz_crd_j] > 0.0) continue;
      double sqer = pow(qcano_i-qcano_j,2);
      //double sqer = pow(qcano_i-qcano_j,2);
      if (sqer < 0) cout << "dbg nega error " << st_i << "  " << st_j << " " << sqer << endl;
      sqer_sum += sqer;
      
    }
  }
  return sqer_sum;
}

int VirtualStateCoupling::mode_test(){
  parse_params(fname_i_params);
  init_transition_table();
  return 0;
}

int VirtualStateCoupling::mode_subzonebased_mc(){
  cout << "mode_subzonebased_mc(): " << fname_i_qraw_is << endl;
  parse_qraw_is(fname_i_qraw_is);
  init_transition_table();
  init_data();
  enum_overlapping_subzones();
  //print_overlapping_subzones();
  //cout << "Error : " << err << endl;
  opt_error = calc_qraw_error_all();
  for (int i=0; i < nstates; i++){
    state_adj_qw_opt[i] = state_adj_qw[i];
  }

  cout << "Error init : " << calc_qraw_error_all() << endl;

  if(greedy_pivot>0){
    cout << "greedy" << endl;
    greedy_search(greedy_pivot-1);
  }

  double gr_error = calc_qraw_error_all();
  cout << "Error gr   : " << gr_error << endl;

  if (gr_error < opt_error){
    opt_error = gr_error;
    for (int i=0; i < nstates; i++){
      state_adj_qw_opt[i] = state_adj_qw[i];
    }
  }else{
    for (int i=0; i < nstates; i++){
      state_adj_qw[i] = state_adj_qw_opt[i];
    }
  }

  mc_acc = 0;
  mc_loop();
  cout << "Error mc1  : " << calc_qraw_error_all() << endl;
  cout << "opt err : " << opt_error << endl;
  cout << "mc acc  : " << mc_acc << endl;

  for (int i=0; i < nstates; i++){
    state_adj_qw[i] = state_adj_qw_opt[i];
  }
  
  mc_acc = 0;
  //mc_target_acc_ratio = 0.0;
  mc_steps *= 0.1;
  mc_delta_x *= 0.1;
  mc_temp *= 0.1;
  mc_loop();
  

  cout << "Error mc2  : " << calc_qraw_error_all() << endl;
  cout << "opt err : " << opt_error << endl;
  cout << "mc acc  : " << mc_acc << endl;

  write_qweight(fname_o_qweight_opt, state_adj_qw_opt, true, qweight_write_mode);
  return 0;
}




int VirtualStateCoupling::mc_loop(){
  double cur_temp = mc_temp;

  double delta_x = mc_delta_x;
  double err = calc_qraw_error_all();
  
  size_t cur_step = 0;

  random_device rnd;
  mt19937 mt(rnd());
  uniform_int_distribution<> ri(0,nstates-1);
  uniform_real_distribution<> rd(0,1);

  //annealing
  vector<double> mc_window_error_sum = vector<double>(mc_n_window_trend);
  int current_anneal_trend = -1;
  int prev_anneal_trend = -1;
  size_t step_anneal_on = mc_error_ave_window_size * mc_n_window_trend * 10;
  
  for (int i=0; i<mc_n_window_trend; i++) mc_window_error_sum[i] = 0.0;
  vector<int> acc_rec;

  //for (int i = 0 ; i < nstates; i++){
  //if(state_qraw[i] == 0) cout << "state_qraw[" << i <<  "] == 0"<<endl;
  //}
  

  for (cur_step=0; cur_step < mc_steps; cur_step++){
    size_t v_id_mig = ri(mt);
    bool flg = true;
    int n_loops = 0;
    while( state_qraw[v_id_mig] > 0 ){
      v_id_mig = ri(mt);
      n_loops++;
      if(n_loops > 1e5){
	cerr << "ERROR. state_qraw[v_id_mig] > 0 did not found."  << endl;
	exit(1);
      }
    }
    double delta_qw = (rd(mt)*2-1) * delta_x;
    double new_qw = state_adj_qw[v_id_mig]+delta_qw;
    if ( new_qw > 0.0 ) new_qw = 0.0;
    //double err_cur = total_error;
    //double err_att = calc_qraw_error_all();
    double err_cur = calc_qraw_error(v_id_mig, state_adj_qw[v_id_mig]);
    double err_att = calc_qraw_error(v_id_mig, new_qw);
    //cout << "dbg0711 " <<cur_step << " "<< err_cur << " " << err_att << endl;
    double delta_error = err_att - err_cur;


    bool acc=false;
    if( delta_error <= 0 ){
      acc=true;
    }else{
      acc = ( rd(mt) < exp(-delta_error/cur_temp) );
    }
    if(acc){
      mc_acc += 1;

      state_adj_qw[v_id_mig] = new_qw;
      total_error += delta_error;

      //double delta_qw_acc =  exp(new_qw) - exp(state_adj_qw[v_id_mig]);
      //double dq = log(1.0+delta_qw_acc);
      //double nrm_fact = (1.0+delta_qw_acc);
      double nrm_fact_test = 0.0;
      double sum_qw_test = 0.0;
      for(int st_i=0; st_i < nstates; st_i++){
	if(state_adj_qw[st_i] <= 0)
	  nrm_fact_test += exp(state_adj_qw[st_i]);
      }
      for(int st_i=0; st_i < nstates; st_i++){
	if(state_adj_qw[st_i] <= 0){
	  //state_adj_qw[st_i] -= dq;
	  state_adj_qw[st_i] = state_adj_qw[st_i] - log(nrm_fact_test);
	  sum_qw_test += exp(state_adj_qw[st_i]);
	}
      }
      //cout << "dbg sum_qw_test "<< cur_step << " : " <<  sum_qw_test << " / " << log(nrm_fact_test) << "   " << delta_x << " " << delta_qw << "  " << new_qw << endl;
      if (total_error < opt_error){
	opt_error = total_error;
	for (int i=0; i < nstates; i++){
	  state_adj_qw_opt[i] = state_adj_qw[i];
	}
      }
      //}else{
      //total_err = err_cur;
    }

    // Acceptance ratio
    if(acc) acc_rec.push_back(1);
    else    acc_rec.push_back(0);
    if(acc_rec.size() > mc_acc_duration){
      acc_rec.erase(acc_rec.begin());
    }
    int acc_count = 0;
    for(int i = 0; i < acc_rec.size(); i++){
      if (acc_rec[i] == 1) acc_count+=1;
      
    }
    double acc_ratio = (double)acc_count/(double)acc_rec.size();

    // Updating dX for each attempt
    double factor = acc_ratio/mc_target_acc_ratio;
    //cout << " dbg " << acc_ratio << " " << mc_target_acc_ratio << " " << factor << endl;
    if ( acc_ratio == 0 ) factor = 0.9;
    if ( factor > 1.01 ) factor = 1.01;
    if ( factor < 0.99 ) factor = 0.99;
    //cout << acc_rec.size() << "-" << mc_acc_duration << "-" << mc_target_acc_ratio << endl;
    if(acc_rec.size() == mc_acc_duration ){ 
      if (mc_target_acc_ratio > 0){
	delta_x *= factor;
	if(delta_x > mc_delta_x_max) delta_x = mc_delta_x_max;
      }
    }
    //    cout << cur_step << " " << cur_step % mc_log_interval << endl;

    
    //annealing

    size_t window_id = (cur_step/mc_error_ave_window_size)%mc_n_window_trend;
    mc_window_error_sum[window_id] += total_error;
    if (mc_delta_temp != 0){
      //cout << "dbg anneal chg "  << cur_temp << " " << current_anneal_trend << " "  << mc_delta_temp << " " << 1.0 + current_anneal_trend * mc_delta_temp << endl;;
      cur_temp *= 1.0 + current_anneal_trend * mc_delta_temp;

      if(cur_temp < mc_min_temp) cur_temp = mc_min_temp;
      if(cur_temp > mc_max_temp) cur_temp = mc_max_temp;
    }
	
    if (mc_delta_temp != 0 &&
	cur_step%mc_error_ave_window_size == 0){
      int trend_count = 0;

      for(int i=1; i < mc_n_window_trend; i++){
	int cur_win = window_id + i;
	if ( cur_win >= mc_n_window_trend ) cur_win -= mc_n_window_trend;
	double diff = mc_window_error_sum[cur_win] -  mc_window_error_sum[window_id];
	if (diff > 0) trend_count += 1;
	else          trend_count -= 1;
      }

      if(trend_count != mc_n_window_trend &&
	 trend_count != -mc_n_window_trend &&
	 current_anneal_trend != 0 &&
	 step_anneal_on == 0){
	prev_anneal_trend = current_anneal_trend;
	current_anneal_trend = 0;
	step_anneal_on = cur_step + mc_error_ave_window_size * mc_n_window_trend;
      }
      if (cur_step == step_anneal_on){
	if(current_anneal_trend == 0){
	  current_anneal_trend = -prev_anneal_trend;
	  step_anneal_on = cur_step + mc_error_ave_window_size * mc_n_window_trend * 10;
	}else{
	  step_anneal_on = 0;
	}
      }
      mc_window_error_sum[window_id] = 0.0;
    }
    
    if (cur_step % mc_log_interval ==0){
      cout << "  " << cur_step << " " << total_error << " " << acc_ratio << " " << delta_x << " " << factor << " " << cur_temp << " " << current_anneal_trend << endl;
    }

    if (total_error < target_error ){
      return 1;
    }
  }
  
  return 0;
}
int VirtualStateCoupling::greedy_search(int in_pivot){
  int pivot = 0;
  n_state_flg2 = 0;
  for (size_t st_i=0; st_i < nstates; st_i++){
    state_flg[st_i] = 0;
  }

  while(n_state_flg2 < nstates){
    size_t max_qraw = 0;
    for (size_t st_i=0; st_i < nstates; st_i++){
      if ( state_qraw[st_i] >= max_qraw && state_flg[st_i] == 0) {
	max_qraw = state_qraw[st_i];
	pivot = st_i;
      }
    }
    greedy_search_pivot(pivot);
  }
}

int VirtualStateCoupling::greedy_search_pivot(int pivot){

  state_flg[pivot] = 2;
  n_state_flg2++;
  size_t cur_step = 0;
  state_adj_qw[pivot] = -1.0;
  cout << " greedy [pivot] " << pivot << " " << state_adj_qw[pivot] << endl;
  vector<size_t> todo_st;
  for (const auto& [st_j, subzones] : overlapping_subzones[pivot]){
    todo_st.push_back(st_j);
    state_flg[st_j] = 1;
  }
  while(!todo_st.empty()){
    size_t st_i = todo_st[0];
    todo_st.erase(todo_st.begin());
    if(state_qraw[st_i] > 0){
      state_flg[st_i] = 2;
      n_state_flg2++;
      continue;
    }
    double pi_factors = 1.0;
    int n_factors = 0;
    for (const auto& [st_j, subzones] : overlapping_subzones[st_i]){
      //cout << "gr " << st_j << " " << state_flg[st_j] << endl;
      if(state_flg[st_j] == 0 ){
	todo_st.push_back(st_j);
	state_flg[st_j] = 1;
      }else if(state_flg[st_j] == 2){
	for (const auto& sz : subzones ) {
	  std::vector<size_t> sz_crd_i(n_dim+1), sz_crd_j(n_dim+1);
	  sz_crd_i[0]=st_i;
	  sz_crd_j[0]=st_j;
	  for ( int i = 0; i < n_dim; i++){
	    sz_crd_i[i+1] = sz[i];
	    sz_crd_j[i+1] = sz[n_dim+i];
	  }
	  double qcano_j = state_adj_qw[st_j] + state_qraw[st_j] + state_qraw_is[sz_crd_j];
	  double qcano_i_uni = state_qraw[st_i] + state_qraw_is[sz_crd_i];
	  
	  //if (qcano_j >= 0 || qcano_i_uni >= 0){
	  if (state_adj_qw[st_j] > 0 || state_qraw[st_j] > 0 || state_qraw_is[sz_crd_j] > 0 ||
	      state_qraw[st_i] > 0 || state_qraw_is[sz_crd_i] > 0){
	    //cout << "i " << st_i << " " << st_j << " " << state_adj_qw[st_j] <<" " << state_qraw[st_j] << " "  << state_qraw_is[sz_crd_j] << " " << state_qraw[st_i] << " "  << state_qraw_is[sz_crd_i] <<endl;;
	    continue;
	  }

	  double factor = qcano_j - qcano_i_uni;
	  //cout << "scale " << st_i << " " << st_j << " ";
	  //for ( int i = 0; i < n_dim; i++)
	  //cout << sz_crd_i[i+1] << " " ;
	  //for ( int i = 0; i < n_dim; i++)
	  //cout << sz_crd_j[i+1] << " " ;
	  //cout << qcano_j << " " << qcano_i_uni << " "<< factor << " " 
	  //<< state_adj_qw[st_j] << " " <<  state_qraw[st_j] << " " <<  state_qraw_is[sz_crd_j]
	  //<< " " <<  state_qraw[st_i] << " " <<  state_qraw_is[sz_crd_i] << endl;
	  pi_factors += factor;
	  n_factors++;
	}
      }
    }
    if ( n_factors > 0 ) 
      state_adj_qw[st_i] += pi_factors/(double)n_factors;
    state_flg[st_i] = 2;
    n_state_flg2++;

    // normalize
    double sum_qw = 0.0;
    for ( int st = 0; st < nstates; st++){
      sum_qw += exp(state_adj_qw[st]);
    }
    for ( int st = 0; st < nstates; st++){
      state_adj_qw[st] = log(exp(state_adj_qw[st])/sum_qw);
    }
    cout << "greedy qw " << st_i << " " << pi_factors << " " << n_factors << " " << state_adj_qw[st_i] << endl;

  }

  return 0;
}

