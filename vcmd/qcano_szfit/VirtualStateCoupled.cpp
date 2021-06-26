#include "VirtualStateCoupled.h"
  
VirtualStateCoupling::VirtualStateCoupling()
{
  default_weight = 1.0;

}  
  
VirtualStateCoupling::~VirtualStateCoupling(){

}

std::string VirtualStateCoupling::get_str_state_definition(){
  std::stringstream ss;
  ss << exchange_interval << std::endl;
  ss << n_dim << std::endl;
  for ( int d = 0; d < n_dim; d ++){
    ss << lowers_vaxis[d].size() << std::endl;
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
    printf("%s\n", buf.c_str());
    unsigned int pos1 = buf.find_first_of("#;");
    if (pos1 != string::npos) { buf = buf.substr(0, pos1); }
    if (buf.size() == 0) continue;
    stringstream ss(buf);
    while (ss >> buf) { args.push_back(buf);
      printf("%s\n", buf.c_str());
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
  parse_params_state_definition(&ifs);
  parse_params_qweight(&ifs);
  ifs.close();
}

void VirtualStateCoupling::parse_params_state_definition(ifstream* ifs){
  vector<string> args;
  string         buf;
  string         cur, cur1, cur2;
  getline(*ifs, buf);
  exchange_interval = atoi(buf.c_str());
  printf("exchange_interval: %d \n", exchange_interval);
  getline(*ifs, buf);
  n_dim = atoi(buf.c_str());
  printf("n_dim: %d \n", n_dim);
  // read lower and upper bounds for each state in each axis
  for (size_t c_dim=0; c_dim < n_dim; c_dim++){
    int cur_nstates; 
    getline(*ifs, buf);
    stringstream ss(buf);
    ss >> cur;  cur_nstates = atoi(cur.c_str());
    printf("%d-th dim : %d states\n", c_dim, cur_nstates);
    vector<double> buf_lowers;
    vector<double> buf_uppers;
    for(int c_state = 0; c_state < cur_nstates; c_state++){
      getline(*ifs, buf);
      stringstream ss2(buf);
      //ss2 << buf;
      ss2 >> cur1 >> cur2;
      printf("dbg param read: %s %s\n", cur1.c_str(), cur2.c_str());
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
  state_qraw = vector<size_t>(nstates);
  state_qraw_sum = 0;
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
    printf("dbg read param: v_id %ld %d ... [%lf ... %lf]\n",v_id, v_crd[0], uppers[v_id][0], lowers[v_id][0]);
    state_weights[v_id] = -1.0;
    printf("v_id %d (", v_id);
    for(int d=0; d<n_dim; d++){
      printf("%d ", v_crd[d]);
    }
    printf(")\n");
  }
  // read information about parameters
}
void VirtualStateCoupling::parse_params_qweight(ifstream* ifs){
  vector<string> args;
  string         buf;
  string         cur, cur1, cur2;
  
  while(*ifs && getline(*ifs, buf)){
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
    if(flg_default){
      default_weight = cur_param;
      continue;
    }
    
    size_t v_id = conv_vstate_crd2id(v_crd);
    state_weights[v_id] = cur_param;
  }
  //
  for(size_t v_id=0; v_id < nstates; v_id++){
    if(state_weights[v_id] < 0.0){
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
void VirtualStateCoupling::parse_params_qraw_is(ifstream* ifs){
  vector<string> args;
  string         buf;
  string         cur, cur1, cur2;
  
  while(*ifs && getline(*ifs, buf)){
    size_t pos1 = buf.find_first_of("#;");
    if (pos1 != string::npos) { buf = buf.substr(0, pos1); }
    if (buf.size() == 0) continue;

    if(buf=="END") break;
    vector<int> v_crd(n_dim);

    stringstream ss(buf);

    bool flg_default = false;
    for (size_t c_dim=0; c_dim < n_dim*2; c_dim++){
      ss >> cur;
      v_crd[c_dim] = atoi(cur.c_str())-1;
      if(atoi(cur.c_str())-1 < 0) flg_default=true;
    }
    ss >> cur;
    double cur_param = atof(cur.c_str()) ;
    if(flg_default){
      default_weight = cur_param;
      continue;
    }
    
    size_t v_id = conv_vstate_crd2id(v_crd);
    state_weights[v_id] = cur_param;
  }
  //
  for(size_t v_id=0; v_id < nstates; v_id++){
    if(state_weights[v_id] < 0.0){
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


void VirtualStateCoupling::init_transition_table()
{
  // FIXME TODO: this loop requires O(nstate^2), may be a problem for nstate >= 10^6
  // One workaround is to input transition information as well from the external script.
  transition_candidates.clear();
  for(size_t i = 0; i < nstates; ++i) {
    printf("dbg: transition from state %ld : ", i);
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
	printf("%ld ", j);
      }

    }
    printf("\n");
    for(size_t d = 0; d < n_dim; ++d) {
      printf("[%lf ... %lf] \n", lowers[i][d], uppers[i][d]);
    }
  }
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
/*
void VirtualStateCoupling::write_qrawstat(){
  std::string state_defs = get_str_state_definition();
  ofstream ofs_qrawstat(fname_o_qrawstat);
  ofs_qrawstat << state_defs;
  for ( int l = 0; l < nstates; l++){
    std::vector<int> vs_crd = conv_vstate_id2crd(l);
    for ( const auto vsc : vs_crd )
      // (+1) convert to 1-origin integer
      ofs_qrawstat << vsc+1 << " ";
    ofs_qrawstat << state_qraw[l] << std::endl;
  }
  ofs_qrawstat << "END" << std::endl;
  ofs_qrawstat.close();
  return;
}
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
  fname_i_params = cfg.fname_i_params;
  fname_o_qcano = cfg.fname_o_qcano;
  parse_params(fname_i_params);
  init_transition_table();
  {
    size_t tottransition = 0;
    for(size_t i = 0; i < nstates; ++i) {
      tottransition += transition_candidates[i].size();
    }
    printf("  Total No. of possible transitions: %ld\n", tottransition);
  }

  return 0;
}

int VirtualStateCoupling::subzonebased_mc(){
  return 0;
}
