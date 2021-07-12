#ifndef __VIRTUAL_STATE_COUPLED_H__
#define __VIRTUAL_STATE_COUPLED_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <algorithm>
#include "Config.h"

using namespace std;
//const double VirtualSystemCoupling::EPS = 1e-10;

class VirtualStateCoupling {
private:
  /////  Config //////
  string fname_i_params;
  string fname_i_qraw_is;
  string fname_o_qcano;
  string fname_o_logstate;
  string fname_o_qrawstat;
  string fname_o_qweight_opt;
  string fname_o_dbg;
  string fname_o_transit_count;
  int flg_transit_count;

  double target_error;

  double mc_temp;
  double mc_delta_x;
  double mc_delta_x_max;
  size_t mc_steps;
  int mc_log_interval;
  double mc_target_acc_ratio;
  int mc_acc_duration;

  int greedy_max_steps;
  size_t greedy_pivot;

  //anneal
  int mc_n_window_trend;
  int mc_error_ave_window_size;
  double mc_max_temp;
  double mc_min_temp;
  double mc_delta_temp;


  size_t nstates;
  //std::vector<double> kappa;
  double kappa;
  double default_weight;

  size_t n_dim;
  // the number of dimensions in the virtual system

  ofstream ofs_logstate;
  ofstream ofs_dbg;
  // output filename specified with LOGSTATE key

  int exchange_interval;
  // interval steps for virtual-state transitions

  float drift_mode;
  // 0  ... Usual VcMD, fulfilling the detailed balance.
  // >0 ... Drifting VcMD, without the detailed balance.

  int verbose;
  int qweight_write_mode;
  // QW_FILE_MODE_RAW = 0 ... raw value
  // QW_FILE_MODE_LOG = 1 ... log value
  int process_unsampled_zone;
  // UNSAMPLED_OMIT = 0
  // UNSAMPLED_MIN = 1


  // properties for each state

  std::vector< std::vector<std::string> > aus_groups;
  std::vector< std::vector<double> > lowers;
  std::vector< std::vector<double> > uppers;
  // lower[v_id][d] ... lower boundary of d-th reaction coordinates of v_id-th state
  //    v_id is in {0, 1, 2, ... nstates-1}

  std::vector<double> state_weights;
  // state_weights[v_id] ... weigth of v_id-th state.
  std::vector< std::vector<size_t> > transition_candidates;
  // transition_candidates[v_id] = [x,y,z, ...]
  //    x, y, z, ... are in {1, 2, ..., nstates}
  //       which can be transitioned from the state v_id.

  // properties for each virtual axis
  std::vector<int> nstates_vaxis;
  // nstates_vaxis[d] ... the number of virtual states in d-th dimension
  std::vector<std::vector<double> > lowers_vaxis;
  std::vector<std::vector<double> > uppers_vaxis;
  // lowers_vaxis[d][i] ... lower boundary of d-th reaction coordinate of i-th state in d-th dimension
  // upper_vaxis[d][i] ... upper boundary ...

  //std::random_device rnd;
//std::default_random_engine generator;
  // current state
  size_t current_state;
  std::vector<size_t> cur_vsis;
  // cur_vsis is the current state in vs coordinate
  // cur_vsis[(vs1, vs2, ..., is1, is2, ...)]

  std::vector<double> state_qraw;
  // state_qraw[d] = qraw
  std::map< std::vector<size_t>, double > state_qraw_is;
  std::map< std::vector<size_t>, size_t > state_transit_count;
  // state_qraw_is[(vs_id, is1, is2, ...)] = qraw
  //   vs_id ... virtual state Id.
  //   is1, is2, ... -1, 0, 1, or 2 corresponding to the lower or higher region from the state boundary
  //                    in each axis

  std::vector<double> state_adj_qw;
  // Scaling factor for each state
  // state_adj_qw[v_id] = q_weight
  // 
  std::vector<double> state_adj_qw_opt;
  
  std::vector<double> state_adj_qw_error;
  double total_error;
  double opt_error;
  size_t mc_acc;

  // 

  size_t state_qraw_sum;

  vector< map<size_t, vector< vector<int> > > > overlapping_subzones;
  // overlapping_subzones[vs_id][vs_id, [[0,0], [0,1], ...]]

  int n_state_flg2;
  std::vector<int> state_flg;
  // flag for greedy search


  void parse_vstate(const string& fname);
  void parse_params(const string& fname);
  void parse_qraw_is(const string& fname);  
  int parse_params_state_definition(ifstream &ifs);
  void parse_params_qweight(ifstream &ifs, int param_mode);
  void parse_params_qraw_is(ifstream &ifs);
  void init_transition_table();
  void init_data();
  size_t conv_vstate_crd2id(std::vector<int> vcrd);
  std::vector<int> conv_vstate_id2crd(size_t v_id);
  bool is_in_range(size_t state, std::vector<double> arg);
  void write_qweight(std::string fname, std::vector<double> in_qw, bool def_val, int param_mode);
  int trial_transition();
  std::string get_str_state_definition(int param_mode);
  void set_cur_vsis(std::vector<int> v_crd, std::vector<double> arg);

  //i/o
  void write_qrawstat();
  void write_qrawstat_is();
  void write_transit_count();

public:
  explicit VirtualStateCoupling();
  ~VirtualStateCoupling();

  int setup(Config cfg);
  int enum_overlapping_subzones();
  int enum_overlapping_subzones_sub(std::vector<int>& overlap_positions,
				    std::vector<int>& cur_vec,
				    std::vector< std::vector<int> >& subzones);
  int print_overlapping_subzones();
  double calc_qraw_error_all();
  double calc_qraw_error(size_t st_i, double qw_i, bool skip_lower_id);
  
  int mode_test();
  int mode_subzonebased_mc();
  int mc_loop();
  int greedy_search(int in_pivot);
  int greedy_search_pivot(int in_pivot);

};



#endif
