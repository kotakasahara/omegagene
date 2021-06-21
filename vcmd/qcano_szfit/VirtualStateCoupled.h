#ifndef __VIRTUAL_STATE_COUPLED_H__
#define __VIRTUAL_STATE_COUPLED_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include "Config.h"

using namespace std;

class VirtualStateCoupling {
private:
  Config cfg;

  size_t nstates;
  //std::vector<double> kappa;
  double kappa;
  double default_weight;

  size_t n_dim;
  // the number of dimensions in the virtual system
  string fname_i_params;
  string fname_o_qcano;
  string fname_o_logstate;
  string fname_o_qrawstat;
  string fname_o_qrawstat_is;
  string fname_o_dbg;
  string fname_o_transit_count;
  int flg_transit_count;
  ofstream ofs_logstate;
  ofstream ofs_dbg;
  // output filename specified with LOGSTATE key

  int exchange_interval;
  // interval steps for virtual-state transitions

  float drift_mode;
  // 0  ... Usual VcMD, fulfilling the detailed balance.
  // >0 ... Drifting VcMD, without the detailed balance.

  int verbose;

  // properties for each state
  std::vector<std::vector<double> > lowers;
  std::vector<std::vector<double> > uppers;
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

  std::vector<size_t> state_qraw;
  // state_qraw[d] = qraw
  std::map< std::vector<size_t>, size_t > state_qraw_is;
  std::map< std::vector<size_t>, size_t > state_transit_count;
  // state_qraw_is[(vs1, vs2, ..., is1, is2, ...)] = qraw
  //   vs1, vs2, ... virtual state coordinate in 1-st, 2-nd, ... -th dimension.
  //   is1, is2, ... -1, 0, 1, or 2 corresponding to the lower or higher region from the state boundary
  //                    in each axis

  size_t state_qraw_sum;
  void parse_vstate(const string& fname);
  void parse_params(const string& fname);
  void init_transition_table();
  size_t conv_vstate_crd2id(std::vector<int> vcrd);
  std::vector<int> conv_vstate_id2crd(size_t v_id);
  bool is_in_range(size_t state, std::vector<double> arg);
  int trial_transition();
  std::string get_str_state_definition();
void set_cur_vsis(std::vector<int> v_crd, std::vector<double> arg);

  //i/o
  void write_qrawstat();
  void write_qrawstat_is();
  void write_transit_count();

public:
  explicit VirtualStateCoupling();
  ~VirtualStateCoupling();

  int setup(Config cfg);
  int subzonebased_mc();
};



#endif
