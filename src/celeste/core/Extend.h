#ifndef __EXTENDED_H__
#define __EXTENDED_H__

#include "CelesteObject.h"

#include "PBC.h"
#include "Write.h"
#include "WriteTrr.h"

#include "celeste/random/Random.h"

#include <cmath>

class Extended : public CelesteObject {
 private:
 protected:
  int write_lambda_interval;
  int write_com_interval;
 public:
  Extended();
  ~Extended();
  void set_lambda_interval(int in_lambda_interval);
  void set_com_interval(int in_com_interval);
};


class VirtualState : public CelesteObject {
 private:
 protected:
  real  lambda_range[2];
  real  trans_prob[2];
  int   poly_order;
  real *poly_params;
  real  alpha[2];
  
 public:
  VirtualState();
  ~VirtualState();
  int set_order(int in_order);
  int get_order() { return poly_order; };
  int set_poly_param(int ord, real param);
  real get_poly_param(int ord) { return poly_params[ord]; };
  int set_params(real in_lambda_low, real in_lambda_high, real in_prob_low, real in_prob_high);
  real get_lambda_low() { return lambda_range[0]; };
  real get_lambda_high() { return lambda_range[1]; };
  int set_alpha(real in_alhpa_low, real in_alpha_high);
  bool is_in_range(real lambda);
  real get_trans_prob(int up_down) { return trans_prob[up_down]; }
};

class ExtendedVMcMD : public Extended {
 private:
 protected:
  unsigned long n_steps;
  
  int  n_vstates;
    int  trans_interval;
    real temperature;
    real const_k;

    VirtualState *vstates;
    int           init_vs;
    int           random_seed;
    bool          flg_vs_transition;

    int cur_vs;

    WriteTTPVMcMDLog writer_vslog;

    WriteTableLog *writer_lambda;
    WriteTableLog *writer_com;

    real_pw *mass;
    real_pw  mass_sum;
    real_pw *mass_groups;
    real_pw *mass_groups_inv;
    real     sigma;
    real     recov_coef;
    // real sigma_half;
    // real sigma_sq_inv;
    //
    int *            n_atoms_in_groups;
    int **           atom_groups;
    int              n_enhance_groups;
    std::vector<int> enhance_groups;
    int              n_enhance_group_pairs;
    int **           enhance_group_pairs;

    real ***crd_groups;
    // crd_groups[group][atom][xyz]
    real **crd_centers;
    // crd_centers[group][xyz]
    real **unit_vec;
    // unit_vec[group][xyz]

    int                      aus_type;
    celeste::random::Random *random_mt;
    // uniform_real_distribution<float> random_gen;

  public:
    ExtendedVMcMD();
    ~ExtendedVMcMD();

    int alloc_crd_centers();
    int free_crd_centers();
    int set_n_vstates(int in_n_vstates);
    void set_trans_interval(int in_trans_interval);
    void set_temperature(real in_tmp);
    int get_trans_interval();
    int get_temperature();
    int apply_bias(unsigned long cur_step, real in_lambda, real_fc *work, int n_atoms_box);

    VirtualState &get_vstate(int vs_id) { return vstates[vs_id]; };

    int  get_init_vs() { return init_vs; };
    void set_init_vs(int in_init_vs) {
        init_vs = in_init_vs;
        cur_vs  = init_vs;
    };
    int  get_random_seed() { return random_seed; };
    void set_random_seed(int in_seed) { random_seed = in_seed; };
    int trial_transition(int source, int rel_dest, real lambda);
    virtual real cal_struct_parameters(real *crd, PBC *pbc);
    int set_current_vstate(real lambda);
    virtual int scale_force(real lambda, real_fc *work, int n_atoms);
    // files
    int set_files(std::string fn_vslog, std::string fn_lambda, int format_lambda,
		  std::string fn_com);
    int close_files();
    int write_vslog(int cur_steps);
    int write_lambda(real lambda);
    int write_com();

    void enable_vs_transition() { flg_vs_transition = true; }
    int set_vs_order(int vs_id, int ord);
    int set_vs_params(int  vs_id,
                      real lambda_low,
                      real lambda_high,
                      real prob_low,
                      real prob_high,
                      real alpha_low,
                      real alpha_high);
    int set_vs_poly_param(int vs_id, int ord, real param);
    int print_info();

    virtual int set_enhance_groups(int *            in_n_atoms_in_groups,
				   int **           in_atom_groups,
				   int              in_n_enhance_groups,
				   std::vector<int> in_enhance_groups);
    int set_mass(real_pw *in_mass, real_pw *in_mass_groups, real_pw *in_mass_groups_inv);
    int set_params(celeste::random::Random *in_mt, real in_sigma, real in_recov_coef, int in_n_steps);
    void set_aus_type(int in_aus_type) { aus_type = in_aus_type; };

    int write_aus_restart(std::string fn_out);
    real ***get_crd_groups() { return crd_groups; };
};

class ExtendedVAUS : public ExtendedVMcMD {
  private:
  protected:
    // int n_vstates;
    // int trans_interval;
    // real temperature;
    // real const_k;

    // VirtualState *vstates;
    // int init_vs;
    // int random_seed;
    // bool flg_vs_transition;

    // int cur_vs;

    // WriteTTPVMcMDLog writer_vslog;

    // WriteTableLog* writer_lambda;

    // from MmSystem

  public:
    ExtendedVAUS();
    ~ExtendedVAUS();

    real set_crd_centers(real *crd, PBC *pbc);
    int set_init_crd_groups(real *crd);
    virtual real cal_struct_parameters(real *crd, PBC *pbc);

    // int set_n_vstates(int in_n_vstates);
    // void set_trans_interval(int in_trans_interval);
    // void set_temperature(real in_tmp);
    // int get_trans_interval();
    // int get_temperature();
    // int apply_bias(unsigned long cur_step,
    // real in_lambda,
    // real_fc* work,
    // int n_atoms_box);

    // VirtualState& get_vstate(int vs_id){ return vstates[vs_id]; };

    // int get_init_vs(){ return init_vs; };
    // void set_init_vs(int in_init_vs){ init_vs = in_init_vs; cur_vs = init_vs;};
    // int get_random_seed(){ return random_seed; };
    // void set_random_seed(int in_seed){ random_seed = in_seed; };
    // int trial_transition(int source, int rel_dest,
    // real lambda);

    // int set_current_vstate(real lambda);
    virtual int scale_force(real lambda, real_fc *work, int n_atoms);

    // files
    // int set_files(string fn_vslog, string fn_lambda, int format_lambda);
    // int close_files();
    // int write_vslog(int cur_steps);
    // int write_lambda(real lambda);

    // void enable_vs_transition(){flg_vs_transition = true;}
    // int set_vs_order(int vs_id, int ord);
    // int set_vs_params(int vs_id,
    // real lambda_low, real lambda_high,
    // real prob_low, real prob_high,
    // real alpha_low, real alpha_high);
    // int set_vs_poly_param(int vs_id, int ord, real param);
    // int print_info();
};

class ExtendedVcMD : public Extended {
 private:
 protected:
  unsigned long  n_steps;
  int  random_seed;
  int  trans_interval;
  bool          flg_vs_transition;
  real temperature;
  real const_k;
  
    WriteTTPVMcMDLog writer_vslog;
    WriteTableLog *writer_lambda;
    //WriteVcMDParam writer_qcano;
    WriteVcMDParam writer_qraw;
    WriteTableLogAscii writer_start;

    unsigned long begin_count_q_raw;

    int reactcrd_type;
    real_pw *mass;
    real_pw  mass_sum;
    real_pw *mass_groups;
    real_pw *mass_groups_inv;
    real     sigma;
    real     recov_coef;
    // real sigma_half;
    // real sigma_sq_inv;
    //
    int *            n_atoms_in_groups;
    int **           atom_groups;

    real ***crd_groups;
    // crd_groups[group][atom][xyz]
    real **crd_centers;
    // crd_centers[group][xyz]
    real **unit_vec;
    // unit_vec[group][xyz]

    int              n_enhance_groups;
    std::vector<int> enhance_groups;
    std::vector< std::vector<int> > enhance_group_pairs;

    int                      aus_type;
    celeste::random::Random *random_mt;
    // uniform_real_distribution<float> random_gen;

    int  n_dim;
    
    std::vector< std::vector<real> > vc_range_min;
    std::vector< std::vector<real> > vc_range_max;
    // range_min[dim][vs]
    // range_max[dim][vs]
    std::vector<int> init_vs;
    // init_vs[dim] = vs
    std::vector< std::vector<int> > grp_ids;
    std::vector< std::vector<std::string> > grp_names;
    std::map< std::vector<int>, real > q_cano;
    std::map< std::vector<int>, real > q_raw;
    real default_q_cano;
    
    std::vector<int> cur_vs;
    
    std::vector<real> lambda;

    std::vector< std::vector<int> > vs_next;
    // vs_next[0,1,2,3] = [x1,y1],[x1,y2],[x2,y1],[x2,y2]
    // the list of destination state for virtual state transition trials
    
    
 public:
    ExtendedVcMD();
    ~ExtendedVcMD();

    int alloc_crd_centers();
    int free_crd_centers();

    //void set_default_q_raw(real in_qr) { default_q_raw = in_qr; };
    int set_default_q_cano();
    void set_reactcrd_type(int in_type) { reactcrd_type = in_type; };
    void set_trans_interval(int in_trans_interval);
    void set_temperature(real in_tmp);
    int  get_random_seed() { return random_seed; };
    void set_random_seed(int in_seed) { random_seed = in_seed; };
    void set_n_dim(int in_n_dim);
    int close_files();
    int get_n_dim(){return n_dim;}
    void enable_vs_transition() { flg_vs_transition = true; }
    int push_vs_range(std::vector<real> new_min,
		      std::vector<real> new_max);
    void set_q_cano(std::map< std::vector<int>, real > in_q);
    void push_grp_ids_name(std::vector<int> in_ids,
			   std::vector<std::string> in_names){
     grp_ids.push_back(in_ids); 
     grp_names.push_back(in_names);
    }
    int set_params(celeste::random::Random *in_mt, real in_sigma,
		   real in_recov_coef, int in_n_steps,
		   int in_begin_count_q_raw);
    real set_crd_centers(real *crd, PBC *pbc);
    bool is_in_range();
    int apply_bias(unsigned long cur_step,
		   real_fc *work, int n_atoms_box);
    int scale_force(real_fc *work, int n_atoms);
    int set_files(std::string fn_vslog, std::string fn_lambda, int format_lambda,
		  std::string fn_qcano, std::string fn_qraw);
    
    std::vector<int> get_init_vs() { return init_vs; };
    void set_init_vs(std::vector<int> in_vs) {
      init_vs = in_vs;
      cur_vs = in_vs;
    };
    int set_vs_next();
    int set_vs_next_sub(const std::vector<int> tmp_vs_next,
			const std::vector< std::vector<int> > vs_next_crd);
    int trial_transition();


    int write_vslog(int cur_steps);
    int write_lambda();
    int write_q();

    int print_info();
    int set_struct_parameters(real *crd, PBC *pbc);
    int set_struct_parameters_mass_center(real *crd, PBC *pbc);
    int set_struct_parameters_min(real *crd, PBC *pbc);
    int set_enhance_groups(int *            in_n_atoms_in_groups,
			   int **           in_atom_groups,
			   int              in_n_enhance_groups,
			   std::vector<int> in_enhance_groups);

    real ***get_crd_groups() { return crd_groups; };
    int set_mass(real_pw *in_mass, real_pw *in_mass_groups, real_pw *in_mass_groups_inv);
    int write_aus_restart(std::string fn_out);
};

#endif
