#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "CelesteObject.h"
#include <array>
#include <vector>

struct Config : public CelesteObject {
    int         mode   = M_TEST;
    std::string fn_cfg = "md_i.cfg";
    std::string fn_inp = "md_i.inp";

    int processor     = PRCS_SINGLE;
    int gpu_device_id = -1;

    int  integrator_type        = INTGRTR_LEAPFROG_PRESTO;
    int  constraint_type        = CONST_NONE;
    real constraint_tolerance   = 0.000001;
    int  constraint_max_loops   = 1000;
    int  thermostat_type        = THMSTT_NONE;
    real temperature            = 300.0;
    real temperature_init       = -1.0;
    int  heating_steps          = 0;
    real thermo_const_tolerance = 0.000001;
    int  thermo_const_max_loops = 1000;

    real_pw     cutoff = 12.0;
    real_pw     cutoff_buf;
    int         n_steps             = 1;
    real        time_step           = 0.0005;
    int         electrostatic       = ELCTRST_ZERODIPOLE;
    real        ele_alpha           = 0.0;
    int         com_motion          = COM_NONE;
    int         n_com_cancel_groups = 0;
    int         com_cancel_groups[MAX_N_COM_GROUPS];
    int         n_com_cancel_groups_name = 0;
    std::string com_cancel_groups_name[MAX_N_COM_GROUPS];
    int         n_enhance_groups_name = 0;
    std::string enhance_groups_name[MAX_N_COM_GROUPS];
    int         random_seed       = -1;
    int         extended_ensemble = EXTENDED_NONE;

    std::array<int, 3> box_div = {{1, 1, 1}};

    int print_intvl_crd = 10000;
    int print_intvl_vel = 0;
    int print_intvl_log;
    int print_intvl_force = 0;
    int print_intvl_energy;
    int print_intvl_energyflow;
    int print_intvl_extended_lambda;

    std::string fn_o_restart     = "md_o.restart";
    std::string fn_o_crd         = "md_o.trr";
    std::string group_o_crd_name = "";
    std::string fn_o_log         = "md_o.log";
    std::string fn_o_energy      = "md_o.erg";
    std::string fn_o_vmcmd_log;
    std::string fn_o_extended_lambda;
    std::string fn_o_energyflow          = "md_o.efl";
    int         format_o_crd             = CRDOUT_GROMACS;
    int         format_o_extended_lambda = LAMBDAOUT_BIN;

    int init_vel_just = 0;
    // 0: initial velocity is 0-dt
    // 1: initial velocity is 0

    real_pw nsgrid_cutoff       = cutoff + 1.0;
    int     nsgrid_update_intvl = 1;
    // real nsgrid_min_width = cutoff * 0.5;
    // real nsgrid_max_n_atoms = 100;

    int         dist_restraint_type   = DISTREST_NONE;
    real        dist_restraint_weight = 0.0;
    int         pos_restraint_type    = POSREST_NONE;
    real        pos_restraint_weight  = 0.0;
    real        enhance_sigma         = 0.2;
    real        enhance_recov_coef    = 50;
    std::string fn_o_aus_restart      = "aus_restart_out.dat";
    // int  aus_type = AUSTYPE_MASSCENTER;

    Config() = default;
    Config(std::vector<std::string> &&arg);
    Config(const std::string &filepath) : Config(extract_args_from_file(filepath)) {}
    Config(int argc, char **argv) : Config(std::vector<std::string>(argv + 1, argv + argc)){};
    std::vector<std::string> extract_args_from_file(const std::string &filepath);
    void set_arguments(std::vector<std::string> &&arg);
};

#endif
