#include "Config.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

Config::Config(vector<string> &&arg) {
    set_arguments(std::forward<vector<string>>(arg));
    if (not fn_cfg.empty()) { set_arguments(extract_args_from_file(fn_cfg)); }
}

vector<string> Config::extract_args_from_file(const string &filepath) {
    ifstream ifs(filepath.c_str());
    if (not ifs) {
        cerr << "Cannot open " << filepath << "!\n";
        std::exit(1);
    }

    vector<string> args;
    string         buf;

    cout << "-----------------------------------\nReading from configuration file:  " << filepath << "\n";
    while (ifs && getline(ifs, buf)) {
      unsigned int pos1 = buf.find_first_of("#;");
      if (pos1 != string::npos) { buf = buf.substr(0, pos1); }
      if (buf.size() == 0) continue;
      cout << buf << endl;
      stringstream ss(buf);
      while (ss >> buf) { args.push_back(buf); }
    }
    cout << "-----------------------------------" << endl;
    return args;
}

void Config::set_arguments(std::vector<std::string> &&arg) {
    string type, val;
    for (auto itr = arg.begin(); itr != arg.end(); itr++) {
        if (*itr == "--mode") {
            itr++;
            if (*itr == "test") {
                mode = M_TEST;
            } else if (*itr == "md") {
                mode = M_DYNAMICS;
            } else {
                cerr << "invalid mode [" << (*itr) << "]\n";
                std::exit(1);
            }

        } else if (*itr == "--cfg") {
            fn_cfg = *++itr;

        } else if (*itr == "--inp") {
            fn_inp = *++itr;

        } else if (*itr == "--processor") {
            itr++;
            if (*itr == "single") {
                processor = PRCS_SINGLE;
            } else if (*itr == "mpi") {
                processor = PRCS_MPI;
            } else if (*itr == "cuda") {
                processor = PRCS_CUDA;
            } else if (*itr == "mpi-cuda") {
                processor = PRCS_MPI_CUDA;
            } else {
                processor = PRCS_DUMMY;
            }

        } else if (*itr == "--integrator") {
            itr++;
            if (*itr == "zhang") {
                integrator_type = INTGRTR_ZHANG;
            } else if (*itr == "leapfrog-presto") {
                integrator_type = INTGRTR_LEAPFROG_PRESTO;
            } else if (*itr == "velocity-verlet") {
                integrator_type = INTGRTR_VELOCITY_VERLET;
            } else if (*itr == "langevin") {
	      integrator_type = INTGRTR_LANGEVIN;
            } else if (*itr == "langevin-vv") {
	      integrator_type = INTGRTR_LANGEVIN_VV;
            } else if (*itr == "mc") {
	      integrator_type = INTGRTR_MC;
            } else {
                integrator_type = INTGRTR_DUMMY;
            }

        } else if (*itr == "--constraint") {
            itr++;
            if (*itr == "none") {
                constraint_type = CONST_NONE;
            } else if (*itr == "shake") {
                constraint_type = CONST_SHAKE;
            } else if (*itr == "shake-settle") {
                constraint_type = CONST_SHAKE_SETTLE;
            } else {
                constraint_type = CONST_DUMMY;
            }
        } else if (*itr == "--const-max-loops") {
            constraint_max_loops = atoi((*++itr).c_str());
        } else if (*itr == "--const-tolerance") {
            constraint_tolerance = atof((*++itr).c_str());
        } else if (*itr == "--thermo-const-max-loops") {
            thermo_const_max_loops = atoi((*++itr).c_str());
        } else if (*itr == "--thermo-const-tolerance") {
            thermo_const_tolerance = atof((*++itr).c_str());
        } else if (*itr == "--cutoff") {
            cutoff = atof((*++itr).c_str());
        } else if (*itr == "--n-steps") {
            n_steps = atoi((*++itr).c_str());
        } else if (*itr == "--time-step") {
            time_step = atof((*++itr).c_str());
        } else if (*itr == "--electrostatic") {
            itr++;
            if (*itr == "zero-dipole") {
                electrostatic = ELCTRST_ZERODIPOLE;
            } else if (*itr == "zero-quadrupole") {
                electrostatic = ELCTRST_ZEROQUADRUPOLE;
            } else if (*itr == "zero-octupole") {
	      electrostatic = ELCTRST_ZEROOCTUPOLE;
            } else if (*itr == "zero-hexadecapole") {
	      electrostatic = ELCTRST_ZEROHEXADECAPOLE;
            } else if (*itr == "debye-huckel") {
	      electrostatic = ELCTRST_DEBYE_HUCKEL;
            } else {
	      electrostatic = ELCTRST_DUMMY;
            }

        } else if (*itr == "--ele-alpha") {
            ele_alpha = atof((*++itr).c_str());

        } else if (*itr == "--thermostat") {
            itr++;
            if (*itr == "none") {
                thermostat_type = THMSTT_NONE;
            } else if (*itr == "scaling") {
                thermostat_type = THMSTT_SCALING;
            } else if (*itr == "hoover-evans") {
                thermostat_type = THMSTT_HOOVER_EVANS;
            } else {
                thermostat_type = THMSTT_NONE;
            }

        } else if (*itr == "--extended-ensemble") {
            itr++;
            if (*itr == "none") {
	      extended_ensemble = EXTENDED_NONE;
            } else if (*itr == "v-mcmd") {
	      extended_ensemble = EXTENDED_VMCMD;
            } else if (*itr == "v-aus") {
	      extended_ensemble = EXTENDED_VAUS;
            } else if (*itr == "vcmd") {
	      extended_ensemble = EXTENDED_VCMD;
            } else {
	      extended_ensemble = EXTENDED_NONE;
            }
	    
        } else if (*itr == "--temperature") {
	  temperature = atof((*++itr).c_str());

        } else if (*itr == "--berendsen-tau") {
	  berendsen_tau = atof((*++itr).c_str());
	  berendsen_tau *= 1000;
	  
        } else if (*itr == "--temperature-init") {
	  temperature_init = atof((*++itr).c_str());

        } else if (*itr == "--heating-steps") {
	  heating_steps = atoi((*++itr).c_str());
	  
        } else if (*itr == "--com-motion") {
	  itr++;
            if (*itr == "none") {
	      com_motion = COM_NONE;
            } else if (*itr == "cancel") {
	      com_motion = COM_CANCEL;
            }

        } else if (*itr == "--com-cancel-group-name") {
	  com_cancel_groups_name[n_com_cancel_groups_name] = ((*++itr).c_str());
	  // cout << "name " << com_cancel_groups_name[n_com_cancel_groups_name] << endl;
            n_com_cancel_groups_name++;

        } else if (*itr == "--com-cancel-group-id") {
            com_cancel_groups[n_com_cancel_groups] = atoi((*++itr).c_str());
            n_com_cancel_groups++;

        } else if (*itr == "--random-seed") {
            random_seed = atoi((*++itr).c_str());

        } else if (*itr == "--box-division") {
            for (auto &b : box_div) b = atoi((*++itr).c_str());

        } else if (*itr == "--nsgrid-cutoff") {
            nsgrid_cutoff = atof((*++itr).c_str());

        }

        // else if (*itr == "--nsgrid-min-width") { nsgrid_min_width= atof((*++itr).c_str()); }
        // else if (*itr == "--nsgrid-max-n-atoms") { nsgrid_max_n_atoms = atof((*++itr).c_str()); }
        else if (*itr == "--nsgrid-update-intvl") {
            nsgrid_update_intvl = atoi((*++itr).c_str());
        } else if (*itr == "--print-interval-coord") {
	  print_intvl_crd.push_back(atoi((*++itr).c_str()));

        } else if (*itr == "--print-interval-velo") {
            print_intvl_vel = atoi((*++itr).c_str());

        } else if (*itr == "--print-interval-force") {
            print_intvl_force = atoi((*++itr).c_str());

        } else if (*itr == "--print-interval-log") {
            print_intvl_log = atoi((*++itr).c_str());

        } else if (*itr == "--print-interval-energy") {
            print_intvl_energy = atoi((*++itr).c_str());

        } else if (*itr == "--print-interval-energyflow") {
            print_intvl_energyflow = atoi((*++itr).c_str());

        } else if (*itr == "--print-interval-extended-lambda") {
            print_intvl_extended_lambda = atoi((*++itr).c_str());
        } else if (*itr == "--print-interval-restart") {
	  print_intvl_restart = atoi((*++itr).c_str());
        } else if (*itr == "--fn-o-restart") {
            fn_o_restart = *++itr;
        } else if (*itr == "--fn-o-coord") {
	  fn_o_crd.push_back(*++itr);
	} else if (*itr == "--group-o-coord") {
	  //group_o_crd_name = ((*++itr).c_str());
	  group_o_crd_name.push_back(*++itr);
        } else if (*itr == "--format-o-coord") {
            itr++;
            if (*itr == "gromacs") {
                format_o_crd = CRDOUT_GROMACS;
            } else if (*itr == "presto") {
                format_o_crd = CRDOUT_PRESTO;
            } else {
                format_o_crd = CRDOUT_DUMMY;
            }

        } else if (*itr == "--fn-o-log") {
            fn_o_log = *++itr;

        } else if (*itr == "--fn-o-energy") {
            fn_o_energy = *++itr;

        } else if (*itr == "--fn-o-vmcmd-log") {
            fn_o_vmcmd_log = *++itr;

        } else if (*itr == "--fn-o-extended-lambda") {
            fn_o_extended_lambda = *++itr;

        } else if (*itr == "--format-o-extended-lambda") {
            itr++;
            if (*itr == "binary") {
                format_o_extended_lambda = LAMBDAOUT_BIN;
            } else if (*itr == "ascii") {
                format_o_extended_lambda = LAMBDAOUT_ASC;
            } else {
                format_o_extended_lambda = LAMBDAOUT_DUMMY;
            }
        }

        else if (*itr == "--fn-o-energyflow") {
            fn_o_energyflow = *++itr;

        } else if (*itr == "--dist-restraint") {
            itr++;
            if (*itr == "none") {
                dist_restraint_type = DISTREST_NONE;
            } else if (*itr == "harmonic") {
                dist_restraint_type = DISTREST_HARMONIC;
            } else {
                dist_restraint_type = DISTREST_DUMMY;
            }

        } else if (*itr == "--dist-restraint-weight") {
            dist_restraint_weight = atof((*++itr).c_str());

        } else if (*itr == "--position-restraint") {
            itr++;
            if (*itr == "none") {
	      pos_restraint_type = POSREST_NONE;
            } else if (*itr == "harmonic") {
	      pos_restraint_type = POSREST_HARMONIC;
            } else {
                pos_restraint_type = POSREST_DUMMY;
            }

        } else if (*itr == "--position-restraint-weight") {
            pos_restraint_weight = atof((*++itr).c_str());

        } else if (*itr == "--gpu-device-id") {
            gpu_device_id = atoi((*++itr).c_str());

        } else if (*itr == "--enhance-group-name") {
            enhance_groups_name[n_enhance_groups_name] = ((*++itr).c_str());
            n_enhance_groups_name++;

        } else if (*itr == "--enhance-sigma") {
            enhance_sigma = atof((*++itr).c_str());

        } else if (*itr == "--enhance-recovery-coef") {
            enhance_recov_coef = atof((*++itr).c_str());
        }
        //    else if (*itr == "--aus-type") {  aus_type = atoi((*++itr).c_str()); }
        else if (*itr == "--fn-o-aus-restart") {
            fn_o_aus_restart = *++itr;
	}else if (*itr == "--fn-o-vcmd-start") {
	  fn_o_vcmd_start = *++itr;
	}else if (*itr == "--fn-o-vcmd-q-raw") {
	  fn_o_vcmd_qraw = *++itr;
	}else if (*itr == "--fn-o-vcmd-q-raw-is") {
	  fn_o_vcmd_qraw_is = *++itr;
	}else if (*itr == "--begin-count-q-raw") {
	  begin_count_qraw = atoi((*++itr).c_str());
	  //}else if (*itr == "--default-q-raw") {
	  //extend_default_q_raw = atof((*++itr).c_str());
        } else if (*itr == "--vcmd-drift"){
	  vcmd_drift = atoi((*++itr).c_str());
        } else if (*itr == "--print-interval-group-com") {
	  print_intvl_group_com = atoi((*++itr).c_str());
        } else if (*itr == "--fn-o-group-com") {
	  fn_o_group_com = *++itr;
        } else if (*itr == "--debye-huckel-dielectric") {
	  dh_dielectric = atof((*++itr).c_str());
        } else if (*itr == "--debye-huckel-ionic-strength") {
	  dh_ionic_strength = atof((*++itr).c_str());
        } else if (*itr == "--debye-huckel-temperature") {
	  dh_temperature = atof((*++itr).c_str());
        } else if (*itr == "--nonbond") {
            itr++;
            if (*itr == "lennard-jones") {
	      nonbond = NONBOND_LJ;
            } else if (*itr == "hydrophobicity-scale-lj") {
	      nonbond = NONBOND_HPS;
	    }
        } else if (*itr == "--hydrophobicity-scale-epsiron") {
	  hps_epsiron = atof((*++itr).c_str());
	  cout << "  HPS EPS " << hps_epsiron << endl;
        } else if (*itr == "--expected-num-density") {
	  expected_num_density = atof((*++itr).c_str());
        } else if (*itr == "--langevin-gamma") {
	  langevin_gamma = atof((*++itr).c_str());
	  langevin_gamma /= 1000;
	  // convert the unit from ps^-1 to fs^-1
        } else if (*itr == "--testmc-delta-x") {
	  testmc_delta_x = atof((*++itr).c_str());
        } else {
	  stringstream ss;
	  error_exit(ss.str(), "1A00001");
        }
    }
    if (temperature_init < 0) temperature_init = temperature;
}
