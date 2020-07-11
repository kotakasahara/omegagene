#include "DynamicsMode.h"

DynamicsMode::DynamicsMode() : RunMode() {}

DynamicsMode::~DynamicsMode() {
    if (DBG >= 1) cout << "DBG1 DynamicsMode::~DynamicsMode()" << endl;
}

int DynamicsMode::test(Config *in_cfg) {
    cfg = in_cfg;
    // RunMode::set_config_parameters(cfg);
    return 0;
}
int DynamicsMode::set_config_parameters(Config *in_cfg) {
    /*
      Setting parameters from Config object.
     */
    // if(DBG>=1)
    // cout << "DBG1: DynamicsMode::set_config_parameters()"<<endl;
    // enecal = new EnergyCalc(&mmsys, &subbox);
    cfg = in_cfg;
    RunMode::set_config_parameters(cfg);
    mmsys.extended_mode = cfg->extended_ensemble;
    if (cfg->extended_ensemble == EXTENDED_VMCMD) {
      mmsys.vmcmd = new ExtendedVMcMD();
      cout << "vmcmd new" << endl;
    } else if (cfg->extended_ensemble == EXTENDED_VAUS) {
      mmsys.vmcmd = new ExtendedVAUS();
      cout << "vaus new" << endl;
    } else if (cfg->extended_ensemble == EXTENDED_VCMD) {
      mmsys.vcmd = new ExtendedVcMD();
      cout << "vcmd new" << endl;
    }
    return 0;
}
int DynamicsMode::initial_preprocess() {
  /*
    Preprocess. Setting up several constant parameters;
  */
  // if (DBG >= 1)
  cout << "DBG1: DynamicsMode::initial_preprocess()" << endl;
  
  time_step      = cfg->time_step;
  time_step_half = cfg->time_step * 0.5;
  
  mmsys.ff_setup(cfg);
  int i=0;
  for ( auto itr = cfg->fn_o_crd.begin();
	itr != cfg->fn_o_crd.end(); itr++){
    if (cfg->format_o_crd == CRDOUT_GROMACS) {
      writer_trr.push_back(new WriteTrrGromacs());
      cout << "writer trr new " << i << endl;
    } else if (cfg->format_o_crd == CRDOUT_PRESTO) {
      writer_trr.push_back(new WriteTrrPresto());
      cout << "writer presto new " << i << endl;
    }
    writer_trr[i]->set_fn(*itr);
    writer_trr[i]->open();
    ++i;
  }
  if (cfg->constraint_type != CONST_NONE) {
    int n_shake_dist =
      mmsys.constraint.get_n_pair() + 3 * mmsys.constraint.get_n_trio() + 6 * mmsys.constraint.get_n_quad();
    mmsys.d_free -= n_shake_dist;
  }
  
  temperature_coeff = 1.0 / (GAS_CONST * (real)mmsys.d_free) * JOULE_CAL * 1e3 * 2.0;
  // enecal->initial_preprocess();
  // set atom coordinates into PBC
  mmsys.revise_crd_inbox();
  mmsys.set_atom_group_info(cfg);
  cout << "subbox_setup() " << cfg->box_div[0] << " " << cfg->box_div[1] << " " << cfg->box_div[2] << endl;
  // grid
  subbox_setup();
  // cout << "nsgrid_setup" << endl;
  // subbox.nsgrid_set(nsgrid_cutoff);
  
  // Random
  if (cfg->random_seed < 0) {
    random_device rnd;
    mmsys.set_random(rnd());
  } else {
    mmsys.set_random(cfg->random_seed);
  }
  if (cfg->extended_ensemble == EXTENDED_VMCMD ||
      cfg->extended_ensemble == EXTENDED_VAUS ) {
    cout << "V_McMD: " << endl;
    cout << "  VS log output ... " << cfg->fn_o_vmcmd_log << endl;
    cout << "  Lambda output ... " << cfg->fn_o_extended_lambda << endl;
    mmsys.vmcmd->set_files(cfg->fn_o_vmcmd_log, cfg->fn_o_extended_lambda, cfg->format_o_extended_lambda, cfg->fn_o_group_com);
    mmsys.vmcmd->set_lambda_interval(cfg->print_intvl_extended_lambda);
    mmsys.vmcmd->set_com_interval(cfg->print_intvl_group_com);
    mmsys.vmcmd->print_info();
    
    mmsys.vmcmd->set_params(&mmsys.random_mt, cfg->enhance_sigma, cfg->enhance_recov_coef,
			    cfg->n_steps); //, cfg->aus_type);
    
    // for(int i_grp=0; i_grp < mmsys.n_groups; i_grp++){
    // cout << "dbg1130 massDM " << i_grp << " " << mmsys.mass_inv_groups[i_grp]<<endl;
    //}
    mmsys.vmcmd->set_mass(subbox.get_mass(), mmsys.mass_groups, mmsys.mass_inv_groups);
    
  }else if(cfg->extended_ensemble == EXTENDED_VCMD){
    cout << "VcMD: " << endl;
    cout << "  VS log output ... " << cfg->fn_o_vmcmd_log << endl;
    cout << "  Lambda output ... " << cfg->fn_o_extended_lambda << endl;
    mmsys.vcmd->set_files(cfg->fn_o_vmcmd_log, cfg->fn_o_extended_lambda, cfg->format_o_extended_lambda,
			  cfg->fn_o_vcmd_qraw, cfg->fn_o_vcmd_start, cfg->fn_o_vcmd_qraw_is);
    mmsys.vcmd->set_lambda_interval(cfg->print_intvl_extended_lambda);
    mmsys.vcmd->print_info();
    cout << "print_info" << endl;
    mmsys.vcmd->set_params(&mmsys.random_mt, cfg->enhance_sigma, cfg->enhance_recov_coef,
			   cfg->n_steps, cfg->begin_count_qraw, cfg->vcmd_drift); //, cfg->aus_type);
    mmsys.vcmd->set_temperature(cfg->temperature);
    
    mmsys.vcmd->set_default_q_cano();
    
    // for(int i_grp=0; i_grp < mmsys.n_groups; i_grp++){
    // cout << "dbg1130 massDM " << i_grp << " " << mmsys.mass_inv_groups[i_grp]<<endl;
    //}
    cout << "set_mass"<< endl;
    mmsys.vcmd->set_mass(subbox.get_mass(), mmsys.mass_groups, mmsys.mass_inv_groups);
    cout << "vcmd"<<endl;
  }      
  
  // cout << "DBG MmSystem.n_bonds: " << mmsys.n_bonds << endl;
  
  //subbox.copy_vel(mmsys.vel_just);
  subbox.copy_vel_next(mmsys.vel_just);
  subbox.set_com_motion(mmsys.n_com_cancel_groups, mmsys.com_cancel_groups, mmsys.n_atoms_in_groups,
			mmsys.atom_groups, mmsys.mass_inv_groups);
  
  mmsys.print_com_cancel_groups();
  // mmsys.print_enhance_groups();
  mmsys.print_out_group();
  
  // subbox.set_com_motion(cfg->n_com_cancel_groups,
  // cfg->com_cancel_groups,
  // mmsys.n_atoms_in_groups,
  // mmsys.atom_groups,
  // mmsys.mass_inv_groups);
  
  cal_kinetic_energy((const real **)mmsys.vel_just);
  cout << "Initial kinetic energy : " << mmsys.kinetic_e << endl;
  cout << "Initial temperature : " << mmsys.temperature << endl;
  cout << "Degree of freedom : " << mmsys.d_free << endl;

  // for velocity-Verlet
  calc_energy_force();
  cout << "test_dynamicsmode" << endl;
  
  return 0;
}

int DynamicsMode::terminal_process() {
  cout << "DynamicsMode::terminal_process()" << endl;

  int i=0;
  for ( auto itr = writer_trr.begin();
	itr != writer_trr.end(); itr++){
    (*itr)->close();
    delete writer_trr[i];
    cout << "delete writer_trr " << i << endl;
    i++;
  }
  
  if (cfg->extended_ensemble == EXTENDED_VMCMD ||
      cfg->extended_ensemble == EXTENDED_VAUS ) {
    mmsys.vmcmd->close_files();
    delete mmsys.vmcmd;
  }else if(cfg->extended_ensemble == EXTENDED_VCMD){
    mmsys.vcmd->close_files();
    delete mmsys.vcmd;
  }

  cout << "term" << endl;
  return 0;
}

int DynamicsMode::main_stream() {
  mmsys.cur_step = 0;
  // while(mmsys.cur_step <= cfg->n_steps){
  for (mmsys.cur_step = 0; mmsys.cur_step < cfg->n_steps; mmsys.cur_step++) {
    sub_output();
    calc_in_each_step();
    if(cfg->print_intvl_restart > 0 && mmsys.cur_step % cfg->print_intvl_restart == 0){
      cout << "output restart " << mmsys.cur_step << endl;
      output_restart();
    }
    if (((cfg->print_intvl_log > 0 && mmsys.cur_step % cfg->print_intvl_log == 0) || mmsys.cur_step == 0
	 || mmsys.cur_step == cfg->n_steps - 1) && cfg->integrator_type != INTGRTR_MC) {
      sub_output_log();
    }
    mmsys.cur_time += cfg->time_step;
  }
  sub_output();
  output_restart();
  cout << "== The last step ==" << endl;
  calc_in_each_step();
  if (cfg->integrator_type == INTGRTR_LANGEVIN){
    output_restart();    
  }
  sub_output_log();
  return 0;
}

int DynamicsMode::output_restart() {
  if (cfg->integrator_type == INTGRTR_LANGEVIN){
    subbox.copy_crd_prev(mmsys.crd);
    //subbox.copy_vel_next(mmsys.vel_just);
  }else{
    subbox.copy_crd(mmsys.crd);
    subbox.copy_vel_next(mmsys.vel_just);
  }
  writer_restart.set_fn(cfg->fn_o_restart);
  writer_restart.write_restart(mmsys.n_atoms, (int)mmsys.cur_step, (double)mmsys.cur_time,
			       (double)(mmsys.pote_bond + mmsys.pote_angle + mmsys.pote_torsion + mmsys.pote_impro
					+ mmsys.pote_14vdw + mmsys.pote_14ele + mmsys.pote_vdw + mmsys.pote_ele),
			       (double)mmsys.kinetic_e, mmsys.crd, mmsys.vel_just);
  
  if (cfg->extended_ensemble == EXTENDED_VAUS){
    subbox.extended_write_aus_restart(cfg->fn_o_aus_restart, EXTENDED_VAUS);
  }else if(cfg->extended_ensemble == EXTENDED_VCMD) { 
    subbox.extended_write_aus_restart(cfg->fn_o_aus_restart, EXTENDED_VCMD);
  }
  
  return 0;
}

int DynamicsMode::calc_in_each_step() {
  return 0;
}
int DynamicsMode::apply_constraint() {

    // if(mmsys.leapfrog_coeff == 1.0){

    if (cfg->thermostat_type == THMSTT_SCALING) {
        subbox.apply_thermostat_with_shake(cfg->thermo_const_max_loops, cfg->thermo_const_tolerance);
        // mmsys.leapfrog_coef = 1.0 ;
    } else {
        subbox.apply_constraint();
    }

    //}
    return 0;
}

int DynamicsMode::apply_dist_restraint() {
  mmsys.pote_dist_rest = mmsys.dist_restraint->apply_restraint(mmsys.n_atoms, subbox.get_crds(), mmsys.pbc, mmsys.force);
  subbox.add_force_from_mmsys(mmsys.force);
  return 0;
}
int DynamicsMode::apply_pos_restraint() {
  mmsys.pote_pos_rest = mmsys.pos_restraint->apply_restraint(mmsys.n_atoms, subbox.get_crds(), mmsys.pbc, mmsys.force);
  subbox.add_force_from_mmsys(mmsys.force);
  return 0;
}

int DynamicsMode::sub_output() {
  // Output
  //cout << mmsys.cur_step % cfg->print_intvl_crd << endl;
  
  //bool out_vel = cfg->print_intvl_vel > 0 && mmsys.cur_step != 0 && (mmsys.cur_step) % cfg->print_intvl_vel == 0;
  //bool out_force =
  //cfg->print_intvl_force > 0 && mmsys.cur_step != 0 && (mmsys.cur_step) % cfg->print_intvl_force == 0;
  //if (out_vel) subbox.copy_vel(mmsys.vel_just);
  //if (out_crd || out_vel || out_force) {
  real total_e = mmsys.set_potential_e() + mmsys.kinetic_e;
  bool cp=false;
  int i_itr=-1;
  for ( auto itr = cfg->print_intvl_crd.begin();
	itr != cfg->print_intvl_crd.end(); itr++){
    i_itr++;
    if ((*itr) < 0 || mmsys.cur_step == 0 || (mmsys.cur_step) % (*itr) != 0)
      continue;
    if (!cp){
      subbox.copy_crd_prev(mmsys.crd);
      cp = true;
    }
    int outgrp = 0;
    if ( mmsys.out_group.size() > i_itr )
      outgrp = mmsys.out_group[i_itr];
    writer_trr[i_itr]->write_trr(mmsys.n_atoms,
				 (int)mmsys.cur_step,
				 mmsys.cur_time,
				 mmsys.pbc.L[0], mmsys.pbc.L[1], mmsys.pbc.L[2],
				 mmsys.crd, mmsys.vel_just, mmsys.force,
				 (float)mmsys.ctime_per_step / (float)CLOCKS_PER_SEC,
				 total_e, mmsys.kinetic_e,
				 mmsys.temperature,
				 mmsys.potential_e,
				 mmsys.pote_vdw, true, true, false, false,
				 mmsys.n_atoms_in_groups[outgrp],
				 mmsys.atom_groups[outgrp]);
				 //mmsys.n_atoms_in_groups[mmsys.out_group[i_itr]],
				 //mmsys.atom_groups[mmsys.out_group[i_itr]]);
  }
  return 0;
}

int DynamicsMode::sub_output_log() {
  stringstream ss;
  string       strbuf;
  char         buf[1024];
  sprintf(buf, "Step: %8lu    Time: %10.4f\n", mmsys.cur_step, mmsys.cur_time);
  ss << string(buf);
  real total_e = mmsys.set_potential_e() + mmsys.kinetic_e;
  sprintf(buf, "Total:     %14.10e\n", total_e);
  ss << string(buf);
  sprintf(buf, "Potential: %14.10e    Kinetic:  %14.10e\n", mmsys.potential_e, mmsys.kinetic_e);
  ss << string(buf);
  sprintf(buf, "Bond:      %14.10e    Angle:    %14.10e\n", mmsys.pote_bond, mmsys.pote_angle);
  ss << string(buf);
  sprintf(buf, "Torsion:   %14.10e    Improper: %14.10e\n", mmsys.pote_torsion, mmsys.pote_impro);
  ss << string(buf);
  sprintf(buf, "14-VDW:    %14.10e    14-Ele:   %14.10e\n", mmsys.pote_14vdw, mmsys.pote_14ele);
  ss << string(buf);
  sprintf(buf, "VDW:       %14.10e    Ele:      %14.10e\n", mmsys.pote_vdw, mmsys.pote_ele);
  ss << string(buf);
  if (cfg->dist_restraint_type != DISTREST_NONE) {
    sprintf(buf, "Distance restraint: %14.10e\n", mmsys.pote_dist_rest);
    ss << string(buf);
  }
  if (cfg->pos_restraint_type != POSREST_NONE) {
    sprintf(buf, "Position restraint: %14.10e\n", mmsys.pote_pos_rest);
    ss << string(buf);
  }
  if (cfg->extended_ensemble != EXTENDED_NONE){
    sprintf(buf, "Extended:           %14.10e\n", mmsys.pote_pos_rest);
    ss << string(buf);
  }
  sprintf(buf, "Temperature:       %14.10e\n", mmsys.temperature);
  ss << string(buf);
  sprintf(buf, "Comput Time:       %14.10e\n", (float)mmsys.ctime_per_step / (float)CLOCKS_PER_SEC);
  ss << string(buf);
  /*
    ss << "Step: " << mmsys.cur_step  << "\t";
    ss << "Time: " << mmsys.cur_time << " [ps]" << endl;
    ss << "Bond:  " << mmsys.pote_bond << "\t";
    ss << "Angle: " << mmsys.pote_angle << "\t";
    ss << "Torsion: " << mmsys.pote_torsion << "\t";
    ss << "Improper: " << mmsys.pote_impro << endl;
    ss << "14vdw:" << mmsys.pote_14vdw <<  "\t";
    ss << "14ele:" << mmsys.pote_14ele <<  "\t";
    ss << "Vdw:" << mmsys.pote_vdw <<  "\t";
    ss << "Ele:" << mmsys.pote_ele << endl;
    ss << "Kine:" << mmsys.kinetic_e << "\t";
    ss << "Temp:" << mmsys.temperature << endl;
    */
    cout << ss.str();
    return 0;
}
int DynamicsMode::cal_kinetic_energy(const real **vel) {
    real_fc kine_pre = 0.0;
    for (int atomid = 0; atomid < mmsys.n_atoms; atomid++) {
        real kine_atom = 0.0;
        for (int d = 0; d < 3; d++) kine_atom += vel[atomid][d] * vel[atomid][d];
        kine_pre += kine_atom * mmsys.mass[atomid];
	//cout << "dbg_kine: " << mmsys.mass[atomid] << " " << vel[atomid][0];
	//cout << " " << vel[atomid][1] <<" " << vel[atomid][2] << endl;
    }
    //cout << "ke_pre " << kine_pre << endl;
    mmsys.kinetic_e = kine_pre * KINETIC_COEFF;
    //cout << "ke " << mmsys.kinetic_e << " " << temperature_coeff << endl;;
    mmsys.temperature = mmsys.kinetic_e * temperature_coeff;
    return 0;
}

int DynamicsMode::subbox_setup() {
    // cout << "subbox.set_parameters" << endl;
  subbox.set_parameters(mmsys.n_atoms, &(mmsys.pbc), cfg, cfg->nsgrid_cutoff, cfg->box_div[0], cfg->box_div[1],
                          cfg->box_div[2]);
    subbox.set_lj_param(mmsys.n_lj_types, mmsys.lj_6term, mmsys.lj_12term,
			mmsys.lj_hps_cutoff, mmsys.lj_hps_lambda);
    // subbox.set_max_n_atoms_region();
    // cout << "alloc_variables" << endl;
    subbox.alloc_variables();

    

    subbox.alloc_variables_for_bonds(mmsys.n_bonds);
    subbox.alloc_variables_for_angles(mmsys.n_angles);
    subbox.alloc_variables_for_torsions(mmsys.n_torsions);
    subbox.alloc_variables_for_impros(mmsys.n_impros);
    subbox.alloc_variables_for_nb14(mmsys.n_nb14);
    subbox.alloc_variables_for_excess(mmsys.n_excess);
    subbox.alloc_variables_for_nb15off(mmsys.max_n_nb15off);
    subbox.initial_division(mmsys.crd, mmsys.vel_just, mmsys.charge, mmsys.mass, mmsys.atom_type);
    subbox_set_bonding_potentials();

    if (cfg->constraint_type != CONST_NONE) {
        subbox.init_constraint(cfg->constraint_type, cfg->constraint_max_loops, cfg->constraint_tolerance,
                               mmsys.constraint.get_n_pair(), mmsys.constraint.get_n_trio(),
                               mmsys.constraint.get_n_quad(), mmsys.settle.get_n_trio());

        subbox.set_subset_constraint(mmsys.constraint, mmsys.settle);
    }
    subbox.init_thermostat(cfg->thermostat_type, cfg->temperature_init, mmsys.d_free);

    if (cfg->extended_ensemble == EXTENDED_VCMD) {
      subbox.set_vcmd(cfg->extended_ensemble, mmsys.vcmd);
    }else if(cfg->extended_ensemble != EXTENDED_NONE) {
      subbox.set_extended(cfg->extended_ensemble, mmsys.vmcmd);
    }

    // cout << "set_nsgrid" << endl;
    subbox.revise_coordinates_pbc();

#ifndef F_WO_NS
    subbox.set_nsgrid();
#endif
    // subbox.set_ff(&ff);

    if(cfg->integrator_type == INTGRTR_LANGEVIN_VV ||
       cfg->integrator_type == INTGRTR_LANGEVIN || 
       cfg->integrator_type == INTGRTR_MC) {
      subbox.set_params_langevin(&mmsys.random_mt, cfg->langevin_gamma);
    }
      //cfg->langevin_gamma,
	//time_step,
	//cfg->temperature);
    //}

    return 0;
}
int DynamicsMode::subbox_set_bonding_potentials() {
    subbox.set_bond_potentials(mmsys.bond_atomid_pairs, mmsys.bond_epsiron, mmsys.bond_r0);
    subbox.set_angle_potentials(mmsys.angle_atomid_triads, mmsys.angle_epsiron, mmsys.angle_theta0);
    subbox.set_torsion_potentials(mmsys.torsion_atomid_quads, mmsys.torsion_energy, mmsys.torsion_overlaps,
                                  mmsys.torsion_symmetry, mmsys.torsion_phase, mmsys.torsion_nb14);
    subbox.set_impro_potentials(mmsys.impro_atomid_quads, mmsys.impro_energy, mmsys.impro_overlaps,
                                mmsys.impro_symmetry, mmsys.impro_phase, mmsys.impro_nb14);
    subbox.set_nb14_potentials(mmsys.nb14_atomid_pairs, mmsys.nb14_atomtype_pairs, mmsys.nb14_coeff_vdw,
                               mmsys.nb14_coeff_ele);
    subbox.set_ele_excess(mmsys.excess_pairs);
    subbox.set_nb15off(mmsys.nb15off);
    return 0;
}
int DynamicsMode::gather_energies() {
    mmsys.pote_bond    = subbox.get_pote_bond();
    mmsys.pote_angle   = subbox.get_pote_angle();
    mmsys.pote_torsion = subbox.get_pote_torsion();
    mmsys.pote_impro   = subbox.get_pote_impro();
    mmsys.pote_14vdw   = subbox.get_pote_14vdw();
    mmsys.pote_14ele   = subbox.get_pote_14ele();
    mmsys.pote_vdw     = subbox.get_pote_vdw();
    mmsys.pote_ele     = subbox.get_pote_ele();
    //real tmp           = mmsys.pote_ele;
    mmsys.pote_ele += mmsys.energy_self_sum;
    // cout << "tmp ele: " << tmp << " " << mmsys.energy_self_sum << " " << mmsys.pote_ele << endl;
    
    return 0;
}

//////////////////////

DynamicsModePresto::DynamicsModePresto() : DynamicsMode() {}

DynamicsModePresto::~DynamicsModePresto() {
  if (DBG >= 1) cout << "DBG1 DynamicsModePresto::~DynamicsModePresto()" << endl;  
}

int DynamicsModePresto::calc_in_each_step() {
  
    const clock_t startTimeStep  = clock();
    const clock_t startTimeReset = clock();

    mmsys.reset_energy();

    const clock_t endTimeReset = clock();
    mmsys.ctime_cuda_reset_work_ene += endTimeReset - startTimeReset;

#ifndef F_WO_NS
    const clock_t startTimeHtod = clock();
    if (mmsys.cur_step % cfg->nsgrid_update_intvl == 0) {
        subbox.nsgrid_update();
    } else {
#if defined(F_CUDA)
        subbox.nsgrid_crd_to_gpu();
#endif
    }
    const clock_t endTimeHtod = clock();
    mmsys.ctime_cuda_htod_atomids += endTimeHtod - startTimeHtod;
#endif
    const clock_t startTimeEne = clock();
    //cout <<"calc_energy()" <<endl;
    subbox.calc_energy(mmsys.cur_step);
     //cout << "gather_energies()"<<endl;
    gather_energies();

    if (cfg->dist_restraint_type != DISTREST_NONE || cfg->pos_restraint_type != POSREST_NONE) {
        subbox.copy_crd(mmsys.crd);
        if (cfg->dist_restraint_type != DISTREST_NONE) apply_dist_restraint();
        if (cfg->pos_restraint_type != POSREST_NONE) apply_pos_restraint();
    }

    const clock_t endTimeEne = clock();
    mmsys.ctime_calc_energy += endTimeEne - startTimeEne;
    if (cfg->extended_ensemble == EXTENDED_VMCMD) {
      subbox.extended_apply_bias(mmsys.cur_step, mmsys.set_potential_e());
    } else if (cfg->extended_ensemble == EXTENDED_VAUS) {
      subbox.extended_apply_bias_struct_param(mmsys.cur_step);
    } else if (cfg->extended_ensemble == EXTENDED_VCMD) {
      subbox.vcmd_apply_bias(mmsys.cur_step);
    }

    const clock_t startTimeVel = clock();
    // cout << "update_velocities"<<endl;
    subbox.cpy_vel_prev();
    subbox.update_velocities(cfg->time_step);
    const clock_t endTimeVel = clock();
    mmsys.ctime_update_velo += endTimeVel - startTimeVel;

    const clock_t startTimeCoord = clock();

    subbox.cancel_com_motion();

    // if(mmsys.leapfrog_coef == 1.0){
    if (cfg->thermostat_type == THMSTT_SCALING) {
        subbox.update_thermostat(mmsys.cur_step);
        if (cfg->constraint_type == CONST_NONE) {
            // mmsys.leapfrog_coef == 1.0){
            subbox.apply_thermostat();
        }
    }
    //}
    // cout << "update_coordinates"<<endl;
    subbox.cpy_crd_prev();
    subbox.update_coordinates_cur(cfg->time_step);

    if (cfg->constraint_type != CONST_NONE) { apply_constraint(); }
// cout << "revise_coordinates"<<endl;
#ifndef F_WO_NS
    subbox.update_coordinates_nsgrid();
#endif
    subbox.revise_coordinates_pbc();

    const clock_t endTimeCoord = clock();
    mmsys.ctime_update_coord += endTimeCoord - startTimeCoord;

    const clock_t startTimeKine = clock();
    // subbox.velocity_average();

    subbox.copy_vel_just(mmsys.vel_just);
    cal_kinetic_energy((const real **)mmsys.vel_just);
    const clock_t endTimeKine = clock();
    mmsys.ctime_calc_kinetic += endTimeKine - startTimeKine;

    const clock_t endTimeStep = clock();
    mmsys.ctime_per_step += endTimeStep - startTimeStep;

    return 0;
}

int DynamicsModePresto::apply_constraint() {

    // if(mmsys.leapfrog_coeff == 1.0){

    if (cfg->thermostat_type == THMSTT_SCALING) {
        subbox.apply_thermostat_with_shake(cfg->thermo_const_max_loops, cfg->thermo_const_tolerance);
        // mmsys.leapfrog_coef = 1.0 ;
    } else {
        subbox.apply_constraint();
    }

    //}
    return 0;
}

//////////////////////

DynamicsModeZhang::DynamicsModeZhang() : DynamicsMode() {}

DynamicsModeZhang::~DynamicsModeZhang() {}

int DynamicsModeZhang::calc_in_each_step() {
    const clock_t startTimeStep = clock();

    const clock_t startTimeReset = clock();
    mmsys.cur_time               = mmsys.cur_step * cfg->time_step;
    mmsys.reset_energy();
    const clock_t endTimeReset = clock();
    mmsys.ctime_cuda_reset_work_ene += endTimeReset - startTimeReset;

    subbox.update_coordinates_cur(time_step_half);
    subbox.cpy_vel_prev();

#ifndef F_WO_NS
    const clock_t startTimeHtod = clock();
    if (mmsys.cur_step % cfg->nsgrid_update_intvl == 0) {
        // cout << "nsgrid_update"<<endl;
        subbox.nsgrid_update();
    } else {
#if defined(F_CUDA)
        subbox.nsgrid_crd_to_gpu();
#endif
    }
    const clock_t endTimeHtod = clock();
    mmsys.ctime_cuda_htod_atomids += endTimeHtod - startTimeHtod;
#endif
    const clock_t startTimeEne = clock();
    // cout << "calc_energy()" << endl;
    subbox.calc_energy(mmsys.cur_step);
    // cout << "gather_energies()"<<endl;
    gather_energies();
    const clock_t endTimeEne = clock();
    mmsys.ctime_calc_energy += endTimeEne - startTimeEne;
    if (cfg->extended_ensemble == EXTENDED_VMCMD) {
      subbox.extended_apply_bias(mmsys.cur_step, mmsys.set_potential_e());
    } else if (cfg->extended_ensemble == EXTENDED_VAUS) {
      subbox.extended_apply_bias_struct_param(mmsys.cur_step);
    } else if (cfg->extended_ensemble == EXTENDED_VCMD) {
      subbox.vcmd_apply_bias(mmsys.cur_step);
    }
    
    if (cfg->dist_restraint_type != DISTREST_NONE) { apply_dist_restraint(); }
    if (cfg->pos_restraint_type != POSREST_NONE) { apply_pos_restraint(); }

    const clock_t startTimeVel = clock();
    subbox.cpy_crd_prev();
    // subbox.apply_thermostat();

    if (cfg->constraint_type != CONST_NONE) {
        // subbox.update_velocities(cfg->time_step);
        // vel_next
        // subbox.update_coordinates_cur(cfg->time_step);
        // apply_constraint();
        // subbox.set_force_from_velocity(cfg->time_step);
        // subbox.cpy_crd_from_prev();
        // subbox.update_velocities(cfg->time_step);
        // subbox.update_coordinates_cur(cfg->time_step);
        // subbox.apply_thermostat();
    }

    subbox.apply_thermostat();

    const clock_t endTimeVel = clock();
    mmsys.ctime_update_velo += endTimeVel - startTimeVel;

    const clock_t startTimeCoord = clock();
    subbox.update_coordinates_cur(time_step_half);
// cout << "revise_coordinates"<<endl;
#ifndef F_WO_NS
    subbox.update_coordinates_nsgrid();
#endif
    subbox.revise_coordinates_pbc();
    const clock_t endTimeCoord = clock();
    mmsys.ctime_update_coord += endTimeCoord - startTimeCoord;
    const clock_t startTimeKine = clock();
    // subbox.velocity_average();

    subbox.copy_vel_next(mmsys.vel_just);
    cal_kinetic_energy((const real **)mmsys.vel_just);
    const clock_t endTimeKine = clock();
    mmsys.ctime_calc_kinetic += endTimeKine - startTimeKine;

    const clock_t endTimeStep = clock();
    mmsys.ctime_per_step += endTimeStep - startTimeStep;
    // test output 0707
    if ((cfg->print_intvl_log > 0 && mmsys.cur_step % cfg->print_intvl_log == 0) || mmsys.cur_step == 0){
      mmsys.set_potential_e();
      cout << "DBG0707c " << mmsys.cur_step  << " " 
	   << subbox.get_crds()[0]  << " " 
	   << subbox.get_crds()[1]  << " " 
	   << subbox.get_crds()[2]  << " " 
	   << mmsys.potential_e  <<  endl;	
    }

    return 0;
}

int DynamicsModeZhang::apply_constraint() {
    subbox.apply_constraint();
    return 0;
}
//////////////////////

DynamicsModeVelocityVerlet::DynamicsModeVelocityVerlet() : DynamicsMode() {}
DynamicsModeVelocityVerlet::~DynamicsModeVelocityVerlet() {}

int DynamicsMode::calc_energy_force() {
    const clock_t startTimeReset = clock();

    mmsys.reset_energy();

    const clock_t endTimeReset = clock();
    mmsys.ctime_cuda_reset_work_ene += endTimeReset - startTimeReset;
#ifndef F_WO_NS
    const clock_t startTimeHtod = clock();
    if (mmsys.cur_step % cfg->nsgrid_update_intvl == 0) {
        subbox.nsgrid_update();
    } else {
#if defined(F_CUDA)
        subbox.nsgrid_crd_to_gpu();
#endif
    }
    const clock_t endTimeHtod = clock();
    mmsys.ctime_cuda_htod_atomids += endTimeHtod - startTimeHtod;
#endif
    const clock_t startTimeEne = clock();
    subbox.calc_energy(mmsys.cur_step);
    // cout << "gather_energies()"<<endl;
    gather_energies();
    if (cfg->dist_restraint_type != DISTREST_NONE || cfg->pos_restraint_type != POSREST_NONE) {
        subbox.copy_crd(mmsys.crd);
        if (cfg->dist_restraint_type != DISTREST_NONE) apply_dist_restraint();
        if (cfg->pos_restraint_type != POSREST_NONE) apply_pos_restraint();
    }
    return 0;
}
int DynamicsModeVelocityVerlet::calc_in_each_step() {
    const clock_t startTimeStep  = clock();

    const clock_t startTimeCoord = clock();
    // if(mmsys.leapfrog_coef == 1.0){
    //if (cfg->thermostat_type == THMSTT_SCALING) {
    //subbox.update_thermostat(mmsys.cur_step);
    //if (cfg->constraint_type == CONST_NONE) {
    // mmsys.leapfrog_coef == 1.0){
    //subbox.apply_thermostat();
    //}
    //}

    //}
    // cout << "update_coordinates"<<endl;
    
    subbox.cpy_crd_prev();
    subbox.update_coordinates_vv(cfg->time_step);
    subbox.cpy_work_prev();
    calc_energy_force();

    if (cfg->constraint_type != CONST_NONE) { apply_constraint(); }
    cout << "revise_coordinates"<<endl;
#ifndef F_WO_NS
    subbox.update_coordinates_nsgrid();
#endif
    subbox.revise_coordinates_pbc();

    const clock_t endTimeCoord = clock();
    mmsys.ctime_update_coord += endTimeCoord - startTimeCoord;

    const clock_t startTimeVel = clock();
     cout << "update_velocities"<<endl;
    subbox.cpy_vel_prev();
    subbox.update_velocities_vv(cfg->time_step);
    const clock_t endTimeVel = clock();
    mmsys.ctime_update_velo += endTimeVel - startTimeVel;

    subbox.cancel_com_motion();

    const clock_t startTimeKine = clock();
    subbox.copy_vel(mmsys.vel_just);
    cal_kinetic_energy((const real **)mmsys.vel_just);
    const clock_t endTimeKine = clock();
    mmsys.ctime_calc_kinetic += endTimeKine - startTimeKine;

    const clock_t endTimeStep = clock();
    mmsys.ctime_per_step += endTimeStep - startTimeStep;

    return 0;
}

int DynamicsModeVelocityVerlet::apply_constraint() {
    // if(mmsys.leapfrog_coeff == 1.0){

    if (cfg->thermostat_type == THMSTT_SCALING) {
        subbox.apply_thermostat_with_shake(cfg->thermo_const_max_loops, cfg->thermo_const_tolerance);
        // mmsys.leapfrog_coef = 1.0 ;
    } else {
        subbox.apply_constraint();
    }

    //}
    return 0;
}

///////

DynamicsModeLangevin::DynamicsModeLangevin() : DynamicsMode() {}

DynamicsModeLangevin::~DynamicsModeLangevin() {
  if (DBG >= 1) cout << "DBG1 DynamicsModeLangevin::~DynamicsModeLangevin()" << endl;  
}

int DynamicsModeLangevin::calc_in_each_step() {

    const clock_t startTimeStep  = clock();
    mmsys.reset_energy();

#ifndef F_WO_NS
    const clock_t startTimeHtod = clock();
    if (mmsys.cur_step % cfg->nsgrid_update_intvl == 0) {
        subbox.nsgrid_update();
    } else {
#if defined(F_CUDA)
        subbox.nsgrid_crd_to_gpu();
#endif
    }
#endif

    //cout <<"calc_energy()" <<endl;
    subbox.calc_energy(mmsys.cur_step);
     //cout << "gather_energies()"<<endl;
    gather_energies();

    if (cfg->dist_restraint_type != DISTREST_NONE || cfg->pos_restraint_type != POSREST_NONE) {
      //subbox.copy_crd(mmsys.crd);
      if (cfg->dist_restraint_type != DISTREST_NONE) apply_dist_restraint();
      if (cfg->pos_restraint_type != POSREST_NONE) apply_pos_restraint();
    }
    
    if (cfg->extended_ensemble == EXTENDED_VMCMD) {
      subbox.extended_apply_bias(mmsys.cur_step, mmsys.set_potential_e());
    } else if (cfg->extended_ensemble == EXTENDED_VAUS) {
      subbox.extended_apply_bias_struct_param(mmsys.cur_step);
    } else if (cfg->extended_ensemble == EXTENDED_VCMD) {
      subbox.vcmd_apply_bias(mmsys.cur_step);
    }
    
    subbox.cpy_crd_prev2();
    if (mmsys.cur_step == 0){
      subbox.cpy_vel_prev();
      //subbox.update_coordinates_from_vel(time_step);
      subbox.update_coordinates_cur(time_step);
    }else{
      subbox.update_coordinates_langevin(time_step_half, cfg->langevin_gamma, cfg->temperature, mmsys.cur_step);
      //subbox.cpy_vel_prev();
      //subbox.update_velocities(time_step_half, cfg->langevin_gamma, cfg->temperature);
      subbox.set_velocities_just_langevin(time_step);
      subbox.copy_vel(mmsys.vel_just);
      cal_kinetic_energy((const real **)mmsys.vel_just);
    }
    //subbox.cancel_com_motion();

    // if(mmsys.leapfrog_coef == 1.0){
    if (cfg->thermostat_type == THMSTT_SCALING) {
        subbox.update_thermostat(mmsys.cur_step);
        if (cfg->constraint_type == CONST_NONE) {
            // mmsys.leapfrog_coef == 1.0){
            subbox.apply_thermostat();
        }
    }
    //if (cfg->thermostat_type == THMSTT_SCALING) {
    //subbox.update_thermostat(mmsys.cur_step);
    //if (cfg->constraint_type == CONST_NONE) {
    //subbox.apply_thermostat();
    //}
    //}
    
    //if (cfg->constraint_type != CONST_NONE) { apply_constraint(); }
// cout << "revise_coordinates"<<endl;

#ifndef F_WO_NS
    subbox.update_coordinates_nsgrid();
#endif
    subbox.revise_coordinates_pbc();

    const clock_t endTimeStep = clock();
    mmsys.ctime_per_step += endTimeStep - startTimeStep;

    // test output 0707
    if ((cfg->print_intvl_log > 0 && mmsys.cur_step % cfg->print_intvl_log == 0) || mmsys.cur_step == 0){
      mmsys.set_potential_e();
      cout << "DBG0707 " << mmsys.cur_step  << " " 
	   << subbox.get_crds()[0]  << " " 
	   << subbox.get_crds()[1]  << " " 
	   << subbox.get_crds()[2]  << " " 
	   << mmsys.potential_e  <<  endl;	
    }

    return 0;
}

int DynamicsModeLangevin::apply_constraint() {

    // if(mmsys.leapfrog_coeff == 1.0){

    if (cfg->thermostat_type == THMSTT_SCALING) {
        subbox.apply_thermostat_with_shake(cfg->thermo_const_max_loops, cfg->thermo_const_tolerance);
        // mmsys.leapfrog_coef = 1.0 ;
    } else {
        subbox.apply_constraint();
    }

    //}
    return 0;
}

/////////////////

DynamicsModeLangevinVV::DynamicsModeLangevinVV() : DynamicsMode() {}

DynamicsModeLangevinVV::~DynamicsModeLangevinVV() {
  if (DBG >= 1) cout << "DBG1 DynamicsModeLangevinVV::~DynamicsModeLangevinVV()" << endl;  
}

int DynamicsModeLangevinVV::calc_in_each_step() {
    const clock_t startTimeStep  = clock();
    mmsys.reset_energy();
    
#ifndef F_WO_NS
    const clock_t startTimeHtod = clock();
    if (mmsys.cur_step % cfg->nsgrid_update_intvl == 0) {
        subbox.nsgrid_update();
    } else {
#if defined(F_CUDA)
        subbox.nsgrid_crd_to_gpu();
#endif
    }
#endif

    subbox.calc_energy(mmsys.cur_step);
     //cout << "gather_energies()"<<endl;
    gather_energies();

    if (cfg->dist_restraint_type != DISTREST_NONE || cfg->pos_restraint_type != POSREST_NONE) {
      //subbox.copy_crd(mmsys.crd);
      if (cfg->dist_restraint_type != DISTREST_NONE) apply_dist_restraint();
      if (cfg->pos_restraint_type != POSREST_NONE) apply_pos_restraint();
    }

    if (cfg->extended_ensemble == EXTENDED_VMCMD) {
      subbox.extended_apply_bias(mmsys.cur_step, mmsys.set_potential_e());
    } else if (cfg->extended_ensemble == EXTENDED_VAUS) {
      subbox.extended_apply_bias_struct_param(mmsys.cur_step);
    } else if (cfg->extended_ensemble == EXTENDED_VCMD) {
      subbox.vcmd_apply_bias(mmsys.cur_step);
    }

    if(mmsys.cur_step > 0){
      subbox.cpy_vel_prev();
      subbox.update_velocities_langevin_vv_second(time_step_half, cfg->langevin_gamma, cfg->temperature);
    }
    subbox.copy_vel_next(mmsys.vel_just);
    cal_kinetic_energy((const real **)mmsys.vel_just);
    subbox.cpy_vel_prev();
    subbox.update_velocities_langevin_vv_first(time_step_half, cfg->langevin_gamma, cfg->temperature);

    subbox.cancel_com_motion();
    // if(mmsys.leapfrog_coef == 1.0){
    if (cfg->thermostat_type == THMSTT_SCALING) {
        subbox.update_thermostat(mmsys.cur_step);
        if (cfg->constraint_type == CONST_NONE) {
            // mmsys.leapfrog_coef == 1.0){
            subbox.apply_thermostat();
        }
    }
    //if (cfg->thermostat_type == THMSTT_SCALING) {
    //subbox.update_thermostat(mmsys.cur_step);
    //if (cfg->constraint_type == CONST_NONE) {
    //subbox.apply_thermostat();
    //}
    //}
    
    //if (cfg->constraint_type != CONST_NONE) { apply_constraint(); }
// cout << "revise_coordinates"<<endl;
    subbox.cpy_crd_prev();
    subbox.update_coordinates_cur(time_step);

#ifndef F_WO_NS
    subbox.update_coordinates_nsgrid();
#endif
    subbox.revise_coordinates_pbc();

    const clock_t endTimeStep = clock();
    mmsys.ctime_per_step += endTimeStep - startTimeStep;

    return 0;
}

int DynamicsModeLangevinVV::apply_constraint() {

    // if(mmsys.leapfrog_coeff == 1.0){

    if (cfg->thermostat_type == THMSTT_SCALING) {
        subbox.apply_thermostat_with_shake(cfg->thermo_const_max_loops, cfg->thermo_const_tolerance);
        // mmsys.leapfrog_coef = 1.0 ;
    } else {
        subbox.apply_constraint();
    }

    //}
    return 0;
}

////

DynamicsModeMC::DynamicsModeMC() : DynamicsMode() {}

DynamicsModeMC::~DynamicsModeMC() {
  if (DBG >= 1) cout << "DBG1 DynamicsModeMC::~DynamicsModeMC()" << endl;  
  cout << "Accepted : " << mmsys.n_acc << endl;
}

int DynamicsModeMC::calc_in_each_step() {
  const clock_t startTimeStep = clock();

  if(mmsys.cur_step > 0){
    subbox.cpy_crd_prev();
    subbox.testmc_trial_move(cfg->testmc_delta_x);
#ifndef F_WO_NS
    subbox.update_coordinates_nsgrid();
#endif
    subbox.revise_coordinates_pbc();
  }    

  mmsys.cpy_energy_to_prev();
  mmsys.reset_energy();

  subbox.calc_energy(mmsys.cur_step);
  gather_energies();

  if (cfg->dist_restraint_type != DISTREST_NONE || cfg->pos_restraint_type != POSREST_NONE) {
    subbox.copy_crd(mmsys.crd);
    if (cfg->dist_restraint_type != DISTREST_NONE) apply_dist_restraint();
    if (cfg->pos_restraint_type != POSREST_NONE) apply_pos_restraint();
  }

  if (cfg->extended_ensemble == EXTENDED_VCMD) {
    mmsys.pote_extend  = subbox.vcmd_apply_bias(mmsys.cur_step);
  }
  mmsys.set_potential_e();
  //// Metropolis
  bool flg_accept = true;
  real delta_e = mmsys.potential_e + mmsys.pote_extend - (mmsys.potential_e_prev + mmsys.pote_extend_prev);
  real rnd = 0;
  real prob = 0;
  if(delta_e > 0) {
    rnd = mmsys.random_mt();
    prob = exp(-delta_e/(GAS_CONST/JOULE_CAL * 1e-3 * mmsys.temperature));
    if (rnd  > prob ) flg_accept = false;
  }
  if(!flg_accept && mmsys.cur_step > 0){
    subbox.cpy_crd_from_prev();    
    mmsys.cpy_energy_from_prev();
#ifndef F_WO_NS
    subbox.update_coordinates_nsgrid();
#endif
  } else{
    mmsys.n_acc ++;
  }

  if ((cfg->print_intvl_log > 0 && mmsys.cur_step % cfg->print_intvl_log == 0) || mmsys.cur_step == 0){  
    cout << "DBG0707b " << mmsys.cur_step  << " " 
	 << subbox.get_crds()[0]  << " " 
	 << subbox.get_crds()[1]  << " " 
      	 << subbox.get_crds()[2]  << " " 
	 << mmsys.potential_e  <<  " " 
	 << mmsys.pote_extend  <<  " "
	 << delta_e << " " << rnd << " " << prob;

  if(!flg_accept && mmsys.cur_step > 0){
    cout << " rej";
  }else{
    cout << " acc";
  }
  cout << endl;	
  
  const clock_t endTimeStep = clock();
  mmsys.ctime_per_step += endTimeStep - startTimeStep;

  return 0;
}
