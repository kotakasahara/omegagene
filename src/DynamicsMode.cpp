#include "DynamicsMode.h"
DynamicsMode::DynamicsMode()
 : RunMode(){
} 

DynamicsMode::~DynamicsMode(){
} 

int DynamicsMode::test(Config* in_cfg){
  cfg = in_cfg;
  //RunMode::set_config_parameters(cfg);
  return 0;
}
int DynamicsMode::set_config_parameters(Config* in_cfg){
  /*
    Setting parameters from Config object.
   */
  //if(DBG>=1)
  //cout << "DBG1: DynamicsMode::set_config_parameters()"<<endl;
  //enecal = new EnergyCalc(&mmsys, &subbox);
  cfg = in_cfg;
  RunMode::set_config_parameters(cfg);
  thermostat = cfg->thermostat;
  temperature = cfg->temperature;
  integrator = cfg->integrator;
  time_step = cfg->time_step;
  nsgrid_cutoff = cfg->nsgrid_cutoff;
  return 0;
}
int DynamicsMode::initial_preprocess(){
  /*
    Preprocess. Setting up several constant parameters;
  */
  //if (DBG >= 1)
  cout << "DBG1: DynamicsMode::initial_preprocess()" << endl;
  mmsys.ff_setup(cfg);
  writer_trr.set_fn(cfg->fn_o_crd);
  writer_trr.open();
  temperature_coeff = 1.0 / (GAS_CONST * (real)mmsys.n_free) * JOULE_CAL * 1e3 * 2.0;
  //enecal->initial_preprocess();
  // set atom coordinates into PBC
  mmsys.revise_crd_inbox();
  cout << "subbox_setup() "<<cfg->box_div[0]<<" "
       << cfg->box_div[1] << " " << cfg->box_div[2] << endl;
  //grid
  subbox_setup();
  //cout << "nsgrid_setup" << endl;
  //subbox.nsgrid_set(nsgrid_cutoff);
  
  cout << "DBG MmSystem.n_bonds: " << mmsys.n_bonds << endl;

  return 0;
}

int DynamicsMode::terminal_process(){
  cout << "DynamicsMode::terminal_process()"<<endl;  
  writer_trr.close();
  return 0;
}

int DynamicsMode::main_stream(){
  //mmsys.nsgrid.move_crd_in_cell(0,0,0.0);
  for(mmsys.cur_step = 0;
      mmsys.cur_step <= n_steps;
      mmsys.cur_step++){
    //cout << "DBG 1 " << mmsys.cur_step << endl;
    calc_in_each_step();
  }
  //mmsys.output_ctimes();
  return 0;
}

int DynamicsMode::calc_in_each_step(){
  const clock_t startTimeStep = clock();

  const clock_t startTimeReset = clock();
  mmsys.cur_time = mmsys.cur_step * time_step;
  mmsys.reset_energy();
  const clock_t endTimeReset = clock();
  mmsys.ctime_cuda_reset_work_ene += endTimeReset - startTimeReset;

  const clock_t startTimeEne = clock();
  //cout << "calc_energy()" << endl;
  subbox.calc_energy();
  //cout << "gather_energies()"<<endl;
  gather_energies();
  const clock_t endTimeEne = clock();
  mmsys.ctime_calc_energy += endTimeEne - startTimeEne;

  if(cfg->expanded_ensemble == EXPAND_VMCMD){
    real potential_e = mmsys.pote_bond + mmsys.pote_angle
      + mmsys.pote_torsion + mmsys.pote_impro
      + mmsys.pote_14vdw + mmsys.pote_14ele
      + mmsys.pote_vdw + mmsys.pote_ele;
    subbox.expand_apply_bias(mmsys.cur_step, potential_e);
  }

  const clock_t startTimeVel = clock();
  //cout << "update_velocities"<<endl;
  subbox.update_velocities(1.0,
			   time_step);
  const clock_t endTimeVel = clock();
  mmsys.ctime_update_velo += endTimeVel - startTimeVel;

  if(cfg->thermostat==THMSTT_HOOVER_EVANS){
    subbox.thermo_hoover_evans(cfg->time_step,
			       mmsys.n_free,
			       cfg->temperature);
  }

  const clock_t startTimeCoord = clock();
  //cout << "update_coordinates"<<endl;  
  subbox.cpy_crd_prev();
  subbox.update_coordinates(time_step);

  if(cfg->constraint != CONST_NONE){
    if(mmsys.leapfrog_coef != 1.0){
      subbox.apply_constraint();
      mmsys.leapfrog_coef = 1.0;
    }else{
      subbox.apply_constraint();
      subbox.set_velocity_from_crd(time_step);
    }
  }

  //cout << "revise_coordinates"<<endl;  
  subbox.revise_coordinates_pbc();
  const clock_t endTimeCoord = clock();
  mmsys.ctime_update_coord += endTimeCoord - startTimeCoord;

  const clock_t startTimeKine = clock();
  subbox.velocity_average();
  subbox.set_vel_just(mmsys.vel_just);
  cal_kinetic_energy((const real**)mmsys.vel_just);
  const clock_t endTimeKine = clock();
  mmsys.ctime_calc_kinetic += endTimeKine - startTimeKine;

#ifndef F_WO_NS
  const clock_t startTimeHtod = clock();
  if(mmsys.cur_step%cfg->nsgrid_update_intvl==0){
    //cout << "nsgrid_update"<<endl;
    subbox.nsgrid_update();
  }else{
    subbox.nsgrid_crd_update();
    //revise_coordinates_pbc();
  }
  const clock_t endTimeHtod = clock();
  mmsys.ctime_cuda_htod_atomids += endTimeHtod - startTimeHtod;
#endif

  const clock_t endTimeStep = clock();
  mmsys.ctime_per_step += endTimeStep - startTimeStep;

  sub_output();
  sub_output_log();

  return 0;


}

int DynamicsMode::sub_output(){
  // Output
  //cout << "cur_step: " << mmsys.cur_step << " ";
  //cout << mmsys.cur_step % cfg->print_intvl_crd << endl;

  bool out_crd = cfg->print_intvl_crd > 0 && mmsys.cur_step % cfg->print_intvl_crd == 0;
  bool out_vel = cfg->print_intvl_vel > 0 && mmsys.cur_step % cfg->print_intvl_vel == 0;
  bool out_force = cfg->print_intvl_force > 0 && mmsys.cur_step % cfg->print_intvl_force == 0;
  
  writer_trr.write_trr(mmsys.n_atoms,
		       mmsys.cur_step, mmsys.cur_time,
		       mmsys.pbc.L[0], mmsys.pbc.L[1], mmsys.pbc.L[2],
		       mmsys.crd, mmsys.vel_just, mmsys.force,
		       true,
		       out_crd, out_vel, out_force);
  return 0;
}

int DynamicsMode::sub_output_log(){
  bool out_log = cfg->print_intvl_log > 0 && mmsys.cur_step % cfg->print_intvl_log == 0;
  stringstream ss;

  if(out_log){
    string strbuf;
    char buf[1024];
    sprintf(buf, "Step: %8lu    Time: %10.4f\n", mmsys.cur_step, mmsys.cur_time);
    ss << string(buf);
    real potential_e = mmsys.pote_bond + mmsys.pote_angle
      + mmsys.pote_torsion + mmsys.pote_impro
      + mmsys.pote_14vdw + mmsys.pote_14ele
      + mmsys.pote_vdw + mmsys.pote_ele;
    real total_e = potential_e + mmsys.kinetic_e;
    sprintf(buf, "Total:     %14.10e\n", total_e);
    ss << string(buf);
    sprintf(buf, "Potential: %14.10e    Kinetic:  %14.10e\n", potential_e, mmsys.kinetic_e);
    ss << string(buf);    
    sprintf(buf, "Bond:      %14.10e    Angle:    %14.10e\n", mmsys.pote_bond, mmsys.pote_angle);
    ss << string(buf);    
    sprintf(buf, "Torsion:   %14.10e    Improper: %14.10e\n", mmsys.pote_torsion, mmsys.pote_impro);
    ss << string(buf);    
    sprintf(buf, "14-VDW:    %14.10e    14-Ele:   %14.10e\n", mmsys.pote_14vdw, mmsys.pote_14ele);
    ss << string(buf);    
    sprintf(buf, "VDW:       %14.10e    Ele:      %14.10e\n", mmsys.pote_vdw, mmsys.pote_ele);
    ss << string(buf);    
    sprintf(buf, "Temperature:       %14.10e\n", mmsys.temperature);
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
  }
  cout << ss.str();
  return 0;
}
/*
int DynamicsMode::update_velocities(real firstcoeff){
  mmsys.velocity_swap();
  // real firstcoeff ...
  //   1.0 for the first step ??
  for(int atomid=0; atomid < mmsys.n_atoms; atomid++){
    for(int d=0; d < 3; d++){
      //cout << "atomid:" << atomid << " dim:"<<d<<endl;
      //cout << "time_step: " << time_step << " force_val: " << FORCE_VEL <<endl;
      //      cout << "mass: " << mmsys.mass[atomid] << endl;
      mmsys.vel_next[atomid][d] = mmsys.vel[atomid][d] - 
	(firstcoeff * time_step *
	 FORCE_VEL * mmsys.force[atomid][d] / mmsys.mass[atomid]);
    }
    //cout << "vel1 " << atomid << " " << mmsys.vel[atomid][0] << " , ";
    //cout << mmsys.vel[atomid][1] << " , ";
    //cout << mmsys.vel[atomid][2] << endl;
    //cout << "vel_next1 " << atomid << " " << mmsys.vel_next[atomid][0] << " , ";
    //cout << mmsys.vel_next[atomid][1] << " , ";
    //cout << mmsys.vel_next[atomid][2] << endl;
  }
  mmsys.velocity_average();
  return 0;
}
*/
int DynamicsMode::cal_kinetic_energy(const real** vel){
  real kine_pre = 0.0;
  for(int atomid=0; atomid < mmsys.n_atoms; atomid++){
    real kine_atom = 0.0;
    for(int d=0; d < 3; d++)
      kine_atom += vel[atomid][d] * vel[atomid][d];
    kine_pre += kine_atom * mmsys.mass[atomid]; 
    //cout << "dbg_kine: " << mmsys.mass[atomid] << " " << vel[atomid][0];
    //cout << " " << vel[atomid][1] <<" " << vel[atomid][2] << endl;
  }
  //cout << "ke_pre " << kine_pre << endl;
  mmsys.kinetic_e = kine_pre * KINETIC_COEFF;
  //cout << "ke " << mmsys.kinetic_e << endl;;
  mmsys.temperature = mmsys.kinetic_e * temperature_coeff;
  return 0;
}


int DynamicsMode::subbox_setup(){
  cout << "subbox.set_parameters" << endl;
  subbox.set_parameters(mmsys.n_atoms, &(mmsys.pbc), cfg,
			nsgrid_cutoff,
			cfg->box_div[0],
			cfg->box_div[1],
			cfg->box_div[2]);
  subbox.set_lj_param(mmsys.n_lj_types,
		      mmsys.lj_6term,
		      mmsys.lj_12term);
  //subbox.set_max_n_atoms_region();
  cout << "alloc_variables" << endl;
  subbox.alloc_variables();
  subbox.alloc_variables_for_bonds(mmsys.n_bonds);
  subbox.alloc_variables_for_angles(mmsys.n_angles);
  subbox.alloc_variables_for_torsions(mmsys.n_torsions);
  subbox.alloc_variables_for_impros(mmsys.n_impros);
  subbox.alloc_variables_for_nb14(mmsys.n_nb14);
  subbox.alloc_variables_for_excess(mmsys.n_excess);
  subbox.alloc_variables_for_nb15off(mmsys.max_n_nb15off);
  cout << "initial_division" << endl;
  subbox.initial_division((const real**)mmsys.crd,
			  (const real**)mmsys.vel_just,
			  (const real*)mmsys.charge,
			  (const real*)mmsys.mass,
			  (const int*)mmsys.atom_type);
  subbox_set_bonding_potentials();

  if(cfg->constraint != CONST_NONE){
    subbox.init_constraint(cfg->constraint,
			   cfg->constraint_max_loops,
			   cfg->constraint_tolerance,
			   mmsys.constraint.get_n_pair(),
			   mmsys.constraint.get_n_trio(),
			   mmsys.constraint.get_n_quad());
    
    int n_shake_dist = mmsys.constraint.get_n_pair() + 
      3 * mmsys.constraint.get_n_trio() + 
      6 * mmsys.constraint.get_n_quad();
    mmsys.n_free -= n_shake_dist;
    
    subbox.set_subset_constraint(mmsys.constraint);
  }
  if(cfg->expanded_ensemble == EXPAND_VMCMD){
    subbox.expand_init();
  }

  //  cout << "set_nsgrid" << endl;
#ifndef F_WO_NS
  subbox.set_nsgrid();
#endif



  //subbox.set_ff(&ff);

  return 0;
}
int DynamicsMode::subbox_set_bonding_potentials(){
  subbox.set_bond_potentials((const int**)mmsys.bond_atomid_pairs,
			     (const real*)mmsys.bond_epsiron,
			     (const real*)mmsys.bond_r0);
  subbox.set_angle_potentials((const int**)mmsys.angle_atomid_triads,
			      (const real*)mmsys.angle_epsiron,
			      (const real*)mmsys.angle_theta0);
  subbox.set_torsion_potentials((const int**)mmsys.torsion_atomid_quads,
				(const real*)mmsys.torsion_energy,
				(const int*)mmsys.torsion_overlaps,
				(const int*)mmsys.torsion_symmetry,
				(const real*)mmsys.torsion_phase,
				(const int*)mmsys.torsion_nb14);
  subbox.set_impro_potentials((const int**)mmsys.impro_atomid_quads,
			      (const real*)mmsys.impro_energy,
			      (const int*)mmsys.impro_overlaps,
			      (const int*)mmsys.impro_symmetry,
			      (const real*)mmsys.impro_phase,
			      (const int*)mmsys.impro_nb14);
  subbox.set_nb14_potentials((const int**)mmsys.nb14_atomid_pairs,
			     (const int**)mmsys.nb14_atomtype_pairs,
			     (const real*)mmsys.nb14_coeff_vdw,
			     (const real*)mmsys.nb14_coeff_ele);
  subbox.set_ele_excess((const int**)mmsys.excess_pairs);
  subbox.set_nb15off((const int*)mmsys.nb15off);
  return 0;
}
int DynamicsMode::gather_energies(){
  mmsys.pote_bond = subbox.get_pote_bond();
  mmsys.pote_angle = subbox.get_pote_angle();
  mmsys.pote_torsion = subbox.get_pote_torsion();
  mmsys.pote_impro = subbox.get_pote_impro();
  mmsys.pote_14vdw = subbox.get_pote_14vdw();
  mmsys.pote_14ele = subbox.get_pote_14ele();
  mmsys.pote_vdw = subbox.get_pote_vdw();
  mmsys.pote_ele = subbox.get_pote_ele();
  mmsys.pote_ele += mmsys.energy_self_sum;
  //cout << "mmsys.energy_self_sum : " << mmsys.energy_self_sum << endl;
  return 0;
}
