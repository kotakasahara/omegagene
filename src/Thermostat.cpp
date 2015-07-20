#include "Thermostat.h"

ThermostatObject::ThermostatObject()
  : CelesteObject(){
}
ThermostatObject::~ThermostatObject(){
}
int ThermostatObject::set_time_step(real in_time_step){
  time_step = in_time_step;
  time_step_inv_sq = 1.0/(time_step*time_step);
  return 0;
}
int ThermostatObject::set_temperature_coeff(int in_d_free){
  d_free = in_d_free;
  temperature_coeff = (2.0 * JOULE_CAL * 1e+3) / (GAS_CONST * (real)d_free) * KINETIC_COEFF;
  return 0;
}
int ThermostatObject::set_constant(int n_atoms, real_pw* mass_inv, real* vel, real* force){
  return 0;
}
int ThermostatObject::apply_thermostat(const int n_atoms,
				       real_fc* work,
				       real* vel, real* vel_next,
				       real_pw* mass,
				       real_pw* mass_inv){

  return 0;
}
int ThermostatObject::apply_thermostat_with_shake(int n_atoms,
						  real_fc* work,
						  real* crd, real* crd_prev,
						  real* vel, real* vel_next,
						  real_pw* mass,
						  real_pw* mass_inv,
						  ConstraintObject* constraint,
						  PBC* pbc,
						  real* buf_crd,
						  const int max_loops,
						  const real tolerance,
						  COMMotion* commotion,
						  int* atomids_rev){

  return 0;
}


////////////////////////////////////////////////////////////


ThermostatScaling::ThermostatScaling() 
  : ThermostatObject(){
}
ThermostatScaling::~ThermostatScaling(){
}
int ThermostatScaling::set_constant(int n_atoms, real_pw* mass_inv, real* vel, real* force){
  return 0;
}

int ThermostatScaling::apply_thermostat(int n_atoms,
					real_fc* work,
					real* vel, real* vel_next,
					real_pw* mass,
					real_pw* mass_inv){
				       
  real_fc kine_pre = 0.0;
  
  for(int i=0, i_3=0; i < n_atoms; i++, i_3+=3){
    real vel_norm = 0;
    for(int d=0; d < 3; d++){
      real vel_tmp = vel[i_3+d] + 0.5 * 
	(time_step * work[i_3+d] * mass_inv[i]);
      vel_norm += vel_tmp * vel_tmp;
    }
    kine_pre += mass[i] * vel_norm;
  }
  
  real dt_temperature = kine_pre * temperature_coeff; 
  real scale = sqrt(temperature / dt_temperature);
  for(int i=0, i_3=0; i < n_atoms; i++, i_3+=3){
    for(int d=0; d < 3; d++){
      real vel_diff = time_step * work[i_3+d] * mass_inv[i];
      vel_next[i_3+d] = (2.0 * scale - 1.0) * vel[i_3+d] + scale * vel_diff;
    }
  }
  return 0;
}

int ThermostatScaling::apply_thermostat_with_shake(int n_atoms,
						   real_fc* work,
						   real* crd, real* crd_prev,
						   real* vel, real* vel_next,
						   real_pw* mass,
						   real_pw* mass_inv,
						   ConstraintObject* constraint,
						   PBC* pbc,
						   real* buf_crd,
						   const int max_loops,
						   const real tolerance,
						   COMMotion* commotion,
						   int* atomids_rev){
  //DBG
  bool converge = false;
  for(int idx = 0; idx < n_atoms*3; idx++)
    buf_crd[idx] =  crd[idx];
  for(int i_loop=0; i_loop < max_loops; i_loop++){

    //cout << "DBG 01 : " << i_loop << endl;
    constraint->apply_constraint(crd, crd_prev, mass_inv, pbc);

    real_fc kine_pre = 0.0;
    for(int i_atom = 0, i_atom_3 = 0;
	i_atom < n_atoms; i_atom++, i_atom_3+=3){
      real vel_norm = 0.0;
      for(int d=0; d < 3; d++){
	//vel_next[i_atom_3+d] = (crd[i_atom_3+d] - crd_prev[i_atom_3+d]) / cfg->time_step;
	real vel_tmp = (vel_next[i_atom_3+d] + vel[i_atom_3+d]) * 0.5;
	vel_norm += (vel_tmp * vel_tmp);
      }
      kine_pre += vel_norm * mass[i_atom];
    }
    real cur_temperature = kine_pre * temperature_coeff;
    
    if(i_loop == max_loops-1) break;
    real diff = fabs(cur_temperature - temperature) / temperature;
    if(diff < tolerance){
      converge = true;
      break;
    }


    kine_pre = 0.0;
    for(int i_atom = 0, i_atom_3 = 0;
	i_atom < n_atoms; i_atom++, i_atom_3+=3){
      real vel_norm = 0.0;
      for(int d=0; d < 3; d++){
	//buf_crd2[i_atom_3+d] += (crd[i_atom_3+d] - buf_crd1[i_atom_3+d]) * time_step_inv_sq;
	//real vel_diff = -FORCE_VEL * work[i_atom_3+d] * mass_inv[i_atom] + buf_crd2[i_atom_3+d];
	real vel_diff = work[i_atom_3+d] * mass_inv[i_atom] + (crd[i_atom_3+d]-buf_crd[i_atom_3+d])*time_step_inv_sq;
	//real vel_tmp = vel[i_atom_3+d] + 0.5 * cfg->time_step * vel_diff;
	vel_next[i_atom_3+d] = vel[i_atom_3+d] + 0.5 * time_step * vel_diff;
	real vel_tmp = vel_next[i_atom_3+d];
	vel_norm += vel_tmp * vel_tmp;
      }
      kine_pre += mass[i_atom] * vel_norm;
    }

    cur_temperature = kine_pre * temperature_coeff;

    real scale = sqrt(temperature / cur_temperature);

    kine_pre = 0.0;
    for(int i_atom = 0, i_atom_3 = 0;
	i_atom < n_atoms; i_atom++, i_atom_3+=3){
      real vel_norm = 0.0;
      for(int d=0; d < 3; d++){
	//real vel_diff = -FORCE_VEL * work[i_atom_3+d] * mass_inv[i_atom] + buf_crd[i_atom_3+d];
	real vel_diff = work[i_atom_3+d] * mass_inv[i_atom] + (crd[i_atom_3+d]-buf_crd[i_atom_3+d])*time_step_inv_sq;
	vel_next[i_atom_3+d] = (2.0 * scale - 1.0) * vel[i_atom_3+d]
	  + scale * time_step * vel_diff;
	real tmp_vel = (vel_next[i_atom_3+d] + vel[i_atom_3+d]) * 0.5;
	vel_norm += tmp_vel*tmp_vel;
      }
      kine_pre += vel_norm * mass[i_atom];
    }
    commotion->cancel_translation(atomids_rev, vel_next);    
    for(int i_atom = 0, i_atom_3 = 0;
	i_atom < n_atoms; i_atom++, i_atom_3+=3){
      for(int d=0; d < 3; d++){
	crd[i_atom_3+d] = crd_prev[i_atom_3+d] + time_step * vel_next[i_atom_3+d];
      }
    }
  }  
  if(!converge){
    cout << "Thermostat was not converged." << endl;
  }
  return 0;
}
////////////////////////////////////////////////////////////


ThermostatHooverEvans::ThermostatHooverEvans() 
  : ThermostatObject(){
}
ThermostatHooverEvans::~ThermostatHooverEvans(){
}

int ThermostatHooverEvans::set_constant(int n_atoms, real_pw* mass_inv, real* vel, real* force){
  real k0=0.0;
  for(int i=0, i_3=0;
      i < n_atoms; i++){
    real pre=0.0;
    for (int d=0; d < 3; d++, i_3++){
      pre += vel[i_3] * vel[i_3];
    }
    k0 += pre * mass_inv[i];
  }
  const_k0_inv = 1.0/k0;
  cout << "Gaussian constraint constant K0 : " << k0 << endl;
  return 0;
}

int ThermostatHooverEvans::apply_thermostat(int n_atoms,
					    real_fc* work,
					    real* vel, real* vel_next,
					    real_pw* mass,
					    real_pw* mass_inv){
  real vf = 0.0;
  real ff = 0.0;
  for (int i = 0, i_3 = 0;
       i < n_atoms; i++){
    real vf_i = 0.0;
    real ff_i = 0.0;
    for(int d = 0; d < 3; d++, i_3++){
      vf_i += vel[i_3]*work[i_3];
      ff_i += work[i_3]*work[i_3];
    }
    vf += vf_i * mass_inv[i];
    ff += ff_i * mass_inv[i];
  }
  real xi = -const_k0_inv * vf;
  real alpha = sqrt(const_k0_inv * ff);
  real beta = exp(-alpha * time_step);
  real gamma = (xi - alpha)/(xi + alpha);
  real gamma_beta = gamma/beta;
  for (int i = 0, i_3 = 0;
       i < n_atoms*3; i++){
    vel_next[i_3] = (1.0 - gamma)/(beta - gamma_beta) * 
      (vel[i_3] + work[i_3] * 
       (1.0 + gamma - beta - gamma_beta) / (alpha - gamma * alpha)
       );
  }

  return 0;
}
int ThermostatHooverEvans::apply_thermostat_with_shake(int n_atoms,
						   real_fc* work,
						   real* crd, real* crd_prev,
						   real* vel, real* vel_next,
						   real_pw* mass,
						   real_pw* mass_inv,
						   ConstraintObject* constraint,
						   PBC* pbc,
						   real* buf_crd,
						   const int max_loops,
						       const real tolerance,
						       COMMotion* commotion,
						       int* atomids_rev){
  cout << "Warning: This function does nothing. ThermostatHooverEvans::apply_thermostat_with_shake " << endl;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////


ThermostatNoseHoover::ThermostatNoseHoover() 
  : ThermostatObject(){
}
ThermostatNoseHoover::~ThermostatNoseHoover(){
}

int ThermostatNoseHoover::set_constant(int n_atoms, real* mass_inv, real* vel, real* force){
  return 0;
}
int ThermostatNoseHoover::apply_thermostat(int n_atoms,
					    real_fc* work,
					    real* vel, real* vel_next,
					    real_pw* mass,
					    real_pw* mass_inv){
  return 0;
}
int ThermostatNoseHoover::apply_thermostat_with_shake(int n_atoms,
						   real_fc* work,
						   real* crd, real* crd_prev,
						   real* vel, real* vel_next,
						   real_pw* mass,
						   real_pw* mass_inv,
						   ConstraintObject* constraint,
						   PBC* pbc,
						   real* buf_crd,
						   const int max_loops,
						      const real tolerance,
						      COMMotion* commotion,
						      int* atomids_rev){
  cout << "Warning: This function does nothing. ThermostatNoseHoover::apply_thermostat_with_shake " << endl;
  return 0;
}
