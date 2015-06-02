#ifndef __THERMOSTAT_H__
#define __THERMOSTAT_H__

#include "CelesteObject.h"
#include <cmath>
#include "Constraint.h"
using namespace std;

class ThermostatObject : public CelesteObject{
 private:
 protected:
  int d_free;
  real temperature;
  real temperature_coeff;
  real time_step;
  real time_step_inv_sq;
  real const_k0_inv;
 public:
  ThermostatObject();
  ~ThermostatObject();
  int set_time_step(real in_time_step);
  int set_temperature_coeff(int in_d_free);
  virtual int set_constant(int n_atoms, real_pw* mass, real* vel, real* force);
  virtual int apply_thermostat(const int n_atoms,
			       real_fc* work,
			       real* vel, real* vel_next,
			       real_pw* mass,
			       real_pw* mass_inv);


  virtual int apply_thermostat_with_shake(int n_atoms,
					  real_fc* work,
					  real* crd, real* crd_prev,
					  real* vel, real* vel_next,
					  real_pw* mass,
					  real_pw* mass_inv,
					  ConstraintObject* constraint,
					  PBC* pbc,
					  real* buf_crd,
					  const int max_loops,
					  const real tolerance);

  inline void set_temperature(real in_t){temperature=in_t;};
  inline real get_temperature(){return temperature;};
};

class ThermostatScaling : public ThermostatObject {
 private:
 protected:
 public:
  ThermostatScaling();
  ~ThermostatScaling();
  virtual int set_constant(int n_atoms, real_pw* mass, real* vel, real* force);
  virtual int apply_thermostat(const int n_atoms,
			       real_fc* work,
			       real* vel, real* vel_next,
			       real_pw* mass,
			       real_pw* mass_inv);
  virtual int apply_thermostat_with_shake(int n_atoms,
					  real_fc* work,
					  real* crd, real* crd_prev,
					  real* vel, real* vel_next,
					  real_pw* mass,
					  real_pw* mass_inv,
					  ConstraintObject* constraint,
					  PBC* pbc,
					  real* buf_crd,
					  const int max_loops,
					  const real tolerance);

};

class ThermostatHooverEvans : public ThermostatObject {
 private:
 protected:
 public:
  ThermostatHooverEvans();
  ~ThermostatHooverEvans();
  virtual int set_constant(int n_atoms, real_pw* mass, real* vel, real* force);
  virtual int apply_thermostat(const int n_atoms,
			       real_fc* work,
			       real* vel, real* vel_next,
			       real_pw* mass,
			       real_pw* mass_inv);
  virtual int apply_thermostat_with_shake(int n_atoms,
					  real_fc* work,
					  real* crd, real* crd_prev,
					  real* vel, real* vel_next,
					  real_pw* mass,
					  real_pw* mass_inv,
					  ConstraintObject* constraint,
					  PBC* pbc,
					  real* buf_crd,
					  const int max_loops,
					  const real tolerance);

};

class ThermostatNoseHoover : public ThermostatObject {
 private:
 protected:
 public:
  ThermostatNoseHoover();
  ~ThermostatNoseHoover();
  virtual int set_constant(int n_atoms, real* mass, real* vel, real* force);
  virtual int apply_thermostat(const int n_atoms,
			       real_fc* work,
			       real* vel, real* vel_next,
			       real_pw* mass,
			       real_pw* mass_inv);
  virtual int apply_thermostat_with_shake(int n_atoms,
					  real_fc* work,
					  real* crd, real* crd_prev,
					  real* vel, real* vel_next,
					  real_pw* mass,
					  real_pw* mass_inv,
					  ConstraintObject* constraint,
					  PBC* pbc,
					  real* buf_crd,
					  const int max_loops,
					  const real tolerance);

};

#endif