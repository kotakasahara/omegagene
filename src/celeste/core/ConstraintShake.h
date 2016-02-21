#ifndef __CONSTRAINT_SHAKE_H__
#define __CONSTRAINT_SHAKE_H__

#include "Constraint.h"

class ConstraintShake : public ConstraintObject {
 private:
 protected:
 public:
  ConstraintShake();
  ~ConstraintShake();
  virtual int apply_constraint(real* in_crd, real* in_crd_prev, real_pw* in_mass_inv,
		       PBC* pbc);
  int shake_pair(real* in_crd, real* in_crd_prev, real_pw* in_mass_inv,
		 PBC* pbc);
  int shake_trio(real* in_crd, real* in_crd_prev, real_pw* in_mass_inv,
		 PBC* pbc);
  int shake_quad(real* in_crd, real* in_crd_prev, real_pw* in_mass_inv,
		 PBC* pbc);


};

class ConstraintSettle: public ConstraintObject {
 private:
 protected:
 public:
  ConstraintSettle();
  ~ConstraintSettle();
  virtual int apply_constraint(real* in_crd, real* in_crd_prev, real_pw* in_mass_inv,
			       PBC* pbc);
};

#endif
