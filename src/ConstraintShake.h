#ifndef __CONSTRAINT_SHAKE_H__
#define __CONSTRAINT_SHAKE_H__

#include "Constraint.h"

class ConstraintShake : public Constraint {
 private:
 protected:
 public:
  ConstraintShake();
  ~ConstraintShake();
  int apply_constraint(real* in_crd, real* in_crd_prev, real* mass,
		       PBC* pbc);
  int shake_pair(real* in_crd, real* in_crd_prev, real* mass,
		 PBC* pbc);
  int shake_trio(real* in_crd, real* in_crd_prev, real* mass,
		 PBC* pbc);
  int shake_quad(real* in_crd, real* in_crd_prev, real* mass,
		 PBC* pbc);


};

#endif
