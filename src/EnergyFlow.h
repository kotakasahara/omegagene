#ifndef __ENERGY_FRLOW_H__
#define __ENERGY_FRLOW_H__

#include "CelesteObject.h"
#include <vector>
using namespase std;

class EnergyFlow : public CelesteObject{
 private:
  static int n_pair_terms;
  static int n_trio_terms;
  static int n_quad_terms;
 public:
  int n_atoms;
  vector< pair<int,int> > pair_terms_atoms;
  vector<real> pair_terms_ene;
  vector<real> pair_terms_work;
  vector<real> pair_ene_i;
  vector<real> pair_eneflo_kine;
  vector<real> pair_eneflo_pote;

  vector< trio<int,int,int> >  trio_terms_atoms;
  vector<real> trio_terms_ene;
  vector<real> trio_terms_work;
  vector<real> trio_ene_i;
  vector<real> trio_eneflo_kine;
  vector<real> trio_eneflo_pote;

  vector<real> quad_terms_atoms;
  vector<real> quad_terms_ene;
  vector<real> quad_terms_work;
  vector<real> quad_ene_i;
  vector<real> quad_eneflo_kine;
  vector<real> quad_eneflo_pote;

  EnergyFlow();
  int initial(int in_n_atoms);
  int reserve_values();
  int reset_values();

};

#endif
