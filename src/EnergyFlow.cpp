#include "EnergyFlow.h"

EnergyFlow::EnergyFlow(){
}

int EnergyFlow::initial(int in_n_atoms){
  n_atoms = in_n_atoms;
  reserve_values();
}

int EnergyFlow::reserve_values(){
  pair_ene_i.reserve(n_atoms);
  pair_eneflo_kine.reserve(n_atoms);
  pair_eneflo_pote.reserve(n_atoms);
  trio_ene_i.reserve(n_atoms);
  trio_eneflo_kine.reserve(n_atoms);
  trio_eneflo_pote.reserve(n_atoms);
  quad_ene_i.reserve(n_atoms);
  quad_eneflo_kine.reserve(n_atoms);
  quad_eneflo_pote.reserve(n_atoms);

}

int EnergyFlow::reset_values(){
}

int EnergyFlow::append_pair_terms(atoms, ene, work){
  pair_terms_atoms.push_back(atoms);
  pair_term_ene.push_back(ene);
  pair_terms_work.push_back(work);
}
