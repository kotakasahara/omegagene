#include "EnergyCalc.h"
#include <bitset>

EnergyCalc::EnergyCalc(MmSystem* in_mmsys,
		       SubBox* in_subbox)
  : EnergyCalcObject(in_mmsys, in_subbox){
}

int EnergyCalc::set_config_parameters(const Config* in_cfg){
  if(DBG>=1)
    cout << "DBG1: EnergyCalc::set_config_parameters()"<<endl;
  mmsys->pbc.print_pbc();

  ff = new ForceField();
  ff->set_config_parameters(in_cfg);
  cutoff = in_cfg->cutoff;
  
  return 0;
}
int EnergyCalc::set_dev_pointers(real_pw *in_d_crd,
				 int *in_d_grid_atom_index,
				 int *in_d_grid_n_atoms){
  return 0;
}
int EnergyCalc::initial_preprocess(){
  if (DBG >= 1)
    cout << "DBG1: EnergyCalc::initial_preprocess"<<endl;
  ff->initial_preprocess(&mmsys->pbc);
  ff->cal_self_energy((const int&)mmsys->n_atoms,
		      (const int&)mmsys->n_bonds,
		      (const int**&)mmsys->bond_atomid_pairs,
		      (const int&)mmsys->n_angles,
		      (const int**&)mmsys->angle_atomid_triads,
		      (const int&)mmsys->n_torsions,
		      (const int**&)mmsys->torsion_atomid_quads,
		      (const int*&)mmsys->torsion_nb14,
		      (const real_pw*&)mmsys->charge,
		      mmsys->energy_self,
		      mmsys->energy_self_sum);
  return 0;
}

int EnergyCalc::calc_energy(){
  //  if (DBG>=1){
  cout << "EnergyCalc::calc_energy()" <<endl;
  //  }
  //calc_energy_pairwise_nogrid();
  calc_energy_pairwise();
  calc_energy_bonds();
  calc_energy_angles();
  calc_energy_torsions();
  calc_energy_impros();
  calc_energy_14nb();
  calc_energy_ele_excess();
  //mmsys->nsgrid.get_ene_forces(mmsys->pote_vdw, mmsys->pote_ele, mmsys->force);
  //mmsys->subbox.set_pair_forces_from_nsgrid();
  mmsys->pote_ele += mmsys->energy_self_sum;
  //subbox->add_work_from_minicell(mmsys->nsgrid.get_work(),
  //mmsys->nsgrid.get_atomids_rev());
  subbox->add_work_from_minicell();
  mmsys->pote_vdw += mmsys->nsgrid.get_energy()[0];
  mmsys->pote_ele += mmsys->nsgrid.get_energy()[1];
  return 0;
}

int EnergyCalc::calc_energy_bonds(){
  for (int i=0; i < mmsys->n_bonds; i++){
    real_pw ene;
    real_pw work[3];
    int atomid1 = mmsys->bond_atomid_pairs[i][0];
    int atomid2 = mmsys->bond_atomid_pairs[i][1];
    //cout << "EnergyCalc::calc_energy_bonds() " << atomid1 << " - " << atomid2 << endl;
    ff->calc_bond(ene, work,
		  mmsys->crd[atomid1],
		  mmsys->crd[atomid2],
		  mmsys->bond_epsiron[i],
		  mmsys->bond_r0[i]);
    if(ene >= 1e3){
      cout << "HIGH BOND ENERGY: atoms " << atomid1 << "-" << atomid2 
	   << "  " << ene << "kcal/mol." << endl;
      cout << "  atom1: " << mmsys->crd[atomid1][0] <<", "
	   << mmsys->crd[atomid1][1] << ", "
	   << mmsys->crd[atomid1][2] <<endl;
      cout << "  atom2: " << mmsys->crd[atomid2][0] <<", "
	   << mmsys->crd[atomid2][1] << ", "
	   << mmsys->crd[atomid2][2] <<endl;
    }

    mmsys->pote_bond += ene;
    //    cout << setprecision(15);
    //    cout << "DBG0203bw ";
    for(int d=0; d<3; d++){
      mmsys->force[atomid1][d] += work[d];
      mmsys->force[atomid2][d] -= work[d];
    }
    /*
      cout << "frc1 " << atomid1 << " " << mmsys->force[atomid1][0] << " , ";
      cout << mmsys->force[atomid1][1] << " , ";
    cout << mmsys->force[atomid1][2] << endl;
    cout << "frc2 " << atomid2 << " " << mmsys->force[atomid2][0] << " , ";
    cout << mmsys->force[atomid2][1] << " , ";
    cout << mmsys->force[atomid2][2] << endl;
    */
    
    //    if(DBG >= 1)
    //cout << "D! bond " << atomid1 << " " << atomid2 <<" "<<ene <<" "  <<work[0] << " "  << work[1] << " " << work[2] << endl;
  }
  return 0;
}

int EnergyCalc::calc_energy_angles(){
  for (int i=0; i < mmsys->n_angles; i++){
    real_pw ene;
    real_pw work1[3], work2[3];
    int atomid1 = mmsys->angle_atomid_triads[i][0];
    int atomid2 = mmsys->angle_atomid_triads[i][1];
    int atomid3 = mmsys->angle_atomid_triads[i][2];
    ff->calc_angle(ene, work1, work2,
		   mmsys->crd[atomid1],
		   mmsys->crd[atomid2],
		   mmsys->crd[atomid3],
		   mmsys->angle_epsiron[i],
		   mmsys->angle_theta0[i]);
    mmsys->pote_angle += ene;
    for(int d=0; d<3; d++){
      mmsys->force[atomid1][d] += work1[d];
      mmsys->force[atomid2][d] -= work1[d] + work2[d];
      mmsys->force[atomid3][d] += work2[d];
    }
    /*
      cout << "D! angle ";
    cout << atomid1 << " " << atomid2 << " " << atomid3 <<" ";
    cout << ene << " ";
    for(int d=0; d<3; d++) cout << work1[d] << " ";
    for(int d=0; d<3; d++) cout << work2[d] << " ";
    cout << endl;
    */
  }
  return 0;
}
int EnergyCalc::calc_energy_torsions(){
  for (int i=0; i < mmsys->n_torsions; i++){
    real_pw ene;
    real_pw work1[3], work2[3], work3[3];
    int atomid1 = mmsys->torsion_atomid_quads[i][0];
    int atomid2 = mmsys->torsion_atomid_quads[i][1];
    int atomid3 = mmsys->torsion_atomid_quads[i][2];
    int atomid4 = mmsys->torsion_atomid_quads[i][3];

    ff->calc_torsion(ene, work1, work2, work3,
		     mmsys->crd[atomid1],
		     mmsys->crd[atomid2],
		     mmsys->crd[atomid3],
		     mmsys->crd[atomid4],
		     mmsys->torsion_energy[i],
		     mmsys->torsion_overlaps[i],
		     mmsys->torsion_symmetry[i],
		     mmsys->torsion_phase[i]);

    mmsys->pote_torsion += ene;
    for(int d=0; d<3; d++){
      mmsys->force[atomid1][d] += work1[d];
      mmsys->force[atomid2][d] += work2[d];
      mmsys->force[atomid3][d] -= work1[d] + work2[d] + work3[d];
      mmsys->force[atomid4][d] += work3[d];      
    }
    /*
    cout << "D! torsion " << atomid1 << " ";
    cout << atomid2 << " " << atomid3 << " " << atomid4 << " ";
    cout  << ene << " ";
    for(int d=0; d<3; d++) cout << work1[d] << " ";
    for(int d=0; d<3; d++) cout << work2[d] << " ";
    for(int d=0; d<3; d++) cout << work3[d] << " ";
    cout << endl;
    */
  }
  return 0;
}
int EnergyCalc::calc_energy_impros(){
  for (int i=0; i < mmsys->n_impros; i++){
    real_pw ene;
    real_pw work1[3], work2[3], work3[3];
    int atomid1 = mmsys->impro_atomid_quads[i][0];
    int atomid2 = mmsys->impro_atomid_quads[i][1];
    int atomid3 = mmsys->impro_atomid_quads[i][2];
    int atomid4 = mmsys->impro_atomid_quads[i][3];

    ff->calc_torsion(ene, work1, work2, work3,
		     mmsys->crd[atomid1],
		     mmsys->crd[atomid2],
		     mmsys->crd[atomid3],
		     mmsys->crd[atomid4],
		     mmsys->impro_energy[i],
		     mmsys->impro_overlaps[i],
		     mmsys->impro_symmetry[i],
		     mmsys->impro_phase[i]);

    mmsys->pote_impro += ene;
    for(int d=0; d<3; d++){
      mmsys->force[atomid1][d] += work1[d];
      mmsys->force[atomid2][d] += work2[d];
      mmsys->force[atomid3][d] -= work1[d] + work2[d] + work3[d];
      mmsys->force[atomid4][d] += work3[d];      
    }

    /*
      cout << "D! impro " << atomid1 << " ";
      cout << atomid2 << " " << atomid3 << " " << atomid4 << " ";
      cout <<ene << " " ;
      for(int d=0; d<3; d++) cout << work1[d] << " ";
      for(int d=0; d<3; d++) cout << work2[d] << " ";
      for(int d=0; d<3; d++) cout << work3[d] << " ";
      cout << endl;
    */
  }
return 0;
}
int EnergyCalc::calc_energy_14nb(){
  for (int i=0; i < mmsys->n_nb14; i++){
    real_pw ene_vdw;
    real_pw ene_ele;
    real_fc work[3];
    int atomid1 = mmsys->nb14_atomid_pairs[i][0];
    int atomid2 = mmsys->nb14_atomid_pairs[i][1];
    int atomtype1 = mmsys->nb14_atomtype_pairs[i][0];
    int atomtype2 = mmsys->nb14_atomtype_pairs[i][1];

    ff->calc_14pair(ene_vdw, ene_ele, work,
		    mmsys->crd[atomid1],
		    mmsys->crd[atomid2],
		    mmsys->lj_6term[atomtype1 * mmsys->n_lj_types + atomtype2],
		    mmsys->lj_12term[atomtype1 * mmsys->n_lj_types + atomtype2],
		    mmsys->charge[atomid1],
		    mmsys->charge[atomid2],
		    mmsys->nb14_coeff_vdw[i],
		    mmsys->nb14_coeff_ele[i]);
    //cout <<  "nb14:"<<atomid1 << "-" << atomid2 <<" " << ene_vdw << " " << ene_ele << endl;;
    mmsys->pote_14vdw += ene_vdw;
    mmsys->pote_14ele += ene_ele;
    for(int d=0; d<3; d++){
      mmsys->force[atomid1][d] += work[d];
      mmsys->force[atomid2][d] -= work[d];
    }
    


    //      cout << "D! 14interaction " << atomid1 << " " << atomid2 << " "
    //    << ene_vdw + ene_ele << " "
    //    << work[0] << " "  << work[1] << " " << work[2] << endl;

  }
  return 0;
}

int EnergyCalc::calc_energy_pairwise_nogrid(){

  cout << " E : " << mmsys->pote_vdw << ", " << mmsys->pote_ele << endl;
  double sum_dist = 0.0;
  double sum_dist_incut = 0.0;
  int atomid1sum = 0;
  int atomid2sum = 0;
  int atomid12mult = 0;
  int n_pairs = 0;
  int n_pairs_incutoff = 0;
  int n_pairs_15off = 0;
  double lj6mult = 0.0;
  double  lj12mult = 0.0;
  double chgmult = 0.0;
  double p_vdw = 0.0;
  double p_ele = 0.0;

  for(int atomid1=0; atomid1 < mmsys->n_atoms; atomid1++){
    int atomtype1 = mmsys->atom_type[atomid1];
    for(int atomid2=atomid1+1; atomid2 < mmsys->n_atoms; atomid2++){
      //pair<int,int> atompair(atomid1, atomid2);
      //set< pair<int,int> >::iterator itr_pair;
      //itr_pair = find(mmsys->nb15off_atomid_pairs.begin(),
      //mmsys->nb15off_atomid_pairs.end(),
      //atompair);
      //if(itr_pair != mmsys->nb15off_atomid_pairs.end())
      //continue;
      n_pairs++;

      if(mmsys->search_nb15off(atomid1, atomid2)){ n_pairs_15off++;  continue;	}

      //cout << "dbg0204 pair: " << atomid1 << "-" << atomid2 << endl;
      real ene_vdw;
      real ene_ele;
      real_fc work[3] = {0.0, 0.0, 0.0};
      int atomtype2 = mmsys->atom_type[atomid2];
      real param_6term = mmsys->lj_6term[atomtype1 * mmsys->n_lj_types + atomtype2];
      real param_12term = mmsys->lj_12term[atomtype1 * mmsys->n_lj_types + atomtype2];      
      //      cout << "DBG1: pair " << atomid1 << " - " << atomid2 << endl;
      //for(int d=0; d<3;d++) cout<<mmsys->crd[atomid1][d]<<" ";
      //cout << endl;
      //for(int d=0; d<3;d++) cout<<mmsys->crd[atomid2][d]<<" ";
      //cout << endl;
      
      real_pw crd1[3] = {mmsys->crd[atomid1][0],
			 mmsys->crd[atomid1][1],
			 mmsys->crd[atomid1][2]};
      real_pw crd2[3] = {mmsys->crd[atomid2][0],
			 mmsys->crd[atomid2][1],
			 mmsys->crd[atomid2][2]};

      //test 
      double r12 = sqrt(pow(crd2[0]-crd1[0],2)+pow(crd2[1]-crd1[1],2)+pow(crd2[2]-crd1[2],2));
      sum_dist += r12;
      if(sum_dist > 100000) sum_dist -= 100000;
      if(ff->calc_pairwise(ene_vdw, ene_ele, work,
			   crd1, crd2,
			   param_6term, param_12term,
			   mmsys->charge[atomid1],
			   mmsys->charge[atomid2]))
	n_pairs_incutoff++;	
      if (ene_vdw != 0.0 || ene_ele != 0.0){
	sum_dist_incut += r12;
	if(sum_dist_incut > 100000) sum_dist_incut -= 100000;	
	lj6mult += param_6term;
	while(lj6mult > 100000) lj6mult -= 100000;

	lj12mult += param_12term;
	while(lj12mult > 100000) lj12mult -= 100000;

	chgmult += mmsys->charge[atomid1] * mmsys->charge[atomid2];
	if(chgmult > 100000) chgmult -= 100000;
	atomid1sum+=atomid1;
	atomid2sum+=atomid2;
	atomid12mult+=atomid2*atomid1;
	atomid1sum = atomid1sum%100000;
	atomid2sum = atomid2sum%100000;
	atomid12mult = atomid12mult%100000;

      }
      p_vdw += ene_vdw;
      p_ele += ene_ele;
      //mmsys->pote_vdw += ene_vdw;
      //mmsys->pote_ele += ene_ele;
      for (int d=0; d<3; d++){
	mmsys->force[atomid1][d] += work[d];
	mmsys->force[atomid2][d] -= work[d];
      }
      //if(DBG >= 1){
	      //cout << "ene : " << ene_vdw << " " << ene_ele << endl;
	//}
    }
  }
  mmsys->pote_vdw = p_vdw;
  mmsys->pote_ele = p_ele;
  cout << "n_pairs_15off " << n_pairs_15off << " / " << mmsys->n_nb15off << endl;
  cout << "15 pairs: " << n_pairs_incutoff << " / " << n_pairs << endl;
  cout << " E : " << mmsys->pote_vdw << ", " << mmsys->pote_ele << endl;
  cout << " E2 : " << p_vdw << ", " << p_ele << endl;
  cout << " atomidsum : " << atomid1sum << " " << atomid2sum << " " << atomid1sum + atomid2sum << " " << atomid12mult <<  endl;
  cout << " sum_dist: " <<  sum_dist << " - " << sum_dist_incut << endl;
  cout << " lj6: " << lj6mult << " lj12: "<< lj12mult <<endl;
  cout << " chg: " << chgmult << endl;
  
  return 0;
}

int EnergyCalc::calc_energy_pairwise(){

  double sum_dist = 0.0;
  double sum_dist_incut = 0.0;
  int atomid1sum = 0;
  int atomid2sum = 0;
  int atomid12mult = 0;
  double lj6mult = 0.0;
  double lj12mult = 0.0;
  double chgmult = 0.0;
  int n_pairs=0;
  int n_pairs_incutoff=0;
  int n_pairs_15off = 0;
  double p_vdw = 0.0;
  double p_ele = 0.0;
  mmsys->nsgrid.init_energy_work();
  //  cout << "calc_energy_pairwise_grid()" << mmsys->nsgrid.n_grid_pairs << endl;
  for(int cp=0; cp < mmsys->nsgrid.get_n_cell_pairs(); cp++){
    CellPair cellpair = mmsys->nsgrid.get_cell_pair(cp);
    int c1 = cellpair.cell_id1;

    //if(c1 != mmsys->nsgrid.get_n_cells()-1) continue;

    int c2 = cellpair.cell_id2;
    //int c2img[3];
    //c2img[0] = cellpair.cx;
    //c2img[1] = cellpair.cy;
    //c2img[2] = cellpair.cz;
    int n_atoms_c1 = mmsys->nsgrid.get_n_atoms_in_cell(c1);
    int n_atoms_c2 = mmsys->nsgrid.get_n_atoms_in_cell(c2);
    int atoms_index_c1 = mmsys->nsgrid.get_idx_cell_head_atom(c1);
    int atoms_index_c2 = mmsys->nsgrid.get_idx_cell_head_atom(c2);
    // grid_atoms1[i]  = atom_id
    //  i = 0, 1, ... n_atoms_c1 ;  atom id in each grid
    //int* grid_atoms1 = mmsys->nsgrid.grid_atom[c1];
    //int* grid_atoms2 = mmsys->nsgrid.grid_atom[c2];

    //cout << "GPTEST c1[" << c1 << "] " << n_atoms_c1 << " ; "
    //<< "c2[" << c2 << "] " << n_atoms_c2 << endl;

    //int tmp_start=-1;
    //int* a1_begin;
    //a1_begin = &tmp_start;
    int a2 = 0;
    //if(c1==c2) a1_begin = &a2;

    for (a2=0; a2 < N_ATOM_CELL; a2++){
      int atomid_grid2 = atoms_index_c2 + a2;
      //int atomid2 = mmsys->nsgrid.grid_atom[g2][a2];
      //if(mmsys->nsgrid.is_dummy(atomid_grid2)==1) continue;
      int atomid2 = mmsys->nsgrid.get_atomid_from_gridorder(atomid_grid2);
      //if(atomid2 == -1) continue;

      //mmsys->nsgrid.get_crd_from_gridorder(atomid_grid2, crd2);      
      //AtomInfo ai2 = 

      //for (int a1=(*a1_begin)+1; a1 < n_atoms_c1; a1++){
      //cout << "DBG1 " << " " << atomid2 << endl;
      if(atomid2 < 0) continue;
      //cout << "DBG2 " << " " << atomid2 << endl;
      real crd2[3];
      mmsys->nsgrid.get_crd(atomid_grid2, crd2[0], crd2[1], crd2[2]);

      //real crd2[3] = {mmsys->crd[atomid2][0],
      //mmsys->crd[atomid2][1],
      //mmsys->crd[atomid2][2]};
      
      //real tc1 = crd2[0] + mmsys->pbc.L[0] * (real)c2img[0];
      //real tc2 = crd2[1] + mmsys->pbc.L[1] * (real)c2img[1];
      //real tc3 = crd2[2] + mmsys->pbc.L[2] * (real)c2img[2];
      fix_pbc_image(crd2, cellpair.image);
      //if(ret == 0){
      //cout << crd2[0] << " " << crd2[1] << " " << crd2[2] << endl;
      //cout << tc1 << " " << tc2 << " " << tc3 << endl;
      //exit(1);
      //}
      
      for (int a1=0; a1 < N_ATOM_CELL; a1++){
	//cout << "DBG3 " << " " << atomid2 << " "  << endl;
	int atomid_grid1 = atoms_index_c1 + a1;
	//if(mmsys->nsgrid.is_dummy(atomid_grid1)==1) continue;
	//int atomtype1 = mmsys->atom_type[atomid1];
	//int atomid1 = mmsys->nsgrid.grid_atom[c1][a1];
	int atomid1 = mmsys->nsgrid.get_atomid_from_gridorder(atomid_grid1);
	//if(atomid1 == -1) continue;
	//AtomInfo ai1 = mmsys->nsgrid.atominfo[atomid_grid1];

	//mmsys->nsgrid.get_crd_from_gridorder(atomid_grid1, crd1);

	//cout << "pair " << atomid1 << "-" << atomid2 << " (grid:" << atomid_grid1 << "-" << atomid_grid2<<")" << endl;
	
	//if(atomid1 >= atomid2 && c1 == c2 ){continue;}
	n_pairs ++;
	//if(mmsys->search_nb15off(atomid1, atomid2)){  n_pairs_15off++; continue; }

	if (check_nb15off(a1, a2, cellpair.pair_mask) ){ 
	  //cout << "nb15:" << a1 << "-" << a2 << "(" << atomid1 <<"-" << atomid2 
	    //<< ") aidingrid(" << atomid_grid1 <<"-" << atomid_grid2
	  //<<") cell:" << c1 << "-" << c2<< " " 
	  //<<":" << bitset<32>(cellpair.pair_mask[0]) << " "
	  //<< bitset<32>(cellpair.pair_mask[1]) << endl;
	  n_pairs_15off++; continue; }
	real crd1[3];
	mmsys->nsgrid.get_crd(atomid_grid1, crd1[0], crd1[1], crd1[2]);

	real ene_vdw = 0.0;
	real ene_ele = 0.0;
	real_fc work[3] = {0.0, 0.0, 0.0};
	real param_6term  = mmsys->lj_6term[mmsys->atom_type[atomid1]  * mmsys->n_lj_types + mmsys->atom_type[atomid2]];
	real param_12term = mmsys->lj_12term[mmsys->atom_type[atomid1] * mmsys->n_lj_types + mmsys->atom_type[atomid2]];
	double r12 = sqrt(pow(crd2[0]-crd1[0],2)+pow(crd2[1]-crd1[1],2)+pow(crd2[2]-crd1[2],2));
	sum_dist += r12;
	if(sum_dist > 100000) sum_dist -= 100000;
	
	if(ff->calc_pairwise(ene_vdw, ene_ele, work,
			     crd1,
			     crd2,
			     param_6term, param_12term,
			     mmsys->charge[atomid1],
			     mmsys->charge[atomid2])==0){
	  mmsys->nsgrid.add_energy(ene_vdw, ene_ele);
	  if(isnan(ene_vdw)){
	    cout << "dbg D! nonbond " << atomid1 << " "  << atomid2 <<" "
		 << " g:" << atomid_grid1 << " "  << atomid_grid2 <<" ";
	    cout << ene_vdw << " " <<  ene_ele << " " << r12 << " "
		 << crd1[0] << "," << crd1[1] << "," << crd1[2] <<" " 
		 << crd2[0] << "," << crd2[1] << "," << crd2[2] <<" " 
		 <<endl;
	  }
	  mmsys->nsgrid.add_work(atomid_grid1, work[0], work[1], work[2]);
	  mmsys->nsgrid.add_work(atomid_grid2, -work[0], -work[1], -work[2]);
	  n_pairs_incutoff++;
	  
	  //printf("DEBUG_W1 %d-%d (%15.12e, %15.12e, %15.12e)\n",
	  //atomid1, atomid2,
	  //work[0], work[1], work[2]);
	  //printf("DEBUG_E1 %d-%d %10e %10e %d %d %10e %10e\n", atomid1, atomid2,
	  //crd1[0], crd1[1], crd1[2],
	  //crd2[0], crd2[1], crd2[2],
	    //	       mmsys->charge[atomid1],
	    //	       mmsys->charge[atomid2],
	    //mmsys->atom_type[atomid1], mmsys->atom_type[atomid2],
	    //ene_vdw, ene_ele);
	}
	p_vdw += ene_vdw;
	p_ele += ene_ele;
	if (ene_vdw != 0.0 || ene_ele != 0.0){
	sum_dist_incut += r12;
	if(sum_dist_incut > 100000) sum_dist_incut -= 100000;	
	  if(atomid1 > atomid2){
	    atomid1sum+=atomid2;
	    atomid2sum+=atomid1;
	  }else{
	    atomid1sum+=atomid1;
	    atomid2sum+=atomid2;
	  }
	lj6mult += param_6term;
	while(lj6mult > 100000) lj6mult -= 100000;
	lj12mult += param_12term;
	while(lj12mult > 100000) lj12mult -= 100000;
	chgmult += mmsys->charge[atomid1] * mmsys->charge[atomid2];
	if(chgmult > 100000) chgmult -= 100000;
	  atomid12mult+=atomid1*atomid2;
	  atomid1sum = atomid1sum%100000;
	  atomid2sum = atomid2sum%100000;
	  atomid12mult = atomid12mult%100000;
	}
	//mmsys->pote_vdw += ene_vdw;
	//mmsys->pote_ele += ene_ele;
	//for (int d=0; d<3; d++){
	//mmsys->force[atomid1][d] += work[d];
	//mmsys->force[atomid2][d] -= work[d];
	//cout << work[d] << " ";
	//}
	//cout <<  endl;
      }
    }
  }
  //mmsys->pote_vdw += p_vdw;
  //mmsys->pote_ele += p_ele;
  //

  cout << "nb15off pairs " << n_pairs_15off << " / "  <<mmsys->n_nb15off << endl;
  cout << "15 pairs: " << n_pairs_incutoff << " / " << n_pairs << endl;
  cout << " E : " << mmsys->pote_vdw << ", " << mmsys->pote_ele << endl;
  cout << " E2 : " << p_vdw << ", " << p_ele << endl;
  cout << " atomidsum : " << atomid1sum << " " << atomid2sum << " " << atomid1sum + atomid2sum << " " << atomid12mult << endl;
  cout << " sum_dist: " <<  sum_dist << " - " << sum_dist_incut << endl;
  cout << " lj6: " << lj6mult << " lj12: "<<lj12mult <<endl;
  cout << " chg: " << chgmult << endl;

  return 0;
}

int EnergyCalc::calc_energy_ele_excess(){
  real_fc ele = 0.0;
  for (int i=0; i < mmsys->n_excess; i++){
    real_pw ene;
    real_pw work[3];
    int atomid1 = mmsys->excess_pairs[i][0];
    int atomid2 = mmsys->excess_pairs[i][1];
    //cout << "excess " << atomid1 << " - " << atomid2 << endl;
    ff->calc_zms_excess(ene, work,
			mmsys->crd[atomid1],
			mmsys->crd[atomid2],
			mmsys->charge[atomid1],
			mmsys->charge[atomid2]);
    ele += ene;
    for(int d=0; d<3; d++){
      mmsys->force[atomid1][d] += work[d];
      mmsys->force[atomid2][d] -= work[d];      
    }
    //cout << "excess " << i << " " << mmsys->pote_ele << endl;;
  }
  //mmsys->pote_ele += ele;
  cout << "EnergyCalc::calc_energy_ele_excess ene:" << ele <<endl;
  return 0;
}

bool EnergyCalc::check_nb15off(const int& a1, const int& a2,
			       const int* bitmask){
  int bit_pos      = a2 * N_ATOM_CELL + a1;
  int mask_id      = bit_pos / 32;
  int interact_pos = bit_pos % 32;
  int interact = 1 << interact_pos;
  return (bitmask[mask_id] & interact) == interact;
}

int EnergyCalc::fix_pbc_image(real* crd, const int image){
  if(image == 0) return 1;
  
  //cout << "fix_... " << crd[0] << " " << crd[1] << " " << crd[2] << endl;
  if     ((image & 1)        == 1)
    crd[0] -= mmsys->pbc.L[0];
  else if((image & 2) == 2) crd[0] += mmsys->pbc.L[0];
  if     ((image & 4) == 4) crd[1] -= mmsys->pbc.L[1];
  else if((image & 8) == 8) crd[1] += mmsys->pbc.L[1];
  if     ((image & 16) == 16) crd[2] -= mmsys->pbc.L[2];
  else if((image & 32) == 32) crd[2] += mmsys->pbc.L[2];
  //cout << "fix_pbc_image " << image << " " << (image & (1<<4)) << " " << (1<<4) << " "
  //<< mmsys->pbc.L[2] << " "
  //<< crd[0] << " " << crd[1] << " " << crd[2] << endl;
  //exit(1);
  return 0;
}
