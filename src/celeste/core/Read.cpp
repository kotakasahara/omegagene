#include "Read.h"
using namespace std;

const int Read::MAX_LEN_NAME = 256;

Read::Read(string inFn) : CelesteObject() {
    op       = false;
    filename = inFn;
}

int Read::open() {
    ifs.open(filename.c_str());
    if (!ifs) {
        cerr << "Cannot open " << filename << "." << endl;
        return 1;
    }
    op = true;
    set_conv_endian_false();
    return 0;
}
int Read::close() {
    ifs.close();
    op = false;
    return 0;
}

vector<int> Read::load_integers() {
    vector<int> intvec;
    open();
    string buf;
    while (getline(ifs, buf)) { intvec.push_back(atoi(buf.c_str())); }
    close();
    return intvec;
}

vector<string> Read::load_strings() {
    vector<string> strvec;
    open();
    string buf;
    while (getline(ifs, buf)) {
        stringstream ss(buf);
        string       str;
        ss >> str;
        strvec.push_back(str);
    }
    close();
    return strvec;
}

int Read::load_launch_set(MmSystem &mmsys) {
    if (open() == 1) { exit(1); }
    cout << "-- Load Celeste input package." << endl;
    load_ls_header(mmsys);
    if (size_box > 0) {
        cout << "--- Load PBC definition   : " << size_box << " bytes." << endl;
        load_ls_box(mmsys);
    }
    if (size_crd > 0) {
        cout << "--- Load atom coordinates : " << size_crd << " bytes." << endl;
        load_ls_crd(mmsys);
    }
    if (size_vel > 0) {
        cout << "--- Load atom velocities  : " << size_vel << " bytes." << endl;
        load_ls_vel(mmsys);
    }
    if (size_topol > 0) {
        cout << "--- Load topology data    : " << size_topol << " bytes." << endl;
        load_ls_tpl(mmsys);
    }
    if (size_constraint > 0) {
        cout << "--- Load constraint definition : " << size_constraint << " bytes." << endl;
        load_ls_constraint(&mmsys.constraint);
    }
    if (size_settle > 0) {
        cout << "--- Load SETTLE definition : " << size_settle << " bytes." << endl;
        load_ls_constraint(&mmsys.settle);
    }
    if (size_extended > 0) {
        cout << "--- Load extendeded ensemble definition : " << size_extended << " bytes." << endl;
        load_ls_vmcmd(mmsys);
    }
    if (size_groups > 0) {
        cout << "--- Load atom group definition : " << size_groups << " bytes." << endl;
        load_ls_atom_groups(mmsys);
    }
    if (size_dist_restraint > 0) {
        cout << "--- Load distance restraint definition : " << size_dist_restraint << " bytes." << endl;
        load_ls_dist_restraint(mmsys.dist_restraint);
    }
    if (size_pos_restraint > 0) {
        cout << "--- Load position restraint definition : " << size_pos_restraint << " bytes." << endl;
        mmsys.n_pos_restraints = load_ls_pos_restraint(mmsys.pos_restraint);
	mmsys.d_free += 3;
	
    }
    if (size_extended_vcmd > 0) {
        cout << "--- Load VcMD: " << size_extended_vcmd << " bytes." << endl;
        load_ls_vcmd(mmsys);
    }
    if (size_group_coord > 0) {
        cout << "--- Load group coordinates for restarting V-AUS: " << size_group_coord << " bytes." << endl;
        load_ls_group_coord(mmsys);
    }
    
    if (size_pote > 0) {
        cout << "--- Load potential: " << size_pote << " bytes." << endl;
        load_ls_pote(mmsys);
    }
    
    // cout << "load_ls_pcluster()" << endl;
    // load_ls_pcluster(mmsys);
    close();
    return 0;
}

int Read::load_ls_header(MmSystem &mmsys) {
    cout << "--- Load file header." << endl;

    int magic;
    ifs.read((char *)&magic, sizeof(int));
    if (magic != MAGIC_NUMBER) {
        set_conv_endian_true();
        cerr << magic << endl;
        magic = reverse_endian(magic);
        if (magic != MAGIC_NUMBER) {
            stringstream ss;
            ss << "ERROR: " << filename << " : the first 4 bytes were not [" << MAGIC_NUMBER << "] but [" << magic
               << "]" << endl;
            error_exit(ss.str(), "1A00004");
            exit(1);
        }
    }

    int  buf_int;
    char version_c[MAX_LEN_NAME];
    read_bin_values(&buf_int, 1); // length of string
    ifs.read(version_c, buf_int);
    mmsys.launchset_version = string(version_c);
    cout << "---- Input file format : version " << mmsys.launchset_version << endl;
    if (mmsys.launchset_version != LS_VERSION) {
        stringstream ss;
        ss << "ERROR: " << filename << " : the version is incompatible." << endl;
        ss << "[" << LS_VERSION << "] is required." << endl;
        ss << "The input file is [" << mmsys.launchset_version << "]" << endl;
        error_exit(ss.str(), "1A00003");
        exit(1);
    }

    read_bin_values(&size_box, 1);
    read_bin_values(&size_crd, 1);
    read_bin_values(&size_vel, 1);
    read_bin_values(&size_topol, 1);
    read_bin_values(&size_constraint, 1);
    read_bin_values(&size_settle, 1);
    read_bin_values(&size_extended, 1);
    read_bin_values(&size_groups, 1);
    read_bin_values(&size_dist_restraint, 1);
    read_bin_values(&size_pos_restraint, 1);
    read_bin_values(&size_extended_vcmd, 1);
    read_bin_values(&size_group_coord, 1);
    read_bin_values(&size_pote, 1);

    if (DBG == 1) {
        cout << "size_box:            " << size_box << endl;
        cout << "size_crd:            " << size_crd << endl;
        cout << "size_vel:            " << size_vel << endl;
        cout << "size_topol:          " << size_topol << endl;
        cout << "size_constraint:     " << size_constraint << endl;
        cout << "size_settle:         " << size_settle << endl;
        cout << "size_extended:       " << size_extended << endl;
        cout << "size_groups:         " << size_groups << endl;
        cout << "size_dist_restraint: " << size_dist_restraint << endl;
        cout << "size_pos_restraint:  " << size_pos_restraint << endl;
        cout << "size_group_coord:    " << size_group_coord << endl;
        cout << "size_extended_vcmd:  " << size_extended_vcmd << endl;
        cout << "size_pote:           " << size_pote << endl;
    }

    return 0;
}

int Read::load_ls_box(MmSystem &mmsys) {
  // read_bin_values(&size_box, 1);
  double pbc_val[12];
  cout << "load_ls_box : ";
  for (int i = 0; i < 12; i++) {
    read_bin_values(&pbc_val[i], 1);
    mmsys.pbc_val[i] = (real)pbc_val[i];
        cout << pbc_val[i] << " ";
  }
  mmsys.pbc.set_pbc(mmsys.pbc_val);
  return 0;
}

int Read::load_ls_crd(MmSystem &mmsys) {
    // COORDINATES
    // int size_crd;
    //  read_bin_values(&size_crd, 1);
    // if(size_crd <= 0){
    // ERROR: size of coordinates data is zero
    // return 1;
    //}
    read_bin_values(&mmsys.n_atoms, 1);
    if (DBG == 1) { cout << "n_atoms: " << mmsys.n_atoms << endl; }
    if (mmsys.n_atoms <= 0) {
        // ERROR: the number of atoms is zero
        return 1;
    }
    mmsys.d_free = mmsys.n_atoms * 3 - 3;
    mmsys.alloc_atom_vars();

    for (int i = 0; i < mmsys.n_atoms; i++) {
        double x, y, z;
        read_bin_values(&x, 1);
        read_bin_values(&y, 1);
        read_bin_values(&z, 1);
        mmsys.crd[i][0] = (real)x;
        mmsys.crd[i][1] = (real)y;
        mmsys.crd[i][2] = (real)z;

        // cout << "load_ls_crd " << i << " : ";
        // cout << mmsys.crd[i][0] << " ";
        // cout << mmsys.crd[i][1] << " ";
        // cout << mmsys.crd[i][2] << endl;
    }
    return 0;
}
int Read::load_ls_vel(MmSystem &mmsys) {
    // VELOCITIES
    // int size_crd;
    // read_bin_values(&size_crd, 1);
    // if(size_crd <= 0){
    // ERROR: size of coordinates data is zero
    // return 1;
    //}
    int n_atoms;
    read_bin_values(&n_atoms, 1);

    if (DBG == 1) { cout << "read_ls_vel ... n_atoms: " << n_atoms << endl; }
    if (mmsys.n_atoms != n_atoms) {
        cerr << "ERROR" << endl;
        // ERROR: the number of atoms is inconsistent to the coordinates field
        return 2;
    }

    for (int i = 0; i < mmsys.n_atoms; i++) {
        double x, y, z;
        read_bin_values(&x, 1);
        read_bin_values(&y, 1);
        read_bin_values(&z, 1);
        mmsys.vel_just[i][0] = (real)x;
        mmsys.vel_just[i][1] = (real)y;
        mmsys.vel_just[i][2] = (real)z;
        /*
        cout << "load_ls_vel " << i << " : ";
        cout << mmsys.vel_next[i][0] << " ";
        cout << mmsys.vel_next[i][1] << " ";
        cout << mmsys.vel_next[i][2] << endl;
        */
    }
    return 0;
}

int Read::load_ls_tpl(MmSystem &mmsys) {
    // TOPOLOGIES
    // int size_tpl;
    // read_bin_values(&size_tpl, 1);
    // if(size_tpl <= 0){
    // ERROR: size of topology data is zero
    // return 1;
    //}
    int n_atoms;
    read_bin_values(&n_atoms, 1);
    if (mmsys.n_atoms != n_atoms) {
        // ERROR: the number of atoms is inconsistent to the coordinates field
        return 2;
    }
    // charge
    for (int i = 0; i < mmsys.n_atoms; i++) {
        double charge;
        read_bin_values(&charge, 1);
        mmsys.charge[i] = (real_pw)charge;
        // if (DBG==1){ cout << "charge " << i << " : " << mmsys.charge[i] << endl; }
    }
    // mass
    for (int i = 0; i < mmsys.n_atoms; i++) {
        double mass;
        read_bin_values(&mass, 1);
        mmsys.mass[i] = (real_pw)mass;
        //    if (DBG==1){ cout << "mass " << i << " : " << mmsys.mass[i] << endl; }
    }
    // atom_type
    for (int i = 0; i < mmsys.n_atoms; i++) {
        int atom_type;
        read_bin_values(&atom_type, 1);
        mmsys.atom_type[i] = atom_type - 1;
        //    if (DBG>=1){ cout << "atom_type " << i << " : " << mmsys.atom_type[i] << endl; }
    }

    // nbpair
    int size_lj;
    read_bin_values(&size_lj, 1);
    read_bin_values(&mmsys.n_lj_types, 1);
    read_bin_values(&mmsys.n_lj_type_pairs, 1);
    // if(DBG==1)
    // cout << "ljpair " << mmsys.n_lj_types << " " << mmsys.n_lj_type_pairs << endl;

    mmsys.alloc_lj_params();

    for (int i = 0; i < mmsys.n_lj_type_pairs; i++) {
        int    type1, type2;
        double lj6, lj12;
        read_bin_values(&type1, 1);
        read_bin_values(&type2, 1);
        read_bin_values(&lj6, 1);
        read_bin_values(&lj12, 1);
        // if(DBG==1)
        // cout << "nbpair " << type1 << " " << type2 << " " << lj6 << " " << lj12 << endl;
        mmsys.set_lj_pair_param(type1 - 1, type2 - 1, (real_pw)lj6, (real_pw)lj12);
    }

    //nbpair hps
    int size_lj_hps;
    read_bin_values(&size_lj_hps, 1);
    for (int i = 0; i < mmsys.n_lj_type_pairs; i++) {
        int    type1, type2;
        double cutoff, lambda;
        read_bin_values(&type1, 1);
        read_bin_values(&type2, 1);
        read_bin_values(&cutoff, 1);
        read_bin_values(&lambda, 1);
        mmsys.set_lj_pair_hps_param(type1 - 1, type2 - 1, (real_pw)cutoff, (real_pw)lambda);
	//cout << "readHPS " << type1 << "-" << type2 << " " << cutoff <<" " <<lambda<<endl; 
    }
    
    // bond
    int size_bond;
    read_bin_values(&size_bond, 1);
    read_bin_values(&mmsys.n_bonds, 1);
    if (DBG == 1) cout << "n_bonds " << mmsys.n_bonds << endl;

    mmsys.alloc_bonds();
    for (int i = 0; i < mmsys.n_bonds; i++) {
        int    atomid1, atomid2;
        double eps, r0;
        read_bin_values(&atomid1, 1);
        read_bin_values(&atomid2, 1);
        read_bin_values(&eps, 1);
        read_bin_values(&r0, 1);
        mmsys.set_bond_param(i, atomid1, atomid2, (real)eps, (real)r0);
        // if (DBG>=1){ cout << "bond " << i << " : " << atomid1 << " - " << atomid2 << endl;}
    }

    // angle
    int size_angle;
    read_bin_values(&size_angle, 1);
    read_bin_values(&mmsys.n_angles, 1);
    if (DBG == 1) cout << "n_angles " << mmsys.n_angles << endl;

    mmsys.alloc_angles();
    for (int i = 0; i < mmsys.n_angles; i++) {
        int    atomid1, atomid2, atomid3;
        double eps, theta0;
        read_bin_values(&atomid1, 1);
        read_bin_values(&atomid2, 1);
        read_bin_values(&atomid3, 1);
        read_bin_values(&eps, 1);
        read_bin_values(&theta0, 1);
        mmsys.set_angle_param(i, atomid1, atomid2, atomid3, (real)eps, (real)theta0);
    }

    // torsion
    int size_torsion;
    read_bin_values(&size_torsion, 1);
    read_bin_values(&mmsys.n_torsions, 1);
    if (DBG == 1) cout << "n_torsions " << mmsys.n_torsions << endl;

    mmsys.alloc_torsions();
    for (int i = 0; i < mmsys.n_torsions; i++) {
        int    atomid1, atomid2, atomid3, atomid4;
        double ene, phase;
        int    overlaps, symmetry, flag_14nb;
        read_bin_values(&atomid1, 1);
        read_bin_values(&atomid2, 1);
        read_bin_values(&atomid3, 1);
        read_bin_values(&atomid4, 1);
        read_bin_values(&ene, 1);
        read_bin_values(&overlaps, 1);
        read_bin_values(&symmetry, 1);
        read_bin_values(&phase, 1);
        read_bin_values(&flag_14nb, 1);
        mmsys.set_torsion_param(i, atomid1, atomid2, atomid3, atomid4, (real)ene, overlaps, symmetry, (real)phase,
                                flag_14nb);
        // cout << "torsion: " << atomid1 << "-" << atomid4 << " " << flag_14nb << endl;
    }

    // improper
    int size_impro;
    read_bin_values(&size_impro, 1);
    read_bin_values(&mmsys.n_impros, 1);
    if (DBG == 1) cout << "n_impros " << mmsys.n_impros << endl;

    mmsys.alloc_impros();
    for (int i = 0; i < mmsys.n_impros; i++) {
        int    atomid1, atomid2, atomid3, atomid4;
        double ene, phase;
        int    overlaps, symmetry, flag_14nb;
        read_bin_values(&atomid1, 1);
        read_bin_values(&atomid2, 1);
        read_bin_values(&atomid3, 1);
        read_bin_values(&atomid4, 1);
        read_bin_values(&ene, 1);
        read_bin_values(&overlaps, 1);
        read_bin_values(&symmetry, 1);
        read_bin_values(&phase, 1);
        read_bin_values(&flag_14nb, 1);
        mmsys.set_impro_param(i, atomid1, atomid2, atomid3, atomid4, (real)ene, overlaps, symmetry, (real)phase,
                              flag_14nb);
    }

    // 14 nonbond
    int size_nb14;
    read_bin_values(&size_nb14, 1);
    read_bin_values(&mmsys.n_nb14, 1);
    if (DBG == 1) cout << "n_nb14 " << mmsys.n_nb14 << endl;

    mmsys.alloc_nb14();
    for (int i = 0; i < mmsys.n_nb14; i++) {
        int    atomid1, atomid2;
        int    atomtype1, atomtype2;
        double coeff_vdw, coeff_ele;
        read_bin_values(&atomid1, 1);
        read_bin_values(&atomid2, 1);
        read_bin_values(&atomtype1, 1);
        read_bin_values(&atomtype2, 1);
        read_bin_values(&coeff_vdw, 1);
        read_bin_values(&coeff_ele, 1);
        mmsys.set_nb14_param(i, atomid1, atomid2, atomtype1 - 1, atomtype2 - 1, (real)coeff_ele, (real)coeff_vdw);
        // cout << "nb14:" << atomid1 <<"-"<< atomid2 << " " << coeff_vdw << " " << coeff_ele << endl;
    }

    // without 15
    int size_nb15off;
    int n_nb15off;
    read_bin_values(&size_nb15off, 1);
    read_bin_values(&n_nb15off, 1);
    if (DBG == 1) cout << "n_nb15off " << n_nb15off << endl;
    mmsys.alloc_nb15off();
    for (int i = 0; i < n_nb15off; i++) {
        int atomid1, atomid2;
        read_bin_values(&atomid1, 1);
        read_bin_values(&atomid2, 1);
        mmsys.set_nb15off(atomid1, atomid2);
        // mmsys.set_nb15off(atomid2, atomid1);
    }
    mmsys.alloc_excess_pairs();
    mmsys.set_excess_pairs();
    return 0;
}

int Read::load_ls_constraint(ConstraintObject *cst) {
    int n_const_2;
    int n_const_3;
    int n_const_4;

    read_bin_values(&n_const_2, 1);
    read_bin_values(&n_const_3, 1);
    read_bin_values(&n_const_4, 1);

    cst->set_max_n_constraints(n_const_2, n_const_3, n_const_4);
    cst->alloc_constraint();

    int    atomid1, atomid2, atomid3, atomid4;
    double dist1, dist2, dist3, dist4, dist5, dist6;
    // 2 atoms
    for (int i = 0; i < n_const_2; i++) {
        read_bin_values(&atomid1, 1);
        read_bin_values(&atomid2, 1);
        read_bin_values(&dist1, 1);
        dist1 = dist1 * dist1;
        cst->add_pair(atomid1, atomid2, (real_cst)dist1);
    }
    // 3 atoms
    for (int i = 0; i < n_const_3; i++) {
        read_bin_values(&atomid1, 1);
        read_bin_values(&atomid2, 1);
        read_bin_values(&atomid3, 1);
        read_bin_values(&dist1, 1);
        read_bin_values(&dist2, 1);
        read_bin_values(&dist3, 1);
        dist1 = dist1 * dist1;
        dist2 = dist2 * dist2;
        dist3 = dist3 * dist3;
        cst->add_trio(atomid1, atomid2, atomid3, (real_cst)dist1, (real_cst)dist2, (real_cst)dist3);
    }
    // 4 atoms
    for (int i = 0; i < n_const_4; i++) {
        read_bin_values(&atomid1, 1);
        read_bin_values(&atomid2, 1);
        read_bin_values(&atomid3, 1);
        read_bin_values(&atomid4, 1);
        read_bin_values(&dist1, 1);
        read_bin_values(&dist2, 1);
        read_bin_values(&dist3, 1);
        read_bin_values(&dist4, 1);
        read_bin_values(&dist5, 1);
        read_bin_values(&dist6, 1);
        dist1 = dist1 * dist1;
        dist2 = dist2 * dist2;
        dist3 = dist3 * dist3;
        dist4 = dist4 * dist4;
        dist5 = dist5 * dist5;
        dist6 = dist6 * dist6;
        cst->add_quad(atomid1, atomid2, atomid3, atomid4, (real_cst)dist1, (real_cst)dist2, (real_cst)dist3,
                      (real_cst)dist4, (real_cst)dist5, (real_cst)dist6);
    }

    return 0;
}

int Read::load_ls_vmcmd(MmSystem &mmsys) {
    int n_vs;
    read_bin_values(&n_vs, 1);
    int interval;
    read_bin_values(&interval, 1);
    double temperature;
    read_bin_values(&temperature, 1);

    mmsys.vmcmd->set_n_vstates(n_vs);
    mmsys.vmcmd->set_trans_interval(interval);
    mmsys.vmcmd->set_temperature((real)temperature);

    for (int i = 0; i < n_vs; i++) {
        int ord;
        read_bin_values(&ord, 1);
        mmsys.vmcmd->set_vs_order(i, ord);
        double lambda_low, lambda_high;
        double prob_low, prob_high;
        read_bin_values(&lambda_low, 1);
        read_bin_values(&lambda_high, 1);
        read_bin_values(&prob_low, 1);
        read_bin_values(&prob_high, 1);
        for (int j = 0; j < ord + 1; j++) {
            double buf;
            read_bin_values(&buf, 1);
            mmsys.vmcmd->set_vs_poly_param(i, j, (real)buf);
        }
        double alpha_low, alpha_high;
        read_bin_values(&alpha_low, 1);
        read_bin_values(&alpha_high, 1);
        mmsys.vmcmd->set_vs_params(i, (real)lambda_low, (real)lambda_high, (real)prob_low, (real)prob_high,
                                   (real)alpha_low, (real)alpha_high);
    }
    int init, seed;
    read_bin_values(&init, 1);
    read_bin_values(&seed, 1);
    mmsys.vmcmd->set_init_vs(init - 1);
    mmsys.vmcmd->set_random_seed(seed);

    return 0;
}

int Read::load_ls_atom_groups(MmSystem &mmsys) {
    int  n_groups;
    int *n_atoms_in_group;
    read_bin_values(&n_groups, 1);
    n_groups++;
    n_atoms_in_group    = new int[n_groups];
    n_atoms_in_group[0] = 0;
    cout << "n_groups " << n_groups << endl;
    mmsys.atom_group_names.push_back(string("null"));
    for (int i = 1; i < n_groups; i++) {
        int  len_name;
        char name[MAX_LEN_NAME];
        read_bin_values(&len_name, 1);
        ifs.read(name, len_name);
        read_bin_values(&n_atoms_in_group[i], 1);
        mmsys.atom_group_names.push_back(string(name));
        cout << "read atom groups : " << i << " " << n_atoms_in_group[i] << " " << name << " "
             << mmsys.atom_group_names[i] << endl;
    }
    mmsys.alloc_atom_groups(n_groups, n_atoms_in_group);
    int buf;
    for (int i = 0; i < n_groups; i++) {
        for (int j = 0; j < n_atoms_in_group[i]; j++) {
            read_bin_values(&buf, 1);
            mmsys.atom_groups[i][j] = buf - 1;
        }
    }
    delete[] n_atoms_in_group;
    return 0;
}
int Read::load_ls_dist_restraint(DistRestraintObject *dr) {
    int n_drunits;
    read_bin_values(&n_drunits, 1);
    dr->alloc_drunits(n_drunits);
    for (int i = 0; i < n_drunits; i++) {
        int   aid1, aid2;
        float coef_low, coef_high;
        float dist_low, dist_high;
        read_bin_values(&aid1, 1);
        read_bin_values(&aid2, 1);
        read_bin_values(&coef_low, 1);
        read_bin_values(&coef_high, 1);
        read_bin_values(&dist_low, 1);
        read_bin_values(&dist_high, 1);
        dr->add_drunit(aid1, aid2, coef_low, coef_high, dist_low, dist_high);
    }
    return 0;
}
int Read::load_ls_pos_restraint(PosRestraintObject *pr) {
  int n_prunits;
  read_bin_values(&n_prunits, 1);
  pr->alloc_prunits(n_prunits);
  for (int i = 0; i < n_prunits; i++) {
    int   aid;
    float crd_x, crd_y, crd_z;
    float dist_margin, coef;
    int rest_type;
    read_bin_values(&aid, 1);
    read_bin_values(&crd_x, 1);
    read_bin_values(&crd_y, 1);
    read_bin_values(&crd_z, 1);
    read_bin_values(&dist_margin, 1);
    read_bin_values(&coef, 1);
    read_bin_values(&rest_type, 1);
    cout << "dbg0708 read c " << coef << endl;
    int n_params; 
    float buf;
    real params[MAX_N_POSRES_PARAMS];
    read_bin_values(&n_params, 1);
    cout << "dbg0708 read n " << n_params << endl;
    for (int j = 0; j < n_params; j++){
      read_bin_values(&buf, 1);
      params[j] = buf;
      cout << "dbg0708 read " << buf << " " << params[j] << endl;
    }
    pr->add_prunit(aid, crd_x, crd_y, crd_z, dist_margin, coef, rest_type, n_params, params);
    
  }
  return n_prunits;
}

int Read::load_ls_group_coord(MmSystem &mmsys) {
  //cout << "dbg1130 test0 : " << endl;
  int buf = 0;
  read_bin_values(&buf, 1);
  char header[MAX_LEN_NAME];
  ifs.read(header, buf);
  //cout << "dbg1130 group_coord : " << string(header) << endl;
  int aus_type = 0;
  read_bin_values(&aus_type, 1);
  //cout << "dbg1130 aus_type: " << aus_type << endl;
  int n_groups = 0;
  read_bin_values(&n_groups, 1);
  //cout << "dbg1130 n_groups: " << n_groups << endl;
  
  vector<int> enhance_groups;
  for (int i = 0; i < n_groups; i++) {
    read_bin_values(&buf, 1);
    enhance_groups.push_back(buf);
  }
  for (int i = 0; i < n_groups; i++) {
    read_bin_values(&buf, 1);
    // if(buf != mmsys.n_atoms_in_groups[enhance_groups[i]]){
    // stringstream ss;
    // ss << "Information in the V-AUS restart file is inconsistent"<<endl;
    // ss << "Enhanced group " << i << " (atom group " << enhance_groups[i] << ") " << endl;
    // ss << mmsys.n_atoms_in_groups[enhance_groups[i]] << " atoms in the group definition." << endl;
    // ss << buf << " atoms in the V-AUS restart file." << endl;
    // error_exit(ss.str(), "1A00005");
    //}
  }
  //cout << "dbg1130 test1: " << endl;  
  if(mmsys.extended_mode == EXTENDED_VAUS){
    mmsys.vmcmd->set_aus_type(aus_type);
    mmsys.vmcmd->set_enhance_groups(mmsys.n_atoms_in_groups,
				    mmsys.atom_groups, 
				    n_groups, enhance_groups);
    
    for (int i = 0; i < n_groups; i++) {
      for (int j = 0; j < mmsys.n_atoms_in_groups[enhance_groups[i]]; j++) {
	read_bin_values(&(mmsys.vmcmd->get_crd_groups()[i][j][0]), 1);
	read_bin_values(&(mmsys.vmcmd->get_crd_groups()[i][j][1]), 1);
	read_bin_values(&(mmsys.vmcmd->get_crd_groups()[i][j][2]), 1);
      }
    }
  }else if(mmsys.extended_mode == EXTENDED_VCMD){
    mmsys.vcmd->set_reactcrd_type(aus_type);
    mmsys.vcmd->set_enhance_groups(mmsys.n_atoms_in_groups,
				   mmsys.atom_groups, 
				   n_groups, enhance_groups);
    for (int i = 0; i < n_groups; i++) {
      for (int j = 0; j < mmsys.n_atoms_in_groups[enhance_groups[i]]; j++) {
	read_bin_values(&(mmsys.vcmd->get_crd_groups()[i][j][0]), 1);
	read_bin_values(&(mmsys.vcmd->get_crd_groups()[i][j][1]), 1);
	read_bin_values(&(mmsys.vcmd->get_crd_groups()[i][j][2]), 1);
      }
    }
  }
  //cout << "dbg1130 test2: " << n_groups << endl;  
  return 0;
}

int Read::load_ls_vcmd(MmSystem &mmsys) {
  int interval;
  read_bin_values(&interval, 1);
  int dim;
  read_bin_values(&dim, 1);
  
  //cout << "dbg 0303 intrv: " << interval << " dim: " << dim << endl;;
  //cout << "dbg 0304 read a inter: " << interval << " dim:" << dim << endl;
  mmsys.vcmd->set_n_dim(dim);
  mmsys.vcmd->set_trans_interval(interval);
  //int n_states = 1;
  std::vector< std::vector<int> > grp_id;
  for (int d = 0; d < dim; d++) {
    int n_vs;
    read_bin_values(&n_vs, 1);
    int n_grp;
    read_bin_values(&n_grp, 1);
    std::vector<int> grp_dim;
    std::vector<string> grp_name;
    for(int i = 0; i < n_grp; i++){
      int grp_id;
      read_bin_values(&grp_id, 1);
      grp_dim.push_back(grp_id);
      int len;
      read_bin_values(&len, 1);
      char name[MAX_LEN_NAME];
      ifs.read(name, len);
      grp_name.push_back(string(name));
      //cout << "read 0304 d:" <<d << " grp:"<< grp_id << ": " << name << endl;
    }
    mmsys.vcmd->push_grp_ids_name(grp_dim, grp_name);
    //n_states *= n_vs;
    std::vector<real> range_min;
    std::vector<real> range_max;
    for (int i = 0; i < n_vs; i++) {
      double buf_min;
      double buf_max;
      read_bin_values(&buf_min, 1);
      read_bin_values(&buf_max, 1);
      range_min.push_back(buf_min);
      range_max.push_back(buf_max);
    }
    mmsys.vcmd->push_vs_range(range_min, range_max);
  }
  //cout << "dbg 0304 read b" << endl;
  std::vector<int> init_vs;
  for (int d = 0; d < dim; d++) {
    int vs;
    read_bin_values(&vs, 1);
    init_vs.push_back(vs - 1);
  }
  mmsys.vcmd->set_init_vs(init_vs);
  int seed;
  read_bin_values(&seed, 1);
  mmsys.vcmd->set_random_seed(seed);
  //cout << "dbg 0304 read c" << endl;
  int n_states;
  read_bin_values(&n_states, 1);
  std::map< std::vector<int>, real > q_cano;
  for (int i = 0; i < n_states; i++){
    std::vector<int> state;
    for(int d = 0; d < dim; d++){
      int vs;
      read_bin_values(&vs, 1);      
      state.push_back(vs-1);
    }
    double param;
    read_bin_values(&param, 1);
    q_cano[state] = param;
  }
  //cout << "dbg 0304 read d" << endl;
  mmsys.vcmd->set_q_cano(q_cano);
  //cout << "dbg 0304 read e" << endl;
  return 0;
}
int Read::load_ls_pote(MmSystem &mmsys) {
  double pote;
  read_bin_values(&pote, 1);
  mmsys.potential_e = pote;
  cout << "dbg0803a " << mmsys.potential_e << endl;
  return 0;
}
template <typename TYPE>
int Read::read_bin_values(TYPE *recept, int len) {
    ifs.read((char *)recept, sizeof(TYPE) * len);
    if (is_conv_endian()) {
        for (int i = 0; i < len; i++) { recept[i] = reverse_endian(recept[i]); }
    }
    return 0;
}
