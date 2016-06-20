#include "WriteTrr.h"
using namespace std;

WriteTrr::WriteTrr() : Write() {}

WriteTrr::~WriteTrr() {}

int WriteTrr::write_trr(int       n_atoms,
                        int       cur_step,
                        real      cur_time,
                        real      lx,
                        real      ly,
                        real      lz,
                        real **   crd,
                        real **   vel_just,
                        real_fc **force,
                        real      cpu_time,
                        real      total_e,
                        real      kinetic_e,
                        real      temperature,
                        real      potential_e,
                        real      vdw_e,
                        bool      out_box,
                        bool      out_crd,
                        bool      out_vel,
                        bool      out_force,
                        int       n_atoms_group,
                        int *     atom_group) {
    return 0;
}

WriteTrrGromacs::WriteTrrGromacs() : WriteTrr() {}

WriteTrrGromacs::~WriteTrrGromacs() {}

int WriteTrrGromacs::write_trr(int       n_atoms,
                               int       cur_step,
                               real      cur_time,
                               real      lx,
                               real      ly,
                               real      lz,
                               real **   crd,
                               real **   vel_just,
                               real_fc **force,
                               real      cpu_time,
                               real      total_e,
                               real      kinetic_e,
                               real      temperature,
                               real      potential_e,
                               real      vdw_e,
                               bool      out_box,
                               bool      out_crd,
                               bool      out_vel,
                               bool      out_force,
                               int       n_atoms_group,
                               int *     atom_group) {

    int box_size          = 0;
    if (out_box) box_size = 9 * sizeof(real);
    int x_size            = 0;
    if (out_crd) x_size   = n_atoms * 3 * sizeof(real);
    int v_size            = 0;
    if (out_vel) v_size   = n_atoms * 3 * sizeof(real);
    int f_size            = 0;
    if (out_force) f_size = n_atoms * 3 * sizeof(real);

    int  magic      = 1993;
    int  nchar1     = 13;
    int  nchar2     = 12;
    real dummy_real = 0.0;
    int  dummy      = 0;

    if (!out_crd && !out_vel && !out_force) { return 1; }
    ofs.write((const char *)&magic, sizeof magic);
    ofs.write((const char *)&nchar1, sizeof(int));
    ofs.write((const char *)&nchar2, sizeof(int));
    ofs.write("GMX_trn_file", 12);
    ofs.write((const char *)&dummy, sizeof(int));
    ofs.write((const char *)&dummy, sizeof(int));
    ofs.write((const char *)&box_size, sizeof box_size);
    ofs.write((const char *)&dummy, sizeof(int));
    ofs.write((const char *)&dummy, sizeof(int));
    ofs.write((const char *)&dummy, sizeof(int));
    ofs.write((const char *)&dummy, sizeof(int));
    ofs.write((const char *)&x_size, sizeof x_size);
    ofs.write((const char *)&v_size, sizeof v_size);
    ofs.write((const char *)&f_size, sizeof f_size);
    ofs.write((const char *)&n_atoms, sizeof(int));
    ofs.write((const char *)&cur_step, sizeof(int));
    ofs.write((const char *)&dummy, sizeof(int));
    ofs.write((const char *)&cur_time, sizeof(real));
    ofs.write((const char *)&dummy_real, sizeof(real));

    if (out_box) {
        real x = lx * 0.1;
        real y = ly * 0.1;
        real z = lz * 0.1;

        ofs.write((const char *)&x, sizeof x);
        ofs.write((const char *)&dummy_real, sizeof(real));
        ofs.write((const char *)&dummy_real, sizeof(real));
        ofs.write((const char *)&dummy_real, sizeof(real));
        ofs.write((const char *)&y, sizeof y);
        ofs.write((const char *)&dummy_real, sizeof(real));
        ofs.write((const char *)&dummy_real, sizeof(real));
        ofs.write((const char *)&dummy_real, sizeof(real));
        ofs.write((const char *)&z, sizeof z);
    }
    if (out_crd) {
        if (n_atoms_group == 0) {
            for (int atomid = 0; atomid < n_atoms; atomid++) {
                real x = crd[atomid][0] * 0.1;
                real y = crd[atomid][1] * 0.1;
                real z = crd[atomid][2] * 0.1;
                ofs.write((const char *)&x, sizeof(real));
                ofs.write((const char *)&y, sizeof(real));
                ofs.write((const char *)&z, sizeof(real));
            }
        } else {
            for (int i = 0; i < n_atoms_group; i++) {
                real x = crd[atom_group[i]][0] * 0.1;
                real y = crd[atom_group[i]][1] * 0.1;
                real z = crd[atom_group[i]][2] * 0.1;
                ofs.write((const char *)&x, sizeof(real));
                ofs.write((const char *)&y, sizeof(real));
                ofs.write((const char *)&z, sizeof(real));
            }
        }
    }
    if (out_vel) {
        if (n_atoms_group == 0) {
            for (int atomid = 0; atomid < n_atoms; atomid++) {
                real x = vel_just[atomid][0] * 0.1;
                real y = vel_just[atomid][1] * 0.1;
                real z = vel_just[atomid][2] * 0.1;
                ofs.write((const char *)&x, sizeof(real));
                ofs.write((const char *)&y, sizeof(real));
                ofs.write((const char *)&z, sizeof(real));
            }
        } else {
            for (int i = 0; i < n_atoms_group; i++) {
                real x = vel_just[atom_group[i]][0] * 0.1;
                real y = vel_just[atom_group[i]][1] * 0.1;
                real z = vel_just[atom_group[i]][2] * 0.1;
                ofs.write((const char *)&x, sizeof(real));
                ofs.write((const char *)&y, sizeof(real));
                ofs.write((const char *)&z, sizeof(real));
            }
        }
    }
    if (out_force) {
        if (n_atoms_group == 0) {
            for (int atomid = 0; atomid < n_atoms; atomid++) {
                real x = force[atomid][0] * 0.1;
                real y = force[atomid][1] * 0.1;
                real z = force[atomid][2] * 0.1;
                ofs.write((const char *)&x, sizeof(real));
                ofs.write((const char *)&y, sizeof(real));
                ofs.write((const char *)&z, sizeof(real));
            }
        } else {
            for (int i = 0; i < n_atoms_group; i++) {
                real x = force[atom_group[i]][0] * 0.1;
                real y = force[atom_group[i]][1] * 0.1;
                real z = force[atom_group[i]][2] * 0.1;
                ofs.write((const char *)&x, sizeof(real));
                ofs.write((const char *)&y, sizeof(real));
                ofs.write((const char *)&z, sizeof(real));
            }
        }
    }
    return 0;
}

WriteTrrPresto::WriteTrrPresto() : WriteTrr() {}

WriteTrrPresto::~WriteTrrPresto() {}

int WriteTrrPresto::write_trr(int       n_atoms,
                              int       cur_step,
                              real      cur_time,
                              real      lx,
                              real      ly,
                              real      lz,
                              real **   crd,
                              real **   vel_just,
                              real_fc **force,
                              real      cpu_time,
                              real      total_e,
                              real      kinetic_e,
                              real      temperature,
                              real      potential_e,
                              real      vdw_e,
                              bool      out_box,
                              bool      out_crd,
                              bool      out_vel,
                              bool      out_force,
                              int       n_atoms_group,
                              int *     atom_group) {
    int buf = 44;
    ofs.write((const char *)&buf, sizeof(int));
    ofs.write((const char *)&cur_step, sizeof(int));
    float buf_cur_time = cur_time;
    ofs.write((const char *)&buf_cur_time, sizeof(float));
    float buf_f = 0.0;
    buf_f       = (float)cpu_time;
    ofs.write((const char *)&buf_f, sizeof(float)); // cpu_time
    buf_f = (float)total_e;
    ofs.write((const char *)&buf_f, sizeof(float)); // total e
    buf_f = (float)kinetic_e;
    ofs.write((const char *)&buf_f, sizeof(float)); // kinetic e
    buf_f = (float)temperature;
    ofs.write((const char *)&buf_f, sizeof(float)); // temperature
    buf_f = (float)potential_e;
    ofs.write((const char *)&buf_f, sizeof(float)); // potential e
    buf_f = 0.0;
    ofs.write((const char *)&buf_f, sizeof(float)); // rmsf
    buf_f = (float)vdw_e;
    ofs.write((const char *)&buf_f, sizeof(float)); // vdw
    buf_f = 0.0;
    ofs.write((const char *)&buf_f, sizeof(float)); // hyd
    buf_f = 0.0;
    ofs.write((const char *)&buf_f, sizeof(float)); // rmsd
    buf = 0;
    ofs.write((const char *)&buf, sizeof(int));

    if (n_atoms_group == 0) {
        buf = n_atoms * 3 * 4;
        ofs.write((const char *)&buf, sizeof(int));
        for (int i = 0; i < n_atoms; i++) {
            float x = crd[i][0];
            float y = crd[i][1];
            float z = crd[i][2];
            ofs.write((const char *)&x, sizeof(float));
            ofs.write((const char *)&y, sizeof(float));
            ofs.write((const char *)&z, sizeof(float));
        }
        ofs.write((const char *)&buf, sizeof(int));
    } else {
        buf = n_atoms_group * 3 * 4;
        ofs.write((const char *)&buf, sizeof(int));
        for (int i = 0; i < n_atoms_group; i++) {
            float x = crd[atom_group[i]][0];
            float y = crd[atom_group[i]][1];
            float z = crd[atom_group[i]][2];
            ofs.write((const char *)&x, sizeof(float));
            ofs.write((const char *)&y, sizeof(float));
            ofs.write((const char *)&z, sizeof(float));
        }
        ofs.write((const char *)&buf, sizeof(int));
    }

    return 0;
}

WriteRestart::WriteRestart() : Write() {}
WriteRestart::~WriteRestart() {}
int WriteRestart::write_restart(int    n_atoms,
                                int    n_steps,
                                double time,
                                double e_potential,
                                double e_kinetic,
                                real **crd,
                                real **vel) {
    int    buf_int;
    double buf_dbl;
    open();
    // title
    buf_int = 80;
    char title[80];
    strcpy(title, ABOUT_ME.c_str());

    ofs.write((const char *)&buf_int, sizeof(int));
    ofs.write(title, 80);
    ofs.write((const char *)&buf_int, sizeof(int));
    // #atoms, #velocities
    buf_int = 8;
    ofs.write((const char *)&buf_int, sizeof(int));
    ofs.write((const char *)&n_atoms, sizeof(int));
    ofs.write((const char *)&n_atoms, sizeof(int));
    ofs.write((const char *)&buf_int, sizeof(int));
    //
    buf_int = 36;
    ofs.write((const char *)&buf_int, sizeof(int));
    ofs.write((const char *)&n_steps, sizeof(int));
    ofs.write((const char *)&time, sizeof(double));
    buf_dbl = e_kinetic + e_potential;
    ofs.write((const char *)&buf_dbl, sizeof(double));
    ofs.write((const char *)&e_kinetic, sizeof(double));
    ofs.write((const char *)&e_potential, sizeof(double));
    ofs.write((const char *)&buf_int, sizeof(int));

    buf_int = n_atoms * 3 * 8;
    ofs.write((const char *)&buf_int, sizeof(int));
    for (int i = 0; i < n_atoms; i++) {
        for (int d = 0; d < 3; d++) {
            buf_dbl = (double)crd[i][d];
            ofs.write((const char *)&buf_dbl, sizeof(double));
        }
    }
    ofs.write((const char *)&buf_int, sizeof(int));

    ofs.write((const char *)&buf_int, sizeof(int));
    for (int i = 0; i < n_atoms; i++) {
        for (int d = 0; d < 3; d++) {
            buf_dbl = (double)vel[i][d];
            ofs.write((const char *)&buf_dbl, sizeof(double));
        }
    }
    ofs.write((const char *)&buf_int, sizeof(int));
    close();
    return 0;
}

WriteGroupCoord::WriteGroupCoord() {}
WriteGroupCoord::~WriteGroupCoord() {}

int WriteGroupCoord::write_aus_restart(const int   aus_type,
                                       int         n_enhance_groups,
                                       vector<int> enhance_groups,
                                       int *       n_atoms_in_groups,
                                       real ***    crd_groups) {

    int buf = 5;
    ofs.write((const char *)&buf, sizeof(int));
    ofs.write("V-AUS", 5);

    ofs.write((const char *)&aus_type, sizeof(int));
    ofs.write((const char *)&n_enhance_groups, sizeof(int));
    for (int k = 0; k < n_enhance_groups; k++) { ofs.write((const char *)&enhance_groups[k], sizeof(int)); }
    for (int k = 0; k < n_enhance_groups; k++) {
        ofs.write((const char *)&n_atoms_in_groups[enhance_groups[k]], sizeof(int));
    }

    for (int k = 0; k < n_enhance_groups; k++) {
        for (int i = 0; i < n_atoms_in_groups[enhance_groups[k]]; i++) {
            double x = crd_groups[k][i][0];
            double y = crd_groups[k][i][1];
            double z = crd_groups[k][i][2];
            ofs.write((const char *)&x, sizeof(double));
            ofs.write((const char *)&y, sizeof(double));
            ofs.write((const char *)&z, sizeof(double));
        }
    }
    return 0;
}
