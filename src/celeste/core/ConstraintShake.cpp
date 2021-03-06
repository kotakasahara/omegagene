#include "ConstraintShake.h"
using namespace std;

ConstraintShake::ConstraintShake() : ConstraintObject() {}

ConstraintShake::~ConstraintShake() {
    free_constraint();
}

int ConstraintShake::apply_constraint(real *in_crd, real *in_crd_prev, real_pw *in_mass_inv, PBC *pbc) {
    shake_pair(in_crd, in_crd_prev, in_mass_inv, pbc);
    shake_trio(in_crd, in_crd_prev, in_mass_inv, pbc);
    shake_quad(in_crd, in_crd_prev, in_mass_inv, pbc);

    return 0;
}

int ConstraintShake::shake_pair(real *in_crd, real *in_crd_prev, real_pw *in_mass_inv, PBC *pbc) {
    for (int i_cst = 0; i_cst < n_pair; i_cst++) {
        int      atomid1_3    = pair_atomids[i_cst][0] * 3;
        int      atomid2_3    = pair_atomids[i_cst][1] * 3;
        real     crd1[3]      = {in_crd[atomid1_3], in_crd[atomid1_3 + 1], in_crd[atomid1_3 + 2]};
        real     crd2[3]      = {in_crd[atomid2_3], in_crd[atomid2_3 + 1], in_crd[atomid2_3 + 2]};
        real     crd1_prev[3] = {in_crd_prev[atomid1_3], in_crd_prev[atomid1_3 + 1], in_crd_prev[atomid1_3 + 2]};
        real     crd2_prev[3] = {in_crd_prev[atomid2_3], in_crd_prev[atomid2_3 + 1], in_crd_prev[atomid2_3 + 2]};
        real_cst d12[3];
        real_cst d12_prev[3];
        pbc->diff_crd_minim_image(d12, crd2, crd1);
        pbc->diff_crd_minim_image(d12_prev, crd2_prev, crd1_prev);
        real_cst mass1_inv    = in_mass_inv[pair_atomids[i_cst][0]];
        real_cst mass2_inv    = in_mass_inv[pair_atomids[i_cst][1]];
        real_cst mass12_inv   = mass1_inv + mass2_inv;
        real_cst diff_sq      = d12[0] * d12[0] + d12[1] * d12[1] + d12[2] * d12[2];
        real_cst diff_sq_prev = d12_prev[0] * d12_prev[0] + d12_prev[1] * d12_prev[1] + d12_prev[2] * d12_prev[2];
        real_cst diff_ip      = d12[0] * d12_prev[0] + d12[1] * d12_prev[1] + d12[2] * d12_prev[2];

        real_cst val_a = mass12_inv * mass12_inv * diff_sq_prev;
        real_cst val_b = mass12_inv * diff_ip;
        real_cst val_c = diff_sq - pair_dist[i_cst];
        real_cst val_g = (val_b - sqrt(max(val_b * val_b - val_c * val_a, (real_cst)0.0))) / val_a;

        for (int d = 0; d < 3; d++) {
            in_crd[atomid1_3 + d] += mass1_inv * val_g * d12_prev[d];
            in_crd[atomid2_3 + d] -= mass2_inv * val_g * d12_prev[d];
        }
    }

    return 0;
}

int ConstraintShake::shake_trio(real *in_crd, real *in_crd_prev, real_pw *in_mass_inv, PBC *pbc) {
    for (int i_cst = 0; i_cst < n_trio; i_cst++) {
        int      atomid1_3    = trio_atomids[i_cst][0] * 3;
        int      atomid2_3    = trio_atomids[i_cst][1] * 3;
        int      atomid3_3    = trio_atomids[i_cst][2] * 3;
        real     crd1[3]      = {in_crd[atomid1_3], in_crd[atomid1_3 + 1], in_crd[atomid1_3 + 2]};
        real     crd2[3]      = {in_crd[atomid2_3], in_crd[atomid2_3 + 1], in_crd[atomid2_3 + 2]};
        real     crd3[3]      = {in_crd[atomid3_3], in_crd[atomid3_3 + 1], in_crd[atomid3_3 + 2]};
        real     crd1_prev[3] = {in_crd_prev[atomid1_3], in_crd_prev[atomid1_3 + 1], in_crd_prev[atomid1_3 + 2]};
        real     crd2_prev[3] = {in_crd_prev[atomid2_3], in_crd_prev[atomid2_3 + 1], in_crd_prev[atomid2_3 + 2]};
        real     crd3_prev[3] = {in_crd_prev[atomid3_3], in_crd_prev[atomid3_3 + 1], in_crd_prev[atomid3_3 + 2]};
        real_cst d12[3], d23[3], d31[3];
        real_cst d12_prev[3], d23_prev[3], d31_prev[3];

        real_cst g[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        pbc->diff_crd_minim_image(d12, crd1, crd2);
        pbc->diff_crd_minim_image(d12_prev, crd1_prev, crd2_prev);
        pbc->diff_crd_minim_image(d23, crd2, crd3);
        pbc->diff_crd_minim_image(d23_prev, crd2_prev, crd3_prev);
        pbc->diff_crd_minim_image(d31, crd3, crd1);
        pbc->diff_crd_minim_image(d31_prev, crd3_prev, crd1_prev);

        real_cst mass1_inv = in_mass_inv[trio_atomids[i_cst][0]];
        real_cst mass2_inv = in_mass_inv[trio_atomids[i_cst][1]];
        real_cst mass3_inv = in_mass_inv[trio_atomids[i_cst][2]];

        real_cst coef11 = mass1_inv + mass2_inv;
        real_cst coef12 = -mass2_inv;
        real_cst coef13 = -mass1_inv;
        real_cst coef21 = -mass2_inv;
        real_cst coef22 = mass2_inv + mass3_inv;
        real_cst coef23 = -mass3_inv;
        real_cst coef31 = -mass1_inv;
        real_cst coef32 = -mass3_inv;
        real_cst coef33 = mass3_inv + mass1_inv;

        bool converge = false;
        for (int i_loop = 0; i_loop < max_loops; i_loop++) {
            for (int d = 0; d < 3; d++) {
                d12[d] += (coef11 * d12_prev[d] * g[0] + coef12 * d23_prev[d] * g[1] + coef13 * d31_prev[d] * g[2]);
                d23[d] += (coef21 * d12_prev[d] * g[0] + coef22 * d23_prev[d] * g[1] + coef23 * d31_prev[d] * g[2]);
                d31[d] += (coef31 * d12_prev[d] * g[0] + coef32 * d23_prev[d] * g[1] + coef33 * d31_prev[d] * g[2]);
            }
            real_cst gradf11 = 2.0 * coef11 * (d12[0] * d12_prev[0] + d12[1] * d12_prev[1] + d12[2] * d12_prev[2]);
            real_cst gradf12 = 2.0 * coef12 * (d12[0] * d23_prev[0] + d12[1] * d23_prev[1] + d12[2] * d23_prev[2]);
            real_cst gradf13 = 2.0 * coef13 * (d12[0] * d31_prev[0] + d12[1] * d31_prev[1] + d12[2] * d31_prev[2]);
            real_cst gradf21 = 2.0 * coef21 * (d23[0] * d12_prev[0] + d23[1] * d12_prev[1] + d23[2] * d12_prev[2]);
            real_cst gradf22 = 2.0 * coef22 * (d23[0] * d23_prev[0] + d23[1] * d23_prev[1] + d23[2] * d23_prev[2]);
            real_cst gradf23 = 2.0 * coef23 * (d23[0] * d31_prev[0] + d23[1] * d31_prev[1] + d23[2] * d31_prev[2]);
            real_cst gradf31 = 2.0 * coef31 * (d31[0] * d12_prev[0] + d31[1] * d12_prev[1] + d31[2] * d12_prev[2]);
            real_cst gradf32 = 2.0 * coef32 * (d31[0] * d23_prev[0] + d31[1] * d23_prev[1] + d31[2] * d23_prev[2]);
            real_cst gradf33 = 2.0 * coef33 * (d31[0] * d31_prev[0] + d31[1] * d31_prev[1] + d31[2] * d31_prev[2]);

            real_cst diff[3] = {trio_dist[i_cst][0] - (d12[0] * d12[0] + d12[1] * d12[1] + d12[2] * d12[2]),
                                trio_dist[i_cst][1] - (d23[0] * d23[0] + d23[1] * d23[1] + d23[2] * d23[2]),
                                trio_dist[i_cst][2] - (d31[0] * d31[0] + d31[1] * d31[1] + d31[2] * d31[2])};

            if (fabs(diff[0] / trio_dist[i_cst][0]) < tolerance && fabs(diff[1] / trio_dist[i_cst][1]) < tolerance
                && fabs(diff[2] / trio_dist[i_cst][2]) < tolerance) {
                converge = true;
                break;
            }
            real_cst determin =
                1.0 / (gradf11 * gradf22 * gradf33 + gradf12 * gradf23 * gradf31 + gradf13 * gradf21 * gradf32
                       - gradf11 * gradf23 * gradf32 - gradf12 * gradf21 * gradf33 - gradf13 * gradf22 * gradf31);
            real_cst detb11 = gradf22 * gradf33 - gradf23 * gradf32;
            real_cst detb12 = gradf32 * gradf13 - gradf33 * gradf12;
            real_cst detb13 = gradf12 * gradf23 - gradf13 * gradf22;
            real_cst detb21 = gradf23 * gradf31 - gradf21 * gradf33;
            real_cst detb22 = gradf33 * gradf11 - gradf31 * gradf13;
            real_cst detb23 = gradf13 * gradf21 - gradf11 * gradf23;
            real_cst detb31 = gradf21 * gradf32 - gradf22 * gradf31;
            real_cst detb32 = gradf31 * gradf12 - gradf32 * gradf11;
            real_cst detb33 = gradf11 * gradf22 - gradf12 * gradf21;

            g[0] = (detb11 * diff[0] + detb12 * diff[1] + detb13 * diff[2]) * determin;
            g[1] = (detb21 * diff[0] + detb22 * diff[1] + detb23 * diff[2]) * determin;
            g[2] = (detb31 * diff[0] + detb32 * diff[1] + detb33 * diff[2]) * determin;
            g[3] += g[0];
            g[4] += g[1];
            g[5] += g[2];
        }

        if (!converge) {
            cout << "Shake not converged, tri: " << i_cst << endl;
            cout << trio_atomids[i_cst][0] << " - " << trio_atomids[i_cst][1] << " - " << trio_atomids[i_cst][2]
                 << endl;
        }

        for (int d = 0; d < 3; d++) {
            in_crd[atomid1_3 + d] += mass1_inv * (g[3] * d12_prev[d] - g[5] * d31_prev[d]);
            in_crd[atomid2_3 + d] += mass2_inv * (-g[3] * d12_prev[d] + g[4] * d23_prev[d]);
            in_crd[atomid3_3 + d] += mass3_inv * (-g[4] * d23_prev[d] + g[5] * d31_prev[d]);
        }
    }

    return 0;
}

int ConstraintShake::shake_quad(real *in_crd, real *in_crd_prev, real_pw *in_mass_inv, PBC *pbc) {
    for (int i_cst = 0; i_cst < n_quad; i_cst++) {
        int atomid_3[4] = {quad_atomids[i_cst][0] * 3, quad_atomids[i_cst][1] * 3, quad_atomids[i_cst][2] * 3,
                           quad_atomids[i_cst][3] * 3};

        real_cst mass_inv[4] = {in_mass_inv[quad_atomids[i_cst][0]], in_mass_inv[quad_atomids[i_cst][1]],
                                in_mass_inv[quad_atomids[i_cst][2]], in_mass_inv[quad_atomids[i_cst][3]]};

        real crd_cur[4][3];
        real crd_prev[4][3];
        for (int i = 0; i < 4; i++) {
            for (int d = 0; d < 3; d++) {
                crd_cur[i][d]  = in_crd[atomid_3[i] + d];
                crd_prev[i][d] = in_crd_prev[atomid_3[i] + d];
            }
        }
        real_cst weight[6][6] = {
            {mass_inv[0] + mass_inv[1], -mass_inv[1], -mass_inv[0], -mass_inv[0], 0.0, -mass_inv[1]},
            {-mass_inv[1], mass_inv[1] + mass_inv[2], -mass_inv[2], 0.0, -mass_inv[2], mass_inv[1]},
            {-mass_inv[0], -mass_inv[2], mass_inv[2] + mass_inv[0], mass_inv[0], mass_inv[2], 0.0},
            {-mass_inv[0], 0.0, mass_inv[0], mass_inv[3] + mass_inv[0], -mass_inv[3], -mass_inv[3]},
            {0.0, -mass_inv[2], mass_inv[2], -mass_inv[3], mass_inv[2] + mass_inv[3], mass_inv[3]},
            {-mass_inv[1], mass_inv[1], 0.0, -mass_inv[3], mass_inv[3], mass_inv[1] + mass_inv[3]}};

        // 1-2 calculate distances between atoms
        real_cst d_cur[6][3];
        real_cst d_prev[6][3];
        for (int i_pair = 0; i_pair < 6; i_pair++) {
            // int id1 = quad_atomids[i_cst][Pairs_idx[i_pair][0]];
            // int id2 = quad_atomids[i_cst][Pairs_idx[i_pair][1]];
            // int id1 = Pairs_idx[i_pair][0];
            // int id2 = Pairs_idx[i_pair][1];
            /*
            for (int d=0; d < 3; d++)
            d_cur[i_pair][d] = crd_cur[Pairs_idx[i_pair][0]][d] -
          crd_cur[Pairs_idx[i_pair][1]][d];
            for (int d=0; d < 3; d++)
            d_prev[i_pair][d] = crd_prev[Pairs_idx[i_pair][0]][d] -
          crd_prev[Pairs_idx[i_pair][1]][d];
            */
            pbc->diff_crd_minim_image(d_cur[i_pair], crd_cur[Pairs_idx[i_pair][0]], crd_cur[Pairs_idx[i_pair][1]]);

            pbc->diff_crd_minim_image(d_prev[i_pair], crd_prev[Pairs_idx[i_pair][0]], crd_prev[Pairs_idx[i_pair][1]]);
        }
        real_cst coef[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        real_cst coef_post[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        // 1-3 convergence loop
        bool     converge = false;
        real_cst d_virtual[6][3];
        // int dbgatom = 151;
        for (int i_loop = 0; i_loop < max_loops; i_loop++) {
            // 1-3-1 calculate virtual vector
            // d_virtual = vec + sum(weight * d_prev * coef)
            for (int i = 0; i < 6; i++) coef[i] = coef_post[i];

            // if(atomid_3[0] == dbgatom*3){
            // cout << "DBG1 coef: ";
            // for(int i=0; i<6; i++) cout << coef[i] << " ";
            // cout << endl;
            //}

            for (int i_pair = 0; i_pair < 6; i_pair++) {
                for (int d = 0; d < 3; d++) d_virtual[i_pair][d] = d_cur[i_pair][d];
                for (int j_pair = 0; j_pair < 6; j_pair++) {
                    for (int d = 0; d < 3; d++)
                        d_virtual[i_pair][d] += weight[j_pair][i_pair] * d_prev[j_pair][d] * coef[j_pair];
                }
            }
            real_cst grad[6][6];
            for (int i_pair = 0; i_pair < 6; i_pair++) {
                for (int j_pair = 0; j_pair < 6; j_pair++) {
                    grad[i_pair][j_pair] = 2.0 * weight[j_pair][i_pair] * (d_virtual[j_pair][0] * d_prev[i_pair][0]
                                                                           + d_virtual[j_pair][1] * d_prev[i_pair][1]
                                                                           + d_virtual[j_pair][2] * d_prev[i_pair][2]);
                    //	  if(atomid_3[0] == dbgatom*3)
                    // cout << "grad["<<i_pair<<"]["<<j_pair<<"] " << grad[i_pair][j_pair]<<endl;
                }
            }
            converge = true;
            real_cst sqrt_dist[6];
            for (int i_pair = 0; i_pair < 6; i_pair++) {
                real_cst norm =
                    (d_virtual[i_pair][0] * d_virtual[i_pair][0] + d_virtual[i_pair][1] * d_virtual[i_pair][1]
                     + d_virtual[i_pair][2] * d_virtual[i_pair][2]);
                sqrt_dist[i_pair] = quad_dist[i_cst][i_pair] - norm;
                // cout << "shkdist " << i_cst << " "<<i_pair << " " << quad_dist[i_cst][i_pair] << " " << norm << endl;
                real_cst diff = fabs(sqrt_dist[i_pair] / quad_dist[i_cst][i_pair]);
                converge      = converge && (diff < tolerance);
            }
            if (converge) break;
            // 1-3-4 calculate right hand vector
            real_cst result[6];
            //	if(atomid_3[0] == dbgatom*3)
            // cout << "DBG result: ";
            for (int i_pair = 0; i_pair < 6; i_pair++) {
                result[i_pair] = sqrt_dist[i_pair];
                for (int j_pair = 0; j_pair < 6; j_pair++) { result[i_pair] += grad[j_pair][i_pair] * coef[j_pair]; }
                // if(atomid_3[0] == dbgatom*3)
                // cout << result[i_pair] << " ";
            }
            // if(atomid_3[0] == dbgatom*3)
            // cout << endl;
            // 1-3-5 solve linear equation (grad * coef_post = result)
            calc_linear_eq(grad, coef_post, result, 6);
        }

        if (!converge) {
            cout << "Shake not converged, quad: " << i_cst << endl;
            cout << quad_atomids[i_cst][0] << " - " << quad_atomids[i_cst][1] << " - " << quad_atomids[i_cst][2]
                 << " - " << quad_atomids[i_cst][3] << endl;
        }
        for (int d = 0; d < 3; d++) {
            in_crd[atomid_3[0] + d] +=
                mass_inv[0] * (coef[0] * d_prev[0][d] - coef[2] * d_prev[2][d] - coef[3] * d_prev[3][d]);
            in_crd[atomid_3[1] + d] +=
                mass_inv[1] * (-coef[0] * d_prev[0][d] + coef[1] * d_prev[1][d] + coef[5] * d_prev[5][d]);
            in_crd[atomid_3[2] + d] +=
                mass_inv[2] * (-coef[1] * d_prev[1][d] + coef[2] * d_prev[2][d] + coef[4] * d_prev[4][d]);
            in_crd[atomid_3[3] + d] +=
                mass_inv[3] * (coef[3] * d_prev[3][d] - coef[4] * d_prev[4][d] - coef[5] * d_prev[5][d]);
        }
    }

    return 0;
}

ConstraintSettle::ConstraintSettle() : ConstraintObject() {}

ConstraintSettle::~ConstraintSettle() {
    free_constraint();
}

int ConstraintSettle::apply_constraint(real *in_crd, real *in_crd_prev, real_pw *in_mass_inv, PBC *pbc) {
    return 0;
}
ConstraintRattle::ConstraintRattle() : ConstraintObject() {}

ConstraintRattle::~ConstraintRattle() {
    free_constraint();
}

int ConstraintRattle::apply_rattle(real *in_crd, real *in_crd_prev, real *in_vel, real *in_work, real_pw *in_mass_inv, PBC *pbc) {
  rattle_pair(in_crd, in_crd_prev, in_vel, in_work, in_mass_inv, pbc);
  rattle_trio(in_crd, in_crd_prev, in_vel, in_work, in_mass_inv, pbc);
  rattle_quad(in_crd, in_crd_prev, in_vel, in_work, in_mass_inv, pbc);
  return 0;
}

int ConstraintRattle::rattle_pair(real *in_crd, real *in_crd_prev, real *in_vel, real *in_work, real_pw *in_mass_inv, PBC *pbc) {
  for (int i_cst = 0; i_cst < n_pair; i_cst++) {
    int      atomid1_3    = pair_atomids[i_cst][0] * 3;
    int      atomid2_3    = pair_atomids[i_cst][1] * 3;
    real     crd1[3]      = {in_crd[atomid1_3], in_crd[atomid1_3 + 1], in_crd[atomid1_3 + 2]};
    real     crd2[3]      = {in_crd[atomid2_3], in_crd[atomid2_3 + 1], in_crd[atomid2_3 + 2]};
    real     crd1_prev[3] = {in_crd_prev[atomid1_3], in_crd_prev[atomid1_3 + 1], in_crd_prev[atomid1_3 + 2]};
    real     crd2_prev[3] = {in_crd_prev[atomid2_3], in_crd_prev[atomid2_3 + 1], in_crd_prev[atomid2_3 + 2]};
    real_cst d12[3];
    real_cst d12_prev[3];
    pbc->diff_crd_minim_image(d12, crd2, crd1);
    pbc->diff_crd_minim_image(d12_prev, crd2_prev, crd1_prev);
    real_cst mass1_inv    = in_mass_inv[pair_atomids[i_cst][0]];
    real_cst mass2_inv    = in_mass_inv[pair_atomids[i_cst][1]];
    real_cst mass12_inv   = mass1_inv + mass2_inv;
    real_cst diff_sq      = d12[0] * d12[0] + d12[1] * d12[1] + d12[2] * d12[2];
    real_cst diff_sq_prev = d12_prev[0] * d12_prev[0] + d12_prev[1] * d12_prev[1] + d12_prev[2] * d12_prev[2];
    real_cst diff_ip      = d12[0] * d12_prev[0] + d12[1] * d12_prev[1] + d12[2] * d12_prev[2];
    
    real q1[3], q2[3];
    real s[3];
    for(int d=0; d<3; d++){
      q1[d] = in_vel[atomid1_3+d] + time_step * mass1_inv * 0.5 * -FORCE_VEL * in_work[atomid1_3+d];
      q2[d] = in_vel[atomid2_3+d] + time_step * mass2_inv * 0.5 * -FORCE_VEL * in_work[atomid2_3+d];
      s[d] = in_crd[atomid1_3+d] + time_step * q1[d] - in_crd[atomid2_3+d] - time_step * q2[d];
    }
    real s2 = s[0]*s[0] + s[1]*s[1] + s[2]*s[2];
    while(s2 >  pair_dist[i_cst]){
      

    }
   real_cst val_a = mass12_inv * mass12_inv * diff_sq_prev;
    real_cst val_b = mass12_inv * diff_ip;
    real_cst val_c = diff_sq - pair_dist[i_cst];
    real_cst val_g = (val_b - sqrt(max(val_b * val_b - val_c * val_a, (real_cst)0.0))) / val_a;
    
    for (int d = 0; d < 3; d++) {
      in_crd[atomid1_3 + d] += mass1_inv * val_g * d12_prev[d];
      in_crd[atomid2_3 + d] -= mass2_inv * val_g * d12_prev[d];
    }
  }
  return 0;
}

int ConstraintRattle::rattle_trio(real *in_crd, real *in_crd_prev, real *in_vel, real *in_work, real_pw *in_mass_inv, PBC *pbc){
    return 0;
}

int ConstraintRattle::rattle_quad(real *in_crd, real *in_crd_prev, real *in_vel, real *in_work, real_pw *in_mass_inv, PBC *pbc){
    return 0;
}

