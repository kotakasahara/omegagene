#ifndef __ZERO_MULTIPOLE_SUM_H__
#define __ZERO_MULTIPOLE_SUM_H__

#include "ElectrostaticObject.h"
#include <cmath>
#include <iostream>

class ZeroMultipoleSum : public ElectrostaticObject {
  private:
  protected:
    const int zms_mode;
    // zms_mode ...  1: dipole
    //               2: quadrupole
    //               3: octapole
    //               4: hexadeca

    const real_pw cutoff;
    const real_pw ewald_alpha;
    const real_pw dielectric_inv;      // for Debye-Huckel
    const real_pw ionic_strength;  // for Debye-Huckel
    const real_pw temperature;
    //const real_pw dh_const; // AVOGADRO/(1000*JOULE_CAL*4*pi*permittivity*dielectric)
    const real_pw debye_length_inv;

    real_pw       bcoeff;
    real_pw       fcoeff;
    real_pw       scoeff;

    real_pw hzqrc;
    real_pw zqcoeff2;
    real_pw zqcoeff4;
    real_pw zqcoeff22;
    real_pw zqcoeff44;

    real_pw zocoeff2;
    real_pw zocoeff4;
    real_pw zocoeff6;
    real_pw zocoeff22;
    real_pw zocoeff44;
    real_pw zocoeff66;

    real_pw zhcoeff2;
    real_pw zhcoeff4;
    real_pw zhcoeff6;
    real_pw zhcoeff8;
    real_pw zhcoeff22;
    real_pw zhcoeff44;
    real_pw zhcoeff66;
    real_pw zhcoeff88;

    real_pw coeffd1;
    real_pw coeffd2;
    real_pw coeffd3;
    real_pw coeffd4;

    real_pw zcore;
    real_pw zqcore;
    real_pw zocore;
    real_pw zhcore;

    real_pw piewald;
    real_pw cutoff_2;
    real_pw cutoff_3;
    real_pw cutoff_4;
    real_pw cutoff_5;
    real_pw cutoff_6;
    real_pw cutoff_7;

    real_pw d_self;
    real_pw d_self_mon;

  public:
    typedef int (ZeroMultipoleSum::*ZmsCalc)(real_pw &,
                                             real_pw &,
                                             const real_pw &,
                                             const real_pw &,
                                             const real_pw &,
                                             const real_pw &,
                                             const real_pw &,
                                             const real_pw &);
    typedef int (ZeroMultipoleSum::*ZmsCalcExcess)(real_pw &,
                                                   real_pw &,
                                                   const real_pw &,
                                                   const real_pw &,
                                                   const real_pw &,
                                                   const real_pw &,
                                                   const real_pw &,
                                                   const real_pw &);

    ZeroMultipoleSum(const int           in_zms_mode,
                     const real          in_alpha,
                     const real_pw       in_cutoff,
                     const ZmsCalc       in_func,
                     const ZmsCalcExcess in_func_excess,
		     const real          in_dielectric,
		     const real          in_ionic_strength,
		     const real          in_temperature);
    // int (* const func)(real_pw&, real_pw&, const real_pw&, const real_pw&, const real_pw&, const real_pw&, const
    // real_pw&, const real_pw&),
    // int (* const func_excess)(real&, real&, const real&, const real&, const real&, const real&, const real&, const
    // real&)
    //);
    const ZmsCalc       func_calc_zms;
    const ZmsCalcExcess func_calc_zms_excess;

    virtual int set_config_parameters(const Config *in_cfg);
    virtual int initial_preprocess();
    int         set_zms_params();
    int cal_self_energy(const int &  n_atoms,
                        const int &  n_excess,
                        const int **&excess_pairs,
                        /*
                        const int& n_bonds,
                        const int**& bond_atomid_pairs,
                        const int& n_angles,
                        const int**& angle_atomid_triads,
                        const int& n_torsions,
                        const int**& torsion_atomid_quads,
                        const int*& torsion_nb14,*/
                        real_pw *&charge,
                        real *&   energy_self,
                        real &    energy_self_sum);
    inline real_pw get_bcoeff() { return bcoeff; }
    inline real_pw get_fcoeff() { return fcoeff; }
    inline real_pw get_scoeff() { return scoeff; }
    inline real_pw get_zcore() { return zcore; }
    int calc_zero02pole_excess_alpha0(real_pw &      ene_ele,
                                      real_pw &      grad_coeff,
                                      const real_pw &r12,
                                      const real_pw &r12_2,
                                      const real_pw &r12_inv,
                                      const real_pw &r12_2_inv,
                                      const real_pw &r12_3_inv,
                                      const real_pw &cc);
    int calc_zero02pole_excess(real_pw &      ene_ele,
                               real_pw &      grad_coeff,
                               const real_pw &r12,
                               const real_pw &r12_2,
                               const real_pw &r12_inv,
                               const real_pw &r12_2_inv,
                               const real_pw &r12_3_inv,
                               const real_pw &cc);
    int calc_zero04pole_excess_alpha0(real_pw &      ene_ele,
                                      real_pw &      grad_coeff,
                                      const real_pw &r12,
                                      const real_pw &r12_2,
                                      const real_pw &r12_inv,
                                      const real_pw &r12_2_inv,
                                      const real_pw &r12_3_inv,
                                      const real_pw &cc);
    int calc_zero04pole_excess(real_pw &      ene_ele,
                               real_pw &      grad_coeff,
                               const real_pw &r12,
                               const real_pw &r12_2,
                               const real_pw &r12_inv,
                               const real_pw &r12_2_inv,
                               const real_pw &r12_3_inv,
                               const real_pw &cc);
    int calc_zero08pole_excess_alpha0(real_pw &      ene_ele,
                                      real_pw &      grad_coeff,
                                      const real_pw &r12,
                                      const real_pw &r12_2,
                                      const real_pw &r12_inv,
                                      const real_pw &r12_2_inv,
                                      const real_pw &r12_3_inv,
                                      const real_pw &cc);
    int calc_zero08pole_excess(real_pw &      ene_ele,
                               real_pw &      grad_coeff,
                               const real_pw &r12,
                               const real_pw &r12_2,
                               const real_pw &r12_inv,
                               const real_pw &r12_2_inv,
                               const real_pw &r12_3_inv,
                               const real_pw &cc);
    int calc_zero16pole_excess_alpha0(real_pw &      ene_ele,
                                      real_pw &      grad_coeff,
                                      const real_pw &r12,
                                      const real_pw &r12_2,
                                      const real_pw &r12_inv,
                                      const real_pw &r12_2_inv,
                                      const real_pw &r12_3_inv,
                                      const real_pw &cc);
    int calc_zero16pole_excess(real_pw &      ene_ele,
                               real_pw &      grad_coeff,
                               const real_pw &r12,
                               const real_pw &r12_2,
                               const real_pw &r12_inv,
                               const real_pw &r12_2_inv,
                               const real_pw &r12_3_inv,
                               const real_pw &cc);
    int calc_zero02pole_alpha0(real_pw &      ene_ele,
                               real_pw &      grad_coeff,
                               const real_pw &r12,
                               const real_pw &r12_2,
                               const real_pw &r12_inv,
                               const real_pw &r12_2_inv,
                               const real_pw &r12_3_inv,
                               const real_pw &cc);
    int calc_zero02pole(real_pw &      ene_ele,
                        real_pw &      grad_coeff,
                        const real_pw &r12,
                        const real_pw &r12_2,
                        const real_pw &r12_inv,
                        const real_pw &r12_2_inv,
                        const real_pw &r12_3_inv,
                        const real_pw &cc);
    int calc_zero04pole_alpha0(real_pw &      ene_ele,
                               real_pw &      grad_coeff,
                               const real_pw &r12,
                               const real_pw &r12_2,
                               const real_pw &r12_inv,
                               const real_pw &r12_2_inv,
                               const real_pw &r12_3_inv,
                               const real_pw &cc);
    int calc_zero04pole(real_pw &      ene_ele,
                        real_pw &      grad_coeff,
                        const real_pw &r12,
                        const real_pw &r12_2,
                        const real_pw &r12_inv,
                        const real_pw &r12_2_inv,
                        const real_pw &r12_3_inv,
                        const real_pw &cc);
    int calc_zero08pole_alpha0(real_pw &      ene_ele,
                               real_pw &      grad_coeff,
                               const real_pw &r12,
                               const real_pw &r12_2,
                               const real_pw &r12_inv,
                               const real_pw &r12_2_inv,
                               const real_pw &r12_3_inv,
                               const real_pw &cc);
    int calc_zero08pole(real_pw &      ene_ele,
                        real_pw &      grad_coeff,
                        const real_pw &r12,
                        const real_pw &r12_2,
                        const real_pw &r12_inv,
                        const real_pw &r12_2_inv,
                        const real_pw &r12_3_inv,
                        const real_pw &cc);
    int calc_zero16pole_alpha0(real_pw &      ene_ele,
                               real_pw &      grad_coeff,
                               const real_pw &r12,
                               const real_pw &r12_2,
                               const real_pw &r12_inv,
                               const real_pw &r12_2_inv,
                               const real_pw &r12_3_inv,
                               const real_pw &cc);
    int calc_zero16pole(real_pw &      ene_ele,
                        real_pw &      grad_coeff,
                        const real_pw &r12,
                        const real_pw &r12_2,
                        const real_pw &r12_inv,
                        const real_pw &r12_2_inv,
                        const real_pw &r12_3_inv,
                        const real_pw &cc);
    int calc_debye_huckel(real_pw &      ene_ele,
			  real_pw &      grad_coeff,
			  const real_pw &r12,
			  const real_pw &r12_2,
			  const real_pw &r12_inv,
			  const real_pw &r12_2_inv,
			  const real_pw &r12_3_inv,
			  const real_pw &cc);
    int calc_null(real_pw &      ene_ele,
		  real_pw &      grad_coeff,
		  const real_pw &r12,
		  const real_pw &r12_2,
		  const real_pw &r12_inv,
		  const real_pw &r12_2_inv,
		  const real_pw &r12_3_inv,
		  const real_pw &cc);
};

#endif
