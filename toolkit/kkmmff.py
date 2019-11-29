import kkpbc
import numpy as np
import scipy
import scipy.special
import kkmmconst as cst

class FF(object):
    DEBUG = True
    def __init__(self):
        self.cutoff = 12.0
        self.ele_cutoff = 10.0
        self.ewaldalpha = 0.0
        self.alphapi = 0.0
        self.func_calc_zd = None
        self.func_calc_excess = None
        self.bcoeff = 0.0
        self.fcoeff = 0.0
        self.zcore = 0.0   
        self.ewaldalpha = 0
        return
    def set_params(self, cutoff):
        self.cutoff = cutoff
        self.vdw_cutoff = cutoff
        self.ele_cutoff = cutoff
        return 
    def set_zd_params(self, ewaldalpha,
                      charges,
                      atom_id_12, atom_id_13, atom_id_14nb):
        self.bcoeff = 0.0
        self.fcoeff = 0.0
        self.zcore = 0.0
        self.ewaldalpha = ewaldalpha
        self.piewald = 0.0

        cutoff_2 = self.ele_cutoff ** 2
        cutoff_3 = self.ele_cutoff * cutoff_2
        tmp = self.ele_cutoff * ewaldalpha
        errorfc = scipy.special.erfc(tmp)
        expterm = np.exp(-tmp*tmp)

        if self.ewaldalpha < cst.Const().eps:
            self.func_calc_zd = self.calc_zd_alpha0
            self.func_calc_zd_excess = self.calc_zd_excess_alpha0
            self.fcoeff = 1.0 / cutoff_3
            self.bcoeff = 0.5 * self.fcoeff
            self.scoeff = 1.5 / cutoff_3
            self.zcore = 1.5 / self.ele_cutoff
        else:
            self.piewald = 2.0 / (np.sqrt(np.pi) * self.ewaldalpha)
            self.func_calc_zd = self.calc_zd_alpha
            self.func_calc_zd_excess = self.calc_zd_excess_alpha
            self.fcoeff = errorfc / cutoff_3 + piewald * expterm / cutoff_2
            self.bcoeff = 0.5 * fcoeff
            self.scoeff = bcoeff + errorfc / cutoff_3
            self.zcore = errorfc / self.ele_cutoff + bcoeff * cutoff_2
        ##self energy
        ##print "self energy"
        dwork = np.zeros(len(charges))
        for pair, params in atom_id_12:
            dwork[pair[0]] += charges[pair[1]]
            dwork[pair[1]] += charges[pair[0]]
        for trio, params in atom_id_13:
            dwork[trio[0]] += charges[trio[2]]
            dwork[trio[2]] += charges[trio[0]]
        for quad, params in atom_id_14nb:
            #if params[4] != 0: continue
            dwork[quad[0]] += charges[quad[1]]
            dwork[quad[1]] += charges[quad[0]]


        chgsum = 0.0
        chgsmm = 0.0
        self.energy_self = np.zeros(len(charges))
        for i, charge_i in enumerate(charges):
            self.energy_self[i] = (-self.piewald * charge_i ** 2 - self.zcore * charge_i * (charge_i+dwork[i])) * 0.5 * cst.Const().charge_coeff
            chgsum += charge_i * charge_i
            chgsmm += charge_i * (charge_i + dwork[i])
            #print "self: " + str(charge_i) + str(dwork[i])
        self.d_self = (-self.zcore * chgsmm - self.piewald * chgsum) * 0.5 * cst.Const().charge_coeff

        self.d_self_mon = np.sum(self.energy_self)

        if FF.DEBUG:
            print("fcoeff : " + str(self.fcoeff))
            print("bcoeff : " + str(self.bcoeff))
            print("scoeff : " + str(self.scoeff))
            print("zcore : " + str(self.zcore))
            print("chgsum : " + str(chgsum))
            print("chgsmm : " + str(chgsmm))
            print("selfEnergy_dcore : " + str(self.zcore*chgsmm*0.5*cst.Const().charge_coeff))
            print("d_self : " + str(self.d_self))
            print("sum of energy_self : " + str(self.d_self_mon))
        return 

    def calc_zd_alpha0(self, r12, r12sq, r12inv, r12inv_2, r12inv_3, cc):
        ene_ele = cc * (r12inv - self.zcore + self.bcoeff * r12sq)
        grad_coef = - cc * (r12inv_3 - self.fcoeff)
        return ene_ele, grad_coef
    def calc_zd_alpha(self, r12, r12sq, r12inv, r12inv_2, r12inv_3, cc):
        tmp = r12 * alphaewald
        errorfc = scipy.special.erfc(tmp)
        ene_ele = cc * (r12inv * errorfc - self.zcore + bcoeff * r12sq)
        grad_coef = - cc * (errorfc * r12inv_3 + self.piewald * np.exp(tmp*tmp) * r12inv_2 - self.fcoeff)
        return ene_ele, grad_coef

    def calc_zd_excess_alpha0(self, r12, r12sq, r12inv_e, cc):
        ene = cc * self.bcoeff * r12sq
        #print "Exess_alpha0: " + str(cc) + " " + str(self.bcoeff) + " " + str(r12sq)
        grad_coef = cc * self.fcoeff
        return ene, grad_coef
    def calc_zd_excess_alpha(self, r12, r12sq, r12inv_3, cc):
        tmp = self.alphaewald * r12
        errfunc = 1.0 - scipy.special.erfc(tmp)
        ene = cc * ( -errfunc * r121 + self.bcoeff * r12sq )
        grad_coef = -cc * ( errfunc * r12inv_3 + self.piewald * np.exp(-tmp*tmp) * r12 - self.fcoeff - r12inv_3)
        return ene, grad_coef

    def calc_bond(self, crd1, crd2, param_e, param_r0, pbc):
        d = pbc.diff_crd_minim_image(crd1, crd2)
        r12 = np.sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2])
        dif = r12 - param_r0
        coef = param_e * dif
        ene = coef * dif
        #gradient
        grad_coef = 2.0 * coef / r12
        work = d * grad_coef
        #if FF.DEBUG:
        #    print "d: "+str(d[0])+" "+str(d[1])+" "+str(d[2])
        #    print "r12: "+str(r12)
        #    print "dif: " +str(dif) + " coef:"+str(coef)
        #    print "ene: " + str(ene)
        #    print "work: " + str(work[0]) + " " + str(work[1]) + " " + str(work[2])
            
        return ene, work
    
    def calc_angle(self, crd1, crd2, crd3, param1, param2, pbc):
    ### eps is a very small value
    ### to avoid angle == 0.0    d12 = pbc.diff_crd_minim_image(crd1, crd2)
        d12 = pbc.diff_crd_minim_image(crd1, crd2)
        d32 = pbc.diff_crd_minim_image(crd3, crd2)
        r12sq = np.sum(d12*d12)
        r32sq = np.sum(d32*d32)
        r12r32 = np.sqrt(r12sq * r32sq)
        r12r32_inv = 1.0 / r12r32
        cost = np.sum(d12 * d32) * r12r32_inv
        
        cost = np.max([-1.0, np.min([1.0, cost])])
        theta = np.arccos(cost)
        
        dif = theta - param2
        coef = param1 * dif
        ene = coef * dif
        
        r12sq_inv = 1.0 / r12sq
        r32sq_inv = 1.0 / r32sq
        
        sin_theta = np.sin(theta)
        sin_theta = max(cst.Const().eps, sin_theta)
        
        coef_grad = - (2.0 * coef) / sin_theta
        work1 = coef_grad * ( d32 * r12r32_inv - d12 * cost * r12sq_inv )
        work2 = coef_grad * ( d12 * r12r32_inv - d32 * cost * r32sq_inv )
        return ene, work1, work2

# torsion energy
    def calc_torsion(self, crd1, crd2, crd3, crd4,
                     param1, param2, param3, param4, param5,
                     pbc):
        ## param1 = FFTOR
        ## param3 = FROT
        ## param4 = FPHS
        if param1 <= cst.Const().eps:
            return [0.0, np.zeros(3), np.zeros(3), np.zeros(3)]
        
        d21 = pbc.diff_crd_minim_image(crd2, crd1)
        d32 = pbc.diff_crd_minim_image(crd3, crd2)
        d43 = pbc.diff_crd_minim_image(crd4, crd3)
    #### VALIDATE THIS ROUTINE LATOR ! 
        p12 = np.cross(d21,d32)
        #print "p12:" + str(p12)
    ## p_x = a_y * b_z - b_y * a_z
    ## p_y = a_z * b_x - b_z * a_x
    ## p_z = a_x * b_y - b_x * a_y
    
        p23 = np.cross(d32,d43)
        #print "p23:" + str(p23)        
        r12sq = np.sum(p12 * p12)
        r23sq = np.sum(p23 * p23)
        #print "r12sq " + str(r12sq)
        #print "r23sq " + str(r23sq)
        r12r23 = np.sqrt(r12sq * r23sq)
    ####
        #print "r12r23 " + str(r12r23)
        if r12r23 <= 0.0:
            return 0.0, np.array([0.0, 0.0, 0.0])
        r12r23_inv = 1.0 / r12r23
        s1223 = np.sum(p12*p23)
        cosp = s1223 * r12r23_inv
        cosp = max(-1.0, min(1.0, cosp))
        p123 = np.cross(p12, p23)
        s32123 = np.sum(d32*p123)
        phi = np.arccos(cosp)
        if s32123 < 0.0:
            phi = 2.0 * np.pi - phi
        phimod = param3 * phi - param4
        div = 1.0 / float(param2)
        ene = param1 * (1.0 + np.cos(phimod)) * div

    ##
        sinp = np.sin(phi)
        nrot = int(param3)

        sinp_sq = sinp*sinp
        cosp_sq = cosp*cosp

        sin_fphs = np.sin(param4)
        cos_fphs = np.cos(param4)
        sin_rot = np.sin(param3 * phi)
        cos_rot = np.cos(param3 * phi)

    ## lim = limit ?? yoku wakannai
        dums = sinp
        gmul = [0.0, 0.0, 2.0, 0.0, 4.0,  0.0,
                6.0, 0.0, 8.0, 0.0, 10.0]
        dflim = cos_fphs * (param3 - gmul[nrot] + gmul[nrot] * cosp)
        #print "cos_fphs: " + str(cos_fphs)
        #print "param3: " + str(param3)
        #print "cosp: " + str(cosp)
        #print "nrot: " + str(nrot)
        #print "gmul[nrot]: " + str(gmul[nrot])
        df0 = cos_fphs * sin_rot - sin_fphs * cos_rot
        #print "dflim: " + str(dflim)
        #print "df0: " + str(df0)

        df1 = dflim
        if cst.Const().eps <= np.fabs(dums):
            #df0 = cos_fphs = np.cos(param4)
            df1 = df0 / dums
    ## kokorahen maji de wake wakaran
    ## ato de zikkuri kangaeyo 
    ## ima ha tsukareta 
        #print "df1: " + str(df1)
        work_coeff = param3 * param1 * div * df1 * r12r23_inv
        #print "frot: " + str(param3)
        #print "fftor: " + str(param1)
        #print "div: " + str(div)
        #print "r12r23: " + str(r12r23_inv)
        #print "work_coeff: " + str(work_coeff)
        r12 = 1.0 / r12sq * s1223
        r23 = 1.0 / r23sq * s1223

        op_p23d32 = np.cross(p23, d32)
        op_p12d32 = np.cross(p12, d32)

        work1 = work_coeff * (op_p23d32 - r12 * op_p12d32)

        d2132 = d21 + d32
        op_d2132_p23 = np.cross(d2132, p23)
        op_d2132_p12 = np.cross(d2132, p12)
        op_p12_d43 = np.cross(p12, d43)
        #op_p12_d2132 = np.cross(p12, d2132)
        op_p23_d43 = np.cross(p23, d43)
        work2 = work_coeff * (op_d2132_p23 + op_p12_d43 \
            - r12 * op_d2132_p12 - r23 * op_p23_d43)

        op_p12_d32 = np.cross(p12, d32)
        op_p23_d32 = np.cross(p23, d32)
        work3 = work_coeff * (op_p12_d32 - r23 * (op_p23_d32))
        return ene, work1, work2, work3

    def calc_14pair(self, crd1, crd4,
                    param_nb1, param_nb2, param_141, param_142,
                    charge1, charge4,
                    pbc):
        d = pbc.diff_crd_minim_image(crd1, crd4)
        r14sq = d[0]*d[0] + d[1]*d[1] + d[2]*d[2]
        r14sq_inv = 1.0 / r14sq
        r14inv = np.sqrt(r14sq_inv)
        
        ## vdw
        r14inv_6 = r14sq_inv * r14sq_inv * r14sq_inv 
        r14inv_12 = r14inv_6 * r14inv_6
        co6 = param_nb1 # fynbpp
        co12 = param_nb2 # fynbpp
        term6 = co6 * r14inv_6
        term12 = co12 * r14inv_12
        #print "co6: " + str(co6)
        #print "r14inv_6: " + str(r14inv_6)
        #print "co12: " + str(co12)
        #print "r14inv_12: " + str(r14inv_12)
        #print "param_142: " + str(param_142)
        ene_vdw = (-term6 + term12) * param_142 # fytvws
        #print "ene_vdw: " + str(ene_vdw)        
        ## vdw grad
        coef_work = r14sq_inv * param_142 * ( -12.0 * term12 + 6.0 * term6 )
        work_vdw = coef_work * d
        #print "work_vdw"
        #print work_vdw
        ## ele
        cc = charge1 * charge4 * cst.Const().charge_coeff * param_141 # fytess
        ene_ele = cc * r14inv

        ## ele grad
        coef_ele = -cc * r14sq_inv * r14inv
        work_ele = coef_ele * d
        #print "work_ele"
        #print work_ele

        return ene_vdw, work_vdw, ene_ele, work_ele
    
    def calc_pairwise(self, crd1, crd2,
                      param1, param2,
                      charge1, charge2, pbc):
        d = pbc.diff_crd_minim_image(crd1, crd2)
        r12sq = d[0]*d[0] + d[1]*d[1] + d[2]*d[2]
        r12 = np.sqrt(r12sq)

        ene_vdw = 0.0
        work_vdw = np.array([0.0, 0.0, 0.0])
        ene_ele = 0.0
        work_ele = np.array([0.0, 0.0, 0.0])
        
        if r12 >= self.cutoff:
            return ene_vdw, work_vdw, ene_ele, work_ele

        ## vdw
        r12inv = 1.0 / r12
        r12inv_2 = 1.0 / r12sq
        r12inv_3 = r12inv * r12inv_2
        r12inv_6 = r12inv_3 * r12inv_3
        r12inv_12 = r12inv_6 * r12inv_6
        
        #print "param1"
        #print param1
        #print "param2"
        #print param2


        term6 = param1 * r12inv_6
        term12 = param2 * r12inv_12
        ene_vdw = -term6 + term12
        ## grad
        work_coef = r12inv_2 * (-12.0*term12 + 6.0*term6)
        work_vdw = work_coef * d
        #print r12inv_12
        ### 
        cc = charge1 * charge2 * cst.Const().charge_coeff
        coef = 0.0
        ene_ele, coef = self.func_calc_zd(r12, r12sq, r12inv, r12inv_2, r12inv_3, cc)
        
        work_ele = coef * d
        return ene_vdw, work_vdw, ene_ele, work_ele
    
    def calc_zd_excess(self, crd1, crd2,
                       charge1, charge2, pbc):
        ene = 0.0
        work = np.array([0.0, 0.0, 0.0])
        d = pbc.diff_crd_minim_image(crd1, crd2)
        r12sq = d[0]*d[0] + d[1]*d[1] + d[2]*d[2]
        r12 = np.sqrt(r12sq)
        r12inv = 1.0 / r12
        r12inv_3 = r12inv ** 3

        cc = charge1 * charge2 * cst.Const().charge_coeff
        ene, grad_coef = self.func_calc_zd_excess(r12, r12sq, r12inv_3, cc)
        work = d * grad_coef
        return  ene, work
    
