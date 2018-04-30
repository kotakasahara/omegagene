#!/usr/bin/python2.7

from optparse import OptionParser
import sys
import os
import numpy as np

import kkmmsystem
import kkmmconfig 
import kkatomgroup

import kkpresto
import kkpresto_restart
import kkpresto_shake
import kkpresto_crd
import presto_generate_velocities as genvelo
import kkpdb 
import kkmm_vcmd
import kkceleste_ausrestart

MAX_RANDINT = 99999999

def get_options():
    p = OptionParser()
    p.add_option('--config', dest='fn_config',
                 # action="append",
                 help="Configuration file for omegagene")
    p.add_option('--crd-list', dest='fn_crd_list',
                 # action="append",
                 help="List of coordinate files")
    p.add_option('--o-vs-frm', dest='fn_o_vs_frm',
                 help="Output file name")
    p.add_option('--i-vs-frm', dest='fn_i_vs_frm',
                 help="Input file name")
    p.add_option("--init-type", dest="init_type",
                 type="choice",
                 choices = ["least", "random", "none"],
                 default="none",
                 help = "How to choose the initials. 'least' or 'random'.")
    p.add_option('--n-cal', dest='n_cal',
                 type="int",
                 help="the number of parallel runs")
    p.add_option('-t', '--temperature', dest='temperature',
                 type="float", default=300,
                 help="Temperature")
    p.add_option('--pref-cal', dest='pref_cal',
                 default="n",
                 help="prefix for the calculation directory")
    p.add_option('--pref-file', dest='pref_file',
                 default="md",
                 help="prefix for the calculation files")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts,args

class VcManager(object):
    def __init__(self, fn_config):
        self.fn_config = fn_config
        self.config = kkmmconfig.ConfigReader(self.fn_config).read_config()        
        self.prob = {}
        self.vs_frm = {}
        self.unmet_neighbor = False
        return
    def load_files(self):
        print "read_tpl"
        self.tpl = kkpresto.TPLReader(self.config.get_val("fn-i-tpl")).read_tpl()
        self.tpl.enumerate_12_13_14()
        self.tpl.convert_units_to_si()        
        print "read initial pdb"
        self.structure = kkpdb.PDBReader(self.config.get_val("fn-i-initial-pdb")).read_model()
        self.system = kkmmsystem.MmSystem(self.structure,
                                          self.config.get_val("cell-x"),
                                          self.config.get_val("cell-y"),
                                          self.config.get_val("cell-z"))
        self.system.set_atom_info_from_tpl(self.tpl)
        
        self.shkr = kkpresto_shake.SHKReader(self.config.get_val("fn-i-shake"))
        self.shkr.read_shk()
        self.shkr.expand_shake_info(self.tpl)

        self.read_extended_vcmd()        
        atom_groups_reader = kkatomgroup.AtomGroupsReader(self.config.get_val("fn-i-atom-groups"))
        print self.config.get_val("fn-i-atom-groups")
        print atom_groups_reader.fn
        self.atom_groups, self.atom_group_names = atom_groups_reader.read_groups()


        return
    def read_extended_vcmd(self):
        self.extended_vcmd = None
        if self.config.get_val("fn-i-vcmd-inp"):
            self.extended_vcmd = kkmm_vcmd.VcMDConf()
            self.extended_vcmd.read_params(self.config.get_val("fn-i-vcmd-inp"))
            # print "group names"
            # print self.extended_vcmd.group_names
        self.aus_restart = kkceleste_ausrestart.CelesteAUSRestart()
        self.aus_restart.set_aus_type(self.config.get_val("aus-type"))

        return
    def read_list(self, fn_list):
        f = open(fn_list)
        txt_list = []
        for line in f:
            tmp = line.strip().split()[0]
            txt_list.append(tmp)
        f.close()
        return txt_list
    def set_crd_files(self, fn_crd_list):
        self.fn_crd_list = self.read_list(fn_crd_list)
        return 
    def get_vs_candidates(self, lmb):
        vs_cand = []
        for dim_tmp, val_lmb in enumerate(lmb):
            dim = dim_tmp+1
            vs_cand_dim = []
            for vsid_tmp, rg in enumerate(self.extended_vcmd.lambda_ranges[dim][1:]):
                vsid = vsid_tmp + 1
                if val_lmb > rg[0] and val_lmb <= rg[1]:
                    vs_cand_dim.append(vsid)

            if len(vs_cand_dim) == 0:
                return []
            elif len(vs_cand_dim) == 1:
                vs_cand.append(vs_cand_dim[0])
            elif len(vs_cand_dim) == 2:
                val1 = self.extended_vcmd.lambda_ranges[dim][vs_cand_dim[0]][1] - val_lmb
                val2 = val_lmb - self.extended_vcmd.lambda_ranges[dim][vs_cand_dim[1]][0]
                if val1 > val2:
                    vs_cand.append(vs_cand_dim[0])
                else:
                    vs_cand.append(vs_cand_dim[1])
        return tuple(vs_cand)
    def cal_lambda(self, frame):
        if self.aus_restart.aus_type == self.aus_restart.AUS_TYPE["dist-mass-centers"]:
            lmb = self.cal_lambda_dist_mass_centers(frame)
        elif self.aus_restart.aus_type == self.aus_restart.AUS_TYPE["dist-min"]:
            lmb = self.cal_lambda_dist_min(frame)
        return lmb
    def cal_lambda_dist_mass_centers(self, frame):
        lmb = []
        def cal_center(grp_id):
            grp_mass = 0.0
            cent = np.zeros(3)
            for at1 in self.atom_groups[grp_id]:
                atid = at1-1
                mass = self.tpl.get_tplatom_from_atom_id(atid).mass
                grp_mass += mass
                cent += mass * frame.crds[atid]
            cent /= grp_mass
            return cent

        for grp_pair in self.extended_vcmd.group_names[1:]:
            grp_id1 = self.atom_group_names.index(grp_pair[0])
            grp_id2 = self.atom_group_names.index(grp_pair[1])
            cent1 = cal_center(grp_id1)
            cent2 = cal_center(grp_id2)
            dx = self.system.pbc.diff_crd_minim_image(cent1, cent2)
            dist = np.sqrt(np.sum(dx*dx))
            lmb.append(dist)
            
        return np.array(lmb)

    def cal_lambda_dist_min(self, frame):
        lmb = []
        for grp_pair in self.extended_vcmd.group_names[1:]:
            grp_id1 = self.atom_group_names.index(grp_pair[0])
            grp_id2 = self.atom_group_names.index(grp_pair[1])
            min_dist = 1e10
            for at1 in self.atom_groups[grp_id1]:
                atid1 = at1-1
                crd1 = frame.crds[atid1]
                for at2 in self.atom_groups[grp_id2]:
                    atid2 = at2-1
                    crd2 = frame.crds[atid2]
                    dx = self.system.pbc.diff_crd_minim_image(crd1, crd2)
                    dist = np.sqrt(np.sum(dx*dx))
                    if dist <= min_dist: min_dist = dist
            lmb.append(min_dist)
        return np.array(lmb)
    def count_vs_frm(self, processed_crd=set()):
        #self.vs_frm = {}
        ## vs_frm[(vs1,vs2,...)] = [(i_crd, i_frame), (i_crd, i_frame), ...]
        for i_crd, fn_crd in enumerate(self.fn_crd_list):
            if i_crd in processed_crd: continue
            print str(i_crd) + " " + fn_crd
            reader = kkpresto_crd.PrestoCrdReader(fn_crd)
            reader.read_information()
            reader.open()
            i_frame = 0
            frame = reader.read_next_frame()
            while frame:
                lmb = self.cal_lambda(frame)
                # print "lambda " + " ".join([str(x) for x in lmb])
                vs = self.get_vs_candidates(lmb)
                # print vs
                if vs:
                    if not vs in self.vs_frm:
                        self.vs_frm[vs] = []
                    self.vs_frm[vs].append((i_crd, i_frame))
                #---------------------------------
                frame = reader.read_next_frame()
                i_frame += 1
            reader.close()
        return
    def print_vs_frm(self, fn_out):
        fo = open(fn_out, "w")
        for vs, frms in self.vs_frm.items():
            line = ",".join([ str(x) for x in vs] )
            for crd, frame in frms:
                line += " " + str(crd)+":"+str(frame) 
            fo.write(line+"\n")
        fo.close()
        return
    def read_vs_frm(self, fn_vsfrm):
        f = open(fn_vsfrm)
        processed_crd = set()
        # self.vs_frm
        for line in f:
            terms = line.strip().split()
            vs = tuple([ int(x) for x in terms[0].split(",") ])
            print line
            print terms[0]
            print "aaa"
            for crdfrm in terms[1:]:
                print crdfrm
                tmp = crdfrm.split(":")
                crd = int(tmp[0])
                frm = int(tmp[1])
                if not vs in self.vs_frm:
                    self.vs_frm[vs] = []
                self.vs_frm[vs].append((crd, frm))
                processed_crd.add(crd)
        return processed_crd
    def pick_init_frames(self, init_type, n_cal):
        init_frames = []
        dict_frames = {}
        if init_type=="least":
            vs_frm_srt = sorted(self.vs_frm.items(), key=lambda x:len(x[1]))
            self.prob = {}
            prob_acc = 0

            for vs, frms in vs_frm_srt:
                prob_acc +=  1/float(len(frms))
                self.prob[vs] = prob_acc

            vs_frm_new = []
            prob_acc_new = 0
            if self.unmet_neighbor:
                for vs, frms in vs_frm_srt:
                    for i, it_vs in enumerate(vs):
                        new_vs = list(vs)
                        new_vs[i] -= 1
                        if new_vs[i] >= 1 and not tuple(new_vs) in self.prob:
                            prob_acc_new += 1
                            self.prob[tuple(new_vs)] = prob_acc_new
                            vs_frm_new.append((tuple(new_vs),frms))
                            print "add " + ",".join([str(x) for x in new_vs])
                        new_vs[i] += 2
                        if new_vs[i] <= self.extended_vcmd.n_vs[i] and not tuple(new_vs) in self.prob:
                            prob_acc_new += 1
                            self.prob[tuple(new_vs)] = prob_acc_new
                            vs_frm_new.append((tuple(new_vs),frms))
                            print "add " + ",".join([str(x) for x in new_vs])

            for vs, frms in vs_frm_srt:
                self.prob[vs] += prob_acc_new

            for vs in self.prob.keys():
                self.prob[vs] /= ( prob_acc + prob_acc_new)
                
            def add_frame(vs, frms):
                return

            vs_frm_new.extend(vs_frm_srt)

            print "n cand"
            print len(vs_frm_new)
            
            for vs, frms in vs_frm_new:
                print ",".join([ str(x) for x in vs ]) + " " + str(self.prob[vs])

            cnt = 0
            while len(init_frames) < n_cal and cnt <= 500:
                rnd = np.random.rand()
                flg = False
                for vs, frms in vs_frm_new:
                    if rnd <= self.prob[vs]:
                        tmp = frms[np.random.randint(len(frms))]
                        buf = (tmp[0], tmp[1], vs)
                        ## buf [ file-id in self.fn_crd_list, 
                        ##       frame-id,
                        ##       vs ]
                        init_frames.append(buf)
                        if not tuple(vs) in dict_frames:
                            dict_frames[tuple(vs)] = 0
                        dict_frames[tuple(vs)] += 1
                        break
                cnt+=1
                #if not flg:
                #print "Sampling is failed."

        elif init_type=="random":
            vs_frm_srt = sorted(self.vs_frm.items(), key=lambda x:len(x[1]))
            for i in range(n_cal):
                rnd = np.random.rand(len(vs_frm_srt))
                frms = vs_frm_srt[rnd]
                tmp = frms[np.random.randint(len(frms))]
                buf = (tmp[0], tmp[1], vs)
                init_frames.append(buf)
                if not tuple(vs) in dict_frames:
                    dict_frames[tuple(vs)] = 0
                dict_frames[tuple(vs)] += 1
                break

        if init_type!="none":
            pass

        print "picked frames:"
        for vs, crdfrm in dict_frames.items():
            print ",".join([str(x) for x in vs]) + " : " + str(crdfrm)
        return init_frames
    def prepare_initials(self, init_type, pref_cal, pref_file, n_cal, temperature):
        init_frames = self.pick_init_frames(init_type, n_cal)
        for i, crdfrmvs in enumerate(init_frames):
            cal_dir = pref_cal + str(i+1)
            try:    os.mkdir(cal_dir)
            except: pass
            #self.gen_pdb_frm(fn_pdb, crdfrmvs[0], crdfrmvs[1])
            fn_restart = os.path.join(cal_dir, pref_file+"_0.restart")
            print self.fn_crd_list[crdfrmvs[0]]+" "+str(crdfrmvs[1])+" " + ",".join([str(x) for x in crdfrmvs[2]])
            self.gen_restart_frm(fn_restart, crdfrmvs[0], crdfrmvs[1],
                                 temperature)
            fn_startvirt = os.path.join(cal_dir, "start.virt")
            
            self.gen_startvirt_frm(fn_startvirt, crdfrmvs[2])
        return
    #def gen_pdb_frm(self, fn, crd, frm):
    #return 
    def gen_restart_frm(self, fn, i_codfile, i_frm, temperature):
        rest = kkpresto_restart.PrestoRestart()        
        rest.n_atoms = len(self.structure.atoms)
        reader = kkpresto_crd.PrestoCrdReader(self.fn_crd_list[i_codfile])
        # print self.fn_crd_list[i_codfile] + " - "  + str(i_frm)
        
        reader.read_information()
        reader.open()
        frame = reader.read_next_frame()
        for i_frame in range(1, i_frm):
            #print i_frame
            frame = reader.read_next_frame()
        reader.close()
        rest.crd = frame.crds
        #self.tpl.convert_units_to_si()
        rest.n_vel = rest.n_atoms
        rest.vel = genvelo.generate_velocities(self.structure,
                                               self.tpl, self.shkr,
                                               temperature, True,
                                               False)
        rest.vel = genvelo.convert_velocity_si_to_ang_fs(rest.vel)
        kkpresto_restart.PrestoRestartWriter(fn).write_restart(rest)
        return 
    def gen_startvirt_frm(self, fn, vs):
        fo = open(fn, "w")
        fo.write(str(len(vs))+"\n")
        fo.write("".join([ str(x) + "\n" for x in vs ]))
        fo.write(str(np.random.randint(MAX_RANDINT))+"\n")
        fo.close()
        return 

def _main():
    opts, args = get_options()
    manager = VcManager(opts.fn_config)
    manager.load_files()
    manager.set_crd_files(opts.fn_crd_list)
    processed_crd = set()
    if opts.fn_i_vs_frm:
        processed_crd = manager.read_vs_frm(opts.fn_i_vs_frm)
    manager.count_vs_frm(processed_crd)
    manager.print_vs_frm(opts.fn_o_vs_frm)
    manager.prepare_initials(opts.init_type, opts.pref_cal,
                             opts.pref_file, opts.n_cal,
                             opts.temperature)
    return

if __name__ == "__main__":
    _main()

