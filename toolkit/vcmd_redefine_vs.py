#!/usr/bin/python2.6

from optparse import OptionParser
import sys
import os
import numpy as np
import re
import kkmm_vcmd
import copy
def opt_parse():
    p = OptionParser()
    p.add_option('--i-new-vs', dest='fn_new_vs',
                 help="q_cano file")
    p.add_option('-o', dest='fn_out',
                 help="filename for output")    
    p.add_option('--o-cano', dest='fn_o_canonical',
                 help="filename for output")    
    p.add_option('-b','--bin-width', dest='bin_width',
                 type="float",
                 help="bin width for lambda")
    p.add_option('--stages', dest='stages',
                 help="numbers separated with 'i,j'. The range can be specified by 'i-j'")
    p.add_option('--series', dest='series',
                 help="numbers separated with 'i,j'. The range can be specified by 'i-j'")
    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

class VcMDData(object):
    def __init__(self, bin_width):
        self.bin_width = bin_width
        # self.distrib[(vs1,vs2...)] = prob
        self.distrib = {}
        self.sum_prob = 0.0
        self.dim = 0

        self.path_cal = "."
        self.fn_lambda = "lambda.out"
        self.fn_vslog = "ttp_vcmd.out"
        self.fn_qcano = "vcmd_next.inp"
        self.stages = []
        self.series = []

        ## self.files_lambda[stage][series] = "fn"
        self.files_lambda = {} 
        self.files_vslog = {} 
        self.files_qcano = {}


        return
    def parse_unit(self, term):
        at = re.compile("^(\d+)$")
        at_range = re.compile("^(\d+)\-(\d+)$")
        m1 = at.match(term)
        m2 = at_range.match(term)
        val = []
        if m1:
            val = [int(m1.group())]
        elif m2:
            val = [x for x in range(int(m2.group(1)), int(m2.group(2))+1)]
        return val

    def read_fnlist(self, fn):
        f = open(fn)
        fnlist = []
        for line in f:
            tmp = line.strip().split()[0]
            fnlist.append(tmp)
        f.close()
        return fnlist

    def set_trajectory_files(self, stages_str, series_str):
        self.stages = self.parse_unit(stages_str)
        self.series = self.parse_unit(series_str)
        for st in self.stages:
            fn_qcano = os.path.join(self.path_cal, str(st-1), self.fn_qcano)
            if not os.path.exists(fn_qcano):
                sys.stderr.write("File not found; qcano for the stage "+str(st)+"\n")
                sys.stderr.write(fn_qcano)
                continue
            self.files_qcano[st] = fn_qcano
            self.files_lambda[st] = {}
            self.files_vslog[st] = {}
            for se in self.series:
                dir_cal = os.path.join(self.path_cal, str(st), "n"+str(se))
                fn_lambda = os.path.join(dir_cal, self.fn_lambda)
                fn_vslog = os.path.join(dir_cal, self.fn_vslog)
                if not os.path.exists(fn_lambda):
                    sys.stderr.write("File not found; \n")
                    sys.stderr.write(fn_lambda+"\n")
                    continue
                if not os.path.exists(fn_vslog):
                    sys.stderr.write("File not found; \n")
                    sys.stderr.write(fn_vslog+"\n")
                    continue
                self.files_lambda[st][se] = fn_lambda
                self.files_vslog[st][se] = fn_vslog
        #print self.files_lambda
        return
        
    def read_trj(self, fn_trj):
        trj = []
        f = open(fn_trj)
        for line in f:
            term = [float(x) for x in line.strip().split()]
            trj.append(term)
        f.close()
        return trj
    def add_distrib(self, vcconf, trj_lambda, trj_vslog):
        for i, vs_l in enumerate(trj_vslog):
            vs = tuple(vs_l)
            
            bin_id = tuple([int(x/self.bin_width) for x in trj_lambda[i]])

            if not bin_id in self.distrib:
                self.distrib[bin_id] = 0.0
                
            self.distrib[bin_id] += vcconf.params[vs][0]
            self.sum_prob +=  vcconf.params[vs][0]
            
        return
    def normalize_distrib(self):
        for bin_id, val in self.distrib.items():
            self.distrib[bin_id] = val/self.sum_prob
        return
    def calc_canonical(self):
        self.distrib = {}
        self.sum_prob = 0.0
        for st, fn_qcano in self.files_qcano.items():
            vcconf = kkmm_vcmd.VcMDConf()
            vcconf.read_params(fn_qcano)
            if self.dim == 0:
                self.dim = vcconf.dim
            if self.dim != vcconf.dim:
                sys.stderr.write("Inconsistency in the VS dimension\n")
                sys.stderr.write(fn_qcano+"\n")
                sys.exit(1)

            for se, fn_lambda in self.files_lambda[st].items():
                print "Stage:"+str(st) + " - Series:" + str(se)
                fn_vslog = self.files_vslog[st][se]
                trj_lambda = self.read_trj(fn_lambda)
                trj_vslog = self.read_trj(fn_vslog)
                self.add_distrib(vcconf, trj_lambda, trj_vslog)
        self.normalize_distrib()
        return
    def write_canonical(self, fn_cano):
        f = open(fn_cano, "w")
        for bin_id, val in self.distrib.items():
            lmb = [ str((x+0.5)*self.bin_width) for x in bin_id ]
            line = " ".join(lmb)+" "+str(val)+"\n"
            f.write(line)
        f.close()
        return

    def enum_vs(self, buf_vs, cur_vs, all_vs):
        if len(buf_vs) < 1: return all_vs
        for i, b0 in enumerate(buf_vs[0]):
            tmp_buf = copy.deepcopy(buf_vs[1:])
            tmp_cur = copy.deepcopy(cur_vs)
            tmp_cur.append(b0)

            if len(tmp_cur) == self.dim:
                all_vs.append(tuple(tmp_cur))
            else:
                all_vs = self.enum_vs(tmp_buf, tmp_cur, all_vs)
        return all_vs
    def find_vs(self, bin_id, vcnew):
        vs = []
        for i_dim, vs_rg in enumerate(vcnew.lambda_ranges[1:]):
            tmp_vs = []
            for tmp, rg in enumerate(vs_rg[1:]):
                vs_id = tmp+1
                if (bin_id[i_dim]+0.5)*self.bin_width > rg[0] and (bin_id[i_dim]+0.5)*self.bin_width <= rg[1]:
                    tmp_vs.append(vs_id)
            vs.append(tmp_vs)
        all_vs = self.enum_vs(vs, [], [])
        return all_vs
    def set_qcano_to_vs(self, vcnew):
        for bin_id, val in self.distrib.items():
            vs = self.find_vs(bin_id, vcnew)
            for i_vs in vs:
                if not tuple(i_vs) in vcnew.params:
                    vcnew.params[i_vs] = [0.0]
                vcnew.params[i_vs][0] += self.distrib[bin_id]/float(len(vs))

        vcnew.normalize_params()
        return vcnew
def _main():
    opts, args = opt_parse()

    vcdat = VcMDData(opts.bin_width)
    vcdat.set_trajectory_files(opts.stages, opts.series)

    vcdat.calc_canonical()
    
    if opts.fn_o_canonical:
        vcdat.write_canonical(opts.fn_o_canonical)
    if opts.fn_out:
        vcnew = kkmm_vcmd.VcMDConf()
        vcnew.read_params(opts.fn_new_vs)
        vcnew = vcdat.set_qcano_to_vs(vcnew)
        kkmm_vcmd.VcMDParamsWriter(opts.fn_out).write(vcnew)
    return

if __name__ == "__main__":
    _main()
