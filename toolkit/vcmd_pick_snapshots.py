#!/usr/bin/python3.7

import argparse
import sys
import os
import numpy as np
import re
import kkmm_vcmd
import copy
import random
import kkgro_trr
import kkpdb
import vcmd_get_state

# from collections import defaultdict

def argparser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i_weight_list', help='')
    parser.add_argument('--i_trr_list', help='')
    parser.add_argument('--i_pdb', help='')
    parser.add_argument('--i_param_temp', help='')
    parser.add_argument('--n_struct', type=int,  help="")
    parser.add_argument('--dir_out', help="")
    parser.add_argument('--n_bins', nargs="*", type=int, help='')
    parser.add_argument('--max_prob', type=float, help='')
    parser.add_argument('--sample_min_lambda', type=float, nargs="*", help='')
    parser.add_argument('--sample_max_lambda', type=float, nargs="*", help='')
    parser.add_argument('--uniform', action="store_true", help='')
    args = parser.parse_args()
    return args

class LambdaGrid:
    def __init__(self, vcmd, n_bins):
        self.vcmd = vcmd
        self.n_bins = n_bins
        self.lambda_min = []
        self.lambda_max = []
        self.lambda_min = np.zeros(self.vcmd.dim)
        self.lambda_max = np.zeros(self.vcmd.dim)

        self.grid_distrib = {}
        ## grid_distrib[(i_1, i_2, ..., i_dim)] = snapshot count
        ##   i_1 ... the
        self.grid_ss = {}
        ## grid_ss[(i_1, i_2, ..., i_dim)] = [(file_id, snapshot_id), (), (), ..]
        self.n_samples = 0
        return

    def mesh_distribution(self, weight_list, sample_min_lambda, sample_max_lambda):
        print("mesh_distribution")
        for dim in range(1, self.vcmd.dim+1):
            print(dim, self.vcmd.n_vs[dim-1])
            self.lambda_min[dim-1] = np.array(self.vcmd.lambda_ranges[dim][1][0])
            self.lambda_max[dim-1] = np.array(self.vcmd.lambda_ranges[dim][self.vcmd.n_vs[dim-1]][1])

        print("- lambda range")
        print(self.lambda_min)
        print(self.lambda_max)

        for file_id, i_weight in enumerate(weight_list):
            print(file_id)
            f = open(i_weight)
            for line_id, line in enumerate(f):
                terms = line.strip().split()
                weight = float(terms[1])
                lmb = np.array([float(x) for x in terms[2:2+self.vcmd.dim]])
                lmb01 = (lmb - self.lambda_min)/(self.lambda_max - self.lambda_min)
                lmb_bin = tuple([ int(x) for  x in lmb01 / (1.0/np.array(self.n_bins)) ])

                flg = True
                for dim, l in enumerate(lmb):
                    if dim < len(sample_min_lambda) and l < sample_min_lambda[dim]:
                        flg = False
                    if dim < len(sample_max_lambda) and l >= sample_max_lambda[dim]:
                        flg = False
                if flg:
                    if not lmb_bin in self.grid_distrib:
                        self.grid_distrib[lmb_bin] = 0
                        self.grid_ss[lmb_bin] = []
                    self.grid_distrib[lmb_bin] += 1
                    self.grid_ss[lmb_bin].append((file_id, line_id, lmb))
                    self.n_samples += 1
        return

    def pick_structures(self, n_struct, max_prob, uniform):
        print("pick_structures")
        picked = []
        distrib_tmp = sorted(self.grid_distrib.items(), key=lambda x:x[1], reverse=False)
        #ratio = self.n_samples/distrib_tmp[0][1]
        #pseudo = 0
        #if ratio > max_prob:
        #    pseudo = distrib_tmp[0][1] * ratio/max_prob
        #n_samples_p = self.n_samples + pseudo*len(distrib_tmp)
        #print([(x[1]+pseudo)/n_samples_p for x in distrib_tmp ])
        distrib = []
        if uniform:
            distrib = np.array([ 1.0 for x in distrib_tmp ])
        else:
            distrib = np.array([ self.n_samples/(x[1]) for x in distrib_tmp ])
        distrib /= np.sum(distrib)
                
        distrib_lmb = [ x[0] for x in distrib_tmp ]

        print("distrib")
        print(distrib_tmp)
        print(distrib)
        for i_strunct in range(n_struct):
            val = random.random()
            weight_acc = 0
            for lmb_bin, weight in zip(distrib_lmb, distrib):
                weight_acc += weight
                if val < weight_acc:
                    rnd = random.randrange(len(self.grid_ss[lmb_bin]))
                    picked.append(self.grid_ss[lmb_bin][rnd])
                    print(val, weight)
                    break
        picked_re = {}
        for file_id, ss_id, lmb in picked:
            if not file_id in picked_re:
                picked_re[file_id] = []
            picked_re[file_id].append((ss_id, lmb))
        return picked_re


def read_fnlist(fn):
    f = open(fn)
    fnlist = []
    for line in f:
        tmp = line.strip().split()[0]
        fnlist.append(tmp)
    f.close()
    return fnlist

def add_cano_weight(vcmd, i_weight):
    f = open(i_weight)
    vcmd.params = {}    
    for line in f:
        terms = line.strip().split()
        weight = float(terms[1])
        lmb = [float(x) for x in terms[2:2+vcmd.dim]]
        vs_set = vcmd.get_states_for_lambda(lmb)
        for vs in vs_set:
            if not vs in vcmd.params:
                vcmd.params[vs] = [0.0]
            vcmd.params[vs][0] += weight
    return vcmd

def gen_pdb_files(model, trr_list,
                  picked_re, dir_out, vcmd):
    print("gen_pdb_files")
    str_id = 0
    for file_id, ss_ids_lmb in picked_re.items():
        trr_reader = kkgro_trr.GroTrrReader(trr_list[file_id])
        trr_reader.open()
        for frame, lmb in ss_ids_lmb:
            str_id+=1

            frame_crd = trr_reader.read_nth_frame(frame)
            frame_crd.crds *= 10.0
            for i, at in enumerate(model.atoms):
                at.crd = np.array(frame_crd.crds[i])
            dir_str = dir_out+"/n"+str(str_id)
            os.makedirs(dir_str)
            fn_out = os.path.join(dir_str,"run.pdb")
            model.pbc_box[0] = frame_crd.box[0]
            model.pbc_box[1] = frame_crd.box[4]
            model.pbc_box[2] = frame_crd.box[8]
            model.pbc_box[3] = 90
            model.pbc_box[4] = 90
            model.pbc_box[5] = 90
            kkpdb.PDBWriter(fn_out).write_model(model)
            
            vs = vcmd_get_state.get_vs_candidates(lmb, vcmd)
            fn_vs = os.path.join(dir_str,"state.dat")
            vcmd_get_state.write_vstate(fn_vs, vs)

            print(dir_str, file_id, frame, lmb, vs)

        trr_reader.close()    
    return 
        
def _main():
    args = argparser()
    weight_list = read_fnlist(args.i_weight_list)
    vcmd = kkmm_vcmd.VcMDConf()
    vcmd.read_params(args.i_param_temp)

    grid = LambdaGrid(vcmd, args.n_bins)
    print(vcmd.lambda_ranges)


    vcmd.read_params(args.i_param_temp)
    grid.mesh_distribution(weight_list, args.sample_min_lambda, args.sample_max_lambda)
    picked_re = grid.pick_structures(args.n_struct, args.max_prob, args.uniform)
    trr_list = read_fnlist(args.i_trr_list)

    model = kkpdb.PDBReader(args.i_pdb).read_model()
    gen_pdb_files(model, trr_list, picked_re, args.dir_out, vcmd)

    return

if __name__ == "__main__":
    _main()
