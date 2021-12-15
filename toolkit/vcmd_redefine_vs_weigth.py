#!/usr/bin/python3.7

import argparse
import sys
import os
import numpy as np
import re
import kkmm_vcmd
import copy
import kkmmconfig
# from collections import defaultdict

def argparser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i_weight_list', help='')
    parser.add_argument('--i_param_temp',  help='')
    parser.add_argument('--o_param_new', help='')
    args = parser.parse_args()
    return args


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
    for line in f:
        terms = line.strip().split()
        weight = float(terms[1])
        lmb = [float(x) for x in terms[2:2+vcmd.dim]]
        vs_set = vcmd.get_states_for_lambda(lmb)
        print(lmb, vs_set)
        for vs in vs_set:
            if not vs in vcmd.params:
                vcmd.params[vs] = [0.0]
            vcmd.params[vs][0] += weight
        
    return vcmd

def _main():
    args = argparser()

    #t = VcMDData(fn_weight_list,
    #fn_param_template)

    weight_list = read_fnlist(args.i_weight_list)
    vcmd = kkmm_vcmd.VcMDConf()
    vcmd.read_params(args.i_param_temp)
    vcmd.params = {}    
    for i_weight in weight_list:
        print(i_weight)
        vcmd = add_cano_weight(vcmd, i_weight)

    writer = kkmm_vcmd.VcMDParamsWriter(args.o_param_new)

    writer.write(vcmd)

    return

if __name__ == "__main__":
    _main()
