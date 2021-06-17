#!/usr/bin/python3.7

import argparse
import kkmm_vcmd

def argparser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i_lambda', help='')
    parser.add_argument('--out', help='')
    parser.add_argument('--i_param', type=str, help='')
    parser.add_argument('--i_vs', type=str, help='')
#    parser.add_argument('--i_vs', action=store_true, help='')
    args = parser.parse_args()
    return args

def read_lambda(fn_lambda):
    f = open(fn_lambda)
    terms = []
    f.readline()
    lmd = []
    for line in f:
        terms = line.strip().split()
        lmd.append([ float(x) for x in  terms ][1:])
    f.close()
    return lmd

def read_vs(fn_vs):
    f = open(fn_vs)
    terms = []
    vs = []
    f.readline()
    for line in f:
        terms = line.strip().split()
        vs.append(tuple([ int(x) for x in  terms ]))
    f.close()
    return vs
            
def write_vstate(fn, state):
    fo = open(fn, "w")
    fo.write(str(len(state))+"\n")
    for st in state:
        fo.write(str(st)+"\n")    
    fo.close()
    return

def check_in_vs(lmb, vs, vcparams):
    flg = True
    #print(vcparams.lambda_ranges)
    for dim_tmp, c_vs in enumerate(vs):
        dim = dim_tmp + 1
        if lmb[dim-1] <=  vcparams.lambda_ranges[dim][c_vs][0] or \
           lmb[dim-1] >   vcparams.lambda_ranges[dim][c_vs][1]:
            flg = False
            break
    return flg

def get_vs_candidates(lmb, vcparams):
    vs_cand = []
    for dim_tmp, val_lmb in enumerate(lmb):
        dim = dim_tmp+1
        vs_cand_dim = []
        for vsid_tmp, rg in enumerate(vcparams.lambda_ranges[dim][1:]):
            vsid = vsid_tmp + 1
            #print(lmb, dim, rg)
            if val_lmb > rg[0] and val_lmb <= rg[1]:
                vs_cand_dim.append(vsid)

        if len(vs_cand_dim) == 0:
            return []
        elif len(vs_cand_dim) == 1:
            vs_cand.append(vs_cand_dim[0])
        elif len(vs_cand_dim) == 2:
            val1 = vcparams.lambda_ranges[dim][vs_cand_dim[0]][1] - val_lmb
            val2 = val_lmb - vcparams.lambda_ranges[dim][vs_cand_dim[1]][0]
            if val1 > val2:
                vs_cand.append(vs_cand_dim[0])
            else:
                vs_cand.append(vs_cand_dim[1])
    return tuple(vs_cand)


def main():
    args = argparser()
    lmd = read_lambda(args.i_lambda)[-1]
    vs = read_vs(args.i_vs)[-1]
    #print(lmd)
    #print(vs)
    vcparams = kkmm_vcmd.VcMDConf()
    vcparams.read_params(args.i_param)
    new_vs = []
    if check_in_vs(lmd, vs, vcparams):
        new_vs = vs
    else:
        print("state mod: ")
        new_vs = get_vs_candidates(lmd, vcparams)
    print(new_vs)
    write_vstate(args.out, new_vs)

            
if __name__ == '__main__':
    main()

