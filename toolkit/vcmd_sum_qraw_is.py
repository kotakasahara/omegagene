#!/usr/bin/python3.7

import argparse
import kkmm_vcmd

def argparser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i_qraw_list', help='')
    parser.add_argument('--i_weight', help='')
    parser.add_argument('--out', help='')

    args = parser.parse_args()
    return args


def read_qraw_list(fn):
    f = open(fn)
    files = []
    for line in f:
        terms = line.strip().split()
        files.append(terms[0])
    f.close()
    return files

def read_weight(fn):
    f = open(fn)
    w = []
    for line in f:
        terms = line.strip().split()
        w.append(float(terms[0]))
    f.close()
    return w


def main():
    args = argparser()
    qraw_files = read_qraw_list(args.i_qraw_list)
    weight = []
    if args.i_weight:
        weight = read_weight(args.i_weight)
        if len(qraw_files) != len(weight):
            stderr.write("The lengths of", args.i_qraw_list, "and", args.i_weight, "differ.")
    
    
    vc = kkmm_vcmd.VcMDConf()
    vc.read_qraw_is(qraw_files[0])
    if weight != []:
        vc.scale_qraw_is(weight[0])
    

    for i, i_fn in enumerate(qraw_files[1:]):
        vc_tmp = kkmm_vcmd.VcMDConf()        
        vc_tmp.read_qraw_is(i_fn)
        #print(vc_tmp.qraw_is)
        if weight != []:
            vc_tmp.scale_qraw_is(weight[i+1])
        vc.sum_qraw_is(vc_tmp)

    kkmm_vcmd.VcMDParamsWriter(args.out).write(vc, mode_param=1)
    
    return
    
            
if __name__ == '__main__':
    main()
