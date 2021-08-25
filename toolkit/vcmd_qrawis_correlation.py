#!/usr/bin/python3.7

import argparse
import kkmm_vcmd

def argparser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i_qrawis_list', nargs='*', help='')
    parser.add_argument('--i_qrawis', nargs='*', help='')
    parser.add_argument('--i_weight', nargs='+', help='')
    parser.add_argument('--out', help='')

    args = parser.parse_args()
    return args


def read_qrawis_list(fn):
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

def read_qrawis(qrawis_files, weight):
    vc = kkmm_vcmd.VcMDConf()
    vc.read_qraw_is(qrawis_files[0])
    if weight != []:
        vc.scale_qrawis(weight[0])

    for i, i_fn in enumerate(qrawis_files[1:]):
        vc_tmp = kkmm_vcmd.VcMDConf()        
        vc_tmp.read_qrawis(i_fn)
        #print(vc_tmp.qrawis)
        if weight != []:
            vc_tmp.scale_qrawis(weight[i+1])
        vc.sum_qrawis(vc_tmp)
    vc.gen_qraw_is_state()
    vc.sum_qraw_is_state()
    return vc

def cal_correlation_zone(vc_a, vc_b):
    sum_ab = 0.0
    sum_a = 0.0
    sum_b = 0.0
    sum_a2 = 0.0
    sum_b2 = 0.0
    n_dat = 0
    for k_state, qraw in vc_a.qraw.items():
        if not k_state in vc_b.qraw: continue
        n_dat += 1
        sum_a += vc_a.qraw[k_state]
        sum_b += vc_b.qraw[k_state]
        sum_ab += vc_a.qraw[k_state] * vc_b.qraw[k_state]
        sum_a2 += vc_a.qraw[k_state] * vc_a.qraw[k_state]
        sum_b2 += vc_b.qraw[k_state] * vc_b.qraw[k_state]
    ave_ab = sum_ab/n_dat
    ave_a = sum_a/n_dat
    ave_b = sum_b/n_dat
    ave_a2 = sum_a2/n_dat
    ave_b2 = sum_b2/n_dat
    sd_a = (ave_a2 - ave_a*ave_a)**0.5
    sd_b = (ave_b2 - ave_b*ave_b)**0.5
    cor = (ave_ab-ave_a*ave_b)/(sd_a*sd_b)
    return cor
        
def cal_correlation_subzone(vc_a, vc_b):
    cor_is = {}
    sum_cor = 0.0
    sum_cor2 = 0.0
    n_cor = 0
    n_complete_zones = 0
    for k_state, sub_dict in vc_a.qraw_is_state.items():
        sum_ab = 0.0
        sum_a = 0.0
        sum_b = 0.0
        sum_a2 = 0.0
        sum_b2 = 0.0
        n_dat = 0
        for k_is, qraw_val in sub_dict.items():
            if not k_state in vc_b.qraw_is_state or \
               not k_is in vc_b.qraw_is_state[k_state]:
                continue
            n_dat += 1
            rho_a = vc_a.qraw_is_state[k_state][k_is][0]/vc_a.qraw[k_state]
            rho_b = vc_b.qraw_is_state[k_state][k_is][0]/vc_b.qraw[k_state]
            sum_ab += rho_a * rho_b
            sum_a += rho_a
            sum_b += rho_b
            sum_a2 += rho_a * rho_a
            sum_b2 += rho_b * rho_b
        if n_dat != 2**vc_a.dim: continue
        n_complete_zones += 1
        ave_ab = sum_ab/n_dat
        ave_a = sum_a/n_dat
        ave_b = sum_b/n_dat
        ave_a2 = sum_a2/n_dat
        ave_b2 = sum_b2/n_dat
        sd_a = (ave_a2 - ave_a*ave_a)**0.5
        sd_b = (ave_b2 - ave_b*ave_b)**0.5
        cor = (ave_ab - ave_a*ave_b)/(sd_a*sd_b)
        sum_cor += cor
        sum_cor2 += cor * cor
        n_cor += 1
        cor_is[k_state] = cor
        
    ave_cor = sum_cor/n_cor
    ave_cor2 = sum_cor2/n_cor
    sd_cor = (ave_cor2 - ave_cor*ave_cor)**0.5
    return ave_cor, sd_cor, n_complete_zones

def main():
    args = argparser()
    if (args.i_qrawis and len(args.i_qrawis) != 2) or\
       (args.i_qrawis_list and len(args.i_qrawis_list) != 2) or\
       (not args.i_qrawis and not args.i_qrawis_list):
        sys.stderr.write("Two lists of qrawis files or two qrawis files are required for --i_qrawis_list or --i_qrawis options, repectively.")
        sys.exit(1)

    qrawis_files_a = []
    qrawis_files_b = []
    weight_a = []
    weight_b = []
    if len(args.i_qrawis) == 2:
        qrawis_files_a.append(args.i_qrawis[0])
        qrawis_files_b.append(args.i_qrawis[1])
    else:        
        qrawis_files_a = read_qrawis_list(args.i_qrawis_list[0])
        qrawis_files_b = read_qrawis_list(args.i_qrawis_list[1])
        if args.i_weight:
            weight_a = read_weight(i_weight[0])
            weight_b = read_weight(i_weight[1])
            if len(qrawis_files_a) != len(weight_a):
                stderr.write("The lengths of", args.i_qrawis_list[0], "and", i_weight[0], "differ.")
            if len(qrawis_files_b) != len(weight_b):
                stderr.write("The lengths of", args.i_qrawis_list[1], "and", i_weight[1], "differ.")


    vc_a = read_qrawis(qrawis_files_a,  weight_a)
    vc_b = read_qrawis(qrawis_files_b,  weight_b)

    sz_cor_ave, sz_cor_sd, n_complete_zones = cal_correlation_subzone(vc_a, vc_b)
    z_cor = cal_correlation_zone(vc_a, vc_b)
    print(z_cor, sz_cor_ave, sz_cor_sd, n_complete_zones)

    return

if __name__ == '__main__':
    main()
