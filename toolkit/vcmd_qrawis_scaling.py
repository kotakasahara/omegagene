#!/usr/bin/python3.7

import argparse
import kkmm_vcmd

def argparser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i_qrawis', help='')
    parser.add_argument('--i_scale',  help='')
    parser.add_argument('--i_ref',  help='')
    parser.add_argument('--out', help='')
    parser.add_argument('--pseudocount', type=float, default=0, help='')

    args = parser.parse_args()
    return args

def main():
    args = argparser()
    if not args.i_qrawis or not args.i_scale or not args.out:
        sys.stderr.write("--i_qrawis: input qraw_is file")
        sys.stderr.write("--i_scale: input q_w file")
        sys.stderr.write("--out: output file")
        sys.exit(1)

    vc = kkmm_vcmd.VcMDConf()
    vc.read_qraw_is(args.i_qrawis)
    vc.read_params(args.i_scale)

    min_qrawis = 1e10
    for vsis, qrawis_v in vc.qraw_is.items():
        if min_qrawis > qrawis_v[0]: min_qrawis = qrawis_v[0]

    if args.i_ref:
        vcref = kkmm_vcmd.VcMDConf()
        vcref.read_qraw_is(args.i_ref)
        for vsis, qrawis_v in vcref.qraw_is.items():
            if vsis in vc.qraw_is:
                vc.qraw_is[vsis][0] += args.pseudocount
            elif args.pseudocount > 0:
                #vc.qraw_is[vsis] = [min_qrawis]
                vc.qraw_is[vsis] = [args.pseudocount]

    min_param = 1e10
    for vs, param in vc.params.items():
        if min_param > param[0]: min_param = param[0]

    for vsis, qrawis_v in vc.qraw_is.items():
        vs = tuple(vsis[:vc.dim])
        if vs in vc.params:
            vc.qraw_is[vsis][0] = qrawis_v[0] * vc.params[vs][0]
        else:
            vc.qraw_is[vsis][0] = qrawis_v[0] * min_param
            #sys.stdeff.write("The scaling parameter is missing for ", vs)

    kkmm_vcmd.VcMDParamsWriter(args.out).write(vc, 1, 0)
    return

if __name__ == '__main__':
    main()
