#!/usr/bin/python2.6

from optparse import OptionParser
import sys
import kkmm_vcmd

def opt_parse():
    p = OptionParser()
    
    p.add_option('-i', dest='fn_qcano',
                 help="q_cano file in the previous run")
    p.add_option('-o', dest='fn_out',
                 help="filename for output")
    p.add_option('--o-qraw', dest='fn_out_qraw',
                 help="filename for output")
    p.add_option('--i-qraw', dest='fn_qraw',
                 action="append",
                 help="q_raw files")
    p.add_option('--i-qraw-list', dest='fn_qraw_list',
                 help="q_raw files")
    p.add_option('--p-count', dest='pseudo_count',
                 type="int", default = 0,
                 help="pseudo_count")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

def read_fnlist(fn):
    f = open(fn)
    fnlist = []
    for line in f:
        tmp = line.strip().split()[0]
        fnlist.append(tmp)
    f.close()
    return fnlist

def _main():
    opts, args = opt_parse()
    vc = kkmm_vcmd.VcMDConf()
    fn_list = []
    if opts.fn_qraw:
        fn_list = opts.fn_qraw
    if opts.fn_qraw_list:
        fn_list.extend(read_fnlist(opts.fn_qraw_list))
    vc.read_params(fn_list[0])
    for fn_qraw in fn_list[1:]:
        vc_sub = kkmm_vcmd.VcMDConf()
        vc_sub.read_params(fn_qraw)
        vc.add_params(vc_sub)
    vc.add_const(opts.pseudo_count)
    vc.normalize_params()

    if opts.fn_out_qraw:
        kkmm_vcmd.VcMDParamsWriter(opts.fn_out_qraw).write(vc)

    vc_prev = kkmm_vcmd.VcMDConf()
    vc_prev.read_params(opts.fn_qcano)

    vc.multiply_params(vc_prev)
    vc.normalize_params()

    if opts.fn_out:
        kkmm_vcmd.VcMDParamsWriter(opts.fn_out).write(vc)


if __name__ == "__main__":
    _main()
