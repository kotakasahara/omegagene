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
    p.add_option('--i-qraw', dest='fn_list_qraw',
                 action="append",
                 help="q_raw files")
    p.add_option('--p-count', dest='pseudo_count',
                 type="int", default = 0,
                 help="pseudo_count")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

def _main():
    opts, args = opt_parse()
    vc = kkmm_vcmd.VcMDConf()
    vc.read_params(opts.fn_list_qraw[0])
    for fn_qraw in opts.fn_list_qraw[1:]:
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

    if opts.fn_out:
        kkmm_vcmd.VcMDParamsWriter(opts.fn_out).write(vc)


if __name__ == "__main__":
    _main()
