#!/usr/bin/python2.6

from optparse import OptionParser
import sys
import os
import numpy as np

import kkmm_vcmd

def opt_parse():
    p = OptionParser()
    
    p.add_option('--i-qraw', dest='fn_qraw',
                 default="vcmd_q_raw.txt",
                 help="q_raw file")
    p.add_option('--i-qcano', dest='fn_qcano',
                 default="vcmd_next.inp",
                 help="q_cano file")
    p.add_option('-o', dest='fn_out',
                 help="filename for output")    
    p.add_option('--stages', dest='stages',
                 type="int",
                 help="")
    p.add_option('--series', dest='series',
                 type="int",
                 help="")
    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

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
    fo = open(opts.fn_out,"w")
    for st in range(2, opts.stages+1):
        fn_cano = os.path.join("..", str(st-1), opts.fn_qcano)
        cano = kkmm_vcmd.VcMDConf()
        cano.read_params(fn_cano)
        for se in range(1, opts.series+1):
            fn_entire = os.path.join("..", str(st), "n"+str(se), opts.fn_qraw)
            print fn_entire
            if not os.path.exists(fn_entire):
                continue
            entire = kkmm_vcmd.VcMDConf()
            entire.read_params(fn_entire)
            for vs, param in entire.params.items():
                if param[0] > 0 and \
                        (vs[0], vs[1]+1) in entire.params and \
                        entire.params[(vs[0], vs[1]+1)][0] > 0 and \
                        (vs[0]+1, vs[1]) in entire.params and \
                        entire.params[(vs[0]+1, vs[1])][0] > 0 and \
                        (vs[0]+1, vs[1]+1) in entire.params and \
                        entire.params[(vs[0]+1, vs[1]+1)][0] > 0:
                    entire_p = np.zeros(4)
                    entire_p[0] = param[0]
                    entire_p[1] = entire.params[(vs[0], vs[1]+1)][0]
                    entire_p[2] = entire.params[(vs[0]+1, vs[1])][0]
                    entire_p[3] = entire.params[(vs[0]+1, vs[1]+1)][0]
                    cano_p = np.zeros(4)
                    cano_p[0] = cano.params[(vs[0], vs[1])][0]
                    cano_p[1] = cano.params[(vs[0], vs[1]+1)][0]
                    cano_p[2] = cano.params[(vs[0]+1, vs[1])][0]
                    cano_p[3] = cano.params[(vs[0]+1, vs[1]+1)][0]
                    sum_cano_p = cano_p.sum()
                    sum_entire_p = entire_p.sum()
                    qcano_pr = cano_p/sum_cano_p * entire_p/sum_entire_p
                    qcano_p = qcano_pr/qcano_pr.sum()

                    line = str(st) + "\t" + str(se) + "\t" 
                    line += str(vs[0]) + "\t" + str(vs[1]) + "\t" 
                    line += "\t".join([ str(x) for x in qcano_p])
                    line += "\n"
                    fo.write(line)
    fo.close()

if __name__ == "__main__":
    _main()
