#!/usr/bin/python2.6

from optparse import OptionParser
import sys
import kkmm_vcmd

def option_parse():
    p = OptionParser()
    
    
    p.add_option('--i-qcano', dest='fn_qcano',
                 help="Omegagene vcmd parameter file")
    p.add_option('--i-cod', dest='fn_cod',
                 help="Omegagene coordinate file")
    p.add_option('--i-lmb', dest='fn_lmb',
                 help="Omegagene lambda file")
    p.add_option('--i-vs', dest='fn_vs',
                 help="Omegagene vs file")
    p.add_option('--interval-vs', dest='itv_vs',
                 type="int",
                 help="VS transition interval")
    p.add_option('--interval-lmb', dest='itv_lmb',
                 type="int",
                 help="Lambda output interval")
    p.add_option('--interval-cod', dest='itv_cod',
                 type="int",
                 help="Coordinate output interval")
    p.add_option('-o', dest='fn_out',
                 help="output file")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"

    return opts, args

def cal_prob(cano, vs, lmb):
    prob = {}
    frames = sorted(vs.keys())
    for frame in frames:
        #print frame
        cur_vs = vs[frame]
        if not frame in lmb:
            print "Lambda value is missing: " + str(frame)
            break
        cur_lmb = lmb[frame]
        cur_prob = cano.params[cur_vs][0]
        if not cano.is_in_range(cur_vs, cur_lmb):
            cur_prob = 0
        prob[frame] = cur_prob
    return prob

def read_dat(fn, itv_dat, itv_cod, valtype="int"):
    dat={}
    f=open(fn)

    for i, line in enumerate(f):
        buf = []
        if valtype=="int":
            buf = [int(x) for x in line.strip().split()]
        elif valtype=="float":
            buf = [float(x) for x in line.strip().split()]
        else:
            sys.stderr.write("Variable type error.")
            sys.exit(1)

        frame = (i+1) * itv_dat
        if frame % itv_cod == 0:
            dat[frame] = tuple(buf)

    f.close()
    return dat

def write_dat(fn_out, vs, lmb, prob):
    fo = open(fn_out,"w")
    frames = sorted(prob.keys())
    for frm in frames:
        #for c_vs, c_lmb, c_prob in zip(vs.values(),lmb.values(), prob):
        line=str(frm) + "\t" + str(prob[frm])+"\t"
        line+="\t".join([str(x) for x in lmb[frm]]) + "\t"
        line+="\t".join([str(x) for x in vs[frm]]) + "\n"
        fo.write(line)
    fo.close()
    return

def _main():
    opts, args = option_parse()

    cano = kkmm_vcmd.VcMDConf()
    cano.read_params(opts.fn_qcano)
    
    vs = read_dat(opts.fn_vs, opts.itv_vs, opts.itv_cod, "int")
    lmb = read_dat(opts.fn_lmb, opts.itv_lmb, opts.itv_cod, "float")
    prob = cal_prob(cano, vs, lmb)

    write_dat(opts.fn_out, vs, lmb, prob)
    
    

if __name__ == "__main__":
    _main()

