#!/usr/bin/python2.6

from optparse import OptionParser
import sys
import kkmm_vcmd

def option_parse():
    p = OptionParser()
    
    
    p.add_option('--i-qcano', dest='fn_qcano',
                 help="Omegagene vcmd parameter file")
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
    p.add_option('--gmx', dest='flg_gmx',
                 action="store_true",
                 help="anayze gromacs output")
    p.add_option('--v57', dest='flg_v57',
                 action="store_true",
                 help="For assertion of the version.")

    opts, args = p.parse_args()
    p.print_help()

    return opts, args

def cal_prob(cano, vs, lmb):


    prob = {}
    frames = sorted(vs.keys())
    for frame in frames:
        #print frame
        cur_vs = vs[frame]
        if not frame in lmb:
            print("Lambda value is missing: " + str(frame))
            break
        cur_lmb = lmb[frame]
        cur_prob = 0.0
        if not cur_vs in cano.params:
            def_vs = tuple([ 0 for x in cur_vs ])
            cur_prob = cano.params[def_vs][0]
        else:
            cur_prob = cano.params[cur_vs][0]
        #print(cur_vs)
        n_overlapping_states = 0
        if not cano.is_in_range(cur_vs, cur_lmb):
            cur_prob = 0
        else:
            n_overlapping_states = cano.count_overlapping_states(cur_vs, cur_lmb)
        prob[frame] = cur_prob/float(n_overlapping_states)
    return prob

def read_dat(fn, itv_dat, itv_cod, first_step=0, skip=0, valtype="int"):
    dat={}
    f=open(fn)
    line_num = first_step
    for line in f:
        if line[0] == "#": continue
        buf = []
        if valtype=="int":
            buf = [int(x) for x in line.strip().split()]
        elif valtype=="float":
            buf = [float(x) for x in line.strip().split()]
        else:
            sys.stderr.write("Variable type error.")
            sys.exit(1)
        buf = buf[skip:]
        frame = line_num * itv_dat
        if frame % itv_cod == 0:
            dat[frame] = tuple(buf)

        line_num += 1

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

    if not opts.flg_v57:
        sys.stderr.write("--v57 option is required.\n")
        sys.exit()

    cano = kkmm_vcmd.VcMDConf()
    cano.read_params(opts.fn_qcano)
    
    vs = read_dat(opts.fn_vs, opts.itv_vs, opts.itv_cod, 1, 0, "int")
    lmb = []
    if opts.flg_gmx:
        lmb = read_dat(opts.fn_lmb, opts.itv_lmb, opts.itv_cod, 0, 1,"float")
    else:
        lmb = read_dat(opts.fn_lmb, opts.itv_lmb, opts.itv_cod, 1, 0, "float")

    prob = cal_prob(cano, vs, lmb)

    write_dat(opts.fn_out, vs, lmb, prob)
    
    

if __name__ == "__main__":
    _main()

