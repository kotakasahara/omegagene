#!/usr/bin/python2.7
import os
import numpy as np
from optparse import OptionParser
import re
import struct as st
import sys
def get_options():
    p = OptionParser()
    
    p.add_option('-i', dest='dir_distrib',
                 default="dat_some",
                 help="Directory including trajectories")
    p.add_option('-o', dest='fn_out',
                 default="inp_c1_all",
                 help="Output file")
    p.add_option('--o-files-ene', dest="fn_o_files_ene",
                 default="fil.all",
                 help="Output file")
    p.add_option('--o-files-vs', dest="fn_o_files_vs",
                 default="filv.all",
                 help="Output file")
    p.add_option('-d', dest='omit_steps',
                 type="int", default=10000,
                 help="Number of steps to be omitted in each trajectory")
    p.add_option('-c', dest='power_coef',
                 type="float", default=2,
                 help="Coefficient for weighting values")
    p.add_option('-b', dest='binwidth',
                 type="float", default=100.0,
                 help="Bin width")
    p.add_option('-n', dest='n_run',
                 type="int", 
                 help="Number of runs")
    p.add_option('--celeste-bin', dest='flg_celeste_bin',
                 action="store_true",
                 help="Read Celeste binary data")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args
    

def read_multene_bin(fn):
    f = open(fn, "rb")
    buf = f.read(4)
    read_sign = ""
    if st.unpack(">i", buf)[0] == 66261:
        read_sign = ">"
    elif st.unpack("<i", buf)[0] == 66261:
        read_sign = "<"
    else:
        sys.stderr.write("File read error : " + fn + "\n")
        sys.stderr.write("The first four bytes are not 66261\n")
        return None
    buf = f.read(4)
    val_len = st.unpack(read_sign + "i", buf)[0]
    #print val_len
    read_val = read_sign
    if val_len==4:   read_val += "f"
    elif val_len==8: read_val += "d"
    buf = f.read(4)    
    n_col = st.unpack(read_sign + "i", buf)[0]

    frame = 0
    val = []
    while True:
        ene = None
        try:
            buf = f.read(val_len)
            #buf = f.read(8)
            #print buf
            ene = st.unpack(read_val, buf)[0]
            #print ene
        except:
            break
        if not ene: break
        val.append(ene)
    f.close()
    return np.array(val)

def read_multene(fn):
    try:
        f = open(fn)
        val = []
        for line in f:
            val.append(float(line))
        f.close()
    except:
        print "cannot open the file: " + str(fn)
        return None
    return np.array(val)

def get_files(dir_distrib):
    files_ene = {}
    files_vs = {}
    re_ene = re.compile("^e(\d+)$")
    re_vs = re.compile("^v(\d+)$")
    for fn in os.listdir(dir_distrib):
        m_ene = re_ene.match(fn)
        m_vs = re_vs.match(fn)
        if m_ene:
            files_ene[int(m_ene.group(1))] = os.path.join(dir_distrib, fn)
        elif m_vs:
            files_vs[int(m_vs.group(1))] = os.path.join(dir_distrib, fn)
    return files_ene, files_vs
                                                 
def output_filelist(files, fn, invalid):
    f_o = open(fn,"w")
    for f_key in sorted(files.keys()):
        if f_key in invalid: continue
        f_o.write(files[f_key]+"\n")
    f_o.close()
    return 

def get_all_summary(files, omit_steps, n_run, flg_bin):
    minene = 1e10
    maxene = -1e10
    summary = {}
    invalid_files = []
    stage_steps = {}
    traj_range = {}
    for num, fn in files.items():
        stage = (num-1) / n_run
        ene = []
        if flg_bin:
            ene = read_multene_bin(fn)        
        else:
            ene = read_multene(fn)        

        if len(ene)==0:
            invalid_files.append(num)
            continue
        stage_steps[stage] = len(ene)
        first = omit_steps
        for i_st in range(stage):
            first -= stage_steps[i_st]
        if first < 0: first = 0
        traj_range[num] = (first, len(ene))
        ene = ene[first:]
        summary[num] = (np.mean(ene), np.var(ene))
        if len(ene) == 0:
            print "ERROR: file does not contain anything"
            print fn
            continue
        tmpmin = np.min(ene)
        tmpmax = np.max(ene)
        if minene >= tmpmin: minene = tmpmin
        if maxene <= tmpmax: maxene = tmpmax
        print "Traj %d (stage:%d ,series%d) steps:%d-%d mean:%f min:%f max:%f var:%f"%(num,stage,(num-1)%n_run,first,len(ene)+first, summary[num][0], tmpmin, tmpmax, summary[num][1])
    return summary, minene, maxene, invalid_files, traj_range

def generate_inp(fn_out, summary, omit_steps, binwidth,
                 minene, maxene, power_coef, n_run, traj_range):

    bottom = (int(minene / binwidth) - 10) * binwidth
    n_bins = int((maxene - bottom)/binwidth + 10)
    n_files = len(summary.keys())
    f = open(fn_out, "w")
    f.write("%6.2f %5d\n"%(binwidth, n_bins))
    f.write("%10.2f\n"%bottom)
    f.write("END\n")
    f.write("%d\n"%n_files)
    max_steps = 990000000

    weight = {}
    print summary.keys()

    for num in summary.keys():
        if traj_range[num][1] - traj_range[num][0] > 0:
            weight[num] = np.power(summary[num][1], power_coef)

    print "min weight = " + str(np.min(weight.values()))

    min_weight = min(weight.values())
    for num in summary.keys():
        std_w = 1
        if traj_range[num][1] - traj_range[num][0] > 0:
            std_w = weight[num]/min_weight
        f.write("%4d %9d %10d %10d\n"%(
                num, traj_range[num][0], traj_range[num][1],
                int(std_w)))

    f.write("END\n")
    f.close()
    return 

def _main():
    opts, args = get_options()    
    files_ene, files_vs = get_files(opts.dir_distrib)
    summary, minene, maxene, invalid, traj_range = get_all_summary(files_ene, opts.omit_steps, opts.n_run, opts.flg_celeste_bin)

    for num, meanval in summary.items():
        print "%4d %10.3f %10.3f"%(num, meanval[0], meanval[1])

    output_filelist(files_ene, opts.fn_o_files_ene, invalid)
    output_filelist(files_vs, opts.fn_o_files_vs, invalid)

    generate_inp(opts.fn_out, summary,
                 opts.omit_steps, opts.binwidth,
                 minene, maxene,
                 opts.power_coef, opts.n_run, traj_range)

    
    return

if __name__ == "__main__":
    _main()

    
