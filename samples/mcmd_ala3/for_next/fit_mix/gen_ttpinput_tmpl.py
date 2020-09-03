#!/local/bin/python

import sys
from optparse import OptionParser
#sys.path.append("/home/kasahara/local/kktools/kkio")
sys.path.append("/home/usr8/14IAW655/local/kktools/kkio")
import kkpresto_mcmd as kkpm

def get_options():
    p = OptionParser()
    
    p.add_option('-i', dest='fn_template',
                 help="file name for template ttp_v_mcmd.inp file")
    p.add_option('-o', dest='fn_out',
                 default="ttp_v_mcmd.inp",
                 help="file name for output ttp_v_mcmd.inp file")
    p.add_option('-s', dest='n_steps_trans',
                 type="int",
                 help="")
    p.add_option('-t', dest='temperature',
                 help="")
    p.add_option('-n', dest='n_vs',
                 type="int",
                 help="The number of virtual state")
    p.add_option('--range-template', dest='flg_range_template',
                 action="store_true",
                 help="Energy range written in the template file will be used.")
    p.add_option('--pref-inp', dest='pref_inp',
                 default=".fort.20",
                 help="Prefix of input files")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

class VirtualState:
    def __init__(self, deg, param, ene_range):
        self.deg = deg
        self.param = param
        self.ene_range = ene_range
        return

def gen_input(fn, vsobjects, n_steps_trans, temperature):
    f = open(fn, "w")
    f.write(";\n")
    f.write(str(len(vsobjects.keys())) + "\n")
    f.write(str(n_steps_trans) + "\n")
    f.write(";\n")
    for vsid, vs in vsobjects.items():
        f.write(str(vs.ene_range[0]) + "  ")
        f.write(str(vs.ene_range[1]) + "\n")
        if vsid == 1:
            f.write("0.0  1.0\n")
        elif vsid == len(vsobjects.keys()):
            f.write("1.0  0.0\n")
        else:
            f.write("1.0  1.0\n")
        f.write(";\n")
    for vsid, vs in vsobjects.items():    
        f.write(str(vs.deg) + "\n")
        for p in vs.param:
            f.write(str(p) + "\n")
    f.write(str(temperature) + "\n")
    f.close()
    return 0

def _main():
    opts, args = get_options()
    mcmd_tmpl = kkpm.PrestoMcmd
    if opts.fn_template:
        mcmd_tmpl = kkpm.PrestoMcmdReader(opts.fn_template).read()
    elif opts.n_vs:
        for i_vs in range(0, opts.n_vs):
            mcmd_tmpl.add_vs(0, [], (0.0, 0.0), (1.0, 1.0))
    if opts.temperature:
        mcmd_tmpl.temperature = opts.temperature
    if opts.n_steps_trans:
        mcmd_tmpl.n_steps_trans = opts.n_steps_trans

    states = {}
    for vs in range(0, len(mcmd_tmpl.vs)):
        if vs >= opts.n_vs: break
        filename = "s" + str(vs+1) + opts.pref_inp
        print filename
        f = open(filename)
        order = int(f.readline().strip())
        param = []
        for i_param in range(order + 3):
            param.append(f.readline().strip())
        ene_range = tuple(f.readline().strip().split())
        if not opts.flg_range_template:
            mcmd_tmpl.vs[vs].range = ene_range
        mcmd_tmpl.vs[vs].params = param
        mcmd_tmpl.vs[vs].order = order
        f.close()
    #gen_input(FN_OUTPUT,states, n_steps_trans, temperature)
    kkpm.PrestoMcmdWriter(opts.fn_out).write(mcmd_tmpl)
    return 0

if __name__ == "__main__":
    _main()
