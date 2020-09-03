#!/local/bin/python

import kkmcmdconf
from optparse import OptionParser

def get_options():
    p = OptionParser()
    
    p.add_option('-i', dest='fn_rangeinfo',
                 help="Input file")
    p.add_option('-o', dest='fn_out',
                 help="Output file")
    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

def generate_ttp_inp_mcmd(fn_range, fn_out, steps_trans=10000, temperature=300):
    ri = kkmcmdconf.RangeInfo(fn_range)
    ri.read()
    f = open(fn_out, "w") 
    n_vs = len(ri.vs_range.keys())
    f.write(";\n")
    f.write(str(n_vs) + "\n")  
    f.write(str(steps_trans) + "\n")
    f.write(";\n")
    for vsid in range(1, n_vs+1):
        f.write(str(ri.vs_range[vsid][0]) + "  " + str(ri.vs_range[vsid][1]) + "\n")
        f.write(str(ri.vs_range[vsid][2]) + "  " + str(ri.vs_range[vsid][3]) + "\n")
        f.write(";\n")
    #for params in self.mcmd_vs_fitfunc[phase]:
    #    f.write(str(len(params)-4) + "\n")
    #    for p in params[1:]:
    #        f.write(str(p)+"\n")
    f.write(str(temperature)+"\n")
    f.close()
    return

def divide_pdf_vs(fn_ttp_inp, fn_pdfs, pref_out="s", suff_out=".pdf"):
    ri = kkmcmdconf.RangeInfo(fn_ttp_inp)
    ri.read_ttp_inp()
    fn_out = {}
    f_out = {}
    dat = {}
    for i_vs, e_range in ri.vs_range.items():
        fn_o = pref_out+str(i_vs)+suff_out
        fn_out[i_vs] = fn_o
        f_out[i_vs] = open(fn_o, "w")
        dat[i_vs] = {}
    for f_pdf in fn_pdfs:
        f = open(f_pdf)
        for line in f:
            terms = line.strip().split()
            if len(terms) != 2: continue
            for i_vs, e_range in ri.vs_range.items():
                e = float(terms[0])
                if e >= e_range[0] and e < e_range[1]:
                    dat[i_vs][e] = line


    for i_vs, f_o  in f_out.items():
        for e in sorted(dat[i_vs].keys()):
            f_out[i_vs].write(dat[i_vs][e])
        f_o.close()
        
def _main():
    opts, args = get_options()    
    generate_ttp_inp_mcmd(opts.fn_rangeinfo,
                          opts.fn_out)

if __name__ == "__main__":
    _main()
