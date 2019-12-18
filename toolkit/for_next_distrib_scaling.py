#!/usr/bin/python2.7
import sys
import os
import numpy as np

from optparse import OptionParser
import re
import sys
sys.path.append(os.environ["HOME"] + "/local/kktools/mdtools")
import kkmcmdconf

def get_options():
    p = OptionParser()
    
    p.add_option('-i', dest='dir_pdf',
                 default="v_pdf",
                 help="Directory including distributions")
    p.add_option('-o', dest='dir_out',
                 default="v_pdf_",
                 help="Output directory")
    p.add_option('--i-ttpv', dest='fn_ttpv',
                 help="ttp_v_mcmd.inp")
    p.add_option('-f', dest='factor',
                 type="float", 
                 help="scaling factor")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args

def read_pdf(fn_pdf):
    f = open(fn_pdf)
    pdf = []
    for line in f:
        terms = line.split()
        pdf.append([float(x) for x in terms])
    f.close()
    return pdf

def write_pdf(fn_out, exp_pdf):
    f = open(fn_out, "w")
    for terms in exp_pdf:
        f.write("%10.8e %10.8e\n"%(terms[0], np.log(terms[1])))
    f.close()
    return 

def cal_mean(exp_pdf, vs_range):
    exp_pdf_val = []
    for i_vs in vs_range.keys():
        exp_pdf_val.extend([ x[1] for  x in exp_pdf[i_vs] if x[0] >= vs_range[i_vs][0] and x[0] <= vs_range[i_vs][1]])
    return np.mean(exp_pdf_val)

def _main():
    opts, args = get_options()    
    ri = kkmcmdconf.RangeInfo(opts.fn_ttpv)
    ri.read_ttp_inp()

    pdfs = {}
    exp_pdf = {}
    for i_vs in ri.vs_range.keys():
        print "vs: "  + str(i_vs)
        fn_pdf = os.path.join(opts.dir_pdf, "s"+str(i_vs)+".pdf")
        pdfs[i_vs] = read_pdf(fn_pdf)
        exp_pdf[i_vs] = [[x[0], np.exp(x[1])] for x in pdfs[i_vs]]

    exp_pdf_mean = cal_mean(exp_pdf, ri.vs_range)

    print "mean : " + str(exp_pdf_mean)
    
    exp_pdf_sum = 0.0
    exp_pdf_scale = {}
    for i_vs in ri.vs_range.keys():
        exp_pdf_scale[i_vs] = [(x[0], exp_pdf_mean + (x[1]-exp_pdf_mean)*opts.factor )
                               for x in exp_pdf[i_vs] if x[0] >= ri.vs_range[i_vs][0] and x[0] <= ri.vs_range[i_vs][1]]
        exp_pdf_sum += np.sum([x[1] for x in exp_pdf_scale[i_vs]])
    
    print "exp_pdf_sum = " + str(exp_pdf_sum)

    for i_vs in ri.vs_range.keys():
        fn_out = os.path.join(opts.dir_out, "s"+str(i_vs)+".pdf")
        s_pdf = [ (x[0], x[1]/exp_pdf_sum) for x in exp_pdf_scale[i_vs] ]
        write_pdf(fn_out, s_pdf)
        
    return

if __name__ == "__main__":
    _main()

    
