#!/usr/bin/python2.7

## CONSTANTS

from optparse import OptionParser
import numpy as np
import re
import os
import sys
import struct as st
## sys.path.append(os.environ["HOME"] + "/local/kktools/mdtools")
import kkmcmdconf

def get_options():
    p = OptionParser()
    p.add_option('-i', dest='fn_inp',
                 default='inp_c1_all',
                 help="file name for input file")
    p.add_option('--i-ene-files', dest='fn_ene_files',
                 default='fil.all',
                 help="file name for the list of energy trajectories")
    p.add_option('--i-vs-files', dest='fn_vs_files',
                 default='filv.all',
                 help="file name for the list of vs trajectories")
    p.add_option('--i-ttpv', dest='fn_ttpvinp',
                 default='ttp_v_mcmd.inp',
                 help="file name for the psygene input file describing ttp v mcmd configurations")
    p.add_option('--pref-out', dest='pref_out',
                 default='v_pdf/s',
                 help="prefix for output dir")
    p.add_option('--ene-margin', dest='ene_margin',
                 type="float", default=0.0,
                 help="margin for ene distribution")
    p.add_option('--ene-interval', dest='ene_interval',
                 type="int", default=1,
                 help="interval frame number for reading energies")
    p.add_option('--ene-interval-file', dest='ene_interval_file',
                 type="int", default=1,
                 help="Interval frame number for recording in the MD runs.")
    p.add_option('--weight-vs', dest='flg_weight_vs',
                 action="store_true",
                 help="Weighting as equally sample for all vs")
    p.add_option('--celeste-bin', dest='flg_celeste_bin',
                 action="store_true",
                 help="Read Celeste binary data")

    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args
    
class DistribData(object):
    MAX_N_VS = 100   # the number of virtual systems (VS)
    N_WIN = 3000    # the number of energy bins in each VS
    N_FILMX = 5000  #
    V_CHANGE = 100000 #

    def __init__(self, ene_interval_file):
        self.resp = re.compile("\s+")
        self.binwidth = 50.0
        self.n_bins = 1
        self.min_ene = -1000000.0
        ## self.trajs[] = (num, first_step, last_step, weight)
        self.trajs = [] 
        
        self.filelist_ene = []
        self.filelist_vs = []
        self.max_vs = 0

        self.n_samples_vs = np.zeros(DistribData.MAX_N_VS, np.int)

        ## number of samples for each VS

        self.n_samples = 0
        self.n_weighted_samples = 0.0

        self.pdf = np.zeros((DistribData.MAX_N_VS, DistribData.N_WIN), np.float)
        
        self.ene_min_vs = np.zeros(DistribData.MAX_N_VS, np.float)
        self.ene_max_vs = np.zeros(DistribData.MAX_N_VS, np.float)
        self.n_bins_vs = np.zeros(DistribData.MAX_N_VS, np.int)

        self.ene_interval_file = ene_interval_file
        return

    def read_inp(self, fn):
        f = open(fn)

        line = f.readline().strip()
        terms = self.resp.split(line)
        self.binwidth = float(terms[0])
        self.n_bins = int(terms[1])

        line = f.readline().strip()
        self.min_ene = float(line)
        self.min_ene -= self.binwidth * 0.5

        line = f.readline().strip()
        assert(line=="END")

        line = f.readline().strip()
        n_trajs = int(line)
        
        for i in range(n_trajs):
            line = f.readline().strip()            
            terms = self.resp.split(line)
            self.trajs.append(([int(x) for x in terms]))
        
        line = f.readline().strip()
        assert(line=="END")

        f.close()
        return

    def read_lists(self, fn):
        strs = []
        f = open(fn)
        for line in f:
            strs.append(line.strip())
        f.close()
        return strs

    def read_ene_list(self, fn):
        self.filelist_ene = self.read_lists(fn)
    def read_vs_list(self, fn):
        self.filelist_vs = self.read_lists(fn)

    def read_vs_traj(self, i_traj, offset):
        frame_vs = []
        ## self.ipos = np.zeros(DistribData.V_CHANGE, np.int) ## 
        ## self.ivert = np.zeros(DistribData.V_CHANGE, np.int) ## 
        ## ipp = 0
        f = open(self.filelist_vs[i_traj])
        line1 = f.readline()
        line2 = f.readline()
        print "%5d %s %s"%(i_traj, self.filelist_ene[i_traj], self.filelist_vs[i_traj])
        print "    %8d %8d"%(self.trajs[i_traj][1], self.trajs[i_traj][2])
        cur_vs = -1;
        f.seek(0)
        for line in f:
            terms = self.resp.split(line.strip())
            cur_vs = int(terms[1])+offset
            if self.max_vs < cur_vs: self.max_vs = cur_vs
            frame_vs.append((int(terms[0]), cur_vs))
            ##ipp+=1
            ##self.ipos[ipp] = int(terms[0])
            ##self.ivert[ipp] = int(terms[1])
            assert(cur_vs <= DistribData.MAX_N_VS)
        print "N of data for virtual system = %d"%len(frame_vs)
        print "Virtual states : 1-" + str(self.max_vs)

        if len(frame_vs) == 1:
            print "  NOTE: A complete interval is not input"
            print "        Convert: Last step --> last of interval"
            frame_vs.append((self.trajs[i_traj][2], cur_vs))
        print frame_vs
        return frame_vs

    def read_ene_traj_bin(self, i_traj, frame_vs, ene_interval):
        f = open(self.filelist_ene[i_traj])
        buf = f.read(4)
        read_sign = ""
        if st.unpack(">i", buf)[0] == 66261:
            read_sign = ">"
        elif st.unpack("<i", buf)[0] == 66261:
            read_sign = "<"
        else:
            sys.stderr.write("File read error : " + self.filelist_ene[i_traj])
            sys.stderr.write("The first four bytes are not 66261")
            return
        buf = f.read(4)
        val_len = st.unpack(read_sign + "i", buf)[0]
        read_val = read_sign
        if val_len==4:   read_val += "f"
        elif val_len==8: read_val += "d"
        buf = f.read(4)
        n_col = st.unpack(read_sign + "i", buf)[0]
        frame = 0
        i_frame_vs = 0
        cur_vs = frame_vs[0][1]
        while True:
            ene = None
            try:
                buf = f.read(val_len)
                ene = st.unpack(read_val, buf)[0]
            except:
                break
            if len(buf) == 0: break
            frame += self.ene_interval_file
            if frame < self.trajs[i_traj][1]: continue
            if frame >= self.trajs[i_traj][2]: break
            ene = st.unpack(read_val, buf)[0]

            if frame >= frame_vs[i_frame_vs+1][0]:
                cur_vs = frame_vs[i_frame_vs+1][1]
                i_frame_vs += 1

            if frame % ene_interval != 0: continue

            self.n_samples += 1
            self.n_weighted_samples += self.trajs[i_traj][3]

            self.n_samples_vs[cur_vs] += 1
            if ene < self.min_ene: continue
            e1 = ene - self.min_ene
            jj = int(e1 / self.binwidth)
                #print "dbg " + str(jj)  + " " + str(self.n_bins)
                #print cur_vs
                #print self.trajs[i_traj][3]
                #sys.exit(1)
            if jj >= self.n_bins: continue
            self.pdf[cur_vs, jj] += self.trajs[i_traj][3]

        f.close()
        return 
    def read_ene_traj(self, i_traj, frame_vs, ene_interval):
        f = open(self.filelist_ene[i_traj])
        i_frame_vs = 0
        cur_vs = frame_vs[0][1]
        for frame, line in enumerate(f):
            ene = float(line.strip())
            if frame < self.trajs[i_traj][1]: continue
            if frame >= self.trajs[i_traj][2]: break

            #if len(frame_vs) != 1:
            #    for i_fvs, fvs in enumerate(frame_vs):
            #        if frame >= fvs[0]: cur_vs = fvs[1]
            #if cur_vs == -1: continue

            #if not i_frame_vs+1 in frame_vs: break
            if frame >= frame_vs[i_frame_vs+1][0]:
                cur_vs = frame_vs[i_frame_vs+1][1]
                i_frame_vs += 1
                #print frame_vs[i_frame_vs]

            if frame % ene_interval != 0: continue

            self.n_samples += 1
            self.n_weighted_samples += self.trajs[i_traj][3]

            self.n_samples_vs[cur_vs] += 1
                
            if ene < self.min_ene: continue
                
            e1 = ene - self.min_ene
            jj = int(e1 / self.binwidth)
                
            if jj >= self.n_bins: continue
            self.pdf[cur_vs, jj] += self.trajs[i_traj][3]
        f.close()
        return

    def set_min_max_vs(self):
        for vs in range(1, self.max_vs+1):
            self.ene_min_vs[vs] = 1e10
            self.ene_max_vs[vs] = -1e10
            print "dbg2 " + str(self.pdf.shape[1])
            for i_bin in range(self.pdf.shape[1]): 
                if self.pdf[vs, i_bin] == 0.0: continue
                self.n_bins_vs[vs] += 1
                ene = self.min_ene + (i_bin * self.binwidth) + self.binwidth*0.5
                if ene < self.ene_min_vs[vs]:
                    self.ene_min_vs[vs] = ene
                if ene > self.ene_max_vs[vs]:
                    self.ene_max_vs[vs] = ene
            print "VS %d : %f - %f"%(vs, self.ene_min_vs[vs], self.ene_max_vs[vs])
            
        return

    def output_distrib(self, pref_out, ene_margin, ri):
        print "the total number of samples = %d"%(self.n_samples)

        for i_vs in range(1, self.max_vs+1):
            val = self.n_samples_vs[i_vs]
            #enumerate(self.n_samples_vs):
            #i_vs = i_vs0+1
            #print "DBG!! " + str(i_vs) + " " + str(self.max_vs)
            if i_vs > self.max_vs: continue
            f = open(pref_out+str(i_vs)+".pdf", "w")
            print "  the number of samples for VS=%d: %d"%(i_vs, val)
            for i_bin in range(self.n_bins):
                ene = self.min_ene + i_bin*self.binwidth + self.binwidth*0.5
                print "%f %d %f %f"%(ene, self.pdf[i_vs, i_bin], self.ene_min_vs[i_vs], self.ene_max_vs[i_vs])

                if self.pdf[i_vs, i_bin] > 0.0 \
                        and ene >= self.ene_min_vs[i_vs]+ene_margin \
                        and ene < self.ene_max_vs[i_vs]-ene_margin:
                    rel_freq = np.log(self.pdf[i_vs, i_bin] / self.n_weighted_samples)
                    f.write("%15.7e  %15.7e\n"%(ene, rel_freq))
                        ## and ene >= ri.vs_range[i_vs][0] 
                        ## and ene < ri.vs_range[i_vs][1]:
            f.close()
        return

    def read_traj(self, flg_bin, ene_interval):
        if flg_bin:
            for i_traj in range(len(self.trajs)):
                print "traj : " + str(i_traj)
                ipp = self.read_vs_traj(i_traj, 0)
                self.read_ene_traj_bin(i_traj, ipp, ene_interval)
        else:
            for i_traj in range(len(self.trajs)):
                print "traj : " + str(i_traj)
                ipp = self.read_vs_traj(i_traj, 0)
                self.read_ene_traj(i_traj, ipp, ene_interval)

        return 
    def read_ttpv(self, fn_ttpv):
        ri = kkmcmdconf.RangeInfo(fn_ttpv)
        ri.read_ttp_inp()
        return ri
    def weight_vs(self, ri):

        #sum_pdf = np.sum(self.pdf, axis=1)
        #print "weight_vs"
        #print np.sum(self.pdf)
        #print sum_pdf
        #print self.pdf.shape
        #print sum_pdf.shape
        self.n_weighted_samples = 0.0

        for i_vs in range(1, self.max_vs+1):
            sum_pdf = 0.0
            n_bins_inrange = 0
            for i_bin in range(self.pdf.shape[1]): 
                ene = self.min_ene + (i_bin * self.binwidth) + self.binwidth*0.5
                if ene >= ri.vs_range[i_vs][0] and ene < ri.vs_range[i_vs][1]:
                    sum_pdf += self.pdf[i_vs, i_bin]
                    n_bins_inrange += 1
            coef = float(n_bins_inrange) / float(sum_pdf)
            print "vs %d n_bins:%d sum_pdf:%d"%(i_vs, n_bins_inrange, sum_pdf)
            for i_bin in range(self.pdf.shape[1]): 
                self.pdf[i_vs, i_bin] *= coef
                ene = self.min_ene + (i_bin * self.binwidth) + self.binwidth*0.5
                if ene >= ri.vs_range[i_vs][0] and ene < ri.vs_range[i_vs][1]:
                    self.n_weighted_samples += self.pdf[i_vs, i_bin]
        return 

def _main():
    opts, args = get_options()
    data = DistribData(opts.ene_interval_file)
    print "data.read_inp"
    data.read_inp(opts.fn_inp)
    print "data.read_ene_list"
    data.read_ene_list(opts.fn_ene_files)
    print "data.read_vs_list"
    data.read_vs_list(opts.fn_vs_files)
    print "data.read_traj"
    data.read_traj(opts.flg_celeste_bin, opts.ene_interval)
    print "set_min_max_vs"
    data.set_min_max_vs()
    ri = data.read_ttpv(opts.fn_ttpvinp)
    if opts.flg_weight_vs:
        data.weight_vs(ri)
    print "output_distrib"
    data.output_distrib(opts.pref_out, opts.ene_margin, ri)
        
if __name__ == "__main__":
    _main()
