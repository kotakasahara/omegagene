#!/usr/bin/python2.7

import sys
import kkkit
import re
import numpy as np

class VcMDConf():
    def __init__(self):
        self.interval = 0
        #self.n_params = 0
        self.dim = 0
        self.group_names = []
        # group_names[dim] = [name, name]
        # self.group_names[0] is must be empty ("")

        # lambda_ranges[dim] = (group A, group B)
        self.lambda_ranges = []
        # lambda_ranges[dim][vsid] = (min, max)
        self.params = {}
        # params[(vs1, vs2, ...)] = (param1, param2)

        self.init_vs = []
        # init_vs[dim] = vsid

        self.seed = -1
    def read_params(self, fn):
        self.interval, self.dim, self.group_names, \
            self.lambda_ranges, self.params = VcMDParamsReader(fn).read()
    def read_init(self, fn):
        self.init_vs, self.seed = VcMDInitReader(fn).read(self.dim)
    def add_params(self, conf):
        for vs, param in conf.params.items():
            for i, p in enumerate(conf.params[vs]):
                self.params[vs][i] += p
        return
    def multiply_params(self, conf):
        for vs, param in conf.params.items():
            for i, p in enumerate(conf.params[vs]):
                self.params[vs][i] *= p
        return
    def normalize_params(self):
        p_sum = np.zeros(len(self.params.values()[0]), dtype=np.float)
        for vs, param in self.params.items():
            p_sum += np.array(param)
        for vs, param in self.params.items():
            for i, q in enumerate(param):
                self.params[vs][i] /= p_sum[i]
        return
    def add_const(self, const):
        for vs, param in self.params.items():
            for i, p in enumerate(param):
                self.params[vs][i] += const
        return
        

class VcMDInitReader(kkkit.FileI):
    def __init__(self, fn):
        super(VcMDInitReader, self).__init__(fn)
        return
    def read(self, in_dim):
        self.open()
        init_vs = [0]
        dim = int(self.f.readline().strip())
        for i in range(dim):
            tmp = int(self.f.readline().strip())
            init_vs.append(tmp)
        seed = int(self.f.readline().strip())
        if not dim == in_dim:
            sys.stderr.write("Inconsistency in the definition of dimensions.\n")
            sys.stderr.write("VcMD paramter file:      " + str(in_dim) + "\n")
            sys.stderr.write("VcMD initial state file: " + str(dim) + "\n")
            sys.exit(1)
        #print "dbg kkmm_vcmd : seed " + str(seed) + "  dim " + str(dim)
        return init_vs, seed

class VcMDParamsWriter(kkkit.FileO):
    def __init__(self, fn):
        super(VcMDParamsWriter, self).__init__(fn)
    def write(self, vc):
        self.open()
        self.f.write(str(vc.interval)+"\n")
        self.f.write(str(vc.dim)+"\n")        
        for d in range(1, vc.dim+1):
            buf = ""
            buf += str(len(vc.lambda_ranges[d])-1)
            for name in vc.group_names[d]:
                buf += " " + name
            self.f.write(buf+"\n")
            for lmbd in vc.lambda_ranges[d][1:]:
                self.f.write(str(lmbd[0]) + " " + str(lmbd[1]) + "\n")
        keys = vc.params.keys()
        keys.sort()
        for vs in keys:
            param = vc.params[vs]
            # for vs, param in vc.params.items():
            buf = " ".join([str(x) for x in vs]) 
            for x in param:
                buf += " " + str(x)
            self.f.write(buf+"\n")
        self.f.write("END")
        self.close()

class VcMDParamsReader(kkkit.FileI):
    def __init__(self, fn):
        super(VcMDParamsReader, self).__init__(fn)
    def read(self):
        self.open()
        params = {}

        interval = int(self.readline_comment().strip().split()[0])
        dim = int(self.readline_comment().strip().split()[0])
        lambda_ranges = [(0.0, 0.0)]
        n_states = 1
        group_names = [[""]]
        #print "dbg kkmm_vcmd : dim " + str(dim) 

        for i in range(dim):
            cur_dim = i+1
            terms = self.readline_comment().strip().split()
            n_vs = int(terms[0])
            group_names.append([])
            for tm in terms[1:]:
                group_names[-1].append(tm)
            n_states *= n_vs
            cur_ranges = [(0,0)]
            for j in range(n_vs):
                terms = self.readline_comment().strip().split()
                cur_ranges.append((float(terms[0]), float(terms[1])))
            lambda_ranges.append(cur_ranges)
            #print "dim " + str(cur_dim)
            #print group_names[-1]
            #print n_vs
            #print cur_ranges

        for i in range(n_states):
            try:
                terms = self.readline_comment().strip().split()
            except:
                sys.stderr.write("An read error was occurred in VcMD param file\n")
                sys.stderr.write(fn+"\n")
                sys.stderr.write("The number of parameters may not be enough.\n")
                sys.stderr.write("The number of combinations of VS was "+pustr(n_states) +"\n")
                sys.exit(0)
            crd = tuple([int(x) for x in terms[:dim]])
            assert(not crd in params)
            param = [float(x) for x in terms[dim:]]
            params[crd] = param
        

        buf = self.readline_comment().strip()
        if not re.match("end", buf, re.IGNORECASE):
            sys.stderr.write("Missing the END keyword.\n")
            sys.exit(0)

        return interval, dim, group_names, lambda_ranges, params

            
