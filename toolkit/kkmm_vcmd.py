#!/usr/bin/python2.7

import sys
import kkkit
import re
import numpy as np
import copy

class VcMDConf():
    def __init__(self):
        self.interval = 0
        #self.n_params = 0
        self.dim = 0
        self.group_names = []
        # group_names[dim] = [name, name]
        # self.group_names[0] is must be empty ("")

        self.n_vs = []
        self.lambda_ranges = []
        # lambda_ranges[dim][vsid] = (min, max)
        self.params = {}
        # params[(vs1, vs2, ...)] = (param1, param2)

        self.init_vs = []
        # init_vs[dim] = vsid

        self.seed = -1
    def sum_params(self, param_od=0):
        s = 0.0;
        for k,v in self.params.items():
            s += v[param_od]
        return s;
    def read_params(self, fn, chk4gen=True):
        self.interval, self.dim, self.group_names, \
            self.lambda_ranges, self.params, self.n_vs = VcMDParamsReader(fn).read(chk4gen)
    def read_init(self, fn):
        self.init_vs, self.seed = VcMDInitReader(fn).read(self.dim)
    def add_params(self, conf):
        key = tuple([ 0 for x in range(self.dim)])
        for vs, param in conf.params.items():
            if vs == key: continue
            if not vs in self.params:
                self.params[vs]  = []
                for i, p in enumerate(conf.params[vs]):
                    self.params[vs].append(p)
            else:
                for i, p in enumerate(conf.params[vs]):
                    self.params[vs][i] += p
        return
    def scale_params(self, factor):
        key_def = tuple([ 0 for x in range(self.dim)])
        for vs, param in self.params.items():
            if vs == key_def: continue
            self.params[vs][0] *= factor
        return
    def multiply_params(self, conf):
        key_def = tuple([ 0 for x in range(self.dim)])
        #print "test"
        #print self.params[key]

        # add default value for VSs in conf
        for vs, param in self.params.items():
            if vs == key_def: continue
            if not vs in conf.params:
                conf.params[vs] = conf.params[key_def]

        conf.normalize_params()
        #for vs, param in conf.params.items():
        #    if vs == key_def: continue
        #    if not vs in self.params:
        #        self.params[vs] = self.params[key_def]
        #self.normalize_params()

        for vs, param in conf.params.items():
            if not vs in self.params: continue
            if vs == key_def: continue
            for i, p in enumerate(conf.params[vs]):
                if self.params[vs][i] > 0:
                    if p > 0:
                        self.params[vs][i] *= p
                    else:
                        self.params[vs][i] *= self.params[vs][i]
            #print self.params[vs]
        return
    def normalize_params(self):
        p_sum = np.zeros(len(self.params.values()[0]), dtype=np.float)
        key_def = tuple([ 0 for x in range(self.dim)])
        for vs, param in self.params.items():
            if vs == key_def: continue
            p_sum += np.array(param)
        print("p_sum")
        print( p_sum)
        p_sum_t = np.zeros(len(self.params.values()[0]), dtype=np.float)
        for vs, param in self.params.items():
            if vs == key_def: continue
            #for i, q in enumerate(param):
            self.params[vs] /= p_sum
            p_sum_t += param
        print(p_sum_t)
        return
    def pop_zero_vs(self):
        key_def = tuple([ 0 for x in range(self.dim)])
        for vs, param in self.params.items():
            if vs == key_def: continue
            if self.params[vs][0] == 0:
                self.params.pop(vs)
        return

    def add_const(self, const):
        for vs, param in self.params.items():
            for i, p in enumerate(param):
                self.params[vs][i] += const
        return
    def set_default_param(self):
        key_def = tuple([ 0 for x in range(self.dim)])
        #if key_def in self.params: return
        min_param = 1e10
        print(len(self.params.items()))
        for k, v in self.params.items():
            if k == key_def: continue
            if v[0] < min_param and v[0] > 0:
                min_param = v[0]
        if min_param == 1e10:
            min_param = 1
        self.params[key_def] = [min_param]
        #self.params[key_def].append(min_param)
        # print min_param
        return
    def symmetrize(self):
        vs_param01 = {}
        vs_num = {}
        vs_vs = {}
        for key, val in self.params.items():
            uvs = list(copy.deepcopy(key))
            uvs.sort()
            uvs = tuple(uvs)
            if not uvs in vs_num:
                vs_num[uvs] = 0
                vs_param01[uvs] = 0.0
                vs_vs[uvs] = []
            vs_num[uvs] += 1
            vs_param01[uvs] += val[0]
            vs_vs[uvs].append(key)
        for uvs, vss in vs_vs.items():
            sym_param01 = vs_param01[uvs]/float(vs_num[uvs])
            for i_vs in vss:
                self.params[i_vs][0] = sym_param01
        return
    def statistics(self):
        buf = []
        key_def = tuple([ 0 for x in range(self.dim)])
        for k, v in self.params.items():
            if k == key_def: continue
            buf.append(v[0])
        param = np.array(buf)
        print("min: " + str(param.min()))
        print("max: " + str(param.max()))
        print("mean: " + str(param.mean()))
        print("sd: " + str(param.std()))
        return
    def is_in_range(self, vs, lmb):
        assert(len(vs)==self.dim and len(lmb)==self.dim)
        flg = True
        for d in range(self.dim):
            if lmb[d] < self.lambda_ranges[d+1][vs[d]][0] or \
               lmb[d] >= self.lambda_ranges[d+1][vs[d]][1]:
                flg=False; break
        return flg

class VcMDInitReader(kkkit.FileI):
    def __init__(self, fn):
        super(VcMDInitReader, self).__init__(fn)
        return
    def read(self, in_dim):
        self.open()
        init_vs = [0]
        ## The fist line:  The number of dimension
        dim = int(self.readline_comment().strip())
        ## The initial VS for each dimension, in each line
        for i in range(dim):
            tmp = int(self.readline_comment().strip())
            init_vs.append(tmp)
        # Random Seed
        seed = int(self.readline_comment().strip())
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
    def read(self, chk4gen=True):
        self.open()
        params = {}

        # The first line: VS transition interval (step)
        interval = int(self.readline_comment().strip().split()[0])
        # The second line: the number of dimensions
        dim = int(self.readline_comment().strip().split()[0])
        lambda_ranges = [(0.0, 0.0)]
        n_states = 1
        group_names = [[""]]
        n_vs_l = []
        # Definitions of VS ranges in each dimension
        for i in range(dim):
            cur_dim = i+1
            # [The number of VS] [Group name] [Group name]
            terms = self.readline_comment().strip().split()
            n_vs = int(terms[0])
            n_vs_l.append(n_vs)
            group_names.append([])
            for tm in terms[1:]:
                group_names[-1].append(tm)
            cur_ranges = [(0,0)]
            for j in range(n_vs):
                # [Min lambda] [Max lambda] 
                terms = self.readline_comment().strip().split()
                assert(len(terms) == 2)
                try:
                    lmin = float(terms[0])
                    lmax = float(terms[1])
                except:
                    sys.stderr.write("Format error in the VS definition.\n")
                    sys.stderr.write("\t".join(terms))
                    sys.stderr.write("The minimum and maximum values of lambda in each VS should be specified in float values.\n")
                    sys.exit(1)
                cur_ranges.append((float(terms[0]), float(terms[1])))
            lambda_ranges.append(cur_ranges)
            #print "dim " + str(cur_dim)
            #print group_names[-1]
            #print n_vs
            #print cur_ranges
        #print "n_states : " + str(n_states)
        while 1:
            terms = self.readline_comment().strip().split()
            if not terms or re.match("end", terms[0], re.IGNORECASE):
                break
            # elif re.match("default", terms[0], re.IGNORECASE):
            # default_q = float(terms[1])
            crd = tuple([int(x) for x in terms[:dim]])
            assert(not crd in params)
            param = [float(x) for x in terms[dim:]]
            try:
                assert(len(param) < 2)
            except:
                sys.stderr.write("WARNING: In the current version, only one parameter for each VS is allowed. The parameters except for the first are ignored.\n")
                sys.stderr.write("\t".join(terms)+"\n")

            if chk4gen:
                try:
                    assert(0 in crd or param[0] > 0)
                except:
                    sys.stderr.write("WARNING: The parameter shoud be larger than zero.\n")
                    sys.stderr.write("\t".join(terms)+"\n")

            params[crd] = param

        return interval, dim, group_names, lambda_ranges, params, n_vs_l

            
