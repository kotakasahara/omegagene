#!/usr/bin/python2.7

import kkkit

class ExpandConf():
    def __init__(self):
        self.vmcmd_range = {}
        self.vmcmd_params = {}
        self.init_vs = -1
        self.seed = -1
        self.temerature = 0.0
        self.interval = 0
        return
    def read_mcmdparams(self, fn):
        self.interval, self.vmcmd_range, \
            self.vmcmd_params, self.temperature = \
            McMDParamsReader(fn).read()
        return
    def read_init(self, fn):
        self.init_vs, self.seed = \
            McMDInitialDef(fn).read()
        return 

class McMDParamsReader(kkkit.FileI):
    def __init__(self, fn):
        super(McMDParamsReader, self).__init__(fn)
    def read(self):
        self.open() 
        #line = self.f.readline().strip()
        #terms = line.split()
        
        #self.ene_range = (float(terms[0]), float(terms[1]))
        #self.ene_width = self.ene_range[1] - self.ene_range[0]
        
        line = self.readline_comment().strip()

        n_vs = int(line.split()[0])
        intvl = int(self.readline_comment())
        ## vs_range[vs_id] = (min_ene, max_ene, tpro1, tpro2)
        vs_range = {}
        vs_params = {}        
        #vsid = 0
        for vsid in range(1,n_vs+1):
            line = self.readline_comment().strip()
            terms = line.split()
            ene_min = float(terms[0])
            ene_max = float(terms[1])
            line = self.readline_comment().strip()
            terms = line.split()
            p1 = float(terms[0])
            p2 = float(terms[1])
            vs_range[vsid] = (ene_min, ene_max, p1, p2)
        for vsid in range(1,n_vs+1):
            vs_params[vsid] = []
            order = int(self.readline_comment().strip())
            for i in range(order+3):
                vs_params[vsid].append(self.readline_comment().strip())
        temperature = float(self.readline_comment())    
        #print "DBG!: " + str(temperature)
        self.close()
        return intvl, vs_range, vs_params, temperature

class McMDInitialDef(kkkit.FileI):
    def __init__(self, fn):
        super(McMDInitialDef, self).__init__(fn)
        return
    def read(self):
        self.open()     
        init_vs = int(self.f.readline().strip())
        seed = int(self.f.readline().strip())
        return init_vs, seed
