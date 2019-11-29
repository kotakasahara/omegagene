#!/usr/bin/python2.7

import kkkit
import sys
import re
import os
import math
import datetime
import shutil
import random

DEBUG = True

class McMDConf(object):
    ## domain  job_status[phase]
    # Phase_prerun = "pre"
    # Phase_mcmd = "mc"
    ## domain of job_status[phase][stage][series]
    State_wait = "waiting"
    State_run = "running"
    State_ready = "ready"
    State_finish = "finished"
    State_error = "error"
    ## type of fitting functions
    Func_poly = "poly"

    Rand_max = 99999999

    def __init__(self):
        self.fn_stem = "md"
        self.project_name = ""
        self.project_home = os.environ["HOME"]
        self.project_log = "kkmcmd.log"
        self.n_digit_run_id = 5

        d = datetime.datetime.now()
        self.random_seed = [d.microsecond]
        self.cell_size = [0.0, 0.0, 0.0]
        self.cell_division = [2, 2, 2]

        self.inp_dist_restraint = ""
        self.inp_topology = ""
        self.inp_shake = "system.shk"

        self.inp_md = "md.inp"
        self.inp_ttp = "ttp_v_mcmd.inp"
        self.inp_vert = "start.vert"

        self.out_mc_ene = "mult.ene"

        self.job_script = ""

        self.start_temp = 10
        self.heat_steps = 1000

        self.prerun_init = []
        self.prerun_force_coef = []
        self.prerun_ene_range = (-1e10, -1e3)
        self.prerun_md_inp_template = ""
        self.prerun_ttp_inp_template = ""

        self.cal_dir = []
        # cal_dir_stage[phase_id] = [dir_stage]
        self.cal_dir_stages = {}
        # cal_dir_series[phase_id] = [dir_series]
        self.cal_dir_series = {}

        self.dir_01_setting = "cal01_init"
        self.dir_02_prerun = "cal02_prerun"
        self.dir_03_mcmd = "cal03_mcmd"
        
        self.mcmd_stages = {}
        self.mcmd_inp_ttp = {}
        self.mcmd_vs_range = {}
        #self.mcmd_init_prev_phase = {}
        #self.mcmd_init_prev_stages = {}
        self.mcmd_init = {}
        self.mcmd_start_vs = {}
        self.mcmd_md_inp_template = {}
        self.mcmd_vs_fitfunc = {}
        self.mcmd_n_steps_transition = {}
        self.mcmd_temperature = {}
        self.dir_mcmd_stages = {}
        self.dir_mcmd_series = {}

        ## job_status[phase][stage][series] = state
        ##  state is one of ("runing", "finished" "waiting")
        ##  value is (jobid, state, path)
        self.job_status = {}

        self.max_running = 16
        self.n_running = 0
        return 

    def get_run_path(self, phase, stage, series):
        d = ""
        #print "get_run_path " + str(phase) + " " + str(stage) + " " + str(series)
        d = os.path.join(self.cal_dir_stages[phase][stage],
                         self.cal_dir_series[phase][series])
        return d

    def set_paths(self):
        self.job_script = self.job_script.replace("${PRJ_HOME}", self.project_home)
        self.job_script_fn = self.job_script.split("/")[-1]
        self.project_log = self.project_log.replace("${PRJ_HOME}", self.project_home)
        self.inp_dist_restraint = self.inp_dist_restraint.replace("${PRJ_HOME}", self.project_home)
        self.inp_topology = self.inp_topology.replace("${PRJ_HOME}", self.project_home)
        self.inp_shake = self.inp_shake.replace("${PRJ_HOME}", self.project_home)

        for ph, cal_dir in enumerate(self.cal_dir):
            self.cal_dir[ph] = cal_dir.replace("${PRJ_HOME}", self.project_home)

        self.dir_01_setting = self.dir_01_setting.replace("${PRJ_HOME}", self.project_home)
        self.dir_02_prerun = self.dir_02_prerun.replace("${PRJ_HOME}", self.project_home)
        self.dir_03_mcmd = self.dir_03_mcmd.replace("${PRJ_HOME}", self.project_home)

        # set prerun stages
        self.prerun_md_inp_template = self.prerun_md_inp_template.replace("${PRJ_HOME}", self.project_home)
        self.prerun_ttp_inp_template = self.prerun_ttp_inp_template.replace("${PRJ_HOME}", self.project_home)
        for i, fn in enumerate(self.prerun_init):
            self.prerun_init[i] = fn.replace("${PRJ_HOME}", self.project_home)


        for phase, vsrange in self.mcmd_vs_range.items():
            self.mcmd_vs_range[phase] = vsrange.replace("${PRJ_HOME}", self.project_home)
        for phase, mdtemp in self.mcmd_md_inp_template.items():
            self.mcmd_md_inp_template[phase] = mdtemp.replace("${PRJ_HOME}", self.project_home)
        for phase, mdtemp in self.mcmd_inp_ttp.items():
            self.mcmd_inp_ttp[phase] = mdtemp.replace("${PRJ_HOME}", self.project_home)

        ## prerun directories
        self.cal_dir_stages[1] = []
        for i, fc in enumerate(self.prerun_force_coef):
            st = "%"+str(self.n_digit_run_id)+"."+str(self.n_digit_run_id-2)+"f"
            dir_stage = str(i+1)
            self.cal_dir_stages[1].append(os.path.join(self.cal_dir[1], dir_stage))
        # set prerun runs
        self.cal_dir_series[1] = []
        for series_id in range(len(self.prerun_init)):
            st = "%0"+str(self.n_digit_run_id)+"d"
            dir_run = "n" + str(series_id+1)
            self.cal_dir_series[1].append(dir_run)

        # mcmd directories
            
        for phase, stg in self.mcmd_stages.items():
            self.cal_dir_stages[phase] = []
            #for i, pre in init_stages:
            #st = "%"+str(self.n_digit_run_id)+"."+str(self.n_digit_run_id-2)+"f"
            #dir_stage = "force_" + str(i) + "_" + st%fc
            for i in range(0,stg):
                dir_stage = str(i+1)
                self.cal_dir_stages[phase].append(os.path.join(self.cal_dir[phase], dir_stage))
        for phase, init in self.mcmd_init.items():
            self.cal_dir_series[phase] = []
            for series_id, ser in enumerate(init):
                dir_run = "n" + str(series_id+1)
                self.cal_dir_series[phase].append(dir_run)
        #for phase, init_stages in self.mcmd_init_prev_stages.items():
        #    self.cal_dir_series[phase] = []
        #    for series_id, ser in enumerate(init_stages):
        #    #st = "%0"+str(self.n_digit_run_id)+"d"
        #    #dir_run = "series_"+st%series_id
        #        dir_run = "n" + str(series_id+1)
        #    #if dir_run[0] == "~":
        #    #dir_run = os.path.join(self.project_home, dir_run)
        #        self.cal_dir_series[phase].append(dir_run)
        return
    def generate_ttp_inp(self, fn_out, force_coef):
        fi = open(self.prerun_ttp_inp_template)
        text = fi.read()
        fi.close()
        text = text.replace("#{FORCE_COEF}", str(force_coef))

        self.file_backup(fn_out)
        fo = open(fn_out,"w")
        fo.write(text)
        fo.close()
        return
    def generate_ttp_inp_mcmd(self, fn_out, phase, stage, series):
        ri = RangeInfo(self.mcmd_vs_range[phase])
        ri.read()
        d_run = self.get_run_path(phase, stage, series)

        self.file_backup(fn_out)
        f = open(fn_out, "w") 
        n_vs = len(self.mcmd_vs_fitfunc[phase])
        f.write(";\n")
        f.write(str(n_vs) + "\n")  
        f.write(str(self.mcmd_n_steps_transition[phase]) + "\n")
        f.write(";\n")
        for vsid in range(1, n_vs+1):
            f.write(str(ri.vs_range[vsid][0]) + "  " + str(ri.vs_range[vsid][1]) + "\n")
            f.write(str(ri.vs_range[vsid][2]) + "  " + str(ri.vs_range[vsid][3]) + "\n")
            f.write(";\n")

        for params in self.mcmd_vs_fitfunc[phase]:
            #if params[0] == self.Func_poly:
            f.write(str(len(params)-4) + "\n")
            for p in params[1:]:
                f.write(str(p)+"\n")
        f.write(str(self.mcmd_temperature[phase])+"\n")
        f.close()
        return
    def replace_md_inp(self, text):
        text = text.replace("#{CELL_SIZE_X}", str(self.cell_size[0]) + "d0")
        text = text.replace("#{CELL_SIZE_Y}", str(self.cell_size[1]) + "d0")
        text = text.replace("#{CELL_SIZE_Z}", str(self.cell_size[2]) + "d0")
        text = text.replace("#{CELL_DIVISION_X}", str(self.cell_division[0]))
        text = text.replace("#{CELL_DIVISION_Y}", str(self.cell_division[1]))
        text = text.replace("#{CELL_DIVISION_Z}", str(self.cell_division[2]))
        text = text.replace("#{INP_TOPOLOGY}", self.inp_topology)
        text = text.replace("#{INP_SHAKE}", self.inp_shake)
        text = text.replace("#{INP_DIST_RESTRAINT}", self.inp_dist_restraint)
        text = text.replace("#{INP_TTP}", self.inp_ttp)
        text = text.replace("#{INP_VERT}", self.inp_vert)
        text = text.replace("#{RANDOM_SEED}", str(random.randint(0, McMDConf.Rand_max)))
        text = text.replace("#{FILENAME_STEM}", self.fn_stem)
        text = text.replace("#{OUT_MC_ENE}", self.out_mc_ene)
        return text

    def generate_md_inp_prerun(self, fn_out, stage, series):
        phase = 1
        fi = open(self.prerun_md_inp_template)
        text = fi.read()
        fi.close()
        text = self.replace_md_inp(text)
        
        inp_pdb = ""
        if stage == 0:
            inp_pdb = self.prerun_init[series]
            inp_restart = self.prerun_init[series][:-3] + "restart"
        else:
            d_prev = self.job_status[phase][stage-1][series][2]
            inp_pdb = os.path.join(self.cal_dir_stages[1][stage-1],
                                   self.cal_dir_series[1][series],
                                   self.fn_stem+".pdb")
            inp_restart = os.path.join(self.cal_dir_stages[1][stage-1],
                                       self.cal_dir_series[1][series],
                                       self.fn_stem+".restart")
        text = text.replace("#{INP_PDB}", inp_pdb)
        text = text.replace("#{INP_RESTART}", inp_restart)

        self.file_backup(fn_out)
        fo = open(fn_out,"w")
        fo.write(text)
        fo.close()
        self.generate_vert_inp(phase,
                               stage, series, 1)
        return
    def read_final_ene(self, cal_dir):
        """
        Read energy at final step from mult.ene
        """
        fn = os.path.join(cal_dir, self.out_mc_ene)
        f = open(fn)
        ene = 0.0
        for line in f:  ene = float(line.strip())
        f.close()
        return ene
    def generate_vert_inp(self, phase, stage, series, vsid):
        dir_cal = self.get_run_path(phase, stage, series)
        fn_vert = os.path.join(dir_cal, self.inp_vert)

        self.file_backup(fn_vert)
        f = open(fn_vert,"w")
        f.write("    " + str(vsid) + "\n")
        f.write(str(random.randint(0, McMDConf.Rand_max))+"\n")
        f.close()
        return

    def generate_md_inp_mcmd(self, fn_out, phase, stage, series):
        #print "gen md inp " + fn_out
        fi = open(self.mcmd_md_inp_template[phase])
        text = fi.read()
        fi.close()
        text = self.replace_md_inp(text)
        
        inp_pdb = ""
        inp_restart = ""
        prev_dir = ""
        if stage == 0:
            prev_dir = os.path.join(self.cal_dir_stages[self.mcmd_init[phase][series][0]]\
                                        [self.mcmd_init[phase][series][1]-1],
                                    self.cal_dir_series[self.mcmd_init[phase][series][0]]\
                                        [self.mcmd_init[phase][series][2]-1])
            
            #prev_dir = os.path.join(self.cal_dir_stages[self.mcmd_init_prev_phase[phase]][self.mcmd_init_prev_stages[phase][series]-1],
                                    #self.cal_dir_series[self.mcmd_init_prev_phase[phase]][series])
            inp_pdb = os.path.join(prev_dir, self.fn_stem+".pdb")
            inp_restart = os.path.join(prev_dir, self.fn_stem+".restart")
            text = text.replace("#{INP_RESTART_FLG}", "NO")                
            heating =  "        INITIA= SET\n" 
            heating += "        STARTT= "+str(self.start_temp) + "\n"
            heating += "        HEATLO= "+str(self.heat_steps) + "\n"
            text = text.replace("#{HEAT_OPTIONS}", heating)                
        else:
            prev_dir = self.job_status[phase][stage-1][series][2]
            inp_pdb = os.path.join(prev_dir,
                                   self.fn_stem+".pdb")
            inp_restart = os.path.join(prev_dir,
                                       self.fn_stem+".restart")
            text = text.replace("#{HEAT_OPTIONS}", "")
            text = text.replace("#{INP_RESTART_FLG}", "YES")

        text = text.replace("#{INP_PDB}", inp_pdb)
        text = text.replace("#{INP_RESTART}", inp_restart)

        self.file_backup(fn_out)
        fo = open(fn_out,"w")
        fo.write(text)
        fo.close()

        ene = self.read_final_ene(prev_dir)
        init_vs = None
        if(stage > 0 or
           not phase in self.mcmd_start_vs or
           self.mcmd_start_vs[phase][series]=='*'):
            ri = None
            if phase in self.mcmd_vs_range:
                ri = RangeInfo(self.mcmd_vs_range[phase])
                ri.read()
            else:
                ri = RangeInfo(self.mcmd_inp_ttp[phase])
                ri.read_ttp_inp()
            init_vs = ri.get_near_state(ene)
        else:
            init_vs = int(self.mcmd_start_vs[phase][series]);
            
        print "prepare mcmd " + str(stage) + " " + str(series)
        print ene
        print init_vs
        self.generate_vert_inp(phase,
                               stage, series, init_vs)
        return

    def set_job_status(self):
        for phase, cdir in enumerate(self.cal_dir):
            if phase==0: continue 
            self.job_status[phase] = {}
            n_stages = 0
            n_series = 0
            if phase == 1:
                n_stages = len(self.prerun_force_coef)
                n_series = len(self.prerun_init)
            else:
                n_stages = self.mcmd_stages[phase]
                n_series = len(self.mcmd_init[phase])
            for stage in range(n_stages):
                self.job_status[phase][stage] = {}
                for series in range(n_series):
                    d_run = self.get_run_path(phase, stage, series)
                    self.job_status[phase][stage][series] = (-1, McMDConf.State_wait, d_run)
        return 
    def prepare_run(self, phase, stage, series):
        if phase == 1:
            self.prepare_run_prerun(stage, series)
        else:
            self.prepare_run_mcmd(phase, stage, series)
        return 
    def prepare_run_prerun(self, stage, series):
        phase = 1

        d_run = self.get_run_path(phase, stage, series)
        if not os.path.exists(d_run):
            os.makedirs(d_run)

        self.generate_ttp_inp(os.path.join(d_run, self.inp_ttp),
                              self.prerun_force_coef[stage])
        self.generate_md_inp_prerun(os.path.join(d_run, self.inp_md),
                                    stage, series)
        shutil.copyfile(self.job_script, os.path.join(d_run, self.job_script_fn))
        self.job_status[phase][stage][series] = (-1, McMDConf.State_ready, d_run)
        if DEBUG:
            print self.job_status[phase][stage][series]            
        return 

    def prepare_run_mcmd(self, phase, stage, series):
        d_run = self.get_run_path(phase, stage, series)
        if not os.path.exists(d_run): os.makedirs(d_run)
        fn_md_inp = os.path.join(d_run, self.inp_md)
        self.generate_md_inp_mcmd(fn_md_inp, phase, stage, series)

        fn_ttp_inp = os.path.join(d_run, self.inp_ttp)
        if phase in self.mcmd_inp_ttp:
            shutil.copyfile(self.mcmd_inp_ttp[phase], fn_ttp_inp)
        else:
            self.generate_ttp_inp_mcmd(fn_ttp_inp, phase, stage, series)

        shutil.copyfile(self.job_script, os.path.join(d_run, self.job_script_fn))
        self.job_status[phase][stage][series] = (-1, McMDConf.State_ready, d_run)        

    def dump_config(self):
        text  = ""
        text += "--fn-stem " + self.fn_stem + "\n"
        text += "--project-name " + self.project_name + "\n"
        text += "--project-home " + self.project_home + "\n"
        text += "--project-log " + self.project_log + "\n"
        text += "--n-digit-run-id " + str(self.n_digit_run_id) + "\n"
        text += "--random-seed " + str(self.random_seed) + "\n"
        text += "--cell-size " + " ".join([str(x) for x in self.cell_size]) + "\n"
        text += "--cell-division " + " ".join([str(x) for x in self.cell_division]) + "\n"
        text += "--inp-dist-restraint " + str(self.inp_dist_restraint)
        text += "; replace #{INP_DIST_RESTRAINT}"+ "\n"
        text += "--inp-topology " + str(self.inp_topology)
        text += "; replace #{INP_TOPOLOTY}"+ "\n"
        text += "--inp-shake " + str(self.inp_shake)
        text += "; replace #{INP_SHAKE}"+ "\n"
        text += "--inp-md " + str(self.inp_md) + "\n"
        text += "--inp-ttp " + str(self.inp_ttp)
        text += "; replace #{INP_TTP}"+ "\n"
        text += "--inp-vert " + str(self.inp_vert)
        text += "; replace #{INP_VERT}"+ "\n"
        text += "--out-mc-ene " + str(self.out_mc_ene)
        text += "; replace #{OUT_MC_ENE}"+ "\n"
        
        text += "--job-script " + self.job_script + "\n"
        for fn in self.prerun_init:
            text += "--prerun-init " + fn + "\n"

        for fc in self.prerun_force_coef:
            text += "--prerun-force-coef " + str(fc) + "\n"
        text += "--prerun-md-inp-template " + self.prerun_md_inp_template + "\n"
        text += "--prerun-ttp-inp-template " + self.prerun_ttp_inp_template + "\n"
        for d in self.dir_prerun_stages:
            text += "; --dir-prerun-stages " + d + "\n"
        for d in self.dir_prerun_series:
            text += "; --dir-prerun-series " + d + "\n"
        text += "--mcmd-stages " + self.mcmd_stages + "\n"
        for cal_dir in self.cal_dir:
            text += "--cal-dir " + cal_dir + "\n"
        text += "--dir-01-setting " + self.dir_01_setting + "\n"
        text += "--dir-02-prerun " + self.dir_02_prerun + "\n"
        text += "--dir-03-mcmd " + self.dir_03_mcmd + "\n"
        text += "--mcmd-vs-range " + self.mcmd_vs_range + "\n"
        text += "--mcmd-md-inp-template " + self.mcmd_md_inp_template + "\n"
        #text += "--mcmd-init-prev-phase " + self.mcmd_init_prev_phase + "\n"
        text += "--mcmd-init-prev-stages " + " ".join([str(x) for x in self.mcmd_init]) + "\n"

        for func in self.dir_vs_fitfunc:
            text += "--mcmd-vs-fitfunc " + " ".join([str(x) for x in func]) + "\n"
        text += "--mcmd-n-steps-transition " + str(self.mcmd_n_steps_transition) + "\n"
        text += "--mcmd-temperature " + str(self.mcmd_temperature) + "\n"
        #text += "--max-running " + str(self.mcmd_temperature) + "\n"
        text += ";job_status \n"
        for phase, stage_series_state in self.job_status.items():
            for stage, series_state in stage_series_state.items():
                for series, jobid_state in series_state.items():
                    jobid = jobid_state[0]
                    state = jobid_state[1]
                    d_run = jobid_state[2]
                    text += ";["+phase+"]-["+str(stage)+"]-["+str(series)+"] "
                    text += str(jobid) + " " + state + " " + d_run + "\n"
                    
        return text
    def file_backup(self, file_path):
        file_path_terms = file_path.split("/")
        filename = file_path_terms[-1]
        file_dir = "/".join(file_path_terms[:-1])
        max_bk_num = -1
        ## -1 : no file
        ## 0  : file exists but no backup files
        for fn in os.listdir(file_dir):
            m1 = re.compile("#"+filename+".(\d+)#").match(fn)
            if m1:
                bk_num = int(m1.group(1))
                if bk_num >= max_bk_num: max_bk_num = bk_num
            if fn == filename and max_bk_num == -1:
                max_bk_num = 0
        new_bk_num = max_bk_num + 1
        new_bk_fn = "#"+filename+"."+str(new_bk_num)+"#"
        if max_bk_num >= 0:
            shutil.copyfile(file_path, os.path.join(file_dir, new_bk_fn))
        return new_bk_fn


class McMDConfReader(kkkit.FileI):
    def __init__(self, fn):
        super(McMDConfReader, self).__init__(fn)
    def read_conf(self):
        conf = McMDConf()
        self.open()
        current_field = ""
        current_text = ""
        for orig_line in self.f:
            line = kkkit.eliminate_comment(orig_line).strip()
            terms = re.compile("\s+").split(line.strip())
            if len(line) == 0:
                continue
            elif terms[0][0:2] != "--":
                sys.stderr.write("Syntax error: Field tag must begin with '--':")
                sys.stderr.write(line + "\n")
            elif terms[0][2:] == "filename-stem":
                conf.fn_stem = terms[1]
            elif terms[0][2:] == "job-script":
                conf.job_script = terms[1]
            elif terms[0][2:] == "start-temp":
                conf.start_temp = float(terms[1])
            elif terms[0][2:] == "heat-steps":
                conf.heat_steps = int(terms[1])
            elif terms[0][2:] == "project-name":
                conf.project_name = terms[1]
            elif terms[0][2:] == "project-log":
                conf.project_log = terms[1]
            elif terms[0][2:] == "project-home":
                conf.project_home = terms[1]
                if terms[1][0:2] == "~/":
                    conf.project_home = os.environ["HOME"] + terms[1][1:]
            elif terms[0][2:] == "n-digit-run-id":
                conf.n_digit_run_id = int(terms[1])
            elif terms[0][2:] == "cell-size":  
                conf.cell_size = [float(terms[1]), float(terms[2]), float(terms[3])]
            elif terms[0][2:] == "cell-division":            
                conf.cell_size = [int(terms[1]), int(terms[2]), int(terms[3])]
            elif terms[0][2:] == "prerun-init":
                if len(terms) == 5:
                    if terms[1][0:2] == "~/":
                        terms[1] = os.environ["HOME"] + terms[1][1:]
                    n_begin = int(terms[2])
                    n_step = int(terms[3])
                    n_cal = int(terms[4])
                    for i in range(n_cal):
                        i_cal = n_begin + i*n_step
                        st = "%0"+str(conf.n_digit_run_id)+"d"
                        str_i_cal = st%i_cal
                        conf.prerun_init.append(terms[1]+str_i_cal+".pdb")
                else:
                    conf.prerun_init.append(terms[1])
            elif terms[0][2:] == "prerun-force-coef":
                if len(terms) == 4:
                    f_begin = float(terms[1])
                    f_step = float(terms[2])
                    n_stage = int(terms[3])
                    for i_stage in range(n_stage):
                        force_coef = f_begin + f_step * i_stage
                        conf.prerun_force_coef.append(force_coef)
                else:
                    conf.prerun_force_coef.append(terms[1])
            elif terms[0][2:] == "prerun-ene-range":
                conf.prerun_ene_range = (float(terms[1]), float(terms[2]))
            elif terms[0][2:] == "random-seed":
                conf.random_seed.append(int(terms[1]))
            elif terms[0][2:] == "prerun-md-inp-template":
                conf.prerun_md_inp_template = terms[1]
            elif terms[0][2:] == "prerun-ttp-inp-template":
                conf.prerun_ttp_inp_template = terms[1]
            elif terms[0][2:] == "inp-dist-restraint":
                conf.inp_dist_restraint = terms[1]
            elif terms[0][2:] == "inp-topology":
                conf.inp_topology = terms[1]
            elif terms[0][2:] == "inp-shake":
                conf.inp_shake = terms[1]
            elif terms[0][2:] == "inp-ttp":
                conf.inp_ttp = terms[1]
            elif terms[0][2:] == "inp-vert":
                conf.inp_vert = terms[1]
            elif terms[0][2:] == "out-mc-ene":
                conf.out_mc_ene = terms[1]
            elif terms[0][2:] == "cal-dir":
                conf.cal_dir.append(terms[1])
            elif terms[0][2:] == "dir-01-setting":
                conf.dir_01_setting = terms[1]
            elif terms[0][2:] == "dir-02-prerun":
                conf.dir_02_prerun = terms[1]
            elif terms[0][2:] == "dir-03-mcmd":
                conf.dir_03_mcmd = terms[1]
            elif terms[0][2:] == "mcmd-stages":
                conf.mcmd_stages[int(terms[1])] = int(terms[2])
            elif terms[0][2:] == "mcmd-inp-ttp":
                conf.mcmd_inp_ttp[int(terms[1])] = terms[2]
            elif terms[0][2:] == "mcmd-vs-range":
                conf.mcmd_vs_range[int(terms[1])] = terms[2]
            elif terms[0][2:] == "mcmd-md-inp-template":
                conf.mcmd_md_inp_template[int(terms[1])] = terms[2]
            #elif terms[0][2:] == "mcmd-init-prev-phase":
            #    conf.mcmd_init_prev_phase[int(terms[1])] = int(terms[2])
            #elif terms[0][2:] == "mcmd-init-prev-stages":
            #    conf.mcmd_init_prev_stages[int(terms[1])] = [int(x) for x in terms[2:]]
            elif terms[0][2:] == "mcmd-start-vs":
                phase = int(terms[1])
                conf.mcmd_start_vs[phase] = terms[2:]
            elif terms[0][2:] == "mcmd-init":
                phase = int(terms[1])
                if not phase in conf.mcmd_init:
                    conf.mcmd_init[phase] = []
                tmp = [x.split(":") for x in terms[2:]]
                for series, elem in enumerate(tmp):
                    if len(elem) < 2:
                        sys.stderr.write("syntax error: --mcmd-init; number of elements must be 2 or 3")
                        sys.exit(1)
                    elif len(elem) == 2:
                        conf.mcmd_init[phase].append((int(elem[0]), int(elem[1]), int(series+1)))
                    else:
                        conf.mcmd_init[phase].append((int(elem[0]), int(elem[1]), int(elem[2])))
            elif terms[0][2:] == "mcmd-vs-fitfunc":
                params = [terms[2]]
                params.extend([x for x in terms[3:]])
                if not int(terms[1]) in conf.mcmd_vs_fitfunc:
                    conf.mcmd_vs_fitfunc[int(terms[1])] = []
                conf.mcmd_vs_fitfunc[int(terms[1])].append(params)
            elif terms[0][2:] == "mcmd-n-steps-transition":
                conf.mcmd_n_steps_transition[int(terms[1])] = int(terms[2])
            elif terms[0][2:] == "mcmd-temperature":
                conf.mcmd_temperature[int(terms[1])] = float(terms[2])
            elif terms[0][2:] == "max-running":
                conf.max_running = int(terms[1])
            else:
                sys.stderr.write("kkmcmdconf.py McMDConfReader::read_conf()")
                sys.stderr.write("Unknown key word: " + terms[0] + "\n")
        conf.set_paths()
        return conf

class RangeInfo(kkkit.FileI):
    def __init__(self, fn):
        super(RangeInfo, self).__init__(fn)
        ## ene_range ... min, max
        self.ene_range = (0.0, 0.0)
        self.ene_width = 0.0
        ## vs_range[vs_id] = (min_ene, max_ene, tpro1, tpro2)
        self.vs_range = {}
        self.vs_params = {}
    def read_ttp_inp(self):
        self.open() 
        #line = self.f.readline().strip()
        #terms = line.split()
        
        #self.ene_range = (float(terms[0]), float(terms[1]))
        #self.ene_width = self.ene_range[1] - self.ene_range[0]
        
        line = self.readline_comment().strip()
        n_vs = int(line.split()[0])
        self.readline_comment()
        self.vs_range = {}
        #vsid = 0
        for vsid in range(1,n_vs+1):
            line = self.readline_comment().strip()
            terms = line.split()
            ene_min = float(terms[0])
            ene_max = float(terms[1])
            line = self.readline_comment().strip()
            p2 = float(terms[0])
            p1 = float(terms[1])
            self.vs_range[vsid] = (ene_min, ene_max, p1, p2)
        for vsid in range(1,n_vs+1):
            self.vs_params[vsid] = []
            order = int(self.readline_comment().strip())
            for i in range(order+3):
                self.vs_params[vsid].append(self.readline_comment().strip())
        self.close()
        return 0
    def read(self):
        self.open()

        line = self.readline_comment().strip()
        terms = line.split()
        self.ene_range = (float(terms[0]), float(terms[1]))
        self.ene_width = self.ene_range[1] - self.ene_range[0]
        
        line = self.readline_comment().strip()
        n_vs = int(line.split()[0])
        
        self.vs_range = {}
        for i in range(0,n_vs):
            line = self.readline_comment().strip()
            terms = line.split()
            vsid = int(terms[0])
            ratio_min = float(terms[1])
            ratio_max = float(terms[2])
            ene_min = self.ene_range[0] + self.ene_width * ratio_min
            ene_max = self.ene_range[0] + self.ene_width * ratio_max
            self.vs_range[vsid] = (ene_min, ene_max, float(terms[3]), float(terms[4]))
        self.close()
        return 0

    def get_near_state(self, ene):
        near_vsid = -1
        near_vsid_diff = 1e10
        for vsid, range in self.vs_range.items():
            #mid = (range[0] + (range[1] - range[0])*0.5)
            #diff = math.fabs(mid - ene)
            #if diff < near_vsid_diff:
            #    near_vsid_diff = diff
            #    near_vsid = vsid
            if range[0] <= ene:
                near_vsid = vsid
        return near_vsid

class McMDInitialDef(kkkit.FileI):
    def __init__(self, fn):
        super(McMDInitialDef, self).__init__(fn)
        self.init_vs = 0
        self.seed = 0
        return
    def read(self):
        self.open()     
        self.init_vs = int(self.readline_comment().strip())
        self.seed = int(self.readline_comment().strip())
        return
