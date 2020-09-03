#!/usr/bin/python2.7

import sys
import re
import os
import random
import datetime
#sys.path.append(os.path.join(os.environ["HOME"],"local","kktools"))
#sys.path.append(os.path.join(os.environ["HOME"],"local","kktools","kkio"))
#sys.path.append(os.path.join(os.environ["HOME"],"local","kktools","mdtools"))
import kkkit
import kkmcmdconf as cf

from optparse import OptionParser
import subprocess

#DEBUG = True
DEBUG = False

def get_options():
    p = OptionParser()
    
    p.add_option('-i', dest='fn_conf',
                 help="file name for config file")
    p.add_option('-s', dest='queue_system',
                 type="choice", default="none",
                 choices=["tsubame", "uge", "tsubame-mpi", "tsubame3", "none"],
                 help="queue system")
    opts, args = p.parse_args()
    print "----------------------------"
    p.print_help()
    print "----------------------------"
    return opts, args
    
def _main():
    opts, args = get_options()
    confreader = cf.McMDConfReader(opts.fn_conf)
    conf = confreader.read_conf()
    if opts.queue_system == "tsubame-mpi" or opts.queue_system == "tsubame3":
        control = McMDControlTpMpi(conf, opts.queue_system)
    else:
        control = McMDControl(conf, opts.queue_system)
    control.main()
    return

class McMDControl:
    def __init__(self, conf, queue_system):
        ## conf ... instance of McMDConf
        self.conf = conf
        self.log = McMDLog(conf)
        self.jobcontrol = None
        self.set_queue_system(queue_system)
        return
    def set_queue_system(self, queue_system):
        if queue_system == "uge":
            self.jobcontrol = JobControlUGE()
        elif queue_system == "tsubame":
            self.jobcontrol = JobControlTsubame()
        else:
            self.jobcontrol = JobControl()
        return 
    def main(self):
        print "conf.set_job_status()"
        self.conf.set_job_status()
        print "log.parse_log()"
        self.conf.n_running = self.log.parse_log()
        print "n_running = " + str(self.conf.n_running)
        print "update_log()"
        self.update_log()
        print "prepare_run()"
        self.prepare_run()
        if self.jobcontrol:
            print "submit_run()"
            self.submit_run()
        return
    def update_log(self):
        jobs_info = self.jobcontrol.check_running_tasks()
        for phase, stage_series_state in self.conf.job_status.items():
            for stage, series_state in stage_series_state.items():
                for series, jobid_state in series_state.items():
                    if series in self.conf.ignore_serieses: continue

                    jobid = jobid_state[0]
                    state = jobid_state[1]
                    d_run = jobid_state[2]
                    if jobid in jobs_info and state != cf.McMDConf.State_run:
                        self.conf.job_status[phase][stage][series] = (jobid, cf.McMDConf.State_run, d_run)
                    elif (not jobid in jobs_info) and state == cf.McMDConf.State_run :
                        flg_end = self.check_end_of_job(phase, stage, series)
                        if flg_end:
                            self.log.append_log("[END]", jobid, phase, stage, series)
                            if self.check_end_of_phase(phase):
                                self.log.append_log("[PHASE_END]", -1, phase, -1, -1)
                        else:
                            msg = "McMDControl::check_end_of_job = False"
                            self.log.append_log("[ERROR]", jobid, phase, stage, series, msg)
        return
    def get_previous_state(self, phase, stage, series):
        if stage == 0:
            return (-1, cf.McMDConf.State_finish, "")
        return self.conf.job_status[phase][stage-1][series]
    def check_end_of_job(self, phase, stage, series):
        #print "dbg1"
        run_dir = self.conf.get_run_path(phase, stage, series)
        #fn_output = os.path.join(run_dir, self.conf.fn_stem+".pdb")
        fn_output = os.path.join(run_dir, self.conf.fn_stem+".restart")
        #print "DEBUG CHECK END OF JOB " + str(phase) + " " + str(stage) + " " + str(series)
        #print fn_output
        ret = os.path.exists(fn_output)
        return ret
    def check_end_of_phase(self, phase):
        for stage, series_state in self.conf.job_status[phase].items():
            for series, jobid_state in series_state.items():
                if series in self.conf.ignore_serieses: continue
                if jobid_state[1] != cf.McMDConf.State_finish:
                    return False
        return True
    def prepare_run(self):
        random.seed(self.conf.random_seed[1])
        prev_phase = 0
        for phase, stage_series_state in self.conf.job_status.items():
            for stage, series_state in stage_series_state.items():
                for series, jobid_state in series_state.items():
                    if series in self.conf.ignore_serieses: continue
                    jobid = jobid_state[0]
                    state = jobid_state[1]
                    d_run = jobid_state[2]
                    if DEBUG:
                        print "McMDControl::prepare_run() [" + str(phase) + "][" + str(stage) + "]["+str(series)+"] "
                        print str(jobid) + " - " + state + " - " + d_run
                        print self.get_previous_state(phase, stage, series)
                        
                    if (phase==1 or self.check_end_of_phase(phase-1)) and \
                            (state == cf.McMDConf.State_wait) and \
                            self.get_previous_state(phase, stage, series)[1] == cf.McMDConf.State_finish:
                        self.conf.prepare_run(phase, stage, series)

            prev_phase = phase
        return
    def submit_run(self):
        for phase, stage_series_state in self.conf.job_status.items():
            for stage, series_state in stage_series_state.items():
                for series, jobid_state in series_state.items():
                    if series in self.conf.ignore_serieses: continue
                    jobid = jobid_state[0]
                    state = jobid_state[1]
                    d_run = jobid_state[2]
                    name = "m" + str(phase) + "." + str(stage) + "." + str(series)
                    #print "RUNNING: " + str(self.conf.n_running)
                    if jobid_state[1] == cf.McMDConf.State_ready and \
                            self.conf.n_running < self.conf.max_running:
                        path = os.path.join(d_run, self.conf.job_script)
                        print "Submit " + str(phase) + " - " + str(stage) + " - " + str(series)
                        #print path
                        jobid = self.jobcontrol.submit_job(path, name, d_run, self.conf.submit_options)
                        self.log.append_log("[START]", jobid, phase, stage, series, note="")
                        self.conf.n_running += 1
                    self.conf.job_status[phase][stage][series] = (jobid, cf.McMDConf.State_run, d_run)
                    
        return

class McMDControlTpMpi(McMDControl):
    def __init__(self, conf, queue_system):
        McMDControl.__init__(self, conf, queue_system)
        return 
    def set_queue_system(self, queue_system):
        #if queue_system == "uge":
        #self.jobcontrol = JobControlUGE()
        #el
        if queue_system == "tsubame-mpi":
            self.jobcontrol = JobControlTsubameTpMpi()
        elif queue_system == "tsubame3":
            self.jobcontrol = JobControlTsubame3()
        else:
            self.jobcontrol = JobControl()
        return 
    def submit_run(self):
        run_info = []
        name = ""
        for phase, stage_series_state in self.conf.job_status.items():
            for stage, series_state in stage_series_state.items():
                for series, jobid_state in series_state.items():
                    if series in self.conf.ignore_serieses: continue
                    jobid = jobid_state[0]
                    state = jobid_state[1]
                    d_run = jobid_state[2]
                    
                    #print "RUNNING: " + str(self.conf.n_running)
                    
                    if jobid_state[1] == cf.McMDConf.State_ready and \
                            self.conf.n_running < self.conf.max_running:
                        #path = os.path.join(d_run, self.conf.job_script)
                        name = "m" + str(phase) + "." + str(stage) + "." + str(series)                        
                        run_info.append((phase, stage, series, jobid_state))
                        if len(run_info) == self.conf.n_mpi:
                            jobid = self.jobcontrol.submit_job(self.conf.job_script_mpi,
                                                               self.conf.job_script,
                                                               name,
                                                               self.conf.cal_dir_stages[phase][stage],
                                                               run_info,
                                                               len(run_info),
                                                               self.conf.submit_options)
                            for ri in run_info:
                                self.log.append_log("[START]", jobid, ri[0], ri[1],
                                                    ri[2], note="")
                                self.conf.n_running += 1
                            run_info = []
                            
                    self.conf.job_status[phase][stage][series] = (jobid, cf.McMDConf.State_run, d_run)
        if len(run_info) > 0:
            jobid = self.jobcontrol.submit_job(self.conf.job_script_mpi,
                                               self.conf.job_script,
                                               name,
                                               self.conf.cal_dir_stages[phase][stage],
                                               run_info,
                                               len(run_info),
                                               self.conf.submit_options)
            for ri in run_info:
                self.log.append_log("[START]", jobid, ri[0], ri[1],
                                    ri[2], note="")
                self.conf.n_running += 1
                    
        return
    
class McMDLog:
    def __init__(self, conf):
        self.conf = conf
        path = ""
        return
    def append_log(self, header, jobid, phase, stage, series, note=""):
        d = datetime.datetime.today()
        str_d = d.strftime("%Y-%m-%d-%H:%M:%S")
        terms = [header, str_d, str(jobid), str(phase), str(stage), str(series)]
        f = open(self.conf.project_log, "a")
        f.write("\t".join(terms))
        if note != "":
            f.write("\t"+note)
        f.write("\n")
        f.close()
    def parse_log(self):
        n_run = 0
        if not os.path.exists(self.conf.project_log):
            f = open(self.conf.project_log, "w")            
            f.close()
            return n_run
        
        f = open(self.conf.project_log)
        for line in f:
            terms = re.compile("\s+").split(line.strip())
            header = terms[0]
            if header[0] == ";" or header==[1]: continue
            if not header in ("[START]","[END]","[ERROR]","[WARNING]", "[PHASE_END]"):
                sys.stderr.write("Unknown header description: "+header)
                continue
            date = terms[1]
            jobid = terms[2]
            phase = int(terms[3])
            stage = int(terms[4])
            series = int(terms[5])
            #print "dbg2"
            d_run = self.conf.get_run_path(phase, stage, series)

            if header == "[START]":
                self.conf.job_status[phase][stage][series] = \
                    (jobid, cf.McMDConf.State_run, d_run)
                n_run += 1
            elif header == "[END]":
                self.conf.job_status[phase][stage][series] = \
                    (jobid, cf.McMDConf.State_finish, d_run)
                n_run -= 1
            elif header == "[ERROR]":
                self.conf.job_status[phase][stage][series] = \
                    (jobid, cf.McMDConf.State_error, d_run)
                n_run -= 1
        f.close()            
        return n_run

class JobControl:
    def __init__(self):
        self.cmd_stat = ""
        return
    def check_running_tasks(self):
        jobs_info = []
        if self.cmd_stat != "":
            cmd = [self.cmd_stat]
            #print cmd
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT,
                                 close_fds=True)
            jobs_info =  self.parse_stat(p.stdout)
        return jobs_info

    def parse_stat(self, handle):
        return 
    def submit_job(self, job_script, name, sub_options=[], options=[]):
        return 

class JobControlUGE(JobControl):
    def __init__(self):
        self.cmd_stat = "qstat"
        return
    def submit_job(self, job_script, name, d_run, sub_options=[]):
        #cmd = ["qsub", "-N", name, job_script, d_run]
        cmd = ["qsub", job_script]
        #print cmd
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             close_fds=True, cwd=d_run)
        jobid = -1
        for line in p.stdout:
            m = re.compile('Your job (\d+) \("(.*)"\) has been submitted').match(line)
            if m:
                jobid = int(m.group(1))
        return jobid
    def parse_stat(self, handle):
        """
        000000000011111111112222222222333333333344444444445555555555666666666677777777778888888888999999999900000000001111111111
        012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
        job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
        -----------------------------------------------------------------------------------------------------------------
        5232949 0.00000 QLOGIN     kotakasahara r     01/04/2014 16:16:41 login.q@t264i                      1        
        """
        jobs_info = {}
        handle.readline()
        handle.readline()
        for line in handle:
            terms = re.compile("\s+").split(line.strip())
            job_info = {}
            job_id = terms[0]
            job_info["prior"] = float(terms[1])
            job_info["name"] = terms[2]
            job_info["user"] = terms[3]
            job_info["state"] = terms[4]
            job_info["start_day"] = terms[5]
            job_info["start_time"] = terms[6]
            #job_info["queue"] = terms[7]
            #job_info["slot"] = terms[8]
            jobs_info[job_id] = job_info
            if DEBUG: print job_info
        return jobs_info
    

class JobControlTsubame(JobControl):
    def __init__(self):
        self.cmd_stat = "t2stat"
        return
    def submit_job(self, job_script, name, d_run, sub_options=[]):
        #t2sub -q G -l walltime=12:00:00 -W group_list=t2g-hp130061 -l select=3:mpiprocs=3:gpus=3:mem=2gb
        cmd = ["t2sub", "-N", name]
        if len(sub_options) > 0:
            cmd.extend(sub_options)
            #cmd.extend(["-q","G","-et","1","-l","walltime=24:00:00"])
            #cmd.extend(["-W","group_list=t2g-hp130061"])
        cmd.extend(["-v","CAL_DIR="+d_run])
        #n_nodes = 3
        #n_procs_node = 3
        #n_gpus_node = 3
        #mem_node = 3
        #cmd.extend(["-l","select="+str(n_nodes)+":mpiprocs="+str(n_procs_node)+":gpus="+str(n_gpus_node)+":mem="+str(mem_node)+"gb"])
        cmd.extend([job_script])

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             close_fds=True,
                             cwd=d_run)


        jobid = -1
        submit_log = ""
        submit_log += p.stdout.readline()
        submit_log += p.stdout.readline()
        submit_log += p.stdout.readline()
        line = p.stdout.readline()
        submit_log += line
        m = re.compile("(\d+\.t2zpbs\d+)").search(line.strip())
        for line in p.stdout: submit_log += line
        if m:
            jobid = m.group(1)

        print " ".join(cmd)

        return jobid
    def parse_stat(self, handle):
        """
        """
        jobs_info = {}
        for line in handle:
            terms = re.compile("\s+").split(line.strip())
            m = re.compile("^\d+\..+$").match(terms[0])
            if len(terms) == 6 and m:
                job_info = {}
                job_id = terms[0]
                job_info["prior"] = 0.0
                job_info["name"] = terms[1]
                job_info["user"] = terms[2]
                job_info["state"] = terms[4]
                job_info["start_day"] = "-"
                job_info["start_time"] = terms[3]
                job_info["queue"] = terms[5]
                job_info["slot"] = 0
                jobs_info[job_id] = job_info
                if DEBUG: print job_info
        return jobs_info
    
class JobControlTsubameTpMpi(JobControl):
    def __init__(self):
        self.cmd_stat = "t2stat"
        return
    def submit_job(self, job_script_mpi, job_script, name, d_run,
                   run_info, n_slots, sub_options=[]):
        #t2sub -q G -l walltime=12:00:00 -W group_list=t2g-hp130061 -l select=3:mpiprocs=3:gpus=3:mem=2gb
        cmd = ["t2sub", "-N", name]
        if len(sub_options) > 0:
            cmd.extend(sub_options)
            #cmd.extend(["-q","G","-et","1","-l","walltime=24:00:00"])
            #cmd.extend(["-W","group_list=t2g-hp130061"])
        cmd.extend(["-v"])
        val = []
        val.append("CAL_DIR="+d_run)
        val.append("NSLOTS="+str(n_slots))
        val.append("JOB_SCRIPT="+job_script)
        for rank in range(len(run_info)):
            val.append("R"+str(rank+1)+"="+run_info[rank][3][2])
        #n_nodes = 3
        #n_procs_node = 3
        #n_gpus_node = 3
        #mem_node = 3
        #cmd.extend(["-l","select="+str(n_nodes)+":mpiprocs="+str(n_procs_node)+":gpus="+str(n_gpus_node)+":mem="+str(mem_node)+"gb"])
        cmd.append(",".join(val))
        cmd.extend([job_script_mpi])

        #print cmd
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             close_fds=True,
                             cwd=d_run)

        jobid = -1
        submit_log = ""
        submit_log += p.stdout.readline()
        submit_log += p.stdout.readline()
        submit_log += p.stdout.readline()
        line = p.stdout.readline()
        submit_log += line
        m = re.compile("(\d+\.t2zpbs\d+)").search(line.strip())
        for line in p.stdout: submit_log += line
        if m:
            jobid = m.group(1)

        print " ".join(cmd)

        return jobid
    def parse_stat(self, handle):
        """
        """
        jobs_info = {}
        for line in handle:
            terms = re.compile("\s+").split(line.strip())
            m = re.compile("^\d+\..+$").match(terms[0])
            if len(terms) == 6 and m:
                job_info = {}
                job_id = terms[0]
                job_info["prior"] = 0.0
                job_info["name"] = terms[1]
                job_info["user"] = terms[2]
                job_info["state"] = terms[4]
                job_info["start_day"] = "-"
                job_info["start_time"] = terms[3]
                job_info["queue"] = terms[5]
                job_info["slot"] = 0
                jobs_info[job_id] = job_info
                if DEBUG: print job_info
        return jobs_info

class JobControlTsubame3(JobControl):
    def __init__(self):
        self.cmd_stat = "qstat"
        return
    def submit_job(self, job_script_mpi, job_script, name, d_run,
                   run_info, n_slots, sub_options=[]):
        #t2sub -q G -l walltime=12:00:00 -W group_list=t2g-hp130061 -l select=3:mpiprocs=3:gpus=3:mem=2gb
        cmd = ["qsub", "-N", name]
        if len(sub_options) > 0:
            cmd.extend(sub_options)
        cmd.extend(["-v"])
        val = []
        val.append("CAL_DIR="+d_run)
        val.append("NSLOTS="+str(n_slots))
        val.append("JOB_SCRIPT_SUB="+job_script)
        for rank in range(len(run_info)):
            val.append("R"+str(rank+1)+"="+run_info[rank][3][2])
        cmd.append(",".join(val))
        cmd.extend([job_script_mpi])

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             close_fds=True,
                             cwd=d_run)

        jobid = -1
        for line in p.stdout:
            m = re.compile('Your job (\d+) \("(.*)"\) has been submitted').match(line)
            if m:
                jobid = int(m.group(1))
        print " ".join(cmd)
        return jobid
    def parse_stat(self, handle):
        """
        000000000011111111112222222222333333333344444444445555555555666666666677777777778888888888999999999900000000001111111111
        012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
        job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
        -----------------------------------------------------------------------------------------------------------------
        5232949 0.00000 QLOGIN     kotakasahara r     01/04/2014 16:16:41 login.q@t264i                      1        
        """
        jobs_info = {}
        handle.readline()
        handle.readline()
        for line in handle:
            terms = re.compile("\s+").split(line.strip())
            job_info = {}
            job_id = terms[0]
            job_info["prior"] = float(terms[1])
            job_info["name"] = terms[2]
            job_info["user"] = terms[3]
            job_info["state"] = terms[4]
            job_info["start_day"] = terms[5]
            job_info["start_time"] = terms[6]
            #job_info["queue"] = terms[7]
            #job_info["slot"] = terms[8]
            jobs_info[job_id] = job_info
            if DEBUG: print job_info
        return jobs_info

if __name__ == "__main__":
    _main()

