#!/bin/bash

###t2sub -q G -l walltime=24:00:00 -W group_list=t2g-hp130061 -l select=3:mpiprocs=3:gpus=3:mem=2gb  ./job_prerun_tbm.bash

NSLOTS=8
cd ${CAL_DIR}
echo ${CAL_DIR}

#export PATH=/usr/apps/openmpi/1.4.2/pgi/bin:$PATH
#export LD_LIBRARY_PATH=/usr/apps/openmpi/1.4.2/pgi/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/opt/cuda/5.5/lib64:$LD_LIBRARY_PATH

pbs_inp=md.inp
pbs_out=md.out
pbs_exe=$HOME/bin/psygene-G_v0926c6

echo "mpirun -np ${NSLOTS} -machinefile ${TMPDIR}/machines ${pbs_exe} < ${pbs_inp} > ${pbs_out}"
mpirun -n $NSLOTS -machinefile ${PBS_NODEFILE} ${pbs_exe} < ${pbs_inp} > ${pbs_out}
