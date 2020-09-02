#!/bin/bash

. /etc/profile.d/modules.sh
module load cuda python-extension/2.7
OMEGAROOT=/gs/hs0/hp170025/local/omegagene/og0_v052a
OMEGABIN=${OMEGAROOT}/omegagene_gpu
OMEGATK=${OMEGAROOT}/toolkit
LD_LIBRARY_PATH=$OMEGAROOT:${LD_LIBRARY_PATH}

for runid in `seq 1 10`; do
# runid=${1}
    runid_str=`printf "%05d" $runid`
    runid_prev=`expr $runid - 1`
    runid_prev_str=`printf "%05d" $runid_prev`

    #if [ $run_id1 -gt 1 ]
    #then
    cp run.system run${runid_str}.system
    perl -pi -e "s/\#\{RUN_ID_PREV\}/${runid_prev_str}/g" run${runid_str}.system
    cp run.inp run${runid_str}.inp
    perl -pi -e "s/\#\{RUN_ID\}/${runid_str}/g" run${runid_str}.inp

    python2.7 ${OMEGATK}/mdinput_generator.py -i run${runid_str}.system -o run${runid_str}.cls -v v.0.52 > log_inputgen_${runid_str}.txt

    ${OMEGABIN} --cfg run${runid_str}.inp --inp run${runid_str}.cls > run${runid_str}.out
done

wait
