#!/bin/bash

path=`pwd`
echo $path
STAGE=${path##*/}
PREV_STAGE=$(( $STAGE - 1 ))
PREV_STAGE_V=$PREV_STAGE

OMEGABIN=${HOME}/local/omegagene/target/bin/omegagene_wons
OMEGATK=${HOME}/local/omegagene/toolkit

SERIES=1
N_SERIES=10

while [ $SERIES -le $N_SERIES ]
do 
    mkdir n${SERIES}
    cd n${SERIES}
    cp ../md.inp.run .
    cp ../md.inp .
    cp ../run_node.bash  .
    perl -pi -e "s/\#\{GPU_DEVICE_ID\}/${GPU}/g" md.inp.run
    perl -pi -e "s/\#\{PREV_STAGE\}/${PREV_STAGE}/g" md.inp
    perl -pi -e "s/\#\{PREV_STAGE_V\}/${PREV_STAGE_V}/g" md.inp
    perl -pi -e "s/\#\{SERIES\}/${SERIES}/g" md.inp
    perl -pi -e "s/\#\{PREV_STAGE\}/${PREV_STAGE}/g" run_node.bash
    perl -pi -e "s/\#\{SERIES\}/${SERIES}/g" run_node.bash
    
    echo "1" > start.virt
    echo "6" >> start.virt
    echo "${RANDOM}" >> start.virt
    TMP=${RANDOM}
    echo $TMP
    python2.7 ${OMEGATK}/presto_generate_velocities.py   -i ../../cg_q8/inp.pdb   --i-tpl ../../cg_q8/out.tpl   -t 10  -o md_0.restart   -s ${TMP}  --mol --check

    cd ..

    SERIES=$(( $SERIES + 1 ))
done

