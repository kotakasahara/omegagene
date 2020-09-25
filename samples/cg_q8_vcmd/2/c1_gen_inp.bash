#!/bin/bash

path=`pwd`
echo $path
STAGE=${path##*/}
PREV_STAGE=$(( $STAGE - 1 ))
PREV_STAGE_V=$PREV_STAGE

#OMEGABIN=${HOME}/local/omegagene/target/bin/omegagene_wons
#OMEGATK=${HOME}/local/omegagene/toolkit

SERIES=1
N_SERIES=10

while [ $SERIES -le $N_SERIES ]
do 
    mkdir n${SERIES}
    cd n${SERIES}
    cp ../md.inp.run .
    cp ../md.inp .
    perl -pi -e "s/\#\{GPU_DEVICE_ID\}/${GPU}/g" md.inp.run
    perl -pi -e "s/\#\{PREV_STAGE\}/${PREV_STAGE}/g" md.inp
    perl -pi -e "s/\#\{PREV_STAGE_V\}/${PREV_STAGE_V}/g" md.inp
    perl -pi -e "s/\#\{SERIES\}/${SERIES}/g" md.inp

    cd ..

    SERIES=$(( $SERIES + 1 ))
done

