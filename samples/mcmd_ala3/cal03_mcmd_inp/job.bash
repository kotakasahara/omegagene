#!/bin/bash

RANK=${1}
. /etc/profile.d/modules.sh
module load cuda python-extension/2.7
OMEGAROOT=/gs/hs0/hp170025/local/omegagene/og0_v052a
OMEGABIN=${OMEGAROOT}/omegagene_gpu
OMEGATK=${OMEGAROOT}/toolkit
LD_LIBRARY_PATH=$OMEGAROOT:${LD_LIBRARY_PATH}

python ${OMEGATK}/mdinput_generator.py -i md.inp -o md.inp.cls -v v.0.52 > log_inputgen.txt

${OMEGABIN} --cfg md.inp.run --inp md.inp.cls > md.out
