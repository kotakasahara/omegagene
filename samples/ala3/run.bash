#!/bin/bash

OMEGABIN=${HOME}/local/og0/target/bin/omegagene_wons
OMEGATK=${HOME}/local/og0/toolkit

export PYTHONPATH=${OMEGATK}:${PYTHONPATH}

python2.7 ${OMEGATK}/presto_generate_velocities.py -i ala3.pdb --i-tpl ala3.tpl --i-shk ala3.shk -t 300 -o md.restart -s ${RANDOM} --mol --check > log01_genvel.txt

python2.7 ${OMEGATK}/mdinput_generator.py -i md.inp -o md.inp.cls -v v.0.52 > log02_inputgen.txt

${OMEGABIN} --cfg md.inp.run --inp md.inp.cls > log03_md.out
