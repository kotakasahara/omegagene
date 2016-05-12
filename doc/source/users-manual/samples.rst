========================
Samples
========================

:Author: Kota Kasahara

--------------------------
Sample files
--------------------------

As a sample, the files for the simulation of a mini-protein, Trp-cage, are included.

:: 

  python ${CELESTETK}/mdinput_generator.py -i system.cfg -o trpc.cls -v v.0.36.f
  omegagene --inp trpc.cls --cfg md.cfg
