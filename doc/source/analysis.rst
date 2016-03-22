========================
Analysis
========================

:Author: Kota Kasahara

--------------------------
Log
--------------------------

Calculation log is printed to the standard output.
The total energy of the system should be constant, in the microcanonical ensemble.
When the total energy drifts, you should take care of the following points:

1. Increase --nsgrid-cutoff 
2. Decrease --nsgrid-update-intvl
3. Increase --cutoff

::
  
  Step:        0    Time:     0.0000
  Total:     -7.4102976562e+04
  Potential: -8.7917343750e+04    Kinetic:  1.3814368164e+04
  Bond:      2.2974748345e+02    Angle:    9.5691961023e+01
  Torsion:   2.9869403994e+02    Improper: 2.9828235618e+00
  14-VDW:    1.0501907707e+02    14-Ele:   1.8181411150e+03
  VDW:       1.5527498751e+04    Ele:      -1.0599511706e+05
  
  Step:   100000    Time: 50000.0000
  Total:     -7.4026648438e+04
  Potential: -8.7759882812e+04    Kinetic:  1.3733237305e+04
  Bond:      9.4730405123e+03    Angle:    2.4954034353e+02
  Torsion:   3.2906814547e+02    Improper: 1.2255384929e+01
  14-VDW:    1.2528642296e+02    14-Ele:   1.8064233115e+03
  VDW:       1.7495527442e+04    Ele:      -1.1725102074e+05


--------------------------
Trajectory
--------------------------

*omegagene* output the trajectory in the format compatible with myPresto/psygene-G.
If you want to convert the trajectory file into Gromacs .trr format, the script to do it is included in *omega-toolkit*.

python ${OMEGATK}/convert_trajectory_presto.py --i-pdb initia.pdb --i-crd traj.crd -o traj.trr

