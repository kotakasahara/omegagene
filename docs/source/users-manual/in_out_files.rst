========================
In/Out Files
========================

:Author: Kota Kasahara

------------------------------------
Input Files
------------------------------------

1. System configuration file (.cfg)
2. Simulation configuration file (.cfg)
3. Structure file (.pdb)
4. Initial coordinates and velocities file (.restart)
5. Topology file (.tpl)
6. Integrated binary (.cls)
7. SHAKE setting file (.shk)
8. V-McMD (or V-AUS) setting files (.inp, .vert)
9. V-AUS restart file (.dat)
10. Atom group definition file (.inp)
11. Distance restraint file (.inp)
12. Position restraint file (.inp)
13. VcMD parameter file (.inp)
14. VcMD initial state file (.inp)

System configuration file (.cfg)
===================================

The input file describing configurations about a simulation system.
Other some input files are specified in this file, and they are integrated into the integrated binary file (.cls) by using *mdinput_generator.py* program.
In the configuration file, a set of a key and value(s) is specified in each line.

* --fn-i-tpl          md.tpl
  * The file name of the topology file (.tpl).
* --fn-i-initial-pdb  md.pdb
  * The file name of the structure file (.pdb).
* --fn-i-restart      md.restart
  * The file name of the initial coordinates and velocities file (.restart)
* --cell-x            61.2425
* --cell-y            61.2425
* --cell-z            61.2425
  * The lengths of the periodic boundary cell in each axis in angestrome unit.
* --fn-i-shake        system.shk
  * The file name of the shake setting file.
* --fn-i-ttp-v-mcmd-inp      ttp_v_mcmd.inp
* --fn-i-ttp-v-mcmd-initial  start.vert
  * The file name of the V-McMD setting files.
* --fn-i-atom-group          atom_groups.inp
  * The file name of the atom group definition file.
* --fn-i-dist-restraint      dist_rest.inp
  * The file name of the distance restraint setting file.
* --fn-i-aus-restart         aus_rest.dat

Simulation configuration file (.cfg)
===================================

The input file describing configurations about simulation conditions.

Common configurations
------------------------------------

* --mode                  md
  * Only the keyword "md" is valid.
* --gpu-device-id         0
  * The device ID of GPGPU board to be used.
* --integrator            leapfrog-presto
  * Type of integrator
  * leapfrog-presto
    * The leap frog algorithm.
    * For NVE, NVT (rescaling), SHAKE
  * zhang
    * The integrator by Zhang [Zhang_1997]_
    * For the Hoover-Evans thermostat ("hoover-evans")
    * SHAKE cannot be applied.
  * langevin
    * Langevin integrator.
    * The thermostat must be set to "none"
    * "langevin-gamma" is required.
* --thermostat            scaling
  * none
    * NVE
  * scaling
    * Scaling velocities.
  * hoover-evans
    * Hoover-Evans. This can be applied for the integrator "zhang".
* --cutoff                12.0
  * The cutoff length for non-bonded potential (angestrome unit)
* --n-steps               10
  * The number of steps of the simulation.
* --time-step             2.0
  * The integration time step
* --electrostatic         zero-dipole
  * "zero-dipole"
    The Zero-dipole summation (ZD) method developed by Fukuda [Fukuda_2011]_ .
  * "zero-quadrupole"
  * "zero-octupole"
  * "zero-hexadecapole"
    The Zero-multipole summation (ZM) method developed by Fukuda [Fukuda_2014]_ .
  * "debye-huckel"
    For coarse-grained simulation.
    For GPU mode, only zero-dipole with --ele-alpha 0.0 is acceptable.
    To use "debye-huckel" with the HPS model with GPU, omegagene must be built with -DCELESTE_HPSGPU option.
* --ele-alpha             0.0
  * The dumping factor fo ZM method.
  * For GPU mode, only 0.0 is acceptable
* --debye-huckel-dielectric
  * A parameter for Debye-Huckel approx.
* --debye-huckel-ionic-strength  
  * A parameter for Debye-Huckel approx.
* --debye-huckel-temperature  300
  * A parameter for Debye-Huckel approx.
* --nonbond   lennard-jones
  * "lennard-jones"
    * Usual LJ potential.
  * "hydrophobicity-scale-lj"
    LJ potential modified with the hydrophobicity scale model.
    To use this otpion on a GPU, omegagene must be built with -DCELESTE_HPSGPU=1 option for cmake.
* --temperature           300
  * Temperature for the thermostat.
* --temperature-init
  * The initial temperature. Default value is the temperature specified by "--temperature" setting.
  * With this setting, "--heating-steps" should be set.
* --heating-steps
  * The temperature is linearly increased or decreased from the --temperature-init to --temperature during the steps specified in this setting.
* --berendsen tau 
  * Tau parameter for the Berendsen thermostat.
* --langevin-gamma
  * Gamma parameter for the Langevin integrator.
* --print-interval-log    1
  * Output interval for the log (the standard output)
* --print-interval-coord  1
  * Output interval for the trajectory file.
* --fn-o-coord            et1.trr
  * Output file name for the trajectory file.
* --format-o-coord        presto
  * The file format of the trajectory.
  * presto
* --fn-o-restart          md.restart
  * Output restart file name.
#* --fn-o-log              et1.log
#  *
#* --fn-o-energy           et1.ene
#  *
* --nsgrid-cutoff         13.0
  * The cutoff length for the neighbor search (angestrome unit).
* --nsgrid-update-intvl   10
  * The interval steps for execution of the neighbor search.
* --com-motion            cancel
  * Settings for canceling the center of mass
  * none
  * cancel
    Translation of the center of mass for some specified groups are cancelled.
    The groups should be specified in "--com-cancel-group-name"
* --com-cancel-group-name  grpA
  * The name of an atom group COM motions of which to be cancelled.
  * Multiple values can be specified.

Configuration for restraints
----------------------------------------------------

* --dist-restraint harmonic
  * Functions for distance restraints
  * none
  * harmonic
* --dist-restraint-weight
0  * Scaling coefficient for the distance restraints
* --position-restraint
  * none
  * harmonic
* --position-restraint-weight
  * Scaling coefficient for the

Configuration for the extended ensemble methods
----------------------------------------------------

* --extended-ensemble             v-mcmd
  * none
    Extended ensemble is not used
  * v-mcmd
    The V-McMD method [Higo et al. (2013)]_
  * v-aus
    The V-AUS method [Higo et al. (2015)]_
  * vcmd
    The VcMD emthod [Higo et al. (2017a)]_ [Higo et al. (2017b)]_ [Hayami et al. (2018)]_ [Hayami et al. (2019)]_
* --fn-o-vmcmd-log                ttp_v_mcmd.out
  * Output file name of a virtual-system trajecotry.
* --fn-o-extended-lambda            mule.ene
  * Output file name of a log of the lambda value.
* --print-interval-extended-lambda  1
  * Output interval for the log of the lambda value.
* --format-o-extended-lambda        ascii
  * File format of the log of the lambda value.
  * ascii
  * binary
#* --aus-type                        3
#  * Type of AUS reaction coordinate.
#  * 3
#    * The reaction coordinate is defined as the distance between centers of mass of the groups.
* --enhance-sigma
  * A parameter for the recovery force in V-AUS simulation.
  * A margin of the lambda value.
* --enhance-recovery-coef
  * Strength of the recovery force which works when the reaction coordinate is out of the predefined range.
* --fn-o-vcmd-start
  * Output file name describing the virtual state in the last step. This can be used for restarting the VcMD runs.
* --fn-o-vcmd-q-raw
  * Output file name for the populations of each virtual state.
* --fn-o-vcmd-q-raw_is
  * Output file name for the populations of each intersection of virtual states.
* --begin-count-q-raw
  * The step number which begans to count the populations.
* --vcmd-drift
  * 0 or 1.
  * 1 indicates the simulations with drifting in the virtual system without detailed-balance.
* --print-interval-group-com
  * Interval steps for output the reaction coordinate values.

.. [Zhang_1997]_ Zhang "Operator-splitting integrators for constant-temperature molecular dynamics" J. Chem. Phys. 106 (14) 1997


Initial coordinates and velocityies file (.restart)
-----------------------------------------------------

This file is compatible for myPresto/Psygene restart file.
For the first run, this file with random velocities can be generated by *presto_generate_velocities.py* scripts in the toolkit.
At the end of a simulation, final coordinates and velocities will be output at the file specified by *--fn-o-restart* option.

Topology file (.tpl)
------------------------------------

This file is compatible for myPresto/Psygene topology file.
It can be prepared by using myPresto/TopolgeneX program.

* --fn-i-tpl md.tpl

SHAKE setting file.
------------------------------------

This file is compatible for myPresto/Psygene SHAKE file.
It can be prepared by using SHAKEinp program.

It should be specified in the system configuration file:

* --fn-i-shake system.shk

V-McMD (or V-AUS) setting files (.inp, .vert)
-----------------------------------------------

These files are compatible for myPresto/Psygene files.

They should be specified in the system configuration file:

* --fn-i-ttp-v-mcmd-inp      ttp_v_mcmd.inp
* --fn-i-ttp-v-mcmd-initial  start.vert

ttp_v_mcmd.inp
~~~~~~~~~~~~~~~~~

This file describes the definition of virtual states and their bias functions.
::

  ;
  3                       # The number of virtual states (VSs).
  100                     # Interval steps for VS transitions.
  ;
  -37700.0   -35560.0     # The lower and upper bound of the 1st VS
  0.0  1.0                # Transition probability for lower and upper VS
  ;
  -36630.0   -34490.0     # The 2nd VS.
  1.0  1.0                # 
  ;
  -35560.0   -33420.0     # The 3rd VS.
  1.0  1.0                # 
  ;
  6                       # The bias function of 1st VS is 6-th order polynomial
  0.539566614289679E+06   # The parameters of the polynomial
  0.107078199910006E+03
  0.885349382820406E-02
  0.390387519416188E-06
  0.968208493740114E-11
  0.128058948306820E-15
  0.705682546575764E-21
  -0.225508543444448E-01  # The last two values are not used
  0.237465821192018E-01   # The last two valuse are not used
  6                       # The 2nd VS.
  0.153644514042787E+06
  0.313777316547802E+02
  0.266888476716589E-02
  0.121017618153941E-06
  0.308531417224255E-11
  0.419329550527650E-16
  0.237358207968186E-21
  -0.502410851186141E-01
  -0.443300654878840E-01
  6                       # The 2nd VS.
  -0.607303057237411E+06
  -0.129555612396211E+03
  -0.115143655222501E-01
  -0.545716541236728E-06
  -0.145465976138503E-10
  -0.206775621289763E-15
  -0.122453436139314E-20
  0.201226906548982E+00
  -0.190884220501175E-03
  300                     # Temperature

start_vert.inp
~~~~~~~~~~~~~~~~~

This file describes the initial virtual state and a random seed.::
  2         # The initial virtual state
  94265278  # The random seed


V-AUS restart file (.dat)
-----------------------------------------------

This file is required for continuing a finished V-AUS job. It should be specified in the system configuration file:

* --fn-i-aus-restart aus_rest.dat

At the end of a V-AUS job, this file is automatically generated at the file, specified by *--fn-o-aus-restart* option.

Atom group definition file (.inp)
-----------------------------------------------

This file defining groups of atoms. This informations are used for:

* Canceling the center of mass motion
* Defining enhaced groups for V-AUS simulations

In this ascii file, each line defines one atom group. The characters at the head of each line indicate the name of each group. Successive columns specify atoms in this group.

For example::

  group1 1 4 6-9
  group2 3 10-11 13

The group1 is composed of atoms 1, 4, 6, 7, 8, and 9. The group2 is composed of atoms 3, 10, 11, and 13.

The atom-ID are began from 1.

Distance restraint file (.inp)
-----------------------------------------------

This file defining the distance restraints between pairs of atoms. Each line between the keywords  "RDDSTC> LIST" and "RDDSTC> STOP" defines a restraint.
Distances between two atoms were restrained with the flat bottom potential. Each line specifies identities of two atoms, restraint forces for lower and upper bound, and lower and upller borders to apply restraint potential.

* Molecule-ID of the 1st atom
* Residue-ID of the 1st atom
* Residue name of the 1st atom
* Atom name of the 1st atom
* Molecule-ID of the 2nd atom
* Residue-ID of the 2nd atom
* Residue name of the 2nd atom
* Atom name of the 2nd atom
* The force coefficient for the lower bound.
* The force coefficient for the upper bound.
* Lower limit of the distance.
* Upper limit of the distance.

Position restraint file (.inp)
-----------------------------------------------

Each line specifies a position restraint for an atom.

* Atom-ID
* Equilibrium X coordinate
* Equilibrium Y coordinate
* Equilibrium Z coordinate
* Margin distance from the equilibrium position.
* Force coefficient.
* Type of restraint. "normal" or "z". Restraint with "z" omits the X and Y coordinates.

VcMD parameter file (.inp)
-----------------------------------------------

This file describes the definition of the virtual system and potentials for each states.

::
  100   ; The interval steps for virtual state transitions.
  2     ; The number of dimensions
  3   group1 group2  ; The number of virtual states, and names of two groups 
                     ; defining the reaction coordinate, for the 1st axis.
  3.0   5.0    ; The range of the reaction coordinates for the 1st state.
  4.0   6.0    ; That for the 2nd state.
  5.0   7.0    ; That for the 3rd state.
  4   group1 group3  ; The definition for the second axis.    
  3.0   5.0
  4.0   6.0
  5.0   7.0
  6.0   8.0
  END

VcMD initial state file (.inp)
-----------------------------------------------

This file specifies the initial virtual state.::
  2           ; The number of dimension
  5           ; Initial virtual state coordinate of the 1st dimension
  6           ; Initial virtual state coordinate of the 2nd dimension
  46582642    ; Random seed




------------------------------------
Output files
------------------------------------

1. Standard output
2. Trajectory file (.cod)
3. Restart file (.restart)
4. V-McMD (or V-AUS) lambda trajectory
5. V-McMD (or V-AUS) virtual-system trajectory
6. V-AUS restart file

Standard output
------------------

A simulation log will be output for the standard output.
Redirection to a file is recommended.

Trajectory file
------------------

The trajectory file format is compatible to myPresto/Psygene.

This file repeats the two pars: a header of the frame, and atomic coordinates at the frame.

For the header part:

* [4 bytes, INT] The size of header parts in bytes. Always "44".
* [4 bytes, INT] Step number
* [4 bytes, FLOAT] Time
* [4 bytes, FLOAT] CPU time
* [4 bytes, FLOAT] Total energy
* [4 bytes, FLOAT] Kinetic energy
* [4 bytes, FLOAT] Temperature
* [4 bytes, FLOAT] Potential energy
* [4 bytes, FLOAT] Always "0.0"
* [4 bytes, FLOAT] Always "0.0"
* [4 bytes, FLOAT] Always "0.0"
* [4 bytes, FLOAT] Always "0.0"
* [4 bytes, INT] The size of header parts in bytes. Always "44".

For the coordinates part:

* [4 bytes, INT] The size of this part. The nubmer of atom * 3 dimensions * 4 bytes.
* [FLOAT] X, Y, and Z coordinates of each atoms.
* [4 bytes, INT] The size of this part. The nubmer of atom * 3 dimensions * 4 bytes.

Restart file (.restart)
-------------------------

The restart file format is compatible to myPresto/Psygene.

This file is composed of the three pars: a header of the frame, atomic coordinates, and velocities.

For the header part:

* [4 bytes, INT] The length of the following text. Alwasy "80".
* [80 bytes, CHAR] Description of this simulation (version information)
* [4 bytes, INT] Alwasy "80".
* [4 bytes, INT] Alwasy "8".
* [4 bytes, INT] The number of atoms for coordinates.
* [4 bytes, INT] The number of atoms for velocities.
* [4 bytes, INT] Alwasy "8".
* [4 bytes, INT] Alwasy "36".
* [4 bytes, INT] Step number.
* [4 bytes, FLOAT] Time.
* [4 bytes, FLOAT] Total energy.
* [4 bytes, FLOAT] Kinetic energy.
* [4 bytes, FLOAT] Potential energy.
* [4 bytes, INT] Alwasy "36".

For the coordinates part:

* [4 bytes, INT] The size of this part in bytes. The number of atoms * 3 dimenstions * 8 bytes.
* [DOUBLE] X, Y, Z coordinates of each atom.
* [4 bytes, INT] The size of this part in bytes. The number of atoms * 3 dimenstions * 8 bytes.

For the velocity part:

* [4 bytes, INT] The size of this part in bytes. The number of atoms * 3 dimenstions * 8 bytes.
* [DOUBLE] X, Y, Z velocities of each atom.
* [4 bytes, INT] The size of this part in bytes. The number of atoms * 3 dimenstions * 8 bytes.


V-McMD (or V-AUS) lambda trajectory
--------------------------------------

For V-McMD or V-AUS simulations, the lambda values are written on this file.

When *--format-o-extended-lambda  ascii* is specified, a lambda value is recorded in each line of the ascii file.

When *--format-o-extended-lambda  binary*, is specified, the values are dumped as a binary file.

* [1-4 bytes] The magic number
* [5-8 bytes] The precision (4 or 8)
* [9-14 bytes] Always 1.
* [After that] The lambda values

V-McMD (or V-AUS) virtual-system trajectory
--------------------------------------------

The trajectory of virtual-system coordinates is written as a two-columns, tab separated table.

* The first column means the step number.
* The second column means the virtual-system coordinate.

For example, in the case that virtual-system transitions are done in every 1000 steps,::

  1     1
  1001  2
  3001  1
  4001  2
  5001  3

V-AUS restart file
---------------------

A binary file required for restarting V-AUS and VcMD simulations.
