========================================
Tutorial for Coarse-grained simulations
========================================
~~~~~~~~~~~~~~~~
Requirements
~~~~~~~~~~~~~~~~

The following environments are required to use myPresto/omegagene.

- Linux or macos system
- Python with numpy and scipy libraries
- cmake
- c++ 11 or later

In addition, use of bash environment is assumed in this document.
Some trivial changes are needed for csh users.

------------------------------
Download myPresto/omegagene
------------------------------


To download myPresto/omegagene, execute the following command on your terminal:

::

  $ git clone https://github.com/kotakasahara/omegagene.git


We also offer you to use our docker container.

::

  git clone https://github.com/terapizawa/myPresto-omegagene.git

------------------------------------
Installation
------------------------------------

Setting up a target build folder:

::

        # in ${PROJECT_ROOT} directory
	$ mkdir target
        $ cd target

The "${PROJECT_ROOT}" indicates the path to the directory of the downloaded omegagene repository in your environment.
Then, CMake evaluates all the external software dependencies for the selected build variant, and exit with errors if the dependency requirements are not met. CMake must be invoked on the `CMakeLists.txt` file in the **${PROJECT_ROOT}** directory.
Run the following command to configure for building the desired variant of myPresto/omegagene in ${PROJECT_ROOT}/target directory

::

   $ cmake -DCELESTE_WO_NS=1 ..

If you use our gpu-based myPreto/omegagene, use the *-DCELESTE_GPUHPS=1 option* as follows:

::

   $ cmake -DCELESTE_GPUHPS=1 ..

Then, make command compiles and builds the software.

::

   $ make

See also "Installation" and "Build manual" in this documentation.
The executable binary is generated in ${PROJECT_ROOT}/target/bin directory.

-----------------------------------------
MD simulations with coarse grained model
-----------------------------------------

~~~~~~~~~~~~~~~~~~~~~~~
Setting up input Files
~~~~~~~~~~~~~~~~~~~~~~~

1. Set paths

::

  $ export OMEGATK=${PROJECT_ROOT}/toolkit
  $ export OMEGABIN=${PROJECT_ROOT}/target/bin/omegagene_wons

You should change these settings depending on your environment.
If you use csh environment, *setenv* command should be used to set the environmental variables.

2. Seting up the directory to the simulation

::

  $ mkdir ${PATH_TO_YOUR_WORKING_DIRECTORY}

Set the variable ${PATH_TO_YOUR_TARGET_DIRECTORY} as the path to your working directory.
(e.g. /home/user/md/test )


Then, copy the sample files to the working directory.

::

  $ cp -r ${PROJECT_ROOT}/samples/cg_q8 ${PATH_TO_YOUR_WORKING_DIRECTORY}
  $ cd ${PATH_TO_YOUR_WORKING_DIRECTORY}/cg_q8

3. make a topology file

This sample system consists of two molecules of poly-Q octapeptides. See *inp.pdb* file.
To make the topology file, a .pdb file consisting of single molecule is needed.

By using an editor, copy the atoms in the first molecule in *inp.pdb* (the first eight lines) and paste to a new file, and save it as *single.pdb*.

inp.pdb::

  ATOM      1  CA  GLN A   1     -20.000 0.00000 -15.200  1.00  0.00      0.0X
  ATOM      2  CA  GLN A   2     -20.000 0.00000 -11.400  1.00  0.00      0.0X
  ATOM      3  CA  GLN A   3     -20.000 0.00000 -7.6000  1.00  0.00      0.0X
  ATOM      4  CA  GLN A   4     -20.000 0.00000 -3.8000  1.00  0.00      0.0X
  ATOM      5  CA  GLN A   5     -20.000 0.00000 0.00000  1.00  0.00      0.0X
  ATOM      6  CA  GLN A   6     -20.000 0.00000 3.80000  1.00  0.00      0.0X
  ATOM      7  CA  GLN A   7     -20.000 0.00000 7.60000  1.00  0.00      0.0X
  ATOM      8  CA  GLN A   8     -20.000 0.00000 11.4000  1.00  0.00      0.0X
  ATOM      9  CA  GLN B   1     0.00000 0.00000 -15.200  1.00  0.00      0.0X
  ATOM     10  CA  GLN B   2     0.00000 0.00000 -11.400  1.00  0.00      0.0X
  ATOM     11  CA  GLN B   3     0.00000 0.00000 -7.6000  1.00  0.00      0.0X
  ATOM     12  CA  GLN B   4     0.00000 0.00000 -3.8000  1.00  0.00      0.0X
  ATOM     13  CA  GLN B   5     0.00000 0.00000 0.00000  1.00  0.00      0.0X
  ATOM     14  CA  GLN B   6     0.00000 0.00000 3.80000  1.00  0.00      0.0X
  ATOM     15  CA  GLN B   7     0.00000 0.00000 7.60000  1.00  0.00      0.0X
  ATOM     16  CA  GLN B   8     0.00000 0.00000 11.4000  1.00  0.00      0.0X

single.pdb::

  ATOM      1  CA  GLN A   1     -20.000 0.00000 -15.200  1.00  0.00      0.0X
  ATOM      2  CA  GLN A   2     -20.000 0.00000 -11.400  1.00  0.00      0.0X
  ATOM      3  CA  GLN A   3     -20.000 0.00000 -7.6000  1.00  0.00      0.0X
  ATOM      4  CA  GLN A   4     -20.000 0.00000 -3.8000  1.00  0.00      0.0X
  ATOM      5  CA  GLN A   5     -20.000 0.00000 0.00000  1.00  0.00      0.0X
  ATOM      6  CA  GLN A   6     -20.000 0.00000 3.80000  1.00  0.00      0.0X
  ATOM      7  CA  GLN A   7     -20.000 0.00000 7.60000  1.00  0.00      0.0X
  ATOM      8  CA  GLN A   8     -20.000 0.00000 11.4000  1.00  0.00      0.0X
  
Then, conduct the following command to generate the topology file, *md.tpl*. 

::

  $ python2.7 ${OMEGATK}/gen_tpl.py --pdb single.pdb --param param.dat --tpl md.tpl --molname mol1

After that, change the number of molecules in the md.tpl.

md.tpl,::


  TPL> TITLE
  MD

  TPL> MOLECULES
  mol1                          2  ## <== HERE. This value must be "2".

  TPL> ATOMS
  mol1


4. Generate the initial atomic velocities

Execute the following command.

::

  python2.7 ${OMEGATK}/presto_generate_velocities.py   -i inp.pdb   --i-tpl md.tpl   -t 100   -o md.restart   -s ${RANDOM}  --mol --check

*-s* option indicates the random seed.

bash environment automatically generates a random number for the ${RANDOM} variable.
When it does not work, replace ${RANDOM} into an arbitral arbitral number.

5. make a cls file

::

  python2.7 ${OMEGATK}/mdinput_generator.py -i md.inp -o md.inp.cls -v v.0.52 > log_inputgen.txt

*md.inp.cls* file is the input file for myPresto/omegagene.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Set up your simulation conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simulation conditions and systems are configured by the following files.

- atom_groups.inp
- md.inp
- md.inp.run

**atom_groups.inp**

::

  mol1 1-8 # amino No for each molecules
  mol2 9-16
  all 1-16　# all amino acids in the input PDB file


**md.inp**

::

  --fn-i-tpl               md.tpl          # tpl file for the simulations
  --fn-i-initial-pdb       inp.pdb          # input PDB files
  --fn-i-restart           md.restart       # all initial positions for the input PDB file
  --cell-x                 50　           # maximum range of x axis
  --cell-y                 50             # maximum range of x axis
  --cell-z                 50             # maximum range of x axis
  --cell-center-x          0              # center position for x axis
  --cell-center-y          0              # center position for y axis
  --cell-center-z          0              # center position for z axis
  --fn-i-atom-groups       atom_groups.inp  # information for all amino acids and its molecules


**md.inp.run**

::

  --processor                   single        ;    # the numner of processors for conducting MD
  --gpu-device-id               0                  # GPU device ID for conducting MD
  --mode                        md            ;    # simulation mode
  --integrator                  langevin  ;        # the method of integration
  --langevin-gamma              1.0   ;      ;     # the parameter for friction coefficient
  --cutoff                      20.0          ;    # the cut-off distance in angstrome
  --n-steps                     2000       ;    # the simulation steps
  --time-step                   5            ;    # the integration time step (fs)
  --electrostatic               debye-huckel  ;    # the electrostatic interactions
  --debye-huckel-dielectric     85            ;    # the value of relative dielectric constant for debye-huckel equation
  --debye-huckel-temperature    300           ;    # the temperature for debye-huckel equation
  --debye-huckel-ionic-strength 0.00015       ;    # the ionic-strength value for debye-huckel equation
  --ele-alpha                   0             ;    # the alpha parameter for ZMM method
  --thermostat                  none               # options for using thermostat in MD
  --temperature                 300           ;    # simulation temperature
  --com-motion                  cancel      ;      # the option for canceling the motion of center-of-mass (COM)
  --com-cancel-group-name       all                # the name of predefined group for the canceling of COM motion
  --group-o-coord    all                           # the name of predefined group to output the trajectory
  --print-interval-log          100           ;    # the interval steps of making logs
  --print-interval-coord        100          ;   # the interval steps of making cods
  --fn-o-coord                  md.cod        ;    # the name of the trajectory output file
  --format-o-coord              presto             # the file format for the trajectory (only "presto" is supported currently)
  --fn-o-restart                md.restart         # the file contains the final conformation's positions
  --nsgrid-cutoff               23                 # the threshhold distance for neighbor molecules
  --nsgrid-update-intvl         5                 # the update interval for nsgrid
  --hydrophobicity-scale-epsiron 0.2               # a parameter for HPS model
  --nonbond hydrophobicity-scale-lj                # indication of using Lennerd-Jones potential
  --expected-num-density        0.1                #  a parameter to define the memory size. It is not recommended to change this default value.

~~~~~~~~~~~~~~~~~~~~~~
Execute omegagene
~~~~~~~~~~~~~~~~~~~~~~

To run an MD simulation using myPresto/omegagene, execute the following command. then please wait untill the job is done.

::

  ${OMEGABIN}  --cfg md.inp.run --inp md.inp.cls > md.out

Simulation log is given in *md.out*, and the trajectory is *md.cod*.

-----------------------------------
Visualize the resulant trajectory
-----------------------------------

The trajectory file md.cod is written in myPresto format. This can be converted into the Gromacs trajectory .trr format.

::

  python2.7 ${OMEGATK}/trajconv_presto_gro.py --i-pdb inp.pdb --i-crd md.cod -o md.trr --lx 50 --ly 50 --lz 50

The trajectory can be visualized by some standard visualizers (e.g., VMD and PyMOL).

In addition, the final snapshot in the restart file can be converted into .pdb file format

::

  python2.7 ${OMEGATK}/restart_to_pdb.py -i md.restart --i-pdb inp.pdb -o finalstep.pdb


========================================
VcMD
========================================

~~~~~~~~~~~~~~~~
Files
~~~~~~~~~~~~~~~~

The directory *${PROJECT_ROOT}/samples/cg_q8_vcmd* in the repository is used for the VcMD tutorial.
In addition, files in *${PROJECT_ROOT}/samples/cg_q8* are also used.
These two directories should be copied into your working directory.

In *${PROJECT_ROOT}/samples/cg_q8_vcmd* directory,  there are directories named as *1* and *2*. They correspond to the first, and second iterations.
In each directory, 10 parallel simulations will be carried out in directories named "n1", "n2", ..., "n10".

~~~~~~~~~~~~~~~~
The first iteration
~~~~~~~~~~~~~~~~

In *1* directory, execute the following scripts attached to the samples.
Note that modify these scripts to adjust the path ${OMEGABIN} and ${OMEGATK} to your omegagene binary and toolkit directory.

::

  $ bash c1_gen_inp.bash 

This script generates the directory *n1* to *n10*.

::

  $ bash c2_exe.bash 

This script sequentially execute simulations from *n1* to *n10*.

::

  $ bash c3_prep_next.bash

This script performes postprocessing for an iteraction.
- *vcmd_next.inp* and *vcmd_next_qraw.dat" are generated.
- *vcmd_next.inp* describes the canonical probability for each virtual state as an input for the next iteration.
- *vcmd_next_qraw.inp* describes the probability in the entire VcMD ensemble for each virtual state.

vcmd_next_qraw.dat::
  10
  1
  7 mol1 mol2
  3.0 5.0
  4.0 6.0
  5.0 9.0
  6.0 12.0
  9.0 17.0
  12.0 22.0
  17.0 30.0
  0 0.0
  2 0.00152933698463
  3 0.0271170564086
  4 0.0949144766083
  5 0.18308553295
  6 0.316262109242
  7 0.377091487806
  END

- The first line indicates the interval steps for virtual state transition trials.
- The second line indicates the number or reaction coordinates (N_RC). 
- The third line indicates the number of virtual state and name of atom groups to define the reaction coordinate.
- The following seven lines are the range of lambda for each virtual state.
- The following lines describes the sampled probability for each state. 

~~~~~~~~~~~~~~~~
The second iteration and further
~~~~~~~~~~~~~~~~

In the *2* directory, the same three script should be executed.
Afther that, make the directory *3* and copy files from *2* to *3*.
Then, repeat the same protocols.

~~~~~~~~~~~~~~~~
Production run
~~~~~~~~~~~~~~~~

After the convergence of the distribution in *vcmd_next_qraw.dat", execute the production run with the same manner.

~~~~~~~~~~~~~~~~~~~
Post-processing
~~~~~~~~~~~~~~~~~~~

Following script performs re-weighting of the ensemble.
Execute it in the directory of each run (e.g., *n1*, *n2*, ... directories).
Note that ${PREV_STAGE} indicates the number of previous iteration.

:: 

  $ python ${OMEGATK}/assign_traj_vs_lambda.py
     --i-qcano ../../${PREV_STAGE}/vcmd_next.inp     
     --i-cod md.cod             
    --interval-cod 1000   
     --i-lmb lambda.out    
    --interval-lmb 1        
     --i-vs ttp_vcmd.out   
     --interval-vs 10  
     -o prob.dat    



- *--i-qcano vcmd_next.inp* is the parameter file used for the VcMD simulation.
- *--i-cod md.cod* is the trajectory file obtained from the VcMD simulation.
- *--interval-cod 1000* should specifies the value same as the *--print-interval-coord* in *md.inp.run* file.
- *--i-lmb lambda.out* is the trajectory file for lambda values obtained from the VcMD simulation.
- *--interval-lmb 1* should specifies the value same as the *--print-interval-extended-lambda* in *md.inp.run* file.
- *--i-vs ttp_vcmd.out*  is the trajectory for the virtual state obtained from the VcMD simulation.
- *--interval-vs 10* is the interval for virtual state transitions.
- *-o prob.dat* is the output file describing probabilistic weight in the canonical ensemble for each snapshot.


