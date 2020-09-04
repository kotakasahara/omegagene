========================
Samples
========================

*samples* directory includes some sample simulation data.

- ala3
- cg_q8
- cg_q8_vcmd
- mcmd_ala3
- trpc

------
ala3
------

A sample for explicitly-solvated all-atom simulation.
This simulation model includes three molecules of capped-Ala peptide.
The simulation can be carried out by executing *run.bash* script.
Output files for a short simulation (100 steps) with the NVE ensemble also attached.
Users can test the built binary by using this data in terms of consistency with the attached output.
The potential energies in each step are recorded in *log03_md.txt* file.

-------
cg_q8
-------

A sample data for coarse-grained simulations.
See the "Tutorial for Coarse-grained simulations" in this documentation.

-------
cg_q8_vcmd
-------

A sample data for VcMD simulation with the coarse-grained model.
See the "Tutorial for Coarse-grained simulations" in this documentation.

-------
mcmd_ala3
-------

A sample data for the McMD simulation.
See the "Tutorial for multi-canonical MD (McMD) simulations" in this documentation.


-------
trpc
-------

A sample data for all-atom simulation. This was used for older versions.
