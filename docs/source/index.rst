===============================================
myPresto/omegagene Documentation
===============================================

--------------------
About myPresto/omegagene
--------------------

**"myPresto/omegagene"** is the molecular dynamics (MD) simulation software.
myPresto/omegagene has several unique features, such as, the zero-multipole summation method and the virtual-system coupled sampler.

  The current version of omegagene can perform:

1. MD simulations on NVE ensemble
2. MD simulations on NVT ensemble (The velocity rescaling, or Hoover-Evans thermostat)
3. MD simulations on Virtual-system coupled samplers
4. Applying constraints with SHAKE algorithm
5. Calculations of electrostatic potentials based on the zero-dipole summation method.
6. Calculations of pairwise potentials on GPGPU (powered by CUDA 7.0, Computer Capability 3.5 is required).
7. Coarse-grained MD simulations with the hydrophobicity scale model and Debye-Huckel approximation.

------------------------------------
Developers
------------------------------------

* KASAHARA Kota, Ritsumeikan Univ., Japan
* TERAZAWA Hiroki, Ritsumeikan Univ., Japan
* ITAYA Hayato, Ritsumeikan Univ., Japan
* GOTO Satoshi, Ritsumeikan Univ., Japan
* TAKAHASHI Takuya, Ritsumeikan Univ., Japan
* MA Benson, Univ. Illinois, US
* GOTO Kota, Tokyo Tech, Japan
* BHASKAR Dasgupta, Osaka Univ., Japan
* HIGO Junichi, Univ. Hyogo, Japan
* FUKUDA Ikuo, Osaka Univ., Japan
* MASHIMO Tadaaki, AIST, Japan
* AKIYAMA Yutaka, Tokyo Tech, Japan
* NAKAMURA Haruki, Osaka Univ., Japan

==============================
Users' Manual
==============================

.. toctree::
    :numbered:
    :maxdepth: 2

    users-manual/run_md
    users-manual/in_out_files
    users-manual/analysis
    users-manual/samples

==============================
Build Manual
==============================

.. toctree::
    :numbered:
    :maxdepth: 2

    build-manual/installation
    build-manual/build_notes

==============================
Tutorial
==============================

.. toctree::
    :numbered:
    :maxdepth: 2

    users-manual/tutorial_cg 
    users-manual/tutorial_mcmd

