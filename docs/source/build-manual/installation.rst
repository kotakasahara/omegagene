========================
Installation
========================

.. contents::

omegagene is written in C++, and its build system utilizes CMake.  The following is a non-exhaustive description for building omegagene. Currently there are three compile options of omegagene:

1. CPU (without the compile option)
2. CPU without the neighbor search algorithm (-DCELESTE_WO_NS=1) 
3. CPU with GPU (CUDA) (-DCELESTE_GPU=1)
4. CPU with GPU (CUDA) for coarse-grained simulations (-DCELESTE_GPUHPS=1)

For platform-specific details on building omegagene, please refer to the :doc:`build_notes`.


+-------------------+--------+-----------------+------------+----------------+
| option            | GPU    | Neighbor-search | all-atom   | coarse-grained | 
+===================+========+=================+============+================+
| (none)            | -      | enable          | enable     | -              |
| -DCELESTE_WO_NS   | -      | -               | enable     | enable         |
| -DCELESTE_GPU     | enable | enable          | enable     | -              |
| -DCELESTE_GPUHPS  | enable | enable          | -          | enable         |
+-------------------+--------+-----------------+------------+----------------+


------------------------------------
Software Requirements
------------------------------------

* CMake 3.4+
* For the GPU version of omegagene, CUDA 7.0+ is required for C++11 support.
* A C++11 compiler:
    * GCC: 4.8+
    * Clang: 3.6+
    * AppleClang: 5.0+
    * Intel: 15.0+ (minimal version required by CUDA 7.0+)
* OpenMP 3.1+
* Python 2.7.x
* numpy

------------------------------------
Building omegagene - Standard Version
------------------------------------

1. Set up a target build folder:

.. code-block:: bash

        # in <PROJECT_ROOT> directory
        localhost$ mkdir target
        localhost$ cd target

2. Configure the build.  CMake will determine all the external software dependencies for the selected build variant, and exit with errors if the dependency requirements are not met.  CMake must be invoked on the ``CMakeLists.txt`` file in the ``<PROJECT_ROOT>`` directory:

.. code-block:: bash

        # in <PROJECT_ROOT>/target directory
        localhost:target local$ cmake ..

3. Build the software:

.. code-block:: bash

        # The verbose flag is optional
        localhost:target local$ make VERBOSE=1

------------------------------------------------------------------------
Building omegagene Without Neighbor Search Routines
------------------------------------------------------------------------

While neighbor search is effective for fast calculations, the implementation is complicated and may be difficult to debug MD runs.  For this reason, a version of omegagene without the neighbor search routines can be built for debugging or testing.

To build this version of omegagene, simply run the following command instead when configuring the build (Step 2):

.. code-block:: bash

        localhost:target local$ cmake -DCELESTE_WO_NS=1 ..

The compiled executable will be named ``omegagene_wons``.

------------------------------------------------------------------------
Building omegagene With GPU Acceleration
------------------------------------------------------------------------

For building this version of omegagene, CUDA 7.0+ is required.  For running the binary, an NVIDIA GPU with Compute Capability >= 3.5 or later is required.

To build this version of omegagene, simply run the following command instead when configuring the build (Step 2):

.. code-block:: bash

        localhost:target local$ cmake -DCELESTE_GPU=1 ..

CMake will automatically determine the default installation paths for the CUDA libraries and ``nvcc``.  Please refer to the Build Notes if you have installed CUDA to a custom filesystem path.

The compiled executable will be named ``omegagene_gpu``.

------------------------------------------------------------------------
Building omegagene for a coarse-grained method
------------------------------------------------------------------------

For the coase-grained simulation on GPU, *-DCELESTE_GPUHPS=1* option is required.

.. code-block:: bash

        localhost:target local$ cmake -DCELESTE_GPUHPS=1 ..

For CPU, the coarse-grained simulation can be performed by the binary with *-DCELESTE_WO_NS=1* option.

------------------------------------
omegagene Toolkit
------------------------------------

*omegagene tookit* is a library of pre- and post-processing scripts for MD simulations to be used with omegagene  It requires Python 2.7.x and the ``numpy`` library.

This manual assumes that the omegagene toolkit directory specified in the environmental variable ``${OMEGATK}``. This path should be added in ``${PYTHONPATH}``:

.. code-block:: bash

    export OMEGATK="${HOME}/omegagene/toolkit"
    export PYTHONPATH=${OMEGATK}:${PYTHONPATH}

.. code-block:: csh

    setenv OMEGATK "${HOME}/omegagene/toolkit"
    setenv PYTHONPATH ${OMEGATK}:${PYTHONPATH}
