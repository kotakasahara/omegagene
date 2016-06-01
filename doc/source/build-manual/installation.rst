========================
Installation
========================

:Author: Kota Kasahara

.. contents::

Celeste is written in C++, and its build system utilizes CMake.  The following is a non-exhaustive description for building Celeste.  Currently there are three compile options of Celeste:

1. CPU
2. CPU without the neighbor search algorithm (inefficient)
3. CPU with GPU (CUDA)

For platform-specific details on building Celeste, please refer to the :doc:`build_notes`.

------------------------------------
Software Requirements
------------------------------------

* CMake 3.4+
* For the GPU version of Celeste, CUDA 7.0+ is required for C++11 support.
* A C++11 compiler:
    * GCC: 4.8+
    * Clang: 3.6+
    * AppleClang: 5.0+
    * Intel: 15.0+ (minimal version required by CUDA 7.0+)
* OpenMP 3.1+
* Python 2.7.x
* numpy

------------------------------------
Building Celeste - Standard Version
------------------------------------

1. Set up a target build folder:

.. code-block:: bash

        # in <PROJECT_ROOT> directory
        localhost:celeste local$ mkdir target
        localhost:celeste local$ cd target

2. Configure the build.  CMake will determine all the external software dependencies for the selected build variant, and exit with errors if the dependency requirements are not met.  CMake must be invoked on the ``CMakeLists.txt`` file in the ``<PROJECT_ROOT>`` directory:

.. code-block:: bash

        # in <PROJECT_ROOT>/target directory
        localhost:target local$ cmake ..

3. Build the software:

.. code-block:: bash

        # The verbose flag is optional
        localhost:target local$ make VERBOSE=1

------------------------------------------------------------------------
Building Celeste Without Neighbor Search Routines
------------------------------------------------------------------------

While neighbor search is effective for fast calculations, the implementation is complicated and may be difficult to debug MD runs.  For this reason, a version of Celeste without the neighbor search routines can be built for debugging or testing.

To build this version of Celeste, simply run the following command instead when configuring the build (Step 2):

.. code-block:: bash

        localhost:target local$ cmake -DCELESTE_WO_NS=1 ..

The compiled executable will be named ``celeste_wons``.

------------------------------------------------------------------------
Building Celeste With GPU Acceleration
------------------------------------------------------------------------

For building this version of Celeste, CUDA 7.0+ is required.  For running the binary, an NVIDIA GPU with Compute Capability >= 3.5 or later is required.

To build this version of Celeste, simply run the following command instead when configuring the build (Step 2):

.. code-block:: bash

        localhost:target local$ cmake -DCELESTE_GPU=1 ..

CMake will automatically determine the default installation paths for the CUDA libraries and ``nvcc``.  Please refer to the Build Notes if you have installed CUDA to a custom filesystem path.

The compiled executable will be named ``celeste_gpu``.

------------------------------------
Celeste Toolkit
------------------------------------

*CelesteTookit* is a library of pre- and post-processing scripts for MD simulations to be used with Celeste.  It requires Python 2.7.x and the ``numpy`` library.

This manual assumes that the CelesteToolkit directory specified in the environmental variable ``${CELESTETK}``. This path should be added in ``${PYTHONPATH}``:

.. code-block:: bash

    export CELESTETK="${HOME}/celeste/toolkit"
    export PYTHONPATH=${CELESTETK}:${PYTHONPATH}

.. code-block:: csh

    setenv CELESTETK "${HOME}/celeste/toolkit"
    setenv PYTHONPATH ${CELESTETK}:${PYTHONPATH}
