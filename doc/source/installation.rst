========================
Installation
========================


:Author: Kota Kasahara

------------------------------------
Summary
------------------------------------

Celeste is written in C++ language. Build the binary by using *make* command.

There are three compile options:

1. For CPU
2. For CPU, without the neighbor search algorithm (not efficient)
3. For CPU with GPU (Cuda)

They can be switched by changeng the EXECUTABLE variable in Makefile:

::

  EXECUTABLE = celeste
  #EXECUTABLE = celeste_wons  #without neighbor search
  #EXECUTABLE = celeste_gpu   #CUDA

------------------------------------
Without neighbor search
------------------------------------

The neighbor search is effective for fast calculation, but the implementation is too complecated. For debugging or tests, the computation witouht the neighbor search is available.

::
 
  EXECUTABLE = celeste_wons

------------------------------------
Acceleration with GPU
------------------------------------

For GPU mode NVIDIA GPU Computer capability 3.5 or later is required.

The path to the CUDA library should be specified as CUDALOC variable in Makefile.

::

  EXECUTABLE = celeste_gpu
  CUDALOC = /opt/cuda/5.5/lib64

------------------------------------
Celeste Toolkit
------------------------------------

*CelesteTookit* is a series of scripts for pre- and post-processing for MD simulations with Celeste.
It requires *python2.7.x* and *numpy* library.
In this manual, is is supposed that the all files of *CelesteToolkit* is contained in the diretory specified in the environmental variable ${CELESTETK}. This path should be added in ${PYTHONPATH}.

bash::

  export CELESTETK="${HOME}/celeste/toolkit"
  export PYTHONPATH=${CELESTETK}:${PYTHONPATH}

csh::

  setenv CELESTETK "${HOME}/celeste/toolkit"
  setenv PYTHONPATH ${CELESTETK}:${PYTHONPATH}

