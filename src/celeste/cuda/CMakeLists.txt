# It appears that the library must be a SHARED library for linking to work on Linux (tested on Centos 7)
CUDA_ADD_LIBRARY(CelesteCUDA SHARED cuda_common.cu cuda_setup.cu)
