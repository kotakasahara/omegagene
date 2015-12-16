SET(LIB_TYPE SHARED)
CUDA_ADD_LIBRARY(CelesteCUDA ${LIB_TYPE} cuda_common.cu cuda_setup.cu)
