FIND_PACKAGE(CUDA REQUIRED)
SET( MINIMUM_CUDA_VERSION_REQUIRED 5.0 )
if (CUDA_VERSION VERSION_LESS MINIMUM_CUDA_VERSION_REQUIRED)
	MESSAGE( FATAL_ERROR "${BoldRed}Installed CUDA version is less than ${MINIMUM_CUDA_VERSION_REQUIRED}!${ColourReset}" )
endif()

#[[
	There is a known issue with CMake flag propagation that causes errors when compiling with nvcc:
		https://cmake.org/Bug/view.php?id=15240
		https://cmake.org/Bug/view.php?id=13450
		http://stackoverflow.com/questions/29121211/cuda-compilation-issue-with-cmake
	The following hack fixes this for now until a patch lands in CMake
]]
SET(CUDA_PROPAGATE_HOST_FLAGS OFF)

SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11;-arch=sm_35")
