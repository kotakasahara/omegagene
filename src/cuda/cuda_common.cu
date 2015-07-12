#include "cuda_common.h"

extern "C" int cuda_set_device(int device_id){
  cudaSetDevice(device_id);
  return 0;
}
extern "C" int cuda_print_device_info(int myrank=0, bool verbose=false){
  cudaDeviceProp prop;
  int count;
  HANDLE_ERROR( cudaGetDeviceCount( &count) );
  int count_bak = count;
  for(int i=0; i< count_bak; i++){
    printf("GPGPU %s %d-%d\n",prop.name,myrank,i);
    HANDLE_ERROR( cudaGetDeviceProperties( &prop, i ) );
    if (!verbose) continue;
    printf("--- General Information for device %d / %d ---\n",i,count_bak-1 );
    printf("  Name: %s\n", prop.name);
    printf("  Compute capability: %d.%d\n", prop.major, prop.minor);
    printf("  Clock rate: %d\n", prop.clockRate);
    printf("  Device copy overlap: ");
    if (prop.deviceOverlap)
      printf("Enabled \n");
    else
      printf("Disabled\n");
    printf("  Kernel execition timeout : ");
    if(prop.kernelExecTimeoutEnabled)
      printf("Enabled\n");
    else
      printf("Disabled\n");
    printf("   --- Memory Information for device %d ---\n",i);
    printf("  Total global mem: %ld\n", prop.totalGlobalMem);
    printf("  Total constant Mem: %ld\n", prop.totalConstMem);
    printf("  Max mem pitch: %ld\n", prop.memPitch);
    printf("  Texture Alignment: %ld\n", prop.textureAlignment);
    printf("   --- MP Information forr device %d ---\n",i);
    printf("  Multiprocessor count: %d\n", prop.multiProcessorCount);
    printf("  Shared mem per mp: %ld\n", prop.sharedMemPerBlock);
    printf("  Registers per mp: %d\n", prop.regsPerBlock);
    printf("  Thread in warp: %d\n", prop.warpSize);
    printf("  Max threads per block: %d\n",prop.maxThreadsPerBlock);
    printf("  Max threads dimensions: (%d, %d, %d)\n",
	   prop.maxThreadsDim[0], prop.maxThreadsDim[1],
	   prop.maxThreadsDim[2]);
    printf("  Max grid dimensions: (%d, %d, %d)\n",
	   prop.maxGridSize[0], prop.maxGridSize[1],
	   prop.maxGridSize[2]);
    printf("\n");
  }
  
  return 0;
}

