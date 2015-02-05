CU = nvcc

GEN_SM = -arch=sm_35
CUDALOC = -I/usr/apps.sp3/cuda/6.0/include
##TRST_INC = -I/opt/cuda/5.0/include

CUINC = $(TRST_INC) $(B40C_INC)
CUDBG =  --ptxas-options=-v -DNDEBUG # --host-compilation C -g 

## cuda 5.0
CUOPT = -D_THRUST_SORT -D_CUDA_RDCTN -D_CUDA_SYMBL -D_CUDA_EXTRCT #-D_CUDA_TMR -D_CUDA_INFO #-D_CUDA_CHKR2 -D_CUDA_SYMBL

#CUTHREADS = -DTHREADS64
#CUTHREADS = -DTHREADS128
CUTHREADS = -DTHREADS256
#CUTHREADS = -DTHREADS512

CUFLAGS = -O3  $(CUDBG) $(CUINC) $(GEN_SM) $(CUOPT) ${CUTHREADS}
#CUFLAGS = -O3  $(CUINC) $(GEN_SM) $(CUOPT) ${CUTHREADS}
CC = icpc 

CFLAGS = -O3 -DNDEBUG -G

OBJS = cuda_common.o \
  cuda_setup.o 

SRCS = cuda_common.cu \
  cuda_setup.cu 


all: $(OBJS)
	$(CULIBRY)
	ar cur libcelestecuda.a $(OBJS)

cuda_common.o : cuda_common.cu cuda_common.h
	$(CU) $(CUFLAGS) -c cuda_common.cu
cuda_setup.o : cuda_setup.cu cuda_setup.h
	$(CU) $(CUFLAGS) -c cuda_setup.cu

clean:
	rm -rf *.o *.a *.so
