EXECUTABLE = celeste
#EXECUTABLE = celeste_gpu
#EXECUTABLE = celeste_mpigpu
#EXECUTABLE = celeste_mpi

### Environment variables for tsubame
CXX = icpc 
CC = icpc 
MPICXX = mpic++
#CXX = g++
#MPICXX = openmpic++
CUDALOC = /opt/cuda/5.5/lib64
# BOOST_INC = ${HOME}/include
# BOOST_LIB = ${HOME}/lib

ICCFLAGS = -ansi -O3 -w1 # -std=c++11
ICCLIB = 
ifeq ($(EXECUTABLE),celeste)
##	MODE_FLAG =
	MODE_FLAG = -DF_WO_NS
	CUDALIB = 
	CCEXE = $(CXX)
endif
ifeq ($(EXECUTABLE), celeste_gpu)
	MODE_FLAG = -DF_CUDA
	CUDALIB = -L$(CUDALOC) -Lcuda -lcudart -lcelestecuda
	CCEXE = $(CXX)
endif
ifeq ($(EXECUTABLE), celeste_mpigpu)
	MODE_FLAG = -DF_CUDA -DF_MPI
	CUDALIB = -L$(CUDALOC) -Lcuda -lcudart -lcelestecuda
#	MPILIB = -I${BOOST_INC} -L${BOOST_LIB} -lboost_mpi -lboost_serialization 
	MPILIB = 
	CCEXE = $(MPICXX)
endif

### Environment variables for DDBJ
#ICC_LIB_PATH=/opt/intel/composer_xe_2013_sp1.0.080/compiler/lib/intel64
#LD_LIBRARY_PATH=/opt/intel/composer_xe_2013_spl.0.080/compiler/lib/intel64:$LD_LIBRARY_PATH
#LD_PRELOAD=${ICC_LIB_PATH}/libirc.so:${ICC_LIB_PATH}/libsvml.so:${ICC_LIB_PATH}/libimf.so
#CC=/opt/intel/composer_xe_2013_sp1.0.080/bin/intel64/icpc
#ICCFLAGS = -ansi -O3 -w1 -I/home/kotakasahara/include -L/home/kotakasahara/lib  -L${ICC_LIB_PATH}
#ICCLIB = -L${ICC_LIB_PATH} -lirc -lsvml -limf --static 
### / For DDBJ
### For MBA2013
# CC=g++
# ICCFLAGS = -ansi -O3  
# ICCLIB = -lm
### For DDBJ
### / For DDBJ

OBJS = CelesteObject.o Celeste.o     Config.o \
       Read.o        Write.o  WriteTrr.o \
       MmSystem.o    PBC.o \
       RunMode.o     DynamicsMode.o \
       EnergyCalcObject.o     EnergyCalc.o \
       ForceFieldObject.o     ForceField.o \
       ElectrostaticObject.o  ZeroMultipoleSum.o \
       MiniCell.o SubBox.o \
       main.o
ifeq ($(EXECUTABLE), celeste_gpu)
OBJS = CelesteObject.o Celeste.o     Config.o \
       Read.o        Write.o  WriteTrr.o \
       MmSystem.o    PBC.o \
       RunMode.o     DynamicsMode.o \
       EnergyCalcObject.o     EnergyCalc.o \
       ForceFieldObject.o     ForceField.o \
       ElectrostaticObject.o  ZeroMultipoleSum.o \
       GpuDynamicsMode.o GpuEnergyCalc.o \
       MiniCell.o \
       main.o
endif
ifeq ($(EXECUTABLE), celeste_mpigpu)
OBJS = CelesteObject.o Celeste.o     Config.o \
       Read.o        Write.o  WriteTrr.o \
       MmSystem.o    PBC.o \
       RunMode.o     DynamicsMode.o \
       EnergyCalcObject.o     EnergyCalc.o \
       ForceFieldObject.o     ForceField.o \
       ElectrostaticObject.o  ZeroMultipoleSum.o \
       GpuDynamicsMode.o GpuEnergyCalc.o \
       MpiGpuDynamicsMode.o MpiGpuEnergyCalc.o \
       MiniCell.o \
       main.o
endif

#OBJS =$(SRCS:.c=.o)

.SUFFIXES: .o .cpp
# $(EXECUTABLE) : $(OBJS)
# 	rm -f $(EXECUTABLE)
# 	$(CCEXE) -o $(EXECUTABLE) $(OBJS) ${ICCLIB}

.cpp.o:
	$(CCEXE) -c $(ICCFLAGS) $(MODE_FLAG) $(MPILIB) $<

# $(OBJS): define.h

all: cpu gpu lnk

cpu: $(OBJS)
	$(CCEXE) -c $(ICCFLAGS) $(MODE_FLAG) $<

gpu:
ifeq ($(EXECUTABLE), celeste_gpu)
	$(MAKE) -C cuda/ 
endif
ifeq ($(EXECUTABLE), celeste_mpigpu)
	$(MAKE) -C cuda/ 
endif

lnk:
	$(CCEXE) -o $(EXECUTABLE) $(OBJS) $(CUDALIB) $(ICCLIB) $(MPILIB)

clean :
	-rm -f *~ *.o $(EXECUTABLE)
	$(MAKE) clean -C cuda/

