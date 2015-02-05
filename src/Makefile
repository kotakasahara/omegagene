#EXECUTABLE = celeste
#EXECUTABLE = celeste_wons
EXECUTABLE = celeste_gpu

#EXECUTABLE = celeste_mpigpu
#EXECUTABLE = celeste_mpi

### Environment variables for tsubame
CXX = icpc 
CC = icpc 
MPICXX = mpic++
#CXX = g++
#MPICXX = openmpic++
#CUDALOC = /opt/cuda/5.5/lib64
CUDALOC = /usr/apps.sp3/cuda/6.0/lib64
# BOOST_INC = ${HOME}/include
# BOOST_LIB = ${HOME}/lib

ICCFLAGS = -ansi -O3 -w1 # -std=c++11
ICCLIB = 

MODE_FLAG =
CUDALIB = 
CCEXE = $(CXX)
ifeq ($(EXECUTABLE), celeste_wons)
	MODE_FLAG = -DF_WO_NS
endif
ifeq ($(EXECUTABLE), celeste_gpu)
	MODE_FLAG = -DF_CUDA
	CUDALIB = -L$(CUDALOC) -Lcuda -lcudart -lcelestecuda
endif
#ifeq ($(EXECUTABLE), celeste_mpigpu)
#	MODE_FLAG = -DF_CUDA -DF_MPI
#	CUDALIB = -L$(CUDALOC) -Lcuda -lcudart -lcelestecuda
#	MPILIB = -I${BOOST_INC} -L${BOOST_LIB} -lboost_mpi -lboost_serialization 
#	MPILIB = 
#	CCEXE = $(MPICXX)
#endif

OBJS = CelesteObject.o Celeste.o     Config.o \
       Read.o        Write.o  WriteTrr.o \
       MmSystem.o    PBC.o \
       RunMode.o     DynamicsMode.o \
       ForceFieldObject.o     ForceField.o \
       ElectrostaticObject.o  ZeroMultipoleSum.o \
       MiniCell.o SubBox.o \
       main.o

## endif

#OBJS =$(SRCS:.c=.o)

.SUFFIXES: .o .cpp

#$(EXECUTABLE) : $(OBJS)
#	rm -f $(EXECUTABLE)
# 	$(CCEXE) -o $(EXECUTABLE) $(OBJS) ${ICCLIB}

# $(OBJS): define.h

$(EXECUTABLE) : $(OBJS) cuda/*.cu cuda/*.h
ifeq ($(EXECUTABLE), celeste_gpu)
	$(MAKE) -C cuda/ 
endif
	$(CCEXE) -o $(EXECUTABLE) $(OBJS) $(CUDALIB) $(ICCLIB) $(MPILIB)

.cpp.o:
	$(CCEXE) -c $(ICCFLAGS) $(MODE_FLAG) $(MPILIB) $<

clean :
	-rm -f *~ *.o $(EXECUTABLE)
	$(MAKE) clean -C cuda/

