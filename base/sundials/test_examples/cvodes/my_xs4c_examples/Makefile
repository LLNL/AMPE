# File Makefile
#--------------------------------------------------------------------------
# Makefile for the CVODES examples using the serial N_Vector implementation.
# Generated automatically from Makefile.in by configure.
#--------------------------------------------------------------------------

# Shell
SHELL = /bin/sh

# SUNDIALS instalation directory
prefix = /home/radu/CODES/sundials
exec_prefix = ${prefix}
includedir = ${prefix}/include
libdir = ${exec_prefix}/lib

# SUNDIALS library and include directories
SUNDIALSINCLUDE = -I${includedir}
SUNDIALSLIBDIRS = -L${libdir}
SUNDIALSLIBS    = -lcvodes -lshared
NVECSERLIB      = -lnvecserial
NVECPARLIB      = -lnvecparallel

# MPI library and include directories
MPIINCLUDE =  -I/usr/apps/mpich/1.2.4/include
MPILIBDIRS =  -L/usr/apps/mpich/1.2.4/lib
MPILIBS    =  -lmpich
MPIFLAGS   =  

# Include
INCDIRS = ${SUNDIALSINCLUDE}

# Libraries
LIBDIRS = ${SUNDIALSLIBDIRS} ${MPILIBDIRS}
LIBS    = ${SUNDIALSLIBS} ${MPILIBS} -lm

# Compiler and compiler options
CPP = mpiCC
CC     = mpicc
CFLAGS = -g -Wall -ffloat-store

#--------------------------------------------------------------------------------
# Make rules
#--------------------------------------------------------------------------------

all: cvfkxCS cvfdxCS pvfkxCS

#-------

cvfkxCS: cvfkxCS.o cvfkxCSf.o cvfkxCSf.c_Cplx_.o
	$(CPP) ${CFLAGS} -o cvfkxCS cvfkxCS.o cvfkxCSf.o cvfkxCSf.c_Cplx_.o ${LIBDIRS} ${LIBS} ${NVECSERLIB}

cvfdxCS: cvfdxCS.o cvfdxCSf.o cvfdxCSf.c_Cplx_.o
	$(CPP) ${CFLAGS} -o cvfdxCS cvfdxCS.o cvfdxCSf.o cvfdxCSf.c_Cplx_.o ${LIBDIRS} ${LIBS} ${NVECSERLIB}

pvfkxCS: pvfkxCS.o pvfkxCSf.o pvfkxCSf.c_Cplx_.o
	$(CPP) ${CFLAGS} -o pvfkxCS pvfkxCS.o pvfkxCSf.o pvfkxCSf.c_Cplx_.o ${LIBDIRS} ${LIBS} ${NVECPARLIB}

#-------

clean:
	rm -f *.o
	rm -f cvfkxCS cvfdxCS pvfkxCS

.c.o:
	${CC} -c ${CFLAGS} ${INCDIRS} $*.c
.cpp.o:
	$(CPP) -c ${CFLAGS} ${INCDIRS} $*.cpp

#---End of Makefile---