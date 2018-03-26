# Shell
SHELL = /bin/sh

# SUNDIALS instalation directory
prefix      = /home/radu/CODES/sundials/LINUX
exec_prefix = ${prefix}
includedir  = ${prefix}/include
libdir      = ${exec_prefix}/lib

# Path to sources
srcdir = ../../../levmar/examples_ser

# Compiler and compiler options
CC       = gcc
CPPFLAGS = 
CFLAGS   = -O -ffloat-store -g
LDFLAGS  = 
TEMP_LIBS     =  

# SUNDIALS library and include directories
SUNDIALS_INC_DIR = ${includedir}
SUNDIALS_LIB_DIR = ${libdir}

INCDIRS = -I${SUNDIALS_INC_DIR}
LIBDIRS = -L${SUNDIALS_LIB_DIR}
LIBS    = ${TEMP_LIBS} -llevmar -lshared -lnvecserial -lm

#--------------------------------------------------------------------------------
# Make rules
#--------------------------------------------------------------------------------

all: ros_d ros_k mad_d mad_k 

ros_k : ${srcdir}/ros_k.c
	${CC} -c ${CFLAGS} ${CPPFLAGS} ${INCDIRS} ${srcdir}/$@.c
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $@.o ${LIBDIRS} ${LIBS}

ros_d : ${srcdir}/ros_d.c
	${CC} -c ${CFLAGS} ${CPPFLAGS} ${INCDIRS} ${srcdir}/$@.c
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $@.o ${LIBDIRS} ${LIBS}

mad_k : ${srcdir}/mad_k.c
	${CC} -c ${CFLAGS} ${CPPFLAGS} ${INCDIRS} ${srcdir}/$@.c
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $@.o ${LIBDIRS} ${LIBS}

mad_d : ${srcdir}/mad_d.c
	${CC} -c ${CFLAGS} ${CPPFLAGS} ${INCDIRS} ${srcdir}/$@.c
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ $@.o ${LIBDIRS} ${LIBS}

clean:
	rm -f *.o
	rm -f ros_d ros_k mad_d mad_k

#---End of Makefile---
