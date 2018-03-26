#--------------------------------------------------------------------------------
# File        : Makefile
# Programmers : Radu Serban @ LLNL
#--------------------------------------------------------------------------------

# Shell
SHELL = /bin/sh

# SUNDIALS instalation directory
prefix      = /home/radu/CODES/sundials/LINUX
exec_prefix = ${prefix}
includedir  = ${prefix}/include
libdir      = ${exec_prefix}/lib

# Path to sources
srcdir = ../../../levmar/source
VPATH  = ../../../levmar/source

# Archiver program
AR = ar rc

# Compiler and compiler options
CC       = gcc
CPPFLAGS = 
CFLAGS   = -O -ffloat-store -g

# LEVMAR and SHARED include directories
SHARED_INC_DIR = ${srcdir}/../../shared/include
LEVMAR_INC_DIR  = ${srcdir}/../include

INCDIRS = -I${LEVMAR_INC_DIR} -I${SHARED_INC_DIR}

# LEVMAR library name
LEVMAR_LIB = liblevmar.a

# Source and object files
SRC =	${srcdir}/levmar.c \
	${srcdir}/lmdense.c \
	${srcdir}/lmspgmr.c

OBJ = levmar.o lmdense.o lmspgmr.o

#--------------------------------------------------------------------------------
# Make rules
#--------------------------------------------------------------------------------

all: ${OBJ}

install: ${OBJ}
	${AR} ${LEVMAR_LIB} ${OBJ}
	cp ${LEVMAR_LIB} ${libdir}
	cp ${LEVMAR_INC_DIR}/levmar.h    ${includedir}
	cp ${LEVMAR_INC_DIR}/lmdense.h  ${includedir}
	cp ${LEVMAR_INC_DIR}/lmspgmr.h  ${includedir}

clean:
	rm -f ${OBJ}
	rm -f ${LEVMAR_LIB}

uninstall:
	rm -f ${libdir}/${LEVMAR_LIB}
	rm -f ${includedir}/levmar.h
	rm -f ${includedir}/lmdense.h
	rm -f ${includedir}/lmspgmr.h

.c.o:
	${CC} -c ${CFLAGS} ${CPPFLAGS} ${INCDIRS} -o $@ $<

#---End of Makefile---
