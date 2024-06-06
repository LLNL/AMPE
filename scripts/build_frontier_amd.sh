#!/bin/tcsh
module load PrgEnv-amd
module load craype-accel-amd-gfx90a
module load boost
module load cmake
module load cray-python
module load sundials/5.8.0-cpu
module load cray-hdf5/1.12.2.9
module load cray-netcdf/4.9.0.9


rm -rf CMakeCache.txt
rm -rf CMakeFiles/
rm cmake_install.cmake
rm Makefile
rm -f ../source/fortran/3d/*.f

#export HDF5_DIR=/opt/cray/pe/hdf5/1.12.2.9/amd/4.3


cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc \
        -DCMAKE_Fortran_COMPILER=ftn \
	-DCMAKE_CXX_FLAGS="-fopenmp -fopenmp-assume-no-thread-state -fopenmp-cuda-mode -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx90a" \
	-DHDF5_DIR=/opt/cray/pe/hdf5/1.12.2.9/amd/5.0 \
        -DMPIEXEC_EXECUTABLE="/usr/bin/srun" \
        -DMPIEXEC_NUMPROCS_FLAG="-n" \
        -DMPIEXEC_PREFLAGS="-c1;--gpus=1" \
        -DCMAKE_EXE_LINKER_FLAGS="-lm -fopenmp" \
        -DSAMRAI_DIR=$HOME/SAMRAI/amd/4.1.2 \
        -DHYPRE_DIR=$HOME/hypre/amd \
        -DSUNDIALS_DIR=$OLCF_SUNDIALS_ROOT \
        -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF \
        -DNDIM="3" \
        -DTHERMO4PFM_DIR=$HOME/Thermo4PFM/amd \
        ..

make -j6
