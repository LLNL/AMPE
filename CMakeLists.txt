cmake_minimum_required(VERSION 3.1)
project (PFiSM)

set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#set(MPI_CXX_COMPILER CMAKE_CXX_COMPILER)
find_package(MPI)
if(MPI_FOUND)
  message(STATUS
    "MPI was found:\n"
    "  MPI_CXX_COMPILER: ${MPI_CXX_COMPILER}\n"
    "  MPIEXEC_NUMPROC_FLAG: ${MPIEXEC_NUMPROC_FLAG}\n"
    "  MPIEXEC_MAX_NUMPROCS: ${MPIEXEC_MAX_NUMPROCS}\n"
    "  MPIEXEC: ${MPIEXEC}\n"
 )
endif()
if(DEFINED MPI_EXEC)
  set(MPIEXEC ${MPI_EXEC})
endif()

message(STATUS "MPIEXEC: ${MPIEXEC}\n")

set(HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5 COMPONENTS CXX HL REQUIRED)
message(STATUS ${HDF5_VERSION} )
include_directories(${HDF5_INCLUDE_DIR})
#set(_hdf5_libs hdf5 hdf5_cpp)
message(STATUS "hdf5 libraries: ${HDF5_LIBRARIES}")
set(_hdf5_libs ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES})
message(STATUS "hdf5 libraries: ${_hdf5_libs}")

find_package(Boost)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

link_directories( ${SAMRAI_DIR}/lib )
link_directories( ${HYPRE_DIR}/lib )
link_directories( ${SUNDIALS_DIR}/lib )
link_directories( ${NETCDF_DIR}/lib )

set( PROJECT_LINK_LIBS1 libSAMRAI_tbox.a )
set( PROJECT_LINK_LIBS2 libSAMRAI_math.a )
set( PROJECT_LINK_LIBS3 libSAMRAI_hier.a )
set( PROJECT_LINK_LIBS4 libSAMRAI_pdat.a )
set( PROJECT_LINK_LIBS5 libSAMRAI_xfer.a )
set( PROJECT_LINK_LIBS6 libSAMRAI_geom.a )
set( PROJECT_LINK_LIBS7 libSAMRAI_mesh.a )
set( PROJECT_LINK_LIBS8 libSAMRAI_solv.a )
set( PROJECT_LINK_LIBS9 libSAMRAI_appu.a )
set( PROJECT_LINK_LIBS_SAMRAI ${PROJECT_LINK_LIBS9}
                              ${PROJECT_LINK_LIBS8}
                              ${PROJECT_LINK_LIBS7}
                              ${PROJECT_LINK_LIBS6}
                              ${PROJECT_LINK_LIBS5}
                              ${PROJECT_LINK_LIBS4}
                              ${PROJECT_LINK_LIBS3}
                              ${PROJECT_LINK_LIBS2}
                              ${PROJECT_LINK_LIBS1} )

set( PROJECT_LINK_LIBS_HYPRE libHYPRE.so )

set( PROJECT_LINK_LIBS_CVODE libsundials_cvode.a )
set( PROJECT_LINK_LIBS_CPODES libsundials_cpodes.a )

set( PROJECT_LINK_LIBS_NETCDFCXX libnetcdf_c++4.a )
set( PROJECT_LINK_LIBS_NETCDF libnetcdf.a )

include_directories( ${SAMRAI_DIR}/include )
include_directories( ${HYPRE_DIR}/include )
include_directories( ${SUNDIALS_DIR}/include )
include_directories( ${NETCDF_DIR}/include )

enable_testing()

add_definitions(-DNDIM=${NDIM})

enable_language(Fortran C CXX)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -std=c++11")

add_subdirectory(source)

add_subdirectory(tests)



