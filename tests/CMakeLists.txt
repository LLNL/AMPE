find_package(PythonInterp 3.2 REQUIRED)

include_directories( ${CMAKE_SOURCE_DIR}/source )

set( fortran_dir ${CMAKE_SOURCE_DIR}/source/fortran )
set( fortran_dir_dim ${CMAKE_SOURCE_DIR}/source/fortran/${NDIM}d )
set( TEST_DIR ${CMAKE_SOURCE_DIR}/tests )

#Within the scope of this file, cmake doesn't know *.f are generated files
set_source_files_properties( ${fortran_dir}/functions.f 
                             PROPERTIES GENERATED TRUE)
set_source_files_properties( ${fortran_dir_dim}/hyprepoisson.f
                             PROPERTIES GENERATED TRUE)
set_source_files_properties( ${fortran_dir_dim}/phase.f 
                             PROPERTIES GENERATED TRUE)

add_executable(testNetCDFwrite ${CMAKE_SOURCE_DIR}/tests/simple_xy_wr.C)
add_executable(testNetCDFread ${CMAKE_SOURCE_DIR}/tests/simple_xy_rd.C)
add_executable(testMPI ${CMAKE_SOURCE_DIR}/tests/testMPI.C)
add_executable(testMath ${CMAKE_SOURCE_DIR}/tests/testMath.C )
add_executable(testNoise ${CMAKE_SOURCE_DIR}/tests/testNoise.C
                         ${CMAKE_SOURCE_DIR}/source/UniformNoise.C 
                         ${CMAKE_SOURCE_DIR}/source/NormalNoise.C )
add_executable(testInterpolation
               ${CMAKE_SOURCE_DIR}/tests/testInterpolationFunctions.C
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f )
add_executable(testHyprePoisson
               ${CMAKE_SOURCE_DIR}/tests/testHyprePoisson.C
               ${CMAKE_SOURCE_DIR}/source/HyprePoisson.C
               ${fortran_dir_dim}/hyprepoisson.f )
add_executable(testPhaseHypre
               ${CMAKE_SOURCE_DIR}/tests/testPhaseHypre.C
               ${CMAKE_SOURCE_DIR}/source/QuatModelParameters.C
               ${CMAKE_SOURCE_DIR}/source/CellPoissonHypreSolver.C
               ${CMAKE_SOURCE_DIR}/source/tools.C
               ${CMAKE_SOURCE_DIR}/source/PoissonSpecifications.C
               ${CMAKE_SOURCE_DIR}/source/PhaseHypre.C
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               ${fortran_dir_dim}/phase.f
               ${fortran_dir_dim}/hyprepoisson.f )
add_executable(testFACPoisson
               ${CMAKE_SOURCE_DIR}/tests/testFACPoisson.C
               ${CMAKE_SOURCE_DIR}/source/FACPoisson.C
               ${fortran_dir_dim}/hyprepoisson.f )
add_executable(testCALPHADFunctions
               ${CMAKE_SOURCE_DIR}/tests/testCALPHADFunctions.C )
add_executable(testCALPHADJacobianTernary
               ${CMAKE_SOURCE_DIR}/tests/testCALPHADJacobianTernary.C )
add_executable(testKKSdiluteBinaryEquilibrium
               ${CMAKE_SOURCE_DIR}/tests/testKKSdiluteBinaryEquilibrium.C
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               )
add_executable(testKKSdiluteBinaryKKS
               ${CMAKE_SOURCE_DIR}/tests/testKKSdiluteBinaryKKS.C
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               )
add_executable(testCALPHADbinaryKKS
               ${CMAKE_SOURCE_DIR}/tests/testCALPHADbinaryKKS.C
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               )
add_executable(testCALPHADternaryKKS
               ${CMAKE_SOURCE_DIR}/tests/testCALPHADternaryKKS.C
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               )
add_executable(testCALPHADbinaryEquilibrium
               ${CMAKE_SOURCE_DIR}/tests/testCALPHADbinaryEquilibrium.C
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               )
add_executable(testCALPHADternaryEquilibrium
               ${CMAKE_SOURCE_DIR}/tests/testCALPHADternaryEquilibrium.C
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               )

add_test(NAME testNetCDFwrite
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testNetCDFwrite)
add_test(NAME testNetCDFread
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testNetCDFread)
add_test(NAME testMPI
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
                 ${CMAKE_CURRENT_BINARY_DIR}/testMPI)
add_test(NAME testMath
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testMath)
add_test(NAME testNoise
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testNoise)
add_test(NAME testInterpolation
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testInterpolation)
add_test(NAME testHyprePoisson
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
                 ${CMAKE_CURRENT_BINARY_DIR}/testHyprePoisson ${NDIM}d.input
         WORKING_DIRECTORY ${TEST_DIR}/HyprePoisson)
add_test(NAME testPhaseHypre
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
                 ${CMAKE_CURRENT_BINARY_DIR}/testPhaseHypre ${NDIM}d.input
         WORKING_DIRECTORY ${TEST_DIR}/PhaseHypre)
add_test(NAME testFACPoisson
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
                 ${CMAKE_CURRENT_BINARY_DIR}/testFACPoisson ${NDIM}d.input
         WORKING_DIRECTORY ${TEST_DIR}/FACPoisson)
add_test(NAME testCALPHADFunctions
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADFunctions)
add_test(NAME testCALPHADJacobianTernary
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADJacobianTernary)
add_test(NAME testKKSdiluteBinaryEquilibrium
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testKKSdiluteBinaryEquilibrium test.input
         WORKING_DIRECTORY ${TEST_DIR}/KKSdiluteBinaryEquilibrium)
add_test(NAME testKKSdiluteBinaryKKS
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testKKSdiluteBinaryKKS test.input
         WORKING_DIRECTORY ${TEST_DIR}/KKSdiluteBinaryKKS)
add_test(NAME testCALPHADbinaryKKS
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADbinaryKKS test.input
         WORKING_DIRECTORY ${TEST_DIR}/CALPHADbinaryKKS)
add_test(NAME testCALPHADternaryKKS
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADternaryKKS test.input
         WORKING_DIRECTORY ${TEST_DIR}/CALPHADternaryKKS)
add_test(NAME testCALPHADbinaryEquilibrium
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
         ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADbinaryEquilibrium test.input
         WORKING_DIRECTORY ${TEST_DIR}/CALPHADbinaryEquilibrium)
add_test(NAME testCALPHADternaryEquilibrium1
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
         ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADternaryEquilibrium test1.input
         WORKING_DIRECTORY ${TEST_DIR}/CALPHADternaryEquilibrium)
add_test(NAME testCALPHADternaryEquilibrium2
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
         ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADternaryEquilibrium test2.input
         WORKING_DIRECTORY ${TEST_DIR}/CALPHADternaryEquilibrium)

#regression tests
add_test(NAME PlanarFront
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/PlanarFront/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}//PlanarFront/${NDIM}d.input
         )

add_test(NAME PlanarFrontMovingFrame
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/PlanarFrontMovingFrame/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/PlanarFrontMovingFrame/${NDIM}d.input
         )

add_test(NAME AlCuMovingFrame
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/AlCuMovingFrame/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/AlCuMovingFrame/${NDIM}d.input
         )

add_test(NAME ChenShen98
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ChenShen98/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/ChenShen98/${NDIM}d.input
         )
add_test(NAME HeatFlow
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/HeatFlow/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/HeatFlow/${NDIM}d.input
         )
add_test(NAME TimeDependentHeatFlow
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentHeatFlow/test.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentHeatFlow/${NDIM}d.input
         ${CMAKE_CURRENT_SOURCE_DIR}/TimeDependentHeatFlow
         )
add_test(NAME SingleGrainGrowthAuNi
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SingleGrainGrowthAuNi/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/SingleGrainGrowthAuNi/${NDIM}d.input
         ${CMAKE_CURRENT_SOURCE_DIR}/../thermodynamic_data
         )
add_test(NAME AlCu
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/AlCu/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/AlCu/${NDIM}d.input
         )
add_test(NAME DiluteCuNi
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/DiluteCuNi/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/DiluteCuNi/${NDIM}d.input
         )
add_test(NAME GaussianT
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/GaussianT/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/GaussianT/${NDIM}d.input
         )
add_test(NAME TwoBilayers
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/TwoBilayers/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/TwoBilayers/${NDIM}d.input
         )

target_link_libraries(testNetCDFwrite ${PROJECT_LINK_LIBS_NETCDF}
                                 ${PROJECT_LINK_LIBS_HDF5})
target_link_libraries(testNetCDFread ${PROJECT_LINK_LIBS_NETCDF}
                                 ${PROJECT_LINK_LIBS_HDF5})
target_link_libraries(testMPI ${PROJECT_LINK_LIBS1})
target_link_libraries(testHyprePoisson ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HYPRE}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     )
target_link_libraries(testPhaseHypre ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HYPRE}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     ${internal_libs}
                                     )
target_link_libraries(testFACPoisson ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HYPRE}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     )
target_link_libraries(testKKSdiluteBinaryEquilibrium ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     ${internal_libs}
                                     )
if(${USE_THERMO4PFM})
  target_link_libraries(testKKSdiluteBinaryEquilibrium ${PROJECT_LINK_LIBS_THERMO4PFM})
endif()
target_link_libraries(testKKSdiluteBinaryKKS ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     ${internal_libs}
                                     )
target_link_libraries(testCALPHADbinaryKKS ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     ${internal_libs}
                                     )
target_link_libraries(testCALPHADternaryKKS ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     ${internal_libs}
                                     )
target_link_libraries(testCALPHADbinaryEquilibrium ${PROJECT_LINK_LIBS_SAMRAI}
                                                   ${internal_libs}
                     )
target_link_libraries(testCALPHADternaryEquilibrium
                      ${PROJECT_LINK_LIBS_SAMRAI}
                      ${internal_libs}
                     )
target_link_libraries(testCALPHADFunctions
                      ${internal_libs})
target_link_libraries(testCALPHADJacobianTernary
                      ${internal_libs})
target_link_libraries(testMath
                      ${internal_libs})
