include_directories( ${CMAKE_SOURCE_DIR}/source )

set( fortran_dir ${CMAKE_SOURCE_DIR}/source/fortran )
set( fortran_dir_dim ${CMAKE_SOURCE_DIR}/source/fortran/${NDIM}d )

#Within the scope of this file, cmake doesn't know *.f are generated files
set_source_files_properties( ${fortran_dir}/functions.f 
                             PROPERTIES GENERATED TRUE)
set_source_files_properties( ${fortran_dir_dim}/hyprepoisson.f
                             PROPERTIES GENERATED TRUE)
set_source_files_properties( ${fortran_dir_dim}/phase.f 
                             PROPERTIES GENERATED TRUE)

add_executable(testNetCDF ${CMAKE_SOURCE_DIR}/tests/simple_xy_wr.C)
add_executable(testMPI ${CMAKE_SOURCE_DIR}/tests/testMPI.C)
add_executable(testMath ${CMAKE_SOURCE_DIR}/tests/testMath.C
                        ${CMAKE_SOURCE_DIR}/source/math_utilities.C )
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
               ${CMAKE_SOURCE_DIR}/tests/testCALPHADFunctions.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADFunctions.C )
add_executable(testCALPHADJacobianTernary
               ${CMAKE_SOURCE_DIR}/tests/testCALPHADJacobianTernary.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADEqPhaseConcSolverTernary.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADConcSolverTernary.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADFunctions.C
               ${CMAKE_SOURCE_DIR}/source/NewtonSolver.C
               ${CMAKE_SOURCE_DIR}/source/DampedNewtonSolver.C
               ${CMAKE_SOURCE_DIR}/source/math_utilities.C )
add_executable(testCALPHADbinaryKKS
               ${CMAKE_SOURCE_DIR}/tests/testCALPHADbinaryKKS.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyFunctionsBinary.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADSpeciesPhaseGibbsEnergy.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADSpeciesPhaseGibbsEnergyExpansion.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADFunctions.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADConcSolverBinary.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADEqConcSolverBinary.C
               ${CMAKE_SOURCE_DIR}/source/NewtonSolver.C
               ${CMAKE_SOURCE_DIR}/source/DampedNewtonSolver.C
               ${CMAKE_SOURCE_DIR}/source/math_utilities.C
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               )
add_executable(testCALPHADternaryKKS
               ${CMAKE_SOURCE_DIR}/tests/testCALPHADternaryKKS.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyFunctionsTernary.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADSpeciesPhaseGibbsEnergy.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADSpeciesPhaseGibbsEnergyExpansion.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADFunctions.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADConcSolverTernary.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADEqConcSolverTernary.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADEqPhaseConcSolverTernary.C
               ${CMAKE_SOURCE_DIR}/source/NewtonSolver.C
               ${CMAKE_SOURCE_DIR}/source/DampedNewtonSolver.C
               ${CMAKE_SOURCE_DIR}/source/math_utilities.C
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               )

add_test(NAME testNetCDF
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testNetCDF)
add_test(NAME testMPI
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 
                 ${MPIEXEC_MAX_NUMPROCS} 
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
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG}
                 ${MPIEXEC_MAX_NUMPROCS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testHyprePoisson ${NDIM}d.input
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/HyprePoisson)
add_test(NAME testPhaseHypre
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG}
                 ${MPIEXEC_MAX_NUMPROCS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testPhaseHypre ${NDIM}d.input
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/PhaseHypre)
add_test(NAME testFACPoisson
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 
                 ${MPIEXEC_MAX_NUMPROCS}
                 ${CMAKE_CURRENT_BINARY_DIR}/testFACPoisson ${NDIM}d.input
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/FACPoisson)
add_test(NAME testCALPHADFunctions
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADFunctions)
add_test(NAME testCALPHADJacobianTernary
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADJacobianTernary)
add_test(NAME testCALPHADbinaryKKS
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADbinaryKKS test.input
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/CALPHADbinaryKKS)
add_test(NAME testCALPHADternaryKKS
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testCALPHADternaryKKS test.input
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/CALPHADternaryKKS)

target_link_libraries(testNetCDF ${PROJECT_LINK_LIBS_NETCDFCXX}
                                 ${PROJECT_LINK_LIBS_NETCDF}
                                 ${PROJECT_LINK_LIBS_HDF5})
target_link_libraries(testMPI ${PROJECT_LINK_LIBS1})
target_link_libraries(testHyprePoisson ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HYPRE}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     )
target_link_libraries(testPhaseHypre ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HYPRE}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     )
target_link_libraries(testFACPoisson ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HYPRE}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     )
target_link_libraries(testCALPHADbinaryKKS ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     )
target_link_libraries(testCALPHADternaryKKS ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     )

