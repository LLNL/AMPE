find_package(PythonInterp 3.2 REQUIRED)

include_directories( ${CMAKE_SOURCE_DIR}/source )

set( fortran_dir ${CMAKE_SOURCE_DIR}/source/fortran )
set( SOURCE_FORTRAN_DIM_DIR ${CMAKE_SOURCE_DIR}/source/fortran/${NDIM}d )
set( TEST_DIR ${CMAKE_SOURCE_DIR}/tests )

if(${USE_THERMO4PFM})
  set(THERMOLIB ${PROJECT_LINK_LIBS_THERMO4PFM})
else()
  set(THERMOLIB ${internal_libs})
endif(${USE_THERMO4PFM})


set(M4FILESDIM ellipticfacops
               hyprepoisson
               phase
               quatdiffs
               quatgrad
               concentrationrhs)

foreach(file ${M4FILESDIM})
  message(STATUS ".m4 -> .f: ${file}...")
  m4tof(${file} ${SOURCE_FORTRAN_DIM_DIR})
endforeach()

set(SOURCES_FORTRAN_DIM "")
foreach(file ${M4FILESDIM})
  list(APPEND SOURCES_FORTRAN_DIM ${SOURCE_FORTRAN_DIM_DIR}/${file}.f)
endforeach()



#Within the scope of this file, cmake doesn't know *.f are generated files
set_source_files_properties( ${fortran_dir}/functions.f 
                             PROPERTIES GENERATED TRUE)
set_source_files_properties( ${SOURCE_FORTRAN_DIM_DIR}/hyprepoisson.f
                             PROPERTIES GENERATED TRUE)
set_source_files_properties( ${SOURCE_FORTRAN_DIM_DIR}/phase.f
                             PROPERTIES GENERATED TRUE)
set_source_files_properties( ${SOURCE_FORTRAN_DIM_DIR}/quatdiffs.f
                             PROPERTIES GENERATED TRUE)
set_source_files_properties( ${SOURCE_FORTRAN_DIM_DIR}/quatgrad.f
                             PROPERTIES GENERATED TRUE)
set_source_files_properties( ${SOURCE_FORTRAN_DIM_DIR}/concentrationrhs.f
                             PROPERTIES GENERATED TRUE)

add_executable(testNetCDFwrite ${CMAKE_SOURCE_DIR}/tests/simple_xy_wr.cc)
add_executable(testNetCDFread ${CMAKE_SOURCE_DIR}/tests/simple_xy_rd.cc)
add_executable(testHDF5 ${CMAKE_SOURCE_DIR}/tests/testHDF5.cc)
add_executable(testMPI ${CMAKE_SOURCE_DIR}/tests/testMPI.cc)
add_executable(testSAMRAIcommunication
               ${CMAKE_SOURCE_DIR}/tests/SAMRAIcommunication/main.cpp
               ${CMAKE_SOURCE_DIR}/tests/SAMRAIcommunication/PatchDataTestStrategy.cpp
               ${CMAKE_SOURCE_DIR}/tests/SAMRAIcommunication/CommTester.cpp
               ${CMAKE_SOURCE_DIR}/tests/SAMRAIcommunication/CellDataTest.cpp)
add_executable(testMath ${CMAKE_SOURCE_DIR}/tests/testMath.cc )
add_executable(testNoise ${CMAKE_SOURCE_DIR}/tests/testNoise.cc
                         ${CMAKE_SOURCE_DIR}/source/UniformNoise.cc 
                         ${CMAKE_SOURCE_DIR}/source/NormalNoise.cc )
add_executable(testInterpolation
               ${CMAKE_SOURCE_DIR}/tests/testInterpolationFunctions.cc
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f )
add_executable(test3PhasesInterpolation
               ${CMAKE_SOURCE_DIR}/tests/test3PhasesInterpolationFunctions.cc
               ${CMAKE_SOURCE_DIR}/source/TiltingFolchPlapp2005.cc
               ${CMAKE_SOURCE_DIR}/source/TiltingMoelans2011.cc )
add_executable(testWeights3Pairs
               ${CMAKE_SOURCE_DIR}/tests/testWeights3Pairs.cc
               ${CMAKE_SOURCE_DIR}/source/computeWeights3Pairs.cc )
add_executable(testHyprePoisson
               ${CMAKE_SOURCE_DIR}/tests/testHyprePoisson.cc
               ${CMAKE_SOURCE_DIR}/tests/HyprePoisson.cc
               ${SOURCE_FORTRAN_DIM_DIR}/hyprepoisson.f )
add_executable(testPhaseHypre
               ${CMAKE_SOURCE_DIR}/tests/testPhaseHypre.cc
               ${CMAKE_SOURCE_DIR}/source/QuatModelParameters.cc
               ${CMAKE_SOURCE_DIR}/source/samrai/CellPoissonHypreSolver.cc
               ${CMAKE_SOURCE_DIR}/source/tools.cc
               ${CMAKE_SOURCE_DIR}/source/samrai/PoissonSpecifications.cc
               PhaseHypre.cc
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               ${SOURCE_FORTRAN_DIM_DIR}/phase.f
               ${SOURCE_FORTRAN_DIM_DIR}/hyprepoisson.f )
add_executable(testFACPoisson
               ${CMAKE_SOURCE_DIR}/tests/testFACPoisson.cc
               ${CMAKE_SOURCE_DIR}/tests/FACPoisson.cc
               ${SOURCE_FORTRAN_DIM_DIR}/hyprepoisson.f )
add_executable(testPhaseFAC
               ${CMAKE_SOURCE_DIR}/tests/testPhaseFAC.cc
               ${CMAKE_SOURCE_DIR}/source/PhaseFACOps.cc
               ${CMAKE_SOURCE_DIR}/source/PhaseFACSolver.cc
               ${CMAKE_SOURCE_DIR}/source/QuatModelParameters.cc
               ${CMAKE_SOURCE_DIR}/source/EllipticFACOps.cc
               ${CMAKE_SOURCE_DIR}/source/EllipticFACSolver.cc
               ${CMAKE_SOURCE_DIR}/source/CartesianRobinBcHelperWithDepth.cc
               ${CMAKE_SOURCE_DIR}/source/tools.cc
               ${CMAKE_SOURCE_DIR}/source/samrai/FACPreconditioner.cc
               ${CMAKE_SOURCE_DIR}/source/samrai/PoissonSpecifications.cc
               ${CMAKE_SOURCE_DIR}/source/samrai/CellPoissonHypreSolver.cc
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               ${SOURCE_FORTRAN_DIM_DIR}/phase.f
               ${SOURCE_FORTRAN_DIM_DIR}/hyprepoisson.f
               ${SOURCE_FORTRAN_DIM_DIR}/ellipticfacops.f
               PhaseFAC.cc)
add_executable(testGradQ
               ${CMAKE_SOURCE_DIR}/tests/testGradQ.cc
               ${CMAKE_SOURCE_DIR}/source/computeQDiffs.cc
               ${SOURCE_FORTRAN_DIM_DIR}/quatdiffs.f
               ${SOURCE_FORTRAN_DIM_DIR}/quatgrad.f
               ${CMAKE_SOURCE_DIR}/source/fortran/quat.f )
add_executable(testAntitrapping3Phases
               ${CMAKE_SOURCE_DIR}/tests/testAntitrapping3Phases.cc
               ${SOURCE_FORTRAN_DIM_DIR}/concentrationrhs.f
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f )
add_executable(testCALPHADFunctions
               ${CMAKE_SOURCE_DIR}/tests/testCALPHADFunctions.cc )
add_executable(testCALPHADJacobianTernary
               ${CMAKE_SOURCE_DIR}/tests/testCALPHADJacobianTernary.cc )
add_executable(testKKSdiluteBinaryEquilibrium
               ${CMAKE_SOURCE_DIR}/tests/testKKSdiluteBinaryEquilibrium.cc
               ${CMAKE_SOURCE_DIR}/source/Database2JSON.cc
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               )
add_executable(testKKSdiluteBinaryKKS
               ${CMAKE_SOURCE_DIR}/tests/testKKSdiluteBinaryKKS.cc
               ${CMAKE_SOURCE_DIR}/source/Database2JSON.cc
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               )
add_executable(testCALPHADbinaryKKS
               ${CMAKE_SOURCE_DIR}/tests/testCALPHADbinaryKKS.cc
               ${CMAKE_SOURCE_DIR}/source/Database2JSON.cc
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               )
add_executable(testCALPHADternaryKKS
               ${CMAKE_SOURCE_DIR}/tests/testCALPHADternaryKKS.cc
               ${CMAKE_SOURCE_DIR}/source/Database2JSON.cc
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               )
add_executable(testCALPHADbinaryEquilibrium
               ${CMAKE_SOURCE_DIR}/tests/testCALPHADbinaryEquilibrium.cc
               ${CMAKE_SOURCE_DIR}/source/Database2JSON.cc
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               )
add_executable(testCALPHADternaryEquilibrium
               ${CMAKE_SOURCE_DIR}/tests/testCALPHADternaryEquilibrium.cc
               ${CMAKE_SOURCE_DIR}/source/Database2JSON.cc
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               )
add_executable(testTemperatureHistory
               ${CMAKE_SOURCE_DIR}/tests/testTemperatureHistory.cc
               ${CMAKE_SOURCE_DIR}/source/TemperatureHistory.cc
               )

add_test(NAME testNetCDFwrite
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testNetCDFwrite)
add_test(NAME testNetCDFread
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testNetCDFread)
add_test(NAME testHDF5
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testHDF5)
add_test(NAME testMPI
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
                 ${CMAKE_CURRENT_BINARY_DIR}/testMPI)
add_test(NAME testSAMRAIcommunication
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
                 ${CMAKE_CURRENT_BINARY_DIR}/testSAMRAIcommunication
                 cell_refine_a.${NDIM}d.input
         WORKING_DIRECTORY ${TEST_DIR}/SAMRAIcommunication)
add_test(NAME testMath
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testMath)
add_test(NAME testNoise
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testNoise)
add_test(NAME testInterpolation
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testInterpolation)
add_test(NAME test3PhasesInterpolation
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/test3PhasesInterpolation)
add_test(NAME testWeights3Pairs
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testWeights3Pairs)
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
add_test(NAME testPhaseFAC
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
                 ${CMAKE_CURRENT_BINARY_DIR}/testPhaseFAC ${NDIM}d.input
         WORKING_DIRECTORY ${TEST_DIR}/PhaseFAC)
add_test(NAME testGradQ
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
                 ${CMAKE_CURRENT_BINARY_DIR}/testGradQ)
add_test(NAME testAntitrapping3Phases
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
                 ${CMAKE_CURRENT_BINARY_DIR}/testAntitrapping3Phases)
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
add_test(NAME testTemperatureHistory
         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/testTemperatureHistory
         WORKING_DIRECTORY ${TEST_DIR}/TemperatureHistory)

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

add_test(NAME GradientMovingFrame
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/GradientMovingFrame/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/GradientMovingFrame/${NDIM}d.input
         )
add_test(NAME TemperatureHistoryMovingFrame
         COMMAND ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/TemperatureHistoryMovingFrame/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/TemperatureHistoryMovingFrame/${NDIM}d.input
         ${CMAKE_CURRENT_SOURCE_DIR}/TemperatureHistoryMovingFrame
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
         ${CMAKE_CURRENT_SOURCE_DIR}/../thermodynamic_data
         )
add_test(NAME TwoGrainsQuadratic
         COMMAND ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/TwoGrainsQuadratic/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/TwoGrainsQuadratic/${NDIM}d.input
         )
add_test(NAME ConservedVolume
         COMMAND ${PYTHON_EXECUTABLE}
         ${CMAKE_CURRENT_SOURCE_DIR}/ConservedVolume/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/ConservedVolume/${NDIM}d.input
         )
add_test(NAME FourCorners
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FourCorners/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/FourCorners/${NDIM}d.input
         )
add_test(NAME SolidifyQuaternions
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SolidifyQuaternions/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/SolidifyQuaternions/${NDIM}d.input
         )
add_test(NAME Dendrite
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Dendrite/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/Dendrite/${NDIM}d.input
         )
add_test(NAME FeNbNiWithNoise
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FeNbNiWithNoise/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/FeNbNiWithNoise/${NDIM}d.input
         ${CMAKE_CURRENT_SOURCE_DIR}/../thermodynamic_data
         )
add_test(NAME FeNbNiWithFluctuation
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/FeNbNiWithFluctuation/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/FeNbNiWithFluctuation/${NDIM}d.input
         ${CMAKE_CURRENT_SOURCE_DIR}/../thermodynamic_data
         )
if(${USE_THERMO4PFM})
  add_test(NAME ThreePhases
           COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ThreePhases/test${NDIM}d.py
           ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
           ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
           ${CMAKE_CURRENT_SOURCE_DIR}/ThreePhases/${NDIM}d.input
           )
  add_test(NAME ThreePhasesCALPHAD
           COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ThreePhasesCALPHAD/test${NDIM}d.py
           ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
           ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
           ${CMAKE_CURRENT_SOURCE_DIR}/ThreePhasesCALPHAD/${NDIM}d.input
           ${CMAKE_CURRENT_SOURCE_DIR}/ThreePhasesCALPHAD
           )
  add_test(NAME ThreePhasesDiffusion
           COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ThreePhasesDiffusion/test${NDIM}d.py
           ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
           ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
           ${CMAKE_CURRENT_SOURCE_DIR}/ThreePhasesDiffusion/${NDIM}d.input
           ${CMAKE_CURRENT_SOURCE_DIR}/ThreePhasesDiffusion
           )
  add_test(NAME 3Ph2Sl
           COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/3Ph2Sl/test${NDIM}d.py
           ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
           ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
           ${CMAKE_CURRENT_SOURCE_DIR}/3Ph2Sl/${NDIM}d.input
           ${CMAKE_CURRENT_SOURCE_DIR}/3Ph2Sl
           )
  add_test(NAME 2Ph1Sl
           COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/2Ph1Sl/test${NDIM}d.py
           ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
           ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
           ${CMAKE_CURRENT_SOURCE_DIR}/2Ph1Sl/${NDIM}d.input
           ${CMAKE_CURRENT_SOURCE_DIR}/2Ph1Sl
           )
endif()
add_test(NAME DumpNetCDF
         COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/DumpNetCDF/test${NDIM}d.py
         ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4
         ${CMAKE_CURRENT_BINARY_DIR}/../source/ampe${NDIM}d
         ${CMAKE_CURRENT_SOURCE_DIR}/DumpNetCDF/${NDIM}d.input
         )

target_link_libraries(testNetCDFwrite ${PROJECT_LINK_LIBS_NETCDF}
                                 ${PROJECT_LINK_LIBS_HDF5})
target_link_libraries(testNetCDFread ${PROJECT_LINK_LIBS_NETCDF}
                                 ${PROJECT_LINK_LIBS_HDF5})
target_link_libraries(testHDF5 ${PROJECT_LINK_LIBS_HDF5})
target_link_libraries(testMPI ${PROJECT_LINK_LIBS1})
target_link_libraries(testSAMRAIcommunication ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     )
target_link_libraries(testHyprePoisson ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HYPRE}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     )
target_link_libraries(testPhaseHypre ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HYPRE}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     ${THERMOLIB}
                                     )
target_link_libraries(testFACPoisson ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HYPRE}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     )
target_link_libraries(testPhaseFAC ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HYPRE}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     ${THERMOLIB}
                                     )
target_link_libraries(testGradQ ${PROJECT_LINK_LIBS_SAMRAI}
                                ${PROJECT_LINK_LIBS_HDF5}
                                )
target_link_libraries(testAntitrapping3Phases ${PROJECT_LINK_LIBS_SAMRAI}
                                              ${PROJECT_LINK_LIBS_HDF5}
                      )
target_link_libraries(testKKSdiluteBinaryEquilibrium ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     ${THERMOLIB}
                                     )
if(${USE_THERMO4PFM})
  target_link_libraries(testKKSdiluteBinaryEquilibrium ${PROJECT_LINK_LIBS_THERMO4PFM})
endif()
target_link_libraries(testKKSdiluteBinaryKKS ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     ${THERMOLIB}
                                     )
target_link_libraries(testCALPHADbinaryKKS ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     ${THERMOLIB}
                                     )
target_link_libraries(testCALPHADternaryKKS ${PROJECT_LINK_LIBS_SAMRAI}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     ${PROJECT_LINK_LIBS_HDF5}
                                     ${THERMOLIB}
                                     )
target_link_libraries(testCALPHADbinaryEquilibrium ${PROJECT_LINK_LIBS_SAMRAI}
                                                   ${THERMOLIB}
                     )
target_link_libraries(testCALPHADternaryEquilibrium
                      ${PROJECT_LINK_LIBS_SAMRAI}
                      ${THERMOLIB}
                     )
target_link_libraries(testCALPHADFunctions
                      ${THERMOLIB})
target_link_libraries(testCALPHADJacobianTernary
                      ${THERMOLIB})
target_link_libraries(testMath
                      ${THERMOLIB})
target_link_libraries(testTemperatureHistory
                      ${PROJECT_LINK_LIBS_SAMRAI})
