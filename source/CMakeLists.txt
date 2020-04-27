# we need the M4 macro processor
find_program( M4_EXECUTABLE m4 DOC "The M4 macro processor" )
if( NOT M4_EXECUTABLE )
   message( SEND_ERROR "Failed to find the M4 macro processor." )
endif( NOT M4_EXECUTABLE )

# script used to convert m4 files
MACRO(m4tof filename dir)
  IF(NOT EXISTS ${dir}/${filename}.f)
    ADD_CUSTOM_COMMAND(
      OUTPUT ${dir}/${filename}.f
      COMMAND m4 
      ARGS -DSAMRAI_FORTDIR=${CMAKE_SOURCE_DIR}/source/fortran 
          ${dir}/${filename}.m4 
             > ${dir}/${filename}.f
      VERBATIM
      )
  ENDIF()
ENDMACRO(m4tof)

set(SOURCES ${CMAKE_SOURCE_DIR}/source/main.C
            ${CMAKE_SOURCE_DIR}/source/InterpolationType.C
            ${CMAKE_SOURCE_DIR}/source/KKSdiluteEquilibriumPhaseConcentrationsStrategy.C
            ${CMAKE_SOURCE_DIR}/source/xlogx.C
            ${CMAKE_SOURCE_DIR}/source/KKSdiluteBinaryConcentrationSolver.C
            ${CMAKE_SOURCE_DIR}/source/KKSFreeEnergyFunctionDiluteBinary.C
            ${CMAKE_SOURCE_DIR}/source/KKSdiluteBinary.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADFunctions.C
            ${CMAKE_SOURCE_DIR}/source/NewtonSolver.C
            ${CMAKE_SOURCE_DIR}/source/DampedNewtonSolver.C
            ${CMAKE_SOURCE_DIR}/source/math_utilities.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyFunctionsBinary.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADConcSolverTernary.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADEqConcSolverTernary.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyFunctionsTernary.C
            ${CMAKE_SOURCE_DIR}/source/QuatModelParameters.C
            ${CMAKE_SOURCE_DIR}/source/CompositionStrategyMobilities.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADSpeciesPhaseGibbsEnergy.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADMobility.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyStrategyBinary.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyStrategyTernary.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyFunctionsWithPenaltyBinary.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADEqConcSolverBinaryWithPenalty.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADEqConcSolverBinary.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADEqPhaseConcSolverTernary.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADConcSolverBinaryWithPenalty.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADConcSolverBinary.C
            ${CMAKE_SOURCE_DIR}/source/tools.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADSpeciesPhaseGibbsEnergyExpansion.C
            ${CMAKE_SOURCE_DIR}/source/FreeEnergyStrategy.C
            ${CMAKE_SOURCE_DIR}/source/AMPE.C
            ${CMAKE_SOURCE_DIR}/source/Noise.C
            ${CMAKE_SOURCE_DIR}/source/NormalNoise.C
            ${CMAKE_SOURCE_DIR}/source/UniformNoise.C
            ${CMAKE_SOURCE_DIR}/source/KimMobilityStrategyInfMob.C
            ${CMAKE_SOURCE_DIR}/source/KimMobilityStrategyFiniteMob.C
            ${CMAKE_SOURCE_DIR}/source/KimMobilityStrategyFiniteMobAntiTrap.C
            ${CMAKE_SOURCE_DIR}/source/KimMobilityStrategy.C
            ${CMAKE_SOURCE_DIR}/source/FieldsInitializer.C
            ${CMAKE_SOURCE_DIR}/source/GradientTemperatureStrategy.C
            ${CMAKE_SOURCE_DIR}/source/toolsSAMRAI.C
            ${CMAKE_SOURCE_DIR}/source/TbasedCompositionDiffusionStrategy.C
            ${CMAKE_SOURCE_DIR}/source/DiffusionForConcInPhaseStrategy.C
            ${CMAKE_SOURCE_DIR}/source/QuatIntegrator.C
            ${CMAKE_SOURCE_DIR}/source/BeckermannCompositionRHSStrategy.C
            ${CMAKE_SOURCE_DIR}/source/DeltaTemperatureFreeEnergyStrategy.C
            ${CMAKE_SOURCE_DIR}/source/ConstantHeatCapacityStrategy.C
            ${CMAKE_SOURCE_DIR}/source/PFModel.C
            ${CMAKE_SOURCE_DIR}/source/CartesianRobinBcHelperWithDepth.C
            ${CMAKE_SOURCE_DIR}/source/ConstantMeltingTemperatureStrategy.C
            ${CMAKE_SOURCE_DIR}/source/LinearMeltingTemperatureStrategy.C
            ${CMAKE_SOURCE_DIR}/source/PartitionCoeffRefinePatchStrategy.C
            ${CMAKE_SOURCE_DIR}/source/PartitionCoefficientStrategy.C
            ${CMAKE_SOURCE_DIR}/source/PhaseIndependentConcentrationsStrategy.C
            ${CMAKE_SOURCE_DIR}/source/AzizPartitionCoefficientStrategy.C
            ${CMAKE_SOURCE_DIR}/source/UniformPartitionCoefficientStrategy.C
            ${CMAKE_SOURCE_DIR}/source/PartitionPhaseConcentrationsStrategy.C
            ${CMAKE_SOURCE_DIR}/source/NKRHeatCapacityStrategy.C
            ${CMAKE_SOURCE_DIR}/source/PhaseFluxStrategySimple.C
            ${CMAKE_SOURCE_DIR}/source/PhaseFluxStrategyIsotropic.C
            ${CMAKE_SOURCE_DIR}/source/PhaseFluxStrategyAnisotropy.C
            ${CMAKE_SOURCE_DIR}/source/TemperatureStrategyFactory.C
            ${CMAKE_SOURCE_DIR}/source/QuatModel.C
            ${CMAKE_SOURCE_DIR}/source/QuatRefinePatchStrategy.C
            ${CMAKE_SOURCE_DIR}/source/ran2.C
            ${CMAKE_SOURCE_DIR}/source/BiasDoubleWellFreeEnergyStrategy.C
            ${CMAKE_SOURCE_DIR}/source/BiasDoubleWellBeckermannFreeEnergyStrategy.C
            ${CMAKE_SOURCE_DIR}/source/BiasDoubleWellUTRCFreeEnergyStrategy.C
            ${CMAKE_SOURCE_DIR}/source/SimpleGradStrategy.C
            ${CMAKE_SOURCE_DIR}/source/SimpleQuatGradStrategy.C
            ${CMAKE_SOURCE_DIR}/source/SimpleQuatMobilityStrategy.C
            ${CMAKE_SOURCE_DIR}/source/QuatMobilityStrategy.C
            ${CMAKE_SOURCE_DIR}/source/QuatSysSolver.C
            ${CMAKE_SOURCE_DIR}/source/QuatFACOps.C
            ${CMAKE_SOURCE_DIR}/source/QuatLevelSolver.C
            ${CMAKE_SOURCE_DIR}/source/QuatLinearRefine.C
            ${CMAKE_SOURCE_DIR}/source/QuatWeightedAverage.C
            ${CMAKE_SOURCE_DIR}/source/MinIntCoarsen.C
            ${CMAKE_SOURCE_DIR}/source/PhaseFreeEnergyStrategy.C
            ${CMAKE_SOURCE_DIR}/source/TemperatureFreeEnergyStrategy.C
            ${CMAKE_SOURCE_DIR}/source/PoissonSpecifications.C
            ${CMAKE_SOURCE_DIR}/source/CellPoissonHypreSolver.C
            ${CMAKE_SOURCE_DIR}/source/EventInterval.C
            ${CMAKE_SOURCE_DIR}/source/HBSMFreeEnergyStrategy.C
            ${CMAKE_SOURCE_DIR}/source/HBSMequilibriumPhaseConcentrationsStrategy.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyStrategyWithPenalty.C
            ${CMAKE_SOURCE_DIR}/source/CALPHADequilibriumPhaseConcentrationsStrategy.C
            ${CMAKE_SOURCE_DIR}/source/GaussianTemperatureStrategy.C
            ${CMAKE_SOURCE_DIR}/source/ScalarTemperatureStrategy.C
            ${CMAKE_SOURCE_DIR}/source/SteadyStateTemperatureStrategy.C
            ${CMAKE_SOURCE_DIR}/source/SteadyStateTemperatureGaussianSource.C
            ${CMAKE_SOURCE_DIR}/source/SteadyStateTemperatureCompositionSource.C
            ${CMAKE_SOURCE_DIR}/source/EllipticFACSolver.C
            ${CMAKE_SOURCE_DIR}/source/EllipticFACOps.C
            ${CMAKE_SOURCE_DIR}/source/PhaseConcentrationsStrategy.C
            ${CMAKE_SOURCE_DIR}/source/PhaseFACSolver.C
            ${CMAKE_SOURCE_DIR}/source/PhaseFACOps.C
            ${CMAKE_SOURCE_DIR}/source/PhaseTemperatureFACOps.C
            ${CMAKE_SOURCE_DIR}/source/EtaFACSolver.C
            ${CMAKE_SOURCE_DIR}/source/EtaFACOps.C
            ${CMAKE_SOURCE_DIR}/source/ConcFACSolver.C
            ${CMAKE_SOURCE_DIR}/source/ConcFACOps.C
            ${CMAKE_SOURCE_DIR}/source/KKSCompositionRHSStrategy.C
            ${CMAKE_SOURCE_DIR}/source/EBSCompositionRHSStrategy.C
            ${CMAKE_SOURCE_DIR}/source/CompositionRHSStrategy.C
            ${CMAKE_SOURCE_DIR}/source/Grains.C
            ${CMAKE_SOURCE_DIR}/source/TemperatureFACSolver.C
            ${CMAKE_SOURCE_DIR}/source/TemperatureFACOps.C
            ${CMAKE_SOURCE_DIR}/source/FACPreconditioner.C
            ${CMAKE_SOURCE_DIR}/source/TimeLocationIndexRobinBcCoefs.C )

set(SOURCES_FORTRAN
            ${CMAKE_SOURCE_DIR}/source/fortran/quat.f
            ${CMAKE_SOURCE_DIR}/source/fortran/functions.f )

set(M4FILES quat functions)
foreach(file ${M4FILES})
  message(STATUS "${file}...")
  m4tof(${file} ${CMAKE_SOURCE_DIR}/source/fortran)
endforeach()

set(M4FILESDIM concentrationdiffusion
              concentrationrhs
              concentrations
              correctrhsquatforsymmetry
              ellipticfacops
              grad
              hyprepoisson
              mobility
              phase
              quatcoarsen
              quatdiffs
              quatdiffusion
              quatenergy
              quatfacops
              quatgrad
              quatlevelsolver
              quatlinrefine
              quatrhs
              quatrotation
              quattag
              settozero
              smoothquat
              temperature )

set(SOURCE_FORTRAN_DIM_DIR ${CMAKE_SOURCE_DIR}/source/fortran/${NDIM}d)
foreach(file ${M4FILESDIM})
  message(STATUS ".m2 -> .f: ${file}...")
  m4tof(${file} ${SOURCE_FORTRAN_DIM_DIR})
endforeach()

set(SOURCES_FORTRAN_DIM "")
foreach(file ${M4FILESDIM})
  list(APPEND SOURCES_FORTRAN_DIM ${SOURCE_FORTRAN_DIM_DIR}/${file}.f)
endforeach()

set(MAIN ampe${NDIM}d)
add_executable(${MAIN} ${SOURCES} ${SOURCES_FORTRAN} ${SOURCES_FORTRAN_DIM})

target_link_libraries(${MAIN} ${PROJECT_LINK_LIBS_SAMRAI}
                              ${PROJECT_LINK_LIBS_HYPRE}
                              ${PROJECT_LINK_LIBS_CVODE}
                              ${PROJECT_LINK_LIBS_CPODES}
                              ${PROJECT_LINK_LIBS_NETCDF}
                              ${PROJECT_LINK_LIBS_HDF5}
                      )

add_executable(plotCALPHADbinary plotCALPHADbinary.C
                                 InterpolationType.C
                                 CALPHADEqConcSolverBinary.C
                                 CALPHADFreeEnergyFunctionsBinary.C
                                 CALPHADSpeciesPhaseGibbsEnergy.C
                                 CALPHADSpeciesPhaseGibbsEnergyExpansion.C
                                 CALPHADFunctions.C
                                 xlogx.C
                                 CALPHADConcSolverBinary.C
                                 DampedNewtonSolver.C
                                 NewtonSolver.C
                                 math_utilities.C
                                 ${CMAKE_SOURCE_DIR}/source/fortran/functions.f)
target_link_libraries(plotCALPHADbinary ${PROJECT_LINK_LIBS_SAMRAI})

add_executable(plotEquilibriumCompositions
                  plotEquilibriumCompositions.C
                  InterpolationType.C
                  CALPHADEqConcSolverBinary.C
                  CALPHADFreeEnergyFunctionsBinary.C
                  CALPHADSpeciesPhaseGibbsEnergy.C
                  CALPHADSpeciesPhaseGibbsEnergyExpansion.C
                  CALPHADFunctions.C
                  xlogx.C
                  CALPHADConcSolverBinary.C
                  DampedNewtonSolver.C
                  NewtonSolver.C
                  math_utilities.C
                  ${CMAKE_SOURCE_DIR}/source/fortran/functions.f)
target_link_libraries(plotEquilibriumCompositions ${PROJECT_LINK_LIBS_SAMRAI}
                                                  ${PROJECT_LINK_LIBS_HDF5})

add_executable(plotEquilibriumCompositionsTernary
                  plotEquilibriumCompositionsTernary.C
                  InterpolationType.C
                  CALPHADEqConcSolverTernary.C
                  CALPHADEqPhaseConcSolverTernary.C
                  CALPHADFreeEnergyFunctionsTernary.C
                  CALPHADSpeciesPhaseGibbsEnergy.C
                  CALPHADSpeciesPhaseGibbsEnergyExpansion.C
                  CALPHADFunctions.C
                  xlogx.C
                  CALPHADConcSolverTernary.C
                  DampedNewtonSolver.C
                  NewtonSolver.C
                  math_utilities.C
                  ${CMAKE_SOURCE_DIR}/source/fortran/functions.f)
target_link_libraries(plotEquilibriumCompositionsTernary
                         ${PROJECT_LINK_LIBS_SAMRAI}
                         ${PROJECT_LINK_LIBS_HDF5})

add_executable(computeCALPHADbinaryEquilibrium
               ${CMAKE_SOURCE_DIR}/source/xlogx.C
               ${CMAKE_SOURCE_DIR}/source/InterpolationType.C
               computeCALPHADbinaryEquilibrium.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyFunctionsBinary.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADEqConcSolverBinary.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADConcSolverBinary.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADSpeciesPhaseGibbsEnergy.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADSpeciesPhaseGibbsEnergyExpansion.C
               ${CMAKE_SOURCE_DIR}/source/CALPHADFunctions.C
               ${CMAKE_SOURCE_DIR}/source/xlogx.C
               ${CMAKE_SOURCE_DIR}/source/NewtonSolver.C
               ${CMAKE_SOURCE_DIR}/source/DampedNewtonSolver.C
               ${CMAKE_SOURCE_DIR}/source/math_utilities.C
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               )
target_link_libraries(computeCALPHADbinaryEquilibrium
                      ${PROJECT_LINK_LIBS_SAMRAI})
