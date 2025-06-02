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

set(SOURCES_WPENALTY
    ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyFunctionsWithPenaltyBinary.cc
    ${CMAKE_SOURCE_DIR}/source/CALPHADEqConcSolverBinaryWithPenalty.cc
    ${CMAKE_SOURCE_DIR}/source/CALPHADConcSolverBinaryWithPenalty.cc
    ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyStrategyWithPenalty.cc
    )

set(SOURCES ${CMAKE_SOURCE_DIR}/source/KKStools.cc
            ${CMAKE_SOURCE_DIR}/source/ParabolicTools.cc
            ${CMAKE_SOURCE_DIR}/source/MultiOrderBinaryDrivingForce.cc
            ${CMAKE_SOURCE_DIR}/source/MultiOrderBinaryThreePhasesDrivingForce.cc
	    ${CMAKE_SOURCE_DIR}/source/MultiOrderBinaryThreePhasesDrivingForceStochioB.cc
            ${CMAKE_SOURCE_DIR}/source/CahnHilliardDoubleWell.cc
            ${CMAKE_SOURCE_DIR}/source/SinteringUWangRHSStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/WangSinteringCompositionRHSStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/WangSinteringDiffusion.cc
            ${CMAKE_SOURCE_DIR}/source/TemperatureDependentWangSinteringDiffusion.cc
            ${CMAKE_SOURCE_DIR}/source/main.cc
            ${CMAKE_SOURCE_DIR}/source/QuatGradModulusStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/ScalarCompositionDiffusionStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/HierarchyStencilOps.cc
            ${CMAKE_SOURCE_DIR}/source/TemperatureHistory.cc
            ${CMAKE_SOURCE_DIR}/source/QuatFaceCoeff.cc
            ${CMAKE_SOURCE_DIR}/source/AdaptMovingFrame.cc
            ${CMAKE_SOURCE_DIR}/source/computeQDiffs.cc
            ${CMAKE_SOURCE_DIR}/source/ThreeArgsInterpolationType.cc
            ${CMAKE_SOURCE_DIR}/source/computeWeights3Pairs.cc
            ${CMAKE_SOURCE_DIR}/source/KimMobilityStrategyInfMob3Phases.cc
            ${CMAKE_SOURCE_DIR}/source/EnergyEvaluationStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/TwoPhasesEnergyEvaluationStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/MultiOrderPEnergyEvaluationStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/RigidBodyForces.cc
            ${CMAKE_SOURCE_DIR}/source/PhaseRigidBodyForces.cc
            ${CMAKE_SOURCE_DIR}/source/WangRigidBodyForces.cc
            ${CMAKE_SOURCE_DIR}/source/TiltingFolchPlapp2005.cc
            ${CMAKE_SOURCE_DIR}/source/TiltingMoelans2011.cc
            ${CMAKE_SOURCE_DIR}/source/TemperatureRHSStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/SimpleTemperatureRHSStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/FolchPlappRHSStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/MultiOrderRHSStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/MovingFrameRHS.cc
            ${CMAKE_SOURCE_DIR}/source/PhaseRHSStrategyWithQ.cc
            ${CMAKE_SOURCE_DIR}/source/diagnostics.cc
            ${CMAKE_SOURCE_DIR}/source/EquilibriumPhaseConcentrationsBinary.cc
            ${CMAKE_SOURCE_DIR}/source/KKSdiluteEquilibriumPhaseConcentrationsStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/KKSdiluteBinary.cc
            ${CMAKE_SOURCE_DIR}/source/QuatModelParameters.cc
            ${CMAKE_SOURCE_DIR}/source/CompositionStrategyMobilities.cc
            ${CMAKE_SOURCE_DIR}/source/CALPHADMobility.cc
            ${CMAKE_SOURCE_DIR}/source/FreeEnergyStrategyBinary.cc
            ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyStrategyBinary.cc
            ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyBinaryMultiOrder.cc
            ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyBinaryMultiOrderThreePhases.cc
            ${CMAKE_SOURCE_DIR}/source/FreeEnergyStrategyThreePhase.cc
            ${CMAKE_SOURCE_DIR}/source/ParabolicFreeEnergyBinary.cc
            ${CMAKE_SOURCE_DIR}/source/ParabolicFreeEnergyMultiOrderBinary.cc
            ${CMAKE_SOURCE_DIR}/source/QuadraticFreeEnergyBinary.cc
            ${CMAKE_SOURCE_DIR}/source/QuadraticFreeEnergyMultiOrderBinary.cc
            ${CMAKE_SOURCE_DIR}/source/QuadraticFreeEnergyMultiOrderBinaryThreePhase.cc
            ${CMAKE_SOURCE_DIR}/source/QuadraticFreeEnergyMultiOrderTernaryThreePhase.cc
            ${CMAKE_SOURCE_DIR}/source/ParabolicFreeEnergyMultiOrderBinaryThreePhase.cc
            ${CMAKE_SOURCE_DIR}/source/EquilibriumPhaseConcentrationsBinaryMultiOrder.cc
            ${CMAKE_SOURCE_DIR}/source/EquilibriumPhaseConcentrationsBinaryMultiOrderThreePhases.cc
            ${CMAKE_SOURCE_DIR}/source/CALPHADequilibriumPhaseConcentrationsStrategyMultiOrder.cc
            ${CMAKE_SOURCE_DIR}/source/CALPHADequilibriumPhaseConcentrationsStrategyMultiOrderThreePhases.cc
            ${CMAKE_SOURCE_DIR}/source/CALPHADequilibriumPhaseConcentrationsThreePhasesStochioB.cc
            ${CMAKE_SOURCE_DIR}/source/CALPHADequilibriumPhaseConcentrationsMultiOrderThreePhasesStochioB.cc
            ${CMAKE_SOURCE_DIR}/source/QuadraticEquilibriumPhaseConcentrationsBinaryMultiOrder.cc
            ${CMAKE_SOURCE_DIR}/source/QuadraticEquilibriumThreePhasesTernaryMultiOrder.cc
            ${CMAKE_SOURCE_DIR}/source/QuadraticEquilibriumPhaseConcentrationsBinary.cc
            ${CMAKE_SOURCE_DIR}/source/ParabolicEquilibriumPhaseConcentrationsBinary.cc
            ${CMAKE_SOURCE_DIR}/source/ParabolicEquilibriumPhaseConcentrationsBinaryMultiOrder.cc
            ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyStrategyBinaryFolchPlapp.cc
            ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyStrategyBinaryFolchPlappStochioB.cc
            ${CMAKE_SOURCE_DIR}/source/ParabolicEquilibriumThreePhasesBinaryMultiOrder.cc
            ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyBinaryMultiOrderThreePhasesStochioB.cc
            ${CMAKE_SOURCE_DIR}/source/CALPHADFreeEnergyStrategyTernary.cc
            ${CMAKE_SOURCE_DIR}/source/tools.cc
            ${CMAKE_SOURCE_DIR}/source/FreeEnergyStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/ConcFreeEnergyStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/AMPE.cc
            ${CMAKE_SOURCE_DIR}/source/Noise.cc
            ${CMAKE_SOURCE_DIR}/source/NormalNoise.cc
            ${CMAKE_SOURCE_DIR}/source/UniformNoise.cc
            ${CMAKE_SOURCE_DIR}/source/KimMobilityStrategyInfMob.cc
            ${CMAKE_SOURCE_DIR}/source/KimMobilityStrategyFiniteMob.cc
            ${CMAKE_SOURCE_DIR}/source/KimMobilityStrategyFiniteMobAntiTrap.cc
            ${CMAKE_SOURCE_DIR}/source/KimMobilityStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/FieldsInitializer.cc
            ${CMAKE_SOURCE_DIR}/source/FieldsWriter.cc
            ${CMAKE_SOURCE_DIR}/source/GradientTemperatureStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/toolsSAMRAI.cc
            ${CMAKE_SOURCE_DIR}/source/TbasedCompositionDiffusionStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/MobilityCompositionDiffusionStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/DerivDiffusionCoeffForQuat.cc
            ${CMAKE_SOURCE_DIR}/source/QuatIntegrator.cc
            ${CMAKE_SOURCE_DIR}/source/BeckermannCompositionRHSStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/DeltaTemperatureFreeEnergyStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/ConstantHeatCapacityStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/PFModel.cc
            ${CMAKE_SOURCE_DIR}/source/CartesianRobinBcHelperWithDepth.cc
            ${CMAKE_SOURCE_DIR}/source/ConstantMeltingTemperatureStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/LinearMeltingTemperatureStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/PartitionCoeffRefinePatchStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/PartitionCoefficientStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/PhaseIndependentConcentrationsStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/AzizPartitionCoefficientStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/UniformPartitionCoefficientStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/PartitionPhaseConcentrationsStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/NKRHeatCapacityStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/PhaseFluxStrategySimple.cc
            ${CMAKE_SOURCE_DIR}/source/PhaseFluxStrategyIsotropic.cc
            ${CMAKE_SOURCE_DIR}/source/PhaseFluxStrategyAnisotropy.cc
            ${CMAKE_SOURCE_DIR}/source/TemperatureStrategyFactory.cc
            ${CMAKE_SOURCE_DIR}/source/QuatModel.cc
            ${CMAKE_SOURCE_DIR}/source/QuatRefinePatchStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/BiasDoubleWellFreeEnergyStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/BiasDoubleWellBeckermannFreeEnergyStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/BiasDoubleWellUTRCFreeEnergyStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/SimpleGradStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/SimpleQuatGradStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/SimpleQuatMobilityStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/QuatMobilityStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/QuatSysSolver.cc
            ${CMAKE_SOURCE_DIR}/source/QuatFACOps.cc
            ${CMAKE_SOURCE_DIR}/source/QuatLevelSolver.cc
            ${CMAKE_SOURCE_DIR}/source/QuatLinearRefine.cc
            ${CMAKE_SOURCE_DIR}/source/QuatWeightedAverage.cc
            ${CMAKE_SOURCE_DIR}/source/MinIntCoarsen.cc
            ${CMAKE_SOURCE_DIR}/source/PhaseFreeEnergyStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/TemperatureFreeEnergyStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/EventInterval.cc
            ${CMAKE_SOURCE_DIR}/source/EquilibriumPhaseConcentrationsThreePhases.cc
            ${CMAKE_SOURCE_DIR}/source/GaussianTemperatureStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/ScalarTemperatureStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/SteadyStateTemperatureStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/SteadyStateTemperatureGaussianSource.cc
            ${CMAKE_SOURCE_DIR}/source/SteadyStateTemperatureCompositionSource.cc
            ${CMAKE_SOURCE_DIR}/source/EllipticFACSolver.cc
            ${CMAKE_SOURCE_DIR}/source/EllipticFACOps.cc
            ${CMAKE_SOURCE_DIR}/source/PhaseConcentrationsStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/PhaseFACSolver.cc
            ${CMAKE_SOURCE_DIR}/source/PhaseFACOps.cc
            ${CMAKE_SOURCE_DIR}/source/PhaseTemperatureFACOps.cc
            ${CMAKE_SOURCE_DIR}/source/EtaFACSolver.cc
            ${CMAKE_SOURCE_DIR}/source/EtaFACOps.cc
            ${CMAKE_SOURCE_DIR}/source/ConcFACSolver.cc
            ${CMAKE_SOURCE_DIR}/source/ConcFACOps.cc
            ${CMAKE_SOURCE_DIR}/source/KKSCompositionRHSStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/EBSCompositionRHSStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/EBSCompositionRHSStrategyWithGradT.cc
            ${CMAKE_SOURCE_DIR}/source/CompositionRHSStrategy.cc
            ${CMAKE_SOURCE_DIR}/source/Grains.cc
            ${CMAKE_SOURCE_DIR}/source/RigidBodyMotionRHS.cc
            ${CMAKE_SOURCE_DIR}/source/RigidBodyMotionConcRHS.cc
            ${CMAKE_SOURCE_DIR}/source/TemperatureFACSolver.cc
            ${CMAKE_SOURCE_DIR}/source/TemperatureFACOps.cc
            ${CMAKE_SOURCE_DIR}/source/TimeLocationIndexRobinBcCoefs.cc )

set(SOURCES_Database2JSON
    ${CMAKE_SOURCE_DIR}/source/Database2JSON.cc )

set(SOURCES_FORTRAN
            ${CMAKE_SOURCE_DIR}/source/fortran/quat.f
            ${CMAKE_SOURCE_DIR}/source/fortran/functions.f )

set(M4FILESDIM concentrationdiffusion
              flux
              concentrationrhs
              concentrations
              correctrhsquatforsymmetry
              ellipticfacops
              forces
              grad
              hyprepoisson
              stencils
              mobility
              multiorder
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
              rbmotion
              settozero
              side2cell
              sintering
              smoothquat
              temperature )

set(SOURCE_FORTRAN_DIM_DIR ${CMAKE_SOURCE_DIR}/source/fortran/${NDIM}d)
foreach(file ${M4FILESDIM})
  message(STATUS ".m4 -> .f: ${file}...")
  m4tof(${file} ${SOURCE_FORTRAN_DIM_DIR})
endforeach()

set(SOURCES_FORTRAN_DIM "")
foreach(file ${M4FILESDIM})
  list(APPEND SOURCES_FORTRAN_DIM ${SOURCE_FORTRAN_DIM_DIR}/${file}.f)
endforeach()

set(MAIN ampe${NDIM}d)
add_executable(${MAIN} ${SOURCES} ${SOURCES_FORTRAN}
                       ${SOURCES_FORTRAN_DIM}
                       ${SOURCES_Database2JSON})

set(THERMOLIB ${PROJECT_LINK_LIBS_THERMO4PFM})

target_link_libraries(${MAIN} ${THERMOLIB})

add_subdirectory(samrai)
set(SAMRAILIB ${internal_samrai})

target_link_libraries(${MAIN} ${SAMRAILIB})

target_link_libraries(${MAIN} ${PROJECT_LINK_LIBS_SAMRAI}
                              ${PROJECT_LINK_LIBS_HYPRE}
                              ${PROJECT_LINK_LIBS_CVODE}
                              ${PROJECT_LINK_LIBS_ARKODE}
                              ${PROJECT_LINK_LIBS_NETCDF}
                              ${PROJECT_LINK_LIBS_HDF5}
                      )

install(TARGETS ${MAIN} DESTINATION bin)
