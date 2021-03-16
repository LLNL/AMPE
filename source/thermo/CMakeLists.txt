set(SOURCES
    NewtonSolver.C
    DampedNewtonSolver.C
    CALPHADFunctions.C
    CALPHADEqConcSolverBinary.C
    CALPHADEqConcSolverTernary.C
    CALPHADConcSolverBinary.C
    CALPHADConcSolverTernary.C
    CALPHADTieLineConcSolverTernary.C
    CALPHADFreeEnergyFunctionsBinary.C
    CALPHADFreeEnergyFunctionsTernary.C
    CALPHADSpeciesPhaseGibbsEnergy.C
    CALPHADSpeciesPhaseGibbsEnergyExpansion.C
    KKSdiluteBinaryConcentrationSolver.C
    KKSFreeEnergyFunctionDiluteBinary.C
    InterpolationType.C
    xlogx.C
    math_utilities.C
)

add_library(ampe_thermo ${SOURCES})

target_link_libraries(ampe_thermo ${PROJECT_LINK_LIBS_SAMRAI})

install(TARGETS ampe_thermo DESTINATION lib)
