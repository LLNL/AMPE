set(SOURCES
    NewtonSolver.cc
    DampedNewtonSolver.cc
    CALPHADFunctions.cc
    CALPHADEqConcSolverBinary.cc
    CALPHADEqConcSolverTernary.cc
    CALPHADConcSolverBinary.cc
    CALPHADConcSolverTernary.cc
    CALPHADTieLineConcSolverTernary.cc
    CALPHADFreeEnergyFunctionsBinary.cc
    CALPHADFreeEnergyFunctionsTernary.cc
    CALPHADSpeciesPhaseGibbsEnergy.cc
    CALPHADSpeciesPhaseGibbsEnergyExpansion.cc
    KKSdiluteBinaryConcentrationSolver.cc
    KKSFreeEnergyFunctionDiluteBinary.cc
    InterpolationType.cc
    xlogx.cc
    math_utilities.cc
)

add_library(ampe_thermo ${SOURCES})

target_link_libraries(ampe_thermo ${PROJECT_LINK_LIBS_SAMRAI})

install(TARGETS ampe_thermo DESTINATION lib)
