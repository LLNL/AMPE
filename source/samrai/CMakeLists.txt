set(SOURCES
    CellPoissonHypreSolver.cc
    PoissonSpecifications.cc
    SundialsAbstractVector.cc
    Sundials_SAMRAIVector.cc
    CVODEAbstractFunctions.cc
    CVODESolver.cc
    FACPreconditioner.cc)

add_library(ampe_samrai ${SOURCES})

target_link_libraries(ampe_samrai ${PROJECT_LINK_LIBS_SAMRAI} ${PROJECT_LINK_LIBS_CVODE})
