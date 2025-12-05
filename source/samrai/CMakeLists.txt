set(SOURCES
    CellPoissonHypreSolver.cc
    PoissonSpecifications.cc
    CVODEAbstractFunctions.cc
    CVODESolver.cc
    ARKODESolver.cc
    FACPreconditioner.cc)

add_library(ampe_samrai ${SOURCES})

target_link_libraries(ampe_samrai ${PROJECT_LINK_LIBS_SAMRAI} ${PROJECT_LINK_LIBS_CVODE})
