set(SOURCES
    CellPoissonHypreSolver.C
    PoissonSpecifications.C
    SundialsAbstractVector.C
    Sundials_SAMRAIVector.C
    CVODEAbstractFunctions.C
    CVODESolver.C
    CPODESAbstractFunctions.C
    CPODESSolver.C
    FACPreconditioner.C)

add_library(ampe_samrai ${SOURCES})

target_link_libraries(ampe_samrai ${PROJECT_LINK_LIBS_SAMRAI} ${PROJECT_LINK_LIBS_CVODE})

install(TARGETS ampe_samrai DESTINATION lib)
