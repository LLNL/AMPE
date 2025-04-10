include_directories( ${CMAKE_SOURCE_DIR}/source )

set(THERMOLIB ${PROJECT_LINK_LIBS_THERMO4PFM})

add_executable(plotCALPHADbinary
               plotCALPHADbinary.cc
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               ${CMAKE_SOURCE_DIR}/source/Database2JSON.cc)
target_link_libraries(plotCALPHADbinary ${PROJECT_LINK_LIBS_SAMRAI}
                                        ${THERMOLIB})

add_executable(plotEnergyVsCompositionParabolicThreePhase
               ${CMAKE_SOURCE_DIR}/source/ParabolicTools.cc
               ${CMAKE_SOURCE_DIR}/drivers/plotEnergyVsCompositionParabolicThreePhase.cc )
target_link_libraries(plotEnergyVsCompositionParabolicThreePhase ${PROJECT_LINK_LIBS_SAMRAI}
                                        ${THERMOLIB})

add_executable(computeEquilibriumCompositions
               computeEquilibriumCompositions.cc
               ${CMAKE_SOURCE_DIR}/source/Database2JSON.cc
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f)
target_link_libraries(computeEquilibriumCompositions ${PROJECT_LINK_LIBS_SAMRAI}
                                                  ${PROJECT_LINK_LIBS_HDF5}
                                                  ${THERMOLIB})

add_executable(computeEquilibriumCompositionsParabolic
               ${CMAKE_SOURCE_DIR}/source/ParabolicTools.cc
               computeEquilibriumCompositionsParabolic.cc)
target_link_libraries(computeEquilibriumCompositionsParabolic ${PROJECT_LINK_LIBS_SAMRAI}
                                                  ${PROJECT_LINK_LIBS_HDF5}
                                                  ${THERMOLIB})

add_executable(computeEquilibriumCompositionsTernary
               computeEquilibriumCompositionsTernary.cc
               ${CMAKE_SOURCE_DIR}/source/Database2JSON.cc
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f)
target_link_libraries(computeEquilibriumCompositionsTernary
                         ${PROJECT_LINK_LIBS_SAMRAI}
                         ${PROJECT_LINK_LIBS_HDF5}
                         ${THERMOLIB})

add_executable(computeCALPHADbinaryEquilibrium
               computeCALPHADbinaryEquilibrium.cc
               ${CMAKE_SOURCE_DIR}/source/Database2JSON.cc
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f
               )
target_link_libraries(computeCALPHADbinaryEquilibrium
                      ${PROJECT_LINK_LIBS_SAMRAI}
                      ${THERMOLIB})

add_executable(computeParabolic3PhasesEquilibrium
               computeParabolic3PhasesEquilibrium.cc
               ${CMAKE_SOURCE_DIR}/source/Database2JSON.cc
               ${CMAKE_SOURCE_DIR}/source/ParabolicTools.cc
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f)
target_link_libraries(computeParabolic3PhasesEquilibrium
                      ${PROJECT_LINK_LIBS_SAMRAI}
                      ${THERMOLIB})

add_executable(computeParabolic3PhasesKKS
               computeParabolic3PhasesKKS.cc
               ${CMAKE_SOURCE_DIR}/source/Database2JSON.cc
               ${CMAKE_SOURCE_DIR}/source/ParabolicTools.cc
               ${CMAKE_SOURCE_DIR}/source/fortran/functions.f)
target_link_libraries(computeParabolic3PhasesKKS
                      ${PROJECT_LINK_LIBS_SAMRAI}
                      ${THERMOLIB})

add_executable(convertDatabase2JSON
               ${CMAKE_SOURCE_DIR}/source/Database2JSON.cc
               convertDatabase2JSON.cc )
target_link_libraries(convertDatabase2JSON
                      ${PROJECT_LINK_LIBS_SAMRAI})

add_executable(computeMobilityZetaFactor
               ${CMAKE_SOURCE_DIR}/source/Database2JSON.cc
               ${CMAKE_SOURCE_DIR}/source/KKStools.cc
               computeMobilityZetaFactor.cc)
target_link_libraries(computeMobilityZetaFactor
                      ${PROJECT_LINK_LIBS_SAMRAI}
                      ${THERMOLIB})
