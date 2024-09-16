// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.

#ifndef included_PhaseConcentrationsStrategyFactory
#define included_PhaseConcentrationsStrategyFactory

#include "CALPHADequilibriumPhaseConcentrationsStrategyThreePhaseStochioB.h"
#include "CALPHADequilibriumPhaseConcentrationsStrategy.h"
#include "KKSdiluteEquilibriumPhaseConcentrationsStrategy.h"
#include "QuadraticEquilibriumPhaseConcentrationsStrategy.h"
#include "QuadraticEquilibriumPhaseConcentrationsStrategyMultiOrder.h"
#include "QuadraticEquilibriumThreePhasesTernaryStrategyMultiOrder.h"
#include "PartitionPhaseConcentrationsStrategy.h"
#include "PhaseIndependentConcentrationsStrategy.h"
#include "Database2JSON.h"
#include "CALPHADFreeEnergyFunctionsBinaryThreePhase.h"
#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"
#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"
#include "CALPHADequilibriumPhaseConcentrationsStrategyMultiOrder.h"
#include "CALPHADequilibriumPhaseConcentrationsStrategyMultiOrderThreePhases.h"
#include "CALPHADFunctions.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

class PhaseConcentrationsStrategyFactory
{
 public:
   static std::shared_ptr<PhaseConcentrationsStrategy> create(
       QuatModelParameters& model_parameters, int conc_l_scratch_id,
       int conc_a_scratch_id, int conc_b_scratch_id, int conc_l_ref_id,
       int conc_a_ref_id, int conc_b_ref_id, int partition_coeff_id,
       const int ncompositions, std::shared_ptr<tbox::Database> conc_db,
       std::shared_ptr<tbox::Database> newton_db)
   {
      boost::property_tree::ptree calphad_pt;
      std::shared_ptr<tbox::MemoryDatabase> calphad_db;

      tbox::plog << "=== PhaseConcentrationsStrategyFactory ===" << std::endl;
      if (model_parameters.isConcentrationModelCALPHAD()) {
         std::shared_ptr<tbox::Database> db(conc_db->getDatabase("Calphad"));
         std::string calphad_filename = db->getString("filename");
         tbox::plog << "Read file " << calphad_filename << std::endl;
         if (calphad_filename.compare(calphad_filename.size() - 4, 4, "json") ==
             0) {
            boost::property_tree::read_json(calphad_filename, calphad_pt);
         } else {
            calphad_db.reset(new tbox::MemoryDatabase("calphad_db"));
            tbox::InputManager::getManager()->parseInputFile(calphad_filename,
                                                             calphad_db);
            copyDatabase(calphad_db, calphad_pt);
         }
      }

      std::shared_ptr<PhaseConcentrationsStrategy> phase_conc_strategy;

      if (model_parameters.kks_phase_concentration()) {
         tbox::plog << "Phase concentration determined by KKS" << std::endl;
         if (model_parameters.isConcentrationModelCALPHAD()) {
            tbox::plog << "CALPHAD..." << std::endl;
            if (ncompositions == 1) {
               tbox::plog << "Binary alloy..." << std::endl;
               bool subl = Thermo4PFM::checkSublattice(calphad_pt);
               double cB = model_parameters.getStochioB();
               if (conc_b_scratch_id >= 0) {
                  tbox::plog << "Three phases..." << std::endl;
                  // three phases
                  if (model_parameters.withMultipleOrderP()) {
                     tbox::plog << "Multi-order parameters..." << std::endl;
                     // multi-order parameters model
                     if (subl) {
                        phase_conc_strategy.reset(
                            new CALPHADequilibriumPhaseConcentrationsStrategyMultiOrderThreePhases<
                                Thermo4PFM::
                                    CALPHADFreeEnergyFunctionsBinary3Ph2Sl>(
                                model_parameters.norderpA(), conc_l_scratch_id,
                                conc_a_scratch_id, conc_b_scratch_id,
                                model_parameters, conc_db, newton_db));
                     } else {
                        phase_conc_strategy.reset(
                            new CALPHADequilibriumPhaseConcentrationsStrategyMultiOrderThreePhases<
                                Thermo4PFM::
                                    CALPHADFreeEnergyFunctionsBinaryThreePhase>(
                                model_parameters.norderpA(), conc_l_scratch_id,
                                conc_a_scratch_id, conc_b_scratch_id,
                                model_parameters, conc_db, newton_db));
                     }
                  } else {
                     // three phases model
                     if (cB >= 0.) {
                        phase_conc_strategy.reset(
                            new CALPHADequilibriumPhaseConcentrationsStrategyThreePhaseStochioB(
                                cB, conc_l_scratch_id, conc_a_scratch_id,
                                conc_b_scratch_id, conc_l_ref_id, conc_a_ref_id,
                                conc_b_ref_id, calphad_pt, newton_db));
                     } else if (subl) {
                        phase_conc_strategy.reset(
                            new CALPHADequilibriumPhaseConcentrationsStrategy<
                                Thermo4PFM::
                                    CALPHADFreeEnergyFunctionsBinary3Ph2Sl>(
                                conc_l_scratch_id, conc_a_scratch_id,
                                conc_b_scratch_id, conc_l_ref_id, conc_a_ref_id,
                                conc_b_ref_id,
                                model_parameters.energy_interp_func_type(),
                                model_parameters.conc_interp_func_type(),
                                model_parameters.with_third_phase(), calphad_pt,
                                newton_db, ncompositions));
                     } else {
                        phase_conc_strategy.reset(
                            new CALPHADequilibriumPhaseConcentrationsStrategy<
                                Thermo4PFM::
                                    CALPHADFreeEnergyFunctionsBinaryThreePhase>(
                                conc_l_scratch_id, conc_a_scratch_id,
                                conc_b_scratch_id, conc_l_ref_id, conc_a_ref_id,
                                conc_b_ref_id,
                                model_parameters.energy_interp_func_type(),
                                model_parameters.conc_interp_func_type(),
                                model_parameters.with_third_phase(), calphad_pt,
                                newton_db, ncompositions));
                     }
                  }
               } else {  // two phases only
                  tbox::plog << "Two phases..." << std::endl;
                  if (model_parameters.norderp() > 1) {
                     // multiple order parameters
                     tbox::plog << "Multi-order parameters..." << std::endl;
                     phase_conc_strategy.reset(
                         new CALPHADequilibriumPhaseConcentrationsStrategyMultiOrder(
                             conc_l_scratch_id, conc_a_scratch_id,
                             model_parameters, conc_db, newton_db));
                  } else {  // single order parameter
                     if (subl) {
                        phase_conc_strategy.reset(
                            new CALPHADequilibriumPhaseConcentrationsStrategy<
                                Thermo4PFM::
                                    CALPHADFreeEnergyFunctionsBinary2Ph1Sl>(
                                conc_l_scratch_id, conc_a_scratch_id,
                                conc_b_scratch_id, conc_l_ref_id, conc_a_ref_id,
                                conc_b_ref_id,
                                model_parameters.energy_interp_func_type(),
                                model_parameters.conc_interp_func_type(),
                                model_parameters.with_third_phase(), calphad_pt,
                                newton_db, ncompositions));

                     } else {
                        phase_conc_strategy.reset(
                            new CALPHADequilibriumPhaseConcentrationsStrategy<
                                Thermo4PFM::CALPHADFreeEnergyFunctionsBinary>(
                                conc_l_scratch_id, conc_a_scratch_id,
                                conc_b_scratch_id, conc_l_ref_id, conc_a_ref_id,
                                conc_b_ref_id,
                                model_parameters.energy_interp_func_type(),
                                model_parameters.conc_interp_func_type(),
                                model_parameters.with_third_phase(), calphad_pt,
                                newton_db, ncompositions));
                     }
                  }
               }
            } else if (ncompositions == 2) {
               tbox::plog << "Ternary..." << std::endl;
               phase_conc_strategy.reset(
                   new CALPHADequilibriumPhaseConcentrationsStrategy<
                       Thermo4PFM::CALPHADFreeEnergyFunctionsTernary>(
                       conc_l_scratch_id, conc_a_scratch_id, conc_b_scratch_id,
                       conc_l_ref_id, conc_a_ref_id, conc_b_ref_id,
                       model_parameters.energy_interp_func_type(),
                       model_parameters.conc_interp_func_type(),
                       model_parameters.with_third_phase(), calphad_pt,
                       newton_db, ncompositions));
            }
         } else if (model_parameters.isConcentrationModelKKSdilute()) {
            tbox::plog << "Dilute..." << std::endl;
            phase_conc_strategy.reset(
                new KKSdiluteEquilibriumPhaseConcentrationsStrategy(
                    conc_l_scratch_id, conc_a_scratch_id, conc_b_scratch_id,
                    conc_l_ref_id, conc_a_ref_id, conc_b_ref_id,
                    model_parameters.energy_interp_func_type(),
                    model_parameters.conc_interp_func_type(), conc_db));
         } else {
            if (model_parameters.isConcentrationModelQuadratic()) {
               if (model_parameters.norderp() > 1) {
                  tbox::plog << "Quadratic, MultiOrder..." << std::endl;
                  if (conc_b_scratch_id >= 0) {
                     tbox::plog << "Quadratic, MultiOrder, Three phases..."
                                << std::endl;
                     phase_conc_strategy.reset(
                         new QuadraticEquilibriumThreePhasesTernaryStrategyMultiOrder(
                             model_parameters.norderpA(), conc_l_scratch_id,
                             conc_a_scratch_id, conc_b_scratch_id,
                             model_parameters, conc_db));
                  } else {
                     phase_conc_strategy.reset(
                         new QuadraticEquilibriumPhaseConcentrationsStrategyMultiOrder(
                             conc_l_scratch_id, conc_a_scratch_id,
                             model_parameters, conc_db));
                  }
               } else {
                  tbox::plog << "Quadratic..." << std::endl;
                  assert(conc_b_scratch_id == -1);
                  phase_conc_strategy.reset(
                      new QuadraticEquilibriumPhaseConcentrationsStrategy(
                          conc_l_scratch_id, conc_a_scratch_id,
                          model_parameters, conc_db));
               }
            }
         }
      } else {
         if (model_parameters.partition_phase_concentration()) {
            phase_conc_strategy.reset(new PartitionPhaseConcentrationsStrategy(
                conc_l_scratch_id, conc_a_scratch_id, conc_b_scratch_id,
                model_parameters.conc_interp_func_type(), partition_coeff_id));
         } else {  // simply use cl=ca=c
            phase_conc_strategy.reset(
                new PhaseIndependentConcentrationsStrategy(conc_l_scratch_id,
                                                           conc_a_scratch_id,
                                                           conc_b_scratch_id));
         }
      }
      return phase_conc_strategy;
   }
};

#endif
