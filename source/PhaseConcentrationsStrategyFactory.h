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

#include "CALPHADequilibriumPhaseConcentrationsStrategy.h"
#include "KKSdiluteEquilibriumPhaseConcentrationsStrategy.h"
#include "QuadraticEquilibriumPhaseConcentrationsStrategy.h"
#include "QuadraticEquilibriumPhaseConcentrationsStrategyMultiOrder.h"
#include "PartitionPhaseConcentrationsStrategy.h"
#include "PhaseIndependentConcentrationsStrategy.h"
#include "Database2JSON.h"
#include "CALPHADFreeEnergyFunctionsBinaryThreePhase.h"
#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"
#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"
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
               bool subl = checkSublattice(calphad_pt);
               if (conc_b_scratch_id >= 0) {
                  if (subl) {
                     phase_conc_strategy.reset(
                         new CALPHADequilibriumPhaseConcentrationsStrategy<
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
                             CALPHADFreeEnergyFunctionsBinaryThreePhase>(
                             conc_l_scratch_id, conc_a_scratch_id,
                             conc_b_scratch_id, conc_l_ref_id, conc_a_ref_id,
                             conc_b_ref_id,
                             model_parameters.energy_interp_func_type(),
                             model_parameters.conc_interp_func_type(),
                             model_parameters.with_third_phase(), calphad_pt,
                             newton_db, ncompositions));
                  }
               } else if (subl) {
                  phase_conc_strategy.reset(
                      new CALPHADequilibriumPhaseConcentrationsStrategy<
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
                          CALPHADFreeEnergyFunctionsBinary>(
                          conc_l_scratch_id, conc_a_scratch_id,
                          conc_b_scratch_id, conc_l_ref_id, conc_a_ref_id,
                          conc_b_ref_id,
                          model_parameters.energy_interp_func_type(),
                          model_parameters.conc_interp_func_type(),
                          model_parameters.with_third_phase(), calphad_pt,
                          newton_db, ncompositions));
               }
            } else if (ncompositions == 2) {
               tbox::plog << "Ternary..." << std::endl;
               phase_conc_strategy.reset(
                   new CALPHADequilibriumPhaseConcentrationsStrategy<
                       CALPHADFreeEnergyFunctionsTernary>(
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
            if (model_parameters.isConcentrationModelQuadratic())
               if (model_parameters.norderp() > 1) {
                  tbox::plog << "Quadratic, MultiOrder..." << std::endl;
                  phase_conc_strategy.reset(
                      new QuadraticEquilibriumPhaseConcentrationsStrategyMultiOrder(
                          conc_l_scratch_id, conc_a_scratch_id,
                          model_parameters, conc_db));
               } else {
                  tbox::plog << "Quadratic..." << std::endl;
                  assert(conc_b_scratch_id == -1);
                  phase_conc_strategy.reset(
                      new QuadraticEquilibriumPhaseConcentrationsStrategy(
                          conc_l_scratch_id, conc_a_scratch_id,
                          model_parameters, conc_db));
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
