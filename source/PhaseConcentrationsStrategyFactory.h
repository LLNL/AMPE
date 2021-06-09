// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.

#ifndef included_PhaseConcentrationsStrategyFactory
#define included_PhaseConcentrationsStrategyFactory

#include "CALPHADequilibriumPhaseConcentrationsStrategy.h"
#include "KKSdiluteEquilibriumPhaseConcentrationsStrategy.h"
#include "HBSMequilibriumPhaseConcentrationsStrategy.h"
#include "PartitionPhaseConcentrationsStrategy.h"

class PhaseConcentrationsStrategyFactory
{
 public:
   static std::shared_ptr<PhaseConcentrationsStrategy> create(
       QuatModelParameters& model_parameters, int conc_l_scratch_id,
       int conc_a_scratch_id, int conc_b_scratch_id, int conc_l_ref_id,
       int conc_a_ref_id, int conc_b_ref_id, int partition_coeff_id,
       const int ncompositions, std::shared_ptr<tbox::Database> conc_db,
       std::shared_ptr<tbox::MemoryDatabase> newton_db)
   {
      std::shared_ptr<tbox::MemoryDatabase> calphad_db;
      if (model_parameters.isConcentrationModelCALPHAD()) {
         std::shared_ptr<tbox::Database> db(conc_db->getDatabase("Calphad"));
         std::string calphad_filename = db->getString("filename");
         calphad_db.reset(new tbox::MemoryDatabase("calphad_db"));
         tbox::InputManager::getManager()->parseInputFile(calphad_filename,
                                                          calphad_db);
      }

      std::shared_ptr<PhaseConcentrationsStrategy> phase_conc_strategy;

      if (model_parameters.kks_phase_concentration()) {
         tbox::plog << "Phase concentration determined by KKS" << std::endl;
         if (model_parameters.isConcentrationModelCALPHAD()) {
            if (ncompositions == 1)
               phase_conc_strategy.reset(
                   new CALPHADequilibriumPhaseConcentrationsStrategy<
                       CALPHADFreeEnergyFunctionsBinary>(
                       conc_l_scratch_id, conc_a_scratch_id, conc_b_scratch_id,
                       conc_l_ref_id, conc_a_ref_id, conc_b_ref_id,
                       model_parameters.energy_interp_func_type(),
                       model_parameters.conc_interp_func_type(),
                       model_parameters.with_third_phase(), calphad_db,
                       newton_db, ncompositions));
            else if (ncompositions == 2)
               phase_conc_strategy.reset(
                   new CALPHADequilibriumPhaseConcentrationsStrategy<
                       CALPHADFreeEnergyFunctionsTernary>(
                       conc_l_scratch_id, conc_a_scratch_id, conc_b_scratch_id,
                       conc_l_ref_id, conc_a_ref_id, conc_b_ref_id,
                       model_parameters.energy_interp_func_type(),
                       model_parameters.conc_interp_func_type(),
                       model_parameters.with_third_phase(), calphad_db,
                       newton_db, ncompositions));

         } else if (model_parameters.isConcentrationModelKKSdilute()) {
            phase_conc_strategy.reset(
                new KKSdiluteEquilibriumPhaseConcentrationsStrategy(
                    conc_l_scratch_id, conc_a_scratch_id, conc_b_scratch_id,
                    conc_l_ref_id, conc_a_ref_id, conc_b_ref_id,
                    model_parameters.energy_interp_func_type(),
                    model_parameters.conc_interp_func_type(), conc_db));
         } else {
            if (model_parameters.isConcentrationModelHBSM())
               phase_conc_strategy.reset(
                   new HBSMequilibriumPhaseConcentrationsStrategy(
                       conc_l_scratch_id, conc_a_scratch_id, conc_b_scratch_id,
                       model_parameters, conc_db));
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
