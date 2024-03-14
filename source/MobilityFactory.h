// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
#ifndef included_MobilityFactory
#define included_MobilityFactory

#include "KimMobilityStrategyInfMob.h"
#include "KimMobilityStrategyInfMob3Phases.h"
#include "KimMobilityStrategyFiniteMob.h"
#include "KimMobilityStrategyFiniteMobAntiTrap.h"
#include "QuatModelParameters.h"
#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "KKSFreeEnergyFunctionDiluteBinary.h"
#include "CALPHADFreeEnergyFunctionsBinary3Ph2Sl.h"

class QuatModel;

class MobilityFactory
{
 public:
   static std::shared_ptr<QuatMobilityStrategy> create(
       QuatModel* model, QuatModelParameters& model_parameters,
       const int conc_l_scratch_id, const int conc_a_scratch_id,
       const int conc_b_scratch_id, const int temperature_scratch_id,
       const unsigned ncompositions, const bool three_phases,
       std::shared_ptr<tbox::Database> conc_db)
   {
      if (model_parameters.isPhaseMobilityScalar()) {
         std::shared_ptr<QuatMobilityStrategy> mobility_strategy;
         mobility_strategy.reset(new SimpleQuatMobilityStrategy(model));
         return mobility_strategy;
      } else {
         // no support for composition dependent diffusion for now
         if (model_parameters.conDiffusionStrategyIsCTD()) {
            TBOX_ERROR(
                "No support for Kim's mobility for composition dependent "
                "diffusion!!\n");
         }

         std::string conc_model =
             conc_db->getStringWithDefault("model", "undefined");
         if (conc_model[0] == 'c') {
            if (ncompositions > 1)
               return createInternal<
                   Thermo4PFM::CALPHADFreeEnergyFunctionsTernary>(
                   model, model_parameters, conc_l_scratch_id,
                   conc_a_scratch_id, conc_b_scratch_id, temperature_scratch_id,
                   ncompositions, conc_db);
            else {
               if (!three_phases) {
                  tbox::plog << "FreeEnergyFunctionType: "
                                "CALPHADFreeEnergyFunctionsBinary"
                             << std::endl;
                  return createInternal<
                      Thermo4PFM::CALPHADFreeEnergyFunctionsBinary>(
                      model, model_parameters, conc_l_scratch_id,
                      conc_a_scratch_id, conc_b_scratch_id,
                      temperature_scratch_id, ncompositions, conc_db);
               } else {
                  return create3phasesMobility(model, model_parameters,
                                               conc_l_scratch_id,
                                               conc_a_scratch_id,
                                               conc_b_scratch_id,
                                               temperature_scratch_id,
                                               ncompositions, conc_db);
               }
            }
         } else {
            tbox::plog << "FreeEnergyFunctionType: "
                          "KKSFreeEnergyFunctionDiluteBinary"
                       << std::endl;
            return createInternal<
                Thermo4PFM::KKSFreeEnergyFunctionDiluteBinary>(
                model, model_parameters, conc_l_scratch_id, conc_a_scratch_id,
                conc_b_scratch_id, temperature_scratch_id, ncompositions,
                conc_db);
         }
      }
   }

 private:
   template <class FreeEnergyFunctionType>
   static std::shared_ptr<QuatMobilityStrategy> createInternal(
       QuatModel* model, QuatModelParameters& model_parameters,
       const int conc_l_scratch_id, const int conc_a_scratch_id,
       const int conc_b_scratch_id, const int temperature_scratch_id,
       const unsigned ncompositions, std::shared_ptr<tbox::Database> conc_db)
   {
      std::shared_ptr<QuatMobilityStrategy> mobility_strategy;

      if (model_parameters.interfaceMobility() > 0.) {
         if (model_parameters.with_antitrapping()) {
            tbox::plog << "KimMobilityStrategyFiniteMobAntiTrap" << std::endl;
            mobility_strategy.reset(new KimMobilityStrategyFiniteMobAntiTrap<
                                    FreeEnergyFunctionType>(
                model_parameters, model, conc_l_scratch_id, conc_a_scratch_id,
                temperature_scratch_id, model_parameters.interfaceMobility(),
                model_parameters.epsilon_phase(),
                model_parameters.phase_well_scale(),
                model_parameters.energy_interp_func_type(),
                model_parameters.conc_interp_func_type(), conc_db,
                ncompositions, model_parameters.D_liquid(),
                model_parameters.Q0_liquid(),
                model_parameters.molar_volume_liquid()));
         } else {
            tbox::plog << "KimMobilityStrategyFiniteMob" << std::endl;
            mobility_strategy.reset(
                new KimMobilityStrategyFiniteMob<FreeEnergyFunctionType>(
                    model, conc_l_scratch_id, conc_a_scratch_id,
                    temperature_scratch_id,
                    model_parameters.interfaceMobility(),
                    model_parameters.epsilon_phase(),
                    model_parameters.phase_well_scale(),
                    model_parameters.energy_interp_func_type(),
                    model_parameters.conc_interp_func_type(), conc_db,
                    ncompositions));
         }
      } else {
         tbox::plog << "KimMobilityStrategyInfMob" << std::endl;
         mobility_strategy.reset(
             new KimMobilityStrategyInfMob<FreeEnergyFunctionType>(
                 model_parameters, model, conc_l_scratch_id, conc_a_scratch_id,
                 temperature_scratch_id, model_parameters.epsilon_phase(),
                 model_parameters.phase_well_scale(),
                 model_parameters.energy_interp_func_type(),
                 model_parameters.conc_interp_func_type(), conc_db,
                 ncompositions, model_parameters.D_liquid(),
                 model_parameters.Q0_liquid(),
                 model_parameters.molar_volume_liquid()));
      }
      return mobility_strategy;
   }

   static std::shared_ptr<QuatMobilityStrategy> create3phasesMobility(
       QuatModel* model, QuatModelParameters& model_parameters,
       const int conc_l_scratch_id, const int conc_a_scratch_id,
       const int conc_b_scratch_id, const int temperature_scratch_id,
       const unsigned ncompositions, std::shared_ptr<tbox::Database> conc_db)
   {
      std::shared_ptr<QuatMobilityStrategy> mobility_strategy;

      tbox::plog << "KimMobilityStrategyInfMob3Phases" << std::endl;
      mobility_strategy.reset(
          new KimMobilityStrategyInfMob3Phases<
              Thermo4PFM::CALPHADFreeEnergyFunctionsBinary3Ph2Sl>(
              model_parameters, model, conc_l_scratch_id, conc_a_scratch_id,
              conc_b_scratch_id, temperature_scratch_id,
              model_parameters.epsilon_phase(),
              model_parameters.phase_well_scale(),
              model_parameters.energy_three_args_interp_func_type(),
              model_parameters.conc_interp_func_type(), conc_db, ncompositions,
              model_parameters.molar_volume_liquid()));
      return mobility_strategy;
   }
};

#endif
