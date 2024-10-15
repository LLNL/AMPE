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
#ifndef included_CompositionDiffusionStrategyFactory
#define included_CompositionDiffusionStrategyFactory

#include "MobilityCompositionDiffusionStrategy.h"
#include "TbasedCompositionDiffusionStrategy.h"
#include "ScalarCompositionDiffusionStrategy.h"
#include "WangSinteringCompositionDiffusionStrategy.h"

class CompositionDiffusionStrategyFactory
{
 public:
   static std::shared_ptr<CompositionDiffusionStrategy> create(
       QuatModel* model, QuatModelParameters& model_parameters,
       const short ncompositions, const int conc_id,
       const int conc_l_scratch_id, const int conc_a_scratch_id,
       const int conc_b_scratch_id,
       const std::vector<int>& conc_pfm_diffusion_id,
       const int conc_pfm_diffusion_l_id, const int conc_pfm_diffusion_a_id,
       const int conc_pfm_diffusion_b_id, const int conc_diffusion_coeff_l_id,
       const int conc_diffusion_coeff_a_id, const int conc_diffusion_coeff_b_id,
       CompositionStrategyMobilities* composition_strategy_mobilities,
       std::shared_ptr<FreeEnergyStrategy> free_energy_strategy)
   {
      std::shared_ptr<CompositionDiffusionStrategy> strategy;

      if (model_parameters.conDiffusionStrategyIsCTD()) {
         tbox::plog << "Uses composition dependent diffusion" << std::endl;
         strategy.reset(new MobilityCompositionDiffusionStrategy(
             ncompositions, conc_l_scratch_id, conc_a_scratch_id,
             conc_pfm_diffusion_l_id, conc_pfm_diffusion_a_id,
             conc_diffusion_coeff_l_id, conc_diffusion_coeff_a_id,
             model_parameters.avg_func_type(),
             model_parameters.diffusion_interp_func_type(),
             composition_strategy_mobilities, free_energy_strategy));
      } else if (model_parameters.concRHSstrategyIsEBS()) {
         // for T-dependent diffusion, phase fraction weight is
         // included in computation and d_conc_diffusion_coeff_*_id
         // are not set
         //
         tbox::plog << "Uses temperature based composition diffusion"
                    << std::endl;

         const bool three_phases_model = model_parameters.with_three_phases();
         const short norderpA =
             three_phases_model ? 1 : model_parameters.norderpA();
         const short norderpB =
             three_phases_model ? 1 : model_parameters.norderpB();

         strategy.reset(new TbasedCompositionDiffusionStrategy(
             model_parameters.norderp(), norderpA, norderpB, three_phases_model,
             conc_pfm_diffusion_l_id, conc_pfm_diffusion_a_id,
             conc_pfm_diffusion_b_id, model_parameters.D_liquid(),
             model_parameters.Q0_liquid(), model_parameters.D_solid_A(),
             model_parameters.Q0_solid_A(), model_parameters.D_solid_B(),
             model_parameters.Q0_solid_B(), model_parameters.D0_LA(),
             model_parameters.Q0_LA(), model_parameters.D0_LB(),
             model_parameters.Q0_LB(), model_parameters.D0_AA(),
             model_parameters.Q0_AA(), model_parameters.D0_AB(),
             model_parameters.Q0_AB(), model_parameters.D0_BB(),
             model_parameters.Q0_BB(),
             model_parameters.diffusion_interp_func_type(),
             model_parameters.avg_func_type()));
      } else if (model_parameters.isConcentrationModelWangSintering()) {
         strategy.reset(new WangSinteringCompositionDiffusionStrategy(
             conc_id, conc_pfm_diffusion_id[0], model_parameters.D_liquid(),
             model_parameters.D_solid_A(), model_parameters.D0_LA(),
             model_parameters.D0_AA(),
             model_parameters.diffusion_interp_func_type(),
             model_parameters.avg_func_type()));

      } else {
         tbox::plog << "Uses temperature based composition for scalar diffusion"
                    << std::endl;
         const bool three_phases_model = model_parameters.with_three_phases();
         const short norderpA =
             three_phases_model ? 1 : model_parameters.norderpA();
         const short norderpB =
             three_phases_model ? 1 : model_parameters.norderpB();

         strategy.reset(new ScalarCompositionDiffusionStrategy(
             model_parameters.norderp(), norderpA, norderpB, three_phases_model,
             conc_pfm_diffusion_id[0], model_parameters.D_liquid(),
             model_parameters.D_solid_A(), model_parameters.D_solid_B(),
             model_parameters.D0_LA(), model_parameters.D0_LB(),
             model_parameters.D0_AA(), model_parameters.D0_AB(),
             model_parameters.D0_BB(),
             model_parameters.diffusion_interp_func_type(),
             model_parameters.avg_func_type()));
      }

      return strategy;
   }
};

#endif
