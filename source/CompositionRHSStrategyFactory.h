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
#ifndef included_CompositionRHSstrategyFactory
#define included_CompositionRHSstrategyFactory

#include "KKSCompositionRHSStrategy.h"
#include "EBSCompositionRHSStrategy.h"
#include "BeckermannCompositionRHSStrategy.h"
#include "QuatModelParameters.h"
#include "QuatModel.h"
#include "CompositionStrategyMobilities.h"
#include "CompositionDiffusionStrategy.h"
#include "CahnHilliardDoubleWell.h"
#include "WangSinteringCompositionRHSStrategy.h"

class CompositionRHSStrategyFactory
{
 public:
   static std::shared_ptr<CompositionRHSStrategy> create(
       QuatModel* quatmodel,
       std::shared_ptr<FreeEnergyStrategy> free_energy_strategy_for_diffusion,
       QuatModelParameters& model_parameters, int ncompositions,
       int conc_scratch_id, int phase_scratch_id, int temperature_scratch_id,
       std::vector<int>& conc_pfm_diffusion_id, int conc_pfm_diffusion_l_id,
       int conc_pfm_diffusion_a_id, int conc_pfm_diffusion_b_id,
       int conc_phase_coupling_diffusion_id, int eta_scratch_id,
       int conc_eta_coupling_diffusion_id, int conc_l_scratch_id,
       int conc_a_scratch_id, int conc_b_scratch_id,
       int partition_coeff_scratch_id, int conc_Mq_id,
       CompositionStrategyMobilities* composition_strategy_mobilities,
       std::shared_ptr<CompositionDiffusionStrategy>
           diffusion_for_conc_in_phase)
   {
      std::shared_ptr<CompositionRHSStrategy> strategy;

      if (model_parameters.concRHSstrategyIsKKS()) {
         tbox::plog << "Use KKSCompositionRHSStrategy..." << std::endl;
         strategy.reset(new KKSCompositionRHSStrategy(
             conc_scratch_id, phase_scratch_id,
             conc_pfm_diffusion_id[0],  // use 1x1 diffusion matrix
             conc_phase_coupling_diffusion_id, temperature_scratch_id,
             conc_l_scratch_id, conc_a_scratch_id, model_parameters.D_liquid(),
             model_parameters.D_solid_A(), model_parameters.Q0_liquid(),
             model_parameters.Q0_solid_A(),
             model_parameters.energy_interp_func_type(),
             model_parameters.avg_func_type()));
      } else if (model_parameters.concRHSstrategyIsEBS()) {
         tbox::plog << "Use EBSCompositionRHSStrategy..." << std::endl;

         assert(diffusion_for_conc_in_phase);

         std::string stencil_type = model_parameters.useEBSisotropicStencil()
                                        ? "isotropic"
                                        : "regular";
         stencil_type =
             model_parameters.useEBS4thOrderStencil() ? "order4" : stencil_type;

         strategy.reset(new EBSCompositionRHSStrategy(
             phase_scratch_id, static_cast<unsigned short>(ncompositions),
             conc_l_scratch_id, conc_a_scratch_id, conc_b_scratch_id,
             temperature_scratch_id, conc_pfm_diffusion_l_id,
             conc_pfm_diffusion_a_id, conc_pfm_diffusion_b_id,
             conc_pfm_diffusion_id, model_parameters.avg_func_type(),
             free_energy_strategy_for_diffusion, diffusion_for_conc_in_phase,
             stencil_type));
      } else if (model_parameters.concRHSstrategyIsBeckermann()) {
         strategy.reset(new BeckermannCompositionRHSStrategy(
             quatmodel, conc_scratch_id, phase_scratch_id,
             partition_coeff_scratch_id, conc_pfm_diffusion_id[0],
             conc_phase_coupling_diffusion_id, model_parameters.D_liquid(),
             model_parameters.D_solid_A(),
             model_parameters.conc_interp_func_type(),
             model_parameters.avg_func_type()));
      } else if (model_parameters.concRHSstrategyIsCahnHilliard()) {
         // use mobility = 1. since multiplied by mobility outside
         strategy.reset(new CahnHilliardDoubleWell(
             conc_scratch_id, 1., model_parameters.CH_ca(),
             model_parameters.CH_cb(), model_parameters.CH_kappa(),
             model_parameters.CH_well_scale(),
             model_parameters.avg_func_type()));
      } else if (model_parameters.isConcentrationModelWangSintering()) {
         tbox::plog << "Use WangSinteringCompositionRHSStrategy..."
                    << std::endl;
         strategy.reset(new WangSinteringCompositionRHSStrategy(
             conc_scratch_id, phase_scratch_id, temperature_scratch_id,
             conc_pfm_diffusion_id[0], model_parameters.conc_mobility(),
             model_parameters.WangSintering_A(),
             model_parameters.WangSintering_B(),
             model_parameters.WangSintering_beta_rho(),
             model_parameters.avg_func_type(), diffusion_for_conc_in_phase));
      } else {
         TBOX_ERROR("Error: unknown composition RHS Strategy");
      }

      return strategy;
   }
};

#endif
