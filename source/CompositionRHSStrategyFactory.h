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

class CompositionRHSStrategyFactory
{
 public:
   static std::shared_ptr<CompositionRHSStrategy> create(
       QuatModel* quatmodel,
       FreeEnergyStrategy* free_energy_strategy_for_diffusion,
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
         strategy.reset(new KKSCompositionRHSStrategy(
             conc_scratch_id, phase_scratch_id,
             conc_pfm_diffusion_id[0],  // use 1x1 diffusion matrix
             conc_phase_coupling_diffusion_id, temperature_scratch_id,
             eta_scratch_id, conc_eta_coupling_diffusion_id, conc_l_scratch_id,
             conc_a_scratch_id, conc_b_scratch_id, model_parameters.D_liquid(),
             model_parameters.D_solid_A(), model_parameters.D_solid_B(),
             model_parameters.Q0_liquid(), model_parameters.Q0_solid_A(),
             model_parameters.Q0_solid_B(),
             model_parameters.energy_interp_func_type(),
             model_parameters.avg_func_type()));
      } else if (model_parameters.concRHSstrategyIsEBS()) {

         assert(d_diffusion_for_conc_in_phase);

         strategy.reset(new EBSCompositionRHSStrategy(
             phase_scratch_id, eta_scratch_id,
             static_cast<unsigned short>(ncompositions), conc_l_scratch_id,
             conc_a_scratch_id, conc_b_scratch_id, temperature_scratch_id,
             conc_pfm_diffusion_l_id, conc_pfm_diffusion_a_id,
             conc_pfm_diffusion_b_id, conc_Mq_id,
             model_parameters.Q_heat_transport(), conc_pfm_diffusion_id,
             model_parameters.avg_func_type(),
             free_energy_strategy_for_diffusion,
             composition_strategy_mobilities, diffusion_for_conc_in_phase));
      } else if (model_parameters.concRHSstrategyIsBeckermann()) {
         strategy.reset(new BeckermannCompositionRHSStrategy(
             quatmodel, conc_scratch_id, phase_scratch_id,
             partition_coeff_scratch_id, conc_pfm_diffusion_id[0],
             conc_phase_coupling_diffusion_id, model_parameters.D_liquid(),
             model_parameters.D_solid_A(),
             model_parameters.conc_interp_func_type(),
             model_parameters.avg_func_type()));
      } else {
         TBOX_ERROR("Error: unknown composition RHS Strategy");
      }

      return strategy;
   }
};

#endif
