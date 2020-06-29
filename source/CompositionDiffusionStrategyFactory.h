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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
#ifndef included_CompositionDiffusionStrategyFactory
#define included_CompositionDiffusionStrategyFactory

#include "DiffusionForConcInPhaseStrategy.h"
#include "TbasedCompositionDiffusionStrategy.h"

class CompositionDiffusionStrategyFactory
{
 public:
   static std::shared_ptr<CompositionDiffusionStrategy> create(
       QuatModel* model, QuatModelParameters& model_parameters,
       const short ncompositions, const int conc_l_scratch_id,
       const int conc_a_scratch_id, const int conc_b_scratch_id,
       const int conc_pfm_diffusion_l_id, const int conc_pfm_diffusion_a_id,
       const int conc_pfm_diffusion_b_id, const int conc_diffusion_coeff_l_id,
       const int conc_diffusion_coeff_a_id, const int conc_diffusion_coeff_b_id,
       CompositionStrategyMobilities* composition_strategy_mobilities,
       FreeEnergyStrategy* free_energy_strategy)
   {
      std::shared_ptr<CompositionDiffusionStrategy> strategy;

      if (model_parameters.conDiffusionStrategyIsCTD()) {
         strategy.reset(new DiffusionForConcInPhaseStrategy(
             ncompositions, conc_l_scratch_id, conc_a_scratch_id,
             conc_b_scratch_id, conc_pfm_diffusion_l_id,
             conc_pfm_diffusion_a_id, conc_pfm_diffusion_b_id,
             conc_diffusion_coeff_l_id, conc_diffusion_coeff_a_id,
             conc_diffusion_coeff_b_id, model_parameters.avg_func_type(),
             model_parameters.diffusion_interp_func_type(),
             composition_strategy_mobilities, free_energy_strategy));
      } else {
         // for T-dependent diffusion, phase fraction weight is
         // included in computation and d_conc_diffusion_coeff_*_id
         // are not set
         strategy.reset(new TbasedCompositionDiffusionStrategy(
             conc_pfm_diffusion_l_id, conc_pfm_diffusion_a_id,
             conc_diffusion_coeff_l_id, conc_diffusion_coeff_a_id,
             model_parameters.D_liquid(), model_parameters.Q0_liquid(),
             model_parameters.D_solid_A(), model_parameters.Q0_solid_A(),
             model_parameters.diffusion_interp_func_type(),
             model_parameters.avg_func_type()));
      }
      return strategy;
   }
};

#endif
