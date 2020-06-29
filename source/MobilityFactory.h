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
#ifndef included_MobilityFactory
#define included_MobilityFactory

#include "KimMobilityStrategyInfMob.h"
#include "KimMobilityStrategyFiniteMob.h"
#include "KimMobilityStrategyFiniteMobAntiTrap.h"
#include "QuatModelParameters.h"

class QuatModel;

class MobilityFactory
{
 public:
   static std::shared_ptr<QuatMobilityStrategy> create(
       QuatModel* model, QuatModelParameters& model_parameters,
       const int conc_l_scratch_id, const int conc_a_scratch_id,
       const int temperature_scratch_id, const unsigned ncompositions,
       std::shared_ptr<tbox::Database> conc_db)
   {
      std::shared_ptr<QuatMobilityStrategy> mobility_strategy;

      if (model_parameters.isPhaseMobilityScalar())
         mobility_strategy.reset(new SimpleQuatMobilityStrategy(model));
      else {
         // no support for composition dependent diffusion for now
         if (model_parameters.conDiffusionStrategyIsCTD()) {
            TBOX_ERROR(
                "No support for Kim's mobility for composition dependent "
                "diffusion!!\n");
         }

         if (model_parameters.interfaceMobility() > 0.) {
            if (model_parameters.with_antitrapping()) {
               mobility_strategy.reset(new KimMobilityStrategyFiniteMobAntiTrap(
                   model, conc_l_scratch_id, conc_a_scratch_id,
                   temperature_scratch_id, model_parameters.interfaceMobility(),
                   model_parameters.epsilon_phase(),
                   model_parameters.phase_well_scale(),
                   model_parameters.energy_interp_func_type(),
                   model_parameters.conc_interp_func_type(), conc_db,
                   ncompositions, model_parameters.D_liquid(),
                   model_parameters.Q0_liquid(),
                   model_parameters.molar_volume_liquid()));
            } else {
               mobility_strategy.reset(new KimMobilityStrategyFiniteMob(
                   model, conc_l_scratch_id, conc_a_scratch_id,
                   temperature_scratch_id, model_parameters.interfaceMobility(),
                   model_parameters.epsilon_phase(),
                   model_parameters.phase_well_scale(),
                   model_parameters.energy_interp_func_type(),
                   model_parameters.conc_interp_func_type(), conc_db,
                   ncompositions));
            }
         } else {
            mobility_strategy.reset(new KimMobilityStrategyInfMob(
                model, conc_l_scratch_id, conc_a_scratch_id,
                temperature_scratch_id, model_parameters.epsilon_phase(),
                model_parameters.phase_well_scale(),
                model_parameters.energy_interp_func_type(),
                model_parameters.conc_interp_func_type(), conc_db,
                ncompositions, model_parameters.D_liquid(),
                model_parameters.Q0_liquid(),
                model_parameters.molar_volume_liquid()));
         }
      }
      return mobility_strategy;
   }
};

#endif
