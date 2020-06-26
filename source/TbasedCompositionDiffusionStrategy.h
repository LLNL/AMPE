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
#ifndef TbasedCompositionDiffusionStrategy_H
#define TbasedCompositionDiffusionStrategy_H

#include "CompositionDiffusionStrategy.h"

#include <cstring>

class TbasedCompositionDiffusionStrategy : public CompositionDiffusionStrategy
{
 public:
   TbasedCompositionDiffusionStrategy(
       const int pfm_diffusion_l_id, const int pfm_diffusion_a_id,
       const int diffusion_coeff_l_id, const int diffusion_coeff_a_id,
       const double D_liquid, const double Q0_liquid, const double D_solid_A,
       const double Q0_solid_A, DiffusionInterpolationType interp_func_type,
       const std::string& avg_func_type);


   /*
    * compute actual diffusion in each phase by weighting diffusion coefficients
    * in each phase with phase variable
    */
   virtual void setDiffusion(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int phase_id, const int eta_id);

   /*
    * Compute diffusion coefficient in each phase
    */
   virtual void setDiffCoeffInEachPhase(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int eta_scratch_id);

 private:
   /*!
    * holds data for diffusion coefficients in composition equation
    * weighted by phase fraction
    */
   int d_pfm_diffusion_l_id;
   int d_pfm_diffusion_a_id;

   int d_diffusion_coeff_l_id;
   int d_diffusion_coeff_a_id;

   double d_D_liquid;
   double d_Q0_liquid;

   double d_D_solid_A;
   double d_Q0_solid_A;

   std::string d_avg_func_type;
};

#endif
