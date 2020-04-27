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
#ifndef included_SpinodalCompositionRHSStrategy
#define included_SpinodalCompositionRHSStrategy

#include "CompositionRHSStrategyWithMobilities.h"

#include <vector>
#include <string>

class SpinodalCompositionRHSStrategy
    : public CompositionRHSStrategyWithMobilities
{
 public:
   SpinodalCompositionRHSStrategy(
       std::shared_ptr<tbox::Database> input_db, const int conc_scratch_id,
       const int phase_scratch_id, const int eta_scratch_id,
       const unsigned int ncompositions, const int conc_a_scratch_id,
       const int conc_b_scratch_id, const int temperature_scratch_id,
       const int diffusion_id, const double kappa, const int Mq_id,
       const std::vector<double>& Q_heat_transport,
       const std::string& phase_interp_func_type,
       const std::string& avg_func_type,
       FreeEnergyStrategy* free_energy_strategy);

   ~SpinodalCompositionRHSStrategy(){};

   void computeFluxOnPatch(hier::Patch& patch, const int flux_id);

   void setDiffusionCoeff(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const double time);

 private:
   int d_conc_scratch_id;
   int d_eta_scratch_id;
   int d_diffusion_id;
   int d_temperature_scratch_id;
   int d_phase_scratch_id;
   int d_conc_a_scratch_id;
   int d_conc_b_scratch_id;

   unsigned int d_ncompositions;

   double d_kappa;

   void setDiffusionForConc(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy);

   void setDiffusionCoeffForConcOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_c,
       std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::SideData<double> > sd_d_coeff,
       const hier::Box& pbox);
};

#endif
