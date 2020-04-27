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
#ifndef included_KKSCompositionRHSStrategy
#define included_KKSCompositionRHSStrategy

#include "CompositionRHSStrategy.h"
#include "InterpolationType.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"

#include <string>

class KKSCompositionRHSStrategy : public CompositionRHSStrategy
{
 public:
   KKSCompositionRHSStrategy(
       const int conc_scratch_id, const int phase_scratch_id,
       const int diffusion0_id, const int phase_coupling_diffusion_id,
       const int temperature_scratch_id, const int eta_scratch_id,
       const int conc_eta_coupling_diffusion_id, const int conc_l_scratch_id,
       const int conc_a_scratch_id, const int conc_b_scratch_id,
       const double D_liquid, const double D_solid_A, const double D_solid_B,
       const double Q0_liquid, const double Q0_solid_A, const double Q0_solid_B,
       const EnergyInterpolationType phase_interp_func_type,
       const std::string& avg_func_type);
   ~KKSCompositionRHSStrategy(){};

   void computeFluxOnPatch(hier::Patch& patch, const int flux_id);
   void setDiffusionCoeff(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const double time);

 private:
   int d_conc_scratch_id;
   int d_phase_scratch_id;
   int d_eta_scratch_id;

   int d_conc_l_scratch_id;
   int d_conc_a_scratch_id;
   int d_conc_b_scratch_id;

   int d_temperature_scratch_id;

   int d_pfm_diffusion_id;
   int d_conc_phase_coupling_diffusion_id;
   int d_conc_eta_coupling_diffusion_id;

   bool d_with_third_phase;

   EnergyInterpolationType d_phase_interp_func_type;
   std::string d_avg_func_type;

   double d_D_liquid;
   double d_D_solid_A;
   double d_D_solid_B;
   double d_Q0_liquid;
   double d_Q0_solid_A;
   double d_Q0_solid_B;

   // Timers
   std::shared_ptr<tbox::Timer> t_set_diffcoeff_timer;

   void setPFMDiffCoeffForConcentration(
       const std::shared_ptr<hier::PatchHierarchy>, const int temperature_id,
       const int phase_id, const int eta_id, const int conc_pfm_diffusion_id);
   void setDiffCoeffForGradPhi(const std::shared_ptr<hier::PatchHierarchy>,
                               const int temperature_id,
                               const int concentration_id, const int phase_id,
                               const int eta_id,
                               const int conc_pfm_diffusion_id,
                               const int conc_phase_coupling_diffusion_id,
                               const int conc_eta_coupling_diffusion_id);

   void setDiffCoeffForPhaseOnPatch(
       std::shared_ptr<pdat::SideData<double> > sd_phi_diff_coeff,
       std::shared_ptr<pdat::SideData<double> > sd_eta_diff_coeff,
       std::shared_ptr<pdat::SideData<double> > sd_pfmd_coeff,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_eta,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b,
       const hier::Box& pbox);
   void computeDiffusionOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_eta,
       std::shared_ptr<pdat::CellData<double> > cd_concentration,
       std::shared_ptr<pdat::SideData<double> > sd_diffusion0,
       std::shared_ptr<pdat::SideData<double> > sd_diffusion,
       const hier::Box& pbox);
};

#endif
