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
// Ref: Eiken, Boettger, Steinbach, PRE 73, 066122 (2006)
#ifndef included_EBSCompositionRHSStrategy
#define included_EBSCompositionRHSStrategy

#include "CompositionRHSStrategy.h"
#include "CALPHADMobility.h"
#include "CALPHADFreeEnergyStrategyBinary.h"
#include "CompositionDiffusionStrategy.h"

#include <vector>
#include <string>

class CompositionDiffusionStrategy;
class CompositionStrategyMobilities;

class EBSCompositionRHSStrategy : public CompositionRHSStrategy
{
 public:
   EBSCompositionRHSStrategy(
       const int phase_scratch_id, const int eta_scratch_id,
       const unsigned short ncompositions, const int conc_l_scratch_id,
       const int conc_a_scratch_id, const int conc_b_scratch_id,
       const int temperature_scratch_id, const int diffusion_l_id,
       const int diffusion_a_id, const int diffusion_b_id, const int Mq_id,
       const std::vector<double>& Q_heat_transport,
       const std::vector<int> diffusion_precond_id,
       const std::string& avg_func_type,
       FreeEnergyStrategy* free_energy_strategy,
       CompositionStrategyMobilities* mobilities_strategy,
       std::shared_ptr<CompositionDiffusionStrategy>
           diffusion_for_conc_in_phase);

   ~EBSCompositionRHSStrategy(){};

   void addFluxFromAntitrappingonPatch(hier::Patch& patch,
                                       const int phase_scratch_id,
                                       const int dphidt_id, const double alpha,
                                       const int flux_id);

   void computeFluxOnPatch(hier::Patch& patch, const int flux_id);

   void setDiffusionCoeff(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const double time);

   // void printDiagnostics(const std::shared_ptr<hier::PatchHierarchy >
   // hierarchy);

   /*
    * Take sum of diffusion coefficients in each phase
    * to get diffusion used in preconditioner
    * (same diffusion coefficient as used in KKS equations)
    * Assumes coefficients in each pahse includes a phase fraction weight
    */
   void setDiffusionCoeffForPreconditioner(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy);

 private:
   unsigned short d_ncompositions;

   std::vector<double> d_Q_heat_transport;

   int d_phase_scratch_id;
   int d_eta_scratch_id;

   int d_conc_l_scratch_id;
   int d_conc_a_scratch_id;
   int d_conc_b_scratch_id;

   /*!
    * holds data for diffusion coefficients in composition equation
    * weighted by phase fraction
    */
   int d_diffusion_l_id;
   int d_diffusion_a_id;
   int d_diffusion_b_id;

   int d_Mq_id;

   int d_temperature_scratch_id;

   bool d_with_third_phase;

   bool d_with_gradT;

   bool d_with_diffusion_for_preconditioner;
   std::vector<int> d_diffusion_precond_id;

   // free energy needed to compute diffusion in each phase
   FreeEnergyStrategy* d_free_energy_strategy;

   std::shared_ptr<CompositionDiffusionStrategy>
       d_diffusion_for_conc_in_phase;

   CompositionStrategyMobilities* d_mobilities_strategy;

   void addFluxFromGradTonPatch(hier::Patch& patch, const int temperature_id,
                                const int flux_id);

   void setDiffusionCoeffForPreconditionerOnPatch(
       std::shared_ptr<pdat::SideData<double> > sd_d_l,
       std::shared_ptr<pdat::SideData<double> > sd_d_a,
       std::shared_ptr<pdat::SideData<double> > sd_d_b,
       const int depth_in_Dmatrix,
       std::shared_ptr<pdat::SideData<double> > sd_d_coeff,
       const hier::Box& pbox, const int depth);
   void setDiffusionCoeffForTOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_c,
       std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::SideData<double> > sd_d_coeff,
       const hier::Box& pbox);
   void setDiffusionCoeffForTOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b,
       std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_eta,
       std::shared_ptr<pdat::SideData<double> > mq, const hier::Box& pbox);
   void setDiffusionCoeffForT(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int concentration_l_id, const int concentration_a_id,
       const int concentration_b_id, const int temperature_id);
};

#endif
