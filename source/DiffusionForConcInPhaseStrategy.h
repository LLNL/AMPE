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
#ifndef included_DIFFUSIONFORCONCINPHASESTRATEGY
#define included_DIFFUSIONFORCONCINPHASESTRATEGY

#include "FuncFort.h"
#include "CompositionDiffusionStrategy.h"

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

#include <string>
using namespace SAMRAI;

class CompositionStrategyMobilities;
class FreeEnergyStrategy;

class DiffusionForConcInPhaseStrategy : public CompositionDiffusionStrategy
{
 public:
   DiffusionForConcInPhaseStrategy(
       const unsigned short ncompositions, const int conc_l_scratch_id,
       const int conc_a_scratch_id, const int conc_b_scratch_id,
       const int pfm_diffusion_l_id, const int pfm_diffusion_a_id,
       const int pfm_diffusion_b_id, const int diffusion_coeff_l_id,
       const int diffusion_coeff_a_id, const int diffusion_coeff_b_id,
       const std::string& avg_func_type,
       DiffusionInterpolationType diff_interp_type,
       CompositionStrategyMobilities* mobilities_strategy,
       FreeEnergyStrategy* free_energy_strategy);

   ~DiffusionForConcInPhaseStrategy(){};

   /*
    * compute actual diffusion in each phase by weighting diffusion coefficients
    * in each phase with phase variable
    */
   void setDiffusion(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                     const int temperature_id, const int phase_id,
                     const int eta_id);

   /*
    * Compute diffusion coefficient in each phase
    */
   void setDiffCoeffInEachPhase(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int eta_scratch_id);

 private:
   double average(const double a, const double b) const
   {
      return FORT_AVERAGE_FUNC(a, b, d_avg_func_type.c_str());
   }

   /*
    * compute PFM diffusion in each phase based on diffusion coefficients
    * in each phase
    */
   void setPFMDiffOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_eta,
       std::shared_ptr<pdat::SideData<double> > sd_d_coeff_l,
       std::shared_ptr<pdat::SideData<double> > sd_d_coeff_a,
       std::shared_ptr<pdat::SideData<double> > sd_d_coeff_b,
       std::shared_ptr<pdat::SideData<double> > sd_d_l,
       std::shared_ptr<pdat::SideData<double> > sd_d_a,
       std::shared_ptr<pdat::SideData<double> > sd_d_b,
       const hier::Box& pbox);

   void setDiffCoeffInEachPhaseOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b,
       std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::CellData<double> > cd_eta,
       std::shared_ptr<pdat::SideData<double> > sd_d_coeff_l,
       std::shared_ptr<pdat::SideData<double> > sd_d_coeff_a,
       std::shared_ptr<pdat::SideData<double> > sd_d_coeff_b,
       const hier::Box& pbox);

   void computeLocalDiffusionMatrixL(const double temperature,
                                     const std::vector<double>& c);
   void computeLocalDiffusionMatrixA(const double temperature,
                                     const std::vector<double>& c);
   void computeLocalDiffusionMatrixB(const double temperature,
                                     const std::vector<double>& c);

   bool d_same_composition_for_third_phase;

   unsigned short d_ncompositions;

   int d_conc_l_scratch_id;
   int d_conc_a_scratch_id;
   int d_conc_b_scratch_id;

   /*!
    * holds data for diffusion coefficients in composition equation
    * weighted by phase fraction
    */
   int d_pfm_diffusion_l_id;
   int d_pfm_diffusion_a_id;
   int d_pfm_diffusion_b_id;

   /*!
    * holds data for diffusion coefficients in each phase
    */
   int d_diffusion_coeff_l_id;
   int d_diffusion_coeff_a_id;
   int d_diffusion_coeff_b_id;

   bool d_with_third_phase;

   CompositionStrategyMobilities* d_mobilities_strategy;

   // free energy needed to compute diffusion in each phase
   FreeEnergyStrategy* d_free_energy_strategy;

   /*!
    * function to use to take averages between two phase values
    * at cell centers and define quantity at cell boundary
    */
   std::string d_avg_func_type;

   /*!
    * small work arrays
    */
   std::vector<double> d_d2f;
   std::vector<double> d_mobmat;
   std::vector<double> d_local_dmat;
};

#endif
