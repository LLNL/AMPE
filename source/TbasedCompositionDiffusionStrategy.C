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
#include "TbasedCompositionDiffusionStrategy.h"
#include "toolsSAMRAI.h"

#include "ConcFort.h"
#include "PhysicalConstants.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

TbasedCompositionDiffusionStrategy::TbasedCompositionDiffusionStrategy(
    const int pfm_diffusion_l_id, const int pfm_diffusion_a_id,
    const int diffusion_coeff_l_id, const int diffusion_coeff_a_id,
    const double D_liquid, const double Q0_liquid, const double D_solid_A,
    const double Q0_solid_A, const DiffusionInterpolationType interp_func_type,
    const std::string& avg_func_type)
    : CompositionDiffusionStrategy(interp_func_type),
      d_pfm_diffusion_l_id(pfm_diffusion_l_id),
      d_pfm_diffusion_a_id(pfm_diffusion_a_id),
      d_diffusion_coeff_l_id(diffusion_coeff_l_id),
      d_diffusion_coeff_a_id(diffusion_coeff_a_id),
      d_D_liquid(D_liquid),
      d_Q0_liquid(Q0_liquid),
      d_D_solid_A(D_solid_A),
      d_Q0_solid_A(Q0_solid_A),
      d_avg_func_type(avg_func_type)
{
   assert(D_liquid >= 0.);
   assert(Q0_liquid >= 0.);
   assert(Q0_solid_A >= 0.);
   assert(D_solid_A >= 0.);
}

void TbasedCompositionDiffusionStrategy::setDiffusion(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int phase_id, const int eta_id)
{
   (void)eta_id;

   // tbox::pout<<"TbasedCompositionDiffusionStrategy::setDiffusion()"<<endl;
   assert(temperature_id >= 0);
   assert(phase_id >= 0);
   assert(d_pfm_diffusion_l_id >= 0);
   assert(d_pfm_diffusion_a_id >= 0);

   const int maxl = hierarchy->getNumberOfLevels();

   const char interp_func_type = interpChar();

   for (int amr_level = 0; amr_level < maxl; amr_level++) {
      std::shared_ptr<hier::PatchLevel> level =
          hierarchy->getPatchLevel(amr_level);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         const hier::Box& pbox = patch->getBox();
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast = pbox.upper();

         std::shared_ptr<pdat::CellData<double> > phi(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(phase_id)));

         std::shared_ptr<pdat::CellData<double> > temperature(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(temperature_id)));

         std::shared_ptr<pdat::SideData<double> > pfm_diffusionL(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(d_pfm_diffusion_l_id)));
         assert(pfm_diffusionL->getGhostCellWidth()[0] == 0);

         std::shared_ptr<pdat::SideData<double> > pfm_diffusionA(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(d_pfm_diffusion_a_id)));

         assert(pfm_diffusionA->getDepth() == pfm_diffusionL->getDepth());
         assert(pfm_diffusionA->getDepth() == 1 ||  // binary
                pfm_diffusionA->getDepth() == 4);   // ternary

         // compute depth 0 of diffusion variables,
         // including phase fraction weight
         FORT_CONCENTRATION_PFMDIFFUSION_OF_T(
             ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
             ifirst(2), ilast(2),
#endif
             phi->getPointer(), phi->getGhostCellWidth()[0],
             pfm_diffusionL->getPointer(0, 0), pfm_diffusionL->getPointer(1, 0),
#if (NDIM == 3)
             pfm_diffusionL->getPointer(2, 0),
#endif
             pfm_diffusionA->getPointer(0, 0), pfm_diffusionA->getPointer(1, 0),
#if (NDIM == 3)
             pfm_diffusionA->getPointer(2, 0),
#endif
             0,  // assuming no ghosts for diffusion data
             temperature->getPointer(), temperature->getGhostCellWidth()[0],
             d_D_liquid, d_Q0_liquid, d_D_solid_A, d_Q0_solid_A,
             gas_constant_R_JpKpmol, &interp_func_type,
             d_avg_func_type.c_str());

         // fill other diagonal value with same value for ternaries for now
         if (pfm_diffusionL->getDepth() > 1) {
            pfm_diffusionL->copyDepth(3, *pfm_diffusionL, 0);
            pfm_diffusionA->copyDepth(3, *pfm_diffusionA, 0);
         }
      }
   }
}

void TbasedCompositionDiffusionStrategy::setDiffCoeffInEachPhase(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int eta_id)
{
   (void)eta_id;

   assert(temperature_id >= 0);
   assert(d_diffusion_coeff_l_id >= 0);
   assert(d_diffusion_coeff_a_id >= 0);

   const int maxl = hierarchy->getNumberOfLevels();

   const char interp_func_type = interpChar();

   for (int amr_level = 0; amr_level < maxl; amr_level++) {
      std::shared_ptr<hier::PatchLevel> level =
          hierarchy->getPatchLevel(amr_level);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         const hier::Box& pbox = patch->getBox();
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast = pbox.upper();

         std::shared_ptr<pdat::CellData<double> > temperature(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(temperature_id)));

         std::shared_ptr<pdat::SideData<double> > diffcoeffL(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(d_diffusion_coeff_l_id)));
         assert(diffcoeffL->getGhostCellWidth()[0] == 0);

         std::shared_ptr<pdat::SideData<double> > diffcoeffA(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(d_diffusion_coeff_a_id)));

         assert(diffcoeffA->getDepth() == diffcoeffL->getDepth());
         assert(diffcoeffA->getDepth() == 1 ||  // binary
                diffcoeffA->getDepth() == 4);   // ternary

         FORT_CONCENTRATION_DIFFCOEFF_OF_T(
             ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
             ifirst(2), ilast(2),
#endif
             diffcoeffL->getPointer(0, 0), diffcoeffL->getPointer(1, 0),
#if (NDIM == 3)
             diffcoeffL->getPointer(2, 0),
#endif
             diffcoeffA->getPointer(0, 0), diffcoeffA->getPointer(1, 0),
#if (NDIM == 3)
             diffcoeffA->getPointer(2, 0),
#endif
             0,  // assuming no ghosts for diffusion data
             temperature->getPointer(), temperature->getGhostCellWidth()[0],
             d_D_liquid, d_Q0_liquid, d_D_solid_A, d_Q0_solid_A,
             gas_constant_R_JpKpmol);
      }
   }
}
