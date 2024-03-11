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
#include "TbasedCompositionDiffusionStrategy.h"
#include "toolsSAMRAI.h"

#include "ConcFort.h"
#include "PhysicalConstants.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

const double gas_constant_R_JpKpmol = GASCONSTANT_R_JPKPMOL;

TbasedCompositionDiffusionStrategy::TbasedCompositionDiffusionStrategy(
    const short norderp, const short norderpA, const bool with3phases,
    const int pfm_diffusion_l_id, const int pfm_diffusion_a_id,
    const int pfm_diffusion_b_id, const double D_liquid, const double Q0_liquid,
    const double D_solid_A, const double Q0_solid_A, const double D_solid_B,
    const double Q0_solid_B, const DiffusionInterpolationType interp_func_type,
    const std::string& avg_func_type)
    : CompositionDiffusionStrategy(interp_func_type),
      d_norderp(norderp),
      d_norderpA(norderpA),
      d_with3phases(with3phases),
      d_pfm_diffusion_l_id(pfm_diffusion_l_id),
      d_pfm_diffusion_a_id(pfm_diffusion_a_id),
      d_pfm_diffusion_b_id(pfm_diffusion_b_id),
      d_D_liquid(D_liquid),
      d_Q0_liquid(Q0_liquid),
      d_D_solid_A(D_solid_A),
      d_Q0_solid_A(Q0_solid_A),
      d_D_solid_B(D_solid_B),
      d_Q0_solid_B(Q0_solid_B),
      d_avg_func_type(avg_func_type),
      d_with_phaseB(d_pfm_diffusion_b_id >= 0)
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

   //   tbox::plog<<"TbasedCompositionDiffusionStrategy::setDiffusion()"<<std::endl;
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

         std::shared_ptr<pdat::SideData<double> > pfm_diffusionA(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(d_pfm_diffusion_a_id)));

         assert(pfm_diffusionA->getDepth() == pfm_diffusionL->getDepth());
         assert(pfm_diffusionA->getDepth() == 1 ||  // binary
                pfm_diffusionA->getDepth() == 4);   // ternary

         // we need at least as many ghost cells for phi and T as in diffusion
         // to calculate diffusion in ghost cells
         assert(phi->getGhostCellWidth()[0] >=
                pfm_diffusionL->getGhostCellWidth()[0]);
         assert(temperature->getGhostCellWidth()[0] >=
                pfm_diffusionL->getGhostCellWidth()[0]);

         std::shared_ptr<pdat::SideData<double> > pfm_diffusionB;

         // compute depth 0 of diffusion variables,
         // including phase fraction weight
         if (d_with_phaseB) {
            pfm_diffusionB =
                SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                    patch->getPatchData(d_pfm_diffusion_b_id));

            assert(phi->getDepth() == d_norderp);

            // distinguish 3 phases and multiple phases conventions
            const double DL = d_with3phases ? d_D_solid_B : d_D_liquid;
            const double DA = d_with3phases ? d_D_liquid : d_D_solid_A;
            const double DB = d_with3phases ? d_D_solid_A : d_D_solid_B;

            const double QL = d_with3phases ? d_Q0_solid_B : d_Q0_liquid;
            const double QA = d_with3phases ? d_Q0_liquid : d_Q0_solid_A;
            const double QB = d_with3phases ? d_Q0_solid_A : d_Q0_solid_B;

            const double* diffLptr0 = d_with3phases
                                          ? pfm_diffusionB->getPointer(0, 0)
                                          : pfm_diffusionL->getPointer(0, 0);
            const double* diffLptr1 = d_with3phases
                                          ? pfm_diffusionB->getPointer(1, 0)
                                          : pfm_diffusionL->getPointer(1, 0);
#if (NDIM == 3)
            const double* diffLptr2 = d_with3phases
                                          ? pfm_diffusionB->getPointer(2, 0)
                                          : pfm_diffusionL->getPointer(2, 0);
#endif

            const double* diffAptr0 = d_with3phases
                                          ? pfm_diffusionL->getPointer(0, 0)
                                          : pfm_diffusionA->getPointer(0, 0);
            const double* diffAptr1 = d_with3phases
                                          ? pfm_diffusionL->getPointer(1, 0)
                                          : pfm_diffusionA->getPointer(1, 0);
#if (NDIM == 3)
            const double* diffAptr2 = d_with3phases
                                          ? pfm_diffusionL->getPointer(2, 0)
                                          : pfm_diffusionA->getPointer(2, 0);
#endif

            const double* diffBptr0 = d_with3phases
                                          ? pfm_diffusionA->getPointer(0, 0)
                                          : pfm_diffusionB->getPointer(0, 0);
            const double* diffBptr1 = d_with3phases
                                          ? pfm_diffusionA->getPointer(1, 0)
                                          : pfm_diffusionB->getPointer(1, 0);
#if (NDIM == 3)
            const double* diffBptr2 = d_with3phases
                                          ? pfm_diffusionA->getPointer(2, 0)
                                          : pfm_diffusionB->getPointer(2, 0);
#endif

            const short norderpA = d_with3phases ? 1 : d_norderpA;

            CONCENTRATION_PFMDIFFUSION_OF_TEMPERATURE_THREEPHASES(
                ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                ifirst(2), ilast(2),
#endif
                phi->getPointer(), d_norderp, norderpA,
                phi->getGhostCellWidth()[0], diffLptr0, diffLptr1,
#if (NDIM == 3)
                diffLptr2,
#endif
                diffAptr0, diffAptr1,
#if (NDIM == 3)
                diffAptr2,
#endif
                diffBptr0, diffBptr1,
#if (NDIM == 3)
                diffBptr2,
#endif
                pfm_diffusionL->getGhostCellWidth()[0],
                temperature->getPointer(), temperature->getGhostCellWidth()[0],
                DL, QL, DA, QA, DB, QB, gas_constant_R_JpKpmol,
                &interp_func_type, d_avg_func_type.c_str());
         } else if (phi->getDepth() > 1) {
            CONCENTRATION_PFMDIFFUSION_OF_TEMPERATURE_MULTIPHASES(
                ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                ifirst(2), ilast(2),
#endif
                phi->getPointer(), phi->getDepth(), phi->getGhostCellWidth()[0],
                pfm_diffusionL->getPointer(0, 0),
                pfm_diffusionL->getPointer(1, 0),
#if (NDIM == 3)
                pfm_diffusionL->getPointer(2, 0),
#endif
                pfm_diffusionA->getPointer(0, 0),
                pfm_diffusionA->getPointer(1, 0),
#if (NDIM == 3)
                pfm_diffusionA->getPointer(2, 0),
#endif
                pfm_diffusionL->getGhostCellWidth()[0],
                temperature->getPointer(), temperature->getGhostCellWidth()[0],
                d_D_liquid, d_Q0_liquid, d_D_solid_A, d_Q0_solid_A,
                gas_constant_R_JpKpmol, &interp_func_type,
                d_avg_func_type.c_str());
         } else {
            CONCENTRATION_PFMDIFFUSION_OF_TEMPERATURE(
                ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                ifirst(2), ilast(2),
#endif
                phi->getPointer(), phi->getGhostCellWidth()[0],
                pfm_diffusionL->getPointer(0, 0),
                pfm_diffusionL->getPointer(1, 0),
#if (NDIM == 3)
                pfm_diffusionL->getPointer(2, 0),
#endif
                pfm_diffusionA->getPointer(0, 0),
                pfm_diffusionA->getPointer(1, 0),
#if (NDIM == 3)
                pfm_diffusionA->getPointer(2, 0),
#endif
                pfm_diffusionL->getGhostCellWidth()[0],
                temperature->getPointer(), temperature->getGhostCellWidth()[0],
                d_D_liquid, d_Q0_liquid, d_D_solid_A, d_Q0_solid_A,
                gas_constant_R_JpKpmol, &interp_func_type,
                d_avg_func_type.c_str());
         }

         // fill other diagonal value with same value for ternaries for now
         if (pfm_diffusionL->getDepth() > 1) {
            pfm_diffusionL->copyDepth(3, *pfm_diffusionL, 0);
            pfm_diffusionA->copyDepth(3, *pfm_diffusionA, 0);
            if (d_with_phaseB) pfm_diffusionB->copyDepth(3, *pfm_diffusionB, 0);
         }
      }
   }
}
