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
#include "WangSinteringDiffusion.h"
#include "toolsSAMRAI.h"

#include "ConcFort.h"

WangSinteringDiffusion::WangSinteringDiffusion(
    const int conc_id, const int diffusion_id, const double D0_liquid,
    const double D0_solidA, const double D0_LA, const double D0_AA,
    const DiffusionInterpolationType interp_func_type,
    const std::string& avg_func_type)
    : CompositionDiffusionStrategy(interp_func_type),
      d_conc_id(conc_id),
      d_diffusion_id(diffusion_id),
      d_D0_liquid(D0_liquid),
      d_D0_solidA(D0_solidA),
      d_d0_LA(D0_LA),
      d_d0_AA(D0_AA),
      d_avg_func_type(avg_func_type)
{
   assert(d_conc_id >= 0);

   assert(D0_liquid >= 0.);
   assert(D0_solidA >= 0.);

   tbox::plog << "WangSinteringDiffusion..." << std::endl;
}

void WangSinteringDiffusion::setDiffusion(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int temperature_id,
    int phase_id)
{
   (void)temperature_id;
   tbox::plog << "WangSinteringDiffusion::setDiffusion()" << std::endl;
   assert(phase_id >= 0);
   assert(d_diffusion_id >= 0);
   assert(d_conc_id >= 0);

   const int maxl = hierarchy->getNumberOfLevels();

   for (int amr_level = 0; amr_level < maxl; amr_level++) {
      std::shared_ptr<hier::PatchLevel> level =
          hierarchy->getPatchLevel(amr_level);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         std::shared_ptr<pdat::CellData<double>> phi(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(phase_id)));

         std::shared_ptr<pdat::CellData<double>> conc(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_conc_id)));
         assert(conc);

         std::shared_ptr<pdat::SideData<double>> diffusion(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(d_diffusion_id)));

         assert(diffusion->getDepth() == 1 ||  // binary
                diffusion->getDepth() == 4);   // ternary

         setDiffusion(patch, phi, conc, diffusion);
      }
   }
}

void WangSinteringDiffusion::setDiffusionInterfaces(
    std::shared_ptr<hier::Patch> patch,
    std::shared_ptr<pdat::CellData<double>> phi,
    std::shared_ptr<pdat::CellData<double>> conc,
    std::shared_ptr<pdat::SideData<double>> diffusion)
{
   assert(diffusion->getGhostCellWidth()[0] > 0);

   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   /*
    *  grain boundary contribution
    */
   const short norderp = phi->getDepth();
   // to avoid adding interface diffusion twice to an
   // interface between two grains of same phase, flag
   // case with dupl=1
   int same_orderp = 1;
   ADD_AB_DIFFUSION(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                    ifirst(2), ilast(2),
#endif
                    phi->getPointer(0), norderp, phi->getPointer(0), norderp,
                    phi->getGhostCellWidth()[0], diffusion->getPointer(0, 0),
                    diffusion->getPointer(1, 0),
#if (NDIM == 3)
                    diffusion->getPointer(2, 0),
#endif
                    diffusion->getGhostCellWidth()[0], d_d0_AA, same_orderp);

   /*
    * surface contribution
    */
   ADD_AB_DIFFUSION_SINGLE(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                           ifirst(2), ilast(2),
#endif
                           conc->getPointer(), conc->getGhostCellWidth()[0],
                           diffusion->getPointer(0, 0),
                           diffusion->getPointer(1, 0),
#if (NDIM == 3)
                           diffusion->getPointer(2, 0),
#endif
                           diffusion->getGhostCellWidth()[0], d_d0_LA);
}

void WangSinteringDiffusion::setDiffusion(
    std::shared_ptr<hier::Patch> patch,
    std::shared_ptr<pdat::CellData<double>> phi,
    std::shared_ptr<pdat::CellData<double>> conc,
    std::shared_ptr<pdat::SideData<double>> diffusion)
{
   assert(phi);
   assert(conc);
   assert(diffusion);
   assert(conc->getDepth() == 1);

   const char interp_func_type = interpChar();

   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   // we need at least as many ghost cells for phi and T as in diffusion
   // to calculate diffusion in ghost cells
   assert(phi->getGhostCellWidth()[0] >= diffusion->getGhostCellWidth()[0]);

   PFMDIFFUSION_SCALAR(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                       ifirst(2), ilast(2),
#endif
                       conc->getPointer(), conc->getGhostCellWidth()[0],
                       diffusion->getPointer(0, 0), diffusion->getPointer(1, 0),
#if (NDIM == 3)
                       diffusion->getPointer(2, 0),
#endif
                       diffusion->getGhostCellWidth()[0], d_D0_liquid,
                       d_D0_solidA, &interp_func_type, d_avg_func_type.c_str());

   setDiffusionInterfaces(patch, phi, conc, diffusion);

   // fill other diagonal value with same value for ternaries for now
   if (diffusion->getDepth() > 1) {
      diffusion->copyDepth(3, *diffusion, 0);
   }
}
