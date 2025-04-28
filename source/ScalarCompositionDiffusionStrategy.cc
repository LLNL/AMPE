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
#include "ScalarCompositionDiffusionStrategy.h"
#include "toolsSAMRAI.h"

#include "ConcFort.h"
#include "PhysicalConstants.h"

const double gas_constant_R_JpKpmol = GASCONSTANT_R_JPKPMOL;

ScalarCompositionDiffusionStrategy::ScalarCompositionDiffusionStrategy(
    const short norderp, const short norderpA, const short norderpB,
    const bool with3phases, const int pfm_diffusion_id, const double D0_liquid,
    const double D0_solidA, const double D0_solidB, const double D0_LA,
    const double D0_LB, const double D0_AA, const double D0_AB,
    const double D0_BB, const DiffusionInterpolationType interp_func_type,
    const std::string& avg_func_type)
    : CompositionDiffusionStrategy(interp_func_type),
      d_norderp(norderp),
      d_norderpA(norderpA),
      d_norderpB(norderpB),
      d_with3phases(with3phases),
      d_pfm_diffusion_id(pfm_diffusion_id),
      d_D0_liquid(D0_liquid),
      d_D0_solidA(D0_solidA),
      d_D0_solidB(D0_solidB),
      d_d0_LA(D0_LA),
      d_d0_LB(D0_LB),
      d_d0_AA(D0_AA),
      d_d0_AB(D0_AB),
      d_d0_BB(D0_BB),
      d_avg_func_type(avg_func_type),
      d_with_phaseB(D0_solidB > 0.)
{
   assert(D0_liquid >= 0.);
   assert(D0_solidA >= 0.);

   tbox::plog << "ScalarCompositionDiffusionStrategy: D0_AB = " << D0_AB
              << std::endl;
}

void ScalarCompositionDiffusionStrategy::setDiffusion(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int temperature_id,
    int phase_id)
{
   (void)temperature_id;
   // tbox::plog<<"ScalarCompositionDiffusionStrategy::setDiffusion()"<<std::endl;
   assert(phase_id >= 0);
   assert(d_pfm_diffusion_id >= 0);

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

         std::shared_ptr<pdat::SideData<double>> pfm_diffusion(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(d_pfm_diffusion_id)));

         assert(pfm_diffusion->getDepth() == 1 ||  // binary
                pfm_diffusion->getDepth() == 4);   // ternary

         setDiffusion(patch, phi, pfm_diffusion);
      }
   }
}

void ScalarCompositionDiffusionStrategy::setDiffusionInterfaces(
    std::shared_ptr<hier::Patch> patch,
    std::shared_ptr<pdat::CellData<double>> phi,
    std::shared_ptr<pdat::SideData<double>> pfm_diffusion)
{
   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   const short nphases = d_with_phaseB ? 3 : 2;
   // std::cout << "nphases = " << nphases << std::endl;
   // if (d_with3phases) std::cout << "with 3 phases..." << std::endl;
   // std::cout << "d_norderpA = " << d_norderpA << std::endl;

   // distinguish 3 phases and multiple phases conventions
   const double* const phiL =
       d_with3phases ? phi->getPointer(0) : phi->getPointer(d_norderp - 1);

   const double* const phiA =
       d_with3phases ? phi->getPointer(1) : phi->getPointer(0);
   const double* phiB = nullptr;
   if (d_with_phaseB) {
      phiB = d_with3phases ? phi->getPointer(2) : phi->getPointer(d_norderpA);
   }

   // create arrays of pointers and parameters to facilitate loop
   // iterations
   const std::vector<const double*> phi_data = {phiL, phiA, phiB};

   std::vector<std::shared_ptr<pdat::SideData<double>>> pfm_diffusion_sdata;
   pfm_diffusion_sdata.push_back(pfm_diffusion);

   const std::vector<int> norderp = {1, d_norderpA, d_norderpB};
   std::vector<std::vector<double>> d0;
   d0.resize(3);
   d0[0].push_back(0.);
   d0[1].push_back(d_d0_LA);
   d0[1].push_back(d_d0_AA);
   d0[2].push_back(d_d0_LB);
   d0[2].push_back(d_d0_AB);
   d0[2].push_back(d_d0_BB);

   // loop over phases L, A, B
   for (int i = 0; i < nphases; i++)
      for (int j = 0; j <= i; j++) {
         if (d0[i][j] > 0.) {
            assert(norderp[i] > 0);
            assert(norderp[j] > 0);
            // to avoid adding interface diffusion twice to an
            // interface between two grains of same phase, flag
            // case with dupl=1
            const int same_phase = (i == j) ? 1 : 0;
            const double dval = d0[i][j];
            // std::cout << "i=" << i << ", j=" << j << ", d = " << dval
            //          << std::endl;

            ADD_AB_DIFFUSION(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                             ifirst(2), ilast(2),
#endif
                             phi_data[i], norderp[i], phi_data[j], norderp[j],
                             phi->getGhostCellWidth()[0],
                             pfm_diffusion->getPointer(0, 0),
                             pfm_diffusion->getPointer(1, 0),
#if (NDIM == 3)
                             pfm_diffusion->getPointer(2, 0),
#endif
                             pfm_diffusion->getGhostCellWidth()[0], dval,
                             same_phase);
         }
      }
}

void ScalarCompositionDiffusionStrategy::setDiffusion(
    std::shared_ptr<hier::Patch> patch,
    std::shared_ptr<pdat::CellData<double>> phi,
    std::shared_ptr<pdat::SideData<double>> pfm_diffusion)
{
   const char interp_func_type = interpChar();


   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   // we need at least as many ghost cells for phi and T as in diffusion
   // to calculate diffusion in ghost cells
   assert(phi->getGhostCellWidth()[0] >= pfm_diffusion->getGhostCellWidth()[0]);

   // compute depth 0 of diffusion variables,
   // including phase fraction weight
   if (d_with_phaseB) {
      assert(phi->getDepth() == d_norderp);

      {
         // Folch-Plapp three phases model assumes order phiL, phiA, phiB
         double* phiL = d_with3phases ? phi->getPointer(0)
                                      : phi->getPointer(d_norderp - 1);
         double* phiA = d_with3phases ? phi->getPointer(1) : phi->getPointer(0);
         double* phiB =
             d_with3phases ? phi->getPointer(2) : phi->getPointer(d_norderpA);

         const int nphiL = 1;
         const int nphiA = d_norderpA;
         const int nphiB = d_norderpB;

         // this call assumes the order phiA, phiB, phiL
         PFMDIFFUSION_SCALAR_MULTIORDER_3PHASES(
             ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
             ifirst(2), ilast(2),
#endif
             phiL, nphiL, phiA, nphiA, phiB, nphiB, phi->getGhostCellWidth()[0],
             pfm_diffusion->getPointer(0, 0), pfm_diffusion->getPointer(1, 0),
#if (NDIM == 3)
             pfm_diffusion->getPointer(2, 0),
#endif
             pfm_diffusion->getGhostCellWidth()[0], d_D0_liquid, d_D0_solidA,
             d_D0_solidB, &interp_func_type, d_avg_func_type.c_str());
      }

      setDiffusionInterfaces(patch, phi, pfm_diffusion);

   } else if (phi->getDepth() > 1) {
      double* phiL =
          d_with3phases ? phi->getPointer(0) : phi->getPointer(d_norderp - 1);
      double* phiA = d_with3phases ? phi->getPointer(1) : phi->getPointer(0);
      const int nphiA = d_with3phases ? 1 : d_norderpA;

      PFMDIFFUSION_SCALAR_2PHASES(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                                  ifirst(2), ilast(2),
#endif
                                  phiL, phiA, nphiA,
                                  phi->getGhostCellWidth()[0],
                                  pfm_diffusion->getPointer(0, 0),
                                  pfm_diffusion->getPointer(1, 0),
#if (NDIM == 3)
                                  pfm_diffusion->getPointer(2, 0),
#endif
                                  pfm_diffusion->getGhostCellWidth()[0],
                                  d_D0_liquid, d_D0_solidA, &interp_func_type,
                                  d_avg_func_type.c_str());

      setDiffusionInterfaces(patch, phi, pfm_diffusion);

   } else {
      PFMDIFFUSION_SCALAR(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                          ifirst(2), ilast(2),
#endif
                          phi->getPointer(), phi->getGhostCellWidth()[0],
                          pfm_diffusion->getPointer(0, 0),
                          pfm_diffusion->getPointer(1, 0),
#if (NDIM == 3)
                          pfm_diffusion->getPointer(2, 0),
#endif
                          pfm_diffusion->getGhostCellWidth()[0], d_D0_liquid,
                          d_D0_solidA, &interp_func_type,
                          d_avg_func_type.c_str());
   }

   // fill other diagonal value with same value for ternaries for now
   if (pfm_diffusion->getDepth() > 1) {
      pfm_diffusion->copyDepth(3, *pfm_diffusion, 0);
   }
}
