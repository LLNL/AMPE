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

const double gas_constant_R_JpKpmol = GASCONSTANT_R_JPKPMOL;

TbasedCompositionDiffusionStrategy::TbasedCompositionDiffusionStrategy(
    const short norderp, const short norderpA, const short norderpB,
    const bool with3phases, const int pfm_diffusion_l_id,
    const int pfm_diffusion_a_id, const int pfm_diffusion_b_id,
    const double D0_liquid, const double Q0_liquid, const double D0_solidA,
    const double Q0_solidA, const double D0_solidB, const double Q0_solidB,
    const double D0_LA, const double Q0_LA, const double D0_LB,
    const double Q0_LB, const double D0_AA, const double Q0_AA,
    const double D0_AB, const double Q0_AB, const double D0_BB,
    const double Q0_BB, const DiffusionInterpolationType interp_func_type,
    const std::string& avg_func_type)
    : CompositionDiffusionStrategy(interp_func_type),
      d_norderp(norderp),
      d_norderpA(norderpA),
      d_norderpB(norderpB),
      d_with3phases(with3phases),
      d_pfm_diffusion_l_id(pfm_diffusion_l_id),
      d_pfm_diffusion_a_id(pfm_diffusion_a_id),
      d_pfm_diffusion_b_id(pfm_diffusion_b_id),
      d_D0_liquid(D0_liquid),
      d_Q0_liquid(Q0_liquid),
      d_D0_solidA(D0_solidA),
      d_Q0_solidA(Q0_solidA),
      d_D0_solidB(D0_solidB),
      d_Q0_solidB(Q0_solidB),
      d_d0_LA(D0_LA),
      d_q0_LA(Q0_LA),
      d_d0_LB(D0_LB),
      d_q0_LB(Q0_LB),
      d_d0_AA(D0_AA),
      d_q0_AA(Q0_AA),
      d_d0_AB(D0_AB),
      d_q0_AB(Q0_AB),
      d_d0_BB(D0_BB),
      d_q0_BB(Q0_BB),
      d_avg_func_type(avg_func_type),
      d_with_phaseB(d_pfm_diffusion_b_id >= 0)
{
   assert(D0_liquid >= 0.);
   assert(Q0_liquid >= 0.);
   assert(Q0_solidA >= 0.);
   assert(D0_solidA >= 0.);

   tbox::plog << "TbasedCompositionDiffusionStrategy: D0_AB = " << D0_AB
              << std::endl;
}

void TbasedCompositionDiffusionStrategy::setDiffusion(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int phase_id)
{
   // tbox::plog<<"TbasedCompositionDiffusionStrategy::setDiffusion()"<<std::endl;
   assert(temperature_id >= 0);
   assert(phase_id >= 0);
   assert(d_pfm_diffusion_l_id >= 0);
   assert(d_pfm_diffusion_a_id >= 0);

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

         std::shared_ptr<pdat::CellData<double>> temperature(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(temperature_id)));

         std::shared_ptr<pdat::SideData<double>> pfm_diffusionL(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(d_pfm_diffusion_l_id)));

         std::shared_ptr<pdat::SideData<double>> pfm_diffusionA(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(d_pfm_diffusion_a_id)));

         std::shared_ptr<pdat::SideData<double>> pfm_diffusionB;
         if (d_with_phaseB) {
            pfm_diffusionB =
                SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                    patch->getPatchData(d_pfm_diffusion_b_id));
         }
         assert(pfm_diffusionA->getDepth() == pfm_diffusionL->getDepth());
         assert(pfm_diffusionA->getDepth() == 1 ||  // binary
                pfm_diffusionA->getDepth() == 4);   // ternary

         setDiffusion(patch, phi, temperature, pfm_diffusionL, pfm_diffusionA,
                      pfm_diffusionB);
      }
   }
}

void TbasedCompositionDiffusionStrategy::setDiffusionInterfaces(
    std::shared_ptr<hier::Patch> patch,
    std::shared_ptr<pdat::CellData<double>> phi,
    std::shared_ptr<pdat::CellData<double>> temperature,
    std::shared_ptr<pdat::SideData<double>> pfm_diffusionL,
    std::shared_ptr<pdat::SideData<double>> pfm_diffusionA,
    std::shared_ptr<pdat::SideData<double>> pfm_diffusionB)
{
   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   const short nphases = d_with_phaseB ? 3 : 2;

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
   pfm_diffusion_sdata.push_back(pfm_diffusionL);
   pfm_diffusion_sdata.push_back(pfm_diffusionA);
   pfm_diffusion_sdata.push_back(pfm_diffusionB);

   const std::vector<int> norderp = {1, d_norderpA, d_norderpB};
   std::vector<std::vector<double>> d0;
   d0.resize(3);
   d0[0].push_back(0.);
   d0[1].push_back(d_d0_LA);
   d0[1].push_back(d_d0_AA);
   d0[2].push_back(d_d0_LB);
   d0[2].push_back(d_d0_AB);
   d0[2].push_back(d_d0_BB);

   std::vector<std::vector<double>> q0;
   q0.resize(3);
   q0[0].push_back(0.);
   q0[1].push_back(d_q0_LA);
   q0[1].push_back(d_q0_AA);
   q0[2].push_back(d_q0_LB);
   q0[2].push_back(d_q0_AB);
   q0[2].push_back(d_q0_BB);
   for (int i = 0; i < nphases; i++)
      for (int j = 0; j <= i; j++) {
         if (d0[i][j] > 0.) {
            assert(norderp[i] > 0);
            assert(norderp[j] > 0);
            // to avoid adding interface diffusion twice to an
            // interface between two grains of same phase, add factor
            // 0.5
            const double dval = (i == j) ? 0.5 * d0[i][j] : d0[i][j];
            // std::cout << "i=" << i << ", j=" << j << ", d = " << dval
            //          << std::endl;
            AB_DIFFUSION_OF_TEMPERATURE(
                ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                ifirst(2), ilast(2),
#endif
                phi_data[i], norderp[i], phi_data[j], norderp[j],
                phi->getGhostCellWidth()[0],
                pfm_diffusion_sdata[i]->getPointer(0, 0),
                pfm_diffusion_sdata[i]->getPointer(1, 0),
#if (NDIM == 3)
                pfm_diffusion_sdata[i]->getPointer(2, 0),
#endif
                pfm_diffusion_sdata[j]->getPointer(0, 0),
                pfm_diffusion_sdata[j]->getPointer(1, 0),
#if (NDIM == 3)
                pfm_diffusion_sdata[j]->getPointer(2, 0),
#endif
                pfm_diffusionL->getGhostCellWidth()[0],
                temperature->getPointer(), temperature->getGhostCellWidth()[0],
                dval, q0[i][j], gas_constant_R_JpKpmol,
                d_avg_func_type.c_str());
         }
      }
}

void TbasedCompositionDiffusionStrategy::setDiffusion(
    std::shared_ptr<hier::Patch> patch,
    std::shared_ptr<pdat::CellData<double>> phi,
    std::shared_ptr<pdat::CellData<double>> temperature,
    std::shared_ptr<pdat::SideData<double>> pfm_diffusionL,
    std::shared_ptr<pdat::SideData<double>> pfm_diffusionA,
    std::shared_ptr<pdat::SideData<double>> pfm_diffusionB)
{
   const char interp_func_type = interpChar();


   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   // we need at least as many ghost cells for phi and T as in diffusion
   // to calculate diffusion in ghost cells
   assert(phi->getGhostCellWidth()[0] >=
          pfm_diffusionL->getGhostCellWidth()[0]);
   assert(temperature->getGhostCellWidth()[0] >=
          pfm_diffusionL->getGhostCellWidth()[0]);

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
         const int nphiA = d_with3phases ? 1 : d_norderpA;
         const int nphiB = d_with3phases ? 1 : d_norderpB;

         // this call assumes the order phiA, phiB, phiL
         CONCENTRATION_PFMDIFFUSION_OF_TEMPERATURE_THREEPHASES(
             ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
             ifirst(2), ilast(2),
#endif
             phiL, nphiL, phiA, nphiA, phiB, nphiB, phi->getGhostCellWidth()[0],
             pfm_diffusionL->getPointer(0, 0), pfm_diffusionL->getPointer(1, 0),
#if (NDIM == 3)
             pfm_diffusionL->getPointer(2, 0),
#endif
             pfm_diffusionA->getPointer(0, 0), pfm_diffusionA->getPointer(1, 0),
#if (NDIM == 3)
             pfm_diffusionA->getPointer(2, 0),
#endif
             pfm_diffusionB->getPointer(0, 0), pfm_diffusionB->getPointer(1, 0),
#if (NDIM == 3)
             pfm_diffusionB->getPointer(2, 0),
#endif
             pfm_diffusionL->getGhostCellWidth()[0], temperature->getPointer(),
             temperature->getGhostCellWidth()[0], d_D0_liquid, d_Q0_liquid,
             d_D0_solidA, d_Q0_solidA, d_D0_solidB, d_Q0_solidB,
             gas_constant_R_JpKpmol, &interp_func_type,
             d_avg_func_type.c_str());
      }

      setDiffusionInterfaces(patch, phi, temperature, pfm_diffusionL,
                             pfm_diffusionA, pfm_diffusionB);

   } else if (phi->getDepth() > 1) {

      CONCENTRATION_PFMDIFFUSION_OF_TEMPERATURE_MULTIPHASES(
          ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
          ifirst(2), ilast(2),
#endif
          phi->getPointer(), phi->getDepth(), phi->getGhostCellWidth()[0],
          pfm_diffusionL->getPointer(0, 0), pfm_diffusionL->getPointer(1, 0),
#if (NDIM == 3)
          pfm_diffusionL->getPointer(2, 0),
#endif
          pfm_diffusionA->getPointer(0, 0), pfm_diffusionA->getPointer(1, 0),
#if (NDIM == 3)
          pfm_diffusionA->getPointer(2, 0),
#endif
          pfm_diffusionL->getGhostCellWidth()[0], temperature->getPointer(),
          temperature->getGhostCellWidth()[0], d_D0_liquid, d_Q0_liquid,
          d_D0_solidA, d_Q0_solidA, gas_constant_R_JpKpmol, &interp_func_type,
          d_avg_func_type.c_str());

      setDiffusionInterfaces(patch, phi, temperature, pfm_diffusionL,
                             pfm_diffusionA, pfm_diffusionB);

   } else {
      CONCENTRATION_PFMDIFFUSION_OF_TEMPERATURE(
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
          pfm_diffusionL->getGhostCellWidth()[0], temperature->getPointer(),
          temperature->getGhostCellWidth()[0], d_D0_liquid, d_Q0_liquid,
          d_D0_solidA, d_Q0_solidA, gas_constant_R_JpKpmol, &interp_func_type,
          d_avg_func_type.c_str());
   }

   // fill other diagonal value with same value for ternaries for now
   if (pfm_diffusionL->getDepth() > 1) {
      pfm_diffusionL->copyDepth(3, *pfm_diffusionL, 0);
      pfm_diffusionA->copyDepth(3, *pfm_diffusionA, 0);
      if (d_with_phaseB) pfm_diffusionB->copyDepth(3, *pfm_diffusionB, 0);
   }
}
