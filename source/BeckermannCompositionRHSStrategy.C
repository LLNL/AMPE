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
#include "BeckermannCompositionRHSStrategy.h"

#include "QuatFort.h"
#include "QuatParams.h"
#include "ConcFort.h"
#include "FuncFort.h"
#include "QuatModel.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

#include <cassert>


BeckermannCompositionRHSStrategy::BeckermannCompositionRHSStrategy(
    QuatModel* quat_model, const int conc_scratch_id,
    const int phase_scratch_id, const int partition_coeff_scratch_id,
    const int conc_tilde_diffusion_id,
    const int conc_phase_coupling_diffusion_id, const double D_liquid,
    const double D_solid_A, const ConcInterpolationType phase_interp_func_type,
    const std::string& avg_func_type)
    : CompositionRHSStrategy(avg_func_type), d_quat_model(quat_model)
{
   assert(conc_scratch_id >= 0);
   assert(phase_scratch_id >= 0);
   assert(conc_tilde_diffusion_id >= 0);
   assert(conc_phase_coupling_diffusion_id >= 0);
   assert(partition_coeff_scratch_id >= 0);

   assert(D_liquid >= 0.);
   assert(D_solid_A >= 0.);

   d_phase_interp_func_type = phase_interp_func_type;
   d_avg_func_type = avg_func_type;

   d_D_liquid = D_liquid;
   d_D_solid_A = D_solid_A;

   d_conc_scratch_id = conc_scratch_id;
   d_phase_scratch_id = phase_scratch_id;

   d_diffusion0_id = conc_tilde_diffusion_id;
   d_conc_phase_coupling_diffusion_id = conc_phase_coupling_diffusion_id;

   d_partition_coeff_scratch_id = partition_coeff_scratch_id;

   tbox::TimerManager* tman = tbox::TimerManager::getManager();
   t_set_diffcoeff_timer = tman->getTimer(
       "AMPE::BeckermannCompositionRHSStrategy::setDiffusionCoeff()");
}

//-----------------------------------------------------------------------

void BeckermannCompositionRHSStrategy::setDiffusionCoeff(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time)
{
   (void)time;

   assert(hierarchy);

   t_set_diffcoeff_timer->start();

   // tbox::pout<<"BeckermannCompositionRHSStrategy::setDiffusionCoeffForConcentration"<<endl;

   // ghost values of k are needed to get k interpolated on faces
   // where diffusion coefficients are evaluated
   d_quat_model->fillPartitionCoeffGhosts();

   // set diffusion coefficients which is a function of the free energy form
   setDiffusionCoeffForConcentration(hierarchy, d_conc_scratch_id,
                                     d_phase_scratch_id, d_diffusion0_id,
                                     d_conc_phase_coupling_diffusion_id);

   t_set_diffcoeff_timer->stop();
}

//=======================================================================

void BeckermannCompositionRHSStrategy::setDiffusionCoeffForConcentration(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int concentration_id, const int phase_id,
    const int conc_tilde_diffusion_id,
    const int conc_phase_coupling_diffusion_id)
{
   // tbox::pout<<"BeckermannCompositionRHSStrategy::setDiffusionCoeffForConcentration()"<<endl;
   assert(concentration_id >= 0);
   assert(phase_id >= 0);
   assert(conc_tilde_diffusion_id >= 0);
   assert(conc_phase_coupling_diffusion_id >= 0);
   assert(d_partition_coeff_scratch_id >= 0);
   assert(d_D_liquid >= 0.);
   assert(d_D_solid_A >= 0.);

   const int maxl = hierarchy->getNumberOfLevels();

   // first compute diffusion for div*D0*grad*c term in concentration equation

   for (int amr_level = 0; amr_level < maxl; amr_level++) {
      std::shared_ptr<hier::PatchLevel> level =
          hierarchy->getPatchLevel(amr_level);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         const hier::Box& pbox = patch->getBox();
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast = pbox.upper();

         std::shared_ptr<pdat::CellData<double> > cd_phi(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(phase_id)));
         assert(cd_phi);

         std::shared_ptr<pdat::SideData<double> > sd_diffusion0(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(conc_tilde_diffusion_id)));
         assert(sd_diffusion0);

         std::shared_ptr<pdat::CellData<double> > cd_partition_coeff(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_partition_coeff_scratch_id)));
         assert(cd_partition_coeff);

         const char interp_func_type = concInterpChar(d_phase_interp_func_type);
         CONCENTRATIONDIFFUSION_BECKERMANN(
             ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
             ifirst(2), ilast(2),
#endif
             cd_phi->getPointer(), NGHOSTS, sd_diffusion0->getPointer(0),
             sd_diffusion0->getPointer(1),
#if (NDIM == 3)
             sd_diffusion0->getPointer(2),
#endif
             0, cd_partition_coeff->getPointer(),
             cd_partition_coeff->getGhostCellWidth()[0], d_D_liquid,
             d_D_solid_A, &interp_func_type, d_avg_func_type.c_str());
      }
   }

   // now set diffusion coefficients for div*D*grad*phi
   // term in r.h.s. of concentration equation
   for (int amr_level = 0; amr_level < maxl; amr_level++) {
      std::shared_ptr<hier::PatchLevel> level =
          hierarchy->getPatchLevel(amr_level);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         const hier::Box& pbox = patch->getBox();

         std::shared_ptr<pdat::SideData<double> > sd_phi_diff_coeff(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(conc_phase_coupling_diffusion_id)));

         std::shared_ptr<pdat::SideData<double> > sd_d0_coeff(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(conc_tilde_diffusion_id)));

         std::shared_ptr<pdat::CellData<double> > cd_phi(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(phase_id)));

         std::shared_ptr<pdat::CellData<double> > cd_c(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_conc_scratch_id)));

         std::shared_ptr<pdat::CellData<double> > cd_k(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_partition_coeff_scratch_id)));

         setDiffusionCoeffForPhaseOnPatch(sd_phi_diff_coeff, sd_d0_coeff,
                                          cd_phi, cd_c, cd_k, pbox);
      }
   }
}

//-----------------------------------------------------------------------

void BeckermannCompositionRHSStrategy::setDiffusionCoeffForPhaseOnPatch(
    std::shared_ptr<pdat::SideData<double> > sd_phi_diff_coeff,
    std::shared_ptr<pdat::SideData<double> > sd_d0_coeff,
    std::shared_ptr<pdat::CellData<double> > cd_phi,
    std::shared_ptr<pdat::CellData<double> > cd_c,
    std::shared_ptr<pdat::CellData<double> > cd_k, const hier::Box& pbox)
{
   double* ptr_phi_diffx = sd_phi_diff_coeff->getPointer(0);
   double* ptr_phi_diffy = sd_phi_diff_coeff->getPointer(1);
   double* ptr_phi_diffz = nullptr;

   double* ptr_dx0_coeff = sd_d0_coeff->getPointer(0);
   double* ptr_dy0_coeff = sd_d0_coeff->getPointer(1);
   double* ptr_dz0_coeff = nullptr;

   if (NDIM > 2) {
      ptr_phi_diffz = sd_phi_diff_coeff->getPointer(2);
      ptr_dz0_coeff = sd_d0_coeff->getPointer(2);
   }

   double* ptr_phi = cd_phi->getPointer();
   double* ptr_c = cd_c->getPointer();
   double* ptr_k = cd_k->getPointer();

   // Assuming d0_coeff and phi_diff_coeff have same box
   const hier::Box& dcoeff_gbox = sd_phi_diff_coeff->getGhostBox();
   int imin_dcoeff = dcoeff_gbox.lower(0);
   int jmin_dcoeff = dcoeff_gbox.lower(1);
   int jp_dcoeff = dcoeff_gbox.numberCells(0);
   int kmin_dcoeff = 0;
   int kp_dcoeff = 0;
#if (NDIM == 3)
   kmin_dcoeff = dcoeff_gbox.lower(2);
   kp_dcoeff = jp_dcoeff * dcoeff_gbox.numberCells(1);
#endif

   const hier::Box& pf_gbox = cd_phi->getGhostBox();
   int imin_pf = pf_gbox.lower(0);
   int jmin_pf = pf_gbox.lower(1);
   int jp_pf = pf_gbox.numberCells(0);
   int kmin_pf = 0;
   int kp_pf = 0;
#if (NDIM == 3)
   kmin_pf = pf_gbox.lower(2);
   kp_pf = jp_pf * pf_gbox.numberCells(1);
#endif

   const hier::Box& k_gbox = cd_k->getGhostBox();
   int imin_k = k_gbox.lower(0);
   int jmin_k = k_gbox.lower(1);
   int jp_k = k_gbox.numberCells(0);
   int kmin_k = 0;
   int kp_k = 0;
#if (NDIM == 3)
   kmin_k = k_gbox.lower(2);
   kp_k = jp_k * k_gbox.numberCells(1);
#endif

   const hier::Box& c_i_gbox = cd_k->getGhostBox();
   int imin_c_i = c_i_gbox.lower(0);
   int jmin_c_i = c_i_gbox.lower(1);
   int jp_c_i = c_i_gbox.numberCells(0);
   int kmin_c_i = 0;
   int kp_c_i = 0;
#if (NDIM == 3)
   kmin_c_i = c_i_gbox.lower(2);
   kp_c_i = jp_c_i * c_i_gbox.numberCells(1);
#endif

   int imin = pbox.lower(0);
   int imax = pbox.upper(0);
   int jmin = pbox.lower(1);
   int jmax = pbox.upper(1);
   int kmin = 0;
   int kmax = 0;
#if (NDIM == 3)
   kmin = pbox.lower(2);
   kmax = pbox.upper(2);
#endif
   const char interp_func_type = concInterpChar(d_phase_interp_func_type);

   // X-side
   for (int kk = kmin; kk <= kmax; kk++) {
      for (int jj = jmin; jj <= jmax; jj++) {
         for (int ii = imin; ii <= imax + 1; ii++) {

            const int idx_dcoeff =
                (ii - imin_dcoeff) + (jj - jmin_dcoeff) * (jp_dcoeff + 1) +
                (kk - kmin_dcoeff) * (kp_dcoeff + dcoeff_gbox.numberCells(1));

            const int idx_pf = (ii - imin_pf) + (jj - jmin_pf) * jp_pf +
                               (kk - kmin_pf) * kp_pf;
            const int idxm1_pf = idx_pf - 1;

            const int idx_c_i = (ii - imin_c_i) + (jj - jmin_c_i) * jp_c_i +
                                (kk - kmin_c_i) * kp_c_i;
            const int idxm1_c_i = idx_c_i - 1;

            const int idx_k =
                (ii - imin_k) + (jj - jmin_k) * jp_k + (kk - kmin_k) * kp_k;
            const int idxm1_k = idx_k - 1;

            double phi = average(ptr_phi[idx_pf], ptr_phi[idxm1_pf]);
            double hphi = INTERP_FUNC(phi, &interp_func_type);

            double c = average(ptr_c[idx_c_i], ptr_c[idxm1_c_i]);

            double k = average(ptr_k[idx_k], ptr_k[idxm1_k]);
            assert(k == k);
            assert(k >= 0.);
            assert(k <= 1.);

            ptr_phi_diffx[idx_dcoeff] = ptr_dx0_coeff[idx_dcoeff] * (1. - k) *
                                        c / (1. - hphi + k * hphi);
         }
      }
   }

   // Y-side
   for (int kk = kmin; kk <= kmax; kk++) {
      for (int jj = jmin; jj <= jmax + 1; jj++) {
         for (int ii = imin; ii <= imax; ii++) {

            const int idx_dcoeff = (ii - imin_dcoeff) +
                                   (jj - jmin_dcoeff) * jp_dcoeff +
                                   (kk - kmin_dcoeff) * (kp_dcoeff + jp_dcoeff);

            const int idx_pf = (ii - imin_pf) + (jj - jmin_pf) * jp_pf +
                               (kk - kmin_pf) * kp_pf;
            const int idxm1_pf = idx_pf - jp_pf;

            const int idx_c_i = (ii - imin_c_i) + (jj - jmin_c_i) * jp_c_i +
                                (kk - kmin_c_i) * kp_c_i;
            const int idxm1_c_i = idx_c_i - jp_c_i;

            const int idx_k =
                (ii - imin_k) + (jj - jmin_k) * jp_k + (kk - kmin_k) * kp_k;
            const int idxm1_k = idx_k - 1;

            // double phi = 0.5 * ( ptr_phi[idx_pf] + ptr_phi[idxm1_pf] );
            double phi = average(ptr_phi[idx_pf], ptr_phi[idxm1_pf]);

            double hphi = INTERP_FUNC(phi, &interp_func_type);

            double c = average(ptr_c[idx_c_i], ptr_c[idxm1_c_i]);

            double k = average(ptr_k[idx_k], ptr_k[idxm1_k]);
            assert(k == k);

            ptr_phi_diffy[idx_dcoeff] = ptr_dy0_coeff[idx_dcoeff] * (1. - k) *
                                        c / (1. - hphi + k * hphi);
         }
      }
   }

   if (NDIM > 2) {
      // Z-side
      for (int kk = kmin; kk <= kmax + 1; kk++) {
         for (int jj = jmin; jj <= jmax; jj++) {
            for (int ii = imin; ii <= imax; ii++) {

               const int idx_dcoeff = (ii - imin_dcoeff) +
                                      (jj - jmin_dcoeff) * jp_dcoeff +
                                      (kk - kmin_dcoeff) * kp_dcoeff;

               const int idx_pf = (ii - imin_pf) + (jj - jmin_pf) * jp_pf +
                                  (kk - kmin_pf) * kp_pf;
               const int idxm1_pf = idx_pf - kp_pf;

               const int idx_c_i = (ii - imin_c_i) + (jj - jmin_c_i) * jp_c_i +
                                   (kk - kmin_c_i) * kp_c_i;
               const int idxm1_c_i = idx_c_i - kp_c_i;

               const int idx_k =
                   (ii - imin_k) + (jj - jmin_k) * jp_k + (kk - kmin_k) * kp_k;
               const int idxm1_k = idx_k - 1;

               // double phi = 0.5 * ( ptr_phi[idx_pf] + ptr_phi[idxm1_pf] );
               double phi = average(ptr_phi[idx_pf], ptr_phi[idxm1_pf]);
               assert(phi == phi);

               double hphi = INTERP_FUNC(phi, &interp_func_type);

               double c = average(ptr_c[idx_c_i], ptr_c[idxm1_c_i]);

               double k = average(ptr_k[idx_k], ptr_k[idxm1_k]);
               assert(k == k);

               ptr_phi_diffz[idx_dcoeff] = ptr_dz0_coeff[idx_dcoeff] *
                                           (1. - k) * c /
                                           (1. - hphi + k * hphi);
            }
         }
      }
   }  // if ( NDIM > 2 )
}

//-----------------------------------------------------------------------

void BeckermannCompositionRHSStrategy::computeFluxOnPatch(hier::Patch& patch,
                                                          const int flux_id)
{
   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));
   const double* dx = patch_geom->getDx();

   const hier::Box& pbox = patch.getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   std::shared_ptr<pdat::CellData<double> > conc(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_scratch_id)));
   assert(conc);

   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_phase_scratch_id)));
   assert(phase);

   std::shared_ptr<pdat::SideData<double> > sd_conc_diffusion0(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(d_diffusion0_id)));
   assert(sd_conc_diffusion0);

   std::shared_ptr<pdat::SideData<double> > sd_conc_phase_coupling_diffusion(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_phase_coupling_diffusion_id)));
   assert(sd_conc_phase_coupling_diffusion);

   std::shared_ptr<pdat::SideData<double> > flux(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(flux_id)));
   assert(flux);

   int three_phase = 0;

#ifdef DEBUG_CHECK_ASSERTIONS
   SAMRAI::math::PatchSideDataNormOpsReal<double> ops;
   double l2d = ops.L2Norm(sd_conc_diffusion0, pbox);
   assert(l2d == l2d);

   double l2p = ops.L2Norm(sd_conc_phase_coupling_diffusion, pbox);
   assert(l2p == l2p);
#endif

   // now compute concentration flux
   CONCENTRATIONFLUX(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                     ifirst(2), ilast(2),
#endif
                     dx, conc->getPointer(), NGHOSTS, phase->getPointer(),
                     NGHOSTS, 0, NGHOSTS, sd_conc_diffusion0->getPointer(0),
                     sd_conc_diffusion0->getPointer(1),
#if (NDIM == 3)
                     sd_conc_diffusion0->getPointer(2),
#endif
                     0, sd_conc_phase_coupling_diffusion->getPointer(0),
                     sd_conc_phase_coupling_diffusion->getPointer(1),
#if (NDIM == 3)
                     sd_conc_phase_coupling_diffusion->getPointer(2),
#endif
                     0, nullptr, nullptr,
#if (NDIM == 3)
                     nullptr,
#endif
                     0, flux->getPointer(0), flux->getPointer(1),
#if (NDIM == 3)
                     flux->getPointer(2),
#endif
                     flux->getGhostCellWidth()[0], three_phase);

#ifdef DEBUG_CHECK_ASSERTIONS
   double l2f = ops.L2Norm(flux, pbox);
   assert(l2f == l2f);
#endif
}
