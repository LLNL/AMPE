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
#include "KKSCompositionRHSStrategy.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "QuatFort.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include <cassert>
#include "QuatParams.h"
#include "ConcFort.h"
#include "FuncFort.h"
#include "SAMRAI/tbox/TimerManager.h"

KKSCompositionRHSStrategy::KKSCompositionRHSStrategy(
    const int conc_scratch_id, const int phase_scratch_id,
    const int conc_pfm_diffusion_id, const int conc_phase_coupling_diffusion_id,
    const int temperature_scratch_id, const int eta_scratch_id,
    const int conc_eta_coupling_diffusion_id, const int conc_l_scratch_id,
    const int conc_a_scratch_id, const int conc_b_scratch_id,
    const double D_liquid, const double D_solid_A, const double D_solid_B,
    const double Q0_liquid, const double Q0_solid_A, const double Q0_solid_B,
    const EnergyInterpolationType phase_interp_func_type,
    const std::string& avg_func_type)
    : CompositionRHSStrategy(avg_func_type)
{
   assert(conc_scratch_id >= 0);
   assert(phase_scratch_id >= 0);
   assert(conc_pfm_diffusion_id >= 0);
   assert(conc_phase_coupling_diffusion_id >= 0);
   assert(temperature_scratch_id >= 0);

   assert(D_liquid >= 0.);
   assert(Q0_liquid >= 0.);
   assert(Q0_solid_A >= 0.);
   assert(D_solid_A >= 0.);

   if (eta_scratch_id >= 0) {
      assert(Q0_solid_B >= 0.);
      assert(D_solid_B >= 0.);
   }

   d_phase_interp_func_type = phase_interp_func_type;
   d_avg_func_type = avg_func_type;

   d_conc_l_scratch_id = conc_l_scratch_id;
   d_conc_a_scratch_id = conc_a_scratch_id;
   d_conc_b_scratch_id = conc_b_scratch_id;

   d_D_liquid = D_liquid;
   d_D_solid_A = D_solid_A;
   d_D_solid_B = D_solid_B;

   d_Q0_liquid = Q0_liquid;
   d_Q0_solid_A = Q0_solid_A;
   d_Q0_solid_B = Q0_solid_B;

   d_conc_scratch_id = conc_scratch_id;
   d_phase_scratch_id = phase_scratch_id;

   d_pfm_diffusion_id = conc_pfm_diffusion_id;
   d_conc_phase_coupling_diffusion_id = conc_phase_coupling_diffusion_id;

   d_temperature_scratch_id = temperature_scratch_id;

   d_eta_scratch_id = eta_scratch_id;
   d_conc_eta_coupling_diffusion_id = conc_eta_coupling_diffusion_id;

   d_with_third_phase = (eta_scratch_id >= 0);

   tbox::TimerManager* tman = tbox::TimerManager::getManager();
   t_set_diffcoeff_timer =
       tman->getTimer("AMPE::KKSCompositionRHSStrategy::setDiffusionCoeff()");
}

//-----------------------------------------------------------------------

void KKSCompositionRHSStrategy::setDiffusionCoeff(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time)
{
   (void)time;

   t_set_diffcoeff_timer->start();

   // tbox::pout<<"QuatIntegrator::setDiffCoeffForConcentration"<<endl;

   assert(hierarchy);

   // set diffusion coefficients which is a function of the free energy form
   setPFMDiffCoeffForConcentration(hierarchy, d_temperature_scratch_id,
                                   d_phase_scratch_id, d_eta_scratch_id,
                                   d_pfm_diffusion_id);
   setDiffCoeffForGradPhi(hierarchy, d_temperature_scratch_id,
                          d_conc_scratch_id, d_phase_scratch_id,
                          d_eta_scratch_id, d_pfm_diffusion_id,
                          d_conc_phase_coupling_diffusion_id,
                          d_conc_eta_coupling_diffusion_id);

   t_set_diffcoeff_timer->stop();
}

//=======================================================================
// compute diffusion for div*D0*grad*c term in concentration equation
void KKSCompositionRHSStrategy::setPFMDiffCoeffForConcentration(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int phase_id, const int eta_id,
    const int conc_pfm_diffusion_id)
{
   // tbox::pout<<"KKSCompositionRHSStrategy::setDiffCoeffForConcentration()"<<endl;
   assert(temperature_id >= 0);
   assert(phase_id >= 0);
   assert(conc_pfm_diffusion_id >= 0);
   if (d_with_third_phase) {
      assert(eta_id >= 0);
   }

   const int maxl = hierarchy->getNumberOfLevels();

   const char interpf = energyInterpChar(d_phase_interp_func_type);

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

         std::shared_ptr<pdat::SideData<double> > pfm_diffusion(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(conc_pfm_diffusion_id)));
         assert(pfm_diffusion->getDepth() == 1);

         int three_phase = 0;
         double* ptr_eta = nullptr;
         std::shared_ptr<pdat::CellData<double> > eta;
         if (d_with_third_phase) {
            three_phase = 1;
            eta = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(eta_id));
            ptr_eta = eta->getPointer();
         }

         FORT_CONCENTRATION_PFMDIFFUSION(
             ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
             ifirst(2), ilast(2),
#endif
             phi->getPointer(), NGHOSTS, ptr_eta, NGHOSTS,
             pfm_diffusion->getPointer(0), pfm_diffusion->getPointer(1),
#if (NDIM == 3)
             pfm_diffusion->getPointer(2),
#endif
             0, temperature->getPointer(), temperature->getGhostCellWidth()[0],
             d_D_liquid, d_Q0_liquid, d_D_solid_A, d_Q0_solid_A, d_D_solid_B,
             d_Q0_solid_B, gas_constant_R_JpKpmol, &interpf,
             d_avg_func_type.c_str(), three_phase);
      }
   }
}

//=======================================================================

void KKSCompositionRHSStrategy::setDiffCoeffForGradPhi(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int concentration_id, const int phase_id,
    const int eta_id, const int conc_pfm_diffusion_id,
    const int conc_phase_coupling_diffusion_id,
    const int conc_eta_coupling_diffusion_id)
{
   (void)temperature_id;
   (void)concentration_id;

   assert(phase_id >= 0);
   assert(conc_pfm_diffusion_id >= 0);
   assert(conc_phase_coupling_diffusion_id >= 0);

   // now set diffusion coefficients for div*D*grad*phi term in concentration
   // equation

   const int maxl = hierarchy->getNumberOfLevels();

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

         std::shared_ptr<pdat::SideData<double> > sd_eta_diff_coeff;

         std::shared_ptr<pdat::SideData<double> > sd_pfmd_coeff(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(conc_pfm_diffusion_id)));

         std::shared_ptr<pdat::CellData<double> > cd_phi(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(phase_id)));

         std::shared_ptr<pdat::CellData<double> > cd_eta;

         std::shared_ptr<pdat::CellData<double> > cd_c_l(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_conc_l_scratch_id)));

         std::shared_ptr<pdat::CellData<double> > cd_c_a(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_conc_a_scratch_id)));

         std::shared_ptr<pdat::CellData<double> > cd_c_b;

         if (d_with_third_phase) {
            cd_eta = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(eta_id));
            cd_c_b = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(d_conc_b_scratch_id));
            sd_eta_diff_coeff =
                SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                    patch->getPatchData(conc_eta_coupling_diffusion_id));
         }

         setDiffCoeffForPhaseOnPatch(sd_phi_diff_coeff, sd_eta_diff_coeff,
                                     sd_pfmd_coeff, cd_phi, cd_eta, cd_c_l,
                                     cd_c_a, cd_c_b, pbox);
      }
   }
}

//-----------------------------------------------------------------------

void KKSCompositionRHSStrategy::setDiffCoeffForPhaseOnPatch(
    std::shared_ptr<pdat::SideData<double> > sd_phi_diff_coeff,
    std::shared_ptr<pdat::SideData<double> > sd_eta_diff_coeff,
    std::shared_ptr<pdat::SideData<double> > sd_d0_coeff,
    std::shared_ptr<pdat::CellData<double> > cd_phi,
    std::shared_ptr<pdat::CellData<double> > cd_eta,
    std::shared_ptr<pdat::CellData<double> > cd_c_l,
    std::shared_ptr<pdat::CellData<double> > cd_c_a,
    std::shared_ptr<pdat::CellData<double> > cd_c_b, const hier::Box& pbox)
{
   double* ptr_phi_diffx = sd_phi_diff_coeff->getPointer(0);
   double* ptr_phi_diffy = sd_phi_diff_coeff->getPointer(1);
   double* ptr_phi_diffz = nullptr;
   double* ptr_eta_diffx = nullptr;
   double* ptr_eta_diffy = nullptr;
   double* ptr_eta_diffz = nullptr;

   double* ptr_dx0_coeff = sd_d0_coeff->getPointer(0);
   double* ptr_dy0_coeff = sd_d0_coeff->getPointer(1);
   double* ptr_dz0_coeff = nullptr;

   if (NDIM > 2) {
      ptr_phi_diffz = sd_phi_diff_coeff->getPointer(2);
      ptr_dz0_coeff = sd_d0_coeff->getPointer(2);
   }

   double* ptr_phi = cd_phi->getPointer();
   double* ptr_eta = nullptr;
   double* ptr_c_l = cd_c_l->getPointer();
   double* ptr_c_a = cd_c_a->getPointer();
   double* ptr_c_b = nullptr;

   if (d_with_third_phase) {
      ptr_eta_diffx = sd_eta_diff_coeff->getPointer(0);
      ptr_eta_diffy = sd_eta_diff_coeff->getPointer(1);
      if (NDIM > 2) {
         ptr_eta_diffz = sd_eta_diff_coeff->getPointer(2);
      }
      ptr_eta = cd_eta->getPointer();
      ptr_c_b = cd_c_b->getPointer();
   }

   // Assuming d0_coeff, phi_diff_coeff, and eta_diff_coeff have same box
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

   // Assuming phi and eta have same box
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

   // Assuming c_l, c_a, and c_b all have same box
   const hier::Box& c_i_gbox = cd_c_l->getGhostBox();
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

   const char interpf = energyInterpChar(d_phase_interp_func_type);

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

            // double phi = 0.5 * ( ptr_phi[idx_pf] + ptr_phi[idxm1_pf] );
            double phi = average(ptr_phi[idx_pf], ptr_phi[idxm1_pf]);

            double eta = 0.0;
            double c_l = 0.5 * (ptr_c_l[idx_c_i] + ptr_c_l[idxm1_c_i]);
            double c_a = 0.5 * (ptr_c_a[idx_c_i] + ptr_c_a[idxm1_c_i]);
            double c_b = 0.0;

            double hphi_prime = FORT_DERIV_INTERP_FUNC(phi, &interpf);

            double heta = 0.0;

            if (d_with_third_phase) {
               eta = 0.5 * (ptr_eta[idx_pf] + ptr_eta[idxm1_pf]);
               c_b = 0.5 * (ptr_c_b[idx_c_i] + ptr_c_b[idxm1_c_i]);

               heta = FORT_INTERP_FUNC(eta, &interpf);
            }

            ptr_phi_diffx[idx_dcoeff] = ptr_dx0_coeff[idx_dcoeff] * hphi_prime *
                                        (c_l - (1.0 - heta) * c_a - heta * c_b);

            if (d_with_third_phase) {
               double hphi = FORT_INTERP_FUNC(phi, &interpf);

               double heta_prime = FORT_DERIV_INTERP_FUNC(eta, &interpf);

               ptr_eta_diffx[idx_dcoeff] =
                   ptr_dx0_coeff[idx_dcoeff] * hphi * heta_prime * (c_a - c_b);
            }
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

            // double phi = 0.5 * ( ptr_phi[idx_pf] + ptr_phi[idxm1_pf] );
            double phi = average(ptr_phi[idx_pf], ptr_phi[idxm1_pf]);
            double eta = 0.0;
            double c_l = 0.5 * (ptr_c_l[idx_c_i] + ptr_c_l[idxm1_c_i]);
            double c_a = 0.5 * (ptr_c_a[idx_c_i] + ptr_c_a[idxm1_c_i]);
            double c_b = 0.0;

            double hphi_prime = FORT_DERIV_INTERP_FUNC(phi, &interpf);

            double heta = 0.0;

            if (d_with_third_phase) {
               eta = 0.5 * (ptr_eta[idx_pf] + ptr_eta[idxm1_pf]);
               c_b = 0.5 * (ptr_c_b[idx_c_i] + ptr_c_b[idxm1_c_i]);

               heta = FORT_INTERP_FUNC(eta, &interpf);
            }

            ptr_phi_diffy[idx_dcoeff] = ptr_dy0_coeff[idx_dcoeff] * hphi_prime *
                                        (c_l - (1.0 - heta) * c_a - heta * c_b);

            if (d_with_third_phase) {
               double hphi = FORT_INTERP_FUNC(phi, &interpf);

               double heta_prime = FORT_DERIV_INTERP_FUNC(eta, &interpf);

               ptr_eta_diffy[idx_dcoeff] =
                   ptr_dy0_coeff[idx_dcoeff] * hphi * heta_prime * (c_a - c_b);
            }
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

               // double phi = 0.5 * ( ptr_phi[idx_pf] + ptr_phi[idxm1_pf] );
               double phi = average(ptr_phi[idx_pf], ptr_phi[idxm1_pf]);
               double eta = 0.0;
               double c_l = 0.5 * (ptr_c_l[idx_c_i] + ptr_c_l[idxm1_c_i]);
               double c_a = 0.5 * (ptr_c_a[idx_c_i] + ptr_c_a[idxm1_c_i]);
               double c_b = 0.0;

               double hphi_prime = FORT_DERIV_INTERP_FUNC(phi, &interpf);

               double heta = 0.0;

               if (d_with_third_phase) {
                  eta = 0.5 * (ptr_eta[idx_pf] + ptr_eta[idxm1_pf]);
                  c_b = 0.5 * (ptr_c_b[idx_c_i] + ptr_c_b[idxm1_c_i]);

                  heta = FORT_INTERP_FUNC(eta, &interpf);
               }

               ptr_phi_diffz[idx_dcoeff] =
                   ptr_dz0_coeff[idx_dcoeff] * hphi_prime *
                   (c_l - (1.0 - heta) * c_a - heta * c_b);

               if (d_with_third_phase) {
                  double hphi = FORT_INTERP_FUNC(phi, &interpf);

                  double heta_prime = FORT_DERIV_INTERP_FUNC(eta, &interpf);

                  ptr_eta_diffz[idx_dcoeff] = ptr_dz0_coeff[idx_dcoeff] * hphi *
                                              heta_prime * (c_a - c_b);
               }
            }
         }
      }
   }  // if ( NDIM > 2 )
}

//=======================================================================

#if 0
void KKSCompositionRHSStrategy::computeDiffusionOnPatch(
   std::shared_ptr< pdat::CellData<double> > cd_temperature,
   std::shared_ptr< pdat::CellData<double> > cd_phi,
   std::shared_ptr< pdat::CellData<double> > cd_eta,
   std::shared_ptr< pdat::CellData<double> > cd_concentration,
   std::shared_ptr< pdat::SideData<double> > sd_diffusion0,
   std::shared_ptr< pdat::SideData<double> > sd_diffusion,
   const hier::Box& pbox )
{
   assert( cd_temperature );
   assert( cd_phi );
   assert( cd_concentration );
   assert( sd_diffusion0 );
   assert( sd_diffusion );
   
   //tbox::pout<<"CALPHADFreeEnergyStrategy::computeDiffusionOnPatch()"<<endl;

   const double* const ptr_temp = cd_temperature->getPointer();
   const double* const ptr_phi  = cd_phi->getPointer();
   const double* const ptr_conc = cd_concentration->getPointer();

   double* ptr_diffusion0_0 = sd_diffusion0->getPointer(0);
   double* ptr_diffusion0_1 = sd_diffusion0->getPointer(1);
   double* ptr_diffusion_0  = sd_diffusion->getPointer(0);
   double* ptr_diffusion_1  = sd_diffusion->getPointer(1);
#if (NDIM == 3)
   double* ptr_diffusion0_2 = sd_diffusion0->getPointer(2);
   double* ptr_diffusion_2  = sd_diffusion->getPointer(2);
#endif
   double* ptr_eta = nullptr;
   if ( d_with_third_phase ) {
      assert( cd_eta );
      ptr_eta = cd_eta->getPointer();
   }

   const hier::Box& temp_gbox = cd_temperature->getGhostBox();
   int imin_temp = temp_gbox.lower(0);
   int jmin_temp = temp_gbox.lower(1);
   int jp_temp = temp_gbox.numberCells(0);
   int kmin_temp = 0;
   int kp_temp = 0;
#if (NDIM == 3)
   kmin_temp = temp_gbox.lower(2);
   kp_temp = jp_temp * temp_gbox.numberCells(1);
#endif

   // Assuming phi, eta, and concentration all have same box
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

   const hier::Box& dcoeff_gbox = sd_diffusion->getGhostBox();
   int imin_dcoeff = dcoeff_gbox.lower(0);
   int jmin_dcoeff = dcoeff_gbox.lower(1);
   int jp_dcoeff = dcoeff_gbox.numberCells(0);
   int kmin_dcoeff = 0;
   int kp_dcoeff = 0;
#if (NDIM == 3)
   kmin_dcoeff = dcoeff_gbox.lower(2);
   kp_dcoeff = jp_dcoeff * dcoeff_gbox.numberCells(1);
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
   const char interpf = energyInterpChar(d_phase_interp_func_type); 
   // X-side
   for ( int kk = kmin; kk <= kmax; kk++ ) {
      for ( int jj = jmin; jj <= jmax; jj++ ) {
         for ( int ii = imin; ii <= imax+1; ii++ ) {

            const int idx_dcoeff = (ii - imin_dcoeff) +
               (jj - jmin_dcoeff) * (jp_dcoeff + 1) +
               (kk - kmin_dcoeff) * (kp_dcoeff + dcoeff_gbox.numberCells(1) );

            const int idx_pf = (ii - imin_pf) +
               (jj - jmin_pf) * jp_pf + (kk - kmin_pf) * kp_pf;
            const int idxm1_pf = idx_pf - 1;

            const int idx_temp = (ii - imin_temp) +
               (jj - jmin_temp) * jp_temp + (kk - kmin_temp) * kp_temp;
            const int idxm1_temp = idx_temp - 1;

            double phi = 
               average( ptr_phi[idx_pf], ptr_phi[idxm1_pf] );
            double c = 
               average( ptr_conc[idx_pf], ptr_conc[idxm1_pf] );
            double t = 
               average( ptr_temp[idx_temp], ptr_temp[idxm1_temp] );
            
            const double hphi =
               FORT_INTERP_FUNC(phi, &interpf);

            double heta = 0.0;
            if ( d_with_third_phase ) {
               double eta = 
                  average( ptr_eta[idx_pf], ptr_eta[idxm1_pf] );
               heta =
                  FORT_INTERP_FUNC(eta, &interpf);
            }

            const double invd2fdc2 = computeLocalInvD2fDc2(c, hphi, heta, t);
            // D=D0*(d2f/dc2)^-1
            ptr_diffusion_0[idx_dcoeff]
               = ptr_diffusion0_0[idx_dcoeff]*invd2fdc2;
         }
      }
   }

   // Y-side
   for ( int kk = kmin; kk <= kmax; kk++ ) {
      for ( int jj = jmin; jj <= jmax+1; jj++ ) {
         for ( int ii = imin; ii <= imax; ii++ ) {

            const int idx_dcoeff = (ii - imin_dcoeff) +
               (jj - jmin_dcoeff) * jp_dcoeff +
               (kk - kmin_dcoeff) * (kp_dcoeff + jp_dcoeff);

            const int idx_pf = (ii - imin_pf) +
               (jj - jmin_pf) * jp_pf + (kk - kmin_pf) * kp_pf;
            const int idxm1_pf = idx_pf - jp_pf;

            const int idx_temp = (ii - imin_temp) +
               (jj - jmin_temp) * jp_temp + (kk - kmin_temp) * kp_temp;
            const int idxm1_temp = idx_temp - jp_temp;

            double phi = 
               average( ptr_phi[idx_pf], ptr_phi[idxm1_pf] );
            double c = 
               average( ptr_conc[idx_pf], ptr_conc[idxm1_pf] );
            double t = 
               average( ptr_temp[idx_temp], ptr_temp[idxm1_temp] );
            
            const double hphi =
               FORT_INTERP_FUNC(phi, &interpf);

            double heta = 0.0;
            if ( d_with_third_phase ) {
               double eta = 
                  average( ptr_eta[idx_pf], ptr_eta[idxm1_pf] );
               heta =
                  FORT_INTERP_FUNC(eta, &interpf);
            }

            const double invd2fdc2 = computeLocalInvD2fDc2(c, hphi, heta, t);
            // D=D0*(d2f/dc2)^-1
            ptr_diffusion_1[idx_dcoeff]
               = ptr_diffusion0_1[idx_dcoeff]*invd2fdc2;
         }
      }
   }

#if (NDIM == 3)
   // Z-side
   for ( int kk = kmin; kk <= kmax+1; kk++ ) {
      for ( int jj = jmin; jj <= jmax; jj++ ) {
         for ( int ii = imin; ii <= imax; ii++ ) {

            const int idx_dcoeff = (ii - imin_dcoeff) +
               (jj - jmin_dcoeff) * jp_dcoeff +
               (kk - kmin_dcoeff) * kp_dcoeff;

            const int idx_pf = (ii - imin_pf) +
               (jj - jmin_pf) * jp_pf + (kk - kmin_pf) * kp_pf;
            const int idxm1_pf = idx_pf - kp_pf;

            const int idx_temp = (ii - imin_temp) +
               (jj - jmin_temp) * jp_temp + (kk - kmin_temp) * kp_temp;
            const int idxm1_temp = idx_temp - kp_temp;

            double phi = 
               average( ptr_phi[idx_pf], ptr_phi[idxm1_pf] );
            double c = 
               average( ptr_conc[idx_pf], ptr_conc[idxm1_pf] );
            double t = 
               average( ptr_temp[idx_temp], ptr_temp[idxm1_temp] );
            
            const double hphi =
               FORT_INTERP_FUNC(phi, &interpf);

            double heta = 0.0;
            if ( d_with_third_phase ) {
               double eta = 
                  average( ptr_eta[idx_pf], ptr_eta[idxm1_pf] );
               heta =
                  FORT_INTERP_FUNC(eta, &interpf);
            }

            const double invd2fdc2 = computeLocalInvD2fDc2(c, hphi, heta, t);
            // D=D0*(d2f/dc2)^-1
            ptr_diffusion_2[idx_dcoeff]
               = ptr_diffusion0_2[idx_dcoeff]*invd2fdc2;
         }
      }
   }
#endif

}
#endif

//-----------------------------------------------------------------------

void KKSCompositionRHSStrategy::computeFluxOnPatch(hier::Patch& patch,
                                                   const int flux_id)
{
   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
           patch.getPatchGeometry()));
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

   std::shared_ptr<pdat::SideData<double> > conc_pfm_diffusion(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(d_pfm_diffusion_id)));
   assert(conc_pfm_diffusion);
   assert(conc_pfm_diffusion->getDepth() == 1);

   std::shared_ptr<pdat::SideData<double> > conc_phase_coupling_diffusion(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_phase_coupling_diffusion_id)));
   assert(conc_phase_coupling_diffusion);
   assert(conc_phase_coupling_diffusion->getDepth() == 1);

   std::shared_ptr<pdat::SideData<double> > flux(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(flux_id)));
   assert(flux);
   assert(flux->getDepth() == conc->getDepth());

   int three_phase = 0;
   double* ptr_eta = nullptr;
   double* ptr_conc_eta_diff0 = nullptr;
   double* ptr_conc_eta_diff1 = nullptr;
#if (NDIM == 3)
   double* ptr_conc_eta_diff2 = nullptr;
#endif
   if (d_with_third_phase) {
      three_phase = 1;
      std::shared_ptr<pdat::CellData<double> > eta(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(d_eta_scratch_id)));
      assert(eta);
      ptr_eta = eta->getPointer();
   }

   // now compute concentration flux
   for (int ic = 0; ic < conc->getDepth(); ++ic)
      FORT_CONCENTRATION_FLUX(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                              ifirst(2), ilast(2),
#endif
                              dx, conc->getPointer(ic), NGHOSTS,
                              phase->getPointer(), NGHOSTS, ptr_eta, NGHOSTS,
                              conc_pfm_diffusion->getPointer(0),
                              conc_pfm_diffusion->getPointer(1),
#if (NDIM == 3)
                              conc_pfm_diffusion->getPointer(2),
#endif
                              0, conc_phase_coupling_diffusion->getPointer(0),
                              conc_phase_coupling_diffusion->getPointer(1),
#if (NDIM == 3)
                              conc_phase_coupling_diffusion->getPointer(2),
#endif
                              0, ptr_conc_eta_diff0, ptr_conc_eta_diff1,
#if (NDIM == 3)
                              ptr_conc_eta_diff2,
#endif
                              0, flux->getPointer(0, ic),
                              flux->getPointer(1, ic),
#if (NDIM == 3)
                              flux->getPointer(2, ic),
#endif
                              flux->getGhostCellWidth()[0], three_phase);
}
