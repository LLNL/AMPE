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
#include "KKSCompositionRHSStrategy.h"
#include "QuatFort.h"
#include "ConcFort.h"
#include "FuncFort.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <cassert>

const double gas_constant_R_JpKpmol = GASCONSTANT_R_JPKPMOL;

KKSCompositionRHSStrategy::KKSCompositionRHSStrategy(
    const int conc_scratch_id, const int phase_scratch_id,
    const int conc_pfm_diffusion_id, const int conc_phase_coupling_diffusion_id,
    const int temperature_scratch_id, const int conc_l_scratch_id,
    const int conc_a_scratch_id, const double D_liquid, const double D_solid_A,
    const double Q0_liquid, const double Q0_solid_A,
    const Thermo4PFM::EnergyInterpolationType phase_interp_func_type,
    const std::string& avg_func_type)
    : CompositionRHSStrategy(avg_func_type),
      d_conc_scratch_id(conc_scratch_id),
      d_phase_scratch_id(phase_scratch_id)
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

   d_phase_interp_func_type = phase_interp_func_type;
   d_avg_func_type = avg_func_type;

   d_conc_l_scratch_id = conc_l_scratch_id;
   d_conc_a_scratch_id = conc_a_scratch_id;

   d_D_liquid = D_liquid;
   d_D_solid_A = D_solid_A;

   d_Q0_liquid = Q0_liquid;
   d_Q0_solid_A = Q0_solid_A;

   d_pfm_diffusion_id = conc_pfm_diffusion_id;
   d_conc_phase_coupling_diffusion_id = conc_phase_coupling_diffusion_id;

   d_temperature_scratch_id = temperature_scratch_id;

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

   // set diffusion coefficients for concentration
   setPFMDiffCoeffForConcentration(hierarchy, d_temperature_scratch_id,
                                   d_phase_scratch_id, d_pfm_diffusion_id);

   // set diffusion coefficients for phase
   setDiffCoeffForGradPhi(hierarchy, d_temperature_scratch_id,
                          d_conc_scratch_id, d_phase_scratch_id,
                          d_pfm_diffusion_id,
                          d_conc_phase_coupling_diffusion_id);

   t_set_diffcoeff_timer->stop();
}

//=======================================================================
// compute diffusion for div*D0*grad*c term in concentration equation
void KKSCompositionRHSStrategy::setPFMDiffCoeffForConcentration(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int phase_id,
    const int conc_pfm_diffusion_id)
{
   // tbox::pout<<"KKSCompositionRHSStrategy::setDiffCoeffForConcentration()"<<endl;
   assert(temperature_id >= 0);
   assert(phase_id >= 0);
   assert(conc_pfm_diffusion_id >= 0);

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

         PFMDIFFUSION(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                      ifirst(2), ilast(2),
#endif
                      phi->getPointer(), phi->getGhostCellWidth()[0],
                      pfm_diffusion->getPointer(0),
                      pfm_diffusion->getPointer(1),
#if (NDIM == 3)
                      pfm_diffusion->getPointer(2),
#endif
                      0, temperature->getPointer(),
                      temperature->getGhostCellWidth()[0], d_D_liquid,
                      d_Q0_liquid, d_D_solid_A, d_Q0_solid_A,
                      gas_constant_R_JpKpmol, &interpf,
                      d_avg_func_type.c_str());
      }
   }
}

//=======================================================================
//
// set diffusion coefficients for div*D*grad*phi term in concentration
//
void KKSCompositionRHSStrategy::setDiffCoeffForGradPhi(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int concentration_id, const int phase_id,
    const int conc_pfm_diffusion_id, const int conc_phase_coupling_diffusion_id)
{
   (void)temperature_id;
   (void)concentration_id;

   assert(phase_id >= 0);
   assert(conc_pfm_diffusion_id >= 0);
   assert(conc_phase_coupling_diffusion_id >= 0);

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

         std::shared_ptr<pdat::SideData<double> > sd_pfmd_coeff(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(conc_pfm_diffusion_id)));

         std::shared_ptr<pdat::CellData<double> > cd_phi(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(phase_id)));

         std::shared_ptr<pdat::CellData<double> > cd_c_l(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_conc_l_scratch_id)));

         std::shared_ptr<pdat::CellData<double> > cd_c_a(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_conc_a_scratch_id)));

         setDiffCoeffForPhaseOnPatch(sd_phi_diff_coeff, sd_pfmd_coeff, cd_phi,
                                     cd_c_l, cd_c_a, pbox);
      }
   }
}

//-----------------------------------------------------------------------

void KKSCompositionRHSStrategy::setDiffCoeffForPhaseOnPatch(
    std::shared_ptr<pdat::SideData<double> > sd_phi_diff_coeff,
    std::shared_ptr<pdat::SideData<double> > sd_d0_coeff,
    std::shared_ptr<pdat::CellData<double> > cd_phi,
    std::shared_ptr<pdat::CellData<double> > cd_c_l,
    std::shared_ptr<pdat::CellData<double> > cd_c_a, const hier::Box& pbox)
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
   double* ptr_c_l = cd_c_l->getPointer();
   double* ptr_c_a = cd_c_a->getPointer();

   // Assuming d0_coeff, phi_diff_coeff have same box
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

   // Assuming c_l, c_a all have same box
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

            double c_l = 0.5 * (ptr_c_l[idx_c_i] + ptr_c_l[idxm1_c_i]);
            double c_a = 0.5 * (ptr_c_a[idx_c_i] + ptr_c_a[idxm1_c_i]);

            double hphi_prime = DERIV_INTERP_FUNC(phi, &interpf);
            ptr_phi_diffx[idx_dcoeff] =
                ptr_dx0_coeff[idx_dcoeff] * hphi_prime * (c_l - c_a);
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
            double c_l = 0.5 * (ptr_c_l[idx_c_i] + ptr_c_l[idxm1_c_i]);
            double c_a = 0.5 * (ptr_c_a[idx_c_i] + ptr_c_a[idxm1_c_i]);

            double hphi_prime = DERIV_INTERP_FUNC(phi, &interpf);

            ptr_phi_diffy[idx_dcoeff] =
                ptr_dy0_coeff[idx_dcoeff] * hphi_prime * (c_l - c_a);
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
               double c_l = 0.5 * (ptr_c_l[idx_c_i] + ptr_c_l[idxm1_c_i]);
               double c_a = 0.5 * (ptr_c_a[idx_c_i] + ptr_c_a[idxm1_c_i]);

               double hphi_prime = DERIV_INTERP_FUNC(phi, &interpf);

               ptr_phi_diffz[idx_dcoeff] =
                   ptr_dz0_coeff[idx_dcoeff] * hphi_prime * (c_l - c_a);
            }
         }
      }
   }  // if ( NDIM > 2 )
}

//-----------------------------------------------------------------------

void KKSCompositionRHSStrategy::computeFluxOnPatch(hier::Patch& patch,
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

   // now compute concentration flux
   // D_phi*grad phi+D_conc*grad conc
   for (int ic = 0; ic < conc->getDepth(); ++ic)
      CONCENTRATIONFLUX(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                        ifirst(2), ilast(2),
#endif
                        dx, conc->getPointer(ic), conc->getGhostCellWidth()[0],
                        phase->getPointer(), phase->getGhostCellWidth()[0],
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
                        0, flux->getPointer(0, ic), flux->getPointer(1, ic),
#if (NDIM == 3)
                        flux->getPointer(2, ic),
#endif
                        flux->getGhostCellWidth()[0]);
}
