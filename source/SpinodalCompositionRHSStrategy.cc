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
//
#include "SpinodalCompositionRHSStrategy.h"
#include "ConcFort.h"
#include "QuatParams.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"


using namespace SAMRAI;

SpinodalCompositionRHSStrategy::SpinodalCompositionRHSStrategy(
    std::shared_ptr<tbox::Database> input_db, const int conc_scratch_id,
    const int phase_scratch_id, const int eta_scratch_id,
    const unsigned int ncompositions, const int conc_a_scratch_id,
    const int conc_b_scratch_id, const int temperature_scratch_id,
    const int diffusion_id, const double kappa, const int Mq_id,
    const std::vector<double>& Q_heat_transport,
    const std::string& phase_interp_func_type, const std::string& avg_func_type,
    FreeEnergyStrategy* free_energy_strategy)
    : CompositionRHSStrategyWithMobilities(input_db, phase_scratch_id,
                                           eta_scratch_id, ncompositions,
                                           temperature_scratch_id, Mq_id,
                                           Q_heat_transport,
                                           phase_interp_func_type,
                                           avg_func_type, free_energy_strategy)
{
   assert(eta_scratch_id > -1);
   assert(conc_a_scratch_id > -1);
   assert(conc_b_scratch_id > -1);
   assert(diffusion_id > -1);
   assert(ncompositions > 0);

   d_conc_scratch_id = conc_scratch_id;
   d_phase_scratch_id = phase_scratch_id;
   d_eta_scratch_id = phase_scratch_id;

   d_ncompositions = ncompositions;

   d_diffusion_id = diffusion_id;
   d_conc_a_scratch_id = conc_a_scratch_id;
   d_conc_b_scratch_id = conc_b_scratch_id;
   d_temperature_scratch_id = temperature_scratch_id;

   d_kappa = kappa;
}

//-----------------------------------------------------------------------

void SpinodalCompositionRHSStrategy::computeFluxOnPatch(hier::Patch& patch,
                                                        const int flux_id)
{
   // tbox::pout<<"EBSCompositionRHSStrategy::computeFluxOnPatch"<<endl;
   assert(d_conc_scratch_id >= 0);
   assert(d_diffusion_id >= 0);
   assert(d_eta_scratch_id >= 0);

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

   std::shared_ptr<pdat::CellData<double> > conca(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_a_scratch_id)));
   assert(conca);

   std::shared_ptr<pdat::CellData<double> > concb(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_b_scratch_id)));
   assert(concb);

   std::shared_ptr<pdat::CellData<double> > eta(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_eta_scratch_id)));
   assert(eta);

   std::shared_ptr<pdat::SideData<double> > conc_diffusion(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(d_diffusion_id)));
   assert(conc_diffusion);
   assert(conc_diffusion->getDepth() == (1));


   std::shared_ptr<pdat::SideData<double> > flux(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(flux_id)));
   assert(flux);
   assert(flux->getDepth() == 1);

   // now compute concentration flux
   CONCENTRATION_FLUX_SPINODAL(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                               ifirst(2), ilast(2),
#endif
                               dx, conc->getPointer(), NGHOSTS, 1,
                               conca->getPointer(), NGHOSTS,
                               concb->getPointer(), NGHOSTS,
                               conc_diffusion->getPointer(0),
                               conc_diffusion->getPointer(1),
#if (NDIM == 3)
                               conc_diffusion->getPointer(2),
#endif
                               0, eta->getPointer(), NGHOSTS, d_kappa,
                               flux->getPointer(0), flux->getPointer(1),
#if (NDIM == 3)
                               flux->getPointer(2),
#endif
                               flux->getGhostCellWidth()[0]);
}

//-----------------------------------------------------------------------

void SpinodalCompositionRHSStrategy::setDiffusionCoeff(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time)
{
   // tbox::pout<<"EBSCompositionRHSStrategy::setDiffusionCoeff"<<endl;
   assert(hierarchy);
   assert(d_free_energy_strategy != NULL);

   // set coefficient for (grad T)/T term
   // if( d_with_gradT )
   //   setDiffusionCoeffForT(
   //      hierarchy,
   //      d_conc_l_scratch_id,
   //      d_conc_a_scratch_id,
   //      d_conc_b_scratch_id,
   //      d_temperature_scratch_id);

   setDiffusionForConc(hierarchy);
}

//=======================================================================

void SpinodalCompositionRHSStrategy::setDiffusionForConc(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   // tbox::pout<<"EBSCompositionRHSStrategy::setDiffusionForConcInPhase"<<endl;
   assert(d_diffusion_id >= 0);

   const int maxl = hierarchy->getNumberOfLevels();

   for (int amr_level = 0; amr_level < maxl; amr_level++) {
      std::shared_ptr<hier::PatchLevel> level =
          hierarchy->getPatchLevel(amr_level);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           p++) {
         std::shared_ptr<hier::Patch> patch = *p;

         std::shared_ptr<pdat::CellData<double> > conc(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_conc_scratch_id)));

         std::shared_ptr<pdat::CellData<double> > temp(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_temperature_scratch_id)));

         std::shared_ptr<pdat::SideData<double> > diffusion(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(d_diffusion_id)));

         setDiffusionCoeffForConcOnPatch(conc, temp, diffusion,
                                         patch->getBox());
      }
   }
}

//-----------------------------------------------------------------------

void SpinodalCompositionRHSStrategy::setDiffusionCoeffForConcOnPatch(
    std::shared_ptr<pdat::CellData<double> > cd_c,
    std::shared_ptr<pdat::CellData<double> > cd_temp,
    std::shared_ptr<pdat::SideData<double> > sd_d_coeff, const hier::Box& pbox)
{
   std::vector<double*> ptr_c(d_ncompositions);
   for (unsigned int ic = 0; ic < d_ncompositions; ic++) {
      ptr_c[ic] = cd_c->getPointer(ic);
   }

   double* ptr_temp = cd_temp->getPointer();

   std::vector<double*> ptr_dx_coeff(d_ncompositions * d_ncompositions);
   std::vector<double*> ptr_dy_coeff(d_ncompositions * d_ncompositions);
   std::vector<double*> ptr_dz_coeff(d_ncompositions * d_ncompositions, NULL);
   for (unsigned int ic = 0; ic < d_ncompositions; ic++)
      for (unsigned int jc = 0; jc < d_ncompositions; jc++) {
         const unsigned int ijc = ic + jc * d_ncompositions;
         ptr_dx_coeff[ijc] = sd_d_coeff->getPointer(0, ijc);
         ptr_dy_coeff[ijc] = sd_d_coeff->getPointer(1, ijc);
         if (NDIM > 2) {
            ptr_dz_coeff[ijc] = sd_d_coeff->getPointer(2, ijc);
         }
      }

   // Assuming all sd_d_* have same box
   const hier::Box& dcoeff_gbox = sd_d_coeff->getGhostBox();
   int imin_dcoeff = dcoeff_gbox.lower(0);
   int jmin_dcoeff = dcoeff_gbox.lower(1);
   int jp_dcoeff = dcoeff_gbox.numberCells(0);
   int kmin_dcoeff = 0;
   int kp_dcoeff = 0;
#if (NDIM == 3)
   kmin_dcoeff = dcoeff_gbox.lower(2);
   kp_dcoeff = jp_dcoeff * dcoeff_gbox.numberCells(1);
#endif

   // assumes all c have same number of ghosts
   const hier::Box& c_gbox = cd_c->getGhostBox();
   int imin_c = c_gbox.lower(0);
   int jmin_c = c_gbox.lower(1);
   int jp_c = c_gbox.numberCells(0);
   int kmin_c = 0;
   int kp_c = 0;
#if (NDIM == 3)
   kmin_c = c_gbox.lower(2);
   kp_c = jp_c * c_gbox.numberCells(1);
#endif

   const hier::Box& temp_gbox = cd_temp->getGhostBox();
   int imin_temp = temp_gbox.lower(0);
   int jmin_temp = temp_gbox.lower(1);
   int jp_temp = temp_gbox.numberCells(0);
   int kmin_temp = 0;
   int kp_temp = 0;
#if (NDIM == 3)
   kmin_temp = temp_gbox.lower(2);
   kp_temp = jp_temp * temp_gbox.numberCells(1);
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

   std::vector<double> c(d_ncompositions);

   std::vector<double> d2f(d_ncompositions * d_ncompositions);
   std::vector<double> mobmat(d_ncompositions * d_ncompositions);

   // X-side
   for (int kk = kmin; kk <= kmax; kk++) {
      for (int jj = jmin; jj <= jmax; jj++) {
         for (int ii = imin; ii <= imax + 1; ii++) {

            const int idx_dcoeff =
                (ii - imin_dcoeff) + (jj - jmin_dcoeff) * (jp_dcoeff + 1) +
                (kk - kmin_dcoeff) * (kp_dcoeff + dcoeff_gbox.numberCells(1));

            const int idx_c =
                (ii - imin_c) + (jj - jmin_c) * jp_c + (kk - kmin_c) * kp_c;
            const int idxm1_c = idx_c - 1;

            const int idx_temp = (ii - imin_temp) + (jj - jmin_temp) * jp_temp +
                                 (kk - kmin_temp) * kp_temp;
            const int idxm1_temp = idx_temp - 1;

            double temp = 0.5 * (ptr_temp[idx_temp] + ptr_temp[idxm1_temp]);

            for (unsigned int ic = 0; ic < d_ncompositions; ic++) {
               c[ic] = 0.5 * (ptr_c[ic][idx_c] + ptr_c[ic][idxm1_c]);
            }

            d_free_energy_strategy->computeSecondDerivativeEnergyPhaseA(temp, c,
                                                                        d2f,
                                                                        false);
            computeDiffusionMobilityPhaseA(c, temp, mobmat);

            ptr_dx_coeff[0][idx_dcoeff] = mobmat[0] * d2f[0];
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

            const int idx_c =
                (ii - imin_c) + (jj - jmin_c) * jp_c + (kk - kmin_c) * kp_c;
            const int idxm1_c = idx_c - jp_c;

            const int idx_temp = (ii - imin_temp) + (jj - jmin_temp) * jp_temp +
                                 (kk - kmin_temp) * kp_temp;
            const int idxm1_temp = idx_temp - jp_temp;

            for (unsigned int ic = 0; ic < d_ncompositions; ic++) {
               c[ic] = 0.5 * (ptr_c[ic][idx_c] + ptr_c[ic][idxm1_c]);
            }

            double temp = 0.5 * (ptr_temp[idx_temp] + ptr_temp[idxm1_temp]);

            d_free_energy_strategy->computeSecondDerivativeEnergyPhaseA(temp, c,
                                                                        d2f,
                                                                        false);
            computeDiffusionMobilityPhaseA(c, temp, mobmat);

            ptr_dy_coeff[0][idx_dcoeff] = mobmat[0] * d2f[0];
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

               const int idx_c =
                   (ii - imin_c) + (jj - jmin_c) * jp_c + (kk - kmin_c) * kp_c;
               const int idxm1_c = idx_c - kp_c;

               const int idx_temp = (ii - imin_temp) +
                                    (jj - jmin_temp) * jp_temp +
                                    (kk - kmin_temp) * kp_temp;
               const int idxm1_temp = idx_temp - kp_temp;

               for (unsigned int ic = 0; ic < d_ncompositions; ic++) {
                  c[ic] = 0.5 * (ptr_c[ic][idx_c] + ptr_c[ic][idxm1_c]);
               }

               double temp = 0.5 * (ptr_temp[idx_temp] + ptr_temp[idxm1_temp]);

               d_free_energy_strategy->computeSecondDerivativeEnergyPhaseA(
                   temp, c, d2f, false);
               computeDiffusionMobilityPhaseA(c, temp, mobmat);

               ptr_dz_coeff[0][idx_dcoeff] = mobmat[0] * d2f[0];
            }
         }
      }
   }  // if ( NDIM > 2 )
}
