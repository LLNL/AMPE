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
#include "EBSCompositionRHSStrategyWithGradT.h"

#include "ConcFort.h"
#include "QuatFort.h"
#include "FuncFort.h"
#include "CompositionDiffusionStrategy.h"
#include "CompositionStrategyMobilities.h"

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

#include <cassert>

EBSCompositionRHSStrategyWithGradT::EBSCompositionRHSStrategyWithGradT(
    const int phase_scratch_id, const int eta_scratch_id,
    const unsigned short ncompositions, const int conc_l_scratch_id,
    const int conc_a_scratch_id, const int conc_b_scratch_id,
    const int temperature_scratch_id, const int diffusion_l_id,
    const int diffusion_a_id, const int diffusion_b_id, const int Mq_id,
    const std::vector<double>& Q_heat_transport,
    const std::vector<int> diffusion_precond_id,
    const std::string& avg_func_type,
    std::shared_ptr<FreeEnergyStrategy> free_energy_strategy,
    CompositionStrategyMobilities* mobilities_strategy,
    std::shared_ptr<CompositionDiffusionStrategy> diffusion_for_conc_in_phase)
    : EBSCompositionRHSStrategy(phase_scratch_id, ncompositions,
                                conc_l_scratch_id, conc_a_scratch_id,
                                conc_b_scratch_id, temperature_scratch_id,
                                diffusion_l_id, diffusion_a_id, diffusion_b_id,
                                diffusion_precond_id, avg_func_type,
                                free_energy_strategy,
                                diffusion_for_conc_in_phase, "regular")
{
   assert(Mq_id >= 0);
   assert(Q_heat_transport.size() > 0);
   tbox::plog << "EBSCompositionRHSStrategy with thermal diffusion"
              << std::endl;
   d_Mq_id = Mq_id;
   d_Q_heat_transport = Q_heat_transport;
}

void EBSCompositionRHSStrategyWithGradT::setDiffusionCoeff(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time)
{
   // tbox::pout<<"EBSCompositionRHSStrategy::setDiffusionCoeff"<<endl;
   assert(hierarchy);

   // set coefficient for (grad T)/T term
   setDiffusionCoeffForT(hierarchy, d_conc_l_scratch_id, d_conc_a_scratch_id,
                         d_conc_b_scratch_id, d_temperature_scratch_id);


   this->setDiffusionCoeff(hierarchy, time);
}

void EBSCompositionRHSStrategyWithGradT::setDiffusionCoeffForT(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int concentration_l_id, const int concentration_a_id,
    const int concentration_b_id, const int temperature_id)
{
   // tbox::pout<<"EBSCompositionRHSStrategy::setDiffusionCoeffForT"<<endl;
   assert(temperature_id >= 0);
   assert(concentration_l_id >= 0);
   assert(concentration_a_id >= 0);
   assert(d_phase_scratch_id >= 0);
   assert(d_Mq_id >= 0);

   const int maxl = hierarchy->getNumberOfLevels();

   for (int amr_level = 0; amr_level < maxl; amr_level++) {
      std::shared_ptr<hier::PatchLevel> level =
          hierarchy->getPatchLevel(amr_level);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           p++) {
         std::shared_ptr<hier::Patch> patch = *p;

         std::shared_ptr<pdat::CellData<double> > temp(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(temperature_id)));

         std::shared_ptr<pdat::CellData<double> > cl(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(concentration_l_id)));
         std::shared_ptr<pdat::CellData<double> > ca(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(concentration_a_id)));
         std::shared_ptr<pdat::CellData<double> > cb;

         std::shared_ptr<pdat::CellData<double> > phi(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_phase_scratch_id)));

         // data array where result of this operation is stored
         std::shared_ptr<pdat::SideData<double> > mq(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(d_Mq_id)));
         assert(mq);

         std::shared_ptr<pdat::CellData<double> > eta;
         if (concentration_b_id >= 0) {
            cb =
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(concentration_b_id));
         }
         if (d_with_three_phases) {
            eta =
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_eta_scratch_id));
         }

         setDiffusionCoeffForTOnPatch(cl, ca, cb, temp, phi, eta, mq,
                                      patch->getBox());
      }
   }
}

//-----------------------------------------------------------------------

void EBSCompositionRHSStrategyWithGradT::setDiffusionCoeffForTOnPatch(
    std::shared_ptr<pdat::CellData<double> > cd_c_l,
    std::shared_ptr<pdat::CellData<double> > cd_c_a,
    std::shared_ptr<pdat::CellData<double> > cd_c_b,
    std::shared_ptr<pdat::CellData<double> > cd_temp,
    std::shared_ptr<pdat::CellData<double> > cd_phi,
    std::shared_ptr<pdat::CellData<double> > cd_eta,
    std::shared_ptr<pdat::SideData<double> > sd_mq,  // output
    const hier::Box& pbox)
{
   assert(d_ncompositions == 1);

   std::vector<double*> ptr_c_l(d_ncompositions);
   std::vector<double*> ptr_c_a(d_ncompositions);
   std::vector<double*> ptr_c_b(d_ncompositions);
   for (unsigned short ic = 0; ic < d_ncompositions; ic++) {
      ptr_c_l[ic] = cd_c_l->getPointer(ic);
      ptr_c_a[ic] = cd_c_a->getPointer(ic);
   }
   if (d_with_three_phases) {
      for (unsigned short ic = 0; ic < d_ncompositions; ic++)
         ptr_c_b[ic] = cd_c_b->getPointer(ic);
   }

   double* ptr_temp = cd_temp->getPointer();

   double* ptr_phi = cd_phi->getPointer();
   double* ptr_eta = nullptr;
   if (d_with_three_phases) {
      ptr_eta = cd_eta->getPointer();
   }

   std::vector<double*> ptr_dx_mq;
   ptr_dx_mq.resize(d_ncompositions * d_ncompositions);
   std::vector<double*> ptr_dy_mq;
   ptr_dy_mq.resize(d_ncompositions * d_ncompositions);
   std::vector<double*> ptr_dz_mq;
   ptr_dz_mq.resize(d_ncompositions * d_ncompositions, nullptr);
   for (unsigned short ic = 0; ic < d_ncompositions; ic++)
      for (unsigned short jc = 0; jc < d_ncompositions; jc++) {
         const unsigned ijc = ic + jc * d_ncompositions;
         ptr_dx_mq[ijc] = sd_mq->getPointer(0, ijc);
         ptr_dy_mq[ijc] = sd_mq->getPointer(1, ijc);
         if (NDIM > 2) {
            ptr_dz_mq[ijc] = sd_mq->getPointer(2, ijc);
         }
      }
   // Assuming all sd_d_* have same box
   const hier::Box& mq_gbox = sd_mq->getGhostBox();
   int imin_mq = mq_gbox.lower(0);
   int jmin_mq = mq_gbox.lower(1);
   int jp_mq = mq_gbox.numberCells(0);
   int kmin_mq = 0;
   int kp_mq = 0;
#if (NDIM == 3)
   kmin_mq = mq_gbox.lower(2);
   kp_mq = jp_mq * mq_gbox.numberCells(1);
#endif

   // assumes all c have same number of ghosts
   const hier::Box& c_gbox = cd_c_l->getGhostBox();
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

   std::vector<double> c_l(d_ncompositions);
   std::vector<double> c_a(d_ncompositions);
   std::vector<double> c_b(d_ncompositions);

   std::vector<double> mobmat0(d_ncompositions * d_ncompositions);
   std::vector<double> mobmat1(d_ncompositions * d_ncompositions);

   char interp_type = 'l';

   // X-side
   for (int kk = kmin; kk <= kmax; kk++) {
      for (int jj = jmin; jj <= jmax; jj++) {
         for (int ii = imin; ii <= imax + 1; ii++) {

            const int idx_mq =
                (ii - imin_mq) + (jj - jmin_mq) * (jp_mq + 1) +
                (kk - kmin_mq) * (kp_mq + mq_gbox.numberCells(1));

            const int idx_pf = (ii - imin_pf) + (jj - jmin_pf) * jp_pf +
                               (kk - kmin_pf) * kp_pf;
            const int idxm1_pf = idx_pf - 1;

            const int idx_c =
                (ii - imin_c) + (jj - jmin_c) * jp_c + (kk - kmin_c) * kp_c;
            const int idxm1_c = idx_c - 1;

            const int idx_temp = (ii - imin_temp) + (jj - jmin_temp) * jp_temp +
                                 (kk - kmin_temp) * kp_temp;
            const int idxm1_temp = idx_temp - 1;

            double temp = 0.5 * (ptr_temp[idx_temp] + ptr_temp[idxm1_temp]);

            double phi = average(ptr_phi[idx_pf], ptr_phi[idxm1_pf]);
            double hphi = INTERP_FUNC(phi, &interp_type);
            double heta = 0.;
            if (d_with_three_phases) {
               double eta = average(ptr_eta[idx_pf], ptr_eta[idxm1_pf]);
               heta = INTERP_FUNC(eta, &interp_type);
            }

            for (unsigned short ic = 0; ic < d_ncompositions; ic++) {
               c_l[ic] = 0.5 * (ptr_c_l[ic][idx_c] + ptr_c_l[ic][idxm1_c]);
               c_a[ic] = 0.5 * (ptr_c_a[ic][idx_c] + ptr_c_a[ic][idxm1_c]);
            }

            double mobmat0 =
                d_mobilities_strategy->computeDiffusionMobilityBinaryPhaseL(
                    c_l[0], temp);
            double mobmat1 =
                d_mobilities_strategy->computeDiffusionMobilityBinaryPhaseL(
                    1. - c_l[0], temp);

            ptr_dx_mq[0][idx_mq] = (1. - hphi) * c_l[0] * (1. - c_l[0]) *
                                   (mobmat0 * d_Q_heat_transport[0] -
                                    mobmat1 * d_Q_heat_transport[1]);

            mobmat0 = d_mobilities_strategy
                          ->computeDiffusionMobilityBinaryPhaseA(c_a[0], temp);
            mobmat1 =
                d_mobilities_strategy->computeDiffusionMobilityBinaryPhaseA(
                    1. - c_a[0], temp);

            ptr_dx_mq[0][idx_mq] += (1. - heta) * hphi * c_a[0] *
                                    (1. - c_a[0]) *
                                    (mobmat0 * d_Q_heat_transport[0] -
                                     mobmat1 * d_Q_heat_transport[1]);

            if (d_with_three_phases) {
               for (unsigned short ic = 0; ic < d_ncompositions; ic++) {
                  c_b[ic] = 0.5 * (ptr_c_b[ic][idx_c] + ptr_c_b[ic][idxm1_c]);
               }
               mobmat0 =
                   d_mobilities_strategy->computeDiffusionMobilityBinaryPhaseB(
                       c_b[0], temp);
               mobmat1 =
                   d_mobilities_strategy->computeDiffusionMobilityBinaryPhaseB(
                       1. - c_b[0], temp);

               ptr_dx_mq[0][idx_mq] += heta * hphi * c_b[0] * (1. - c_b[0]) *
                                       (mobmat0 * d_Q_heat_transport[0] -
                                        mobmat1 * d_Q_heat_transport[1]);
            }
         }
      }
   }

   // Y-side
   for (int kk = kmin; kk <= kmax; kk++) {
      for (int jj = jmin; jj <= jmax + 1; jj++) {
         for (int ii = imin; ii <= imax; ii++) {

            const int idx_mq = (ii - imin_mq) + (jj - jmin_mq) * jp_mq +
                               (kk - kmin_mq) * (kp_mq + jp_mq);

            const int idx_pf = (ii - imin_pf) + (jj - jmin_pf) * jp_pf +
                               (kk - kmin_pf) * kp_pf;
            const int idxm1_pf = idx_pf - jp_pf;

            const int idx_c =
                (ii - imin_c) + (jj - jmin_c) * jp_c + (kk - kmin_c) * kp_c;
            const int idxm1_c = idx_c - jp_c;

            const int idx_temp = (ii - imin_temp) + (jj - jmin_temp) * jp_temp +
                                 (kk - kmin_temp) * kp_temp;
            const int idxm1_temp = idx_temp - jp_temp;

            for (unsigned short ic = 0; ic < d_ncompositions; ic++) {
               c_l[ic] = 0.5 * (ptr_c_l[ic][idx_c] + ptr_c_l[ic][idxm1_c]);
               c_a[ic] = 0.5 * (ptr_c_a[ic][idx_c] + ptr_c_a[ic][idxm1_c]);
            }

            double temp = 0.5 * (ptr_temp[idx_temp] + ptr_temp[idxm1_temp]);

            double phi = average(ptr_phi[idx_pf], ptr_phi[idxm1_pf]);
            double hphi = INTERP_FUNC(phi, &interp_type);

            double heta = 0.0;
            if (d_with_three_phases) {
               double eta = average(ptr_eta[idx_pf], ptr_eta[idxm1_pf]);
               heta = INTERP_FUNC(eta, &interp_type);
            }

            double mobmat0 =
                d_mobilities_strategy->computeDiffusionMobilityBinaryPhaseL(
                    c_l[0], temp);
            double mobmat1 =
                d_mobilities_strategy->computeDiffusionMobilityBinaryPhaseL(
                    1. - c_l[0], temp);

            ptr_dy_mq[0][idx_mq] = (1. - hphi) * c_l[0] * (1. - c_l[0]) *
                                   (mobmat0 * d_Q_heat_transport[0] -
                                    mobmat1 * d_Q_heat_transport[1]);

            mobmat0 = d_mobilities_strategy
                          ->computeDiffusionMobilityBinaryPhaseA(c_a[0], temp);
            mobmat1 =
                d_mobilities_strategy->computeDiffusionMobilityBinaryPhaseA(
                    1. - c_a[0], temp);

            ptr_dy_mq[0][idx_mq] += (1. - heta) * hphi * c_a[0] *
                                    (1. - c_a[0]) *
                                    (mobmat0 * d_Q_heat_transport[0] -
                                     mobmat1 * d_Q_heat_transport[1]);

            if (d_with_three_phases) {
               for (unsigned short ic = 0; ic < d_ncompositions; ic++)
                  c_b[ic] = 0.5 * (ptr_c_b[ic][idx_c] + ptr_c_b[ic][idxm1_c]);

               mobmat0 =
                   d_mobilities_strategy->computeDiffusionMobilityBinaryPhaseB(
                       c_b[0], temp);
               mobmat1 =
                   d_mobilities_strategy->computeDiffusionMobilityBinaryPhaseB(
                       1. - c_b[0], temp);

               ptr_dy_mq[0][idx_mq] += heta * hphi * c_b[0] * (1. - c_b[0]) *
                                       (mobmat0 * d_Q_heat_transport[0] -
                                        mobmat1 * d_Q_heat_transport[1]);
            }
         }
      }
   }

   if (NDIM > 2) {
      // Z-side
      for (int kk = kmin; kk <= kmax + 1; kk++) {
         for (int jj = jmin; jj <= jmax; jj++) {
            for (int ii = imin; ii <= imax; ii++) {

               const int idx_mq = (ii - imin_mq) + (jj - jmin_mq) * jp_mq +
                                  (kk - kmin_mq) * kp_mq;

               const int idx_pf = (ii - imin_pf) + (jj - jmin_pf) * jp_pf +
                                  (kk - kmin_pf) * kp_pf;
               const int idxm1_pf = idx_pf - kp_pf;

               const int idx_c =
                   (ii - imin_c) + (jj - jmin_c) * jp_c + (kk - kmin_c) * kp_c;
               const int idxm1_c = idx_c - kp_c;

               const int idx_temp = (ii - imin_temp) +
                                    (jj - jmin_temp) * jp_temp +
                                    (kk - kmin_temp) * kp_temp;
               const int idxm1_temp = idx_temp - kp_temp;

               for (unsigned short ic = 0; ic < d_ncompositions; ic++) {
                  c_l[ic] = 0.5 * (ptr_c_l[ic][idx_c] + ptr_c_l[ic][idxm1_c]);
                  c_a[ic] = 0.5 * (ptr_c_a[ic][idx_c] + ptr_c_a[ic][idxm1_c]);
               }

               double temp = 0.5 * (ptr_temp[idx_temp] + ptr_temp[idxm1_temp]);

               double phi = average(ptr_phi[idx_pf], ptr_phi[idxm1_pf]);
               double hphi = INTERP_FUNC(phi, &interp_type);

               double heta = 0.0;
               if (d_with_three_phases) {
                  double eta = average(ptr_eta[idx_pf], ptr_eta[idxm1_pf]);
                  heta = INTERP_FUNC(eta, &interp_type);
               }

               double mobmat0 =
                   d_mobilities_strategy->computeDiffusionMobilityBinaryPhaseL(
                       c_l[0], temp);
               double mobmat1 =
                   d_mobilities_strategy->computeDiffusionMobilityBinaryPhaseL(
                       1. - c_l[0], temp);

               ptr_dz_mq[0][idx_mq] = (1. - hphi) * c_l[0] * (1. - c_l[0]) *
                                      (mobmat0 * d_Q_heat_transport[0] -
                                       mobmat1 * d_Q_heat_transport[1]);

               mobmat0 =
                   d_mobilities_strategy->computeDiffusionMobilityBinaryPhaseA(
                       c_a[0], temp);
               mobmat1 =
                   d_mobilities_strategy->computeDiffusionMobilityBinaryPhaseA(
                       1. - c_a[0], temp);

               ptr_dz_mq[0][idx_mq] += (1. - heta) * hphi * c_a[0] *
                                       (1. - c_a[0]) *
                                       (mobmat0 * d_Q_heat_transport[0] -
                                        mobmat1 * d_Q_heat_transport[1]);

               if (d_with_three_phases) {
                  for (unsigned short ic = 0; ic < d_ncompositions; ic++)
                     c_b[ic] =
                         0.5 * (ptr_c_b[ic][idx_c] + ptr_c_b[ic][idxm1_c]);

                  mobmat0 =
                      d_mobilities_strategy
                          ->computeDiffusionMobilityBinaryPhaseB(c_b[0], temp);
                  mobmat1 =
                      d_mobilities_strategy
                          ->computeDiffusionMobilityBinaryPhaseB(1. - c_b[0],
                                                                 temp);

                  ptr_dz_mq[0][idx_mq] += heta * hphi * c_b[0] * (1. - c_b[0]) *
                                          (mobmat0 * d_Q_heat_transport[0] -
                                           mobmat1 * d_Q_heat_transport[1]);
               }
            }
         }
      }
   }  // if ( NDIM > 2 )
}

//-----------------------------------------------------------------------

void EBSCompositionRHSStrategyWithGradT::addFluxFromGradTonPatch(
    hier::Patch& patch, const int temperature_id, const int flux_id)
{
   assert(d_Mq_id >= 0);

   // tbox::plog<<"EBSCompositionRHSStrategy::addFluxFromGradTonPatch()"<<endl;

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));
   const double* dx = patch_geom->getDx();

   const hier::Box& pbox = patch.getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   std::shared_ptr<pdat::CellData<double> > temperature(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_id)));
   assert(temperature);
   std::shared_ptr<pdat::SideData<double> > flux(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(flux_id)));
   assert(flux);
   std::shared_ptr<pdat::SideData<double> > mq(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(d_Mq_id)));
   assert(mq);

   // now compute concentration flux
   ADDCONCENTRATIONFLUXFROMGRADT(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                                 ifirst(2), ilast(2),
#endif
                                 dx, temperature->getPointer(),
                                 temperature->getGhostCellWidth()[0],
                                 mq->getPointer(0), mq->getPointer(1),
#if (NDIM == 3)
                                 mq->getPointer(2),
#endif
                                 mq->getGhostCellWidth()[0],
                                 flux->getPointer(0), flux->getPointer(1),
#if (NDIM == 3)
                                 flux->getPointer(2),
#endif
                                 flux->getGhostCellWidth()[0], average());
}
