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

// Ref: Eiken, Boettger, Steinbach, PRE 73, 066122 (2006)
#include "EBSCompositionRHSStrategy.h"
#include "QuatParams.h"
#include "ConcFort.h"
#include "QuatFort.h"
#include "CompositionDiffusionStrategy.h"
#include "CompositionStrategyMobilities.h"

#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/PatchSideDataBasicOps.h"

#include <cassert>


EBSCompositionRHSStrategy::EBSCompositionRHSStrategy(
    const int phase_scratch_id, const int eta_scratch_id,
    const unsigned short ncompositions, const int conc_l_scratch_id,
    const int conc_a_scratch_id, const int conc_b_scratch_id,
    const int temperature_scratch_id, const int diffusion_l_id,
    const int diffusion_a_id, const int diffusion_b_id, const int Mq_id,
    const std::vector<double>& Q_heat_transport,
    const std::vector<int> diffusion_precond_id,
    const std::string& avg_func_type, FreeEnergyStrategy* free_energy_strategy,
    CompositionStrategyMobilities* mobilities_strategy,
    std::shared_ptr<CompositionDiffusionStrategy> diffusion_for_conc_in_phase)
    : CompositionRHSStrategy(avg_func_type),
      d_mobilities_strategy(mobilities_strategy),
      d_diffusion_for_conc_in_phase(diffusion_for_conc_in_phase),
      d_free_energy_strategy(free_energy_strategy)
{
   assert(diffusion_l_id >= 0);
   assert(temperature_scratch_id >= 0);
   assert(d_free_energy_strategy != nullptr);

   if (Mq_id >= 0) {
      assert(Q_heat_transport.size() > 0);
      tbox::plog << "EBSCompositionRHSStrategy with thermal diffusion"
                 << std::endl;
   } else {
      tbox::plog << "EBSCompositionRHSStrategy without thermal diffusion"
                 << std::endl;
   }

   d_ncompositions = ncompositions;

   d_phase_scratch_id = phase_scratch_id;
   d_eta_scratch_id = eta_scratch_id;

   d_conc_l_scratch_id = conc_l_scratch_id;
   d_conc_a_scratch_id = conc_a_scratch_id;
   d_conc_b_scratch_id = conc_b_scratch_id;

   d_diffusion_l_id = diffusion_l_id;
   d_diffusion_a_id = diffusion_a_id;
   d_diffusion_b_id = diffusion_b_id;

   d_Mq_id = Mq_id;
   d_Q_heat_transport = Q_heat_transport;

   d_diffusion_precond_id = diffusion_precond_id;

   d_temperature_scratch_id = temperature_scratch_id;

   d_with_third_phase = (conc_b_scratch_id >= 0);

   d_with_diffusion_for_preconditioner = (!diffusion_precond_id.empty());

   d_with_gradT = (Mq_id >= 0);
}


//-----------------------------------------------------------------------

void EBSCompositionRHSStrategy::setDiffusionCoeff(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time)
{
   (void)time;

   // tbox::pout<<"EBSCompositionRHSStrategy::setDiffusionCoeff"<<endl;
   assert(hierarchy);
   assert(d_free_energy_strategy != nullptr);

   // set coefficient for (grad T)/T term
   if (d_with_gradT)
      setDiffusionCoeffForT(hierarchy, d_conc_l_scratch_id, d_conc_a_scratch_id,
                            d_conc_b_scratch_id, d_temperature_scratch_id);

   // compute actual diffusion, including phase fraction weight
   d_diffusion_for_conc_in_phase->setDiffusion(hierarchy,
                                               d_temperature_scratch_id,
                                               d_phase_scratch_id,
                                               d_eta_scratch_id);

   if (d_with_diffusion_for_preconditioner)
      setDiffusionCoeffForPreconditioner(hierarchy);
}

//-----------------------------------------------------------------------

void EBSCompositionRHSStrategy::computeFluxOnPatch(hier::Patch& patch,
                                                   const int flux_id)
{
   // tbox::pout<<"EBSCompositionRHSStrategy::computeFluxOnPatch"<<endl;
   assert(d_conc_l_scratch_id >= 0);
   assert(d_conc_a_scratch_id >= 0);
   assert(d_diffusion_l_id >= 0);
   assert(d_diffusion_a_id >= 0);

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
           patch.getPatchGeometry()));
   const double* dx = patch_geom->getDx();

   const hier::Box& pbox = patch.getBox();
   const hier::Index& ifirst(pbox.lower());
   const hier::Index& ilast(pbox.upper());

   const unsigned nc2 = d_ncompositions * d_ncompositions;

   std::shared_ptr<pdat::CellData<double> > conc_l(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_l_scratch_id)));
   assert(conc_l);
   assert(conc_l->getDepth() == d_ncompositions);

   std::shared_ptr<pdat::CellData<double> > conc_a(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_a_scratch_id)));
   assert(conc_a);
   assert(conc_a->getDepth() == d_ncompositions);

   std::shared_ptr<pdat::SideData<double> > conc_diffusionl(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(d_diffusion_l_id)));
   assert(conc_diffusionl);
   assert(conc_diffusionl->getDepth() == (nc2));

   // math::PatchSideDataBasicOps<double> ops;
   //{
   // double vmax=ops.max( conc_diffusionl, pbox );
   // double vmin=ops.min(conc_diffusionl, pbox );
   // tbox::pout<<"Min-Max. DL on patch="<<vmin<<","<<vmax<<endl;
   //}

   std::shared_ptr<pdat::SideData<double> > conc_diffusiona(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(d_diffusion_a_id)));
   assert(conc_diffusiona);
   assert(conc_diffusiona->getDepth() == (int)nc2);

   //{
   // double vmax=ops.max( conc_diffusiona, pbox );
   // double vmin=ops.min(conc_diffusiona, pbox );
   // tbox::pout<<"Min-Max. DS on patch="<<vmin<<","<<vmax<<endl;
   //}

   std::shared_ptr<pdat::SideData<double> > flux(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(flux_id)));
   assert(flux);
   assert(flux->getDepth() == d_ncompositions);

   flux->fillAll(0.);

   // now add components of concentration flux,
   // one phase at a time
   ADD_CONCENTRATIONFLUX_EBS(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                             ifirst(2), ilast(2),
#endif
                             dx, conc_l->getPointer(), NGHOSTS, d_ncompositions,
                             conc_diffusionl->getPointer(0),
                             conc_diffusionl->getPointer(1),
#if (NDIM == 3)
                             conc_diffusionl->getPointer(2),
#endif
                             0, flux->getPointer(0), flux->getPointer(1),
#if (NDIM == 3)
                             flux->getPointer(2),
#endif
                             flux->getGhostCellWidth()[0]);

   ADD_CONCENTRATIONFLUX_EBS(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                             ifirst(2), ilast(2),
#endif
                             dx, conc_a->getPointer(), NGHOSTS, d_ncompositions,
                             conc_diffusiona->getPointer(0),
                             conc_diffusiona->getPointer(1),
#if (NDIM == 3)
                             conc_diffusiona->getPointer(2),
#endif
                             0, flux->getPointer(0), flux->getPointer(1),
#if (NDIM == 3)
                             flux->getPointer(2),
#endif
                             flux->getGhostCellWidth()[0]);

   if (d_with_third_phase) {
      assert(d_conc_b_scratch_id >= 0);
      std::shared_ptr<pdat::CellData<double> > conc_b(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(d_conc_b_scratch_id)));
      assert(conc_b);

      std::shared_ptr<pdat::SideData<double> > conc_diffusionb(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch.getPatchData(d_diffusion_b_id)));
      assert(conc_diffusionb);

      assert(conc_diffusionb->getPointer(0) != nullptr);
      assert(conc_diffusionb->getPointer(1) != nullptr);
#if (NDIM == 3)
      assert(conc_diffusionb->getPointer(2) != nullptr);
#endif

      ADD_CONCENTRATIONFLUX_EBS(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                                ifirst(2), ilast(2),
#endif
                                dx, conc_b->getPointer(), NGHOSTS,
                                d_ncompositions, conc_diffusionb->getPointer(0),
                                conc_diffusionb->getPointer(1),
#if (NDIM == 3)
                                conc_diffusionb->getPointer(2),
#endif
                                0, flux->getPointer(0), flux->getPointer(1),
#if (NDIM == 3)
                                flux->getPointer(2),
#endif
                                flux->getGhostCellWidth()[0]);
   }
}

//-----------------------------------------------------------------------

void EBSCompositionRHSStrategy::setDiffusionCoeffForPreconditioner(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   // tbox::pout<<"EBSCompositionRHSStrategy::setDiffusionCoeffForPreconditioner"<<endl;

   assert(d_diffusion_l_id >= 0);
   assert(d_diffusion_a_id >= 0);
   if (d_with_third_phase) {
      assert(d_diffusion_b_id >= 0);
   }
   assert(!d_diffusion_precond_id.empty());

   const int maxl = hierarchy->getNumberOfLevels();

   for (int amr_level = 0; amr_level < maxl; amr_level++) {
      std::shared_ptr<hier::PatchLevel> level =
          hierarchy->getPatchLevel(amr_level);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           p++) {
         std::shared_ptr<hier::Patch> patch = *p;

         std::shared_ptr<pdat::SideData<double> > dl(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(d_diffusion_l_id)));
         std::shared_ptr<pdat::SideData<double> > da(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(d_diffusion_a_id)));
         std::shared_ptr<pdat::SideData<double> > db;

         if (d_with_third_phase) {
            db = SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                patch->getPatchData(d_diffusion_b_id));
         }

         assert((int)(d_diffusion_precond_id.size() *
                      d_diffusion_precond_id.size()) == dl->getDepth());

         for (unsigned int ic = 0; ic < d_diffusion_precond_id.size(); ic++) {

            const int depth_in_Dmat = (ic + 1) * (ic + 1) - 1;
            std::shared_ptr<pdat::SideData<double> > diffusion(
                SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                    patch->getPatchData(d_diffusion_precond_id[ic])));
            TBOX_ASSERT(diffusion);

            setDiffusionCoeffForPreconditionerOnPatch(dl, da, db, depth_in_Dmat,
                                                      diffusion,
                                                      patch->getBox(), 0);
         }
      }
   }
}

//=======================================================================
// phase weighted average of phase diffusions
//
void EBSCompositionRHSStrategy::setDiffusionCoeffForPreconditionerOnPatch(
    std::shared_ptr<pdat::SideData<double> > sd_d_l,
    std::shared_ptr<pdat::SideData<double> > sd_d_a,
    std::shared_ptr<pdat::SideData<double> > sd_d_b,
    const int depth_in_Dmatrix,
    std::shared_ptr<pdat::SideData<double> > sd_d_coeff,
    const hier::Box& pbox, const int depth)
{
   assert(depth_in_Dmatrix < sd_d_l->getDepth());

   double* ptr_dx_coeff = sd_d_coeff->getPointer(0, depth);
   double* ptr_dy_coeff = sd_d_coeff->getPointer(1, depth);
   double* ptr_dz_coeff = nullptr;
   if (NDIM > 2) {
      ptr_dz_coeff = sd_d_coeff->getPointer(2, depth);
   }

   double* ptr_dx_l = sd_d_l->getPointer(0, depth_in_Dmatrix);
   double* ptr_dy_l = sd_d_l->getPointer(1, depth_in_Dmatrix);
   double* ptr_dz_l = nullptr;
   if (NDIM > 2) {
      ptr_dz_l = sd_d_l->getPointer(2, depth_in_Dmatrix);
   }

   double* ptr_dx_a = sd_d_a->getPointer(0, depth_in_Dmatrix);
   double* ptr_dy_a = sd_d_a->getPointer(1, depth_in_Dmatrix);
   double* ptr_dz_a = nullptr;
   if (NDIM > 2) {
      ptr_dz_a = sd_d_a->getPointer(2, depth_in_Dmatrix);
   }

   double* ptr_dx_b = nullptr;
   double* ptr_dy_b = nullptr;
   double* ptr_dz_b = nullptr;
   if (d_with_third_phase) {
      ptr_dx_b = sd_d_b->getPointer(0, depth_in_Dmatrix);
      ptr_dy_b = sd_d_b->getPointer(1, depth_in_Dmatrix);
      if (NDIM > 2) ptr_dz_b = sd_d_b->getPointer(2, depth_in_Dmatrix);
   }

   std::vector<double*> ptr_c_b;
   ptr_c_b.resize(d_ncompositions);

   // Assuming all sd_d_coeff_* have same box
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

   // X-side
   for (int kk = kmin; kk <= kmax; kk++) {
      for (int jj = jmin; jj <= jmax; jj++) {
         for (int ii = imin; ii <= imax + 1; ii++) {

            const int idx_dcoeff =
                (ii - imin_dcoeff) + (jj - jmin_dcoeff) * (jp_dcoeff + 1) +
                (kk - kmin_dcoeff) * (kp_dcoeff + dcoeff_gbox.numberCells(1));
            TBOX_ASSERT(idx_dcoeff < sd_d_coeff->getArrayData(0).getOffset());

            ptr_dx_coeff[idx_dcoeff] =
                ptr_dx_l[idx_dcoeff] + ptr_dx_a[idx_dcoeff];

            if (d_with_third_phase) {
               ptr_dx_coeff[idx_dcoeff] += ptr_dx_b[idx_dcoeff];
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
            TBOX_ASSERT(idx_dcoeff < sd_d_coeff->getArrayData(1).getOffset());

            ptr_dy_coeff[idx_dcoeff] =
                ptr_dy_l[idx_dcoeff] + ptr_dy_a[idx_dcoeff];

            if (d_with_third_phase) {
               ptr_dy_coeff[idx_dcoeff] += ptr_dy_b[idx_dcoeff];
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
               TBOX_ASSERT(idx_dcoeff <
                           sd_d_coeff->getArrayData(2).getOffset());

               ptr_dz_coeff[idx_dcoeff] =
                   ptr_dz_l[idx_dcoeff] + ptr_dz_a[idx_dcoeff];

               if (d_with_third_phase) {
                  ptr_dz_coeff[idx_dcoeff] += ptr_dz_b[idx_dcoeff];
               }
            }
         }
      }
   }  // if ( NDIM > 2 )
}

//-----------------------------------------------------------------------

void EBSCompositionRHSStrategy::addFluxFromAntitrappingonPatch(
    hier::Patch& patch, const int phase_scratch_id, const int dphidt_id,
    const double alpha, const int flux_id)
{
   assert(dphidt_id >= 0);

   // tbox::plog<<"EBSCompositionRHSStrategy::addFluxFromGradTonPatch()"<<endl;

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
           patch.getPatchGeometry()));
   const double* dx = patch_geom->getDx();

   const hier::Box& pbox = patch.getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(phase_scratch_id)));
   assert(phase);
   std::shared_ptr<pdat::CellData<double> > cl(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_l_scratch_id)));
   assert(cl);
   std::shared_ptr<pdat::CellData<double> > ca(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_conc_a_scratch_id)));
   assert(ca);
   std::shared_ptr<pdat::CellData<double> > dphidt(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(dphidt_id)));
   assert(dphidt);
   std::shared_ptr<pdat::SideData<double> > flux(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(flux_id)));
   assert(flux);
   assert(flux->getDepth() == d_ncompositions);

   // needs ghost values to define fluxes through cell boundaries
   assert(cl->getGhostCellWidth()[0] > 0);
   assert(dphidt->getGhostCellWidth()[0] > 0);

   // now compute concentration flux
   ADDCONCENTRATIONFLUXFROMANTITRAPPING(
       ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
       ifirst(2), ilast(2),
#endif
       dx, phase->getPointer(), NGHOSTS, cl->getPointer(), ca->getPointer(),
       cl->getGhostCellWidth()[0], d_ncompositions, dphidt->getPointer(),
       dphidt->getGhostCellWidth()[0], alpha, flux->getPointer(0),
       flux->getPointer(1),
#if (NDIM == 3)
       flux->getPointer(2),
#endif
       flux->getGhostCellWidth()[0]);
}

//-----------------------------------------------------------------------

void EBSCompositionRHSStrategy::setDiffusionCoeffForT(
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
         if (d_with_third_phase) {
            cb = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(concentration_b_id));
            eta = SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                patch->getPatchData(d_eta_scratch_id));
         }

         setDiffusionCoeffForTOnPatch(cl, ca, cb, temp, phi, eta, mq,
                                      patch->getBox());
      }
   }
}

//-----------------------------------------------------------------------

void EBSCompositionRHSStrategy::setDiffusionCoeffForTOnPatch(
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
   if (d_with_third_phase) {
      for (unsigned short ic = 0; ic < d_ncompositions; ic++)
         ptr_c_b[ic] = cd_c_b->getPointer(ic);
   }

   double* ptr_temp = cd_temp->getPointer();

   double* ptr_phi = cd_phi->getPointer();
   double* ptr_eta = nullptr;
   if (d_with_third_phase) {
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
            if (d_with_third_phase) {
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

            if (d_with_third_phase) {
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
            if (d_with_third_phase) {
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

            if (d_with_third_phase) {
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
               if (d_with_third_phase) {
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

               if (d_with_third_phase) {
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

void EBSCompositionRHSStrategy::addFluxFromGradTonPatch(
    hier::Patch& patch, const int temperature_id, const int flux_id)
{
   assert(d_Mq_id >= 0);

   // tbox::plog<<"EBSCompositionRHSStrategy::addFluxFromGradTonPatch()"<<endl;

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
           patch.getPatchGeometry()));
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
                                 dx, temperature->getPointer(), NGHOSTS,
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
