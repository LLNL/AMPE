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

// Ref: Eiken, Boettger, Steinbach, PRE 73, 066122 (2006)
#include "EBSCompositionRHSStrategy.h"
#include "QuatParams.h"
#include "ConcFort.h"
#include "QuatFort.h"
#include "FuncFort.h"
#include "CompositionDiffusionStrategy.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
//#include "SAMRAI/math/PatchSideDataBasicOps.h"

#include <cassert>


EBSCompositionRHSStrategy::EBSCompositionRHSStrategy(
    const int phase_scratch_id, const unsigned short ncompositions,
    const int conc_l_scratch_id, const int conc_a_scratch_id,
    const int conc_b_scratch_id, const int temperature_scratch_id,
    const int diffusion_l_id, const int diffusion_a_id,
    const int diffusion_b_id, const std::vector<int> diffusion_precond_id,
    const std::string& avg_func_type,
    std::shared_ptr<FreeEnergyStrategy> free_energy_strategy,
    std::shared_ptr<CompositionDiffusionStrategy> diffusion_for_conc_in_phase)
    : CompositionRHSStrategy(avg_func_type),
      d_diffusion_for_conc_in_phase(diffusion_for_conc_in_phase),
      d_free_energy_strategy(free_energy_strategy)
{
   assert(diffusion_l_id >= 0);
   assert(temperature_scratch_id >= 0);
   assert(d_free_energy_strategy);

   tbox::plog << "EBSCompositionRHSStrategy without thermal diffusion"
              << std::endl;

   d_ncompositions = ncompositions;

   d_phase_scratch_id = phase_scratch_id;

   d_conc_l_scratch_id = conc_l_scratch_id;
   d_conc_a_scratch_id = conc_a_scratch_id;
   d_conc_b_scratch_id = conc_b_scratch_id;

   d_diffusion_l_id = diffusion_l_id;
   d_diffusion_a_id = diffusion_a_id;
   d_diffusion_b_id = diffusion_b_id;

   d_diffusion_precond_id = diffusion_precond_id;

   d_temperature_scratch_id = temperature_scratch_id;

   d_with_three_phases = (conc_b_scratch_id >= 0);

   d_with_diffusion_for_preconditioner = (!diffusion_precond_id.empty());

   if (d_with_three_phases)
      tbox::plog << "EBSCompositionRHSStrategy with three phases" << std::endl;
}


//-----------------------------------------------------------------------

void EBSCompositionRHSStrategy::setDiffusionCoeff(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time)
{
   (void)time;

   // tbox::pout<<"EBSCompositionRHSStrategy::setDiffusionCoeff"<<endl;
   assert(hierarchy);
   assert(d_free_energy_strategy);

   // compute actual diffusion, including phase fraction weight
   d_diffusion_for_conc_in_phase->setDiffusion(hierarchy,
                                               d_temperature_scratch_id,
                                               d_phase_scratch_id, -1);

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
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));
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

   if (d_conc_b_scratch_id >= 0) {
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
   if (d_with_three_phases) {
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

         if (d_with_three_phases) {
            db =
                SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                    patch->getPatchData(d_diffusion_b_id));
            assert(db);
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
    std::shared_ptr<pdat::SideData<double> > sd_d_b, const int depth_in_Dmatrix,
    std::shared_ptr<pdat::SideData<double> > sd_d_coeff, const hier::Box& pbox,
    const int depth)
{
   assert(sd_d_l);
   assert(sd_d_a);
   assert(sd_d_a->getDirectionVector()[0] != 0);
   assert(sd_d_a->getDirectionVector()[1] != 0);
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
   assert(ptr_dy_l != nullptr);

   double* ptr_dx_a = sd_d_a->getPointer(0, depth_in_Dmatrix);
   double* ptr_dy_a = sd_d_a->getPointer(1, depth_in_Dmatrix);
   double* ptr_dz_a = nullptr;
   if (NDIM > 2) {
      ptr_dz_a = sd_d_a->getPointer(2, depth_in_Dmatrix);
   }
   assert(ptr_dx_a != nullptr);
   assert(ptr_dy_a != nullptr);

   double* ptr_dx_b = nullptr;
   double* ptr_dy_b = nullptr;
   double* ptr_dz_b = nullptr;
   if (sd_d_b) {
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

            if (sd_d_b) {
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

            if (sd_d_b) {
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

               if (sd_d_b) {
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
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));
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

   // now compute concentration flux anti-trapping correction
   if (d_with_three_phases) {
      std::shared_ptr<pdat::CellData<double> > cb(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(d_conc_b_scratch_id)));
      assert(cb);
      ADDCONCENTRATIONFLUXFROMANTITRAPPING3PHASES(
          ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
          ifirst(2), ilast(2),
#endif
          dx, phase->getPointer(), NGHOSTS, cl->getPointer(), ca->getPointer(),
          cb->getPointer(), cl->getGhostCellWidth()[0], dphidt->getPointer(),
          dphidt->getGhostCellWidth()[0], alpha, flux->getPointer(0),
          flux->getPointer(1),
#if (NDIM == 3)
          flux->getPointer(2),
#endif
          flux->getGhostCellWidth()[0]);
   } else
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
