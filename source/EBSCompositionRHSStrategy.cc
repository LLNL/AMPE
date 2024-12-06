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
#include "ConcFort.h"
#include "QuatFort.h"
#include "FuncFort.h"
#include "CompositionDiffusionStrategy.h"
#include "ArrayOperation.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

#include <cassert>

//-----------------------------------------------------------------------
#if (NDIM == 3)
void add_flux_ebs(const int& ifirst0, const int& ilast0, const int& ifirst1,
                  const int& ilast1, const int& ifirst2, const int& ilast2,
                  const double* dx, const double* conc, const int& ngconc,
                  const int& ncomp, const double* diffconc0,
                  const double* diffconc1, const double* diffconc2,
                  const int& ngdiffconc, const double* flux0,
                  const double* flux1, const double* flux2, const int& ngflux,
                  const int* physb)
{
   (void)physb;

   ADD_FLUX(ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2, dx, conc, ngconc,
            ncomp, diffconc0, diffconc1, diffconc2, ngdiffconc, flux0, flux1,
            flux2, ngflux);
}

void add_flux_iso(const int& ifirst0, const int& ilast0, const int& ifirst1,
                  const int& ilast1, const int& ifirst2, const int& ilast2,
                  const double* dx, const double* conc, const int& ngconc,
                  const int& ncomp, const double* diffconc0,
                  const double* diffconc1, const double* diffconc2,
                  const int& ngdiffconc, const double* flux0,
                  const double* flux1, const double* flux2, const int& ngflux,
                  const int* physb)
{
   (void)physb;
   assert(ngdiffconc > 0);
   ADD_FLUX_ISO(ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2, dx, conc,
                ngconc, ncomp, diffconc0, diffconc1, diffconc2, ngdiffconc,
                flux0, flux1, flux2, ngflux);
}

void add_flux_4th(const int& ifirst0, const int& ilast0, const int& ifirst1,
                  const int& ilast1, const int& ifirst2, const int& ilast2,
                  const double* dx, const double* conc, const int& ngconc,
                  const int& ncomp, const double* diffconc0,
                  const double* diffconc1, const double* diffconc2,
                  const int& ngdiffconc, const double* flux0,
                  const double* flux1, const double* flux2, const int& ngflux,
                  const int* physb)
{
   //   assert(ngdiffconc > 0);
   ADD_FLUX_4TH(ifirst0, ilast0, ifirst1, ilast1, ifirst2, ilast2, dx, conc,
                ngconc, ncomp, diffconc0, diffconc1, diffconc2, ngdiffconc,
                flux0, flux1, flux2, ngflux, physb);
}
#else
void add_flux_ebs(const int& ifirst0, const int& ilast0, const int& ifirst1,
                  const int& ilast1, const double* dx, const double* conc,
                  const int& ngconc, const int& ncomp, const double* diffconc0,
                  const double* diffconc1, const int& ngdiffconc,
                  const double* flux0, const double* flux1, const int& ngflux,
                  const int* physb)
{
   (void)physb;
   ADD_FLUX(ifirst0, ilast0, ifirst1, ilast1, dx, conc, ngconc, ncomp,
            diffconc0, diffconc1, ngdiffconc, flux0, flux1, ngflux);
}

void add_flux_iso(const int& ifirst0, const int& ilast0, const int& ifirst1,
                  const int& ilast1, const double* dx, const double* conc,
                  const int& ngconc, const int& ncomp, const double* diffconc0,
                  const double* diffconc1, const int& ngdiffconc,
                  const double* flux0, const double* flux1, const int& ngflux,
                  const int* physb)
{
   (void)physb;
   //   assert(ngdiffconc > 0);
   ADD_FLUX_ISO(ifirst0, ilast0, ifirst1, ilast1, dx, conc, ngconc, ncomp,
                diffconc0, diffconc1, ngdiffconc, flux0, flux1, ngflux);
}

void add_flux_4th(const int& ifirst0, const int& ilast0, const int& ifirst1,
                  const int& ilast1, const double* dx, const double* conc,
                  const int& ngconc, const int& ncomp, const double* diffconc0,
                  const double* diffconc1, const int& ngdiffconc,
                  const double* flux0, const double* flux1, const int& ngflux,
                  const int* physb)
{
   //   assert(ngdiffconc > 0);
   ADD_FLUX_4TH(ifirst0, ilast0, ifirst1, ilast1, dx, conc, ngconc, ncomp,
                diffconc0, diffconc1, ngdiffconc, flux0, flux1, ngflux, physb);
}

#endif


EBSCompositionRHSStrategy::EBSCompositionRHSStrategy(
    const int phase_scratch_id, const unsigned short ncompositions,
    const int conc_l_scratch_id, const int conc_a_scratch_id,
    const int conc_b_scratch_id, const int temperature_scratch_id,
    const int diffusion_l_id, const int diffusion_a_id,
    const int diffusion_b_id, const std::vector<int> diffusion_precond_id,
    const std::string& avg_func_type,
    std::shared_ptr<FreeEnergyStrategy> free_energy_strategy,
    std::shared_ptr<CompositionDiffusionStrategy> diffusion_for_conc_in_phase,
    const std::string flux_type)
    : CompositionRHSStrategy(avg_func_type),
      d_ncompositions(ncompositions),
      d_phase_scratch_id(phase_scratch_id),
      d_conc_l_scratch_id(conc_l_scratch_id),
      d_conc_a_scratch_id(conc_a_scratch_id),
      d_conc_b_scratch_id(conc_b_scratch_id),
      d_temperature_scratch_id(temperature_scratch_id),
      d_diffusion_l_id(diffusion_l_id),
      d_diffusion_a_id(diffusion_a_id),
      d_diffusion_b_id(diffusion_b_id),
      d_diffusion_for_conc_in_phase(diffusion_for_conc_in_phase),
      d_free_energy_strategy(free_energy_strategy)
{
   assert(diffusion_l_id >= 0);
   assert(temperature_scratch_id >= 0);

   tbox::plog << "EBSCompositionRHSStrategy" << std::endl;

   d_diffusion_precond_id = diffusion_precond_id;

   d_with_phaseB = (conc_b_scratch_id >= 0);

   d_with_diffusion_for_preconditioner = (!diffusion_precond_id.empty());
   if (flux_type == "isotropic")
      d_add_flux = add_flux_iso;
   else if (flux_type == "order4")
      d_add_flux = add_flux_4th;
   else
      d_add_flux = add_flux_ebs;

   if (d_with_phaseB)
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
                                               d_phase_scratch_id);

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

   std::shared_ptr<pdat::SideData<double> > diffusionl(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(d_diffusion_l_id)));
   assert(diffusionl);
   assert(diffusionl->getDepth() == (d_ncompositions * d_ncompositions));

   std::shared_ptr<pdat::SideData<double> > diffusiona(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(d_diffusion_a_id)));
   assert(diffusiona);
   assert(diffusiona->getDepth() == (int)(d_ncompositions * d_ncompositions));

   std::shared_ptr<pdat::SideData<double> > flux(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(flux_id)));
   assert(flux);
   assert(flux->getDepth() == d_ncompositions);

   flux->fillAll(0.);

   // get CartesianPatchGeometry associated with patch
   std::shared_ptr<geom::CartesianPatchGeometry> pg(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));
   int physbc[6] = {0, 0, 0, 0, 0, 0};
   if (pg->getTouchesRegularBoundary(0, 0)) {
      physbc[0] = 1;
   }
   if (pg->getTouchesRegularBoundary(0, 1)) {
      physbc[1] = 1;
   }
   if (pg->getTouchesRegularBoundary(1, 0)) {
      physbc[2] = 1;
   }
   if (pg->getTouchesRegularBoundary(1, 1)) {
      physbc[3] = 1;
   }
   if (pg->getTouchesRegularBoundary(2, 0)) {
      physbc[4] = 1;
   }
   if (pg->getTouchesRegularBoundary(2, 1)) {
      physbc[5] = 1;
   }


   // now add components of concentration flux,
   // one phase at a time
   d_add_flux(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
              ifirst(2), ilast(2),
#endif
              dx, conc_l->getPointer(), conc_l->getGhostCellWidth()[0],
              d_ncompositions, diffusionl->getPointer(0),
              diffusionl->getPointer(1),
#if (NDIM == 3)
              diffusionl->getPointer(2),
#endif
              diffusionl->getGhostCellWidth()[0], flux->getPointer(0),
              flux->getPointer(1),
#if (NDIM == 3)
              flux->getPointer(2),
#endif
              flux->getGhostCellWidth()[0], physbc);

   d_add_flux(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
              ifirst(2), ilast(2),
#endif
              dx, conc_a->getPointer(), conc_a->getGhostCellWidth()[0],
              d_ncompositions, diffusiona->getPointer(0),
              diffusiona->getPointer(1),
#if (NDIM == 3)
              diffusiona->getPointer(2),
#endif
              diffusiona->getGhostCellWidth()[0], flux->getPointer(0),
              flux->getPointer(1),
#if (NDIM == 3)
              flux->getPointer(2),
#endif
              flux->getGhostCellWidth()[0], physbc);

   if (d_conc_b_scratch_id >= 0) {
      std::shared_ptr<pdat::CellData<double> > conc_b(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(d_conc_b_scratch_id)));
      assert(conc_b);

      std::shared_ptr<pdat::SideData<double> > diffusionb(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch.getPatchData(d_diffusion_b_id)));
      assert(diffusionb);

      assert(diffusionb->getPointer(0) != nullptr);
      assert(diffusionb->getPointer(1) != nullptr);
#if (NDIM == 3)
      assert(diffusionb->getPointer(2) != nullptr);
#endif

      d_add_flux(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                 ifirst(2), ilast(2),
#endif
                 dx, conc_b->getPointer(), conc_b->getGhostCellWidth()[0],
                 d_ncompositions, diffusionb->getPointer(0),
                 diffusionb->getPointer(1),
#if (NDIM == 3)
                 diffusionb->getPointer(2),
#endif
                 diffusionb->getGhostCellWidth()[0], flux->getPointer(0),
                 flux->getPointer(1),
#if (NDIM == 3)
                 flux->getPointer(2),
#endif
                 flux->getGhostCellWidth()[0], physbc);
   }
}

//-----------------------------------------------------------------------

void EBSCompositionRHSStrategy::setDiffusionCoeffForPreconditioner(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   // tbox::pout<<"EBSCompositionRHSStrategy::setDiffusionCoeffForPreconditioner"<<endl;

   assert(d_diffusion_l_id >= 0);
   assert(d_diffusion_a_id >= 0);
   if (d_with_phaseB) {
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

         if (d_with_phaseB) {
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
            assert(diffusion);

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

   const unsigned int num_depth = 1;
   const unsigned int src_depth = depth_in_Dmatrix;
   const unsigned int dst_depth = depth;
   const hier::IntVector src_shift(pbox.getDim(), 0);
   AddOperation<double> addop;

   sd_d_coeff->fill(0., pbox);

   // sd_d_* already include the phase fraction weight, so it is a simple sum
   for (int side_normal = 0; side_normal < NDIM; side_normal++) {
      // create temporary box to pass data indexes to doArrayDataOperationOnBox
      hier::Box sbox(pbox);
      hier::Index up(pbox.upper());
      up[side_normal]++;
      sbox.setUpper(up);

      pdat::ArrayDataOperationUtilities<double, AddOperation<double> >::
          doArrayDataOperationOnBox(sd_d_coeff->getArrayData(side_normal),
                                    sd_d_l->getArrayData(side_normal), sbox,
                                    src_shift, dst_depth, src_depth, num_depth,
                                    addop);
      pdat::ArrayDataOperationUtilities<double, AddOperation<double> >::
          doArrayDataOperationOnBox(sd_d_coeff->getArrayData(side_normal),
                                    sd_d_a->getArrayData(side_normal), sbox,
                                    src_shift, dst_depth, src_depth, num_depth,
                                    addop);
      if (sd_d_b)
         pdat::ArrayDataOperationUtilities<double, AddOperation<double> >::
             doArrayDataOperationOnBox(sd_d_coeff->getArrayData(side_normal),
                                       sd_d_b->getArrayData(side_normal), sbox,
                                       src_shift, dst_depth, src_depth,
                                       num_depth, addop);
   }
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
   const int norderp = phase->getDepth();
   if (d_with_phaseB) {
      std::shared_ptr<pdat::CellData<double> > cb(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(d_conc_b_scratch_id)));
      assert(cb);
      ADDCONCENTRATIONFLUXFROMANTITRAPPING3PHASES(
          ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
          ifirst(2), ilast(2),
#endif
          dx, phase->getPointer(), phase->getGhostCellWidth()[0],
          cl->getPointer(), ca->getPointer(), cb->getPointer(),
          cl->getGhostCellWidth()[0], dphidt->getPointer(),
          dphidt->getGhostCellWidth()[0], alpha, flux->getPointer(0),
          flux->getPointer(1),
#if (NDIM == 3)
          flux->getPointer(2),
#endif
          flux->getGhostCellWidth()[0]);
   } else if (norderp > 1) {
      ADDCONCENTRATIONFLUXFROMANTITRAPPINGMULTIORDERP(
          ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
          ifirst(2), ilast(2),
#endif
          dx, phase->getPointer(), phase->getGhostCellWidth()[0], norderp,
          cl->getPointer(), ca->getPointer(), cl->getGhostCellWidth()[0],
          dphidt->getPointer(), dphidt->getGhostCellWidth()[0], alpha,
          flux->getPointer(0), flux->getPointer(1),
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
          dx, phase->getPointer(), phase->getGhostCellWidth()[0],
          cl->getPointer(), ca->getPointer(), cl->getGhostCellWidth()[0],
          d_ncompositions, dphidt->getPointer(), dphidt->getGhostCellWidth()[0],
          alpha, flux->getPointer(0), flux->getPointer(1),
#if (NDIM == 3)
          flux->getPointer(2),
#endif
          flux->getGhostCellWidth()[0]);
}
