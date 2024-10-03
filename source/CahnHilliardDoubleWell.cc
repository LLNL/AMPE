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
#include "CahnHilliardDoubleWell.h"
#include "QuatFort.h"
#include "ConcFort.h"
#include "FuncFort.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <cassert>

const double gas_constant_R_JpKpmol = GASCONSTANT_R_JPKPMOL;

CahnHilliardDoubleWell::CahnHilliardDoubleWell(
    const int conc_id, const int temperature_id, const int diffusion_id,
    const double mobility, const double ca, const double cb, const double kappa,
    const double well_scale, const std::string& avg_func_type)
    : CompositionRHSStrategy(avg_func_type),
      d_conc_id(conc_id),
      d_temperature_id(temperature_id),
      d_diffusion_id(diffusion_id),
      d_mobility(mobility),
      d_ca(ca),
      d_cb(cb),
      d_kappa(kappa),
      d_well_scale(well_scale)
{
   assert(conc_id >= 0);
   assert(mobility > 0.);
   assert(d_kappa > 0.);
   assert(d_well_scale > 0.);

   tbox::TimerManager* tman = tbox::TimerManager::getManager();
   t_set_diffcoeff_timer =
       tman->getTimer("AMPE::CahnHilliardDoubleWell::setDiffusionCoeff()");
}

//-----------------------------------------------------------------------

void CahnHilliardDoubleWell::setDiffusionCoeff(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time)
{
   (void)time;

   assert(d_diffusion_id >= 0);
   assert(hierarchy);

   t_set_diffcoeff_timer->start();

   // tbox::pout<<"QuatIntegrator::setDiffCoeffForConcentration"<<endl;
   for (int amr_level = 0; amr_level < maxl; amr_level++) {
      std::shared_ptr<hier::PatchLevel> level =
          hierarchy->getPatchLevel(amr_level);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         const hier::Box& pbox = patch->getBox();
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast = pbox.upper();

         std::shared_ptr<pdat::CellData<double> > temperature(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_temperature_id)));

         std::shared_ptr<pdat::SideData<double> > diffusion(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(d_diffusion_id)));
         assert(diffusion->getDepth() == 1);

         diffusion->fillAll(0.);
         addDiffusionCoeff(temperature, diffusion);
      }
   }

   t_set_diffcoeff_timer->stop();
}

//-----------------------------------------------------------------------

void CahnHilliardDoubleWell::addDiffusionCoeff(
    std::shared_ptr<pdat::CellData<double> > temperature,
    std::shared_ptr<pdat::SideData<double> > diffusion, const hier::Box& pbox)
{
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   DIFFUSION_OF_TEMPERATURE(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                            ifirst(2), ilast(2),
#endif
                            temperature->getPointer(),
                            temperature->getGhostCellWidth()[0],
                            diffusion->getPointer(0), diffusion->getPointer(1),
#if (NDIM == 3)
                            diffusion->getPointer(2),
#endif
                            d_d0, d_q0, gas_constant_R_JpKpmol);
}

//-----------------------------------------------------------------------

void CahnHilliardDoubleWell::computeFluxOnPatch(hier::Patch& patch,
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
           patch.getPatchData(d_conc_id)));
   assert(conc);
   assert(conc->getGhostCellWidth()[0] > 1);

   std::shared_ptr<pdat::SideData<double> > flux(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(flux_id)));
   assert(flux);
   assert(flux->getDepth() == conc->getDepth());
   assert(flux->getGhostCellWidth()[0] > 0);

   // now compute concentration flux
   // D_phi*grad phi+D_conc*grad conc
   flux->fillAll(0.);
   for (int ic = 0; ic < conc->getDepth(); ++ic)
      ADD_CAHNHILLIARDDOUBLEWELL_FLUX(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                                      ifirst(2), ilast(2),
#endif
                                      dx, conc->getPointer(ic),
                                      conc->getGhostCellWidth()[0], d_mobility,
                                      d_ca, d_cb, d_well_scale, d_kappa,
                                      flux->getPointer(0, ic),
                                      flux->getPointer(1, ic),
#if (NDIM == 3)
                                      flux->getPointer(2, ic),
#endif
                                      flux->getGhostCellWidth()[0]);
}
