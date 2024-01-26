// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#include "QuatFort.h"
#include "SimpleTemperatureRHSStrategy.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/PatchCellDataBasicOps.h"

using namespace SAMRAI;

SimpleTemperatureRHSStrategy::SimpleTemperatureRHSStrategy(
    const double thermal_diffusivity, const double latent_heat,
    const int temperature_scratch_id, const int cp_id)
    : d_thermal_diffusivity(thermal_diffusivity),
      d_latent_heat(latent_heat),
      d_temperature_scratch_id(temperature_scratch_id),
      d_cp_id(cp_id)
{
   assert(d_temperature_scratch_id >= 0);
   assert(d_cp_id >= 0);
}

void SimpleTemperatureRHSStrategy::evaluateRHS(
    std::shared_ptr<hier::Patch> patch, const int temperature_rhs_id,
    const int dphidt_scratch_id)
{
   assert(temperature_rhs_id >= 0);
   assert(d_cp_id >= 0);
   assert(d_latent_heat > 0.);
   assert(d_latent_heat < 1.e32);
   assert(d_temperature_scratch_id >= 0);

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch->getPatchGeometry()));
   const double* dx = patch_geom->getDx();

   std::shared_ptr<pdat::CellData<double> > temperature(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_temperature_scratch_id)));
   assert(temperature);
   assert(temperature->getGhostCellWidth()[0] > 0);

   std::shared_ptr<pdat::CellData<double> > cp(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_cp_id)));
   assert(cp);

   std::shared_ptr<pdat::CellData<double> > temperature_rhs(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(temperature_rhs_id)));
   assert(temperature_rhs);

   double* phase_rhs_ptr = nullptr;
   int phase_rhs_nghosts = 0;
   const bool with_phase = (dphidt_scratch_id > -1);
   if (with_phase) {
      std::shared_ptr<pdat::CellData<double> > phase_rhs(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(dphidt_scratch_id)));
      assert(phase_rhs);
      phase_rhs_ptr = phase_rhs->getPointer();
      phase_rhs_nghosts = phase_rhs->getGhostCellWidth()[0];
   }

   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

#ifdef DEBUG_CHECK_ASSERTIONS
   math::PatchCellDataBasicOps<double> mathops;
   const double mincp = mathops.min(cp, pbox);
   tbox::plog << "mincp=" << mincp << std::endl;
   assert(mincp > 0.);
#endif

   COMPUTERHSTEMP(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                  ifirst(2), ilast(2),
#endif
                  dx, d_thermal_diffusivity, d_latent_heat,
                  temperature->getPointer(),
                  temperature->getGhostCellWidth()[0], cp->getPointer(),
                  cp->getGhostCellWidth()[0], with_phase, phase_rhs_ptr,
                  phase_rhs_nghosts, temperature_rhs->getPointer(),
                  temperature_rhs->getGhostCellWidth()[0]);
}
