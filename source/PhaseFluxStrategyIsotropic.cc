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
#include "PhaseFluxStrategyIsotropic.h"
#include "QuatFort.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"

void PhaseFluxStrategyIsotropic::computeFluxes(
    const std::shared_ptr<hier::PatchLevel> level, const int phase_id,
    const int quat_id, const int flux_id)
{
   // this strategy is independent of grain orientation
   (void)quat_id;

   //  Compute phase "flux" on patches in level.
   for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
        ++ip) {
      std::shared_ptr<hier::Patch> patch = *ip;

      const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
          SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                                 hier::PatchGeometry>(
              patch->getPatchGeometry()));
      const double* dx = patch_geom->getDx();

      std::shared_ptr<pdat::CellData<double> > phase(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(phase_id)));

      assert(phase->getGhostCellWidth()[0] > 0);

      std::shared_ptr<pdat::SideData<double> > phase_flux(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(flux_id)));

      const hier::Box& pbox = patch->getBox();
      const hier::Index& ifirst = pbox.lower();
      const hier::Index& ilast = pbox.upper();

      for (int i = 0; i < phase->getDepth(); i++)
         COMPUTE_FLUX_ISOTROPIC(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                                ifirst(2), ilast(2),
#endif
                                dx, d_epsilon_phase, phase->getPointer(i),
                                phase->getGhostCellWidth()[0],
                                phase_flux->getPointer(0, i),
                                phase_flux->getPointer(1, i),
#if (NDIM == 3)
                                phase_flux->getPointer(2, i),
#endif
                                phase_flux->getGhostCellWidth()[0]);
   }
}
