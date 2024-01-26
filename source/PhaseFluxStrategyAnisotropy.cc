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
#include "PhaseFluxStrategyAnisotropy.h"
#include "QuatFort.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/math/PatchSideDataNormOpsReal.h"

void PhaseFluxStrategyAnisotropy::computeFluxes(
    const std::shared_ptr<hier::PatchLevel> level, const int phase_id,
    const int quat_id, const int flux_id)
{
   assert(quat_id >= 0);

   //  Compute phase "flux" on patches in level.
   for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
        ++ip) {
      std::shared_ptr<hier::Patch> patch = *ip;

      const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
          SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                                 hier::PatchGeometry>(
              patch->getPatchGeometry()));
      TBOX_ASSERT(patch_geom);

      const double* dx = patch_geom->getDx();

      std::shared_ptr<pdat::CellData<double> > phase(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(phase_id)));
      TBOX_ASSERT(phase);

      std::shared_ptr<pdat::SideData<double> > phase_flux(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(flux_id)));
      TBOX_ASSERT(phase_flux);

      std::shared_ptr<pdat::CellData<double> > quat(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(quat_id)));
      TBOX_ASSERT(quat);

      const hier::Box& pbox = patch->getBox();
      const hier::Index& ifirst = pbox.lower();
      const hier::Index& ilast = pbox.upper();

      for (int i = 0; i < phase->getDepth(); i++)
         ANISOTROPIC_GRADIENT_FLUX(
             ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
             ifirst(2), ilast(2),
#endif
             dx, d_epsilon_phase, d_nu, d_knumber, phase->getPointer(i),
             phase->getGhostCellWidth()[0], quat->getPointer(),
             quat->getGhostCellWidth()[0], quat->getDepth(),
             phase_flux->getPointer(0, i), phase_flux->getPointer(1, i),
#if (NDIM == 3)
             phase_flux->getPointer(2, i),
#endif
             phase_flux->getGhostCellWidth()[0]);

#ifdef DEBUG_CHECK_ASSERTIONS
      SAMRAI::math::PatchSideDataNormOpsReal<double> ops;
      double l2f = ops.L2Norm(phase_flux, pbox);
      assert(l2f == l2f);
#endif
   }
}
