// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE

#include "RigidBodyMotionRHS.h"
#include "QuatFort.h"
#include "toolsSAMRAI.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"


RigidBodyMotionRHS::RigidBodyMotionRHS(const int data_id, const int norderp,
                                       const int weight_id,
                                       const double mobility)
    : d_data_id(data_id),
      d_norderp(norderp),
      d_weight_id(weight_id),
      d_mobility(mobility)
{
   tbox::plog << "RigidBodyMotionRHS with norderp = " << norderp << std::endl;
}

void RigidBodyMotionRHS::computeVolumes(
    std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   d_volumes.resize(d_norderp);
   for (unsigned i = 0; i < d_norderp; i++) {
      d_volumes[i] =
          integralDepthCellData(hierarchy, d_data_id, i, d_weight_id);
   }
}

// add component related to moving frame if moving velocity!=0
void RigidBodyMotionRHS::addRHS(
    std::shared_ptr<hier::PatchHierarchy> hierarchy, const int ydot_id,
    const std::vector<std::array<double, NDIM>>& forces)
{
   assert(forces.size() == d_norderp);

   computeVolumes(hierarchy);

   assert(forces.size() == d_volumes.size());
   for (unsigned i = 0; i < d_volumes.size(); i++)
      assert(d_volumes[i] > 0.);

   std::vector<std::array<double, NDIM>> velocities(forces);
   for (unsigned i = 0; i < forces.size(); i++)
      for (int j = 0; j < NDIM; j++)
         velocities[i][j] *= (d_mobility / d_volumes[i]);

   // for (unsigned i = 0; i < forces.size(); i++)
   //   for (int j = 0; j < NDIM; j++)
   //      std::cout << "velocity = " << velocities[i][j] << std::endl;

   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
           ++ip) {

         std::shared_ptr<hier::Patch> patch = *ip;

         // get CartesianPatchGeometry associated with patch
         const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
             SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                                    hier::PatchGeometry>(
                 patch->getPatchGeometry()));
         const double* dx = patch_geom->getDx();

         const hier::Box& pbox = patch->getBox();
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast = pbox.upper();

         std::shared_ptr<pdat::CellData<double>> data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_data_id)));
         assert(data);

         std::shared_ptr<pdat::CellData<double>> data_rhs(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(ydot_id)));
         assert(data_rhs);
         assert(data->getGhostCellWidth()[0] > 0);
         assert(data->getDepth() == data_rhs->getDepth());
         assert(data->getDepth() >= d_norderp);

         for (short i = 0; i < d_norderp; i++) {
            const double* const velocity = velocities[i].data();
            ADDRBMOTION(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                        ifirst(2), ilast(2),
#endif
                        dx, data->getPointer(i), data->getGhostCellWidth()[0],
                        velocity, data_rhs->getPointer(i),
                        data_rhs->getGhostCellWidth()[0]);
         }
      }
   }
}
