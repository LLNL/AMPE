// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE

#include "RigidBodyMotionConcRHS.h"
#include "QuatFort.h"
#include "toolsSAMRAI.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"


RigidBodyMotionConcRHS::RigidBodyMotionConcRHS(const int data_id,
                                               const int weight_id,
                                               const double mobility)
    : d_data_id(data_id), d_weight_id(weight_id), d_mobility(mobility)
{
}

void RigidBodyMotionConcRHS::computeVolumes(
    std::shared_ptr<hier::PatchHierarchy> hierarchy, const unsigned nbodies)
{
   d_volumes.resize(nbodies);
   for (unsigned i = 0; i < nbodies; i++) {
      d_volumes[i] =
          integralDepthCellData(hierarchy, d_data_id, i, d_weight_id);
   }
}

// add component related to moving frame if moving velocity!=0
void RigidBodyMotionConcRHS::addRHS(
    std::shared_ptr<hier::PatchHierarchy> hierarchy, const int ydot_id,
    const std::vector<std::array<double, NDIM>>& forces)
{
   computeVolumes(hierarchy, forces.size());

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
         assert(data_rhs->getDepth() == 1);

         const short nfield = data->getDepth() - 1;

         // div ( c v ) = div ( sum_i eta_i * v_i )
         for (short i = 0; i < nfield; i++) {
            const double* const velocity = velocities[i].data();
            assert(!std::isnan(velocity[0]));
            ADDRBMOTION(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                        ifirst(2), ilast(2),
#endif
                        dx, data->getPointer(i), data->getGhostCellWidth()[0],
                        velocity, data_rhs->getPointer(0),
                        data_rhs->getGhostCellWidth()[0]);
         }
      }
   }
}
