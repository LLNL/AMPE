// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE

#include "MovingFrameRHS.h"
#include "QuatFort.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"

MovingFrameRHS::MovingFrameRHS(const int data_id) : d_data_id(data_id) {}

void MovingFrameRHS::addRHS(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                            const int ydot_id, const double frame_velocity)
{
   // add component related to moving frame if moving velocity!=0
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
           ++ip) {

         std::shared_ptr<hier::Patch> patch = *ip;

         const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
             SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                                    hier::PatchGeometry>(
                 patch->getPatchGeometry()));
         const double* dx = patch_geom->getDx();

         const hier::Box& pbox = patch->getBox();
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast = pbox.upper();

         std::shared_ptr<pdat::CellData<double> > data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_data_id)));
         assert(data);

         std::shared_ptr<pdat::CellData<double> > data_rhs(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(ydot_id)));
         assert(data_rhs);
         assert(data->getGhostCellWidth()[0] > 0);
         assert(data->getDepth() == data_rhs->getDepth());

         const short nfield = data->getDepth();
         // std::cout<<"ADDVDPHIDX for "<<ndatas<<" datas"<<std::endl;

         for (short i = 0; i < nfield; i++)
            ADDVDPHIDX(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                       ifirst(2), ilast(2),
#endif
                       dx, data->getPointer(i), data->getGhostCellWidth()[0],
                       frame_velocity, data_rhs->getPointer(i),
                       data_rhs->getGhostCellWidth()[0]);
      }
   }
}
