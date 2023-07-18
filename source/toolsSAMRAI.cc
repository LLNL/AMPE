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
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/HierarchySideDataOpsReal.h"

using namespace SAMRAI;


void copyDepthSideData(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
                       const int dst_id, const int dst_depth, const int src_id,
                       const int src_depth)
{
   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end();
           ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::SideData<double> > dst(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 p->getPatchData(dst_id)));
         TBOX_ASSERT(dst);

         std::shared_ptr<pdat::SideData<double> > src(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 p->getPatchData(src_id)));
         TBOX_ASSERT(src);

         dst->copyDepth(dst_depth, *src, src_depth);
      }
   }
}

void copyDepthCellData(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
                       const int dst_id, const int dst_depth, const int src_id,
                       const int src_depth)
{
   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end();
           ++ip) {
         const std::shared_ptr<hier::Patch>& p = *ip;

         std::shared_ptr<pdat::CellData<double> > dst(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 p->getPatchData(dst_id)));
         TBOX_ASSERT(dst);

         std::shared_ptr<pdat::CellData<double> > src(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 p->getPatchData(src_id)));
         TBOX_ASSERT(src);

         dst->copyDepth(dst_depth, *src, src_depth);
      }
   }
}

int checkForNans(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
                 const int data_id)
{
#ifndef NDEBUG
   math::HierarchyCellDataOpsReal<double> mathops(hierarchy);
   const double norm_data = mathops.L1Norm(data_id);

   if (norm_data != norm_data) return 1;
#else
   (void)hierarchy;
   (void)data_id;
#endif
   return 0;
}

int checkSideDataForNans(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
                         const int data_id)
{
#ifndef NDEBUG
   math::HierarchySideDataOpsReal<double> mathops(hierarchy);
   const double norm_data = mathops.L1Norm(data_id);

   if (norm_data != norm_data) return 1;
#else
   (void)hierarchy;
   (void)data_id;
#endif
   return 0;
}
