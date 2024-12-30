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
#include "SAMRAI/math/PatchCellDataNormOpsReal.h"

#include "QuatFort.h"

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

double integralDepthCellData(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, const int data_id,
    const int depth, const int weight_id)
{
   math::PatchCellDataNormOpsReal<double> ops;

   double integral = 0.;
   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end();
           ++ip) {
         const std::shared_ptr<hier::Patch>& patch = *ip;
         const hier::Box& pbox = patch->getBox();

         std::shared_ptr<pdat::CellData<double> > data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(data_id)));
         TBOX_ASSERT(data);

         std::shared_ptr<pdat::CellData<double> > vol(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(weight_id)));
         TBOX_ASSERT(vol);

         std::shared_ptr<pdat::CellData<double> > data_copy;
         data_copy.reset(new pdat::CellData<double>(
             pbox, 1, hier::IntVector(tbox::Dimension(NDIM), 0)));
         data_copy->copyDepth(0, *data, depth);
         integral += ops.integral(data_copy, pbox, vol);
      }
   }

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   mpi.AllReduce(&integral, 1, MPI_SUM);

   return integral;
}

void sideToCell(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
                const int cdata_id, const int cdepth, const int sdata_id,
                const int sdepth)
{
   assert(cdata_id >= 0);
   assert(sdata_id >= 0);

   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end();
           ++ip) {
         const std::shared_ptr<hier::Patch>& patch = *ip;
         const hier::Box& pbox = patch->getBox();
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast = pbox.upper();

         std::shared_ptr<pdat::CellData<double> > cdata(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(cdata_id)));
         assert(cdata);
         assert(cdepth < cdata->getDepth());

         std::shared_ptr<pdat::SideData<double> > sdata(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(sdata_id)));
         assert(sdata);
         assert(sdepth < sdata->getDepth());

         SIDE2CELL(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                   ifirst(2), ilast(2),
#endif
                   sdata->getPointer(0, sdepth), sdata->getPointer(1, sdepth),
#if (NDIM == 3)
                   sdata->getPointer(2, sdepth),
#endif
                   sdata->getGhostCellWidth()[0], cdata->getPointer(cdepth),
                   cdata->getGhostCellWidth()[0]);
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

void checkCellDataFieldValues(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy, const int data_id)
{
   math::HierarchyCellDataOpsReal<double> mathops(hierarchy);

   tbox::plog << "--- Field check ---" << std::endl;
   tbox::plog << "Data id: " << data_id << std::endl;
   const double max = mathops.max(data_id);
   tbox::plog << "Max. field value: " << max << std::endl;
   const double min = mathops.min(data_id);
   tbox::plog << "Min. field value: " << min << std::endl;

   assert(!std::isnan(max));
   assert(!std::isnan(min));
}
