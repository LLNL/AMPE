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
#include "HierarchyStencilOps.h"

#include "QuatFort.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

void HierarchyStencilOps::computeGradSide(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int diffs_id,
    int grad_side_id)
{
   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         computeGradSide(patch, diffs_id, grad_side_id);
      }
   }
}

void HierarchyStencilOps::computeGradSide(std::shared_ptr<hier::Patch> patch,
                                          int diffs_id, int grad_side_id)
{
   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   std::shared_ptr<pdat::SideData<double> > diff_data(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch->getPatchData(diffs_id)));
   assert(diff_data);
   assert(diff_data->getGhostCellWidth()[0] > 0);

   std::shared_ptr<pdat::SideData<double> > grad_side_data(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch->getPatchData(grad_side_id)));
   assert(grad_side_data);
   assert(grad_side_data->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));

   const hier::Box& diff_gbox = diff_data->getGhostBox();
   const hier::Index& d_lower = diff_gbox.lower();
   const hier::Index& d_upper = diff_gbox.upper();

   const hier::Box& grad_gbox = grad_side_data->getGhostBox();
   const hier::Index& g_lower = grad_gbox.lower();
   const hier::Index& g_upper = grad_gbox.upper();

   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch->getPatchGeometry()));
   TBOX_ASSERT(patch_geom);

   const double* dx = patch_geom->getDx();

   // there is a gradient component for each dimension x,y,z
   assert(grad_side_data->getDepth() == NDIM);

   GRAD_SIDE(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
             ifirst(2), ilast(2),
#endif
             diff_data->getPointer(0), diff_data->getPointer(1),
#if (NDIM == 3)
             diff_data->getPointer(2),
#endif
             d_lower[0], d_upper[0], d_lower[1], d_upper[1],
#if (NDIM == 3)
             d_lower[2], d_upper[2],
#endif
             dx,
             grad_side_data->getPointer(0, 0),  // side 0, depth 0 (x component)
             grad_side_data->getPointer(0, 1),  // side 0, depth 1 (y component)
#if (NDIM == 3)
             grad_side_data->getPointer(0, 2),
#endif
             grad_side_data->getPointer(1, 0), grad_side_data->getPointer(1, 1),
#if (NDIM == 3)
             grad_side_data->getPointer(1, 2),
#endif
#if (NDIM == 3)
             grad_side_data->getPointer(2, 0), grad_side_data->getPointer(2, 1),
             grad_side_data->getPointer(2, 2),
#endif
             g_lower[0], g_upper[0], g_lower[1], g_upper[1]
#if (NDIM == 3)
             ,
             g_lower[2], g_upper[2]
#endif
   );
}

void HierarchyStencilOps::computeDiffs(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int var_id,
    int diffs_id)
{
   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
      computeDiffs(level, var_id, diffs_id);
   }
}

void HierarchyStencilOps::computeDiffs(std::shared_ptr<hier::PatchLevel> level,
                                       int var_id, int diffs_id)
{
   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      std::shared_ptr<hier::Patch> patch = *p;

      computeDiffs(patch, var_id, diffs_id);
   }
}


void HierarchyStencilOps::computeDiffs(std::shared_ptr<hier::Patch> patch,
                                       int var_id, int diffs_id)
{
   const hier::Box& box = patch->getBox();
   const hier::Index& ifirst = box.lower();
   const hier::Index& ilast = box.upper();

   std::shared_ptr<pdat::CellData<double> > var_data(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(var_id)));
   assert(var_data);
   assert(var_data->getGhostCellWidth()[0] > 0);

   std::shared_ptr<pdat::SideData<double> > diff_data(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch->getPatchData(diffs_id)));
   assert(diff_data);
   assert(diff_data->getGhostCellWidth()[0] > 0);

   DIFFS(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
         ifirst(2), ilast(2),
#endif
         var_data->getPointer(), var_data->getGhostCellWidth()[0],
         diff_data->getPointer(0), diff_data->getPointer(1),
#if (NDIM == 3)
         diff_data->getPointer(2),
#endif
         diff_data->getGhostCellWidth()[0]);
}

void HierarchyStencilOps::computeGradCell(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int diffs_id,
    int grad_id)
{
   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> patch_level =
          hierarchy->getPatchLevel(ln);

      computeGradCell(patch_level, diffs_id, grad_id);
   }
}

void HierarchyStencilOps::computeGradCell(
    std::shared_ptr<hier::PatchLevel> level, int diffs_id, int grad_id)
{
   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      std::shared_ptr<hier::Patch> patch = *p;

      computeGradCell(patch, diffs_id, grad_id);
   }
}

void HierarchyStencilOps::computeGradCell(std::shared_ptr<hier::Patch> patch,
                                          int diffs_id, int grad_id)
{
   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   std::shared_ptr<pdat::SideData<double> > diff_data(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch->getPatchData(diffs_id)));
   assert(diff_data);

   std::shared_ptr<pdat::CellData<double> > grad_cell_data(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(grad_id)));
   assert(grad_cell_data);

   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch->getPatchGeometry()));
   TBOX_ASSERT(patch_geom);

   const double* dx = patch_geom->getDx();

   assert(grad_cell_data->getDepth() == NDIM);

   GRAD_CELL(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
             ifirst(2), ilast(2),
#endif
             diff_data->getPointer(0), diff_data->getPointer(1),
#if (NDIM == 3)
             diff_data->getPointer(2),
#endif
             diff_data->getGhostCellWidth()[0], dx,
             grad_cell_data->getPointer(0), grad_cell_data->getPointer(1),
#if (NDIM == 3)
             grad_cell_data->getPointer(2),
#endif
             grad_cell_data->getGhostCellWidth()[0]);
}
