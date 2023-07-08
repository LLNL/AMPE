// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/Patch.h"

#include "QuatFort.h"

using namespace SAMRAI;

void computeQGradSide(std::shared_ptr<pdat::SideData<double> > diff_data,
                      std::shared_ptr<pdat::SideData<double> > grad_data,
                      const double dx[NDIM], const bool isSymmetryAware,
                      const bool isotropicStencil,
                      std::shared_ptr<pdat::SideData<int> > rotation_index)
{
   assert(diff_data);
   assert(grad_data);

   int qlen = diff_data->getDepth();
   if (isSymmetryAware) qlen /= 2;

   const hier::Box& box = grad_data->getBox();
   const hier::Index& ifirst = box.lower();
   const hier::Index& ilast = box.upper();

   if (isSymmetryAware) {
      assert(rotation_index);

      if (!isotropicStencil) {

         // tbox::pout<<"compute quat grad at sides..."<<endl;

         QUATGRAD_SIDE_SYMM(
             ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
             ifirst(2), ilast(2),
#endif
             qlen, dx, diff_data->getPointer(0, 0), diff_data->getPointer(1, 0),
#if (NDIM == 3)
             diff_data->getPointer(2, 0),
#endif
             diff_data->getGhostCellWidth()[0],
             grad_data->getPointer(0, 0 * qlen),  // grad_x, xside
             grad_data->getPointer(0, 1 * qlen),  // grad_y, xside
#if (NDIM == 3)
             grad_data->getPointer(0, 2 * qlen),  // grad_z, xside
#endif
             grad_data->getPointer(1, 0 * qlen),  // grad_x, yside
             grad_data->getPointer(1, 1 * qlen),  // grad_y, yside
#if (NDIM == 3)
             grad_data->getPointer(1, 2 * qlen),  // grad_z, yside
#endif
#if (NDIM == 3)
             grad_data->getPointer(2, 0 * qlen),
             grad_data->getPointer(2, 1 * qlen),
             grad_data->getPointer(2, 2 * qlen),
#endif
             grad_data->getGhostCellWidth()[0], rotation_index->getPointer(0),
             rotation_index->getPointer(1),
#if (NDIM == 3)
             rotation_index->getPointer(2),
#endif
             rotation_index->getGhostCellWidth()[0]);
      } else {
         QUATGRAD_SIDE_SYMM_ISOTROPIC(
             ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
             ifirst(2), ilast(2),
#endif
             qlen, dx, diff_data->getPointer(0, 0), diff_data->getPointer(1, 0),
#if (NDIM == 3)
             diff_data->getPointer(2, 0),
#endif
             diff_data->getGhostCellWidth()[0],
             grad_data->getPointer(0, 0 * qlen),  // grad_x, xside
             grad_data->getPointer(0, 1 * qlen),  // grad_y, xside
#if (NDIM == 3)
             grad_data->getPointer(0, 2 * qlen),  // grad_z, xside
#endif
             grad_data->getPointer(1, 0 * qlen),  // grad_x, yside
             grad_data->getPointer(1, 1 * qlen),  // grad_y, yside
#if (NDIM == 3)
             grad_data->getPointer(1, 2 * qlen),  // grad_z, yside
#endif
#if (NDIM == 3)
             grad_data->getPointer(2, 0 * qlen),
             grad_data->getPointer(2, 1 * qlen),
             grad_data->getPointer(2, 2 * qlen),
#endif
             grad_data->getGhostCellWidth()[0], rotation_index->getPointer(0),
             rotation_index->getPointer(1),
#if (NDIM == 3)
             rotation_index->getPointer(2),
#endif
             rotation_index->getGhostCellWidth()[0]);
      }
   } else {
      if (!isotropicStencil) {
         QUATGRAD_SIDE(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                       ifirst(2), ilast(2),
#endif
                       qlen, dx, diff_data->getPointer(0, 0),
                       diff_data->getPointer(1, 0),
#if (NDIM == 3)
                       diff_data->getPointer(2, 0),
#endif
                       diff_data->getGhostCellWidth()[0],
                       // output
                       grad_data->getPointer(0, 0 * qlen),  // grad_x, xside
                       grad_data->getPointer(0, 1 * qlen),  // grad_y, xside
#if (NDIM == 3)
                       grad_data->getPointer(0, 2 * qlen),  // grad_z, xside
#endif
                       grad_data->getPointer(1, 0 * qlen),  // grad_x, yside
                       grad_data->getPointer(1, 1 * qlen),  // grad_y, yside
#if (NDIM == 3)
                       grad_data->getPointer(1, 2 * qlen),  // grad_z, yside
#endif
#if (NDIM == 3)
                       grad_data->getPointer(2, 0 * qlen),
                       grad_data->getPointer(2, 1 * qlen),
                       grad_data->getPointer(2, 2 * qlen),
#endif
                       grad_data->getGhostCellWidth()[0]);
      } else {
         QUATGRAD_SIDE_ISOTROPIC(
             ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
             ifirst(2), ilast(2),
#endif
             qlen, dx, diff_data->getPointer(0, 0), diff_data->getPointer(1, 0),
#if (NDIM == 3)
             diff_data->getPointer(2, 0),
#endif
             diff_data->getGhostCellWidth()[0],
             // output
             grad_data->getPointer(0, 0 * qlen),  // grad_x, xside
             grad_data->getPointer(0, 1 * qlen),  // grad_y, xside
#if (NDIM == 3)
             grad_data->getPointer(0, 2 * qlen),  // grad_z, xside
#endif
             grad_data->getPointer(1, 0 * qlen),  // grad_x, yside
             grad_data->getPointer(1, 1 * qlen),  // grad_y, yside
#if (NDIM == 3)
             grad_data->getPointer(1, 2 * qlen),  // grad_z, yside
#endif
#if (NDIM == 3)
             grad_data->getPointer(2, 0 * qlen),
             grad_data->getPointer(2, 1 * qlen),
             grad_data->getPointer(2, 2 * qlen),
#endif
             grad_data->getGhostCellWidth()[0]);
      }
   }
}

void computeQGrad(std::shared_ptr<pdat::SideData<double> > diff_data,
                  std::shared_ptr<pdat::CellData<double> > grad_data,
                  const double dx[NDIM], const bool isSymmetryAware,
                  std::shared_ptr<pdat::SideData<int> > rotation_index)
{
   assert(diff_data);
   assert(grad_data);

   int qlen = diff_data->getDepth();
   if (isSymmetryAware) qlen /= 2;

   const hier::Box& box = grad_data->getBox();
   const hier::Index& ifirst = box.lower();
   const hier::Index& ilast = box.upper();

   if (isSymmetryAware) {
      assert(rotation_index);

      QUATGRAD_CELL_SYMM(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                         ifirst(2), ilast(2),
#endif
                         qlen, dx, diff_data->getPointer(0, 0),
                         diff_data->getPointer(1, 0),
#if (NDIM == 3)
                         diff_data->getPointer(2, 0),
#endif
                         diff_data->getGhostCellWidth()[0],
                         grad_data->getPointer(0 * qlen),
                         grad_data->getPointer(1 * qlen),
#if (NDIM == 3)
                         grad_data->getPointer(2 * qlen),
#endif
                         grad_data->getGhostCellWidth()[0],
                         rotation_index->getPointer(0),
                         rotation_index->getPointer(1),
#if (NDIM == 3)
                         rotation_index->getPointer(2),
#endif
                         rotation_index->getGhostCellWidth()[0]);

   } else {
      QUATGRAD_CELL(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                    ifirst(2), ilast(2),
#endif
                    qlen, dx, diff_data->getPointer(0, 0),
                    diff_data->getPointer(1, 0),
#if (NDIM == 3)
                    diff_data->getPointer(2, 0),
#endif
                    diff_data->getGhostCellWidth()[0],
                    grad_data->getPointer(0 * qlen),
                    grad_data->getPointer(1 * qlen),
#if (NDIM == 3)
                    grad_data->getPointer(2 * qlen),
#endif
                    grad_data->getGhostCellWidth()[0]);
   }
}

void computeQDiffs(std::shared_ptr<pdat::CellData<double> > quat_data,
                   std::shared_ptr<pdat::SideData<double> > diff_data,
                   const bool isSymmetryAware,
                   std::shared_ptr<pdat::SideData<int> > rotation_index)
{
   assert(quat_data);
   assert(diff_data);
   if (isSymmetryAware)
      assert(2 * quat_data->getDepth() == diff_data->getDepth());
   else
      assert(quat_data->getDepth() == diff_data->getDepth());
   assert(quat_data->getGhostCellWidth()[0] > 0);

   const int qlen = quat_data->getDepth();

   const hier::Box& box = quat_data->getBox();
   const hier::Index& ifirst = box.lower();
   const hier::Index& ilast = box.upper();

   // If symmetry is on, "symmetric diffs" are stored first, with
   // "nonsymmetric diffs" stored offset by qlen.
   const int symm_depth_offset = 0;
   const int nonsymm_depth_offset = isSymmetryAware ? qlen : 0;

   QUATDIFFS(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
             ifirst(2), ilast(2),
#endif
             qlen, quat_data->getPointer(), quat_data->getGhostCellWidth()[0],
             diff_data->getPointer(0, nonsymm_depth_offset),
             diff_data->getPointer(1, nonsymm_depth_offset),
#if (NDIM == 3)
             diff_data->getPointer(2, nonsymm_depth_offset),
#endif
             diff_data->getGhostCellWidth()[0]);

   if (isSymmetryAware) {
      assert(rotation_index);

      QUATDIFFS_SYMM(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                     ifirst(2), ilast(2),
#endif
                     qlen, quat_data->getPointer(),
                     quat_data->getGhostCellWidth()[0],
                     diff_data->getPointer(0, symm_depth_offset),
                     diff_data->getPointer(1, symm_depth_offset),
#if (NDIM == 3)
                     diff_data->getPointer(2, symm_depth_offset),
#endif
                     diff_data->getGhostCellWidth()[0],
                     rotation_index->getPointer(0),
                     rotation_index->getPointer(1),
#if (NDIM == 3)
                     rotation_index->getPointer(2),
#endif
                     rotation_index->getGhostCellWidth()[0]);
   }
}
