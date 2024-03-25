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
#include "QuatGradModulusStrategy.h"

#include "QuatFort.h"

#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"

void QuatGradModulusStrategy::computeQuatGradModulus(
    const std::shared_ptr<hier::PatchLevel> level, int& grad_cell_id,
    int& grad_modulus_id)
{
   assert(d_qlen > 0);
   assert(grad_cell_id >= 0);
   assert(grad_modulus_id >= 0);

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      std::shared_ptr<hier::Patch> patch = *p;
      const hier::Box& pbox = patch->getBox();
      const hier::Index& ifirst = pbox.lower();
      const hier::Index& ilast = pbox.upper();

      std::shared_ptr<pdat::CellData<double> > grad_cell_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(grad_cell_id)));
      assert(grad_cell_data);
      assert(grad_cell_data->getGhostCellWidth() ==
             hier::IntVector(tbox::Dimension(NDIM), 0));

      std::shared_ptr<pdat::CellData<double> > grad_modulus_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(grad_modulus_id)));
      assert(grad_modulus_data);
      assert(grad_modulus_data->getGhostCellWidth() ==
             hier::IntVector(tbox::Dimension(NDIM), 0));

      assert(grad_cell_data->getDepth() == NDIM * d_qlen);
      assert(grad_modulus_data->getDepth() == 1);

      QUATGRAD_MODULUS(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                       ifirst(2), ilast(2),
#endif
                       d_qlen, grad_cell_data->getPointer(0 * d_qlen),
                       grad_cell_data->getPointer(1 * d_qlen),
#if (NDIM == 3)
                       grad_cell_data->getPointer(2 * d_qlen),
#endif
                       grad_cell_data->getGhostCellWidth()[0],
                       grad_modulus_data->getPointer(),
                       grad_modulus_data->getGhostCellWidth()[0]);
   }
}

void QuatGradModulusStrategy::computeQuatGradModulusFromSides(
    const std::shared_ptr<hier::PatchLevel> level, int& grad_side_id,
    int& grad_modulus_id)
{
   assert(grad_side_id >= 0);
   assert(grad_modulus_id >= 0);

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      std::shared_ptr<hier::Patch> patch = *p;
      const hier::Box& pbox = patch->getBox();
      const hier::Index& ifirst = pbox.lower();
      const hier::Index& ilast = pbox.upper();

      std::shared_ptr<pdat::SideData<double> > grad_side_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(grad_side_id)));
      assert(grad_side_data);
      assert(grad_side_data->getGhostCellWidth() ==
             hier::IntVector(tbox::Dimension(NDIM), 0));

      std::shared_ptr<pdat::CellData<double> > grad_modulus_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(grad_modulus_id)));
      assert(grad_modulus_data);
      assert(grad_modulus_data->getGhostCellWidth() ==
             hier::IntVector(tbox::Dimension(NDIM), 0));

      assert(grad_side_data->getDepth() == NDIM * d_qlen);
      assert(grad_modulus_data->getDepth() == 1);

      QUATGRAD_MODULUS_FROM_SIDES_COMPACT(
          ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
          ifirst(2), ilast(2),
#endif
          d_qlen, grad_side_data->getPointer(0), grad_side_data->getPointer(1),
#if (NDIM == 3)
          grad_side_data->getPointer(2),
#endif
          grad_side_data->getGhostCellWidth()[0],
          grad_modulus_data->getPointer(),
          grad_modulus_data->getGhostCellWidth()[0]);
   }
}
