// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.

#include "LimitDifference.h"

#include "FuncFort.h"

#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellData.h"

using namespace SAMRAI;

LimitDifference::LimitDifference(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
    : d_hierarchy(hierarchy)
{
}

void LimitDifference::apply(const int cell_data_id, const int ref_cell_data_id,
                            const int weight_id, const int depth,
                            const double diff)

{
   assert(cell_data_id >= 0);
   assert(ref_cell_data_id >= 0);

   tbox::pout << "Limit difference to " << diff << std::endl;
   const int maxln = d_hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> patch_level =
          d_hierarchy->getPatchLevel(ln);

      apply(patch_level, cell_data_id, ref_cell_data_id, weight_id, depth,
            diff);
   }
}

void LimitDifference::apply(const std::shared_ptr<hier::PatchLevel> level,
                            const int cell_data_id, const int ref_cell_data_id,
                            const int weight_id, const int depth,
                            const double diff)

{
   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      std::shared_ptr<hier::Patch> patch = *p;

      std::shared_ptr<pdat::CellData<double> > data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(cell_data_id)));
      assert(data);
      const hier::Box& gbox = data->getGhostBox();

      std::shared_ptr<pdat::CellData<double> > ref_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(ref_cell_data_id)));
      assert(ref_data);

      std::shared_ptr<pdat::CellData<double> > weight(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(weight_id)));
      assert(weight);

      const double inv1pd = 1. / (1. + diff);
      pdat::CellIterator iend(pdat::CellGeometry::end(gbox));
      for (pdat::CellIterator i(pdat::CellGeometry::begin(gbox)); i != iend;
           ++i) {
         pdat::CellIndex cell = *i;
         assert((*ref_data)(cell) >= 0.);
         assert((*ref_data)(cell) <= 1.);
         const double delta = (*data)(cell) - (*ref_data)(cell);
         const double w = std::max(0., (*weight)(cell, depth));
         const double x = 2. * w;
         const double t = x > 1. ? 1. : x * x * x * (10. - 15. * x + 6. * x);
         const double tt = (t + diff) * inv1pd;
         // const double tt = std::abs(delta)>diff ? diff/delta : 1.;
         (*data)(cell) = (*ref_data)(cell) + tt * delta;
         if ((*data)(cell) > 1.) (*data)(cell) = 1.;
         if ((*data)(cell) < 0.) (*data)(cell) = 0.;
      }
   }
}
