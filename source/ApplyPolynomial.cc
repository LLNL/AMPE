// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.

#include "ApplyPolynomial.h"

#include "FuncFort.h"

#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellData.h"

ApplyPolynomial::ApplyPolynomial(EnergyInterpolationType interp_func_type)
    : d_interp_func_type(interp_func_type)
{
}

void ApplyPolynomial::apply(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int src_cell_data_id, const int dst_cell_data_id)

{
   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> patch_level =
          hierarchy->getPatchLevel(ln);

      apply(patch_level, src_cell_data_id, dst_cell_data_id);
   }
}

void ApplyPolynomial::apply(const std::shared_ptr<hier::PatchLevel> level,
                            const int src_cell_data_id,
                            const int dst_cell_data_id)

{
   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      std::shared_ptr<hier::Patch> patch = *p;

      std::shared_ptr<pdat::CellData<double> > sdata(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(src_cell_data_id)));
      assert(sdata);
      const hier::Box& gbox = sdata->getGhostBox();

      std::shared_ptr<pdat::CellData<double> > ddata(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(dst_cell_data_id)));
      assert(ddata);

      // issue:mew: Potentially replace this loop with fortran kernel.
      const char interp = energyInterpChar(d_interp_func_type);
      pdat::CellIterator iend(pdat::CellGeometry::end(gbox));
      for (pdat::CellIterator i(pdat::CellGeometry::begin(gbox)); i != iend;
           ++i) {
         pdat::CellIndex cell = *i;
         const double phi = (*sdata)(cell);
         const double hphi = INTERP_FUNC(phi, &interp);
         (*ddata)(cell) = hphi;
      }
   }
}
