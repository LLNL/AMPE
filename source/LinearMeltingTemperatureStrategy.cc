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
//
#include "LinearMeltingTemperatureStrategy.h"
#include "ConcFort.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/PatchCellDataNormOpsReal.h"

LinearMeltingTemperatureStrategy::LinearMeltingTemperatureStrategy(
    const double Tref, const double c0, const double liquidus_slope,
    const int concentration_id, const int equilibrium_temperature_id)
    : d_Tref(Tref),
      d_c0(c0),
      d_liquidus_slope(liquidus_slope),
      d_concentration_id(concentration_id),
      d_equilibrium_temperature_id(equilibrium_temperature_id)
{
   tbox::plog << "uses LinearMeltingTemperatureStrategy with Tref = " << Tref
              << "..." << std::endl;

   assert(d_concentration_id >= 0);
   assert(d_equilibrium_temperature_id >= 0);
   assert(d_Tref >= 0.);
   assert(d_Tref < 100000.);
}

void LinearMeltingTemperatureStrategy::evaluate(hier::Patch& patch)
{
   assert(d_equilibrium_temperature_id >= 0);
   assert(d_concentration_id >= 0);

   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   std::shared_ptr<pdat::CellData<double> > temperature(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_equilibrium_temperature_id)));
   std::shared_ptr<pdat::CellData<double> > concentration(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_concentration_id)));

   assert(temperature);
   assert(concentration);

#ifdef DEBUG_CHECK_ASSERTIONS
   SAMRAI::math::PatchCellDataNormOpsReal<double> ops;
   double l2rhs = ops.L2Norm(concentration, patch.getBox());
   assert(l2rhs == l2rhs);
#endif

   LINEARMELTINGLINE(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                     ifirst(2), ilast(2),
#endif
                     temperature->getPointer(),
                     temperature->getGhostCellWidth()[0],
                     concentration->getPointer(),
                     concentration->getGhostCellWidth()[0], d_Tref, d_c0,
                     d_liquidus_slope);

#ifdef DEBUG_CHECK_ASSERTIONS
   l2rhs = ops.L2Norm(temperature, patch.getBox());
   assert(l2rhs == l2rhs);
#endif
}
