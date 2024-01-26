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
#include "ConstantHeatCapacityStrategy.h"

#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

ConstantHeatCapacityStrategy::ConstantHeatCapacityStrategy(const double cp,
                                                           const int cp_id)
    : d_cp(cp), d_cp_id(cp_id)
{
   assert(d_cp_id >= 0);
   assert(cp > 0.);
}

void ConstantHeatCapacityStrategy::setCurrentValue(
    std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   assert(d_cp_id >= 0);

   static bool first_time = true;

   if (first_time) {
      math::HierarchyCellDataOpsReal<double> mathops(hierarchy);
      mathops.setToScalar(d_cp_id, d_cp);
   }

   first_time = false;
}
