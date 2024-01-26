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
#include "QuatMobilityStrategy.h"

#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

void QuatMobilityStrategy::computePhaseTemperatureMobility(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    int& phase_mobility_id, int& cp_id, int& mobility_id)
{
   assert(phase_mobility_id >= 0);
   assert(cp_id >= 0);
   assert(mobility_id >= 0);

   math::HierarchyCellDataOpsReal<double> ops(hierarchy);

   ops.divide(mobility_id, phase_mobility_id, cp_id, false);
}
