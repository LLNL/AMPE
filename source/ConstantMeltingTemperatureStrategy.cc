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
#include "ConstantMeltingTemperatureStrategy.h"
#include "ConcFort.h"

#include "SAMRAI/pdat/CellData.h"

void ConstantMeltingTemperatureStrategy::evaluate(hier::Patch& patch)
{
   assert(d_equilibrium_temperature_id >= 0);

   std::shared_ptr<pdat::CellData<double> > temperature(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_equilibrium_temperature_id)));

   assert(temperature);

   temperature->fill(d_Tref);
}
