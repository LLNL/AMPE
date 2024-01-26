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
#ifndef MeltingTemperatureStrategy_H
#define MeltingTemperatureStrategy_H

#include "SAMRAI/hier/Patch.h"

using namespace SAMRAI;

class MeltingTemperatureStrategy
{
 public:
   MeltingTemperatureStrategy() {}

   virtual ~MeltingTemperatureStrategy() {}

   virtual void evaluate(hier::Patch& patch) = 0;

   virtual int equilibrium_temperature_id() = 0;
};

#endif
