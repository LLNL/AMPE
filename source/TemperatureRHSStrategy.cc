// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#include "TemperatureRHSStrategy.h"


void TemperatureRHSStrategy::evaluateRHS(
    std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_rhs_id, const int dphidt_scratch_id)
{
   assert(temperature_rhs_id >= 0);

   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
           ++ip) {
         std::shared_ptr<hier::Patch> patch = *ip;

         evaluateRHS(patch, temperature_rhs_id, dphidt_scratch_id);
      }
   }
}
