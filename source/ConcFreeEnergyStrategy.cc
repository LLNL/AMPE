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
#include "ConcFreeEnergyStrategy.h"

void ConcFreeEnergyStrategy::computeDerivFreeEnergyLiquid(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int fl_id)
{
   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         computeDerivFreeEnergyLiquid(*patch, temperature_id, fl_id);
      }
   }
}

void ConcFreeEnergyStrategy::computeDerivFreeEnergySolidA(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int fl_id)
{
   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         computeDerivFreeEnergySolidA(*patch, temperature_id, fl_id);
      }
   }
}

void ConcFreeEnergyStrategy::computeDerivFreeEnergySolidB(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int fl_id)
{
   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         computeDerivFreeEnergySolidB(*patch, temperature_id, fl_id);
      }
   }
}
