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
#include "FreeEnergyStrategy.h"

#include "SAMRAI/pdat/CellData.h"

void FreeEnergyStrategy::computeFreeEnergyLiquid(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int fl_id, const bool gp)
{
   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         computeFreeEnergyLiquid(*patch, temperature_id, fl_id, gp);
      }
   }
}

void FreeEnergyStrategy::computeFreeEnergySolidA(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int fa_id, const bool gp)
{
   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         computeFreeEnergySolidA(*patch, temperature_id, fa_id, gp);
      }
   }
}

void FreeEnergyStrategy::computeDrivingForce(
    const double time, const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int phase_id, const int conc_id,
    const int f_l_id, const int f_a_id, const int f_b_id, const int rhs_id)
{
   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;
         computeDrivingForce(time, *patch, temperature_id, phase_id, conc_id,
                             f_l_id, f_a_id, f_b_id, rhs_id);
      }
   }
}

void FreeEnergyStrategy::computeFreeEnergySolidB(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int fb_id, const bool gp)
{
   assert(temperature_id >= 0);
   assert(fb_id >= 0);

   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         computeFreeEnergySolidB(*patch, temperature_id, fb_id, gp);
      }
   }
}
