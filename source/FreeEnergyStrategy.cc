// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
#include "FreeEnergyStrategy.h"

#include "SAMRAI/pdat/CellData.h"

// default implementation: derivative is 0 (energy independent of c)
void FreeEnergyStrategy::computeDerivFreeEnergy(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const int df_id)
{
   assert(df_id >= 0);

   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
           ip++) {
         std::shared_ptr<hier::Patch> patch = *ip;

         std::shared_ptr<pdat::CellData<double> > df(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(df_id)));
         assert(df);

         df->fillAll(0.);
      }
   }
}

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
    const int temperature_id, const int phase_id, const int eta_id,
    const int conc_id, const int f_l_id, const int f_a_id, const int f_b_id,
    const int rhs_id)
{
   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;
         computeDrivingForce(time, *patch, temperature_id, phase_id, eta_id,
                             conc_id, f_l_id, f_a_id, f_b_id, rhs_id);
      }
   }
}

void FreeEnergyStrategy::computeFreeEnergySolidB(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int temperature_id, const int fb_id, const bool gp)
{
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
