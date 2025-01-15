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

ConcFreeEnergyStrategy::ConcFreeEnergyStrategy(const int conc_l_id,
                                               const int conc_a_id,
                                               const int conc_b_id)
    : d_conc_l_id(conc_l_id), d_conc_a_id(conc_a_id), d_conc_b_id(conc_b_id)
{
}

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

void ConcFreeEnergyStrategy::computeDerivFreeEnergyLiquid(
    hier::Patch& patch, const int temperature_id, const int df_id)
{
   assert(temperature_id >= 0);
   assert(df_id >= 0);

   computeDerivFreeEnergy(patch, temperature_id, df_id, d_conc_l_id,
                          Thermo4PFM::PhaseIndex::phaseL);
}

void ConcFreeEnergyStrategy::computeDerivFreeEnergySolidA(
    hier::Patch& patch, const int temperature_id, const int df_id)
{
   assert(temperature_id >= 0);
   assert(df_id >= 0);

   computeDerivFreeEnergy(patch, temperature_id, df_id, d_conc_a_id,
                          Thermo4PFM::PhaseIndex::phaseA);
}

void ConcFreeEnergyStrategy::computeDerivFreeEnergySolidB(
    hier::Patch& patch, const int temperature_id, const int df_id)
{
   assert(temperature_id >= 0);
   assert(df_id >= 0);

   computeDerivFreeEnergy(patch, temperature_id, df_id, d_conc_b_id,
                          Thermo4PFM::PhaseIndex::phaseB);
}

void ConcFreeEnergyStrategy::computeFreeEnergyLiquid(hier::Patch& patch,
                                                     const int temperature_id,
                                                     const int fl_id,
                                                     const bool gp)
{
   assert(temperature_id >= 0);
   assert(fl_id >= 0);

   assert(d_conc_l_id >= 0);

   computeFreeEnergy(patch, temperature_id, fl_id, d_conc_l_id,
                     Thermo4PFM::PhaseIndex::phaseL, gp);
}

void ConcFreeEnergyStrategy::computeFreeEnergySolidA(hier::Patch& patch,
                                                     const int temperature_id,
                                                     const int fa_id,
                                                     const bool gp)
{
   assert(temperature_id >= 0.);
   assert(fa_id >= 0);

   computeFreeEnergy(patch, temperature_id, fa_id, d_conc_a_id,
                     Thermo4PFM::PhaseIndex::phaseA, gp);
}

//=======================================================================

void ConcFreeEnergyStrategy::computeFreeEnergySolidB(hier::Patch& patch,
                                                     const int temperature_id,
                                                     const int fb_id,
                                                     const bool gp)
{
   assert(temperature_id >= 0.);
   assert(fb_id >= 0);

   computeFreeEnergy(patch, temperature_id, fb_id, d_conc_b_id,
                     Thermo4PFM::PhaseIndex::phaseB, gp);
}
