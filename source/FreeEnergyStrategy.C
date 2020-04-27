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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
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
