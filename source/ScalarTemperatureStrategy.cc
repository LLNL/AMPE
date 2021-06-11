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
#include "ScalarTemperatureStrategy.h"

#include "SAMRAI/pdat/CellData.h"

ScalarTemperatureStrategy::ScalarTemperatureStrategy(
    const int temperature_id, const int temperature_scratch_id,
    const double temperature0, std::shared_ptr<tbox::Database> temperature_db)
{
   assert(temperature_id >= 0);
   assert(temperature_scratch_id >= 0);
   assert(temperature0 > 0.);

   d_temperature_id = temperature_id;
   d_temperature_scratch_id = temperature_scratch_id;
   d_temperature0 = temperature0;

   d_dtemperaturedt =
       temperature_db->getDoubleWithDefault("dtemperaturedt", 0.0);
   d_target_temperature =
       temperature_db->getDoubleWithDefault("target_temperature", 0.0);
}

double ScalarTemperatureStrategy::getCurrentMaxTemperature(
    std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
{
   (void)patch_hierarchy;

   return getCurrentTemperature(time);
}

double ScalarTemperatureStrategy::getCurrentMinTemperature(
    std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
{
   (void)patch_hierarchy;

   return getCurrentTemperature(time);
}

double ScalarTemperatureStrategy::getCurrentAverageTemperature(
    std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
{
   (void)patch_hierarchy;

   return getCurrentTemperature(time);
}

double ScalarTemperatureStrategy::getCurrentTemperature(const double time)
{
   double t = d_temperature0 + d_dtemperaturedt * time;

   // limit temperature to target if set.

   if (d_dtemperaturedt < 0. && t < d_target_temperature) {

      t = d_target_temperature;

   } else if (d_target_temperature > 0.0 && d_dtemperaturedt > 0. &&
              t > d_target_temperature) {

      t = d_target_temperature;
   }

   return t;
}

void ScalarTemperatureStrategy::setCurrentTemperature(
    std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
{
   assert(d_temperature_id >= 0);
   assert(d_temperature_scratch_id >= 0);

   if (d_dtemperaturedt != 0.0 || time == 0.0) {

      double t = getCurrentTemperature(time);

      int maxln = patch_hierarchy->getFinestLevelNumber();

      for (int ln = 0; ln <= maxln; ln++) {

         std::shared_ptr<hier::PatchLevel> level =
             patch_hierarchy->getPatchLevel(ln);

         for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
              p++) {

            std::shared_ptr<hier::Patch> patch = *p;

            std::shared_ptr<pdat::CellData<double> > t_data(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_temperature_id)));

            t_data->fillAll(t);

            std::shared_ptr<pdat::CellData<double> > ts_data(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_temperature_scratch_id)));

            ts_data->fillAll(t);
         }
      }
   }
}
