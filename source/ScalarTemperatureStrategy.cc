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
#include "ScalarTemperatureStrategy.h"

#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

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
    std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time)
{
   assert(d_temperature_id >= 0);
   assert(d_temperature_scratch_id >= 0);

   if (d_dtemperaturedt != 0.0 || time == 0.0) {

      double t = getCurrentTemperature(time);

      math::HierarchyCellDataOpsReal<double> mathops(hierarchy);
      mathops.setToScalar(d_temperature_id, t);
      mathops.setToScalar(d_temperature_scratch_id, t, false);
   }
}
