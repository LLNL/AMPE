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
#ifndef included_ScalarTemperatureStrategy
#define included_ScalarTemperatureStrategy

#include "TemperatureStrategy.h"
#include <cassert>

// This strategy fills the temperature field from a single scalar value
// which may depend on time.

class ScalarTemperatureStrategy : public TemperatureStrategy
{
 public:
   ScalarTemperatureStrategy(const int temperature_id,
                             const int temperature_scratch_id,
                             const double temperature0,
                             std::shared_ptr<tbox::Database> temperature_db);

   ~ScalarTemperatureStrategy(){};

   virtual double getCurrentMaxTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time);
   virtual double getCurrentMinTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time);
   virtual double getCurrentAverageTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time);

   virtual void setCurrentTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time);

 private:
   int d_temperature_id;
   int d_temperature_scratch_id;
   int d_weight_id;

   double d_temperature0;
   double d_dtemperaturedt;
   double d_target_temperature;

   double getCurrentTemperature(const double time);
};

#endif
