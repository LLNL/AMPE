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
#ifndef included_GradientTemperatureStrategy
#define included_GradientTemperatureStrategy

#include "TemperatureStrategy.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

#include <cassert>


// This strategy fills the temperature field from a single scalar value
// which may depend on time.

class GradientTemperatureStrategy : public TemperatureStrategy
{
 public:
   GradientTemperatureStrategy(const int temperature_id,
                               const int temperature_scratch_id,
                               const double temperature0,
                               const double velocity_x,
                               std::shared_ptr<tbox::Database> temperature_db);

   ~GradientTemperatureStrategy(){};

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

   double d_gradient[NDIM];
   double d_center[3];
   double d_center_init[3];

   // moving frame velocity in x direction
   double d_velocity_x;
   ;

   void setCurrentTemperature(const double temperature, hier::Patch& patch,
                              std::shared_ptr<pdat::CellData<double> > cd_temp);
};

#endif
