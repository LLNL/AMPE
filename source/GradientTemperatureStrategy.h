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
#include "TemperatureHistory.h"

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
                               const double velocity_x,
                               std::shared_ptr<tbox::Database> temperature_db);

   ~GradientTemperatureStrategy(){};

   void resetVelocity(const double time, const double velocity_x)
   {
      // reset d_ref_center and d_ref_time first according
      // to velocity used between d_ref_time and time
      if (d_temperature_history)
         d_ref_center[0] += (time - d_ref_time) * d_velocity_x;
      else
         d_ref_center[0] -= (time - d_ref_time) * d_velocity_x;

      // now reset d_ref_time and d_velocity_x
      d_ref_time = time;
      d_velocity_x = velocity_x;
   }

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

   // center at first time step
   double d_init_center[NDIM];

   double d_gradient[NDIM];
   double d_center[NDIM];
   double d_ref_center[NDIM];

   // time associated with d_ref_center
   double d_ref_time;

   // moving frame velocity in x direction
   double d_velocity_x;

   std::shared_ptr<TemperatureHistory> d_temperature_history;

   void setCurrentTemperature(const double temperature, hier::Patch& patch,
                              std::shared_ptr<pdat::CellData<double> > cd_temp);
};

#endif
