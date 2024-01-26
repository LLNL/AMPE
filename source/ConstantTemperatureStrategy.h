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
#ifndef included_ConstantTemperatureStrategy
#define included_ConstantTemperatureStrategy

#include "TemperatureStrategy.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

// This strategy leaves the temperature field alone, assuming it will
// be initialized from an initial conditions file and left constant
// thereafter, or set through a heat equation time evolution

class ConstantTemperatureStrategy : public TemperatureStrategy
{
 public:
   ConstantTemperatureStrategy(const int temperature_id, const int weight_id)
   {
      assert(temperature_id >= 0);
      assert(weight_id >= 0);

      d_temperature_id = temperature_id;
      d_weight_id = weight_id;
   }

   ~ConstantTemperatureStrategy(){};

   virtual double getCurrentMinTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
   {
      (void)time;
      math::HierarchyCellDataOpsReal<double> cellops(patch_hierarchy);

      return cellops.min(d_temperature_id);
   }

   virtual double getCurrentMaxTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
   {
      (void)time;
      math::HierarchyCellDataOpsReal<double> cellops(patch_hierarchy);

      return cellops.max(d_temperature_id);
   }

   virtual double getCurrentAverageTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
   {
      (void)time;
      math::HierarchyCellDataOpsReal<double> cellops(patch_hierarchy);

      return cellops.integral(d_temperature_id, d_weight_id) /
             cellops.sumControlVolumes(d_temperature_id, d_weight_id);
   }


   virtual void setCurrentTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
   {
      (void)patch_hierarchy;
      (void)time;
   }

 private:
   int d_temperature_id;
   int d_weight_id;
};

#endif
