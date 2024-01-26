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
#ifndef included_SteadyStateTemperatureStrategy
#define included_SteadyStateTemperatureStrategy

#include "TemperatureStrategy.h"
#include "HeatCapacityStrategy.h"

#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/tbox/Database.h"

#include <cassert>
#include <vector>

// This strategy fills the temperature field from the solution of
// diffusion equation
#include "TemperatureFACSolver.h"

class SteadyStateTemperatureStrategy : public TemperatureStrategy
{
 public:
   SteadyStateTemperatureStrategy(
       const int temperature_scratch_id,
       const int rhs_id,  // used internally only, but allocated outside class
       const int weight_id,
       std::shared_ptr<tbox::Database> temperature_sys_solver_database,
       solv::LocationIndexRobinBcCoefs* bc_coefs);

   virtual ~SteadyStateTemperatureStrategy();

   void initialize(
       const std::shared_ptr<hier::PatchHierarchy>& patch_hierarchy);

   virtual double getCurrentMinTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
   {
      (void)time;
      math::HierarchyCellDataOpsReal<double> cellops(patch_hierarchy);

      return cellops.min(d_temperature_scratch_id);
   }

   virtual double getCurrentMaxTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
   {
      (void)time;
      math::HierarchyCellDataOpsReal<double> cellops(patch_hierarchy);

      return cellops.max(d_temperature_scratch_id);
   }

   virtual double getCurrentAverageTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
   {
      (void)time;
      math::HierarchyCellDataOpsReal<double> cellops(patch_hierarchy);

      return cellops.integral(d_temperature_scratch_id, d_weight_id) /
             cellops.sumControlVolumes(d_temperature_scratch_id, d_weight_id);
   }

   virtual void setCurrentTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time) = 0;

   void resetSolversState(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy);

   void solveSystem()
   {
      d_temperature_sys_solver->solveSystem(d_temperature_scratch_id, d_rhs_id);
   }

 protected:
   TemperatureFACSolver* d_temperature_sys_solver;
   int d_temperature_scratch_id;
   int d_rhs_id;

 private:
   int d_weight_id;

   double getCurrentTemperature(const double time);
};


#endif
