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
#ifndef included_SteadyStateTemperatureCompositionSource
#define included_SteadyStateTemperatureCompositionSource

#include "SteadyStateTemperatureStrategy.h"
#include "HeatCapacityStrategy.h"

#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/tbox/Database.h"

#include <cassert>
#include <vector>

// This strategy fills the temperature field from teh solution of
// diffusion equation
class TemperatureFACSolver;

class SteadyStateTemperatureCompositionSource
    : public SteadyStateTemperatureStrategy
{
 public:
   SteadyStateTemperatureCompositionSource(
       const int temperature_scratch_id,
       const int composition_id,  // source depends on composition
       const int rhs_id,  // used internally only, but allocated outside class
       const int weight_id, const double thermal_diffusivity, const int cp_id,
       const std::vector<double>& T_source,
       std::shared_ptr<tbox::Database> temperature_sys_solver_database,
       HeatCapacityStrategy* heat_capacity_strategy,
       solv::LocationIndexRobinBcCoefs* bc_coefs);

   ~SteadyStateTemperatureCompositionSource(){};

   void setCurrentTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time);

 private:
   const double d_thermal_diffusivity;
   const std::vector<double> d_T_source;

   int d_composition_id;
   int d_cp_id;

   HeatCapacityStrategy* d_heat_capacity_strategy;
};


#endif
