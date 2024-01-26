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
#ifndef included_SteadyStateTemperatureGaussianSource
#define included_SteadyStateTemperatureGaussianSource

#include "SteadyStateTemperatureStrategy.h"
#include "HeatCapacityStrategy.h"
#include "ConcFort.h"

#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

#include <cassert>
#include <vector>

// This strategy fills the temperature field from teh solution of
// diffusion equation
class TemperatureFACSolver;

class SteadyStateTemperatureGaussianSource
    : public SteadyStateTemperatureStrategy
{
 public:
   SteadyStateTemperatureGaussianSource(
       const int temperature_scratch_id,
       const int rhs_id,  // used internally only, but allocated outside class
       const int weight_id, const double thermal_diffusivity, const int cp_id,
       std::shared_ptr<tbox::Database> heat_source_db,
       std::shared_ptr<tbox::Database> temperature_sys_solver_database,
       std::shared_ptr<geom::CartesianGridGeometry> grid_geometry,
       HeatCapacityStrategy* heat_capacity_strategy,
       solv::LocationIndexRobinBcCoefs* bc_coefs);

   ~SteadyStateTemperatureGaussianSource(){};

   void setCurrentTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time);

 private:
   const double d_thermal_diffusivity;

   int d_cp_id;

   HeatCapacityStrategy* d_heat_capacity_strategy;

   // Gaussian parameters
   double d_center[3];

   // Gaussian center at t=0
   double d_center0[3];
   double d_periodic_length[3];

   // velocity of Gaussian center
   double d_velocity[3];

   // Gaussian in time for pulse
   double d_pulse_time;
   double d_pulse_width;

   std::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;
   int d_periodic_flag[3];

   double d_source_peak;
   double d_standard_dev;

   bool d_verbose;
};


#endif
