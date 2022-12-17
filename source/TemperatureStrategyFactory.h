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
#ifndef included_TemperatureStrategyFactory
#define included_TemperatureStrategyFactory

#include "QuatModelParameters.h"
#include "HeatCapacityStrategy.h"

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

class TemperatureStrategy;
using namespace SAMRAI;

class TemperatureStrategyFactory
{
 public:
   TemperatureStrategyFactory(
       const int temperature_id, const int temperature_scratch_id,
       const int conc_id, const int weight_id, const int temperature_rhs_id,
       const int cp_id, const double molar_volume,
       const bool with_concentration,
       std::shared_ptr<geom::CartesianGridGeometry> grid_geometry,
       HeatCapacityStrategy* heat_capacity_strategy);

   std::shared_ptr<TemperatureStrategy> create(
       std::shared_ptr<tbox::Database> model_db,
       std::shared_ptr<tbox::Database> integrator_db,
       const QuatModelParameters&);

 private:
   const int d_temperature_id;
   const int d_temperature_scratch_id;
   const int d_conc_id;
   const int d_weight_id;
   const int d_temperature_rhs_id;
   const int d_cp_id;

   double d_molar_volume;
   bool d_with_concentration;
   solv::LocationIndexRobinBcCoefs* d_temperature_bc_coefs;
   std::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;

   HeatCapacityStrategy* d_heat_capacity_strategy;

   double readTemperature0(
       std::shared_ptr<tbox::Database> temperature_db,
       QuatModelParameters::TemperatureType temperature_type);
};

#endif
