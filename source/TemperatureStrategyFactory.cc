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
#include "TemperatureStrategyFactory.h"

#include "ScalarTemperatureStrategy.h"
#include "GaussianTemperatureStrategy.h"
#include "GradientTemperatureStrategy.h"
#include "ConstantTemperatureStrategy.h"
#include "SteadyStateTemperatureCompositionSource.h"
#include "SteadyStateTemperatureGaussianSource.h"
#include "tools.h"

#include <vector>

using namespace SAMRAI;

TemperatureStrategyFactory::TemperatureStrategyFactory(
    const int temperature_id, const int temperature_scratch_id,
    const int conc_id, const int weight_id, const int temperature_rhs_id,
    const int cp_id, const double molar_volume, const bool with_concentration,
    std::shared_ptr<geom::CartesianGridGeometry> grid_geometry,
    HeatCapacityStrategy* heat_capacity_strategy)
    : d_temperature_id(temperature_id),
      d_temperature_scratch_id(temperature_scratch_id),
      d_conc_id(conc_id),
      d_weight_id(weight_id),
      d_temperature_rhs_id(temperature_rhs_id),
      d_cp_id(cp_id),
      d_molar_volume(molar_volume),
      d_with_concentration(with_concentration),
      d_grid_geometry(grid_geometry),
      d_heat_capacity_strategy(heat_capacity_strategy)
{
   d_temperature_bc_coefs = 0;
}

double TemperatureStrategyFactory::readTemperature0(
    std::shared_ptr<tbox::Database> temperature_db,
    QuatModelParameters::TemperatureType temperature_type)
{
   assert(temperature_db);

   double temperature0 = -1.;
   if (temperature_type == QuatModelParameters::TemperatureType::SCALAR ||
       temperature_type == QuatModelParameters::TemperatureType::GAUSSIAN ||
       temperature_type == QuatModelParameters::TemperatureType::GRADIENT) {
      if (temperature_db->keyExists("temperature0")) {
         temperature0 = temperature_db->getDouble("temperature0");
         printDeprecated("temperature0", "temperature");
      } else if (temperature_db->keyExists("T_parameter")) {
         temperature0 = temperature_db->getDouble("T_parameter");
         printDeprecated("T_parameter", "temperature");
      } else if (temperature_db->keyExists("temperature")) {
         temperature0 = temperature_db->getDouble("temperature");
      } else if (temperature_db->keyExists("initial")) {
         temperature0 = temperature_db->getDouble("initial");
      } else {
         TBOX_ERROR(
             "Error in TemperatureStrategyFactory: temperature not specified");
      }
   }
   return temperature0;
}

std::shared_ptr<TemperatureStrategy> TemperatureStrategyFactory::create(
    std::shared_ptr<tbox::Database> model_db,
    std::shared_ptr<tbox::Database> integrator_db,
    const QuatModelParameters& model_parameters)
{
   std::shared_ptr<TemperatureStrategy> strategy;

   std::shared_ptr<tbox::Database> temperature_db;
   if (model_db->keyExists("Temperature")) {
      temperature_db = model_db->getDatabase("Temperature");
   } else {
      temperature_db = model_db;
   }

   if (model_db->keyExists("BoundaryConditions") &&
       model_parameters.with_steady_temperature()) {
      std::shared_ptr<tbox::Database> bc_db =
          model_db->getDatabase("BoundaryConditions");
      if (bc_db->keyExists("Temperature")) {
         std::shared_ptr<tbox::Database> temperature_bc_db =
             bc_db->getDatabase("Temperature");
         d_temperature_bc_coefs =
             new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
                                                 "TFactoryTemperatureBcCoefs",
                                                 temperature_bc_db);
      }
   }

   // create strategy
   if (model_parameters
           .isTemperatureUniform()) {  // uniform T across spatial domain
      const double temperature0 =
          readTemperature0(temperature_db,
                           QuatModelParameters::TemperatureType::SCALAR);

      strategy.reset(new ScalarTemperatureStrategy(d_temperature_id,
                                                   d_temperature_scratch_id,
                                                   temperature0,
                                                   temperature_db));
   } else if (model_parameters.isTemperatureGaussian()) {

      strategy.reset(
          new GaussianTemperatureStrategy(d_temperature_id,
                                          d_temperature_scratch_id, d_weight_id,
                                          temperature_db, d_grid_geometry));
   } else if (model_parameters.isTemperatureGradient()) {
      const double frame_velocity = model_parameters.movingVelocity();
      strategy.reset(new GradientTemperatureStrategy(d_temperature_id,
                                                     d_temperature_scratch_id,
                                                     frame_velocity,
                                                     temperature_db));
   } else if (model_parameters.with_heat_equation()) {

      if (model_parameters.with_steady_temperature()) {
         assert(d_temperature_scratch_id >= 0);
         assert(d_temperature_rhs_id >= 0);
         assert(d_temperature_bc_coefs != 0);
         assert(d_heat_capacity_strategy);

         std::shared_ptr<tbox::Database> temperature_sys_solver_database;
         if (integrator_db->isDatabase("TemperatureSysSolver")) {
            temperature_sys_solver_database =
                integrator_db->getDatabase("TemperatureSysSolver");
         }
         if (model_parameters.isHeatSourceCompositionDependent())
            strategy.reset(new SteadyStateTemperatureCompositionSource(
                d_temperature_scratch_id, d_conc_id, d_temperature_rhs_id,
                d_weight_id, model_parameters.thermal_diffusivity(), d_cp_id,
                model_parameters.T_source(), temperature_sys_solver_database,
                d_heat_capacity_strategy, d_temperature_bc_coefs));
         else {
            std::shared_ptr<tbox::Database> heat_source_db =
                temperature_db->getDatabase("HeatSource");
            strategy.reset(new SteadyStateTemperatureGaussianSource(
                d_temperature_scratch_id, d_temperature_rhs_id, d_weight_id,
                model_parameters.thermal_diffusivity(), d_cp_id, heat_source_db,
                temperature_sys_solver_database, d_grid_geometry,
                d_heat_capacity_strategy, d_temperature_bc_coefs));
         }
      } else {
         strategy.reset(
             new ConstantTemperatureStrategy(d_temperature_id, d_weight_id));
      }

   } else {
      strategy.reset(
          new ConstantTemperatureStrategy(d_temperature_id, d_weight_id));
   }

   return strategy;
}
