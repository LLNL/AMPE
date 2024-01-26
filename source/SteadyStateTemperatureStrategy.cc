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
#include "SteadyStateTemperatureStrategy.h"

#include "TemperatureFACSolver.h"
#include "QuatFort.h"
#include "TemperatureFACOps.h"

#include "SAMRAI/tbox/Database.h"

SteadyStateTemperatureStrategy::SteadyStateTemperatureStrategy(
    const int temperature_scratch_id, const int rhs_id, const int weight_id,
    std::shared_ptr<tbox::Database> temperature_sys_solver_database,
    solv::LocationIndexRobinBcCoefs* bc_coefs)
{
   assert(temperature_scratch_id >= 0);
   assert(rhs_id >= 0);
   assert(weight_id >= 0);

   d_temperature_scratch_id = temperature_scratch_id;
   d_rhs_id = rhs_id;
   d_weight_id = weight_id;

   std::shared_ptr<TemperatureFACOps> fac_ops(
       new TemperatureFACOps("SteadyStateTemperatureStrategyFACOps",
                             temperature_sys_solver_database));

   d_temperature_sys_solver =
       new TemperatureFACSolver("SteadyStateTemperatureStrategySysSolver",
                                fac_ops, temperature_sys_solver_database);

   if (bc_coefs != 0) {
      d_temperature_sys_solver->setBcObject(bc_coefs);
   } else {
      d_temperature_sys_solver->setBoundaries("Dirichlet");
   }
}

//-----------------------------------------------------------------------

SteadyStateTemperatureStrategy::~SteadyStateTemperatureStrategy()
{
   delete d_temperature_sys_solver;
}

//-----------------------------------------------------------------------

void SteadyStateTemperatureStrategy::initialize(
    const std::shared_ptr<hier::PatchHierarchy>& patch_hierarchy)
{
   int finest = patch_hierarchy->getFinestLevelNumber();

   d_temperature_sys_solver->initializeSolverState(d_temperature_scratch_id,
                                                   d_rhs_id, patch_hierarchy, 0,
                                                   finest);
}

//-----------------------------------------------------------------------

void SteadyStateTemperatureStrategy::resetSolversState(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   assert(d_temperature_sys_solver);
   d_temperature_sys_solver->resetSolverState(d_temperature_scratch_id,
                                              d_rhs_id, hierarchy);
}
