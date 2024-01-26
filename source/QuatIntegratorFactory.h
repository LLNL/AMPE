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
#ifndef QUATINTEGRATORFACTORY_H
#define QUATINTEGRATORFACTORY_H

#include "QuatIntegrator.h"

#include "SAMRAI/tbox/MemoryDatabase.h"

class QuatIntegratorFactory
{
 public:
   static QuatIntegrator* create(
       const std::string integrator_name, QuatModelParameters& model_parameters,
       QuatModel* quat_model,
       std::shared_ptr<geom::CartesianGridGeometry> grid_geometry,
       const int qlen, const int ncompositions,
       std::shared_ptr<tbox::Database> input_db, const bool use_warm_start,
       const bool symmetry_aware, const bool all_periodic)
   {
      hier::VariableDatabase* variable_db =
          hier::VariableDatabase::getDatabase();

      std::shared_ptr<hier::VariableContext> current =
          variable_db->getContext("CURRENT");
      std::shared_ptr<hier::VariableContext> scratch =
          variable_db->getContext("SCRATCH");

      QuatIntegrator* integrator = nullptr;

      std::shared_ptr<tbox::Database> integrator_db =
          input_db->getDatabase("Integrator");

      std::shared_ptr<tbox::Database> model_db =
          input_db->getDatabase("ModelParameters");
      std::shared_ptr<tbox::Database> bc_db;
      if (!all_periodic) {
         bc_db = model_db->getDatabase("BoundaryConditions");
      }

      if (integrator_name == "quatonly") {
         integrator = new QuatIntegrator("QuatOnlyIntegrator", model_parameters,
                                         quat_model, current, scratch, qlen,
                                         0,  // 0 composition fields
                                         input_db, grid_geometry, bc_db,
                                         false,  // without phase
                                         false,  // without concentration
                                         false, false, false, false, false,
                                         false, symmetry_aware, false);
      } else {
         integrator = new QuatIntegrator(
             "Integrator", model_parameters, quat_model, current, scratch, qlen,
             ncompositions, input_db, grid_geometry, bc_db,
             model_parameters.with_phase(),
             model_parameters.with_concentration(),
             model_parameters.with_heat_equation(),
             model_parameters.with_steady_temperature(),
             model_parameters.with_gradT(),
             model_parameters.with_antitrapping(),
             model_parameters.with_partition_coeff(), use_warm_start,
             symmetry_aware,
             model_parameters
                 .use_diffs_to_compute_flux());  // use_gradq_for_flux
      }

      return integrator;
   }
};

#endif
