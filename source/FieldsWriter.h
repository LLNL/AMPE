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
#include "QuatModelParameters.h"
#include "QuatRefinePatchStrategy.h"

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

using namespace SAMRAI;


class FieldsWriter
{
 public:
   FieldsWriter(QuatModelParameters& model_parameters, std::string filename,
                int initial_conditions_level,
                std::shared_ptr<geom::CartesianGridGeometry> grid_geometry,
                int phase_id, int phase_scratch_id, int temperature_id,
                int temperature_scratch_id, int quat_id, int quat_scratch_id,
                int conc_id, int conc_scratch_id, int eta_id,
                int eta_scratch_id, const int ncompositions, const int qlen,
                QuatRefinePatchStrategy* all_refine_patch_strategy);

   void writeInitialConditionsFile(
       const std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time);

 private:
   QuatModelParameters& d_model_parameters;
   std::string d_filename;
   int d_initial_conditions_level;
   std::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;
   int d_phase_id;
   int d_phase_scratch_id;
   int d_temperature_id;
   int d_temperature_scratch_id;
   int d_quat_id;
   int d_quat_scratch_id;
   int d_conc_id;
   int d_conc_scratch_id;
   int d_eta_id;
   int d_eta_scratch_id;
   int d_ncompositions;
   int d_qlen;

   QuatRefinePatchStrategy* d_all_refine_patch_strategy;

   std::shared_ptr<hier::PatchLevel> FlattenHierarchy(
       const std::shared_ptr<hier::PatchHierarchy> src_hierarchy,
       const int level_number, const double time);
};
