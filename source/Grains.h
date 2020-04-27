// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
#ifndef included_Grains
#define included_Grains

// Headers for basic SAMRAI objects
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/IEEE.h"

// Headers for major algorithm/data structure objects
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

using namespace SAMRAI;

class Grains
{
 public:
   Grains(const int qlen, const bool visit_output,
          std::shared_ptr<tbox::Database> input_db);
   virtual ~Grains() {}

   void initialize(std::shared_ptr<tbox::Database> input_db,
                   const bool all_periodic);

   void resetHierarchyConfiguration(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
       const int coarsest_level, const int finest_level);
   void initializeLevelData(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
       const int level_number, const double time, const bool can_be_refined,
       const bool initial_time,
       const std::shared_ptr<hier::PatchLevel>& old_level,
       const bool allocate_data);
   void initializeRefineCoarsenAlgorithms(
       std::shared_ptr<geom::CartesianGridGeometry> grid_geom,
       std::shared_ptr<hier::CoarsenOperator> quat_coarsen_op);

   void registerVariables(void);
   void registerWithVisit(std::shared_ptr<appu::VisItDataWriter>);

   void computeGrainDiagnostics(void);

   void findAndNumberGrains(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int phase_id, const int weight_id, const double time);

   void computeGrainVolumes(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int weight_id);

   void computeGrainConcentrations(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const double time, const int conc_id, const int weight_id);

   void extendGrainOrientation(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const double time, int quat_scratch_id, int phase_id, int quat_id);

   int getNumberOfGrains() const { return d_number_of_grains; }

 private:
   hier::ComponentSelector d_local_data;

   const int d_qlen;
   const bool d_visit_output;
   bool d_grain_diag_isActive;
   bool d_grain_extend_isActive;

   double d_grain_phase_threshold;
   double d_grain_size_minimum;

   void initializeRefineCoarsenAlgorithms();

   std::shared_ptr<pdat::CellVariable<int> > d_grain_number_var;
   int d_grain_number_id;
   int d_grain_number_scr_id;
   std::shared_ptr<pdat::CellVariable<double> > d_grain_volume_var;
   int d_grain_volume_id;

   std::shared_ptr<pdat::CellVariable<int> > d_grain_extend_var;
   int d_grain_extend_id;
   int d_grain_extend_scr_id;
   std::shared_ptr<pdat::CellVariable<double> > d_grain_quat_var;
   int d_grain_quat_id;
   int d_grain_quat_scr_id;

   std::shared_ptr<hier::RefineOperator> d_grain_number_refine_op;
   std::shared_ptr<hier::CoarsenOperator> d_grain_number_coarsen_op;
   std::shared_ptr<hier::RefineOperator> d_grain_extend_refine_op;
   std::shared_ptr<hier::CoarsenOperator> d_grain_extend_coarsen_op;
   std::shared_ptr<hier::RefineOperator> d_grain_quat_refine_op;
   std::shared_ptr<hier::CoarsenOperator> d_grain_quat_coarsen_op;

   std::shared_ptr<xfer::RefineAlgorithm> d_grain_number_refine_alg;
   std::shared_ptr<xfer::CoarsenAlgorithm> d_grain_number_coarsen_alg;
   std::shared_ptr<xfer::RefineAlgorithm> d_grain_extend_refine_alg;
   std::shared_ptr<xfer::CoarsenAlgorithm> d_grain_extend_coarsen_alg;
   std::shared_ptr<xfer::RefineAlgorithm> d_grain_quat_refine_alg;
   std::shared_ptr<xfer::CoarsenAlgorithm> d_grain_quat_coarsen_alg;

   std::vector<std::shared_ptr<xfer::RefineSchedule> >
       d_grain_number_refine_sched;
   std::vector<std::shared_ptr<xfer::CoarsenSchedule> >
       d_grain_number_coarsen_sched;
   std::vector<std::shared_ptr<xfer::RefineSchedule> >
       d_grain_extend_refine_sched;
   std::vector<std::shared_ptr<xfer::CoarsenSchedule> >
       d_grain_extend_coarsen_sched;
   std::vector<std::shared_ptr<xfer::RefineSchedule> >
       d_grain_quat_refine_sched;
   std::vector<std::shared_ptr<xfer::CoarsenSchedule> >
       d_grain_quat_coarsen_sched;

   int d_number_of_grains;

   // Timers
   std::shared_ptr<tbox::Timer> t_findAndNumberGrains_timer;
   std::shared_ptr<tbox::Timer> t_extendGrainOrientation_timer;
};

#endif
