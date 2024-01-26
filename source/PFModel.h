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
#ifndef included_PFModel
#define included_PFModel

#include "EventInterval.h"
#include "Verbosity.h"

// Headers for basic SAMRAI objects
#include "SAMRAI/tbox/MemoryDatabase.h"

// Headers for major algorithm/data structure objects
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/hier/PatchHierarchy.h"

using namespace SAMRAI;

class PFModel : public tbox::Serializable,
                public mesh::StandardTagAndInitStrategy
{
 public:
   PFModel(void);
   virtual ~PFModel();

   virtual void Initialize(std::shared_ptr<tbox::MemoryDatabase>& input_db,
                           const std::string& run_name,
                           const bool is_from_restart,
                           const std::string& restart_read_dirname,
                           const int restore_num);

   virtual void readInitialDatabase(
       std::shared_ptr<tbox::Database> main_input_db);

   void setVerbosity(Verbosity* v) { d_verbosity = v; }

   virtual void getFromInput(std::shared_ptr<tbox::Database> input_db);

   virtual void getFromRestart(std::shared_ptr<tbox::Database> input_db);

   virtual void setupInitialDataLevel(void);

   virtual void setupHierarchy(void);

   virtual void Run(void);

   // pure virtual functions to be implemented by model
   virtual void CreateIntegrator(std::shared_ptr<tbox::Database> input_db) = 0;

   virtual void RegisterVariables(void) = 0;

   virtual void InitializeIntegrator(void) = 0;

   virtual void initializeCoarseRefineOperators(void) = 0;

   virtual void RegisterWithVisit(void) = 0;

   virtual double Advance(void) = 0;

   virtual void postAdvanceDiagnostics(void);

   virtual void preRunDiagnostics(void);

   virtual void postRunDiagnostics(void);

   virtual void writeRestartFile(void);

   virtual void Regrid(const std::shared_ptr<hier::PatchHierarchy> hierarchy);

   virtual void computeGrainDiagnostics(void);

   virtual void WriteInitialConditionsFile(
       std::string initial_conditions_file_name,
       int initial_conditions_level) = 0;

   // deallocate some temporary data to free some memory,
   // for example before some high memory footprint postprocessing
   virtual void DeallocateIntermediateLocalPatchData(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy)
   {
      (void)hierarchy;
   };

   //-----------------------------------------------------------------------
   //
   // Methods inherited from Serializable
   //
   virtual void putToRestart(const std::shared_ptr<tbox::Database>& db) const;

   //-----------------------------------------------------------------------
   //
   // Methods inherited from StandardTagAndInitStrategy
   //
   virtual void initializeLevelData(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
       const int level_number, const double time, const bool can_be_refined,
       const bool initial_time,
       const std::shared_ptr<hier::PatchLevel>& old_level,
       const bool allocate_data) = 0;

   virtual void resetHierarchyConfiguration(
       const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
       const int coarsest_level, const int finest_level) = 0;

   //-----------------------------------------------------------------------

 protected:
   std::shared_ptr<hier::PatchHierarchy> d_patch_hierarchy;
   std::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;
   std::shared_ptr<mesh::GriddingAlgorithm> d_gridding_algorithm;
   std::shared_ptr<appu::VisItDataWriter> d_visit_data_writer;

   bool d_amr_enabled;

   std::vector<int> d_tag_buffer_array;
   bool d_all_periodic;
   bool d_periodic_flag[NDIM];

   double d_time;
   double d_previous_timestep;

   std::shared_ptr<EventInterval> d_time_info_interval;

   std::shared_ptr<EventInterval> d_restart_interval;
   std::string d_restart_write_dirname;

   bool d_write_initial_conditions_file;
   std::string d_initial_conditions_file_name;
   int d_initial_conditions_level;

   std::shared_ptr<EventInterval> d_grain_diag_interval;

   // Initialization variables
   std::string d_init_data_filename;
   int d_level_of_init_data;
   int d_slice_index;
   hier::IntVector d_ratio_of_init_to_coarsest;
   std::shared_ptr<hier::PatchLevel> d_initial_level;
   std::vector<float> d_init_q;
   std::vector<float> d_init_c;
   std::vector<float> d_init_cl;
   std::vector<float> d_init_ca;
   std::vector<float> d_init_cb;
   float d_init_t;

   int d_max_cycles;
   double d_start_time;
   double d_end_time;
   int d_tag_buffer;
   int d_cycle;

   std::shared_ptr<EventInterval> d_regrid_interval;

   std::string d_object_name;

   std::shared_ptr<tbox::Timer> t_run_time;

 private:
   std::shared_ptr<EventInterval> d_visit_dump_interval;

   Verbosity* d_verbosity;
};

#endif
