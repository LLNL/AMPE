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
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/hier/VariableDatabase.h"

#include <cassert>


#include "PFModel.h"

PFModel::PFModel() : d_ratio_of_init_to_coarsest(tbox::Dimension(NDIM))
{
   d_object_name = "PFModel";

   d_verbosity = NULL;

   d_slice_index = -1;

   d_init_t = -1.;

   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);

   tbox::TimerManager* time_man = tbox::TimerManager::getManager();
   t_run_time = time_man->getTimer("AMPE::PFModel::run");
}

//=======================================================================

PFModel::~PFModel()
{
   d_patch_hierarchy.reset();
   d_grid_geometry.reset();
   d_gridding_algorithm.reset();
   d_visit_data_writer.reset();
}

//=======================================================================

void PFModel::Initialize(std::shared_ptr<tbox::MemoryDatabase>& input_db,
                         const std::string& run_name,
                         const bool is_from_restart,
                         const std::string& restart_read_dirname,
                         const int restore_num)
{
   //-----------------------------------------------------------------------
   tbox::plog << "PFModel::Initialize()" << std::endl;

   tbox::Dimension dim(NDIM);

   if (input_db->isDatabase("Verbosity")) {
      std::shared_ptr<tbox::Database> verbosity_db =
          input_db->getDatabase("Verbosity");

      if (verbosity_db->getBoolWithDefault("silent", false)) {
         d_verbosity->setBasicLevel(0);
      } else {
         d_verbosity->setBasicLevel(
             verbosity_db->getIntegerWithDefault("level", 1));
      }

      d_verbosity->setAmrLevel(
          verbosity_db->getIntegerWithDefault("amr_level", 1));
   }

   if (d_verbosity->notSilent()) {
      tbox::pout << "Run name = " << run_name << std::endl;
   }

   //-----------------------------------------------------------------------
   //
   // Restart
   //
   d_restart_write_dirname = "r." + run_name;

   d_restart_interval.reset(
       new EventInterval(input_db, "Restart", 0.0, "step", false, true));

   if (input_db->isDatabase("Restart")) {
      std::shared_ptr<tbox::Database> restart_db =
          input_db->getDatabase("Restart");

      if (restart_db->keyExists("dirname")) {
         d_restart_write_dirname = restart_db->getString("dirname");
      }
   }

   /*
    * Get restart manager and root restart database.  If run is from
    * restart, open the restart file.
    */

   tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   if (is_from_restart) {
      restart_manager->openRestartFile(restart_read_dirname, restore_num,
                                       mpi.getSize());

      getFromRestart(input_db);
   } else {
      getFromInput(input_db);

      d_time = d_start_time;
      d_previous_timestep = 0.0;
      d_cycle = 0;
   }

   //-----------------------------------------------------------------------
   //
   // Time step info
   //
   d_time_info_interval.reset(
       new EventInterval(input_db, "TimestepOutput", 1, "step"));

   //-----------------------------------------------------------------------
   //
   // Grain Diagnostics
   //
   d_grain_diag_interval.reset(
       new EventInterval(input_db, "GrainDiagnostics", 0.0, "step"));

   //-----------------------------------------------------------------------
   //
   // Writing initial conditions file
   //
   d_write_initial_conditions_file = false;
   d_initial_conditions_file_name = "";
   d_initial_conditions_level = 0;

   if (input_db->isDatabase("InitialConditions")) {
      std::shared_ptr<tbox::Database> ic_db =
          input_db->getDatabase("InitialConditions");

      if (ic_db->isDatabase("WriteEndingFile")) {
         std::shared_ptr<tbox::Database> end_ic_db =
             ic_db->getDatabase("WriteEndingFile");

         d_write_initial_conditions_file = true;
         d_initial_conditions_file_name = end_ic_db->getString("filename");

         if (end_ic_db->keyExists("amr_level_of_data")) {
            d_initial_conditions_level =
                end_ic_db->getInteger("amr_level_of_data");
         }
      }
   }

   //-----------------------------------------------------------------------
   //
   // SAMRAI
   //
   //-----------------------------------------------------------------------

   std::shared_ptr<tbox::Database> cart_db(
       new tbox::MemoryDatabase("CartesianGeometry"));
   std::shared_ptr<tbox::Database> geo_db(input_db->getDatabase("Geometry"));


   int periodic_array[NDIM];
   if (geo_db->keyExists("periodic_dimension"))
      geo_db->getIntegerArray("periodic_dimension", periodic_array, NDIM);
   else
      for (int ii = 0; ii < NDIM; ii++) {
         periodic_array[ii] = 1;
      }

   cart_db->putIntegerArray("periodic_dimension", periodic_array, NDIM);

   int tmp_int_array[NDIM];
   geo_db->getIntegerArray("coarsest_level_resolution", tmp_int_array, NDIM);

   // Possibly assert that tmp_int_array has powers of 2, as long as that is
   // a hypre restriction.

   int lo_array[NDIM];
   int up_array[NDIM];

   for (int ii = 0; ii < NDIM; ii++) {
      lo_array[ii] = 0;
      up_array[ii] = tmp_int_array[ii] - 1;
   }
   for (int ii = 0; ii < NDIM; ii++) {
      if (up_array[ii] == lo_array[ii] && periodic_array[ii] == 1)
         TBOX_ERROR("Periodic BC require more than one cell along axis "
                    << ii << std::endl);
   }

   tbox::DatabaseBox domain_box(tbox::Dimension(NDIM), lo_array, up_array);

   cart_db->putDatabaseBox("domain_boxes", domain_box);

   std::vector<double> x_lo = geo_db->getDoubleVector("x_lo");
   std::vector<double> x_up = geo_db->getDoubleVector("x_up");

   cart_db->putDoubleVector("x_lo", x_lo);
   cart_db->putDoubleVector("x_up", x_up);

   d_grid_geometry.reset(new geom::CartesianGridGeometry(tbox::Dimension(NDIM),
                                                         "CartesianGeometry",
                                                         cart_db));


   //-----------------------------------------------------------------------

   std::shared_ptr<tbox::Database> tag_db(
       new tbox::MemoryDatabase("StandardTagAndInitialize"));

   if (d_amr_enabled) {
      std::shared_ptr<tbox::Database> amr_db(input_db->getDatabase("Amr"));

      d_regrid_interval.reset(new EventInterval(amr_db, "Regrid", 5, "step"));

      if (amr_db->isDatabase("StandardTagAndInitialize")) {
         tag_db = amr_db->getDatabase("StandardTagAndInitialize");
      } else {
         tag_db->putString("tagging_method", "GRADIENT_DETECTOR");
      }
   } else {
      d_regrid_interval.reset(
          new EventInterval(std::shared_ptr<tbox::Database>(), "", 0));
   }

   //-----------------------------------------------------------------------

   d_visit_dump_interval.reset(
       new EventInterval(input_db, "Visit", 0.0, "step", true, true));

   //-----------------------------------------------------------------------

   // required: max_levels, largest_patch_size, ratio_to_coarser

   std::shared_ptr<tbox::Database> grid_db(
       new tbox::MemoryDatabase("StandardTagAndInitialize"));

   if (d_amr_enabled) {
      std::shared_ptr<tbox::Database> amr_db(input_db->getDatabase("Amr"));

      grid_db = amr_db->getDatabase("GriddingAlgorithm");
   } else {
      grid_db->putInteger("max_levels", 1);

      std::shared_ptr<tbox::Database> largest_patch_db(
          grid_db->putDatabase("largest_patch_size"));

      int lg_patch_size[NDIM];
      for (int ii = 0; ii < NDIM; ii++) {
         lg_patch_size[ii] = 1000000;
      }

      largest_patch_db->putIntegerArray("level_0", lg_patch_size, NDIM);

      // This database needs to be created, but it can be left empty
      grid_db->putDatabase("ratio_to_coarser");
   }

   std::shared_ptr<mesh::StandardTagAndInitialize> error_detector(
       new mesh::StandardTagAndInitialize("StandardTagAndInitialize", this,
                                          tag_db));

   std::shared_ptr<mesh::BergerRigoutsos> box_generator(
       new mesh::BergerRigoutsos(dim));

   std::shared_ptr<mesh::TreeLoadBalancer> load_balancer(
       new mesh::TreeLoadBalancer(dim, "LoadBalancer"));
   load_balancer->setSAMRAI_MPI(SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld());

   d_patch_hierarchy.reset(
       new hier::PatchHierarchy("PatchHierarchy", d_grid_geometry,
                                input_db->isDatabase("PatchHierarchy")
                                    ? input_db->getDatabase("PatchHierarchy")
                                    : std::shared_ptr<tbox::Database>()));

   d_gridding_algorithm.reset(
       new mesh::GriddingAlgorithm(d_patch_hierarchy, "GriddingAlgorithm",
                                   grid_db, error_detector, box_generator,
                                   load_balancer));

   // Check for periodic coordinates

   const hier::IntVector& one_vec = hier::IntVector::getOne(dim);
   const hier::IntVector periodic = d_grid_geometry->getPeriodicShift(one_vec);

   d_all_periodic = true;
   for (int dd = 0; dd < NDIM; dd++) {
      d_periodic_flag[dd] = periodic[dd] != 0;
      d_all_periodic = d_all_periodic && d_periodic_flag[dd];
   }

   CreateIntegrator(input_db);

   // Both the PFModel and the integrator will define their variables here.
   //
   RegisterVariables();

   // initialize Refine operators and algorithms, coarsen, etc.
   initializeCoarseRefineOperators();

   //-----------------------------------------------------------------------
   //
   // Visit
   //
   if (d_visit_dump_interval->isActive()) {

      std::shared_ptr<tbox::Database> visit_db(input_db->getDatabase("Visit"));

      std::string visit_dump_dirname = "v." + run_name;
      if (visit_db->keyExists("dirname")) {
         visit_dump_dirname = visit_db->getString("dirname");
      }

      int visit_number_procs_per_file =
          visit_db->getIntegerWithDefault("number_procs_per_file", 1);

      d_visit_data_writer.reset(
          new appu::VisItDataWriter(dim, "PFM VisIt Writer", visit_dump_dirname,
                                    visit_number_procs_per_file));
      assert(d_visit_data_writer);

      RegisterWithVisit();
   }

   d_tag_buffer_array =
       std::vector<int>(d_patch_hierarchy->getMaxNumberOfLevels());
}

//=======================================================================

void PFModel::readInitialDatabase(std::shared_ptr<tbox::Database> input_db)
{
   // InitialConditions {
   //    filename = "initial_data.nc"
   //    amr_level_of_data = 0
   //    slice_index = 0
   // }

   std::shared_ptr<tbox::Database> data_db(
       input_db->getDatabase("InitialConditions"));
   if (!data_db) {
      TBOX_ERROR("Input error: no InitialConditions database" << std::endl);
   }

   if (data_db->keyExists("filename")) {
      d_init_data_filename = data_db->getString("filename");
   }
   d_level_of_init_data = 0;
   if (data_db->keyExists("amr_level_of_data")) {
      d_level_of_init_data = data_db->getInteger("amr_level_of_data");
   }

   // Read value for uniform quat if available
   if (data_db->keyExists("init_q")) {
      size_t n = data_db->getArraySize("init_q");
      d_init_q.resize(n);
      data_db->getFloatArray("init_q", &d_init_q[0], n);
   }
   // Read value for uniform concentration if available
   if (data_db->keyExists("init_c")) {
      size_t n = data_db->getArraySize("init_c");
      d_init_c.resize(n);
      data_db->getFloatArray("init_c", &d_init_c[0], n);
   }
   if (data_db->keyExists("init_cl")) {
      size_t n = data_db->getArraySize("init_cl");
      d_init_cl.resize(n);
      data_db->getFloatArray("init_cl", &d_init_cl[0], n);
   }
   if (data_db->keyExists("init_ca")) {
      size_t n = data_db->getArraySize("init_ca");
      d_init_ca.resize(n);
      data_db->getFloatArray("init_ca", &d_init_ca[0], n);
   }
   if (data_db->keyExists("init_cb")) {
      size_t n = data_db->getArraySize("init_cb");
      d_init_cb.resize(n);
      data_db->getFloatArray("init_cb", &d_init_cb[0], n);
   }

   // Read value for uniform temperature if available
   if (data_db->keyExists("init_t")) {
      d_init_t = data_db->getFloat("init_t");
   }
#if (NDIM == 2)
   if (data_db->keyExists("slice_index")) {
      d_slice_index = data_db->getInteger("slice_index");
   }
#endif

   hier::IntVector ratio(tbox::Dimension(NDIM), 1);
   for (int l = d_level_of_init_data; l > 0; l--)
      ratio *= d_patch_hierarchy->getRatioToCoarserLevel(l);
   d_ratio_of_init_to_coarsest = ratio;
}

//=======================================================================

void PFModel::setupInitialDataLevel(void)
{
   assert(d_grid_geometry);
   return;
}

//=======================================================================
// Setup hierarchy

void PFModel::setupHierarchy(void)
{
   tbox::plog << "\nPFModel::setupHierarchy..." << std::endl;
   double init_time = 0.;
   d_gridding_algorithm->makeCoarsestLevel(init_time);

   bool done = false;
   bool initial_time = true;
   for (int ln = 0; d_patch_hierarchy->levelCanBeRefined(ln) && !done; ln++) {
      d_gridding_algorithm->makeFinerLevel(d_tag_buffer_array[ln], initial_time,
                                           d_cycle, init_time);
      done = !(d_patch_hierarchy->finerLevelExists(ln));
   }
}

//=======================================================================

void PFModel::Run(void)
{
   t_run_time->start();

   preRunDiagnostics();

   if (d_visit_dump_interval->includeInitial(d_time)) {
      assert(d_visit_data_writer);

      d_visit_data_writer->writePlotData(d_patch_hierarchy, d_cycle, d_time);
   }

   if (d_verbosity->notSilent() &&
       d_time_info_interval->includeInitial(d_time)) {
      tbox::pout << "cycle # " << d_cycle << " : t = " << d_time << std::endl;
   }

   if (d_regrid_interval->includeInitial(d_time)) {
      writeRestartFile();
   }

   if (d_grain_diag_interval->includeInitial(d_time)) {
      computeGrainDiagnostics();
   }

   bool flag_visit_dump = true;
   while ((d_time < d_end_time) && (d_cycle < d_max_cycles)) {

      d_cycle++;

      const double dt = Advance();

      d_previous_timestep = dt;

      if (d_verbosity->notSilent() &&
          d_time_info_interval->hasIntervalPassed(d_cycle, d_time)) {

         tbox::pout << "cycle # " << d_cycle << " : t = " << d_time
                    << " : dt = " << dt << std::endl;
      }

      if (d_restart_interval->hasIntervalPassed(d_cycle, d_time)) {

         writeRestartFile();
      }

      if (d_grain_diag_interval->hasIntervalPassed(d_cycle, d_time)) {

         computeGrainDiagnostics();
      }

      postAdvanceDiagnostics();

      if (d_visit_dump_interval->hasIntervalPassed(d_cycle, d_time)) {

         assert(d_visit_data_writer);

         d_visit_data_writer->writePlotData(d_patch_hierarchy, d_cycle, d_time);
         flag_visit_dump = true;
      } else {
         flag_visit_dump = false;
      }

      if (d_patch_hierarchy->getMaxNumberOfLevels() > 1 &&
          d_regrid_interval->hasIntervalPassed(d_cycle, d_time)) {

         Regrid(d_patch_hierarchy);
      }
   }

   // Final diagnostics

   postRunDiagnostics();

   if (d_visit_dump_interval->includeFinal(d_time) && !flag_visit_dump) {
      assert(d_visit_data_writer);

      d_visit_data_writer->writePlotData(d_patch_hierarchy, d_cycle, d_time);
   }

   if (d_restart_interval->includeFinal(d_time)) {
      writeRestartFile();
   }

   DeallocateIntermediateLocalPatchData(d_patch_hierarchy);

   if (d_write_initial_conditions_file) {
      WriteInitialConditionsFile(d_initial_conditions_file_name,
                                 d_initial_conditions_level);
      if (d_verbosity->notSilent()) {
         tbox::pout << "Wrote initial conditions file "
                    << d_initial_conditions_file_name << std::endl;
      }
   }

   if (d_verbosity->notSilent()) {
      tbox::pout << "Run complete" << std::endl;
   }

   t_run_time->stop();
}

//=======================================================================

void PFModel::postAdvanceDiagnostics(void) {}

//=======================================================================

void PFModel::preRunDiagnostics(void) {}

//=======================================================================

void PFModel::postRunDiagnostics(void) {}

//=======================================================================

void PFModel::writeRestartFile(void)
{
   tbox::RestartManager::getManager()->writeRestartFile(d_restart_write_dirname,
                                                        d_cycle);
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   mpi.Barrier();
   if (d_verbosity->notSilent()) {
      tbox::pout << "Wrote restart file" << std::endl;
   }
}

//=======================================================================

void PFModel::Regrid(const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   int fine_lev_before = d_patch_hierarchy->getFinestLevelNumber();

   d_gridding_algorithm->regridAllFinerLevels(0, d_tag_buffer_array, d_cycle,
                                              d_time);

   int fine_lev_after = d_patch_hierarchy->getFinestLevelNumber();

   if (d_verbosity->amrLevel() > 0) {
      if (d_verbosity->amrLevel() > 1) {
         tbox::pout << "############################################\n"
                    << "Amr hierarchy regridded\n"
                    << "   Finest level before regrid: " << fine_lev_before
                    << "\n"
                    << "   Finest level after regrid: " << fine_lev_after
                    << "\n"
                    << "############################################"
                    << std::endl;
      } else {
         tbox::pout << "AMR hierarchy regridded" << std::endl;
      }
   }
}

//=======================================================================

void PFModel::computeGrainDiagnostics(void)
{
   // This will be a virtual override when implemented.  Derived
   // classes should not call this base class function.

   if (d_verbosity->notSilent()) {
      tbox::pout << "Grain diagnostics unimplemented" << std::endl;
   }
}

//=======================================================================
//
// Methods inherited from Serializable
//

void PFModel::putToRestart(const std::shared_ptr<tbox::Database>& db) const
{
   assert(db);

   db->putInteger("max_cycles", d_max_cycles);
   db->putDouble("start_time", d_start_time);
   db->putDouble("end_time", d_end_time);
   db->putDouble("time", d_time);
   db->putInteger("cycle", d_cycle);
   db->putDouble("previous_timestep", d_previous_timestep);
   db->putBool("amr_enabled", d_amr_enabled);
}


void PFModel::getFromInput(std::shared_ptr<tbox::Database> input_db)
{
   assert(input_db);

   if (input_db->keyExists("max_cycles")) {
      d_max_cycles = input_db->getInteger("max_cycles");
   } else if (input_db->keyExists("max_timesteps")) {
      d_max_cycles = input_db->getInteger("max_timesteps");
   } else {
      d_max_cycles = INT_MAX;
   }

   if (input_db->keyExists("max_delta_cycles")) {
      d_max_cycles = input_db->getInteger("max_delta_cycles");
   }

   d_start_time = input_db->getDoubleWithDefault("start_time", 0.0);

   if (input_db->keyExists("end_time")) {
      d_end_time = input_db->getDouble("end_time");
   } else {
      TBOX_ERROR("No \"end_time\" set in input file");
   }

   if (input_db->isDatabase("Amr")) {
      std::shared_ptr<tbox::Database> amr_db(input_db->getDatabase("Amr"));
      d_amr_enabled = amr_db->getBoolWithDefault("enabled", true);
   } else {
      d_amr_enabled = false;
   }
}

void PFModel::getFromRestart(std::shared_ptr<tbox::Database> input_db)
{
   std::shared_ptr<tbox::Database> root_db(
       tbox::RestartManager::getManager()->getRootDatabase());

   std::shared_ptr<tbox::Database> restart_db;

   if (root_db->isDatabase(d_object_name)) {
      restart_db = root_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
                 << d_object_name << " not found in the restart file.");
   }

   if (restart_db->keyExists("max_cycles")) {
      d_max_cycles = restart_db->getInteger("max_cycles");
   } else {
      d_max_cycles = restart_db->getInteger("max_timesteps");
   }

   d_start_time = restart_db->getDouble("start_time");
   d_end_time = restart_db->getDouble("end_time");

   if (restart_db->keyExists("time")) {
      d_time = restart_db->getDouble("time");
   } else {
      d_time = restart_db->getDouble("loop_time");
   }

   if (restart_db->keyExists("cycle")) {
      d_cycle = restart_db->getInteger("cycle");
   } else {
      d_cycle = restart_db->getInteger("iteration_number");
   }

   if (restart_db->keyExists("previous_timestep")) {
      d_previous_timestep = restart_db->getDouble("previous_timestep");
   } else {
      d_previous_timestep = 0.0;
   }

   if (restart_db->keyExists("amr_enabled")) {
      d_amr_enabled = restart_db->getBool("amr_enabled");
   } else {
      d_amr_enabled = true;
   }

   // Allow override of these from input deck

   if (input_db->keyExists("max_delta_cycles")) {
      int max_delta_cycles = input_db->getInteger("max_delta_cycles");
      d_max_cycles = d_cycle + max_delta_cycles;
   }

   if (input_db->keyExists("max_cycles")) {
      d_max_cycles = input_db->getInteger("max_cycles");
   } else if (input_db->keyExists("max_timesteps")) {
      d_max_cycles = input_db->getInteger("max_timesteps");
   }

   if (input_db->keyExists("end_time")) {
      d_end_time = input_db->getDouble("end_time");
   }
}

//=======================================================================
//
// Methods inherited from StandardTagAndInitStrategy
//

void PFModel::initializeLevelData(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int level_number, const double time, const bool can_be_refined,
    const bool initial_time, const std::shared_ptr<hier::PatchLevel>& old_level,
    const bool allocate_data)
{
   // Note that this method is pure virtual and MUST be implemented in derived
   // classes However, we might want some default behavior here, so it should be
   // called by the derived method
   (void)hierarchy;
   (void)level_number;
   (void)time;
   (void)can_be_refined;
   (void)initial_time;
   (void)old_level;
   (void)allocate_data;
}

void PFModel::resetHierarchyConfiguration(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int coarsest_level, const int finest_level)
{
   // Note that this method is pure virtual and MUST be implemented in derived
   // classes However, we might want some default behavior here, so it should be
   // called by the derived method
   (void)hierarchy;
   (void)coarsest_level;
   (void)finest_level;
}

//=======================================================================
