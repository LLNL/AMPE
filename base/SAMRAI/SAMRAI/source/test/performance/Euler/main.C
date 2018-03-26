/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Main program for SAMRAI Euler gas dynamics sample application
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
using namespace std;

#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <sys/stat.h>

// Headers for basic SAMRAI objects

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/hier/BoxList.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/hier/VariableDatabase.h"

// Headers for major algorithm/data structure objects from SAMRAI

#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/algs/HyperbolicLevelIntegrator.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/algs/TimeRefinementIntegrator.h"
#include "SAMRAI/algs/TimeRefinementLevelStrategy.h"
#include "SAMRAI/appu/VisItDataWriter.h"

#define RECORD_STATS
//#undef RECORD_STATS
#ifdef RECORD_STATS
#include "SAMRAI/tbox/Statistic.h"
#include "SAMRAI/tbox/Statistician.h"
#endif

// Header for application-specific algorithm/data structure object

#include "Euler.h"

// Classes for autotesting.

#if (TESTING == 1)
#include "AutoTester.h"
#endif

#include <boost/shared_ptr.hpp>

void
outputStats(
   mesh::GriddingAlgorithm& gridding_algorithm,
   algs::HyperbolicLevelIntegrator& hyp_level_integrator);

using namespace SAMRAI;

using namespace tbox;

/*
 ************************************************************************
 *
 * This is the main program for an AMR Euler gas dynamics application
 * built using SAMRAI.   The application program is constructed by
 * composing a variety of algorithm objects found in SAMRAI plus some
 * others that are specific to this application.   The following brief
 * discussion summarizes these objects.
 *
 *    hier::PatchHierarchy - A container for the AMR patch hierarchy and
 *       the data on the grid.
 *
 *    geom::CartesianGridGeometry - Defines and maintains the Cartesian
 *       coordinate system on the grid.  The hier::PatchHierarchy
 *       maintains a reference to this object.
 *
 * A single overarching algorithm object drives the time integration
 * and adaptive gridding processes:
 *
 *    algs::TimeRefinementIntegrator - Coordinates time integration and
 *       adaptive gridding procedures for the various levels
 *       in the AMR patch hierarchy.  Local time refinement is
 *       employed during hierarchy integration; i.e., finer
 *       levels are advanced using smaller time increments than
 *       coarser level.  Thus, this object also invokes data
 *       synchronization procedures which couple the solution on
 *       different patch hierarchy levels.
 *
 * The time refinement integrator is not specific to the numerical
 * methods used and the problem being solved.   It maintains references
 * to two other finer grain algorithmic objects, more specific to
 * the problem at hand, with which it is configured when they are
 * passed into its constructor.   They are:
 *
 *    algs::HyperbolicLevelIntegrator - Defines data management procedures
 *       for level integration, data synchronization between levels,
 *       and tagging cells for refinement.  These operations are
 *       tailored to explicit time integration algorithms used for
 *       hyperbolic systems of conservation laws, such as the Euler
 *       equations.  This integrator manages data for numerical
 *       routines that treat individual patches in the AMR patch
 *       hierarchy.  In this particular application, it maintains a
 *       pointer to the Euler object that defines variables and
 *       provides numerical routines for the Euler model.
 *
 *       Euler - Defines variables and numerical routines for the
 *          discrete Euler equations on each patch in the AMR
 *          hierarchy.
 *
 *    mesh::GriddingAlgorithm - Drives the AMR patch hierarchy generation
 *       and regridding procedures.  This object maintains
 *       references to three other algorithmic objects with
 *       which it is configured when they are passed into its
 *       constructor.   They are:
 *
 *       mesh::BergerRigoutsos - Clusters cells tagged for refinement on a
 *          patch level into a collection of logically-rectangular
 *          box domains.
 *
 *       mesh::TreeLoadBalancer - Processes the boxes generated by the
 *          mesh::BergerRigoutsos algorithm into a configuration from
 *          which patches are contructed.  The algorithm we use in this
 *          class assumes a spatially-uniform workload distribution;
 *          thus, it attempts to produce a collection of boxes
 *          each of which contains the same number of cells.  The
 *          load balancer also assigns patches to processors.
 *
 *       mesh::StandardTagAndInitialize - Couples the gridding algorithm
 *          to the HyperbolicIntegrator. Selects cells for
 *          refinement based on either Gradient detection, Richardson
 *          extrapolation, or pre-defined Refine box region.  The
 *          object maintains a pointer to the algs::HyperbolicLevelIntegrator,
 *          which is passed into its constructor, for this purpose.
 *
 ************************************************************************
 */

/*
 *******************************************************************
 *
 * For each run, the input filename and restart information
 * (if needed) must be given on the command line.
 *
 *      For non-restarted case, command line is:
 *
 *          executable <input file name>
 *
 *      For restarted run, command line is:
 *
 *          executable <input file name> <restart directory> \
 *                     <restart number>
 *
 * Accessory routines used within the main program:
 *
 *   dumpVizData1dPencil - Writes 1d pencil of Euler solution data
 *      to plot files so that it may be viewed in MatLab.  This
 *      routine assumes a single patch level in 2d and 3d.  In
 *      other words, it only plots data on level zero.  It can
 *      handle AMR in 1d.
 *
 *******************************************************************
 */

int main(
   int argc,
   char* argv[])
{

   /*
    * Initialize MPI and SAMRAI, enable logging, and process command line.
    */

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   int num_failures = 0;

   {
      const tbox::Dimension dim(NDIM);

      string input_filename;
      string restart_read_dirname;
      int restore_num = 0;
      string case_name;

      bool is_from_restart = false;

      if ((argc != 2) && (argc != 3) && (argc != 4)) {
         pout << "USAGE:\n"
              << argv[0] << " <input filename> "
              << "or\n"
              << argv[0] << " <input filename> <case name>"
              << "or\n"
              << argv[0] << " <input filename> "
              << "  <restart dir> <restore number> [options]\n"
              << "  options:\n"
              << "  none at this time"
              << endl;
         tbox::SAMRAI_MPI::abort();
         return -1;
      } else {
         input_filename = argv[1];

         if (argc == 3) {
            case_name = argv[2];
         }
         if (argc == 4) {
            restart_read_dirname = argv[2];
            restore_num = atoi(argv[3]);

            is_from_restart = true;
         }
      }

      plog << "input_filename = " << input_filename << endl;
      plog << "restart_read_dirname = " << restart_read_dirname << endl;
      plog << "restore_num = " << restore_num << endl;

      /*
       * Create input database and parse all data in input file.
       */

      boost::shared_ptr<Database> input_db(new InputDatabase("input_db"));
      InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Retrieve "Main" section of the input database.  First, read dump
       * information, which is used for writing plot files.  Second, if
       * proper restart information was given on command line, and the
       * restart interval is non-zero, create a restart database.
       */

      boost::shared_ptr<Database> main_db(input_db->getDatabase("Main"));

      string base_name = main_db->getStringWithDefault("base_name", "unnamed");

      /*
       * Modify basename for this particular run.
       * Add the number of processes and the case name.
       */
      if (!case_name.empty()) {
         base_name = base_name + '-' + case_name;
      }
      base_name = base_name + '-' + tbox::Utilities::intToString(
            SAMRAI_MPI::getNodes(),
            5);
      pout << "Added case name (" << case_name << ") and nprocs ("
           << SAMRAI_MPI::getNodes() << ") to base name -> '"
           << base_name << "'\n";

      /*
       * Logging.
       */
      string log_filename = base_name + ".log";
      log_filename =
         main_db->getStringWithDefault("log_filename", base_name + ".log");

      bool log_all_nodes = false;
      log_all_nodes =
         main_db->getBoolWithDefault("log_all_nodes", log_all_nodes);
      if (log_all_nodes) {
         PIO::logAllNodes(log_filename);
      } else {
         PIO::logOnlyNodeZero(log_filename);
      }

      int viz_dump_interval = 0;
      if (main_db->keyExists("viz_dump_interval")) {
         viz_dump_interval = main_db->getInteger("viz_dump_interval");
      }

      Array<string> viz_writer(1);
      viz_writer[0] = "VisIt";
      string viz_dump_filename;
      string visit_dump_dirname = base_name + ".visit";
      bool uses_visit = false;
      int visit_number_procs_per_file = 1;
      if (viz_dump_interval > 0) {
         if (main_db->keyExists("viz_writer")) {
            viz_writer = main_db->getStringArray("viz_writer");
         }
         if (main_db->keyExists("viz_dump_filename")) {
            viz_dump_filename = main_db->getString("viz_dump_filename");
         }
         string viz_dump_dirname = base_name;
         if (main_db->keyExists("viz_dump_dirname")) {
            viz_dump_dirname = main_db->getString("viz_dump_dirname");
         }
         for (int i = 0; i < viz_writer.getSize(); i++) {
            if (viz_writer[i] == "VisIt") uses_visit = true;
         }
         visit_dump_dirname = viz_dump_dirname + ".visit";
         if (uses_visit) {
            if (viz_dump_dirname.empty()) {
               TBOX_ERROR("main(): "
                  << "\nviz_dump_dirname is null ... "
                  << "\nThis must be specified for use with VisIt"
                  << endl);
            }
            if (main_db->keyExists("visit_number_procs_per_file")) {
               visit_number_procs_per_file =
                  main_db->getInteger("visit_number_procs_per_file");
            }
         }
      }

      int restart_interval = 0;
      if (main_db->keyExists("restart_interval")) {
         restart_interval = main_db->getInteger("restart_interval");
      }

      string restart_write_dirname;
      if (restart_interval > 0) {
         if (main_db->keyExists("restart_write_dirname")) {
            restart_write_dirname = main_db->getString("restart_write_dirname");
         } else {
            TBOX_ERROR(
               "`restart_interval' > 0, but key `restart_write_dirname'"
               << "not found in input file.");
         }
      }

      bool use_refined_timestepping = true;
      if (main_db->keyExists("timestepping")) {
         string timestepping_method = main_db->getString("timestepping");
         if (timestepping_method == "SYNCHRONIZED") {
            use_refined_timestepping = false;
         }
      }

#if (TESTING == 1) && !(HAVE_HDF5)
      /*
       * If we are autotesting on a system w/o HDF5, the read from
       * restart will result in an error.  We want this to happen
       * for users, so they know there is a problem with the restart,
       * but we don't want it to happen when autotesting.
       */
      is_from_restart = false;
      restart_interval = 0;
#endif

      const bool write_restart = (restart_interval > 0)
         && !(restart_write_dirname.empty());

      /*
       * Get restart manager and root restart database.  If run is from
       * restart, open the restart file.
       */

      RestartManager* restart_manager = RestartManager::getManager();

      if (is_from_restart) {
         restart_manager->
         openRestartFile(restart_read_dirname, restore_num,
            tbox::SAMRAI_MPI::getNodes());
      }

      /*
       * Setup the timer manager to trace timing statistics during execution
       * of the code.  The list of timers is given in the TimerManager
       * section of the input file.  Timing information is stored in the
       * restart file.  Timers will automatically be initialized to their
       * previous state if the run is restarted, unless they are explicitly
       * reset using the TimerManager::resetAllTimers() routine.
       */

      TimerManager::createManager(input_db->getDatabase("TimerManager"));

      /*
       * Create major algorithm and data objects which comprise application.
       * Each object is initialized either from input data or restart
       * files, or a combination of both.  Refer to each class constructor
       * for details.  For more information on the composition of objects
       * and the roles they play in this application, see comments at top of file.
       */

      boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
         new geom::CartesianGridGeometry(
            "CartesianGeometry",
            input_db->getDatabase("CartesianGeometry")));

      boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy(
         new hier::PatchHierarchy(
            "PatchHierarchy",
            grid_geometry,
            input_db->getDatabase("PatchHierarchy")));

      Euler* euler_model = new Euler("Euler",
            input_db->getDatabase("Euler"),
            grid_geometry);

      boost::shared_ptr<algs::HyperbolicLevelIntegrator> hyp_level_integrator(
         new algs::HyperbolicLevelIntegrator(
            "HyperbolicLevelIntegrator",
            input_db->getDatabase("HyperbolicLevelIntegrator"),
            euler_model,
            true,
            use_refined_timestepping));

      boost::shared_ptr<mesh::StandardTagAndInitialize> error_detector(
         new mesh::StandardTagAndInitialize(
            "StandardTagAndInitialize",
            hyp_level_integrator,
            input_db->getDatabase("StandardTagAndInitialize")));

      boost::shared_ptr<Database> abr_db(
         input_db->getDatabase("BergerRigoutsos"));
      boost::shared_ptr<mesh::BergerRigoutsos> new_box_generator(
         new mesh::BergerRigoutsos(abr_db));
#if 0
      boost::shared_ptr<mesh::BergerRigoutsos<NDIM> > old_box_generator;
      const char which_br = main_db->getCharWithDefault("which_br", 'o');
      boost::shared_ptr<mesh::BoxGeneratorStrategy> box_generator(
         which_br == 'o'
         ? boost::shared_ptr<mesh::BoxGeneratorStrategy>(old_box_generator)
         : boost::shared_ptr<mesh::BoxGeneratorStrategy>(new_box_generator));
#endif

      boost::shared_ptr<mesh::TreeLoadBalancer> load_balancer(
         new mesh::TreeLoadBalancer(
            dim,
            "TreeLoadBalancer",
            input_db->getDatabase("TreeLoadBalancer")));
      load_balancer->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

      boost::shared_ptr<mesh::GriddingAlgorithm> gridding_algorithm(
         new mesh::GriddingAlgorithm(
            dim,
            "GriddingAlgorithm",
            input_db->getDatabase("GriddingAlgorithm"),
            error_detector,
            new_box_generator,
            load_balancer));

      boost::shared_ptr<algs::TimeRefinementIntegrator> time_integrator(
         new algs::TimeRefinementIntegrator(
            "TimeRefinementIntegrator",
            input_db->getDatabase("TimeRefinementIntegrator"),
            patch_hierarchy,
            hyp_level_integrator,
            gridding_algorithm));

      /*
       * Set up Visualization writer.  Note that the Euler application
       * creates some derived data quantities so we register the Euler model
       * as a derived data writer.  If no derived data is written, this step
       * is not necessary.
       */
#ifdef HAVE_HDF5
      boost::shared_ptr<appu::VisItDataWriter> visit_data_writer(
         new appu::VisItDataWriter(
            "Euler VisIt Writer",
            visit_dump_dirname,
            visit_number_procs_per_file));
      if (uses_visit) {
         euler_model->registerVisItDataWriter(visit_data_writer);
      }
#endif

      /*
       * Initialize hierarchy configuration and data on all patches.
       * Then, close restart file and write initial state for visualization.
       */

      hier::Connector dummy; // Cause communicator set-up before performance timings.
      tbox::SAMRAI_MPI::getSAMRAIWorld().Barrier(); // Synchronize for the sake of accurate timings.
      double dt_now = time_integrator->initializeHierarchy();

      RestartManager::getManager()->closeRestartFile();

#if (TESTING == 1)
      /*
       * Create the autotesting component which will verify correctness
       * of the problem. If no automated testing is done, the object does
       * not get used.
       */
      AutoTester autotester("AutoTester", input_db);
#endif

      /*
       * After creating all objects and initializing their state, we
       * print the input database and variable database contents
       * to the log file.
       */

#if 1
      plog << "\nCheck input data and variables before simulation:" << endl;
      plog << "Input database..." << endl;
      input_db->printClassData(plog);
      plog << "\nVariable database..." << endl;
      hier::VariableDatabase::getDatabase()->printClassData(plog);

#endif
      plog << "\nCheck Euler data... " << endl;
      euler_model->printClassData(plog);

      /*
       * Create timers for measuring I/O.
       */
      boost::shared_ptr<Timer> t_write_viz(
         TimerManager::getManager()->getTimer("apps::main::write_viz"));
      boost::shared_ptr<Timer> t_write_restart(
         TimerManager::getManager()->getTimer("apps::main::write_restart");

      t_write_viz->start();
      if (viz_dump_interval > 0) {
#ifdef HAVE_HDF5
         if (uses_visit) {
            visit_data_writer->writePlotData(
               patch_hierarchy,
               time_integrator->getIntegratorStep(),
               time_integrator->getIntegratorTime());
         }
#endif
      }
      t_write_viz->stop();

      /*
       * Time step loop.  Note that the step count and integration
       * time are maintained by algs::TimeRefinementIntegrator.
       */

      double loop_time = time_integrator->getIntegratorTime();
      double loop_time_end = time_integrator->getEndTime();

#if (TESTING == 1)
      /*
       * If we are doing autotests, check result...
       */
      num_failures += autotester.evalTestData(
            time_integrator->getIntegratorStep(),
            patch_hierarchy,
            time_integrator,
            hyp_level_integrator,
            gridding_algorithm);
#endif

      while ((loop_time < loop_time_end) &&
             time_integrator->stepsRemaining()) {

         int iteration_num = time_integrator->getIntegratorStep() + 1;

         plog << endl << endl;
         pout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;
         pout << "At begining of timestep # " << iteration_num - 1 << endl;
         pout << "Simulation time is " << loop_time << endl;
         pout << "Current dt is " << dt_now << endl;
         pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
         plog << endl << endl;

         double dt_new = time_integrator->advanceHierarchy(dt_now);

         loop_time += dt_now;
         dt_now = dt_new;

#if 0
         int hierarchy_cell_count = 0;
         for (int ln = 0; ln < patch_hierarchy->getNumberLevels(); ++ln) {
            int level_cell_count = 0;
            boost::shared_ptr<hier::PatchLevel> patch_level(
               patch_hierarchy->getPatchLevel(ln));
            for (hier::PatchLevel::iterator pi(patch_level->begin());
                 pi != patch_level->end(); ++pi) {
               boost::shared_ptr<hier::Patch> patch(
                  patch_level->getPatch(*pi));
               level_cell_count += patch->getBox().size();
            }
            cell_count_stat[ln]->recordProcStat(level_cell_count);
            hierarchy_cell_count += level_cell_count;
         }
         for (int ln = patch_hierarchy->getNumberLevels();
              ln < patch_hierarchy->getMaxNumberOfLevels(); ++ln) {
            cell_count_stat[ln]->recordProcStat(0);
         }
         cell_count_stat[patch_hierarchy->getMaxNumberOfLevels()]->
         recordProcStat(hierarchy_cell_count);
         patch_hierarchy->recursivePrint(plog, "", 2);
         sim_time_stat->recordProcStat(dt_now);
#endif

         plog << "Hierarchy summary:\n";
         patch_hierarchy->recursivePrint(plog, "H-> ", 1);
         plog << "PatchHierarchy summary:\n";
         patch_hierarchy->recursivePrint(plog,
            "H-> ",
            1);

         plog << endl << endl;
         pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
         pout << "At end of timestep # " << iteration_num - 1 << endl;
         pout << "Simulation time is " << loop_time << endl;
         pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
         plog << endl << endl;

         /*
          * At specified intervals, write restart files.
          */
         if (write_restart) {

            if ((iteration_num % restart_interval) == 0) {
               t_write_restart->start();
               RestartManager::getManager()->
               writeRestartFile(restart_write_dirname,
                  iteration_num);
               t_write_restart->stop();
            }
         }

         /*
          * At specified intervals, write out data files for plotting.
          */
         t_write_viz->start();
         if ((viz_dump_interval > 0)
             && (iteration_num % viz_dump_interval) == 0) {
#ifdef HAVE_HDF5
            if (uses_visit) {
               visit_data_writer->writePlotData(patch_hierarchy,
                  iteration_num,
                  loop_time);
            }
#endif

         }
         t_write_viz->stop();

#if (TESTING == 1)
         /*
          * If we are doing autotests, check result...
          */
         num_failures += autotester.evalTestData(iteration_num,
               patch_hierarchy,
               time_integrator,
               hyp_level_integrator,
               gridding_algorithm);
#endif

         /*
          * Write byte transfer information to log file.
          */
#if 0
         char num_buf[8];
         sprintf(num_buf, "%02d", iteration_num);
         tbox::plog << "Step " << num_buf
                    << " P" << tbox::SAMRAI_MPI::getRank()
                    << ": " << tbox::SAMRAI_MPI::getIncomingBytes()
                    << " bytes in" << endl;
#endif

      }

      /*
       * Output timer results.
       */
      TimerManager::getManager()->print(plog);

#ifdef RECORD_STATS
      outputStats(*gridding_algorithm, *hyp_level_integrator);
#endif

      /*
       * At conclusion of simulation, deallocate objects.
       */
      patch_hierarchy.reset();
      grid_geometry.reset();

      new_box_generator.reset();
      load_balancer.reset();
      hyp_level_integrator.reset();
      error_detector.reset();
      gridding_algorithm.reset();
      time_integrator.reset();

#ifdef HAVE_HDF5
      visit_data_writer.reset();
#endif

      if (euler_model) delete euler_model;

   }

   if (num_failures == 0) {
      tbox::pout << "\nPASSED:  Euler" << endl;
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return num_failures;
}

#ifdef RECORD_STATS
void outputStats(
   mesh::GriddingAlgorithm& gridding_algorithm,
   algs::HyperbolicLevelIntegrator& hyp_level_integrator)
{
   /*
    * Output statistics.
    */
   tbox::plog << "HyperbolicLevelIntegrator statistics:" << endl;
   hyp_level_integrator.printStatistics(tbox::plog);
   tbox::plog << "\nGriddingAlgorithm statistics:" << endl;
   gridding_algorithm.printStatistics(tbox::plog);
}
#endif
