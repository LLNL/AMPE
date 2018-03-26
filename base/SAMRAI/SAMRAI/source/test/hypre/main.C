/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Main program for Hypre Poisson example
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include <string>
using namespace std;

#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/appu/VisItDataWriter.h"

#include "HyprePoisson.h"

#include <boost/shared_ptr.hpp>

using namespace SAMRAI;

/*
 ************************************************************************
 *
 * This is the driver program to demonstrate
 * how to use the Hypre Poisson solver.
 *
 * We set up the simple problem
 *          u + div(grad(u)) = sin(x)*sin(y)
 * in the domain [0:1]x[0:1], with u=0 on the
 * boundary.
 *
 * HyprePoisson is the primary object used to
 * set up and solve the system.  It maintains
 * the data for the computed solution u, the
 * exact solution, and the right hand side.
 *
 * The hierarchy created to solve this problem
 * has only one level.  (The Hypre Poisson solver
 * is a single-level solver.)
 *
 *************************************************************************
 */

int main(
   int argc,
   char* argv[])
{
   /*
    * Initialize MPI, SAMRAI, and enable logging.
    */

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {
      bool converged = true;

#if !defined(HAVE_HYPRE)
      tbox::pout << "This example requires the package HYPRE"
                 << "\nto work properly.  SAMRAI was not configured"
                 << "\nwith this package."
                 << endl;
#else

      /*
       * Process command line arguments.  For each run, the input
       * filename must be specified.  Usage is:
       *
       *    executable <input file name>
       *
       */
      string input_filename;

      if (argc != 2) {
         TBOX_ERROR("USAGE:  " << argv[0] << " <input file> \n"
                               << "  options:\n"
                               << "  none at this time" << endl);
      } else {
         input_filename = argv[1];
      }

      /*
       * Create input database and parse all data in input file.
       */

      boost::shared_ptr<tbox::InputDatabase> input_db(
         new tbox::InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Retrieve "Main" section from input database.
       * The main database is used only in main().
       * The base_name variable is a base name for
       * all name strings in this program.
       */

      boost::shared_ptr<tbox::Database> main_db(
         input_db->getDatabase("Main"));

      const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

      string base_name = "unnamed";
      base_name = main_db->getStringWithDefault("base_name", base_name);

      /*
       * Start logging.
       */
      const string log_file_name = base_name + ".log";
      bool log_all_nodes = false;
      log_all_nodes = main_db->getBoolWithDefault("log_all_nodes",
            log_all_nodes);
      if (log_all_nodes) {
         tbox::PIO::logAllNodes(log_file_name);
      } else {
         tbox::PIO::logOnlyNodeZero(log_file_name);
      }

      /*
       * Create major algorithm and data objects which comprise application.
       * Each object will be initialized either from input data or restart
       * files, or a combination of both.  Refer to each class constructor
       * for details.  For more information on the composition of objects
       * for this application, see comments at top of file.
       */

      boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
         new geom::CartesianGridGeometry(
            dim,
            base_name + "CartesianGeometry",
            input_db->getDatabase("CartesianGeometry")));
      tbox::plog << "Cartesian Geometry:" << endl;
      grid_geometry->printClassData(tbox::plog);

      boost::shared_ptr<hier::PatchHierarchy> patch_hierarchy(
         new hier::PatchHierarchy(
            base_name + "::PatchHierarchy",
            grid_geometry,
            input_db->getDatabase("PatchHierarchy")));

      /*
       * The HyprePoisson object is the main user object specific to the
       * problem being solved.  It provides the implementations for setting
       * up the grid and plotting data.  It also wraps up the solve
       * process that includes making the initial guess, specifying the
       * boundary conditions and call the solver.
       */
      HyprePoisson hypre_poisson(
         base_name + "::HyprePoisson",
         dim,
         input_db->isDatabase("HyprePoisson") ?
         input_db->getDatabase("HyprePoisson") :
         boost::shared_ptr<tbox::Database>());

      /*
       * Create the tag-and-initializer, box-generator and load-balancer
       * object references required by the gridding_algorithm object.
       */
      boost::shared_ptr<mesh::StandardTagAndInitialize> tag_and_initializer(
         new mesh::StandardTagAndInitialize(
            dim,
            "CellTaggingMethod",
            &hypre_poisson,
            input_db->getDatabase("StandardTagAndInitialize")));
      boost::shared_ptr<mesh::BergerRigoutsos> box_generator(
         new mesh::BergerRigoutsos(
            dim,
            input_db->getDatabase("BergerRigoutsos")));
      boost::shared_ptr<mesh::TreeLoadBalancer> load_balancer(
         new mesh::TreeLoadBalancer(
            dim,
            "load balancer",
            input_db->getDatabase("TreeLoadBalancer")));
      load_balancer->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());

      /*
       * Create the gridding algorithm used to generate the SAMR grid
       * and create the grid.
       */
      boost::shared_ptr<mesh::GriddingAlgorithm> gridding_algorithm(
         new mesh::GriddingAlgorithm(
            patch_hierarchy,
            "DistributedGridding Algorithm",
            input_db->getDatabase("GriddingAlgorithm"),
            tag_and_initializer,
            box_generator,
            load_balancer));
      tbox::plog << "Gridding algorithm:" << endl;
      gridding_algorithm->printClassData(tbox::plog);

      /*
       * Make the coarsest patch level where we will be solving.
       */
      gridding_algorithm->makeCoarsestLevel(0.0);

      /*
       * Set up the plotter for the hierarchy just created.
       * The FACPoisson object handles the data and has the
       * function setupExternalPlotter to register its data
       * with the plotter.
       */
      tbox::Array<string> vis_writer(1);
      vis_writer[0] = "Visit";
      if (main_db->keyExists("vis_writer")) {
         vis_writer = main_db->getStringArray("vis_writer");
      }
      bool use_visit = false;
      for (int i = 0; i < vis_writer.getSize(); i++) {
         if (vis_writer[i] == "VisIt") use_visit = true;
      }
      string vis_filename =
         main_db->getStringWithDefault("vis_filename", base_name);
#ifdef HAVE_HDF5
      boost::shared_ptr<appu::VisItDataWriter> visit_writer;
      if (use_visit) {
         visit_writer.reset(new appu::VisItDataWriter(dim,
               "Visit Writer",
               vis_filename + ".visit"));
         hypre_poisson.registerVariablesWithPlotter(*visit_writer);
      }
#endif

      /*
       * After creating all objects and initializing their state,
       * we print the input database and variable database contents
       * to the log file.
       */
      tbox::plog << "\nCheck input data and variables before simulation:"
                 << endl;
      tbox::plog << "Input database..." << endl;
      input_db->printClassData(tbox::plog);

      /*
       * Solve.
       */
      converged = hypre_poisson.solvePoisson();

      /*
       * Plot.
       */
#ifdef HAVE_HDF5
      if (use_visit) {
         visit_writer->writePlotData(patch_hierarchy, 0);
      }
#endif

      /*
       * Deallocate objects when done.
       */

      tbox::TimerManager::getManager()->print(tbox::plog);

#endif

#ifdef TESTING
      if (converged) {
         tbox::pout << "\nPASSED:  hypre" << endl;
      } else {
         TBOX_WARNING("Hypre test did not converge.");
      }
#endif
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return 0;
}
