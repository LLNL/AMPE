/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Test program for performance of tree search algorithm.
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/geom/GridGeometry.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/HDFDatabase.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <algorithm>
#include <vector>
#include <iomanip>

using namespace SAMRAI;
using namespace tbox;

/*
 * Break up boxes in the given BoxLevel.  This method is meant
 * to create a bunch of small boxes from the user-input boxes in order
 * to set up a non-trivial mesh configuration.
 */
void
breakUpBoxes(
   hier::BoxLevel& mapped_box_level,
   const hier::IntVector& max_box_size);

/*
 * Find overlapping Boxes using an exhaustive search.
 *
 * Intended as an alternate way of finding overlaps, for verifying
 * results from tree search.
 */
void
exhaustiveFindOverlapBoxes(
   hier::Connector& overlap_connector,
   const hier::Box& mapped_box,
   const hier::IntVector& refinement_ratio,
   const boost::shared_ptr<const hier::BaseGridGeometry>& grid_geometry,
   const hier::BoxContainer& search_mapped_boxes);

/*
 ************************************************************************
 *
 * This is accuracy test for the multiblock tree search algorithm
 * in BoxContainer:
 *
 * 1. Generate a set of Boxes.
 *
 * 2. Sort the Boxes into trees using makeTree.
 *
 * 3. Search for overlaps.
 *
 *************************************************************************
 */

int main(
   int argc,
   char* argv[])
{
   /*
    * Initialize MPI, SAMRAI.
    */

   SAMRAI_MPI::init(&argc, &argv);
   SAMRAIManager::initialize();
   SAMRAIManager::startup();
   tbox::SAMRAI_MPI mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   if (mpi.getSize() != 1) {
      TBOX_ERROR("mbtree test is intended to run on just one processor.");
   }

   const int rank = mpi.getRank();
   int fail_count = 0;

   {

      /*
       * Process command line arguments.  For each run, the input
       * filename must be specified.  Usage is:
       *
       * executable <input file name>
       */
      std::string input_filename;

      if (argc != 2) {
         TBOX_ERROR("USAGE:  " << argv[0] << " <input file> \n"
                               << "  options:\n"
                               << "  none at this time" << std::endl);
      } else {
         input_filename = argv[1];
      }

      /*
       * Create input database and parse all data in input file.
       */

      boost::shared_ptr<InputDatabase> input_db(new InputDatabase("input_db"));
      tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

      /*
       * Set up the timer manager.
       */
      if (input_db->isDatabase("TimerManager")) {
         TimerManager::createManager(input_db->getDatabase("TimerManager"));
      }

      /*
       * Retrieve "Main" section from input database.
       * The main database is used only in main().
       * The base_name variable is a base name for
       * all name strings in this program.
       */

      boost::shared_ptr<Database> main_db(input_db->getDatabase("Main"));

      const tbox::Dimension dim(static_cast<unsigned short>(main_db->getInteger("dim")));

      std::string base_name = "unnamed";
      base_name = main_db->getStringWithDefault("base_name", base_name);

      /*
       * Start logging.
       */
      const std::string log_filename = base_name + ".log";
      bool log_all_nodes = false;
      log_all_nodes = main_db->getBoolWithDefault("log_all_nodes",
            log_all_nodes);
      if (log_all_nodes) {
         PIO::logAllNodes(log_filename);
      } else {
         PIO::logOnlyNodeZero(log_filename);
      }

      plog << "Input database after initialization..." << std::endl;
      input_db->printClassData(plog);

      /*
       * Generate the grid geometry.
       */
      if (!main_db->keyExists("GridGeometry")) {
         TBOX_ERROR("Multiblock tree search test: could not find entry GridGeometry"
            << "\nin input.");
      }
      boost::shared_ptr<const hier::BaseGridGeometry> grid_geometry(
         new geom::GridGeometry(
            dim,
            "GridGeometry",
            main_db->getDatabase("GridGeometry")));

      /*
       * Baseline stuff:
       *
       * Whether to generate a baseline or compare against it.
       * If generating a baseline, the tests are NOT checked!
       */

      const std::string baseline_dirname = main_db->getString("baseline_dirname");
      const std::string baseline_filename = baseline_dirname + "/" + base_name + ".baselinedb";
      tbox::HDFDatabase basline_db(baseline_filename);

      const bool generate_baseline =
         main_db->getBoolWithDefault("generate_baseline", false);

      boost::shared_ptr<tbox::HDFDatabase> baseline_db(
         new tbox::HDFDatabase("mbtree baseline"));
      boost::shared_ptr<tbox::Database> mapped_box_level_db;
      boost::shared_ptr<tbox::Database> connector_db;
      if (generate_baseline) {
         baseline_db->create(baseline_filename);
         mapped_box_level_db = baseline_db->putDatabase("MappedBoxLevel");
         connector_db = baseline_db->putDatabase("Connector");
      } else {
         baseline_db->open(baseline_filename);
         mapped_box_level_db = baseline_db->getDatabase("MappedBoxLevel");
         connector_db = baseline_db->getDatabase("Connector");
      }

      /*
       * Print input database again to fully show usage.
       */
      plog << "Input database after running..." << std::endl;
      input_db->printClassData(plog);

      const hier::IntVector& one_vector(hier::IntVector::getOne(dim));

      hier::BoxLevel mapped_box_level(
         one_vector,
         grid_geometry,
         tbox::SAMRAI_MPI::getSAMRAIWorld());
      grid_geometry->computePhysicalDomain(
         mapped_box_level,
         hier::IntVector::getOne(dim));
      mapped_box_level.finalize();

      /*
       * Generate boxes from the multiblock domain description.
       */
      hier::IntVector max_box_size(dim, tbox::MathUtilities<int>::getMax());
      if (main_db->isInteger("max_box_size")) {
         main_db->getIntegerArray("max_box_size", &max_box_size[0], dim.getValue());
      }
      breakUpBoxes(mapped_box_level, max_box_size);

      /*
       * Write the baseline BoxLevel or check to ensure it is
       * the same as the one in the baseline database.  The regression
       * test is invalid of we don't have the same BoxLevel.
       */
      if (generate_baseline) {
         tbox::pout << "\nBoxLevel for review:\n"
                    << mapped_box_level.format("REVIEW: ", 2)
                    << std::endl;
         mapped_box_level.putUnregisteredToDatabase(mapped_box_level_db);
      } else {
         /*
          * Get the baselined BoxLevel and compare.
          */
         hier::BoxLevel baseline_mapped_box_level(dim);
         baseline_mapped_box_level.getFromDatabase(
            *mapped_box_level_db,
            grid_geometry);
         if (mapped_box_level != baseline_mapped_box_level) {
            tbox::perr << "Multiblock Tree test problem:\n"
                       << "the BoxLevel generated is different\n"
                       << "from the one in the database.  Thus the check\n"
                       << "cannot be done.\n";
            ++fail_count;
            tbox::pout << mapped_box_level.format("M: ", 2)
                       << std::endl
                       << baseline_mapped_box_level.format("B: ", 2);
         }
      }

      /*
       * Generate boxes from the multiblock domain description.
       */

      const hier::BoxContainer& level_boxes = mapped_box_level.getBoxes();
      level_boxes.makeTree(grid_geometry.get());

      /*
       * Find overlaps.
       */
      hier::IntVector connector_width(dim, 1);
      if (main_db->isInteger("connector_width")) {
         main_db->getIntegerArray("connector_width", &connector_width[0], dim.getValue());
      }

      hier::Connector connector(
         mapped_box_level,
         mapped_box_level,
         connector_width);

      const hier::IntVector& refinement_ratio(one_vector);

      for (hier::BoxContainer::const_iterator bi = mapped_box_level.getBoxes().begin();
           bi != mapped_box_level.getBoxes().end(); ++bi) {

         const hier::Box& mapped_box(*bi);

         hier::Box grown_box(mapped_box);
         grown_box.grow(connector_width);

         hier::BoxContainer overlap_boxes;
         level_boxes.findOverlapBoxes(
            overlap_boxes,
            grown_box,
            refinement_ratio,
            true);

         connector.insertNeighbors(overlap_boxes, mapped_box.getId());
      }

      /*
       * Write the baseline NeighborhoodSet or check against it.
       */
      if (generate_baseline) {

         /*
          * If writing baseline, verify the results against the
          * exhaustive search method first.
          */
         hier::Connector connector_from_exhaustive_search(
            mapped_box_level,
            mapped_box_level,
            connector_width);
         for (hier::BoxContainer::const_iterator bi =
              mapped_box_level.getBoxes().begin();
              bi != mapped_box_level.getBoxes().end(); ++bi) {

            const hier::Box& mapped_box(*bi);

            hier::Box grown_mapped_box(mapped_box);
            grown_mapped_box.grow(connector_width);

            exhaustiveFindOverlapBoxes(
               connector_from_exhaustive_search,
               grown_mapped_box,
               refinement_ratio,
               grid_geometry,
               mapped_box_level.getBoxes());

         }

         if (!connector.localNeighborhoodsEqual(
                 connector_from_exhaustive_search)) {

            tbox::perr << "Failed verification in baseline generation:\n"
                       << "Neighborhoods from the tree search do not match\n"
                       << "neighborhoods from exhaustive search.\n"
                       << "You should determine the cause and fix the mismatch\n"
                       << "before generating the baseline.\n"
                       << "Connector from tree search:\n"
                       << connector.format("TREE: ", 2)
                       << "Connector from exhautive search:\n"
                       << connector_from_exhaustive_search.format("EXHAUSTIVE: ", 2)
                       << std::endl;

            hier::Connector exhaustive_minus_tree, tree_minus_exhaustive;
            hier::Connector::computeNeighborhoodDifferences(
               exhaustive_minus_tree,
               connector_from_exhaustive_search,
               connector);
            hier::Connector::computeNeighborhoodDifferences(
               tree_minus_exhaustive,
               connector,
               connector_from_exhaustive_search);
            tbox::perr << "What's found by exhaustive search but not by tree search:\n"
                       << exhaustive_minus_tree.format("", 2)
                       << "\nWhat's found by tree search but not by exhaustive search:\n"
                       << tree_minus_exhaustive.format("", 2)
                       << std::endl;

            tbox::perr << "Baseline was NOT generated due to the above problem!"
                       << std::endl;
            ++fail_count;
         } else {

            connector.putNeighborhoodsToDatabase(connector_db);
            tbox::pout << "Connector for review:\n"
                       << connector.format("REVIEW: ", 2)
                       << "This data has been verified by comparing against the results\n"
                       << "of an exhaustive search and a new baseline has been written.\n"
                       << "However, you may wan to manually verify\n"
                       << "it before using the new baseline.\n"
                       << std::endl;
         }

      } else {
         /*
          * Get the baseline Connector NeighborhoodSet and compare.
          */
         hier::Connector baseline_connector;
         baseline_connector.getNeighborhoodsFromDatabase(*connector_db);
         if (!baseline_connector.localNeighborhoodsEqual(connector)) {
            tbox::perr << "Multiblock Tree test problem:\n"
                       << "the NeighborhoodSets generated is different\n"
                       << "from the one in the database.\n"
                       << "computed neighborhood set:\n";
            connector.writeNeighborhoodsToErrorStream("COMPUTED:  ");
            tbox::perr << "baseline neighborhood set\n";
            baseline_connector.writeNeighborhoodsToErrorStream("BASELINE:  ");
            ++fail_count;
         }
      }

      /*
       * Remind user to manually check the new baseline.
       */
      if (generate_baseline) {
         tbox::pout << "NEW BASELINE GENERATED!\n"
                    <<
         "Please manually review the output for accuracy before using the baseline.\n"
                    << "If it is correct, accept it by setting generate_baseline = FALSE\n"
                    << "in input file."
                    << std::endl;
      }

      if (fail_count == 0) {
         tbox::pout << "\nPASSED:  Multiblock tree search" << std::endl;
      }

      input_db.reset();
      main_db.reset();

      /*
       * Exit properly by shutting down services in correct order.
       */
      tbox::plog << "\nShutting down..." << std::endl;

   }

   /*
    * Shut down.
    */
   SAMRAIManager::shutdown();
   SAMRAIManager::finalize();

   if (fail_count == 0) {
      SAMRAI_MPI::finalize();
   } else {
      tbox::pout << "Process " << std::setw(5) << rank << " aborting."
                 << std::endl;
      tbox::Utilities::abort("Aborting due to nonzero fail count",
         __FILE__, __LINE__);
   }

   tbox::plog << "Process " << std::setw(5) << rank << " exiting." << std::endl;
   return fail_count;
}

/*
 * Break up boxes in the given BoxLevel.  This method is meant
 * to create a bunch of small boxes from the user-input boxes in order
 * to set up a non-trivial mesh configuration.
 */
void breakUpBoxes(
   hier::BoxLevel& mapped_box_level,
   const hier::IntVector& max_box_size) {

   const tbox::Dimension& dim(mapped_box_level.getDim());

   hier::BoxLevel domain_mapped_box_level(mapped_box_level);
   domain_mapped_box_level.setParallelState(hier::BoxLevel::GLOBALIZED);

   mesh::TreeLoadBalancer load_balancer(mapped_box_level.getDim());

   hier::Connector dummy_connector;

   const hier::IntVector min_size(dim, 2);
   const hier::IntVector bad_interval(dim, 1);
   const hier::IntVector cut_factor(dim, 1);

   load_balancer.loadBalanceBoxLevel(
      mapped_box_level,
      dummy_connector,
      dummy_connector,
      boost::shared_ptr<hier::PatchHierarchy>(),
      0,
      dummy_connector,
      dummy_connector,
      min_size,
      max_box_size,
      domain_mapped_box_level,
      bad_interval,
      cut_factor);
}

/*
 * Find overlapping Boxes using an exhaustive search.
 *
 * Intended as an alternate way of finding overlaps, for verifying
 * results from tree search.
 */
void exhaustiveFindOverlapBoxes(
   hier::Connector& overlap_connector,
   const hier::Box& mapped_box,
   const hier::IntVector& refinement_ratio,
   const boost::shared_ptr<const hier::BaseGridGeometry>& grid_geometry,
   const hier::BoxContainer& search_mapped_boxes)
{
   const hier::BoxId& box_id = mapped_box.getId();
   hier::Box transformed_box(mapped_box);
   hier::BlockId transformed_block_id(mapped_box.getBlockId());
   hier::Connector::NeighborhoodIterator base_box_itr =
      overlap_connector.findLocal(box_id);
   bool has_base_box = base_box_itr != overlap_connector.end();

   for (hier::BoxContainer::const_iterator bi = search_mapped_boxes.begin();
        bi != search_mapped_boxes.end(); ++bi) {

      const hier::Box& search_mapped_box(*bi);

      /*
       * Get transformed_box in coordinate of search_mapped_box, if
       * the two blocks are neighbors.  If the blocks are not
       * neighbors, there can't be any overlap between mapped_box and
       * search_mapped_box.
       */
      if (transformed_block_id != search_mapped_box.getBlockId()) {
         transformed_box = mapped_box;
         bool transformed = grid_geometry->transformBox(transformed_box,
               refinement_ratio,
               search_mapped_box.getBlockId(),
               mapped_box.getBlockId());
         transformed_block_id = transformed ?
            search_mapped_box.getBlockId() :
            mapped_box.getBlockId();
      }

      if (transformed_block_id == search_mapped_box.getBlockId()) {
         if (transformed_box.intersects(search_mapped_box)) {
            if (!has_base_box) {
               base_box_itr = overlap_connector.makeEmptyLocalNeighborhood(
                  box_id);
               has_base_box = true;
            }
            overlap_connector.insertLocalNeighbor(
               search_mapped_box,
               base_box_itr);
         }
      }
   }
}
