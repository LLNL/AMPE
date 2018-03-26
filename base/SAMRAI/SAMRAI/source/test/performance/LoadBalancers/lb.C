/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Test program for performance and quality of TreeLoadBalancer.
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include <iomanip>

#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/BoxLevelStatistics.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/ConnectorStatistics.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/mesh/BalanceUtilities.h"
#include "SAMRAI/mesh/TreeLoadBalancer.h"
#include "SAMRAI/mesh/TreeLoadBalancerOld.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/appu/VisItDataWriter.h"

#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include <vector>

#include "DerivedVisOwnerData.h"
#include "SinusoidalFrontTagger.h"

using namespace SAMRAI;
using namespace tbox;

/*
 ************************************************************************
 *
 *
 *************************************************************************
 */

void
generatePrebalance(
   hier::BoxLevel &Lb,
   hier::BoxLevel &La,
   hier::Connector &La_to_Lb,
   hier::Connector &Lb_to_La,
   const std::string &Lb_box_gen_method,
   tbox::Database &main_db,
   const boost::shared_ptr<hier::PatchHierarchy> &hierarchy,
   int coarser_ln,
   const hier::IntVector &min_size,
   const hier::IntVector &required_connector_width );

void
generatePrebalanceByUserBoxes(
   hier::BoxLevel& L1,
   hier::Connector& L0_to_L1,
   hier::Connector& balance_to_L0,
   const hier::BoxLevel& L0,
   boost::shared_ptr<tbox::Database> database,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const hier::IntVector& min_size,
   const hier::IntVector& connector_width );

void
generatePrebalanceByUserShells(
   hier::BoxLevel& L1,
   hier::Connector& L0_to_L1,
   hier::Connector& L1_to_L0,
   const hier::BoxLevel& L0,
   const boost::shared_ptr<tbox::Database>& database,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const hier::IntVector& min_size,
   const hier::IntVector& connector_width );

void
generatePrebalanceBySinusoidalFront(
   hier::BoxLevel& L1,
   hier::Connector& L0_to_L1,
   hier::Connector& L1_to_L0,
   const hier::BoxLevel& L0,
   const boost::shared_ptr<tbox::Database>& database,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int coarser_ln,
   const hier::IntVector& min_size,
   const hier::IntVector& connector_width );

void
generatePrebalanceByShrinkingLevel(
   hier::BoxLevel& L2,
   hier::Connector& L1_to_L2,
   hier::Connector& L2_to_L1,
   const hier::BoxLevel& L1,
   const boost::shared_ptr<tbox::Database>& database,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int coarser_ln,
   const hier::IntVector& min_size,
   const hier::IntVector& connector_width );

void
sortNodes(
   hier::BoxLevel& new_box_level,
   hier::Connector& tag_to_new,
   hier::Connector& new_to_tag,
   bool sort_by_corners,
   bool sequentialize_global_indices);

void
refineHead(
   hier::BoxLevel& head,
   hier::Connector& ref_to_head,
   hier::Connector& head_to_ref,
   const hier::IntVector &refinement_ratio );

void outputMetadataL0(
   const hier::BoxLevel &L0,
   const hier::Connector &L0_to_L0,
   int level_output_depth = 1,
   int connector_output_depth = 0 );

void outputMetadataBefore(
   const hier::Connector &La_to_Lb,
   const hier::Connector &Lb_to_La,
   const std::string &La_name,
   const std::string &Lb_name,
   int level_output_depth = 1,
   int connector_output_depth = 0 );

void outputMetadataAfter(
   const hier::Connector &La_to_Lb,
   const hier::Connector &Lb_to_La,
   const hier::Connector &Lb_to_Lb,
   const std::string &La_name,
   const std::string &Lb_name,
   int level_output_depth = 1,
   int connector_output_depth = 0 );

boost::shared_ptr<mesh::LoadBalanceStrategy>
createLoadBalancer(
   boost::shared_ptr<tbox::Database> &input_db,
   const std::string &lb_type,
   int ln,
   const tbox::Dimension &dim );


/*
********************************************************************************
*
* Performance testing for load balancers.
*
* 1. Build "level 0" from the domain description (input parameter
* "domain_boxes").  L0 is for doing the test, not for checking load
* balancer performance.
*
* 2. Build "level 1" and write out performance data for balancing it.
* The prebalance boxes for L1 are user-specified (input parameter
* "L1_box_gen_method").  This configuration tries to mimick real problems
* where the tags occupy a small portion of the tag level, leading to a
* limited number of owners for prebalance boxes.
*
* 3. Build "level 2" and write out performance data for balancing it.
* The prebalance boxes for L2 are generated by clustering tags on L1.
* All L1 cells are tagged except for a small margin by the L1 boundary
* (input parameter "tag_margin".  This configuration tries to mimick
* real problems where the tags occupy a large portion of the tag
* level, leading to a greater number of owners for prebalance boxes.
*
********************************************************************************
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

   const int rank = mpi.getRank();
   int fail_count = 0;

   /*
    * Process command line arguments.  For each run, the input
    * filename must be specified.  Usage is:
    *
    * executable <input file name>
    */
   std::string input_filename;

   if (argc < 2) {
      TBOX_ERROR("USAGE:  " << argv[0] << " <input file> [case name]\n"
                            << "  options:\n"
                            << "  none at this time" << std::endl);
   } else {
      input_filename = argv[1];
   }

   std::string case_name;
   if (argc > 2) {
      case_name = argv[2];
   }

   {
      /*
       * Scope to force destruction of objects that would otherwise
       * leave allocated memory reported by the memory test.
       */

      /*
       * Create input database and parse all data in input file.
       */

      boost::shared_ptr<InputDatabase> input_db(new InputDatabase("input_db"));
      boost::shared_ptr<Database> base_db(input_db);
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

      const tbox::Dimension
         dim(static_cast<unsigned short>(main_db->getInteger("dim")));

      const hier::IntVector &zero_vec = hier::IntVector::getZero(dim);

      std::string base_name = "unnamed";
      base_name = main_db->getStringWithDefault("base_name", base_name);

      /*
       * Modify basename for this particular run.
       * Add the number of processes and the case name.
       */
      if (!case_name.empty()) {
         base_name = base_name + '-' + case_name;
      }
      base_name = base_name + '-' + tbox::Utilities::intToString(
            mpi.getSize(),
            5);
      tbox::plog << "Added case name (" << case_name << ") and nprocs ("
                 << mpi.getSize() << ") to base name -> '"
                 << base_name << "'\n";

      if (!case_name.empty()) {
         tbox::plog << "Added case name (" << case_name << ") and nprocs ("
                    << mpi.getSize() << ") to base name -> '"
                    << base_name << "'\n";
      }

      /*
       * Start logging.
       */
      const std::string log_file_name = base_name + ".log";
      bool log_all_nodes = false;
      log_all_nodes = main_db->getBoolWithDefault("log_all_nodes",
            log_all_nodes);
      if (log_all_nodes) {
         PIO::logAllNodes(log_file_name);
      } else {
         PIO::logOnlyNodeZero(log_file_name);
      }

      plog << "Input database after initialization..." << std::endl;
      input_db->printClassData(plog);



      /*
       * Parameters.  Some of these can be specified by input deck.
       */
      hier::IntVector ghost_cell_width(dim, 2);
      if (main_db->isInteger("ghost_cell_width")) {
         main_db->getIntegerArray("ghost_cell_width", &ghost_cell_width[0], dim.getValue());
      }

      hier::IntVector bad_interval(dim, 2);
      hier::IntVector cut_factor(dim, 1);

      hier::OverlapConnectorAlgorithm oca;



      /*
       * Set up the domain from input.
       */

      hier::BoxContainer input_boxes(main_db->getDatabaseBoxArray("domain_boxes"));

      hier::BoxContainer domain_boxes;
      hier::LocalId local_id(0);
      for (hier::BoxContainer::iterator itr = input_boxes.begin();
           itr != input_boxes.end(); ++itr) {
         itr->setBlockId(hier::BlockId(0));
         domain_boxes.pushBack(hier::Box(*itr, local_id++, 0));
      }

      std::vector<double> xlo(dim.getValue());
      std::vector<double> xhi(dim.getValue());
      for (int i = 0; i < dim.getValue(); ++i) {
         xlo[i] = 0.0;
         xhi[i] = 1.0;
      }
      if (main_db->isDouble("xlo")) {
         main_db->getDoubleArray("xlo", &xlo[0], dim.getValue());
      }
      if (main_db->isDouble("xhi")) {
         main_db->getDoubleArray("xhi", &xhi[0], dim.getValue());
      }

      /*
       * If num_procs_in_tile is given, take the domain_boxes, xlo and xhi
       * to be the size for the (integer) value of num_procs_in_tile.  Scale
       * the problem from there to the number of process running by
       * doubling the dimension starting with the j direction.
       *
       * The number of processes must be a power of 2 times the value
       * of num_procs_in_tile.
       */
      if ( main_db->isInteger("num_procs_in_tile") ) {
         int num_procs_in_tile = main_db->getInteger("num_procs_in_tile");
         int doubling_dir = 1;

         while (num_procs_in_tile < mpi.getSize()) {
            for ( hier::BoxContainer::iterator bi=domain_boxes.begin();
                  bi!=domain_boxes.end(); ++bi ) {
               hier::Box &input_box = *bi;
               input_box.upper()(doubling_dir) += input_box.numberCells(doubling_dir);
            }
            xhi[doubling_dir] += xhi[doubling_dir] - xlo[doubling_dir];
            doubling_dir = (doubling_dir + 1)%dim.getValue();
            num_procs_in_tile *= 2;
            tbox::plog << "num_procs_in_tile = " << num_procs_in_tile << std::endl
                       << domain_boxes.format("IB: ", 2) << std::endl;
         }

         if ( num_procs_in_tile != mpi.getSize() ) {
            TBOX_ERROR("If num_procs_in_tile (" << num_procs_in_tile << ") is given,\n"
                       <<"number of processes (" << mpi.getSize() << ") must be\n"
                       <<"a power-of-2 times the value of num_procs_in_tile.");
         }

      }


      {
         /*
          * Add a dummy PatchData with a big ghost width.
          * GridGeometry forbids increasing the max data ghost width
          * after it starts computing boundary boxes for a patch.  We
          * force it to accept a big ghost width here so that the
          * methods below (particularly those registering tag
          * data) won't crash by asking for more ghost width than what
          * was registered when the first boundary boxes were
          * computed.
          */
      hier::VariableDatabase* vdb =
         hier::VariableDatabase::getDatabase();

      boost::shared_ptr<pdat::CellVariable<int> > dummy_variable(
         new pdat::CellVariable<int>(dim, "DummyVariable"));

      boost::shared_ptr<hier::VariableContext> dummy_context(
         vdb->getContext("DUMMY"));

      vdb->registerVariableAndContext(
         dummy_variable,
         dummy_context,
         hier::IntVector(dim,10));
      }



      /*
       * Create hierarchy.
       */

      boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
         new geom::CartesianGridGeometry(
            "GridGeometry",
            &xlo[0],
            &xhi[0],
            domain_boxes));

      boost::shared_ptr<hier::PatchHierarchy> hierarchy(
         new hier::PatchHierarchy(
            "Hierarchy",
            grid_geometry,
            input_db->getDatabase("PatchHierarchy") ));

      const int max_levels = hierarchy->getMaxNumberOfLevels();

      hier::BoxLevel domain_box_level(
         hier::IntVector(dim, 1),
         grid_geometry,
         tbox::SAMRAI_MPI::getSAMRAIWorld(),
         hier::BoxLevel::GLOBALIZED);
      hier::BoxContainer::iterator domain_boxes_itr(domain_boxes);
      for (int i = 0; i < domain_boxes.size(); ++i, ++domain_boxes_itr) {
         domain_box_level.addBox(hier::Box(*domain_boxes_itr,
               hier::LocalId(i), 0));
      }


      /*
       * Set up the load balancers.
       */



      std::string load_balancer_type =
         main_db->getStringWithDefault("load_balancer_type", "TreeLoadBalancer");




      /*
       * Step 1: Build L0.
       */
      tbox::pout << "\nGenerating L0" << std::endl;

      hier::BoxLevel L0(hier::IntVector(dim, 1), grid_geometry);

      {
         hier::BoxContainer L0_boxes(
            main_db->isDatabase("L0_boxes") ?
            main_db->getDatabaseBoxArray("L0_boxes") : domain_boxes );
         const int boxes_per_proc =
            (L0_boxes.size() + L0.getMPI().getSize()
             - 1) / L0.getMPI().getSize();
         const int my_boxes_start = L0.getMPI().getRank()
            * boxes_per_proc;
         const int my_boxes_stop =
            tbox::MathUtilities<int>::Min(my_boxes_start + boxes_per_proc,
               L0_boxes.size());
         hier::BoxContainer::iterator L0_boxes_itr(L0_boxes);
         for (int i = 0; i < my_boxes_start; ++i) {
            ++L0_boxes_itr;
         }
         for (int i = my_boxes_start; i < my_boxes_stop; ++i, ++L0_boxes_itr) {
            L0.addBox(*L0_boxes_itr, hier::BlockId::zero());
         }
      }

      {
         /*
          * Load balance the L0 BoxLevel, using the domain as its L0.
          *
          * This is not a part of the performance test because does not
          * reflect the load balancer use in real apps.  We just neeed a
          * distributed L0 for the real load balancing performance test.
          */
         hier::Connector L0_to_domain(
            L0,
            domain_box_level,
            hier::IntVector(dim, 2));
         hier::Connector domain_to_L0(
            domain_box_level,
            L0,
            hier::IntVector(dim, 2));
         oca.findOverlaps(L0_to_domain);
         oca.findOverlaps(domain_to_L0);

         boost::shared_ptr<mesh::LoadBalanceStrategy> lb0(
            createLoadBalancer( base_db, load_balancer_type, 0, dim ));

         tbox::plog << "\n\n\ninitial L0 loads:\n";
         mesh::BalanceUtilities::gatherAndReportLoadBalance(
            (double)L0.getLocalNumberOfCells(),
            L0.getMPI());

         tbox::SAMRAI_MPI::getSAMRAIWorld().Barrier();
         lb0->loadBalanceBoxLevel(
            L0,
            L0_to_domain,
            domain_to_L0,
            hierarchy,
            0,
            hier::Connector(),
            hier::Connector(),
            hierarchy->getSmallestPatchSize(0),
            hierarchy->getLargestPatchSize(0),
            domain_box_level,
            bad_interval,
            cut_factor);

         sortNodes(L0,
            domain_to_L0,
            L0_to_domain,
            false,
            true);

         oca.assertOverlapCorrectness(L0_to_domain);
         oca.assertOverlapCorrectness(domain_to_L0);

         L0.cacheGlobalReducedData();

         tbox::plog << "\n\n\nfinal L0 loads:\n";
         mesh::BalanceUtilities::gatherAndReportLoadBalance(
            (double)L0.getLocalNumberOfCells(),
            L0.getMPI());

         outputMetadataL0(
            L0,
            L0.getPersistentOverlapConnectors().
            findOrCreateConnector( L0, ghost_cell_width, true ) );
      }



      hier::BoxLevel L1(dim);
      hier::Connector L1_to_L0;
      hier::Connector L0_to_L1;
      hier::Connector L1_to_L1;

      hier::BoxLevel L2(dim);
      hier::Connector L2_to_L1;
      hier::Connector L1_to_L2;
      hier::Connector L2_to_L2;



      if ( max_levels > 1 ) {
         /*
          * Step 2: Build L1.
          */
         tbox::pout << "\nGenerating L1" << std::endl;


         const int coarser_ln = 0;
         const int finer_ln = coarser_ln + 1;

         // Get the prebalanced L1:
         const hier::IntVector required_connector_width =
            hierarchy->getRequiredConnectorWidth(coarser_ln, finer_ln);
         const hier::IntVector min_size = hier::IntVector::ceilingDivide(
            hierarchy->getSmallestPatchSize(finer_ln), hierarchy->getRatioToCoarserLevel(finer_ln) );

         std::string L1_box_gen_method =
            main_db->getStringWithDefault("L1_box_gen_method", "PrebalanceByUserBoxes");
         generatePrebalance(
            L1,
            L0,
            L0_to_L1,
            L1_to_L0,
            L1_box_gen_method,
            *main_db,
            hierarchy,
            coarser_ln,
            min_size,
            required_connector_width );


         // Output metadata before balancing L1.
         outputMetadataBefore( L0_to_L1, L1_to_L0, "L0", "L1pre" );

         if ( L1.getGlobalNumberOfBoxes() == 0 ) {
            TBOX_ERROR("Level " << finer_ln << " box generator resulted in no boxes.");
         }

         boost::shared_ptr<mesh::LoadBalanceStrategy> lb1(
            createLoadBalancer( base_db, load_balancer_type, 1 , dim));

         mesh::BalanceUtilities::gatherAndReportLoadBalance(
            (double)L1.getLocalNumberOfCells(),
            L1.getMPI());

         tbox::SAMRAI_MPI::getSAMRAIWorld().Barrier();
         // Load balance L1.
         lb1->loadBalanceBoxLevel(
            L1,
            L1_to_L0,
            L0_to_L1,
            hierarchy,
            1,
            hier::Connector(),
            hier::Connector(),
            hier::IntVector::ceilingDivide(hierarchy->getSmallestPatchSize(1), hierarchy->getRatioToCoarserLevel(1)),
            hier::IntVector::ceilingDivide(hierarchy->getLargestPatchSize(1), hierarchy->getRatioToCoarserLevel(1)),
            domain_box_level,
            bad_interval,
            cut_factor);

         oca.assertOverlapCorrectness(L1_to_L0);
         oca.assertOverlapCorrectness(L0_to_L1);

         sortNodes(L1,
                   L0_to_L1,
                   L1_to_L0,
                   false,
                   true);

         // Refine L1.
         if ( hierarchy->getRatioToCoarserLevel(1) != zero_vec ) {
            refineHead(
               L1,
               L0_to_L1,
               L1_to_L0,
               hierarchy->getRatioToCoarserLevel(1) );
         }

         // Get the L1_to_L1 for edge statistics.
         oca.bridge(
            L1_to_L1,
            L1_to_L0,
            L0_to_L1,
            L1_to_L0,
            L0_to_L1);

         // Output metadata after balancing L1.
         outputMetadataAfter( L0_to_L1, L1_to_L0, L1_to_L1, "L0", "L1post" );

         mesh::BalanceUtilities::gatherAndReportLoadBalance(
            (double)L1.getLocalNumberOfCells(),
            L1.getMPI());

      }



      if ( max_levels > 2 ) {
         /*
          * Step 3: Build L2.
          */
         tbox::pout << "\nGenerating L2" << std::endl;

         const int coarser_ln = 1;
         const int finer_ln = coarser_ln + 1;

         // Get the prebalanced L2:
         const hier::IntVector required_connector_width =
            hierarchy->getRequiredConnectorWidth(coarser_ln, finer_ln);
         const hier::IntVector min_size = hier::IntVector::ceilingDivide(
            hierarchy->getSmallestPatchSize(finer_ln), hierarchy->getRatioToCoarserLevel(finer_ln) );

         std::string L2_box_gen_method =
            main_db->getStringWithDefault("L2_box_gen_method", "PrebalanceByShrinkingLevel");
         generatePrebalance(
            L2,
            L1,
            L1_to_L2,
            L2_to_L1,
            L2_box_gen_method,
            *main_db,
            hierarchy,
            coarser_ln,
            min_size,
            required_connector_width );


         // Output metadata before balancing L2.
         outputMetadataBefore( L1_to_L2, L2_to_L1, "L1", "L2pre" );

         if ( L2.getGlobalNumberOfBoxes() == 0 ) {
            TBOX_ERROR("Level " << finer_ln << " box generator resulted in no boxes.");
         }

         boost::shared_ptr<mesh::LoadBalanceStrategy> lb2(
            createLoadBalancer( base_db, load_balancer_type, 2 , dim));

         mesh::BalanceUtilities::gatherAndReportLoadBalance(
            (double)L2.getLocalNumberOfCells(),
            L2.getMPI());

         tbox::SAMRAI_MPI::getSAMRAIWorld().Barrier();
         // Load balance L2.
         lb2->loadBalanceBoxLevel(
            L2,
            L2_to_L1,
            L1_to_L2,
            hierarchy,
            1,
            hier::Connector(),
            hier::Connector(),
            hier::IntVector::ceilingDivide(hierarchy->getSmallestPatchSize(2), hierarchy->getRatioToCoarserLevel(2)),
            hier::IntVector::ceilingDivide(hierarchy->getLargestPatchSize(2), hierarchy->getRatioToCoarserLevel(2)),
            domain_box_level,
            bad_interval,
            cut_factor);

         oca.assertOverlapCorrectness(L2_to_L1);
         oca.assertOverlapCorrectness(L1_to_L2);

         sortNodes(L2,
                   L1_to_L2,
                   L2_to_L1,
                   false,
                   true);

         // Refine L2.
         if ( hierarchy->getRatioToCoarserLevel(2) != zero_vec ) {
            refineHead(
               L2,
               L1_to_L2,
               L2_to_L1,
               hierarchy->getRatioToCoarserLevel(2) );
         }

         // Get the L2_to_L2 for edge statistics.
         oca.bridge(
            L2_to_L2,
            L2_to_L1,
            L1_to_L2,
            L2_to_L1,
            L1_to_L2);

         // Output metadata after balancing L2.
         outputMetadataAfter( L1_to_L2, L2_to_L1, L2_to_L2, "L1", "L2post" );

         mesh::BalanceUtilities::gatherAndReportLoadBalance(
            (double)L2.getLocalNumberOfCells(),
            L2.getMPI());

      }




      bool write_visit =
         main_db->getBoolWithDefault("write_visit", false);
      if ( write_visit ) {
#ifdef HAVE_HDF5
         hierarchy->makeNewPatchLevel(0, L0);
         if ( hierarchy->getMaxNumberOfLevels() > 1 ) {
            hierarchy->makeNewPatchLevel(1, L1);
         }
         if ( hierarchy->getMaxNumberOfLevels() > 2 ) {
            hierarchy->makeNewPatchLevel(2, L2);
         }

         if ((dim == tbox::Dimension(2)) || (dim == tbox::Dimension(3))) {
            /*
             * Create the VisIt data writer.
             * Write the plot file.
             */
            DerivedVisOwnerData owner_writer;
            const std::string visit_filename = base_name + ".visit";
            appu::VisItDataWriter visit_data_writer(dim,
                                                    "VisIt Writer",
                                                    visit_filename);
            visit_data_writer.registerDerivedPlotQuantity("Owner",
                                                          "SCALAR",
                                                          &owner_writer);
            visit_data_writer.writePlotData(hierarchy, 0);
         }
#else
         TBOX_WARNING("main: You set write_visit to TRUE,\n"
                      << "but VisIt dumps are not supported due to\n"
                      << "not having configured with HDF5.\n");
#endif
      }

   }

   /*
    * Output timer results.
    */
   tbox::TimerManager::getManager()->print(tbox::plog);


   /*
    * Print input database again to fully show usage.
    */
   plog << "Input database after running..." << std::endl;
   tbox::InputManager::getManager()->getInputDatabase()->printClassData(plog);

   tbox::pout << "\nPASSED:  treelb" << std::endl;

   /*
    * Exit properly by shutting down services in correct order.
    */
   tbox::plog << "\nShutting down..." << std::endl;

   /*
    * Shut down.
    */
   SAMRAIManager::shutdown();
   SAMRAIManager::finalize();

   if (fail_count == 0) {
      SAMRAI_MPI::finalize();
   } else {
      std::cout << "Process " << std::setw(5) << rank << " aborting."
                << std::endl;
      tbox::Utilities::abort("Aborting due to nonzero fail count",
         __FILE__, __LINE__);
   }

   return fail_count;
}



/*
****************************************************************************
* Output data for L0.
****************************************************************************
*/
void outputMetadataL0(
   const hier::BoxLevel &L0,
   const hier::Connector &L0_to_L0,
   int level_output_depth,
   int connector_output_depth )
{
   const std::string L0_name("L0");
   const std::string L0_to_L0_name(L0_name + "_to_" + L0_name );
   const std::string arrow("-> ");

   L0.cacheGlobalReducedData();

   tbox::plog << "\n\n\n" << L0_name << ":\n";

   hier::BoxLevelStatistics L0_stats(L0);
   tbox::plog << L0_name << " stats:\n";
   L0_stats.printBoxStats(tbox::plog, L0_name + arrow);
   tbox::plog << L0_name << ":\n" << L0.format( L0_name + arrow, level_output_depth );

   hier::ConnectorStatistics L0_L0_stats(L0_to_L0);
   tbox::plog << L0_to_L0_name << " neighbor stats:\n";
   L0_L0_stats.printNeighborStats(tbox::plog, L0_to_L0_name + arrow );
   tbox::plog << L0_to_L0_name << ":\n" << L0_to_L0.format( L0_to_L0_name + arrow, connector_output_depth );

   return;
}



/*
****************************************************************************
* Output "before" data.
****************************************************************************
*/
void outputMetadataBefore(
   const hier::Connector &La_to_Lb,
   const hier::Connector &Lb_to_La,
   const std::string &La_name,
   const std::string &Lb_name,
   int level_output_depth,
   int connector_output_depth )
{
   const hier::BoxLevel &La = La_to_Lb.getBase();
   const hier::BoxLevel &Lb = La_to_Lb.getHead();

   const std::string La_to_Lb_name(La_name + "_to_" + Lb_name );
   const std::string Lb_to_La_name(Lb_name + "_to_" + La_name );
   const std::string arrow("-> ");

   La.cacheGlobalReducedData();
   Lb.cacheGlobalReducedData();

   tbox::plog << "\n\n\nBefore balancing " << Lb_name << ":\n";

   hier::BoxLevelStatistics La_stats(La);
   tbox::plog << La_name << " stats:\n";
   La_stats.printBoxStats(tbox::plog, La_name + arrow);
   tbox::plog << La_name << ":\n" << La.format( La_name + arrow, level_output_depth );

   hier::BoxLevelStatistics Lb_stats(Lb);
   tbox::plog << Lb_name << " stats:\n";
   Lb_stats.printBoxStats(tbox::plog, Lb_name + arrow);
   tbox::plog << Lb_name << ":\n" << Lb.format( Lb_name + arrow, level_output_depth );

   hier::ConnectorStatistics Lb_La_stats(Lb_to_La);
   tbox::plog << Lb_to_La_name << " neighbor stats:\n";
   Lb_La_stats.printNeighborStats(tbox::plog, Lb_to_La_name + arrow );
   tbox::plog << Lb_to_La_name << ":\n" << Lb_to_La.format( Lb_to_La_name + arrow, connector_output_depth );

   hier::ConnectorStatistics La_Lb_stats(La_to_Lb);
   tbox::plog << La_to_Lb_name << " neighbor stats:\n";
   La_Lb_stats.printNeighborStats(tbox::plog, La_to_Lb_name + arrow );
   tbox::plog << La_to_Lb_name << ":\n" << La_to_Lb.format( La_to_Lb_name + arrow, connector_output_depth );

   return;
}



/*
****************************************************************************
* Output "after" data.
****************************************************************************
*/
void outputMetadataAfter(
   const hier::Connector &La_to_Lb,
   const hier::Connector &Lb_to_La,
   const hier::Connector &Lb_to_Lb,
   const std::string &La_name,
   const std::string &Lb_name,
   int level_output_depth,
   int connector_output_depth )
{
   const hier::BoxLevel &Lb = La_to_Lb.getHead();

   const std::string La_to_Lb_name(La_name + "_to_" + Lb_name );
   const std::string Lb_to_La_name(Lb_name + "_to_" + La_name );
   const std::string Lb_to_Lb_name(Lb_name + "_to_" + Lb_name );
   const std::string arrow("-> ");

   Lb.cacheGlobalReducedData();

   tbox::plog << "\n\n\nAfter balancing " << Lb_name << ":\n";

   hier::BoxLevelStatistics Lb_stats(Lb);
   tbox::plog << Lb_name << " stats:\n";
   Lb_stats.printBoxStats(tbox::plog, Lb_name + arrow);
   tbox::plog << Lb_name << ":\n" << Lb.format( Lb_name + arrow, level_output_depth );

   hier::ConnectorStatistics Lb_Lb_stats(Lb_to_Lb);
   tbox::plog << Lb_to_Lb_name << " neighbor stats:\n";
   Lb_Lb_stats.printNeighborStats(tbox::plog, Lb_to_Lb_name + arrow );
   tbox::plog << Lb_to_Lb_name << ":\n" << Lb_to_Lb.format( Lb_to_Lb_name + arrow, connector_output_depth );

   hier::ConnectorStatistics Lb_La_stats(Lb_to_La);
   tbox::plog << Lb_to_La_name << " neighbor stats:\n";
   Lb_La_stats.printNeighborStats(tbox::plog, Lb_to_La_name + arrow );
   tbox::plog << Lb_to_La_name << ":\n" << Lb_to_La.format( Lb_to_La_name + arrow, connector_output_depth );

   hier::ConnectorStatistics La_Lb_stats(La_to_Lb);
   tbox::plog << La_to_Lb_name << " neighbor stats:\n";
   La_Lb_stats.printNeighborStats(tbox::plog, La_to_Lb_name + arrow );
   tbox::plog << La_to_Lb_name << ":\n" << La_to_Lb.format( La_to_Lb_name + arrow, connector_output_depth );

   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void generatePrebalanceByUserShells(
   hier::BoxLevel& L1,
   hier::Connector& L0_to_L1,
   hier::Connector& L1_to_L0,
   const hier::BoxLevel& L0,
   const boost::shared_ptr<tbox::Database>& database,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const hier::IntVector& min_size,
   const hier::IntVector& connector_width )
{

   const tbox::Dimension dim(hierarchy->getDim());
   const hier::IntVector &zero_vec(hier::IntVector::getZero(dim));
   const hier::IntVector &one_vec(hier::IntVector::getOne(dim));

   /*
    * Starting at shell origin, tag cells with centroids
    * at radii[0]<r<radii[1], radii[2]<r<radii[3], and so on.
    */
   tbox::Array<double> radii;
   double efficiency_tol = 0.75;
   double combine_tol = 0.75;

   std::vector<double> r0(dim.getValue());
   for (int d = 0; d < dim.getValue(); ++d) r0[d] = 0;

   boost::shared_ptr<tbox::Database> abr_db;
   if (database) {
      efficiency_tol = database->getDoubleWithDefault("efficiency_tol",
            efficiency_tol);
      combine_tol = database->getDoubleWithDefault("combine_tol", combine_tol);
      if (database->isDouble("r0")) {
         database->getDoubleArray("r0", &r0[0], dim.getValue());
      }
      if (database->isDouble("radii")) {
         radii = database->getDoubleArray("radii");
      }
      abr_db = database->getDatabaseWithDefault("BergerRigoutsos", abr_db);
      TBOX_ASSERT(radii.size() % 2 == 0);
   }

   const int tag_val = 1;

   hier::VariableDatabase* vdb =
      hier::VariableDatabase::getDatabase();
   boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
      hierarchy->getGridGeometry(),
      boost::detail::dynamic_cast_tag());

   boost::shared_ptr<hier::PatchLevel> tag_level(
      new hier::PatchLevel(L0,
         grid_geometry,
         vdb->getPatchDescriptor()));

   boost::shared_ptr<pdat::CellVariable<int> > tag_variable(
      new pdat::CellVariable<int>(
         dim,
         "UserShellsTagVariable"));

   boost::shared_ptr<hier::VariableContext> default_context(
      vdb->getContext("TagVariable"));

   const int tag_id = vdb->registerVariableAndContext(
         tag_variable,
         default_context,
         hier::IntVector::getZero(dim));

   tag_level->allocatePatchData(tag_id);

   const double* xlo = grid_geometry->getXLower();
   const double* h = grid_geometry->getDx();
   std::vector<double> r(dim.getValue());
   for (hier::PatchLevel::iterator pi(tag_level->begin());
        pi != tag_level->end(); ++pi) {
      const boost::shared_ptr<hier::Patch>& patch = *pi;

      pdat::NodeData<double> node_tag_data(patch->getBox(), 1, zero_vec);
      pdat::NodeData<int>::iterator niend(node_tag_data.getGhostBox(), false);
      for (pdat::NodeData<int>::iterator ni(node_tag_data.getGhostBox(), true);
           ni != niend; ++ni) {
         const pdat::NodeIndex& idx = *ni;
         double rr = 0;
         for (int d = 0; d < dim.getValue(); ++d) {
            r[d] = xlo[d] + h[d] * idx(d) - r0[d];
            rr += r[d] * r[d];
         }
         rr = sqrt(rr);
         for (int i = 0; i < radii.size(); i += 2) {
            if (radii[i] < rr && rr < radii[i + 1]) {
               (node_tag_data)(idx) = tag_val;
               break;
            }
         }
      }

      boost::shared_ptr<pdat::CellData<int> > tag_data(
         patch->getPatchData(tag_id),
         boost::detail::dynamic_cast_tag());

      tag_data->getArrayData().fillAll(0);

      const hier::BlockId& block_id = patch->getBox().getBlockId();

      pdat::CellData<int>::iterator ciend(tag_data->getGhostBox(), false);
      for (pdat::CellData<int>::iterator ci(tag_data->getGhostBox(), true);
           ci != ciend; ++ci) {
         const pdat::CellIndex& cid = *ci;

         // Loop through nodes of cell cid.  Tag cell if node is tagged.
         const hier::Box cell_box(cid,cid, block_id);
         pdat::NodeIterator node_itr_end(cell_box, false);
         for ( pdat::NodeIterator node_itr(cell_box, true);
               node_itr != node_itr_end; ++node_itr ) {
            if ( node_tag_data(*node_itr) == tag_val ) {
               (*tag_data)(cid) = tag_val;
               break;
            }
         }
      }

   }

   mesh::BergerRigoutsos abr(dim, abr_db);
   abr.setMPI(L0.getMPI());
   abr.findBoxesContainingTags(
      L1,
      L0_to_L1,
      L1_to_L0,
      tag_level,
      tag_id,
      tag_val,
      L0.getGlobalBoundingBox(0),
      min_size,
      efficiency_tol,
      combine_tol,
      connector_width,
      hier::BlockId::zero(),
      hier::LocalId(0));

   /*
    * The clustering step generated Connectors to/from the temporary
    * tag_level->getBoxLevel(), which is not the same as the
    * L0 BoxLevel.  We need to reset the Connectors to use
    * the L0 instead.
    */
   L0_to_L1.setBase(L0);
   L0_to_L1.setHead(L1, true);
   L1_to_L0.setBase(L1);
   L1_to_L0.setHead(L0, true);


   /*
    * Make L1 nest inside L0 by one cell.
    */
   hier::BoxLevel L1nested(dim);
   hier::Connector L1_to_L1nested;
   hier::BoxLevelConnectorUtils blcu;
   blcu.computeInternalParts( L1nested,
                              L1_to_L1nested,
                              L1_to_L0,
                              -one_vec,
                              grid_geometry->getDomainSearchTree() );
   hier::MappingConnectorAlgorithm mca;
   mca.modify( L0_to_L1,
               L1_to_L0,
               L1_to_L1nested,
               &L1,
               &L1nested );

   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void generatePrebalanceByShrinkingLevel(
   hier::BoxLevel& L2,
   hier::Connector& L1_to_L2,
   hier::Connector& L2_to_L1,
   const hier::BoxLevel& L1,
   const boost::shared_ptr<tbox::Database>& database,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int coarser_ln,
   const hier::IntVector& min_size,
   const hier::IntVector& connector_width )
{

   const tbox::Dimension dim(hierarchy->getDim());

   boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
      hierarchy->getGridGeometry(),
      boost::detail::dynamic_cast_tag());

   // Parameters set by database, with defaults.
   double efficiency_tol = 1.00;
   double combine_tol = 1.00;
   hier::IntVector shrink_width(dim, 2);
   boost::shared_ptr<tbox::Database> abr_db;

   if (database) {

      efficiency_tol = database->getDoubleWithDefault("efficiency_tol",
            efficiency_tol);

      combine_tol = database->getDoubleWithDefault("combine_tol", combine_tol);

      abr_db = database->getDatabaseWithDefault("BergerRigoutsos", abr_db);

      database->getIntegerArray("shrink_width", &shrink_width[0], dim.getValue());
   }


   hier::BoxLevel L1tags(dim);
   hier::Connector L1_to_L1tags;
   const hier::Connector &L1_to_L1 =
      L1.getPersistentOverlapConnectors().findOrCreateConnector(
         L1,
         shrink_width );

   hier::BoxLevelConnectorUtils blcu;
   blcu.computeInternalParts( L1tags,
                              L1_to_L1tags,
                              L1_to_L1,
                              -shrink_width,
                              grid_geometry->getDomainSearchTree() );
   tbox::plog << "L1_to_L1tags:\n" << L1_to_L1tags.format("L1->L1tags: ", 2);


   const int tag_val = 1;

   hier::VariableDatabase* vdb =
      hier::VariableDatabase::getDatabase();

   boost::shared_ptr<hier::PatchLevel> tag_level(
      new hier::PatchLevel(
         L1,
         grid_geometry,
         vdb->getPatchDescriptor()));

   boost::shared_ptr<pdat::CellVariable<int> > tag_variable(
      new pdat::CellVariable<int>(
         dim,
         "ShrinkingLevelTagVariable"));

   boost::shared_ptr<hier::VariableContext> default_context(
      vdb->getContext("TagVariable"));

   const int tag_id = vdb->registerVariableAndContext(
         tag_variable,
         default_context,
         hier::IntVector::getZero(dim));

   tag_level->allocatePatchData(tag_id);

   for (hier::PatchLevel::iterator pi(tag_level->begin());
        pi != tag_level->end(); ++pi) {

      const boost::shared_ptr<hier::Patch>& patch = *pi;
      boost::shared_ptr<pdat::CellData<int> > tag_data(
         patch->getPatchData(tag_id),
         boost::detail::dynamic_cast_tag());

      tag_data->getArrayData().fillAll(0);

      if ( !L1_to_L1tags.hasNeighborSet(patch->getBox().getId()) ) {
         tag_data->getArrayData().fillAll(1);
      }
      else {
         hier::Connector::ConstNeighborhoodIterator ni =
            L1_to_L1tags.find(patch->getBox().getId());

         for ( hier::Connector::ConstNeighborIterator na = L1_to_L1tags.begin(ni);
               na != L1_to_L1tags.end(ni); ++na ) {

            const hier::Box &tag_box = *na;
            tag_data->getArrayData().fillAll(1, tag_box);

         }
      }

   }

   mesh::BergerRigoutsos abr(dim, abr_db);
   abr.setMPI(L1.getMPI());
   abr.findBoxesContainingTags(
      L2,
      L1_to_L2,
      L2_to_L1,
      tag_level,
      tag_id,
      tag_val,
      L1.getGlobalBoundingBox(0),
      min_size,
      efficiency_tol,
      combine_tol,
      connector_width,
      hier::BlockId::zero(),
      hier::LocalId(0));


   /*
    * The clustering step generated Connectors to/from the temporary
    * tag_level->getBoxLevel(), which is not the same as the
    * L1 BoxLevel.  We need to reset the Connectors to use
    * the L1 instead.
    */
   L1_to_L2.setBase(L1);
   L1_to_L2.setHead(L2, true);
   L2_to_L1.setBase(L2);
   L2_to_L1.setHead(L1, true);


   /*
    * Make L2 nest inside L1 by shrink_width.
    */
   const hier::IntVector nesting_width(dim, hierarchy->getProperNestingBuffer(coarser_ln));
   hier::BoxLevel L2nested(dim);
   hier::Connector L2_to_L2nested;
   blcu.computeInternalParts( L2nested,
                              L2_to_L2nested,
                              L2_to_L1,
                              -nesting_width,
                              grid_geometry->getDomainSearchTree() );
   hier::MappingConnectorAlgorithm mca;
   mca.modify( L1_to_L2,
               L2_to_L1,
               L2_to_L2nested,
               &L2,
               &L2nested );

   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void generatePrebalanceBySinusoidalFront(
   hier::BoxLevel& L2,
   hier::Connector& L1_to_L2,
   hier::Connector& L2_to_L1,
   const hier::BoxLevel& L1,
   const boost::shared_ptr<tbox::Database>& database,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int coarser_ln,
   const hier::IntVector& min_size,
   const hier::IntVector& connector_width )
{

   const tbox::Dimension dim(hierarchy->getDim());

   boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
      hierarchy->getGridGeometry(),
      boost::detail::dynamic_cast_tag());

   // Parameters set by database, with defaults.
   double efficiency_tol = 0.70;
   double combine_tol = 0.70;
   hier::IntVector tag_buffer(dim, 2);
   boost::shared_ptr<tbox::Database> abr_db; // BergerRigoutsos database.
   boost::shared_ptr<tbox::Database> sft_db; // SinusoidalFrontTagger database.

   if (database) {

      efficiency_tol = database->getDoubleWithDefault("efficiency_tol",
            efficiency_tol);

      combine_tol = database->getDoubleWithDefault("combine_tol", combine_tol);

      abr_db = database->getDatabaseWithDefault("BergerRigoutsos", abr_db);

      sft_db = database->getDatabaseWithDefault("SinusoidalFrontTagger", sft_db);

      if ( database->isInteger("tag_buffer") ) {
         database->getIntegerArray("tag_buffer", &tag_buffer[0], dim.getValue());
      }
   }




   const int tag_val = 1;

   hier::VariableDatabase* vdb =
      hier::VariableDatabase::getDatabase();

   boost::shared_ptr<hier::PatchLevel> tag_level(
      new hier::PatchLevel(
         L1,
         grid_geometry,
         vdb->getPatchDescriptor()));

   boost::shared_ptr<pdat::CellVariable<int> > tag_variable(
      new pdat::CellVariable<int>(
         dim,
         "SinusoidalFrontTagVariable"));

   boost::shared_ptr<hier::VariableContext> default_context(
      vdb->getContext("TagVariable"));

   const int tag_id = vdb->registerVariableAndContext(
         tag_variable,
         default_context,
         tag_buffer );

   tag_level->allocatePatchData(tag_id);

   SinusoidalFrontTagger sinusoidal_front_tagger(
      "SinusoidalFrontTagger",
      dim,
      sft_db.get() );
   sinusoidal_front_tagger.resetHierarchyConfiguration(hierarchy, 0, 1);

   for (hier::PatchLevel::iterator pi(tag_level->begin());
        pi != tag_level->end(); ++pi) {

      const boost::shared_ptr<hier::Patch>& patch = *pi;

      boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
         patch->getPatchGeometry(),
         boost::detail::dynamic_cast_tag());

      boost::shared_ptr<pdat::CellData<int> > tag_data(
         patch->getPatchData(tag_id),
         boost::detail::dynamic_cast_tag());

      sinusoidal_front_tagger.computeFrontsData(
         NULL /* distance data */,
         tag_data.get(),
         tag_buffer,
         patch_geom->getXLower(),
         patch_geom->getDx(),
         0.0 );

      // tbox::plog << "Tag data for patch " << patch->getBox() << ":\n";
      // tag_data->print(tag_data->getGhostBox(),0,tbox::plog);

   }

   mesh::BergerRigoutsos abr(dim, abr_db);
   abr.setMPI(L1.getMPI());
   abr.findBoxesContainingTags(
      L2,
      L1_to_L2,
      L2_to_L1,
      tag_level,
      tag_id,
      tag_val,
      L1.getGlobalBoundingBox(0),
      min_size,
      efficiency_tol,
      combine_tol,
      connector_width,
      hier::BlockId::zero(),
      hier::LocalId(0));


   /*
    * The clustering step generated Connectors to/from the temporary
    * tag_level->getBoxLevel(), which is not the same as the
    * L1 BoxLevel.  We need to reset the Connectors to use
    * the L1 instead.
    */
   L1_to_L2.setBase(L1);
   L1_to_L2.setHead(L2, true);
   L2_to_L1.setBase(L2);
   L2_to_L1.setHead(L1, true);


   /*
    * Make L2 nest inside L1 by nesting_width.
    */
   const hier::IntVector nesting_width(dim, hierarchy->getProperNestingBuffer(coarser_ln));
   hier::BoxLevel L2nested(dim);
   hier::Connector L2_to_L2nested;
   hier::BoxLevelConnectorUtils blcu;
   blcu.computeInternalParts( L2nested,
                              L2_to_L2nested,
                              L2_to_L1,
                              -nesting_width,
                              grid_geometry->getDomainSearchTree() );
   hier::MappingConnectorAlgorithm mca;
   mca.modify( L1_to_L2,
               L2_to_L1,
               L2_to_L2nested,
               &L2,
               &L2nested );

   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void generatePrebalanceByUserBoxes(
   hier::BoxLevel& L1,
   hier::Connector& L0_to_L1,
   hier::Connector& L1_to_L0,
   const hier::BoxLevel& L0,
   boost::shared_ptr<tbox::Database> database,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const hier::IntVector& min_size,
   const hier::IntVector& connector_width )
{
   NULL_USE(min_size);

   const tbox::Dimension& dim(hierarchy->getDim());

   hier::BoxContainer prebalance_boxes(database->getDatabaseBoxArray("prebalance_boxes"));
   tbox::Array<int> initial_owners(1);
   initial_owners[0] = 0;
   initial_owners = database->getIntegerArray("initial_owners");

   L1.initialize(hier::IntVector(dim, 1),
      hierarchy->getGridGeometry(),
      L0.getMPI());
   hier::BoxContainer::iterator prebalance_boxes_itr(prebalance_boxes);
   for (int i = 0; i < prebalance_boxes.size(); ++i, ++prebalance_boxes_itr) {
      const int owner = i % initial_owners.size();
      if (owner == L1.getMPI().getRank()) {
         prebalance_boxes_itr->setBlockId(hier::BlockId(0));
         L1.addBox(hier::Box(*prebalance_boxes_itr,
               hier::LocalId(i), owner));
      }
   }

   // Generate the balance<===>L0 Connectors.
   L1_to_L0.setBase(L1);
   L1_to_L0.setHead(L0);
   L1_to_L0.setWidth(connector_width, true);
   L0_to_L1.setBase(L0);
   L0_to_L1.setHead(L1);
   L0_to_L1.setWidth(connector_width, true);
   hier::OverlapConnectorAlgorithm oca;
   oca.findOverlaps(L1_to_L0);
   oca.findOverlaps(L0_to_L1);

   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void sortNodes(
   hier::BoxLevel& new_box_level,
   hier::Connector& tag_to_new,
   hier::Connector& new_to_tag,
   bool sort_by_corners,
   bool sequentialize_global_indices)
{
   const hier::MappingConnectorAlgorithm mca;

   hier::Connector sorting_map;
   hier::BoxLevel seq_box_level(new_box_level.getDim());
   hier::BoxLevelConnectorUtils dlbg_edge_utils;
   dlbg_edge_utils.makeSortingMap(
      seq_box_level,
      sorting_map,
      new_box_level,
      sort_by_corners,
      sequentialize_global_indices);

   mca.modify(tag_to_new,
      new_to_tag,
      sorting_map,
      &new_box_level);

   return;
}




/*
 ***********************************************************************
 ***********************************************************************
 */
void refineHead(
   hier::BoxLevel& head,
   hier::Connector& ref_to_head,
   hier::Connector& head_to_ref,
   const hier::IntVector &refinement_ratio )
{
   head.refineBoxes(
      head,
      refinement_ratio,
      head.getRefinementRatio()*refinement_ratio);
   head.finalize();

   const hier::IntVector& head_to_ref_width =
      refinement_ratio * head_to_ref.getConnectorWidth();
   head_to_ref.setBase(head);
   head_to_ref.setWidth(head_to_ref_width, true);

   ref_to_head.setHead(head, true);
   ref_to_head.refineLocalNeighbors(refinement_ratio);

   return;
}




/*
 ***********************************************************************
 ***********************************************************************
 */
void generatePrebalance(
   hier::BoxLevel &Lb,
   hier::BoxLevel &La,
   hier::Connector &La_to_Lb,
   hier::Connector &Lb_to_La,
   const std::string &box_gen_method,
   tbox::Database &main_db,
   const boost::shared_ptr<hier::PatchHierarchy> &hierarchy,
   int coarser_ln,
   const hier::IntVector &min_size,
   const hier::IntVector &required_connector_width )
{
   boost::shared_ptr<tbox::Database> box_gen_db(
      main_db.getDatabaseWithDefault(
         box_gen_method, boost::shared_ptr<Database>()));

   if (box_gen_method == "PrebalanceByUserBoxes") {
      generatePrebalanceByUserBoxes(
         Lb,
         La_to_Lb,
         Lb_to_La,
         La,
         box_gen_db,
         hierarchy,
         min_size,
         required_connector_width );
   } else if (box_gen_method == "PrebalanceByUserShells") {
      generatePrebalanceByUserShells(
         Lb,
         La_to_Lb,
         Lb_to_La,
         La,
         box_gen_db,
         hierarchy,
         min_size,
         required_connector_width );
   } else if (box_gen_method == "PrebalanceBySinusoidalFront") {
      generatePrebalanceBySinusoidalFront(
         Lb,
         La_to_Lb,
         Lb_to_La,
         La,
         box_gen_db,
         hierarchy,
         coarser_ln,
         min_size,
         required_connector_width );
   } else if (box_gen_method == "PrebalanceByShrinkingLevel") {
      generatePrebalanceByShrinkingLevel(
         Lb,
         La_to_Lb,
         Lb_to_La,
         La,
         box_gen_db,
         hierarchy,
         coarser_ln,
         min_size,
         required_connector_width );
   } else {
      TBOX_ERROR("Bad box_gen_method: '" << box_gen_method << "'");
   }

   return;
}

boost::shared_ptr<mesh::LoadBalanceStrategy>
createLoadBalancer(
   boost::shared_ptr<tbox::Database> &input_db,
   const std::string &lb_type,
   int ln,
   const tbox::Dimension &dim )
{

   if (lb_type == "TreeLoadBalancer") {

      boost::shared_ptr<mesh::TreeLoadBalancer> tree_lb(
         new mesh::TreeLoadBalancer(
            dim,
            std::string("mesh::TreeLoadBalancer") + tbox::Utilities::intToString(ln),
            input_db->getDatabaseWithDefault("TreeLoadBalancer",
                                             boost::shared_ptr<tbox::Database>())));
      tree_lb->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());
      return tree_lb;

   }else if (lb_type == "TreeLoadBalancerOld") {

      boost::shared_ptr<mesh::TreeLoadBalancerOld> tree_lb(
         new mesh::TreeLoadBalancerOld(
            dim,
            std::string("mesh::TreeLoadBalancerOld") + tbox::Utilities::intToString(ln),
            input_db->getDatabaseWithDefault("TreeLoadBalancerOld",
                                             boost::shared_ptr<tbox::Database>())));
      tree_lb->setSAMRAI_MPI(tbox::SAMRAI_MPI::getSAMRAIWorld());
      return tree_lb;

   } else if (lb_type == "ChopAndPackLoadBalancer") {

      boost::shared_ptr<mesh::ChopAndPackLoadBalancer> cap_lb(
         new mesh::ChopAndPackLoadBalancer(
            dim,
            std::string("mesh::ChopAndPackLoadBalancer") + tbox::Utilities::intToString(ln),
            input_db->getDatabaseWithDefault("ChopAndPackLoadBalancer",
                                             boost::shared_ptr<tbox::Database>())));
      return cap_lb;

   }
   else {
      TBOX_ERROR(
         "Missing or bad load_balancer specification in Main database.\n"
         << "Specify load_balancer_type = STRING, where STRING can be\n"
         << "\"ChopAndPackLoadBalancer\" or \"TreeLoadBalancer\".");
   }

   return boost::shared_ptr<mesh::LoadBalanceStrategy>();
}
