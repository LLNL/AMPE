/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Asynchronous Berger-Rigoutsos algorithm wrapper
 *
 ************************************************************************/
#ifndef included_mesh_BergerRigoutsos_C
#define included_mesh_BergerRigoutsos_C

#include <stdlib.h>

#include "SAMRAI/mesh/BergerRigoutsos.h"

#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/mesh/BergerRigoutsosNode.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/TimerManager.h"

namespace SAMRAI {
namespace mesh {

boost::shared_ptr<tbox::Timer> BergerRigoutsos::t_barrier_before;
boost::shared_ptr<tbox::Timer> BergerRigoutsos::t_barrier_after;
boost::shared_ptr<tbox::Timer> BergerRigoutsos::t_find_boxes_with_tags;
boost::shared_ptr<tbox::Timer> BergerRigoutsos::t_run_abr;
boost::shared_ptr<tbox::Timer> BergerRigoutsos::t_global_reductions;
boost::shared_ptr<tbox::Timer> BergerRigoutsos::t_sort_output_nodes;

tbox::StartupShutdownManager::Handler
BergerRigoutsos::s_initialize_finalize_handler(
   BergerRigoutsos::initializeCallback,
   0,
   0,
   BergerRigoutsos::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 ************************************************************************
 * Constructor stores parameters of options for ussing
 * the asynchronous Berger-Rigoutsos implementation.
 ************************************************************************
 */
BergerRigoutsos::BergerRigoutsos(
   const tbox::Dimension& dim,
   const boost::shared_ptr<tbox::Database>& database):
   d_dim(dim),
   d_mpi(tbox::SAMRAI_MPI::commNull),
   d_max_box_size(hier::IntVector(d_dim, tbox::MathUtilities<int>::getMax())),
   d_max_lap_cut_from_center(1.0),
   d_laplace_cut_threshold_ar(0.0),
   d_log_node_history(false),
   d_log_cluster_summary(false),
   d_log_cluster(false),
   d_owner_mode("MOST_OVERLAP"),
   d_algo_advance_mode("ADVANCE_SOME"),
   d_sort_output_nodes(false),
   d_check_min_box_size('w'),
   d_barrier_before(false),
   d_barrier_after(false)
{

   /*
    * Set database-dependent parameters or cache them for use
    * when we construct a dendogram root.
    */
   if (database) {
      if (database->isInteger("max_box_size")) {
         database->getIntegerArray("max_box_size", &d_max_box_size[0], d_dim.getValue());
      }
      d_max_lap_cut_from_center =
         database->getDoubleWithDefault("max_lap_cut_from_center",
            d_max_lap_cut_from_center);
      d_laplace_cut_threshold_ar =
         database->getDoubleWithDefault("laplace_cut_threshold_ar",
            d_laplace_cut_threshold_ar);
      d_log_node_history =
         database->getBoolWithDefault("log_node_history",
            d_log_node_history);
      d_log_cluster_summary =
         database->getBoolWithDefault("log_cluster_summary",
            d_log_cluster_summary);
      d_log_cluster =
         database->getBoolWithDefault("log_cluster",
            d_log_cluster);
      d_algo_advance_mode =
         database->getStringWithDefault("algo_advance_mode",
            d_algo_advance_mode);
      d_owner_mode =
         database->getStringWithDefault("owner_mode",
            d_owner_mode);
      d_sort_output_nodes =
         database->getBoolWithDefault("sort_output_nodes",
            d_sort_output_nodes);

      std::string tmp_str;

      tmp_str = database->getStringWithDefault("check_min_box_size",
            std::string("WARN"));
      d_check_min_box_size = char(tolower(*tmp_str.c_str()));
      if (d_check_min_box_size != 'i' &&
          d_check_min_box_size != 'w' &&
          d_check_min_box_size != 'e') {
         TBOX_ERROR("BergerRigoutsos: input parameter check_min_box_size\n"
            << "can only be \"IGNORE\", \"WARN\" or \"ERROR\"");
      }

      d_barrier_before =
         database->getBoolWithDefault("barrier_before",
            d_barrier_before);
      d_barrier_after =
         database->getBoolWithDefault("barrier_after",
            d_barrier_after);
   }

}

BergerRigoutsos::~BergerRigoutsos()
{
   if (d_mpi.getCommunicator() != tbox::SAMRAI_MPI::commNull) {
      // Free the private communicator (if SAMRAI_MPI has not been finalized).
      int flag;
      tbox::SAMRAI_MPI::Finalized(&flag);
      if (!flag) {
         d_mpi.freeCommunicator();
      }
   }
}

/*
 *************************************************************************
 * Set the MPI communicator.
 *************************************************************************
 */
void
BergerRigoutsos::setMPI(
   const tbox::SAMRAI_MPI& mpi)
{
   if (d_mpi.getCommunicator() != tbox::SAMRAI_MPI::commNull) {
      d_mpi.freeCommunicator();
   }
   if (mpi.getCommunicator() != tbox::SAMRAI_MPI::commNull) {
      d_mpi.dupCommunicator(mpi);
   }
}

/*
 ************************************************************************
 *
 * Implement the BoxGeneratorStrategy interface method using
 * the asynchronous Berger-Rigoutsos implementation.
 *
 * Create objects for using the ABR recursion tree, set options for
 * using the ABR implementation, then run it.
 *
 * The output boxes from the dendogram root is in the form of a
 * BoxLevel.  This method postprocess that data to
 * convert the output to the box list form required by the
 * box clustering strategy interface.
 *
 ************************************************************************
 */
void
BergerRigoutsos::findBoxesContainingTags(
   hier::BoxLevel& new_mapped_box_level,
   hier::Connector& tag_to_new,
   hier::Connector& new_to_tag,
   const boost::shared_ptr<hier::PatchLevel>& tag_level,
   const int tag_data_index,
   const int tag_val,
   const hier::Box& bound_box,
   const hier::IntVector& min_box,
   const double efficiency_tol,
   const double combine_tol,
   const hier::IntVector& max_gcw,
   const hier::BlockId& block_id,
   const hier::LocalId& first_local_id) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS5(new_mapped_box_level,
      *tag_level,
      bound_box,
      min_box,
      max_gcw);

   tbox::SAMRAI_MPI mpi(tag_level->getBoxLevel()->getMPI());

   if (!(bound_box.numberCells() >= min_box)) {
      if (d_check_min_box_size == 'e') {
         TBOX_ERROR("BergerRigoutsos::findBoxesContainingTags input error:\n"
            << "Input box " << bound_box << " has size " << bound_box.numberCells()
            << "\nwhich is already smaller than the minimum box size\n"
            << min_box << "\n\n"
            << "To ignore or just issue a warning, see the input parameter\n"
            << "check_min_box_size.\n");
      } else if (d_check_min_box_size == 'w') {
         TBOX_WARNING("BergerRigoutsos::findBoxesContainingTags input warning:\n"
            << "Input box " << bound_box << " has size " << bound_box.numberCells()
            << "\nwhich is already smaller than the minimum box size\n"
            << min_box << "\n\n"
            << "To ignore or issue error, see the input parameter\n"
            << "check_min_box_size.\n");
      }
   }

   if (d_barrier_before) {
      t_barrier_before->start();
      mpi.Barrier();
      t_barrier_before->stop();
   }

   if (bound_box.empty()) {
      TBOX_ERROR("BergerRigoutsos: empty bounding box not allowed.");
   }

   const hier::BoxLevel& tag_mapped_box_level =
      *tag_level->getBoxLevel();

   /*
    * If using a duplicate MPI communicator, check that the duplicate
    * and the communicator of the tag_mapped_box_level are congruent.
    */
   if (d_mpi.getCommunicator() != tbox::SAMRAI_MPI::commNull) {
      tbox::SAMRAI_MPI mpi1(d_mpi);
      tbox::SAMRAI_MPI mpi2(tag_mapped_box_level.getMPI());
      TBOX_ASSERT(mpi1.getSize() == mpi2.getSize());
      TBOX_ASSERT(mpi1.getRank() == mpi2.getRank());
      if (mpi1.getSize() > 1) {
         int compare_result;
         tbox::SAMRAI_MPI::Comm_compare(
            d_mpi.getCommunicator(),
            tag_mapped_box_level.getMPI().getCommunicator(),
            &compare_result);
         if (compare_result != MPI_CONGRUENT) {
            TBOX_ERROR("BergerRigoutsos set-up error: MPI communicator\n"
               << "set by setMPI() (" << d_mpi.getCommunicator()
               << ") and the communicator of the input tag_mapped_box_level ("
               << tag_mapped_box_level.getMPI().getCommunicator() << ") are not congruent.");
         }
      }
   }

   t_find_boxes_with_tags->start();

   BergerRigoutsosNode root_node(d_dim, block_id, first_local_id);

   // Set standard Berger-Rigoutsos clustering parameters.
   root_node.setClusteringParameters(tag_data_index,
      tag_val,
      min_box,
      efficiency_tol,
      combine_tol,
      d_max_box_size,
                                     d_max_lap_cut_from_center,
                                     d_laplace_cut_threshold_ar);

   // Set the parallel algorithm and DLBG parameters.
   root_node.setAlgorithmAdvanceMode(d_algo_advance_mode);
   root_node.setOwnerMode(d_owner_mode);
   root_node.setComputeRelationships("BIDIRECTIONAL", max_gcw);

   // Set debugging/verbosity parameters.
   root_node.setLogNodeHistory(d_log_node_history);

   root_node.clusterAndComputeRelationships(new_mapped_box_level,
      tag_to_new,
      new_to_tag,
      bound_box,
      tag_level,
      d_mpi);

   if (d_sort_output_nodes == true) {
      /*
       * Sorting the node indices is not required.
       * This optional step makes the results order
       * deterministic, which makes the results repeatable.
       * (The natural order of the output of the asynchronous
       * clustering algorithm is non-deterministic because
       * it depends on the order of asynchronous messages.)
       */
      sortOutputBoxes(new_mapped_box_level,
         tag_to_new,
         new_to_tag);
   }

   /*
    * Get some global parameters.  Do it before logging to prevent
    * the logging flag from having an undue side effect on performance.
    */
   t_global_reductions->start();
   new_mapped_box_level.getGlobalNumberOfBoxes();
   new_mapped_box_level.getGlobalNumberOfCells();
   new_mapped_box_level.getGlobalBoundingBox(block_id.getBlockValue());
   t_global_reductions->stop();

   if (d_log_cluster) {
      tbox::plog << "BergerRigoutsos cluster log:\n"
      << "\tNew mapped_box_level clustered by BergerRigoutsos:\n" << new_mapped_box_level.format("",
         2)
      << "\tBergerRigoutsos tag_to_new:\n" << tag_to_new.format("", 2)
      << "\tBergerRigoutsos new_to_tag:\n" << new_to_tag.format("", 2);
   }
   if (d_log_cluster_summary) {
      /*
       * Log summary of clustering and dendogram.
       */
      tbox::plog << "BergerRigoutsos summary:\n"
                 << "\tAsync BR on proc " << mpi.getRank()
                 << " owned "
                 << root_node.getMaxOwnership() << " participating in "
                 << root_node.getMaxNodes() << " nodes ("
                 << (double)root_node.getMaxOwnership() / root_node.getMaxNodes()
                 << ") in " << root_node.getMaxGeneration() << " generations,"
                 << "   " << root_node.getNumBoxesGenerated()
                 << " boxes generated.\n\t"
                 << root_node.getMaxTagsOwned() << " locally owned tags on new BoxLevel.\n\t"
                 << "Initial bounding box = " << bound_box << ", "
                 << bound_box.size() << " cells, "
                 << "final global bounding box = "
                 << new_mapped_box_level.getGlobalBoundingBox(block_id.getBlockValue())
                 << ", "
                 << new_mapped_box_level.getGlobalBoundingBox(block_id.getBlockValue()).size()
                 << " cells\n\t"
                 << "Final output has " << root_node.getNumTags()
                 << " tags in "
                 << new_mapped_box_level.getGlobalNumberOfCells()
                 << " global cells [" << new_mapped_box_level.getMinNumberOfCells()
                 << "-" << new_mapped_box_level.getMaxNumberOfCells() << "], "
                 << new_mapped_box_level.getGlobalNumberOfBoxes()
                 << " global mapped boxes [" << new_mapped_box_level.getMinNumberOfBoxes()
                 << "-" << new_mapped_box_level.getMaxNumberOfBoxes() << "]\n\t"
                 << "Number of continuations: avg = "
                 << root_node.getAvgNumberOfCont()
                 << "   max = " << root_node.getMaxNumberOfCont() << '\n'
                 << "\tBergerRigoutsos new_level summary:\n" << new_mapped_box_level.format("\t\t",0)
                 << "\tBergerRigoutsos new_level statistics:\n" << new_mapped_box_level.formatStatistics("\t\t")
                 << "\tBergerRigoutsos new_to_tag summary:\n" << new_to_tag.format("\t\t",0)
                 << "\tBergerRigoutsos new_to_tag statistics:\n" << new_to_tag.formatStatistics("\t\t")
                 << "\tBergerRigoutsos tag_to_new summary:\n" << tag_to_new.format("\t\t",0)
                 << "\tBergerRigoutsos tag_to_new statistics:\n" << tag_to_new.formatStatistics("\t\t")
                 << "\n";
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   assertNoMessageForPrivateCommunicator();
#endif

   if (d_barrier_after) {
      t_barrier_after->start();
      mpi.Barrier();
      t_barrier_after->stop();
   }

   t_find_boxes_with_tags->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BergerRigoutsos::sortOutputBoxes(
   hier::BoxLevel& new_mapped_box_level,
   hier::Connector& tag_to_new,
   hier::Connector& new_to_tag) const
{

   t_sort_output_nodes->start();

   hier::OverlapConnectorAlgorithm oca;

   if (0) {
      // Check inputs.
      int errs = 0;
      if (oca.checkOverlapCorrectness(tag_to_new, false, true)) {
         ++errs;
         tbox::perr << "Error found in tag_to_new!\n";
      }
      if (oca.checkOverlapCorrectness(new_to_tag, false, true)) {
         ++errs;
         tbox::perr << "Error found in new_to_tag!\n";
      }
      if (errs != 0) {
         TBOX_ERROR(
            "Errors found before sorting nodes."
            << "new_mapped_box_level:\n" << new_mapped_box_level.format("", 2)
            << "tag mapped_box_level:\n" << tag_to_new.getBase().format("", 2)
            << "tag_to_new:\n" << tag_to_new.format("", 2)
            << "new_to_tag:\n" << new_to_tag.format("", 2));
      }
   }

   /*
    * Sort local indices by corners to make the output deterministic.
    */
   hier::Connector sorting_map;
   hier::BoxLevel sorted_mapped_box_level(d_dim);
   hier::BoxLevelConnectorUtils dlbg_edge_utils;
   dlbg_edge_utils.makeSortingMap(
      sorted_mapped_box_level,
      sorting_map,
      new_mapped_box_level,
      true /* sort nodes by corners */,
      false /* don't sequentialize indices globally */);
   if (0) {
      tbox::plog
      << "tag mapped_box_level:\n" << tag_to_new.getBase().format("", 2)
      << "tag_to_new:\n" << tag_to_new.format("", 2)
      << "new_to_tag:\n" << new_to_tag.format("", 2)
      << "Sorting map:\n" << sorting_map.format("", 2);
   }
   if (0) {
      // Check sorting_map before using it.
      int errs = 0;
      if (oca.checkOverlapCorrectness(sorting_map, false, true)) {
         ++errs;
         tbox::perr << "Error found in sorting_map!\n";
      }
      if (errs != 0) {
         TBOX_ERROR(
            "Errors in load balance mapping found."
            << "presorted mapped_box_level:\n" << new_mapped_box_level.format("", 2)
            << "sorted mapped_box_level:\n" << sorted_mapped_box_level.format("", 2)
            << "sorting_map:\n" << sorting_map.format("", 2));
      }
   }
   hier::MappingConnectorAlgorithm mca;
   mca.modify(tag_to_new,
      new_to_tag,
      sorting_map,
      &new_mapped_box_level);
   if (0) {
      // Check result of mapping.
      int errs = 0;
      if (oca.checkOverlapCorrectness(tag_to_new, false, true)) {
         ++errs;
         tbox::perr << "Error found in tag_to_new!\n";
      }
      if (oca.checkOverlapCorrectness(new_to_tag, false, true)) {
         ++errs;
         tbox::perr << "Error found in new_to_tag!\n";
      }
      if (errs != 0) {
         TBOX_ERROR(
            "Errors found after sorting nodes."
            << "new_mapped_box_level:\n" << new_mapped_box_level.format("", 2)
            << "tag mapped_box_level:\n" << tag_to_new.getBase().format("", 2)
            << "tag_to_new:\n" << tag_to_new.format("", 2)
            << "new_to_tag:\n" << new_to_tag.format("", 2));
      }
   }

   t_sort_output_nodes->stop();
}

/*
 ***************************************************************************
 *
 ***************************************************************************
 */
void
BergerRigoutsos::assertNoMessageForPrivateCommunicator() const
{
   /*
    * If using a private communicator, double check to make sure
    * there are no remaining messages.  This is not a guarantee
    * that there is no messages in transit, but it can find
    * messages that have arrived but not received.
    */
   if (d_mpi.getCommunicator() != tbox::SAMRAI_MPI::commNull) {
      int flag;
      tbox::SAMRAI_MPI::Status mpi_status;
      int mpi_err = d_mpi.Iprobe(MPI_ANY_SOURCE,
            MPI_ANY_TAG,
            &flag,
            &mpi_status);
      if (mpi_err != MPI_SUCCESS) {
         TBOX_ERROR("Error probing for possible lost messages.");
      }
      if (flag == true) {
         int count = -1;
         mpi_err = tbox::SAMRAI_MPI::Get_count(&mpi_status, MPI_INT, &count);
         TBOX_ERROR("Library error!\n"
            << "BergerRigoutsos detected before or after\n"
            << "running BergerRigoutsosNode that there\n"
            << "is a message yet to be received.  This is\n"
            << "an error because all messages using the\n"
            << "private communicator should have been\n"
            << "accounted for.  Message status:\n"
            << "source " << mpi_status.MPI_SOURCE << '\n'
            << "tag " << mpi_status.MPI_TAG << '\n'
            << "count " << count << " (assuming integers)\n");
      }
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
BergerRigoutsos::initializeCallback()
{
   TBOX_ASSERT(!t_global_reductions);
   t_run_abr = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsos::run_abr");
   t_find_boxes_with_tags = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsos::find_boxes_with_tags");
   t_global_reductions = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsos::global_reductions");
   t_sort_output_nodes = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsos::sort_output_nodes");
   t_barrier_before = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsos::barrier_before");
   t_barrier_after = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsos::barrier_after");
}

/*
 ***************************************************************************
 *
 * Release static timers.  To be called by shutdown registry to make sure
 * memory for timers does not leak.
 *
 ***************************************************************************
 */
void
BergerRigoutsos::finalizeCallback()
{
   t_barrier_before.reset();
   t_barrier_after.reset();
   t_find_boxes_with_tags.reset();
   t_run_abr.reset();
   t_global_reductions.reset();
   t_sort_output_nodes.reset();
}

}
}
#endif
