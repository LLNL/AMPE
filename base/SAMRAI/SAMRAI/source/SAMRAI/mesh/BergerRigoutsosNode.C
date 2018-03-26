/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Node in asynchronous Berger-Rigoutsos dendogram
 *
 ************************************************************************/
#ifndef included_mesh_BergerRigoutsosNode_C
#define included_mesh_BergerRigoutsosNode_C

#include <cstring>
#include <set>
#include <algorithm>

#include "SAMRAI/mesh/BergerRigoutsosNode.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

const int BergerRigoutsosNode::BAD_INTEGER = -9999999;

tbox::StartupShutdownManager::Handler
BergerRigoutsosNode::s_initialize_handler(
   BergerRigoutsosNode::initializeCallback,
   0,
   0,
   BergerRigoutsosNode::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_cluster;
boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_cluster_and_compute_relationships;
boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_continue_algorithm;
boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_compute;
boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_comm_wait;
boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_MPI_wait;
boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_compute_new_graph_relationships;
boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_share_new_relationships;
boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_share_new_relationships_send;
boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_share_new_relationships_recv;
boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_share_new_relationships_unpack;
boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_local_tasks;
boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_local_histogram;
boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_reduce_histogram;
boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_bcast_acceptability;
boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_gather_grouping_criteria;
boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_bcast_child_groups;
boost::shared_ptr<tbox::Timer> BergerRigoutsosNode::CommonParams::t_bcast_to_dropouts;

/*
 *******************************************************************
 * Construct root node of the dendogram.
 *******************************************************************
 */
BergerRigoutsosNode::BergerRigoutsosNode(
   const tbox::Dimension& dim,
   const hier::BlockId& block_id,
   const hier::LocalId& first_local_id):
   d_dim(dim),
   d_pos(1),
   d_common(new CommonParams(d_dim)),
   d_parent(NULL),
   d_lft_child(NULL),
   d_rht_child(NULL),
   d_box(d_dim),
   d_owner(0),
   d_group(0),
   d_mpi_tag(-1),
   d_overlap(-1),
   d_box_acceptance(undetermined),
   d_mapped_box(d_dim),
   d_mapped_box_iterator(hier::BoxContainer().end()),
   d_wait_phase(to_be_launched),
   d_send_msg(),
   d_recv_msg(),
   d_comm_group(NULL),
   d_block_id(block_id),
   d_first_local_id(first_local_id),
   d_generation(1),
   d_n_cont(0)
{

   ++(d_common->num_nodes_owned);
   ++(d_common->max_nodes_owned);
   ++(d_common->num_nodes_allocated);
   ++(d_common->max_nodes_allocated);
   if (d_common->max_generation < d_generation) {
      d_common->max_generation = d_generation;
   }

   if (d_common->log_node_history) {
      tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                 << d_common->num_nodes_active << "-active  "
                 << d_common->num_nodes_owned << "-owned  "
                 << d_common->num_nodes_completed << "-completed  "
                 << "Construct " << d_generation << ':' << d_pos
                 << ".\n";
   }
}

/*
 *******************************************************************
 * Construct non-root node of the dendogram.  This is private!
 * Public constructors are only for root nodes.
 *******************************************************************
 */
BergerRigoutsosNode::BergerRigoutsosNode(
   CommonParams* common_params,
   BergerRigoutsosNode* parent,
   const int child_number,
   const hier::BlockId& block_id,
   const hier::LocalId& first_local_id):
   d_dim(parent->d_dim),
   d_pos((parent->d_pos > 0 && parent->d_pos <
          tbox::MathUtilities<int>::getMax() / 2) ?
         2 * parent->d_pos + child_number :
         (child_number == 0 ? -1 : -2)),
   d_common(common_params),
   d_parent(parent),
   d_lft_child(NULL),
   d_rht_child(NULL),
   d_box(d_dim),
   d_owner(-1),
   d_group(0),
   d_mpi_tag(-1),
   d_overlap(-1),
   d_box_acceptance(undetermined),
   d_mapped_box(d_dim),
   d_mapped_box_iterator(hier::BoxContainer().end()),
   d_wait_phase(for_data_only),
   d_send_msg(),
   d_recv_msg(),
   d_comm_group(NULL),
   d_block_id(block_id),
   d_first_local_id(first_local_id),
   d_generation(d_parent->d_generation + 1),
   d_n_cont(0)
{

#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parent->d_pos >= 0 && d_pos < 0) {
      TBOX_WARNING("Too many generations for node identification.\n"
         << "The node id cannot be increased any further.\n"
         << "This affects only the node id, which is only\n"
         << "used for analysis and debugging and does not\n"
         << "affect the algorithm.\n"
         << "Last valid node id is " << d_parent->d_pos << '\n');
   }
#endif

#ifdef DEBUG_CHECK_ASSERTIONS
   d_mapped_box_iterator = BoxContainer().end();
#endif

   ++(d_common->num_nodes_allocated);
   ++(d_common->max_nodes_allocated);
   if (d_common->max_generation < d_generation) {
      d_common->max_generation = d_generation;
   }

   if (d_common->log_node_history) {
      tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                 << d_common->num_nodes_active << "-active  "
                 << d_common->num_nodes_owned << "-owned  "
                 << d_common->num_nodes_completed << "-completed  "
                 << "Construct " << d_generation << ':' << d_pos
                 << ", child of "
                 << d_parent->d_generation << ':' << d_parent->d_pos
                 << "   " << d_parent->d_mapped_box
                 << ".\n";
   }
}

BergerRigoutsosNode::~BergerRigoutsosNode()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * Forbid deleting a node that is running because there may
    * be pending communication (by the node or its children).
    * Note that this is NOT an extra restriction over the
    * recursive implementation.
    */
   if (d_wait_phase != for_data_only &&
       d_wait_phase != to_be_launched &&
       d_wait_phase != completed) {
      TBOX_ERROR("Should not delete a node that is currently running\n"
         << "the Berger-Rigoutsos algorithm because there\n"
         << "may be pending communications.");
   }
#endif

   if (d_comm_group != NULL) {
      if (!d_comm_group->isDone()) {
         TBOX_ERROR("Library error: Destructing a node with an unfinished\n"
            << "communication tree is bad because it leaves\n"
            << "pending MPI messages.");
      }
      delete d_comm_group;
      d_comm_group = NULL;
   }

   --(d_common->num_nodes_allocated);

   if (d_parent != NULL && d_common->log_node_history) {
      tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                 << d_common->num_nodes_active << "-active  "
                 << d_common->num_nodes_owned << "-owned  "
                 << d_common->num_nodes_completed << "-completed  "
                 << "Destruct " << d_generation << ':' << d_pos
                 << "  " << d_mapped_box
                 << "  " << d_box
                 << ".\n";
   }

   if (d_parent == NULL) {
      delete d_common;
   }

   d_wait_phase = deallocated;
}

/*
 ********************************************************************
 ********************************************************************
 */
void
BergerRigoutsosNode::setClusteringParameters(
   const int tag_data_index,
   const int tag_val,
   const hier::IntVector min_box,
   const double efficiency_tol,
   const double combine_tol,
   const hier::IntVector& max_box_size,
   const double max_lap_cut_from_center,
   const double laplace_cut_threshold_ar)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(d_dim, min_box, max_box_size);

   d_common->tag_data_index = tag_data_index;
   d_common->tag_val = tag_val;
   d_common->min_box = min_box;
   d_common->efficiency_tol = efficiency_tol;
   d_common->combine_tol = combine_tol;
   d_common->max_box_size = max_box_size;
   d_common->max_lap_cut_from_center = max_lap_cut_from_center;
   d_common->laplace_cut_threshold_ar = laplace_cut_threshold_ar;
}

/*
 ********************************************************************
 ********************************************************************
 */
void
BergerRigoutsosNode::clusterAndComputeRelationships(
   hier::BoxLevel& new_mapped_box_level,
   hier::Connector& tag_to_new,
   hier::Connector& new_to_tag,
   const hier::Box& bound_box,
   const boost::shared_ptr<hier::PatchLevel>& tag_level,
   const tbox::SAMRAI_MPI& mpi_object)
{
   TBOX_ASSERT(d_parent == NULL);
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS3(d_dim,
      new_mapped_box_level,
      bound_box,
      *tag_level);

   d_common->tag_level = tag_level;
   d_common->tag_mapped_box_level = tag_level->getBoxLevel().get();

   setMPI(mpi_object);

   d_box = bound_box;

   /*
    * During the algorithm, we kept the results in primitive
    * containers to avoid the overhead of fine-grain changes to the
    * output objects.  Now initialize the outputs using those
    * primitive containers.
    */

   new_mapped_box_level.initialize(
      d_common->tag_level->getRatioToLevelZero(),
      d_common->tag_level->getGridGeometry(),
      d_common->tag_mapped_box_level->getMPI(),
      hier::BoxLevel::DISTRIBUTED);

   if (d_common->compute_relationships >= 1) {
      tag_to_new.clearNeighborhoods();
      tag_to_new.setBase(*tag_level->getBoxLevel());
      tag_to_new.setHead(new_mapped_box_level);
      tag_to_new.setWidth(d_common->max_gcw, true);
   }
   if (d_common->compute_relationships >= 2) {
      new_to_tag.clearNeighborhoods();
      new_to_tag.setBase(new_mapped_box_level);
      new_to_tag.setHead(*tag_level->getBoxLevel());
      new_to_tag.setWidth(d_common->max_gcw, true);
   }

   d_common->new_mapped_box_level = &new_mapped_box_level;
   d_common->tag_to_new = &tag_to_new;
   d_common->new_to_tag = &new_to_tag;

   clusterAndComputeRelationships();

   new_mapped_box_level.finalize();

   /*
    * Clear temporary parameters that are only used during active
    * clustering.
    */
   d_common->new_mapped_box_level = NULL;
   d_common->tag_to_new = NULL;
   d_common->new_to_tag = NULL;
   d_common->tag_mapped_box_level = NULL;
   d_common->tag_level.reset();
}

/*
 ********************************************************************
 * Preprocess, run the asynchronous algorithm and postprocess.
 ********************************************************************
 */
void
BergerRigoutsosNode::clusterAndComputeRelationships()
{
   TBOX_ASSERT(d_parent == NULL);
   tbox::SAMRAI_MPI mpi(d_common->mpi_object);

   d_common->t_cluster_and_compute_relationships->start();
   d_common->t_cluster->start();

   /*
    * If compute_relationships == 1:
    *   - Compute relationships from tagged level to new levels.
    *     These relationships are organized around the tagged nodes.
    *     They do not need to be shared with the owners of the
    *     new nodes.
    *
    * If compute_relationships == 2:
    *   - Compute relationships as in compute_relationships == 1 case.
    *   - Owners of new relationships send new relationship data to owners
    *     of new nodes.  This creates the neighbor data
    *     organized around the new nodes.
    */

   if (d_common->compute_relationships > 0) {

      /*
       * Create empty neighbor lists for nodes on tagged mapped_box_level.
       * As new nodes are finalized, they will be added to
       * these lists.
       */
      const BoxContainer& tag_mapped_boxes =
         d_common->tag_mapped_box_level->getBoxes();
      for (hier::RealBoxConstIterator ni(tag_mapped_boxes.realBegin());
           ni != tag_mapped_boxes.realEnd(); ++ni) {
         d_common->tag_to_new->makeEmptyLocalNeighborhood(ni->getId());
      }
      TBOX_ASSERT(
         static_cast<int>(d_common->tag_mapped_box_level->getLocalNumberOfBoxes()) ==
         d_common->tag_to_new->getLocalNumberOfNeighborSets());

   }

   TBOX_ASSERT(d_common->algo_advance_mode == ADVANCE_SOME ||
      d_common->algo_advance_mode == ADVANCE_ANY ||
      d_common->algo_advance_mode == SYNCHRONOUS);            // No other supported currently.
   {
      int n_comm_group_completed = 0;
      // d_common->relaunch_queue.appendItem(this);
      d_common->relaunch_queue.push_back(this);


      do {

         d_common->t_compute->start();
         while (!d_common->relaunch_queue.empty()) {
            BergerRigoutsosNode* node_for_relaunch = d_common->relaunch_queue.front();
            d_common->relaunch_queue.pop_front();
            if (0) {
               tbox::plog << "Continuing from queue ";
               node_for_relaunch->printState(tbox::plog);
               tbox::plog << std::endl;
            }
            node_for_relaunch->continueAlgorithm();
            if (0) {
               tbox::plog << "Exiting continueAlgorithm ";
               node_for_relaunch->printState(tbox::plog);
               tbox::plog << std::endl;
            }
         }
         d_common->t_compute->stop();

         d_common->t_comm_wait->start();
         n_comm_group_completed =
            static_cast<int>(d_common->comm_stage.advanceSome());
         d_common->t_comm_wait->stop();

         d_common->t_compute->start();
         while ( d_common->comm_stage.numberOfCompletedMembers() > 0 ) {
            BergerRigoutsosNode* node_for_relaunch =
               (BergerRigoutsosNode *)(d_common->comm_stage.popCompletionQueue()->getHandler());
            if (0) {
               tbox::plog << "Continuing from stage ";
               node_for_relaunch->printState(tbox::plog);
               tbox::plog << std::endl;
            }
            node_for_relaunch->continueAlgorithm();
            if (0) {
               tbox::plog << "Exiting continueAlgorithm ";
               node_for_relaunch->printState(tbox::plog);
               tbox::plog << std::endl;
            }
         }
         d_common->t_compute->stop();

         if (0) {
            tbox::plog << "relaunch_queue size "
                       << d_common->relaunch_queue.size()
                       << "   groups completed: " << n_comm_group_completed
                       << std::endl;
            tbox::plog << "Stage has " << d_common->comm_stage.numberOfMembers()
                       << " members, "
                       << d_common->comm_stage.numberOfPendingMembers()
                       << " pending members, "
                       << d_common->comm_stage.numberOfPendingRequests()
                       << " pending requests." << std::endl;
         }
      } while ( !d_common->relaunch_queue.empty() || d_common->comm_stage.hasPendingRequests() );

   }

#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_wait_phase != completed) {
      printDendogramState(tbox::plog, "ERR->");
      TBOX_ERROR(
         "Root node finished but d_wait_phase is not set to completed.\n"
         << "d_wait_phase=" << d_wait_phase);
   }
   if (d_common->compute_relationships > 2) {
      // Each new node should have its own neighbor list.
      TBOX_ASSERT(d_common->new_mapped_box_level->getBoxes().size() ==
         d_common->new_to_tag->getLocalNumberOfNeighborSets());
   }
#endif

   // Barrier to separate clustering cost from relationship sharing cost.
   mpi.Barrier();

   d_common->t_cluster->stop();

   /*
    * Share relationships with owners, if requested.
    * This is a one-time operation that is not considered a part
    * of continueAlgorithm(), so it lies outside that timimg.
    */
   if (d_common->compute_relationships > 1) {
      shareNewNeighborhoodSetsWithOwners();
   }

   d_common->t_cluster_and_compute_relationships->stop();
}

/*
 **************************************************************************
 * This private method sets multiple internal data and presumes that
 * d_common->tag_level has been set.
 **************************************************************************
 */
void
BergerRigoutsosNode::setMPI(
   const tbox::SAMRAI_MPI& mpi_object)
{
   TBOX_ASSERT(d_parent == NULL);
   TBOX_ASSERT(d_common->tag_mapped_box_level != NULL);

   /*
    * If a valid MPI communicator is given, use it instead of the
    * tag_mapped_box_level's communicator.  It must be congruent with
    * the tag_mapped_box_level's.
    */
   if (mpi_object.getCommunicator() != tbox::SAMRAI_MPI::commNull) {
#if defined(DEBUG_CHECK_ASSERTIONS)
      /*
       * If user supply a communicator to use, make sure it is
       * compatible with the BoxLevel involved.
       */
      tbox::SAMRAI_MPI mpi1(mpi_object);
      tbox::SAMRAI_MPI mpi2(d_common->tag_mapped_box_level->getMPI());
      TBOX_ASSERT(mpi1.getSize() == mpi2.getSize());
      TBOX_ASSERT(mpi1.getRank() == mpi2.getRank());
      if (mpi1.getSize() > 1) {
         int compare_result;
         tbox::SAMRAI_MPI::Comm_compare(
            mpi_object.getCommunicator(),
            d_common->tag_mapped_box_level->getMPI().getCommunicator(),
            &compare_result);
         if (compare_result != MPI_CONGRUENT) {
            TBOX_ERROR("BergerRigoutsosNode::setMPI:\n"
               << "MPI communicator (" << mpi_object.getCommunicator()
               << ") and the communicator of the input tag_mapped_box_level ("
               << d_common->tag_mapped_box_level->getMPI().getCommunicator()
               << ") are not congruent.");
         }
      }
#endif
      d_common->mpi_object = mpi_object;
   } else {
      d_common->mpi_object =
         d_common->tag_mapped_box_level->getMPI();
   }

   tbox::SAMRAI_MPI mpi(d_common->mpi_object);

   /*
    * Reserve the tag upper bound for the relationship-sharing phase.
    * Divide the rest into tag pools divided among all processes.
    */
   if (tbox::SAMRAI_MPI::usingMPI()) {
      /*
       * For some MPI implementations, I cannot get the attribute for
       * any communicator except for MPI_COMM_WORLD.  Assuming the tag
       * upper bound is the same for all communicators, I will try
       * some other communicators to get it.
       */
      int* tag_upper_bound_ptr, flag;
      mpi.Attr_get(MPI_TAG_UB,
         &tag_upper_bound_ptr,
         &flag);
      if (tag_upper_bound_ptr == NULL) {
         tbox::SAMRAI_MPI mpi1(tbox::SAMRAI_MPI::getSAMRAIWorld().getCommunicator());
         mpi1.Attr_get(MPI_TAG_UB,
            &tag_upper_bound_ptr,
            &flag);
      }
      if (tag_upper_bound_ptr == NULL) {
         tbox::SAMRAI_MPI mpi1(tbox::SAMRAI_MPI::commWorld);
         mpi1.Attr_get(MPI_TAG_UB,
            &tag_upper_bound_ptr,
            &flag);
      }
      TBOX_ASSERT(tag_upper_bound_ptr != NULL);
      d_common->tag_upper_bound = *tag_upper_bound_ptr;

   } else {
      d_common->tag_upper_bound = 1000000;
   }

   d_common->nproc = mpi.getSize();
   d_common->rank = mpi.getRank();

   d_common->available_mpi_tag =
      d_common->tag_upper_bound / d_common->nproc * d_common->rank;

   // Divide the rest into tag pools divided among all processes.
   d_common->available_mpi_tag = d_common->tag_upper_bound / d_common->nproc
      * d_common->rank;

   if (d_common->rank == 0) {
      claimMPITag();
   } else {
      d_mpi_tag = 0;
   }
   /*
    * Even though owner has different way of getting MPI tag,
    * all processes should start with the same mpi tag.
    */
   TBOX_ASSERT(d_mpi_tag == 0);

   /*
    * Set the processor group (for the root dendogram node).
    */
   d_group.resize(d_common->nproc, BAD_INTEGER);
   for (unsigned int i = 0; i < d_group.size(); ++i) {
      d_group[i] = i;
   }
}

/*
 ********************************************************************
 * This method looks messy, but it is just the BR agorithm,
 * with multiple pause and continue points implemented by
 * the goto and labels.  Each pause point is accompanied by
 * a line setting d_wait_phase so that the algorithm can
 * continue where it left off when this method is called again.
 * The BR algorithm is not completed until this method returns
 * the WaitPhase value "completed".
 ********************************************************************
 */
BergerRigoutsosNode::WaitPhase
BergerRigoutsosNode::continueAlgorithm()
{
   d_common->t_continue_algorithm->start();
   ++d_n_cont;

   TBOX_ASSERT(d_parent == NULL || d_parent->d_wait_phase != completed);
   // TBOX_ASSERT( ! inRelaunchQueue(this) );
   TBOX_ASSERT(inRelaunchQueue(this) == d_common->relaunch_queue.end());

   /*
    * Skip right to where we left off,
    * which is specified by the wait phase variable.
    */
   switch (d_wait_phase) {
      case for_data_only:
         TBOX_ERROR("Library error: Attempt to execute data-only node.");
      case to_be_launched:
         goto TO_BE_LAUNCHED;
      case reduce_histogram:
         goto REDUCE_HISTOGRAM;
      case bcast_acceptability:
         goto BCAST_ACCEPTABILITY;
      case gather_grouping_criteria:
         goto GATHER_GROUPING_CRITERIA;
      case bcast_child_groups:
         goto BCAST_CHILD_GROUPS;
      case run_children:
         goto RUN_CHILDREN;
      case bcast_to_dropouts:
         goto BCAST_TO_DROPOUTS;
      case completed:
         TBOX_ERROR("Library error: Senseless continuation of completed node.");
      default:
         TBOX_ERROR("Library error: Nonexistent phase.");
   }

   bool sub_completed;

   /*
    * Delegated tasks: Major tasks are delegated to private methods.
    * These methods may check whether the process is the owner or
    * just a contributor and fork appropriately.  The communication
    * checking tasks return whether communication is completed, but
    * they do NOT change the d_wait_phase variable, which is done
    * in this function.
    */

TO_BE_LAUNCHED:

   ++(d_common->num_nodes_active);
   if (d_common->max_nodes_active < d_common->num_nodes_active) {
      d_common->max_nodes_active = d_common->num_nodes_active;
   }

   if (d_common->log_node_history) {
      tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                 << d_common->num_nodes_active << "-active  "
                 << d_common->num_nodes_owned << "-owned  "
                 << d_common->num_nodes_completed << "-completed  "
                 << "Commence " << d_generation << ':' << d_pos
                 << "  " << d_box
                 << "  accept=" << d_box_acceptance
                 << "  ovlap=" << d_overlap
                 << "  owner=" << d_owner
                 << "  gsize=" << d_group.size()
                 << ".\n";
   }

   if (d_parent == NULL || d_overlap > 0 || d_common->rank == d_owner) {

      TBOX_ASSERT(inGroup(d_group));

      // Set up communication group for operations in participating group.
      d_comm_group = new tbox::AsyncCommGroup(
            computeCommunicationTreeDegree(static_cast<int>(d_group.size())),
            &d_common->comm_stage,
            this);
      d_comm_group->setGroupAndRootRank(d_common->mpi_object,
         &d_group[0], static_cast<int>(d_group.size()), d_owner);
      if (d_parent == NULL) {
         /*
          * For the global group, MPI collective functions are presumably
          * faster than the peer-to-peer collective implementation in
          * AsyncCommGroup.
          *
          * Enable this mode only for the root node.  Child nodes are
          * not guaranteed to execute the communication operation at
          * the same point on all processors (even if all proccessors
          * participate).
          */
         d_comm_group->setUseMPICollectiveForFullGroups(true);
      }

      d_common->t_local_tasks->start();
      makeLocalTagHistogram();
      d_common->t_local_tasks->stop();

      if (d_group.size() > 1) {
         d_common->t_reduce_histogram->start();
         reduceHistogram_start();
REDUCE_HISTOGRAM:
         if (!d_common->t_reduce_histogram->isRunning())
            d_common->t_reduce_histogram->start();
         if (d_common->algo_advance_mode == SYNCHRONOUS) {
            d_comm_group->completeCurrentOperation();
         }
         sub_completed = reduceHistogram_check();
         d_common->t_reduce_histogram->stop();
         if (!sub_completed) {
            d_wait_phase = reduce_histogram;
            goto RETURN;
         }
      }

      if (d_common->rank == d_owner) {
         /*
          * The owner node saves the tag count.  Participant nodes get
          * tag count from broadcastAcceptability().  This data is just for
          * analysis (not required) and I expect it to have trivial cost.
          */
         int narrowest_dir = 0;
         for (int d = 0; d < d_dim.getValue(); ++d) {
            if (d_histogram[d].size() < d_histogram[narrowest_dir].size())
               narrowest_dir = d;
         }
         d_num_tags = 0;
         for (size_t i = 0; i < d_histogram[narrowest_dir].size(); ++i) {
            d_num_tags += d_histogram[narrowest_dir][i];
         }

         /*
          * If this is the root node, d_num_tags is the total tag count
          * in all nodes.
          */
         if (d_parent == NULL) {
            d_common->num_tags_in_all_nodes = d_num_tags;
         }
      }

      if (d_common->rank == d_owner) {
         d_common->t_local_tasks->start();
         computeMinimalBoundingBoxForTags();
         acceptOrSplitBox();
         d_common->t_local_tasks->stop();
         TBOX_ASSERT(boxAccepted() || boxRejected() ||
            (boxHasNoTag() && d_parent == 0));
         if (!boxHasNoTag()) {
            /*
             * A mapped_box_level node is created even if box is not acceptable,
             * so that the children can reference its local index in case
             * the box is later accepted based on the combined tolerance
             * of the children.  The node would be erased later if
             * it is not finally accepted.
             */
            createBox();
         }
      }

      if (d_group.size() > 1) {
         d_common->t_bcast_acceptability->start();
         broadcastAcceptability_start();
BCAST_ACCEPTABILITY:
         if (!d_common->t_bcast_acceptability->isRunning())
            d_common->t_bcast_acceptability->start();
         if (d_common->algo_advance_mode == SYNCHRONOUS) {
            d_comm_group->completeCurrentOperation();
         }
         sub_completed = broadcastAcceptability_check();
         d_common->t_bcast_acceptability->stop();
         if (!sub_completed) {
            d_wait_phase = bcast_acceptability;
            goto RETURN;
         }
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_common->rank == d_owner) {
         TBOX_ASSERT(d_box_acceptance == accepted_by_calculation ||
            d_box_acceptance == rejected_by_calculation ||
            d_box_acceptance == hasnotag_by_owner);
      } else {
         TBOX_ASSERT(d_box_acceptance == accepted_by_owner ||
            d_box_acceptance == rejected_by_owner ||
            d_box_acceptance == hasnotag_by_owner);
      }
#endif

      /*
       * If this is the root node, d_num_tags is the total tag count
       * in all nodes.
       */
      if (d_parent == NULL && d_common->rank != d_owner) {
         d_common->num_tags_in_all_nodes = d_num_tags;
      }

      if (boxRejected()) {

         /*
          * Compute children groups and owners without assuming
          * entire mesh structure is known locally.
          */
         d_common->t_local_tasks->start();
         countOverlapWithLocalPatches();
         d_common->t_local_tasks->stop();

         if (d_group.size() > 1) {
            d_common->t_gather_grouping_criteria->start();
            gatherGroupingCriteria_start();
GATHER_GROUPING_CRITERIA:
            if (!d_common->t_gather_grouping_criteria->isRunning())
               d_common->t_gather_grouping_criteria->start();
            if (d_common->algo_advance_mode == SYNCHRONOUS) {
               d_comm_group->completeCurrentOperation();
            }
            sub_completed = gatherGroupingCriteria_check();
            d_common->t_gather_grouping_criteria->stop();
            if (!sub_completed) {
               d_wait_phase = gather_grouping_criteria;
               goto RETURN;
            }
         }

         if (d_common->rank == d_owner) {
            d_common->t_local_tasks->start();
            formChildGroups();
            d_common->t_local_tasks->stop();
         }

         if (d_group.size() > 1) {
            d_common->t_bcast_child_groups->start();
            broadcastChildGroups_start();
BCAST_CHILD_GROUPS:
            if (!d_common->t_bcast_child_groups->isRunning())
               d_common->t_bcast_child_groups->start();
            if (d_common->algo_advance_mode == SYNCHRONOUS) {
               d_comm_group->completeCurrentOperation();
            }
            sub_completed = broadcastChildGroups_check();
            d_common->t_bcast_child_groups->stop();
            if (!sub_completed) {
               d_wait_phase = bcast_child_groups;
               goto RETURN;
            }
         }

         if (d_lft_child->d_owner == d_common->rank) {
            ++(d_common->num_nodes_owned);
            ++(d_common->max_nodes_owned);
         }
         if (d_rht_child->d_owner == d_common->rank) {
            ++(d_common->num_nodes_owned);
            ++(d_common->max_nodes_owned);
         }

         runChildren_start();
RUN_CHILDREN:
         sub_completed = runChildren_check();
         if (!sub_completed) {
            d_wait_phase = run_children;
            goto RETURN;
         }
      } else if (boxAccepted()) {
         if (d_common->rank == d_owner) {
            ++(d_common->num_boxes_generated);
         }
      } else {
         // Box has no tag.
      }

      // All done with communication within participating group.
      delete d_comm_group;
      d_comm_group = NULL;

   } else {
      /*
       * This process is not in the group that decides on the box for
       * this dendogram node.
       */
      TBOX_ASSERT(!inGroup(d_group));
   }

   if (d_parent == NULL) {
      /*
       * Compute relationships and set up relationship sharing data.
       * This is usually done by a node's parent in the
       * runChildren_check() method because only the
       * parent can know if the node's box will be
       * kept or recombined with the sibling.
       * But the root node must do this itself because it has no parent.
       */
      if (d_common->compute_relationships > 0 && boxAccepted()) {
         computeNewNeighborhoodSets();
      }
   }

   TBOX_ASSERT(d_lft_child == NULL);
   TBOX_ASSERT(d_rht_child == NULL);
   // TBOX_ASSERT( ! inRelaunchQueue(this) );
   TBOX_ASSERT(inRelaunchQueue(this) == d_common->relaunch_queue.end());

   /*
    * Broadcast the result to dropouts.
    * Dropout processes are those that participated in the
    * parent but not in this dendogram node.  They need the
    * result to perform combined efficiency check for the
    * parent.
    *
    * Processes that should participate in the dropout broadcast
    * are the dropouts (processes with zero overlap) and the owner.
    *
    * Broadcast to dropouts is only needed if:
    *
    *    - In multi-owner mode and relationship-computing mode.
    *      In single-owner mode, only the original owner needs
    *      the final result, and it participates everywhere,
    *      so there is no need for this phase.
    *      When computing relationships, participant processors must
    *      know results to do recombination check, to determine
    *      if parent box is preferred.
    *
    *    - This is NOT the root dendogram node.  The root node
    *      has no parent and no corresponding dropout group.
    *
    *    - Dropout group is not empty.  Number of dropouts
    *      is the difference between parent group size and this
    *      group size.
    */
   if (d_overlap == 0 || d_common->rank == d_owner) {

      if ((d_common->owner_mode != SINGLE_OWNER ||
           d_common->compute_relationships > 0) &&
          d_parent != NULL &&
          d_parent->d_group.size() > d_group.size()) {

         d_common->t_bcast_to_dropouts->start();
         {
            // Create the communication group for the dropouts.
            VectorOfInts dropouts(0);
            d_common->t_local_tasks->start();
            computeDropoutGroup(d_parent->d_group,
               d_group,
               dropouts,
               d_owner);
            d_comm_group = new tbox::AsyncCommGroup(
                  computeCommunicationTreeDegree(
                     static_cast<int>(d_group.size())),
                  &d_common->comm_stage,
                  this);
            d_comm_group->setGroupAndRootIndex(d_common->mpi_object,
               &dropouts[0], static_cast<int>(dropouts.size()), 0);
            d_common->t_local_tasks->stop();
         }

         broadcastToDropouts_start();
BCAST_TO_DROPOUTS:
         if (!d_common->t_bcast_to_dropouts->isRunning())
            d_common->t_bcast_to_dropouts->start();
         sub_completed = broadcastToDropouts_check();
         d_common->t_bcast_to_dropouts->stop();

         if (!sub_completed) {
            d_wait_phase = bcast_to_dropouts;
            goto RETURN;
         }

         if (d_common->log_node_history && d_common->rank != d_owner) {
            tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                       << d_common->num_nodes_active << "-active  "
                       << d_common->num_nodes_owned << "-owned  "
                       << d_common->num_nodes_completed << "-completed  "
                       << "DO Recv " << d_generation << ':' << d_pos
                       << "  " << d_mapped_box
                       << "  accept=" << d_box_acceptance
                       << ".\n";
         }

         delete d_comm_group;
         d_comm_group = NULL;
      }
   }

   d_wait_phase = completed;

   if (d_comm_group != NULL) {
      // No further communication.  Deallocate the communication group.
      delete d_comm_group;
      d_comm_group = NULL;
   }

   TBOX_ASSERT(d_common->num_nodes_owned >= 0);

   // Adjust counters.
   --(d_common->num_nodes_active);
   ++(d_common->num_nodes_completed);
   if (d_owner == d_common->rank) {
      TBOX_ASSERT(d_common->num_nodes_owned > 0);
      --(d_common->num_nodes_owned);
      TBOX_ASSERT(d_common->num_nodes_owned >= 0);
   }
   d_common->num_conts_to_complete += d_n_cont;
   if (d_common->max_conts_to_complete < d_n_cont) {
      d_common->max_conts_to_complete = d_n_cont;
   }

   TBOX_ASSERT(d_common->num_nodes_owned >= 0);

   if (d_common->log_node_history) {
      tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                 << d_common->num_nodes_active << "-active  "
                 << d_common->num_nodes_owned << "-owned  "
                 << d_common->num_nodes_completed << "-completed  "
                 << "Complete " << d_generation << ':' << d_pos
                 << "  " << d_mapped_box
                 << "  accept=" << d_box_acceptance
                 << ".\n";
   }

   /*
    * Recall that an dendogram node waiting for its children
    * is not placed in the relaunch queue (because it is
    * pointless to relaunch it until the children are completed).
    * Therefore, to eventually continue that node, its last
    * child to complete must put it on the queue.  If this node
    * and its sibling are completed, put the parent on the FRONT
    * queue to be checked immediately (required for synchronous
    * mode).
    */
   if (d_parent != NULL &&
       d_parent->d_lft_child->d_wait_phase == completed &&
       d_parent->d_rht_child->d_wait_phase == completed) {
      TBOX_ASSERT(d_parent->d_wait_phase == run_children);
      // TBOX_ASSERT( ! inRelaunchQueue(d_parent) );
      TBOX_ASSERT(inRelaunchQueue(d_parent) == d_common->relaunch_queue.end());
      // d_common->relaunch_queue.addItem(d_parent);
      d_common->relaunch_queue.push_front(d_parent);
      if (d_common->log_node_history) {
         tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                    << d_common->num_nodes_active << "-active  "
                    << d_common->num_nodes_owned << "-owned  "
                    << d_common->num_nodes_completed << "-completed  "
                    << "Parent " << d_parent->d_generation << ':'
                    << d_parent->d_pos
                    << " awoken by last child of "
                    << d_parent->d_lft_child->d_generation << ':'
                    << d_parent->d_lft_child->d_pos
                    << ", "
                    << d_parent->d_rht_child->d_generation << ':'
                    << d_parent->d_rht_child->d_pos
                    << " queue size " << d_common->relaunch_queue.size()
                    << ".\n";
      }
   }

RETURN:

#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_wait_phase != completed && d_wait_phase != run_children) {
      TBOX_ASSERT(!d_comm_group->isDone());
      TBOX_ASSERT(d_common->comm_stage.hasPendingRequests());
   }
   if (d_wait_phase == run_children) {
      // TBOX_ASSERT( ! d_common->relaunch_queue.isEmpty() );
      TBOX_ASSERT(!d_common->relaunch_queue.empty());
   }
#endif

   d_common->t_continue_algorithm->stop();

   return d_wait_phase;
}

void
BergerRigoutsosNode::runChildren_start()
{
   /*
    * Children were created to store temporary data
    * and determine participation. Now, run them.
    */

   /*
    * Should only be here if box is rejected based on calculation.
    */
   TBOX_ASSERT(d_box_acceptance == rejected_by_calculation ||
      d_box_acceptance == rejected_by_owner);

   d_lft_child->d_wait_phase = to_be_launched;
   d_rht_child->d_wait_phase = to_be_launched;

   /*
    * Queue the children so they get executed.
    * Put them at the front so that in synchronous
    * mode, they can complete first before moving
    * to another task (important in synchronous mode).
    * It also does not hurt to put children at the
    * front of the queue because they have
    * immediate computation (compute histogram)
    * to perform.  Put the left child in front
    * of the right to more closely match the
    * progression of the recursive BR (not essential).
    */
   d_common->relaunch_queue.push_front(d_rht_child);
   d_common->relaunch_queue.push_front(d_lft_child);
}


/*
 ********************************************************************
 * Check for combined tolerance.
 * If both children accepted their boxes without further splitting
 * but their combined efficiency is not good enough to make
 * the splitting worth accepting, use the current box instead
 * of the children boxes.  Otherwise, use the children boxes.
 ********************************************************************
 */
bool
BergerRigoutsosNode::runChildren_check()
{
   if (d_lft_child->d_wait_phase != completed ||
       d_rht_child->d_wait_phase != completed) {
      return false;
   }

   if (d_lft_child->boxAccepted() &&
       d_rht_child->boxAccepted() &&
       d_box.numberCells() <= d_common->max_box_size &&
       (d_lft_child->d_box.size() + d_rht_child->d_box.size() >=
        d_common->combine_tol * d_box.size())) {

      // Discard childrens' graph nodes in favor of recombination.

      d_box_acceptance = accepted_by_recombination;

      if (d_common->log_node_history) {
         tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                    << d_common->num_nodes_active << "-active  "
                    << d_common->num_nodes_owned << "-owned  "
                    << d_common->num_nodes_completed << "-completed  "
                    << "Recombine " << d_generation << ':' << d_pos
                    << "  " << d_mapped_box
                    << " <= " << d_lft_child->d_mapped_box
                    << " + " << d_rht_child->d_mapped_box
                    << "  " << "accept=" << d_box_acceptance
                    << ".\n";
      }

      if (d_lft_child->d_owner == d_common->rank) {
         d_lft_child->eraseBox();
         d_lft_child->d_box_acceptance = rejected_by_recombination;
         --(d_common->num_boxes_generated);
      }

      if (d_rht_child->d_owner == d_common->rank) {
         d_rht_child->eraseBox();
         d_rht_child->d_box_acceptance = rejected_by_recombination;
         --(d_common->num_boxes_generated);
      }

      if (d_owner == d_common->rank) {
         ++(d_common->num_boxes_generated);
      }

   } else {

      // Accept childrens' results, discarding graph node.

      if (d_owner == d_common->rank) {
         eraseBox();
      }
      if (d_common->compute_relationships > 0) {
         if (d_lft_child->boxAccepted() &&
             d_lft_child->d_box_acceptance != accepted_by_dropout_bcast) {
            d_lft_child->computeNewNeighborhoodSets();
         }
         if (d_rht_child->boxAccepted() &&
             d_rht_child->d_box_acceptance != accepted_by_dropout_bcast) {
            d_rht_child->computeNewNeighborhoodSets();
         }
         if (d_common->log_node_history) {
            tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                       << d_common->num_nodes_active << "-active  "
                       << d_common->num_nodes_owned << "-owned  "
                       << d_common->num_nodes_completed << "-completed  "
                       << "Discard " << d_generation << ':' << d_pos
                       << "  " << d_mapped_box
                       << " => " << d_lft_child->d_mapped_box
                       << " + " << d_rht_child->d_mapped_box
                       << "  " << "accept=" << d_box_acceptance
                       << ".\n";
         }
      }

   }

   /*
    * No longer need children nodes after this point.
    */
   delete d_lft_child;
   delete d_rht_child;
   d_lft_child = 0;
   d_rht_child = 0;

   return true;
}

BergerRigoutsosNode::CommonParams::CommonParams(
   const tbox::Dimension& dim):
   d_dim(dim),
   tag_level((hier::PatchLevel*)NULL),
   tag_mapped_box_level(NULL),
   new_mapped_box_level(NULL),
   tag_to_new(NULL),
   new_to_tag(NULL),
   // Parameters not from clustering algorithm interface ...
   max_lap_cut_from_center(1.0),
   laplace_cut_threshold_ar(0.0),
   max_box_size(d_dim, tbox::MathUtilities<int>::getMax()),
   // Parameters from clustering algorithm interface ...
   tag_data_index(-1),
   tag_val(0),
   min_box(d_dim, 1),
   efficiency_tol(tbox::MathUtilities<double>::getMax()),
   combine_tol(tbox::MathUtilities<double>::getMax()),
   // Implementation flags and data...
   compute_relationships(2),
   relationship_senders(),
   relationship_messages(),
   max_gcw(d_dim, 1),
   owner_mode(MOST_OVERLAP),
   // Communication parameters ...
   mpi_object(tbox::SAMRAI_MPI::commNull),
   rank(-1),
   nproc(-1),
   available_mpi_tag(-1),
   // Analysis support ...
   log_node_history(false),
   num_tags_in_all_nodes(0),
   max_tags_owned(0),
   num_nodes_allocated(0),
   max_nodes_allocated(0),
   num_nodes_active(0),
   max_nodes_active(0),
   num_nodes_owned(0),
   max_nodes_owned(0),
   num_nodes_completed(0),
   max_generation(0),
   num_boxes_generated(0),
   num_conts_to_complete(0),
   max_conts_to_complete(0)
{
   // Set the timer for the communication stage's MPI waiting.
   comm_stage.setCommunicationWaitTimer(t_MPI_wait);
}

/*
 ********************************************************************
 *
 * Asynchronous methods: these methods have _start and _check
 * suffices.  They involve initiating some task and checking
 * whether that task is completed.
 *
 ********************************************************************
 */

void
BergerRigoutsosNode::reduceHistogram_start()
{
   if (d_group.size() == 1) {
      return;
   }
   d_comm_group->setMPITag(d_mpi_tag + reduce_histogram_tag);
   const int hist_size = getHistogramBufferSize(d_box);
   if (d_common->rank == d_owner) {
      d_recv_msg.resize(hist_size, BAD_INTEGER);
      putHistogramToBuffer(&d_recv_msg[0]);
      d_comm_group->beginSumReduce(&d_recv_msg[0], hist_size);
   } else {
      d_send_msg.resize(hist_size, BAD_INTEGER);
      putHistogramToBuffer(&d_send_msg[0]);
      d_comm_group->beginSumReduce(&d_send_msg[0], hist_size);
   }
}

bool
BergerRigoutsosNode::reduceHistogram_check()
{
   if (d_group.size() == 1) {
      return true;
   }
   d_comm_group->proceedToNextWait();
   if (d_comm_group->isDone() && d_common->rank == d_owner) {
      getHistogramFromBuffer(&d_recv_msg[0]);
   }
   return d_comm_group->isDone();
}

void
BergerRigoutsosNode::broadcastAcceptability_start()
{
   if (d_group.size() == 1) {
      return;
   }
   d_comm_group->setMPITag(d_mpi_tag + bcast_acceptability_tag);
   /*
    * Items communicated:
    * - local index of node
    * - whether box is accepted
    * - in case box is accepted:
    *   . box (which may have been trimmed to minimal tag bounding box)
    * - in case box is rejected:
    *   . left/right child boxes
    *   . left/right child MPI tags
    */

   const int buffer_size = 1          // Number of tags in candidate
      + 1                             // Acceptability flag.
      + 1                             // Local index of node.
      + d_dim.getValue() * 2       // Box.
      + d_dim.getValue() * 4       // Children boxes.
      + 2                             // Children MPI tags
   ;

   if (d_common->rank == d_owner) {
      TBOX_ASSERT(d_box_acceptance == rejected_by_calculation ||
         d_box_acceptance == accepted_by_calculation ||
         (d_parent == NULL && d_box_acceptance == hasnotag_by_owner));
      d_send_msg.resize(buffer_size, BAD_INTEGER);
      int* ptr = &d_send_msg[0];
      *(ptr++) = d_num_tags;
      *(ptr++) = d_box_acceptance >= 0 ?
         d_box_acceptance + 2 /* indicate remote decision */ :
         d_box_acceptance;
      if (!boxHasNoTag()) {
         *(ptr++) = d_mapped_box.getLocalId().getValue();
         ptr = putBoxToBuffer(d_box, ptr);
         if (boxRejected()) {
            ptr = putBoxToBuffer(d_lft_child->d_box, ptr);
            ptr = putBoxToBuffer(d_rht_child->d_box, ptr);
            *(ptr++) = d_lft_child->d_mpi_tag;
            *(ptr++) = d_rht_child->d_mpi_tag;
         }
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      else {
         // This may not be needed now that the messages are in vector<int>.
         // Suppress memory check warnings about uninitialized data.
         for (size_t c = ptr - (&d_send_msg[0]); c < d_send_msg.size(); ++c) {
            d_send_msg[c] = -1;
         }
      }
#endif
      d_comm_group->beginBcast(&d_send_msg[0], buffer_size);
   } else {
      d_recv_msg.resize(buffer_size, BAD_INTEGER);
      d_comm_group->beginBcast(&d_recv_msg[0], buffer_size);
   }
}

bool
BergerRigoutsosNode::broadcastAcceptability_check()
{
   if (d_group.size() == 1) {
      return true;
   }
   d_comm_group->checkBcast();
   if (d_comm_group->isDone() && d_common->rank != d_owner) {

      int* ptr = &d_recv_msg[0];
      d_num_tags = *(ptr++);
      if (d_parent == NULL) {
         d_common->num_tags_in_all_nodes = d_num_tags;
      }
      d_box_acceptance = intToBoxAcceptance(*(ptr++));
      TBOX_ASSERT(boxAccepted() || boxRejected() ||
         (boxHasNoTag() && d_parent == NULL));
      if (!boxHasNoTag()) {
         const hier::LocalId node_local_id(*(ptr++));
         ptr = getBoxFromBuffer(d_box, ptr);
         d_mapped_box = hier::Box(d_box, node_local_id, d_owner);
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(d_mapped_box.getBlockId() == d_block_id);
         TBOX_ASSERT(d_mapped_box.getLocalId() >= 0);
         if (d_parent != NULL) {
            /*
             * Do not check for min_box violation in root node.  That
             * check should be done outside of this class in order to
             * have flexibility regarding how to handle it.
             */
            TBOX_ASSERT(d_box.numberCells() >= d_common->min_box);
         }
#endif
      }

      if (boxRejected()) {

         /*
          * The owner formed its children earlier so it can
          * use their parameters while determining which to run.
          * Contributors create the children when the receive
          * the d_box_acceptance flag indicates that further
          * branching is required.
          */
         d_lft_child = new BergerRigoutsosNode(d_common,
               this,
               0,
               d_block_id,
               d_first_local_id);
         d_rht_child = new BergerRigoutsosNode(d_common,
               this,
               1,
               d_block_id,
               d_first_local_id);

         ptr = getBoxFromBuffer(d_lft_child->d_box, ptr);
         ptr = getBoxFromBuffer(d_rht_child->d_box, ptr);
         d_lft_child->d_box.setBlockId(d_block_id);
         d_rht_child->d_box.setBlockId(d_block_id);

         d_lft_child->d_mpi_tag = *(ptr++);
         d_rht_child->d_mpi_tag = *(ptr++);

         TBOX_ASSERT(d_lft_child->d_box.numberCells() >= d_common->min_box);
         TBOX_ASSERT(d_rht_child->d_box.numberCells() >= d_common->min_box);
         TBOX_ASSERT(d_lft_child->d_mpi_tag > -1);
         TBOX_ASSERT(d_rht_child->d_mpi_tag > -1);
         if (d_common->log_node_history) {
            tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                       << d_common->num_nodes_active << "-active  "
                       << d_common->num_nodes_owned << "-owned  "
                       << d_common->num_nodes_completed << "-completed  "
                       << "Rm Split " << d_generation << ':' << d_pos
                       << "  " << d_mapped_box
                       << " => " << d_lft_child->d_box
                       << " + " << d_rht_child->d_box
                       << ".\n";
         }

      }
      else {
         if (d_common->log_node_history) {
            tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                       << d_common->num_nodes_active << "-active  "
                       << d_common->num_nodes_owned << "-owned  "
                       << d_common->num_nodes_completed << "-completed  "
                       << "Rm Accepted " << d_generation << ':' << d_pos
                       << "  " << d_box
                       << "  accept=" << d_box_acceptance << ".\n";
         }
      }
   }
   return d_comm_group->isDone();
}

void
BergerRigoutsosNode::gatherGroupingCriteria_start()
{
   if (d_group.size() == 1) {
      return;
   }
   d_comm_group->setMPITag(d_mpi_tag + gather_grouping_criteria_tag);

   if (d_common->rank == d_owner) {
      d_recv_msg.resize(4 * d_group.size(), BAD_INTEGER);
      d_comm_group->beginGather(&d_recv_msg[0], 4);
   } else {
      d_send_msg.resize(4, BAD_INTEGER);
      d_send_msg[0] = d_lft_child->d_overlap;
      d_send_msg[1] = d_rht_child->d_overlap;
      // Use negative burden measures for uniformity of criteria comparison.
      d_send_msg[2] = -d_common->num_nodes_owned;
      d_send_msg[3] = -d_common->num_nodes_active;
      d_comm_group->beginGather(&d_send_msg[0], 4);
   }
}

void
BergerRigoutsosNode::broadcastChildGroups_start()
{
   if (d_group.size() == 1) {
      return;
   }
   /*
    * Items communicated:
    * - left/right owner
    * - left/right group
    */
   d_comm_group->setMPITag(d_mpi_tag + bcast_child_groups_tag);

   if (d_common->rank == d_owner) {

      /*
       * When d_parent == NULL, use d_comm_group's MPI collective call option.
       * The option uses MPI_Bcast, which requires the buffer size is the same
       * on all processors.  When this is not the case, use the child group
       * sizes to save memory and possibly improve performance.
       */
      const int buffer_size = 2                // Left/right owners.
         + 2                                   // Left/right group sizes.
         + (d_parent == NULL ? static_cast<int>(d_group.size())
            : static_cast<int>(d_lft_child->d_group.size()))    // Left group.
         + (d_parent == NULL ? static_cast<int>(d_group.size())
            : static_cast<int>(d_rht_child->d_group.size()))    // Right group.
      ;

      d_send_msg.resize(buffer_size, BAD_INTEGER);
      int* ptr = &d_send_msg[0];

      *(ptr++) = d_lft_child->d_owner;
      *(ptr++) = static_cast<int>(d_lft_child->d_group.size());
      for (size_t i = 0; i < d_lft_child->d_group.size(); ++i) {
         *(ptr++) = d_lft_child->d_group[i];
      }
      *(ptr++) = d_rht_child->d_owner;
      *(ptr++) = static_cast<int>(d_rht_child->d_group.size());
      for (size_t i = 0; i < d_rht_child->d_group.size(); ++i) {
         *(ptr++) = d_rht_child->d_group[i];
      }
      if (d_parent == NULL) {
         // Initialize unused data to avoid warnings and weird numbers.
         for (size_t i =
                 (d_lft_child->d_group.size() + d_rht_child->d_group.size());
              i < 2 * d_group.size(); ++i) {
            *(ptr++) = -1;
         }
      }

      d_comm_group->beginBcast(&d_send_msg[0], buffer_size);
   } else {
      const int buffer_size = 2                // Left/right owners.
         + 2                                   // Left/right group sizes.
         + 2 * static_cast<int>(d_group.size())   // Left/right groups.
      ;
      d_recv_msg.resize(buffer_size, BAD_INTEGER);

      d_comm_group->beginBcast(&d_recv_msg[0], buffer_size);
   }
}

bool
BergerRigoutsosNode::broadcastChildGroups_check()
{
   if (d_group.size() == 1) {
      return true;
   }
   d_comm_group->checkBcast();
   if (d_comm_group->isDone() && d_common->rank != d_owner) {

      int* ptr = &d_recv_msg[0];

      d_lft_child->d_owner = *(ptr++);
      d_lft_child->d_group.resize(*(ptr++), BAD_INTEGER);
      for (size_t i = 0; i < d_lft_child->d_group.size(); ++i) {
         d_lft_child->d_group[i] = *(ptr++);
      }
      d_rht_child->d_owner = *(ptr++);
      d_rht_child->d_group.resize(*(ptr++), BAD_INTEGER);
      for (size_t i = 0; i < d_rht_child->d_group.size(); ++i) {
         d_rht_child->d_group[i] = *(ptr++);
      }

      TBOX_ASSERT(d_lft_child->d_owner >= 0);
      TBOX_ASSERT(d_lft_child->d_group.size() > 0);
      TBOX_ASSERT((d_lft_child->d_overlap > 0) ==
         inGroup(d_lft_child->d_group));
      TBOX_ASSERT(d_rht_child->d_owner >= 0);
      TBOX_ASSERT(d_rht_child->d_group.size() > 0);
      TBOX_ASSERT((d_rht_child->d_overlap > 0) ==
         inGroup(d_rht_child->d_group));

   }

   return d_comm_group->isDone();
}

void
BergerRigoutsosNode::broadcastToDropouts_start()
{
   TBOX_ASSERT(d_common->rank == d_owner || d_overlap == 0);
   d_comm_group->setMPITag(d_mpi_tag + bcast_to_dropouts_tag);

   const int buffer_size = 1      // d_box_acceptance
      + 1                         // local index of graph node
      + d_dim.getValue() * 2   // d_box (in case it got reduced)
   ;
   d_send_msg.clear();
   d_recv_msg.clear();
   if (d_common->rank == d_owner) {
      d_send_msg.resize(buffer_size, BAD_INTEGER);
      d_send_msg[0] = d_box_acceptance;
      d_send_msg[1] = d_mapped_box.getLocalId().getValue();
      putBoxToBuffer(d_box, &d_send_msg[2]);
      d_comm_group->beginBcast(&d_send_msg[0],
         buffer_size);
   } else {
      d_recv_msg.resize(buffer_size, BAD_INTEGER);
      d_comm_group->beginBcast(&d_recv_msg[0],
         buffer_size);
   }
}

bool
BergerRigoutsosNode::broadcastToDropouts_check()
{
   TBOX_ASSERT(d_common->rank == d_owner || d_overlap == 0);
   d_comm_group->checkBcast();
   if (d_comm_group->isDone()) {
      if (d_common->rank != d_owner) {
         /*
          * We check for the case of the box having no tags,
          * to keeps things explicit and help detect bugs.
          * But in fact, having no tags is impossible
          * in the broadcastToDropout step, because it is
          * only possible for the root dendogram node,
          * which has no dropout group.
          */
         TBOX_ASSERT(d_recv_msg[0] >= 0);

         d_box_acceptance = intToBoxAcceptance((d_recv_msg[0] % 2)
               + rejected_by_dropout_bcast);
         const hier::LocalId local_id(d_recv_msg[1]);
         getBoxFromBuffer(d_box, &d_recv_msg[2]);
#ifdef DEBUG_CHECK_ASSERTIONS
         if (d_parent != NULL) {
            /*
             * Do not check for min_box violation in root node.  That
             * check should be done outside of this class in order to
             * have flexibility regarding how to handle it.
             */
            TBOX_ASSERT(d_box.numberCells() >= d_common->min_box);
         }
#endif
         d_mapped_box = hier::Box(d_box, local_id, d_owner);
         TBOX_ASSERT(d_mapped_box.getBlockId() == d_block_id);
      }
   }
   return d_comm_group->isDone();
}

/*
 ********************************************************************
 * Utility computations using local data.
 ********************************************************************
 */

void
BergerRigoutsosNode::makeLocalTagHistogram()
{
   d_common->t_local_histogram->start();

   /*
    * Compute the histogram size and allocate space for it.
    */
   for (int d = 0; d < d_dim.getValue(); ++d) {
      TBOX_ASSERT(d_box.numberCells(d) > 0);
      d_histogram[d].clear();
      d_histogram[d].insert(d_histogram[d].end(), d_box.numberCells(d), 0);
   }

   /*
    * Accumulate tag counts in the histogram variable.
    */
   const hier::PatchLevel& tag_level = *d_common->tag_level;
   for (hier::PatchLevel::iterator ip(tag_level.begin());
        ip != tag_level.end(); ++ip) {
      hier::Patch& patch = **ip;

      const hier::BlockId& block_id = patch.getBox().getBlockId();

      if (block_id == d_block_id) {
         const hier::Box intersection = patch.getBox() * d_box;
         const hier::IntVector& lower = d_box.lower();

         if (!(intersection.empty())) {

            boost::shared_ptr<pdat::CellData<int> > tag_data_(
               patch.getPatchData(d_common->tag_data_index),
               boost::detail::dynamic_cast_tag());

            pdat::CellData<int>& tag_data = *tag_data_;

            pdat::CellIterator ciend(intersection, false);
            for (pdat::CellIterator ci(intersection, true);
                 ci != ciend; ++ci) {
               if (tag_data(*ci) == d_common->tag_val) {
                  const hier::IntVector& idx = *ci;
                  for (int d = 0; d < d_dim.getValue(); ++d) {
                     ++(d_histogram[d][idx(d) - lower(d)]);
                  }
               }
            }
         }
      }
   }
   d_common->t_local_histogram->stop();
}

/*
 ********************************************************************
 * Change d_box to that of the minimal bounding box for tags.
 * If d_box is changed, reduce d_histogram to new d_box.
 ********************************************************************
 */
void
BergerRigoutsosNode::computeMinimalBoundingBoxForTags()
{
   TBOX_ASSERT(!d_box.empty());

   int d;

   hier::Index new_lower = d_box.lower();
   hier::Index new_upper = d_box.upper();

   const hier::IntVector& min_box = d_common->min_box;
   hier::IntVector box_size = d_box.numberCells();

   /*
    * Bring the lower side of the box up past untagged index planes.
    * Bring the upper side of the box down past untagged index planes.
    * Do not make the box smaller than the min_box requirement.
    */
   for (d = 0; d < d_dim.getValue(); ++d) {
      TBOX_ASSERT(d_histogram[d].size() != 0);
      int* histogram_beg = &d_histogram[d][0];
      int* histogram_end = histogram_beg + d_box.numberCells(d) - 1;
      while (*histogram_beg == 0 &&
             box_size(d) > min_box(d)) {
         ++new_lower(d);
         ++histogram_beg;
         --box_size(d);
      }
      while (*histogram_end == 0 &&
             box_size(d) > min_box(d)) {
         --new_upper(d);
         --histogram_end;
         --box_size(d);
      }
   }

   const hier::Box new_box(new_lower, new_upper, d_block_id);
   const hier::IntVector new_size = new_box.numberCells();

   if (!new_box.isSpatiallyEqual(d_box)) {
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_parent != NULL) {
         /*
          * Do not check for min_box violation in root node.  That
          * check should be done outside of this class in order to
          * have flexibility regarding how to handle it.
          */
         TBOX_ASSERT(new_box.numberCells() >= d_common->min_box);
      }
#endif
      /*
       * Save tagged part of the current histogram and reset the box.
       * Is this step really required?  No, we can just keep the
       * shift in a hier::IntVector and adjust.
       */
      for (d = 0; d < d_dim.getValue(); ++d) {
         VectorOfInts& h = d_histogram[d];
         const int shift = new_lower(d) - d_box.lower() (d);
         if (shift > 0) {
            int i;
            for (i = 0; i < new_size(d); ++i) {
               h[i] = h[i + shift];
            }
         }
         h.resize(new_size(d), BAD_INTEGER);
      }
      if (d_common->log_node_history) {
         tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                    << d_common->num_nodes_active << "-active  "
                    << d_common->num_nodes_owned << "-owned  "
                    << d_common->num_nodes_completed << "-completed  "
                    << "Shrunken " << d_generation << ':' << d_pos
                    << "  " << d_box << " -> " << new_box
                    << ".\n";
      }
      d_box = new_box;
   }
}

/*
 *********************************************************************
 * Accept the box or split it, setting d_box_acceptance accordingly.
 *********************************************************************
 */
void
BergerRigoutsosNode::acceptOrSplitBox()
{

#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_owner != d_common->rank) {
      TBOX_ERROR("Only the owner can determine\n"
         "whether to accept or split a box.\n");
   }
   TBOX_ASSERT(d_box_acceptance == undetermined);
#endif

   const hier::IntVector boxdims(d_box.numberCells());
   const hier::IntVector oversize(boxdims - d_common->max_box_size);

   /*
    * Box d_box is acceptable if
    * - it has a high enough fraction of tagged cells, or
    * - it cannot be split without breaking the minimum
    *   box requirement, or
    *
    * If d_box has no tags:
    * - set d_box_acceptance = hasnotag_by_owner
    * If accepting d_box:
    * - set d_box_acceptance = accepted_by_calculation
    * If rejecting d_box:
    * - set d_box_acceptance = rejected_by_calculation
    * - create left and right children
    * - set children boxes
    * - claim MPI tags for communication by children nodes
    *
    * Instead of writing from scratch,
    * the code to find the split plane was copied
    * from mesh::BergerRigoutsos<DIM>::splitTagBoundBox()
    * and modified.
    */

   if (d_box_acceptance == undetermined) {
      if (oversize <= hier::IntVector::getZero(d_dim)) {
         /*
          * See if d_box should be accepted based on efficiency.
          */
         int num_tagged = 0;
         for (size_t i = 0; i < d_histogram[0].size(); ++i) {
            num_tagged += d_histogram[0][i];
         }
         int boxsize = d_box.size();
         double efficiency = (boxsize == 0 ? 1.e0 :
                              ((double)num_tagged) / boxsize);

         if (d_common->max_tags_owned < num_tagged) {
            d_common->max_tags_owned = num_tagged;
         }

         if (efficiency >= d_common->efficiency_tol) {
            d_box_acceptance = accepted_by_calculation;
            if (d_common->log_node_history) {
               tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                          << d_common->num_nodes_active << "-active  "
                          << d_common->num_nodes_owned << "-owned  "
                          << d_common->num_nodes_completed << "-completed  "
                          << "Accepted " << d_generation << ':' << d_pos
                          << "  " << d_box << " by sufficient efficiency of " << efficiency
                          << "  accept=" << d_box_acceptance << ".\n";
            }
         } else if (num_tagged == 0) {
            // No tags!  This should be caught at the dendogram root.
            TBOX_ASSERT(d_parent == NULL);
            d_box_acceptance = hasnotag_by_owner;
            if (d_common->log_node_history) {
               tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                          << d_common->num_nodes_active << "-active  "
                          << d_common->num_nodes_owned << "-owned  "
                          << d_common->num_nodes_completed << "-completed  "
                          << "HasNoTag " << d_generation << ':' << d_pos
                          << "  " << d_box
                          << ".\n";
            }
         }
      }
   }

   /*
    * If cut_margin is negative in any direction, we cannot cut d_box
    * across that direction without violating min_box.
    */
   const hier::IntVector cut_margin = boxdims - (d_common->min_box) * 2;

   if (d_box_acceptance == undetermined) {
      /*
       * If d_box cannot be split without violating min_box, it should
       * be accepted.
       */
      if (cut_margin < hier::IntVector::getZero(d_dim)) {
         d_box_acceptance = accepted_by_calculation;
      }
   }

   hier::IntVector sorted_margins(d_dim);

   if (d_box_acceptance == undetermined) {
      /*
       * Sort the bounding box dimensions from largest to smallest cut
       * margin.  If there are multiple cuttable directions, we will
       * favor the direction with the greatest cut_margin.
       */
      int dim;
      for (dim = 0; dim < d_dim.getValue(); dim++) {
         sorted_margins(dim) = dim;
      }
      for (int d0 = 0; d0 < d_dim.getValue() - 1; d0++) {
         for (int d1 = d0 + 1; d1 < d_dim.getValue(); d1++) {
            if (cut_margin(sorted_margins(d0)) <
                cut_margin(sorted_margins(d1))) {
               int tmp_dim = sorted_margins(d0);
               sorted_margins(d0) = sorted_margins(d1);
               sorted_margins(d1) = tmp_dim;
            }
         }
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      for (dim = 0; dim < d_dim.getValue() - 1; dim++) {
         TBOX_ASSERT(cut_margin(sorted_margins(dim)) >=
            cut_margin(sorted_margins(dim + 1)));
      }
#endif
   }

   const int max_margin_dir = sorted_margins(0);
   const int min_margin_dir = sorted_margins(d_dim.getValue()-1);

   int num_cuttable_dim = 0;

   if (d_box_acceptance == undetermined) {
      /*
       * Determine number of coordinate directions that are cuttable
       * according to the cut_margin.
       */
      for (num_cuttable_dim = 0; num_cuttable_dim < d_dim.getValue();
           num_cuttable_dim++) {
         if (cut_margin(sorted_margins(num_cuttable_dim)) < 0) {
            break;
         }
      }
      TBOX_ASSERT(num_cuttable_dim > 0);   // We already accounted for un-cuttable case before this point.
   }

   if (d_box_acceptance == undetermined) {

      /*
       * Attempt to split box at a zero interior point in the
       * histogram.  Check each cuttable direction, from
       * largest to smallest, until zero point found.
       */

      int cut_lo, cut_hi;
      int cut_pt = -(tbox::MathUtilities<int>::getMax());
      int cut_dir = -1;
      int dir = -1;
      const hier::Index box_lo(d_box.lower());
      const hier::Index box_hi(d_box.upper());
      hier::Index lft_hi(box_hi);
      hier::Index rht_lo(box_lo);

      for (dir = 0; dir < d_dim.getValue(); dir++) {
         cut_dir = sorted_margins(dir);
         if (cut_margin(cut_dir) < 0) {
            continue;  // This direction is too small to cut.
         }
         if (findZeroCutSwath(cut_lo, cut_hi, cut_dir)) {
            // Split bound box at cut_pt; cut_dir is splitting dimension.
            TBOX_ASSERT(cut_hi - cut_lo >= 0);
            lft_hi(cut_dir) = cut_lo - 1;
            rht_lo(cut_dir) = cut_hi + 1;
            if (d_common->log_node_history) {
               tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                          << d_common->num_nodes_active << "-active  "
                          << d_common->num_nodes_owned << "-owned  "
                          << d_common->num_nodes_completed << "-completed  "
                          << "HoleCut " << d_generation << ':' << d_pos
                          << "  " << d_box << " d=" << cut_dir
                          << " at " << cut_lo << '-' << cut_hi
                          << ".\n";
            }
            break;
         }
      }

      /*
       * If no zero point found, try Laplacian cut.
       */

      if (dir == d_dim.getValue()) {

         /*
          * laplace_cut_threshold_ar specifies the mininum box
          * thickness that can be cut, as a ratio to the thinnest box
          * direction.  If the box doesn't have any direction thick
          * enough, then it has a reasonable aspect ratio, so we can
          * cut it in any direction.
          *
          * Degenerate values of laplace_cut_threshold_ar:
          *
          * 1: cut any direction except the thinnest.
          *
          * (0,1) and huge values: cut any direction.
          *
          * 0: Not a degenerate case but a special case meaning cut
          * only the thickest direction.  This leads to more cubic
          * boxes but can miss feature edges aligned across other
          * directions.
          */
         int max_box_length_to_leave = boxdims(max_margin_dir) - 1;
         if ( d_common->laplace_cut_threshold_ar > 0.0 ) {
            max_box_length_to_leave = static_cast<int>(0.5 + boxdims(min_margin_dir)*d_common->laplace_cut_threshold_ar);
            if ( max_box_length_to_leave >= boxdims(max_margin_dir) ) {
               /*
                * Box aspect ratio is not too bad. Disable preference
                * for cutting longer dirs.
                */
               max_box_length_to_leave = 0;
            }
         }

         int diff_laplace = -1;
         for ( int d=0; d<d_dim.getValue(); ++d ) {
            if ( cut_margin(d) < 0 || boxdims(d) <= max_box_length_to_leave ) {
               continue;  // Direction d is too small to cut.
            }
            int try_cut_pt, try_diff_laplace;
            cutAtLaplacian(try_cut_pt, try_diff_laplace, d);
            if ( diff_laplace < try_diff_laplace ||
                 ( diff_laplace == try_diff_laplace && cut_margin(d) > cut_margin(cut_dir) ) ) {
               cut_dir = d;
               cut_pt = try_cut_pt;
               diff_laplace = try_diff_laplace;
            }
         }
         TBOX_ASSERT( cut_dir >= 0 && cut_dir < d_dim.getValue() );

         // Split bound box at cut_pt; cut_dir is splitting dimension.
         lft_hi(cut_dir) = cut_pt - 1;
         rht_lo(cut_dir) = cut_pt;
         if (d_common->log_node_history) {
            tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                       << d_common->num_nodes_active << "-active  "
                       << d_common->num_nodes_owned << "-owned  "
                       << d_common->num_nodes_completed << "-completed  "
                       << "LapCut " << d_generation << ':' << d_pos
                       << "  " << d_box
                       << " d=" << cut_dir << " at " << cut_pt
                       << ".\n";
         }
      }

      /*
       * The owner forms its children now so it can use their
       * parameters while determining which to run.
       * Contributors create the children when they receive
       * the d_box_acceptance flag from the owner.
       */
      d_lft_child = new BergerRigoutsosNode(d_common, this, 0, d_block_id, d_first_local_id);
      d_rht_child = new BergerRigoutsosNode(d_common, this, 1, d_block_id, d_first_local_id);

      d_lft_child->d_box = hier::Box(box_lo, lft_hi, d_block_id);
      d_rht_child->d_box = hier::Box(rht_lo, box_hi, d_block_id);
      TBOX_ASSERT(d_lft_child->d_box.numberCells() >= d_common->min_box);
      TBOX_ASSERT(d_rht_child->d_box.numberCells() >= d_common->min_box);

      d_lft_child->claimMPITag();
      d_rht_child->claimMPITag();

      d_box_acceptance = rejected_by_calculation;

      if (d_common->log_node_history) {
         tbox::plog << d_common->num_nodes_allocated << "-allocated  "
                    << d_common->num_nodes_active << "-active  "
                    << d_common->num_nodes_owned << "-owned  "
                    << d_common->num_nodes_completed << "-completed  "
                    << "Lc Split "
                    << d_generation << ':' << d_pos << "  " << d_box
                    << " => " << d_lft_child->d_generation << ':'
                    << d_lft_child->d_pos << d_lft_child->d_box
                    << " + " << d_rht_child->d_generation << ':'
                    << d_rht_child->d_pos << d_rht_child->d_box
                    << ".\n";
      }

   }
}

/*
 ********************************************************************
 *
 * Attempt to find a range with zero histogram value near the
 * middle of d_box in the given coordinate direction.
 * Note that the hole is kept more than a minimium distance from
 * the endpoints of of the index interval.
 *
 * Note that it is assumed that box indices are cell indices.
 *
 * If a hole is found, cut_lo and cut_hi are set to the
 * range of zero tag cells.
 *
 ********************************************************************
 */

bool
BergerRigoutsosNode::findZeroCutSwath(
   int& cut_lo,
   int& cut_hi,
   const int dim)
{
   const int lo = d_box.lower(dim);
   const int hi = d_box.upper(dim);
   // Compute the limit for the swath.
   const int cut_lo_lim = lo + d_common->min_box(dim);
   const int cut_hi_lim = hi - d_common->min_box(dim);

   /*
    * Start in the middle of the box.
    * Move cut_lo down and cut_hi up until a hole is found.
    * Keep moving in same direction of the hole until the
    * other side of the hole is found.  The two planes form
    * the widest cut possible at the hole.
    */
   cut_lo = cut_hi = (lo + hi) / 2;
   while ((cut_lo >= cut_lo_lim) && (cut_hi <= cut_hi_lim)) {
      if (d_histogram[dim][cut_lo - lo] == 0) {
         /* The narrow cut is at cut_lo.  Initialize the cut swath here
          * and move cut_lo down until the far side the hole is found.
          */
         cut_hi = cut_lo;
         while (((cut_lo > cut_lo_lim)) &&
                (d_histogram[dim][cut_lo - lo - 1] == 0)) {
            --cut_lo;
         }
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(cut_hi >= cut_lo);
         TBOX_ASSERT(cut_lo - lo >= d_common->min_box(dim));
         TBOX_ASSERT(hi - cut_hi >= d_common->min_box(dim));
         for (int i = cut_lo; i <= cut_hi; ++i) {
            TBOX_ASSERT(d_histogram[dim][i - lo] == 0);
         }
#endif
         return true;
      }
      if (d_histogram[dim][cut_hi - lo] == 0) {
         /* The narrow cut is at cut_hi.  Initialize the cut swath here
          * and move cut_hi up until the far side the hole is found.
          */
         cut_lo = cut_hi;
         while (((cut_hi < cut_hi_lim)) &&
                (d_histogram[dim][cut_hi - lo + 1] == 0)) {
            ++cut_hi;
         }
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(cut_hi >= cut_lo);
         TBOX_ASSERT(cut_lo - lo >= d_common->min_box(dim));
         TBOX_ASSERT(hi - cut_hi >= d_common->min_box(dim));
         for (int i = cut_lo; i <= cut_hi; ++i) {
            TBOX_ASSERT(d_histogram[dim][i - lo] == 0);
         }
#endif
         return true;
      }
      --cut_lo;
      ++cut_hi;
   }

   return false;
}

/*
 ***********************************************************************
 *
 * Attempt to find a point in the given coordinate direction near an
 * inflection point in the histogram for that direction. Note that the
 * cut point is kept more than a minimium distance from the endpoints
 * of the index interval (lo, hi).  Also, the box must have at least
 * three cells along a side to apply the Laplacian test.  If no
 * inflection point is found, the mid-point of the interval is
 * returned as the cut point.
 *
 * Note that it is assumed that box indices are cell indices.
 *
 ***********************************************************************
 */

void
BergerRigoutsosNode::cutAtLaplacian(
   int& cut_pt,
   int& diff_laplace,
   const int dim)
{
   /*
    * New implementation prefers and possibly restricts the Laplace cut
    * to the center part of the box.
    *
    * The cuts refer to face indices, not cell indices.
    *
    * Note that we work in the index space centered on the box's lower
    * cell and add the box lower cell index at the end.
    */

   const VectorOfInts& hist = d_histogram[dim];
   const unsigned int hist_size = static_cast<int>(hist.size());
   TBOX_ASSERT(d_box.upper() (dim) - d_box.lower() (dim) + 1 == static_cast<int>(hist_size));
   TBOX_ASSERT(hist_size >= 2);

   /*
    * Laplacian cut requires at least 4 cells of histogram, so it can
    * compare at 2 Laplacians.  Without 4 cells, we just cut across the
    * largest change in the histogram.
    */
   if (hist_size < 4) {
      cut_pt = 1;
      for (unsigned int i = 2; i < hist_size; ++i) {
         if (tbox::MathUtilities<int>::Abs(hist[cut_pt] - hist[cut_pt - 1]) <
             tbox::MathUtilities<int>::Abs(hist[i] - hist[i - 1])) {
            cut_pt = i;
         }
      }
      cut_pt += d_box.lower() (dim);
      diff_laplace = 0;  // Did not use any Laplace values.
      return;
   }

   const int box_lo = 0;
   const int box_hi = hist_size - 1;
   const int max_dist_from_center =
      int(d_common->max_lap_cut_from_center * hist_size / 2);
   const int box_mid = (box_lo + box_hi + 1) / 2;

   const int cut_lo_lim = tbox::MathUtilities<int>::Max(
         box_lo + d_common->min_box(dim),
         box_mid - max_dist_from_center);

   const int cut_hi_lim = tbox::MathUtilities<int>::Min(
         box_hi - d_common->min_box(dim) + 1,
         box_mid + max_dist_from_center);

   /*
    * Initial cut point and differences between the Laplaces on either
    * side of it.  We want to cut where this difference is biggest.
    */
   cut_pt = box_mid;
   diff_laplace =
      (hist[cut_pt - 1] - 2 * hist[cut_pt] + hist[cut_pt + 1])
      - (hist[cut_pt - 2] - 2 * hist[cut_pt - 1] + hist[cut_pt]);
   diff_laplace = tbox::MathUtilities<int>::Abs(diff_laplace);

   int ic = 1;
   int cut_lo = box_mid - ic;
   int cut_hi = box_mid + ic;

   while (cut_lo > cut_lo_lim || cut_hi < cut_hi_lim) {
      if (cut_lo > cut_lo_lim) {
         const int la = (hist[cut_lo - 1] - 2 * hist[cut_lo] + hist[cut_lo + 1]);
         const int lb = (hist[cut_lo - 2] - 2 * hist[cut_lo - 1] + hist[cut_lo]);
         if ( la*lb <= 0 ) {
            const int try_diff_laplace = tbox::MathUtilities<int>::Abs(la-lb);
            if (try_diff_laplace > diff_laplace) {
               cut_pt = cut_lo;
               diff_laplace = try_diff_laplace;
            }
         }
      }
      if (cut_hi < cut_hi_lim) {
         const int la = (hist[cut_hi - 1] - 2 * hist[cut_hi] + hist[cut_hi + 1]);
         const int lb = (hist[cut_hi - 2] - 2 * hist[cut_hi - 1] + hist[cut_hi]);
         if ( la*lb <= 0 ) {
            const int try_diff_laplace = tbox::MathUtilities<int>::Abs(la-lb);
            if (try_diff_laplace > diff_laplace) {
               cut_pt = cut_hi;
               diff_laplace = try_diff_laplace;
            }
         }
      }
      --cut_lo;
      ++cut_hi;
   }

   cut_pt += d_box.lower() (dim);
}

/*
 ********************************************************************
 * Create a DLBG Box in d_common->new_mapped_box_level,
 * where the output boxes of the algorithm is saved.
 *
 * Only the owner should create the mapped_box_level node this way.
 * Other processes build mapped_box_level node using data from owner.
 ********************************************************************
 */
void
BergerRigoutsosNode::createBox()
{
   TBOX_ASSERT(d_common->rank == d_owner);
   hier::LocalId last_index =
      d_common->new_mapped_box_level->getBoxes().isEmpty() ? d_first_local_id :
      d_common->new_mapped_box_level->getBoxes().back().getLocalId();

   hier::Box new_box(d_box, last_index + 1, d_common->rank);
   TBOX_ASSERT(new_box.getBlockId() == d_block_id);
   d_common->new_mapped_box_level->addBoxWithoutUpdate(new_box);
   d_mapped_box_iterator = d_common->new_mapped_box_level->getBox(new_box);

   d_mapped_box = *d_mapped_box_iterator;
}

/*
 ********************************************************************
 * Discard the Box.  On the owner, this Box is a part of
 * d_common->new_mapped_box_level where it must be removed.  On
 * contributors the Box can just be ignored.  To prevent bugs,
 * the node and its iterator are set to unusable values.
 ********************************************************************
 */
void
BergerRigoutsosNode::eraseBox()
{
   if (d_common->rank == d_owner) {
      d_common->new_mapped_box_level->eraseBoxWithoutUpdate(
         *d_mapped_box_iterator);
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   d_mapped_box_iterator = BoxContainer().end();
   d_mapped_box = hier::Box(d_dim);
#endif
}

void
BergerRigoutsosNode::countOverlapWithLocalPatches()
{
   /*
    * Count overlaps for the left and right sides.
    *
    * Remove the child if it has zero overlap.
    */
   hier::Box lft_grown_box = d_lft_child->d_box;
   lft_grown_box.grow(d_common->max_gcw);
   hier::Box rht_grown_box = d_rht_child->d_box;
   rht_grown_box.grow(d_common->max_gcw);
   int& lft_overlap = d_lft_child->d_overlap;
   int& rht_overlap = d_rht_child->d_overlap;
   lft_overlap = rht_overlap = 0;

   const hier::PatchLevel& tag_level = *d_common->tag_level;
   for (hier::PatchLevel::iterator ip(tag_level.begin());
        ip != tag_level.end(); ++ip) {

      const hier::BlockId& block_id = (*ip)->getBox().getBlockId();

      if (block_id == d_block_id) {
         const hier::Box& patch_box = (*ip)->getBox();

         hier::Box lft_intersection = patch_box * lft_grown_box;
         lft_overlap += lft_intersection.size();

         hier::Box rht_intersection = patch_box * rht_grown_box;
         rht_overlap += rht_intersection.size();

      }
   }
}

/*
 *************************************************************************
 * Child groups are subsets of current group.  Each child group
 * includes processes owning patches that overlap the box of that child.
 * The overlap data has been gathered in d_recv_msg.
 * See gatherGroupingCriteria_start() for the format of the message.
 *************************************************************************
 */
void
BergerRigoutsosNode::formChildGroups()
{
   /*
    * Form child groups and determine owners from data gathered
    * in the gather_overlap_counts phase.
    */
   if (d_group.size() == 1) {
      // Short cut for trivial groups.
      d_lft_child->d_group.resize(1, BAD_INTEGER);
      d_rht_child->d_group.resize(1, BAD_INTEGER);
      d_lft_child->d_group[0] = d_group[0];
      d_rht_child->d_group[0] = d_group[0];
      d_lft_child->d_owner = d_owner;
      d_rht_child->d_owner = d_owner;
      return;
   }

   d_lft_child->d_group.resize(d_group.size(), BAD_INTEGER);
   d_rht_child->d_group.resize(d_group.size(), BAD_INTEGER);

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * Only owner process should be here.
    */
   if (d_common->rank != d_owner) {
      TBOX_ERROR("Library error!");
   }
   TBOX_ASSERT(d_recv_msg.size() == 4 * d_group.size());
#endif

   int* lft_overlap = &d_recv_msg[0];
   int* rht_overlap = &d_recv_msg[1];

   const int imyself = findOwnerInGroup(d_common->rank, d_group);
   lft_overlap[imyself * 4] = d_lft_child->d_overlap;
   rht_overlap[imyself * 4] = d_rht_child->d_overlap;

   int* lft_criteria = NULL;
   int* rht_criteria = NULL;
   switch (d_common->owner_mode) {
      case SINGLE_OWNER:
         lft_criteria = &d_recv_msg[0];
         rht_criteria = &d_recv_msg[1];
         lft_criteria[imyself * 4] = tbox::MathUtilities<int>::getMax();
         rht_criteria[imyself * 4] = tbox::MathUtilities<int>::getMax();
         break;
      case MOST_OVERLAP:
         lft_criteria = &d_recv_msg[0];
         rht_criteria = &d_recv_msg[1];
         lft_criteria[imyself * 4] = d_lft_child->d_overlap;
         rht_criteria[imyself * 4] = d_rht_child->d_overlap;
         break;
      case FEWEST_OWNED:
         lft_criteria = &d_recv_msg[2];
         rht_criteria = &d_recv_msg[2];
         lft_criteria[imyself * 4] = -d_common->num_nodes_owned;
         rht_criteria[imyself * 4] = -d_common->num_nodes_owned;
         break;
      case LEAST_ACTIVE:
         lft_criteria = &d_recv_msg[3];
         rht_criteria = &d_recv_msg[3];
         lft_criteria[imyself * 4] = -d_common->num_nodes_active;
         rht_criteria[imyself * 4] = -d_common->num_nodes_active;
         break;
      default:
         TBOX_ERROR("LIBRARY error");
         break;
   }

   int n_lft = 0;
   int n_rht = 0;

   int lft_owner_score = tbox::MathUtilities<int>::getMin();
   int rht_owner_score = tbox::MathUtilities<int>::getMin();

   /*
    * Loop through the group to see which process should participate
    * on the left/right sides.  Also see which process should be the
    * owner of the left/right sides.  For efficiency in some searches
    * through d_groups, make sure that d_group is ordered.
    */
   for (unsigned int i = 0; i < d_group.size(); ++i) {
      int i4 = i * 4;
      if (lft_overlap[i4] != 0) {
         d_lft_child->d_group[n_lft++] = d_group[i];
         if (lft_criteria[i4] > lft_owner_score) {
            d_lft_child->d_owner = d_group[i];
            lft_owner_score = lft_criteria[i4];
         }
      }
      if (rht_overlap[i4] != 0) {
         d_rht_child->d_group[n_rht++] = d_group[i];
         if (rht_criteria[i4] > rht_owner_score) {
            d_rht_child->d_owner = d_group[i];
            rht_owner_score = rht_criteria[i4];
         }
      }
   }

   d_lft_child->d_group.resize(n_lft, BAD_INTEGER);
   d_rht_child->d_group.resize(n_rht, BAD_INTEGER);

#ifdef DEBUG_CHECK_ASSERTIONS
   // Recall that only the owner should execute this code.
   TBOX_ASSERT(d_lft_child->d_owner >= 0);
   TBOX_ASSERT(d_lft_child->d_group.size() > 0);
   TBOX_ASSERT(d_lft_child->d_group.size() <= d_group.size());
   TBOX_ASSERT(d_common->owner_mode == SINGLE_OWNER ||
      ((d_lft_child->d_overlap == 0) !=
       inGroup(d_lft_child->d_group)));
   TBOX_ASSERT(d_rht_child->d_owner >= 0);
   TBOX_ASSERT(d_rht_child->d_group.size() > 0);
   TBOX_ASSERT(d_rht_child->d_group.size() <= d_group.size());
   TBOX_ASSERT(d_common->owner_mode == SINGLE_OWNER ||
      ((d_rht_child->d_overlap == 0) !=
       inGroup(d_rht_child->d_group)));
   if (d_common->owner_mode == SINGLE_OWNER) {
      TBOX_ASSERT(inGroup(d_lft_child->d_group, d_owner));
      TBOX_ASSERT(inGroup(d_rht_child->d_group, d_owner));
   }
   for (size_t i = 0; i < d_group.size(); ++i) {
      TBOX_ASSERT(i == 0 || d_group[i] > d_group[i - 1]);
      TBOX_ASSERT((lft_overlap[i * 4] > 0 ||
                   (d_group[i] == d_lft_child->d_owner))
         == inGroup(d_lft_child->d_group, d_group[i]));
      TBOX_ASSERT((rht_overlap[i * 4] > 0 ||
                   (d_group[i] == d_rht_child->d_owner))
         == inGroup(d_rht_child->d_group, d_group[i]));
   }
#endif
}

/*
 *************************************************************************
 *
 * Compute overlaps between the new graph node and nodes on
 * the tagged level, saving that data in the form of relationships.
 *
 * Note that the relationship data may be duplicated in two objects.
 * - tag_to_new stores the relationships organized around each node
 *   in the tagged level.  For each node on the tagged level,
 *   we store a container of neighbors on the new mapped_box_level.
 * - new_to_tag stores the relationships organized around each NEW node.
 *   For each new node we store a container of neighbors on the
 *   tagged level.
 *
 * If compute_relationships > 0, we store tag_to_new.
 *
 * If compute_relationships > 1, we also compute new_to_tag.
 * The data in new_to_tag are
 * computed by the particant processes but eventually stored on the
 * owners of the new nodes, so their computation requires caching
 * the relationship data in relationship_messages for sending to the appropriate
 * processes later.
 *
 *************************************************************************
 */
void
BergerRigoutsosNode::computeNewNeighborhoodSets()
{
   d_common->t_compute_new_graph_relationships->start();
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(d_common->compute_relationships > 0);
   TBOX_ASSERT(d_mapped_box.getLocalId() >= 0);
   TBOX_ASSERT(boxAccepted());
   TBOX_ASSERT(d_box_acceptance != accepted_by_dropout_bcast);
   if (d_parent != NULL) {
      /*
       * Do not check for min_box violation in root node.  That
       * check should be done outside of this class in order to
       * have flexibility regarding how to handle it.
       */
      TBOX_ASSERT(d_box.numberCells() >= d_common->min_box);
   }
   /*
    * We should not compute nabrs if we got the node
    * by a dropout broadcast because we already know
    * there is no overlap!
    */
   TBOX_ASSERT(d_box_acceptance != accepted_by_dropout_bcast);
#endif

   // Create an expanded box for intersection check.
   hier::Box grown_box = d_box;
   grown_box.grow(d_common->max_gcw);

   /*
    * On the owner process, we store the neighbors of the new node.
    * This data is NOT required on other processes.
    */
   bool on_owner_process = d_common->rank == d_owner;
   if (on_owner_process) {
      d_common->new_to_tag->makeEmptyLocalNeighborhood(d_mapped_box.getId());
   }

   // Data to send to d_owner regarding new relationships found by local process.
   VectorOfInts* relationship_message = NULL;
   if (d_common->compute_relationships > 1 && d_common->rank != d_owner) {
      /*
       * Will have to send to d_owner the relationships found locally for
       * graph node d_mapped_box.
       * Label the id of the new node and the (yet unknown) number
       * of relationship found for it.
       *
       * The message to be sent to owner is appended the following
       * data:
       * - index of new node
       * - number of relationships found for the new node
       * - index of nodes on the tagged level overlapping new node.
       */
      relationship_message = &d_common->relationship_messages[d_owner];
      relationship_message->insert(relationship_message->end(), d_mapped_box.getLocalId().getValue());
      relationship_message->insert(relationship_message->end(), 0);
   }

   const int index_of_counter =
      (relationship_message != NULL ? static_cast<int>(relationship_message->size()) : 0) - 1;
   const int ints_per_node = hier::Box::commBufferSize(d_dim);

   const BoxContainer& tag_mapped_boxes =
      d_common->tag_mapped_box_level->getBoxes();

   for (hier::RealBoxConstIterator ni(tag_mapped_boxes.realBegin());
        ni != tag_mapped_boxes.realEnd(); ++ni) {

      const hier::Box& tag_mapped_box = *ni;

      if (tag_mapped_box.getBlockId() == d_block_id) {
         hier::Box intersection = tag_mapped_box * grown_box;

         if (!intersection.empty()) {

            // Add d_mapped_box as a neighbor of tag_mapped_box.
            d_common->tag_to_new->insertLocalNeighbor(d_mapped_box,
               tag_mapped_box.getId());

            if (on_owner_process) {
               // Owner adds tag_mapped_box as a neighbor of d_mapped_box.
               d_common->new_to_tag->insertLocalNeighbor(tag_mapped_box,
                  d_mapped_box.getId());
            }

            if (relationship_message != NULL) {
               /* Non-owners put found relationship in the message
                * to (eventually) send to d_owner.
                */
               relationship_message->insert(relationship_message->end(), ints_per_node, 0);
               int* ptr = &(*relationship_message)[relationship_message->size() - ints_per_node];
               tag_mapped_box.putToIntBuffer(ptr);
               ++(*relationship_message)[index_of_counter];
            }
         }
      }
   }

   if (d_common->compute_relationships > 1 &&
       d_common->rank == d_owner) {
      /*
       * If box was accepted, the owner should remember
       * which process will be sending relationship data.
       * Update the list of relationship senders to make sure
       * it includes all processes in the group.
       * We use this list in shareNewNeighborhoodSetsWithOwners to
       * tell us which processors are sending us new relationships.
       * The relationship senders are the participants of the group.
       */
      d_common->relationship_senders.insert(d_group.begin(), d_group.end());
   }

   d_common->t_compute_new_graph_relationships->stop();
}

/*
 **********************************************************************
 *
 * Send new relationships found by local process to owners of the new nodes
 * associated with those relationships.  Receive similar data from other
 * processes.
 *
 * Messages to be sent out were placed in d_common->relationship_messages by
 * computeNewNeighborhoodSets().  This method sends out these messages
 * and receives anticipated messages from processes listed in
 * d_common->relationship_senders.  Received messages are unpacked to get
 * data on new relationships.
 *
 **********************************************************************
 */
void
BergerRigoutsosNode::shareNewNeighborhoodSetsWithOwners()
{
   tbox::SAMRAI_MPI mpi(d_common->mpi_object);
   if (mpi.getSize() == 1) {
      return;
   }

   d_common->t_share_new_relationships->start();

   IntSet relationship_senders = d_common->relationship_senders;
   std::map<int, VectorOfInts>& relationship_messages = d_common->relationship_messages;

   const int ints_per_node = hier::Box::commBufferSize(d_dim);

   int ierr;
   tbox::SAMRAI_MPI::Status mpi_status;

   // Nonblocking send of relationship data.
   d_common->t_share_new_relationships_send->start();
   tbox::Array<tbox::SAMRAI_MPI::Request> mpi_request(
      static_cast<int>(relationship_messages.size()));
   std::map<int, VectorOfInts>::iterator send_i;
   int nsend = 0;
   for (send_i = relationship_messages.begin(), nsend = 0;
        send_i != relationship_messages.end();
        ++send_i, ++nsend) {
      const int& owner = (*send_i).first;
      VectorOfInts& msg = (*send_i).second;
      ierr = mpi.Isend(&msg[0],
            static_cast<int>(msg.size()),
            MPI_INT,
            owner,
            d_common->tag_upper_bound,
            &mpi_request[nsend]);
#ifndef DEBUG_CHECK_ASSERTIONS
      NULL_USE(ierr);
#endif
      TBOX_ASSERT(ierr == MPI_SUCCESS);
   }
   d_common->t_share_new_relationships_send->stop();

   {
      /*
       * The rest of this method assumes current process is NOT
       * in relationship_senders, so remove it.  For efficiency, method
       * computeNewNeighborhoodSets() (which created the relationship senders)
       * did not remove it.
       */
      IntSet::iterator local = relationship_senders.find(d_common->rank);
      if (local != relationship_senders.end()) {
         relationship_senders.erase(local);
      }
   }

   /*
    * Create set recved_from which is to contain ranks of
    * processes from which we've received the expected relationship data.
    * The while loop goes until all expected messages have
    * been received from relationship_senders.
    *
    * In the while loop:
    *    - Probe for an incomming message.
    *    - Determine its size allocate memory for receiving the message.
    *    - Receive the message.
    *    - Get relationship data from the message.
    */
   IntSet recved_from;
   while (recved_from.size() < relationship_senders.size()) {

      d_common->t_share_new_relationships_recv->start();
      ierr = mpi.Probe(MPI_ANY_SOURCE,
            d_common->tag_upper_bound,
            &mpi_status);
      TBOX_ASSERT(ierr == MPI_SUCCESS);

      const int sender = mpi_status.MPI_SOURCE;
      int mesg_size = -1;
      mpi.Get_count(&mpi_status, MPI_INT, &mesg_size);
      TBOX_ASSERT(relationship_senders.find(sender) != relationship_senders.end());
      TBOX_ASSERT(recved_from.find(sender) == recved_from.end());
      TBOX_ASSERT(mesg_size >= 0);

      tbox::Array<int> buf(mesg_size);
      int* ptr = buf.getPointer();
      ierr = mpi.Recv(ptr,
            mesg_size,
            MPI_INT,
            sender,
            d_common->tag_upper_bound,
            &mpi_status);
      TBOX_ASSERT(ierr == MPI_SUCCESS);
      d_common->t_share_new_relationships_recv->stop();

      d_common->t_share_new_relationships_unpack->start();
      int consumed = 0;
      while (ptr < buf.getPointer() + buf.size()) {
         const hier::LocalId new_local_id(*(ptr++));
         hier::BoxId box_id(new_local_id, d_common->rank);
         int n_new_relationships = *(ptr++);
         TBOX_ASSERT(d_common->new_to_tag->hasNeighborSet(box_id));
         if (n_new_relationships > 0) {
            hier::Connector::NeighborhoodIterator base_box_itr =
               d_common->new_to_tag->makeEmptyLocalNeighborhood(box_id);
            for (int n = 0; n < n_new_relationships; ++n) {
               hier::Box node(d_dim);
               node.getFromIntBuffer(ptr);
               ptr += ints_per_node;
               d_common->new_to_tag->insertLocalNeighbor(node, base_box_itr);
            }
         }
         consumed += 2 + n_new_relationships * ints_per_node;
      }
      recved_from.insert(sender);
      d_common->t_share_new_relationships_unpack->stop();
   }

   if (nsend > 0) {
      // Make sure all nonblocking sends completed.
      d_common->t_share_new_relationships_send->start();
      tbox::Array<tbox::SAMRAI_MPI::Status> mpi_statuses(
         static_cast<int>(relationship_messages.size()));
      ierr = mpi.Waitall(static_cast<int>(relationship_messages.size()),
            mpi_request.getPointer(),
            mpi_statuses.getPointer());
      TBOX_ASSERT(ierr == MPI_SUCCESS);
      d_common->t_share_new_relationships_send->stop();
   }

   d_common->t_share_new_relationships->stop();

}

/*
 ********************************************************************
 * Utility methods.
 ********************************************************************
 */

int*
BergerRigoutsosNode::putHistogramToBuffer(
      int* buffer)
{
  int dim_val = d_dim.getValue();
   for (int d = 0; d < dim_val; ++d) {
      d_histogram[d].resize(d_box.numberCells(d), BAD_INTEGER);
      memcpy(buffer,
         &d_histogram[d][0],
         d_box.numberCells(d) * sizeof(int));
      buffer += d_box.numberCells(d);
   }
   return buffer;
}

int*
BergerRigoutsosNode::getHistogramFromBuffer(
      int* buffer)
{
   unsigned int dim_val = d_dim.getValue();
   for (unsigned int d = 0; d < dim_val; ++d) {
      TBOX_ASSERT((int)d_histogram[d].size() == d_box.numberCells(d));
      // d_histogram[d].resizeArray( d_box.numberCells(d) );
      memcpy(&d_histogram[d][0],
         buffer,
         d_box.numberCells(d) * sizeof(int));
      buffer += d_box.numberCells(d);
   }
   return buffer;
}

int*
BergerRigoutsosNode:: putBoxToBuffer(
      const hier::Box& box,
      int* buffer) const
{
   const hier::IntVector& l = box.lower();
   const hier::IntVector& u = box.upper();
   int dim_val = d_dim.getValue();
   for (int d = 0; d < dim_val; ++d) {
      *(buffer++) = l(d);
      *(buffer++) = u(d);
   }
   return buffer;
}

int*
BergerRigoutsosNode::getBoxFromBuffer(
   hier::Box& box,
   int* buffer) const
{
   hier::IntVector& l = box.lower();
   hier::IntVector& u = box.upper();
   int dim_val = d_dim.getValue();
   for (int d = 0; d < dim_val; ++d) {
      l(d) = *(buffer++);
      u(d) = *(buffer++);
   }
   return buffer;
}

/*
 ***********************************************************************
 * Put in dropouts things that are in main_group but
 * not in sub_group.
 *
 * Assume that sub_group is a subset of elements in main_group.
 * Assume that sub_group and main_group are sorted in ascending order.
 *
 * Assume add_root is NOT in the dropout and add it anyway.
 ***********************************************************************
 */
void
BergerRigoutsosNode::computeDropoutGroup(
   const VectorOfInts& main_group,
   const VectorOfInts& sub_group,
   VectorOfInts& dropout_group,
   int add_root) const
{
   TBOX_ASSERT(main_group.size() >= sub_group.size());

   dropout_group.resize(main_group.size(), BAD_INTEGER);

   size_t i, j, k = 0;
   dropout_group[k++] = add_root;
   for (i = 0, j = 0; i < main_group.size(); ++i) {
      if (main_group[i] != sub_group[j]) {
         dropout_group[k++] = main_group[i];
      } else {
         ++j;
         if (j == sub_group.size()) {
            // No more in the sub_group so the rest of main_group
            // goes in dropout_group.
            for (i = i + 1; i < main_group.size(); ++i, ++k) {
               dropout_group[k] = main_group[i];
            }
         }
      }
   }

   TBOX_ASSERT(j = sub_group.size());
   dropout_group.resize(k, BAD_INTEGER);
}

/*
 **********************************************************************
 * Claim a unique tag from the processor's available tag pool.
 * Check that the pool is not overused.
 **********************************************************************
 */
void
BergerRigoutsosNode::claimMPITag()
{
   /*
    * Each dendogram node should claim no more than one MPI tag
    * so make sure it does not already have one.
    */
   TBOX_ASSERT(d_mpi_tag < 0);

   d_mpi_tag = d_common->available_mpi_tag;
   d_common->available_mpi_tag = d_mpi_tag + total_phase_tags;
   if (d_mpi_tag + total_phase_tags - 1 >
       d_common->tag_upper_bound / (d_common->nproc) * (d_common->rank + 1)) {
      /*
       * Each process is alloted tag_upper_bound/(d_common->nproc)
       * tag values.  If it needs more than this, it will encroach
       * on the tag pool of the next process and may lead to using
       * non-unique tags.
       */
      TBOX_ERROR("Out of MPI tag values need to ensure that\n"
         << "messages are properly differentiated."
         << "\nd_mpi_tag = " << d_mpi_tag
         << "\ntag_upper_bound = " << d_common->tag_upper_bound
         << "\nmber of nodes = " << d_common->nproc
         << "\nmax tag required = " << d_mpi_tag + total_phase_tags - 1
         << "\nmax tag available = "
         << d_common->tag_upper_bound / (d_common->nproc) * (d_common->rank + 1));
      /*
       * It is probably safe to recycle tags if we run out of MPI tags.
       * This is not implemented because thus far, there is no need for it.
       * Recycling is starting over from the initial tag set aside for the
       * local process.  To make sure that recycled tags are not still
       * in use, we should claim a new (or recycled) tag for the dropout
       * broadcast phase.  This is because descendant nodes may recycle
       * the current claimed tag before this phase starts.  All other
       * phases are not interupted by descendant communications, so we
       * are assured that their tag is not doubly claimed.
       */
   }
}

/*
 **********************************************************************
 * Convert an integer value to BoxAcceptance.
 * This is needed because the compiler cannot
 * cast an integer to an enum type.
 **********************************************************************
 */
BergerRigoutsosNode::BoxAcceptance
BergerRigoutsosNode::intToBoxAcceptance(
   int i) const
{
   switch (i) {
      case undetermined: return undetermined;

      case hasnotag_by_owner: return hasnotag_by_owner;

      case rejected_by_calculation: return rejected_by_calculation;

      case accepted_by_calculation: return accepted_by_calculation;

      case rejected_by_owner: return rejected_by_owner;

      case accepted_by_owner: return accepted_by_owner;

      case rejected_by_recombination: return rejected_by_recombination;

      case accepted_by_recombination: return accepted_by_recombination;

      case rejected_by_dropout_bcast: return rejected_by_dropout_bcast;

      case accepted_by_dropout_bcast: return accepted_by_dropout_bcast;

      default:
         TBOX_ERROR("Library error: bad BoxAcceptance data of " << i << ".\n");
   }
   return undetermined;
}

void
BergerRigoutsosNode::printClassData(
   std::ostream& os,
   int detail_level) const
{
   os << "ID              " << d_pos << " owner=" << d_owner << " box="
      << d_box
   ;
   if (detail_level > 0) {
      os << "\nfamily          " << (d_parent == NULL ? 0 : d_parent->d_pos)
         << ' ' << (d_lft_child ? (d_lft_child->d_pos) : -1)
         << ' ' << (d_rht_child ? (d_rht_child->d_pos) : -1)
      ;
   }
   if (detail_level > 1) {
      os << "\nthis           " << this
         << "\ngeneration     " << d_generation << " place="
         << ((d_pos % 2) ? 'r' : 'l')
         << "\nnode           " << d_mapped_box
         << "\nbox_acceptance " << d_box_acceptance
         << "\noverlap        " << d_overlap
         << "\ngroup          " << d_group.size() << ':'
      ;
      for (size_t i = 0; i < d_group.size(); ++i) {
         os << ' ' << d_group[i];
      }
      os << std::endl;
   }
}

/*
 **********************************************************************
 *
 * Methods for setting algorithm parameters before running.
 *
 **********************************************************************
 */

void
BergerRigoutsosNode::setAlgorithmAdvanceMode(
   const std::string& mode)
{
   if (mode == "ADVANCE_ANY") {
      d_common->algo_advance_mode = ADVANCE_ANY;
   } else if (mode == "ADVANCE_SOME") {
      d_common->algo_advance_mode = ADVANCE_SOME;
   } else if (mode == "SYNCHRONOUS") {
      d_common->algo_advance_mode = SYNCHRONOUS;
   } else {
      TBOX_ERROR("No such algorithm choice: " << mode << "\n");
   }
}

void
BergerRigoutsosNode::setOwnerMode(
   const std::string& mode)
{
   if (mode == "SINGLE_OWNER") {
      d_common->owner_mode = SINGLE_OWNER;
   } else if (mode == "MOST_OVERLAP") {
      d_common->owner_mode = MOST_OVERLAP;
   } else if (mode == "FEWEST_OWNED") {
      d_common->owner_mode = FEWEST_OWNED;
   } else if (mode == "LEAST_ACTIVE") {
      d_common->owner_mode = LEAST_ACTIVE;
   } else {
      TBOX_ERROR("BergerRigoutsosNode: Unrecognized owner mode request: "
         << mode);
   }
}

void
BergerRigoutsosNode::setComputeRelationships(
   const std::string mode,
   const hier::IntVector& ghost_cell_width)
{
   if (mode == "NONE") {
      d_common->compute_relationships = 0;
   } else if (mode == "TAG_TO_NEW") {
      d_common->compute_relationships = 1;
   } else if (mode == "BIDIRECTIONAL") {
      d_common->compute_relationships = 2;
   } else {
      TBOX_ERROR("BergerRigoutsosNode::setComputeRelationships error:\n"
         << "bad mode '" << mode << "' specified.\n"
         << "Should be one of NONE, TAG_TO_NEW, BIDIRECTIONAL");
   }
   TBOX_ASSERT(ghost_cell_width >= hier::IntVector::getZero(d_dim));
   d_common->max_gcw = ghost_cell_width;
}

/*
 **********************************************************************
 *
 * Methods for collecting data on dendogram after it completes.
 *
 **********************************************************************
 */

void
BergerRigoutsosNode::printDendogramState(
   std::ostream& co,
   const std::string& border) const
{
   co << border;
   for (int i = 0; i < d_generation; ++i) {
      co << "  ";
   }
   printState(co);
   co << std::endl;
   if (d_lft_child) {
      d_lft_child->printDendogramState(co, border);
   }
   if (d_rht_child) {
      d_rht_child->printDendogramState(co, border);
   }
}

void
BergerRigoutsosNode::printState(
   std::ostream& co) const
{
   co << d_generation << ':' << d_pos << '=' << d_mapped_box
      << "  o=" << d_owner << ',' << (d_common->rank == d_owner)
      << "  a=" << d_box_acceptance
      << "  w=" << d_wait_phase << '/' << bool(d_comm_group)
      << (d_comm_group ? d_comm_group->isDone() : true)
      << "  t=" << d_num_tags;
   if (d_lft_child) {
      co << "  l=" << d_lft_child->d_generation << ':' << d_lft_child->d_pos
         << '=' << d_lft_child->d_mapped_box;
   }
   if (d_rht_child) {
      co << "  r=" << d_rht_child->d_generation << ':' << d_rht_child->d_pos
         << '=' << d_rht_child->d_mapped_box;
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
void
BergerRigoutsosNode::initializeCallback()
{
   // Timers
   CommonParams::t_cluster = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::cluster");
   CommonParams::t_cluster_and_compute_relationships = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::clusterAndComputeRelationships()");
   CommonParams::t_continue_algorithm = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::continueAlgorithm()");

   CommonParams::t_compute = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::compute");
   CommonParams::t_comm_wait = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::Comm_wait");
   CommonParams::t_MPI_wait = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::MPI_wait");

   CommonParams::t_compute_new_graph_relationships = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::computeNewNeighborhoodSets()");
   CommonParams::t_share_new_relationships = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::shareNewNeighborhoodSetsWithOwners()");
   CommonParams::t_share_new_relationships_send = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::shareNewNeighborhoodSetsWithOwners()_send");
   CommonParams::t_share_new_relationships_recv = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::shareNewNeighborhoodSetsWithOwners()_recv");
   CommonParams::t_share_new_relationships_unpack = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::shareNewNeighborhoodSetsWithOwners()_unpack");

   CommonParams::t_local_histogram = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::makeLocalTagHistogram()");
   CommonParams::t_local_tasks = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::continueAlgorithm()_local_tasks");

   // Multi-stage timers.
   CommonParams::t_reduce_histogram = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::reduce_histogram");
   CommonParams::t_bcast_acceptability = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::bcast_acceptability");
   CommonParams::t_gather_grouping_criteria = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::gather_grouping_criteria");
   CommonParams::t_bcast_child_groups = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::bcast_child_groups");
   CommonParams::t_bcast_to_dropouts = tbox::TimerManager::getManager()->
      getTimer("mesh::BergerRigoutsosNode::bcast_to_dropouts");

}

/*
 *************************************************************************
 *************************************************************************
 */
void
BergerRigoutsosNode::finalizeCallback()
{
   CommonParams::t_cluster.reset();
   CommonParams::t_cluster_and_compute_relationships.reset();
   CommonParams::t_continue_algorithm.reset();
   CommonParams::t_compute.reset();
   CommonParams::t_comm_wait.reset();
   CommonParams::t_MPI_wait.reset();
   CommonParams::t_compute_new_graph_relationships.reset();
   CommonParams::t_share_new_relationships.reset();
   CommonParams::t_share_new_relationships_send.reset();
   CommonParams::t_share_new_relationships_recv.reset();
   CommonParams::t_share_new_relationships_unpack.reset();
   CommonParams::t_local_tasks.reset();
   CommonParams::t_local_histogram.reset();
   CommonParams::t_reduce_histogram.reset();
   CommonParams::t_bcast_acceptability.reset();
   CommonParams::t_gather_grouping_criteria.reset();
   CommonParams::t_bcast_child_groups.reset();
   CommonParams::t_bcast_to_dropouts.reset();
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
