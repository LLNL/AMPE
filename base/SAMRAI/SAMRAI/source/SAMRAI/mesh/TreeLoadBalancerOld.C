/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_TreeLoadBalancerOld_C
#define included_mesh_TreeLoadBalancerOld_C

#include "SAMRAI/mesh/TreeLoadBalancerOld.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"

#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/AsyncCommGroup.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/Statistician.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <cmath>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

const int TreeLoadBalancerOld::TreeLoadBalancerOld_LOADTAG0;
const int TreeLoadBalancerOld::TreeLoadBalancerOld_LOADTAG1;
const int TreeLoadBalancerOld::TreeLoadBalancerOld_EDGETAG0;
const int TreeLoadBalancerOld::TreeLoadBalancerOld_EDGETAG1;
const int TreeLoadBalancerOld::TreeLoadBalancerOld_PREBALANCE0;
const int TreeLoadBalancerOld::TreeLoadBalancerOld_PREBALANCE1;
const int TreeLoadBalancerOld::TreeLoadBalancerOld_FIRSTDATALEN;
const int TreeLoadBalancerOld::TreeLoadBalancerOld_MIN_NPROC_FOR_AUTOMATIC_MULTICYCLE;

const int TreeLoadBalancerOld::d_default_data_id = -1;

/*
 *************************************************************************
 * TreeLoadBalancerOld constructor.
 *************************************************************************
 */

TreeLoadBalancerOld::TreeLoadBalancerOld(
   const tbox::Dimension& dim,
   const std::string& name,
   const boost::shared_ptr<tbox::Database>& input_db):
   d_dim(dim),
   d_object_name(name),
   d_mpi(tbox::SAMRAI_MPI::commNull),
   d_mpi_is_dupe(false),
   d_n_root_cycles(-1),
   d_degree(2),
   d_master_workload_data_id(d_default_data_id),
   d_balance_penalty_wt(1.0),
   d_surface_penalty_wt(1.0),
   d_slender_penalty_wt(1.0),
   d_slender_penalty_threshold(3.0),
   d_precut_penalty_wt(1.0),
   // Data shared during balancing.
   d_min_size(d_dim),
   d_max_size(d_dim),
   d_bad_interval(d_dim),
   d_cut_factor(d_dim),
   // Output control.
   d_report_load_balance(false),
   // Performance evaluation.
   d_barrier_before(false),
   d_barrier_after(false),
   d_print_steps(false),
   d_print_break_steps(false),
   d_print_swap_steps(false),
   d_print_edge_steps(false),
   d_check_connectivity(false),
   d_check_map(false)
{
   TBOX_ASSERT(!name.empty());
   getFromInput(input_db);
   setTimers();
}



/*
 *************************************************************************
 * TreeLoadBalancerOld constructor.
 *************************************************************************
 */

TreeLoadBalancerOld::~TreeLoadBalancerOld()
{
   freeMPICommunicator();
}



/*
 *************************************************************************
 * Accessory functions to get/set load balancing parameters.
 *************************************************************************
 */

bool
TreeLoadBalancerOld::getLoadBalanceDependsOnPatchData(
   int level_number) const
{
   return getWorkloadDataId(level_number) < 0 ? false : true;
}



/*
**************************************************************************
**************************************************************************
*/
void
TreeLoadBalancerOld::setWorkloadPatchDataIndex(
   int data_id,
   int level_number)
{
   boost::shared_ptr<pdat::CellDataFactory<double> > datafact(
      hier::VariableDatabase::getDatabase()->getPatchDescriptor()->
      getPatchDataFactory(data_id),
      boost::detail::dynamic_cast_tag());
   if (!datafact) {
      TBOX_ERROR(
         d_object_name << " error: "
                       << "\n   data_id " << data_id << " passed to "
                       << "setWorkloadPatchDataIndex()"
                       << " does not refer to cell-centered double patch data. " << std::endl);
   }

   if (level_number >= 0) {
      int asize = d_workload_data_id.getSize();
      if (asize < level_number + 1) {
         d_workload_data_id.resizeArray(level_number + 1);
         for (int i = asize; i < level_number - 1; i++) {
            d_workload_data_id[i] =
               d_master_workload_data_id;
         }
         d_workload_data_id[level_number] = data_id;
      }
   } else {
      d_master_workload_data_id = data_id;
      for (int ln = 0; ln < d_workload_data_id.getSize(); ln++) {
         d_workload_data_id[ln] = d_master_workload_data_id;
      }
   }
}



/*
 *************************************************************************
 * This method implements the abstract LoadBalanceStrategy interface,
 * but it is not where the tree load balancer algorithm is implemented.
 *
 * This method does some preliminary setup then calls
 * loadBalanceWithinRankGroup to compute the new balanced
 * BoxLevel and the mapping Connectors between the old and the new.
 * Then it applies the mapping to update the balance<==>anchor
 * Connectors.  It may do this multiple times, as specified by the
 * cycling parameter.
 *
 * After load balancing, it enforces the maximum size restriction
 * by breaking up large boxes and update balance<==>anchor again.
 *************************************************************************
 */
void
TreeLoadBalancerOld::loadBalanceBoxLevel(
   hier::BoxLevel& balance_box_level,
   hier::Connector& balance_to_anchor,
   hier::Connector& anchor_to_balance,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const hier::Connector& balance_to_attractor,
   const hier::Connector& attractor_to_balance,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::BoxLevel& domain_box_level,
   const hier::IntVector& bad_interval,
   const hier::IntVector& cut_factor,
   const tbox::RankGroup& rank_group) const
{
   NULL_USE(balance_to_attractor);
   NULL_USE(attractor_to_balance);
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   TBOX_ASSERT(anchor_to_balance.isFinalized() ==
      balance_to_anchor.isFinalized());
   if (anchor_to_balance.isFinalized()) {
      TBOX_ASSERT(anchor_to_balance.isTransposeOf(balance_to_anchor));
   }
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS6(d_dim,
      balance_box_level,
      min_size,
      max_size,
      domain_box_level,
      bad_interval,
      cut_factor);
   if (hierarchy) {
      TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, *hierarchy);
   }


   if ( d_mpi_is_dupe ) {
      /*
       * If user has set the duplicate communicator, make sure it is
       * compatible with the BoxLevel involved.
       */
      TBOX_ASSERT(d_mpi.getSize() == balance_box_level.getMPI().getSize());
      TBOX_ASSERT(d_mpi.getRank() == balance_box_level.getMPI().getRank());
#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_mpi.getSize() > 1) {
         int compare_result;
         tbox::SAMRAI_MPI::Comm_compare(
            d_mpi.getCommunicator(),
            balance_box_level.getMPI().getCommunicator(),
            &compare_result);
         if (compare_result != MPI_CONGRUENT) {
            TBOX_ERROR("TreeLoadBalancerOld::loadBalanceBoxLevel:\n"
               << "The input balance_box_level has a SAMRAI_MPI that is\n"
               << "not congruent with the one set with setSAMRAI_MPI().\n"
               << "You must use freeMPICommunicator() before balancing\n"
               << "a BoxLevel with an incongruent SAMRAI_MPI.");
         }
      }
#endif
   }
   else {
      d_mpi = balance_box_level.getMPI();
   }

   if (d_print_steps ||
       d_print_break_steps) {
      tbox::plog << "TreeLoadBalancerOld::loadBalanceBoxLevel called with:"
                 << "\n  min_size = " << min_size
                 << "\n  max_size = " << max_size
                 << "\n  bad_interval = " << bad_interval
                 << "\n  cut_factor = " << cut_factor
                 << std::endl << balance_box_level.format("", 2);
   }


   /*
    * Periodic image Box should be ignored during load balancing
    * because they have no real work.  The load-balanced results
    * should contain no periodic images.
    *
    * To avoid need for special logic to skip periodic images while
    * load balancing, we just remove periodic images in the
    * balance_box_level and all periodic edges in
    * anchor<==>balance.
    */

   balance_box_level.removePeriodicImageBoxes();
   if (balance_to_anchor.isFinalized()) {

      anchor_to_balance.removePeriodicRelationships();
      anchor_to_balance.setHead(balance_box_level, true);

      balance_to_anchor.removePeriodicRelationships();
      balance_to_anchor.setBase(balance_box_level, true);

   }


   if (d_barrier_before) {
      t_barrier_before->start();
      d_mpi.Barrier();
      t_barrier_before->stop();
   }

   if (!rank_group.containsAllRanks()) {
      prebalanceBoxLevel(balance_box_level,
         balance_to_anchor,
         anchor_to_balance,
         rank_group);
   }

   t_load_balance_box_level->start();

   d_min_size = min_size;
   d_max_size = max_size;
   d_bad_interval = bad_interval;
   d_cut_factor = cut_factor;
   /*
    * Domain boxes are used by breakOffLoad to determine where
    * the bad cuts are.  Computing domain_boxes from domain_box_level
    * should be moved above the this method.
    */

   /*
    * We expect the domain box_level to be in globalized state.
    */
   TBOX_ASSERT(
      domain_box_level.getParallelState() ==
      hier::BoxLevel::GLOBALIZED);

   d_block_domain_boxes.clear();
   int nblocks =
      domain_box_level.getGridGeometry()->getNumberBlocks();
   d_block_domain_boxes.resize(nblocks);

   if (nblocks == 1) {
      domain_box_level.getGlobalBoxes(d_block_domain_boxes[0]);
      d_block_domain_boxes[0].refine(balance_box_level.getRefinementRatio());
   } else {
      for (int b = 0; b < nblocks; ++b) {
         d_block_domain_boxes[b] = hier::BoxContainer(
            domain_box_level.getGlobalBoxes(), hier::BlockId(b));

         d_block_domain_boxes[b].refine(balance_box_level.getRefinementRatio());
      }
   }

   if (d_print_steps) {
      tbox::plog << "Pre balanced:\n" << balance_box_level.format("", 2);
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   if (balance_to_attractor.isFinalized()) {
      /*
       * If balance_to_attractor is given, sanity-check it.
       */
      if (&balance_box_level != &balance_to_attractor.getBase() &&
          !(balance_box_level == balance_to_attractor.getBase())) {
         TBOX_ERROR(
            "TreeLoadBalancerOld::loadBalanceBoxLevel: balance_box_level\n"
            << "does not match the base of balance_to_attractor.");
      }
      if (!balance_to_attractor.isTransposeOf(attractor_to_balance)) {
         TBOX_ERROR("TreeLoadBalancerOld::loadBalanceBoxLevel:\n"
            << "attractor_to_balance and balance_to_attractor\n"
            << "are not transposes of each other.");
      }
   }
#endif


   t_compute_local_load->start();
   double local_load = computeLocalLoads(balance_box_level);
   t_compute_local_load->stop();

   size_t nproc_with_initial_load =
      balance_box_level.getLocalNumberOfBoxes() > 0;

   double global_sum_load;

   {
      /*
       * Determine the total load and number of processes that has any
       * initial load.
       *
       * TODO: If there's more than one rank group, shouldn't this a
       * global reduction for each rank group instead of a single one
       * for all?
       */
      t_compute_global_load->start();
      if (d_mpi.getSize() > 1) {
         double dtmp[2], dtmp_sum[2];
         dtmp[0] = local_load;
         dtmp[1] = static_cast<double>(nproc_with_initial_load);
         d_mpi.Allreduce(dtmp,
            dtmp_sum,
            2,
            MPI_DOUBLE,
            MPI_SUM);
         global_sum_load = dtmp_sum[0];
         nproc_with_initial_load = (size_t)dtmp_sum[1];
      } else {
         global_sum_load = local_load;
      }
      t_compute_global_load->stop();
      if (d_print_steps) {
         tbox::plog << "TreeLoadBalancerOld::loadBalanceBoxLevel balancing "
                    << global_sum_load << " (initially born on "
                    << nproc_with_initial_load << " procs) across all "
                    << d_mpi.getSize()
                    << " procs, averaging " << global_sum_load / d_mpi.getSize()
                    << " or " << pow(global_sum_load / d_mpi.getSize(), 1.0 / d_dim.getValue())
                    << "^" << d_dim << " per proc." << std::endl;
      }
   }


   d_global_avg_load = global_sum_load / rank_group.size();


   /*
    * User can set the number of cycles to use (see n_root_cycles
    * input parameter), or leave it negative to let this class set it
    * automatically using the following heuristic algorithm:
    *
    * If machine size is small enough to have negligible scaling
    * issues (< TreeLoadBalancerOld_MIN_NPROC_FOR_AUTOMATIC_MULTICYCLE), use 1 cycle.
    *
    * Else if the initial load is sufficiently spread out (across at
    * least sqrt(nproc)) processes, use 1 cycle.
    *
    * Else use 2 cycles.
    */
   int number_of_cycles = d_n_root_cycles;
   if (number_of_cycles < 0) {
      // User requested automatic number of cycles.
      if (balance_box_level.getMPI().getSize() < TreeLoadBalancerOld_MIN_NPROC_FOR_AUTOMATIC_MULTICYCLE ) {
         number_of_cycles = 1;
      }
      else if ( int(nproc_with_initial_load * nproc_with_initial_load) >=
                balance_box_level.getMPI().getSize() ) {
         number_of_cycles = 1;
      } else {
         number_of_cycles = 2;
      }
   }



   /*
    * The icycle loop spreads out the work each time through.  If
    * using more than one cycle, only the last one tries to balance
    * across all of d_mpi.
    */

   for (int icycle = 0; icycle < number_of_cycles; ++icycle) {

      // If not the first cycle, local_load needs updating.
      if (icycle > 0) {
         t_compute_local_load->start();
         local_load = computeLocalLoads(balance_box_level);
         t_compute_local_load->stop();
      }

      if (d_report_load_balance) {
         // Debugging: check overall load balance at intermediate cycles.
         tbox::plog
         << "TreeLoadBalancerOld::loadBalanceBoxLevel results before cycle "
         << icycle << ":" << std::endl;
         BalanceUtilities::gatherAndReportLoadBalance(
            local_load,
            balance_box_level.getMPI());
      }


      const bool last_cycle = (icycle == number_of_cycles-1);

      /*
       * Determine whether to use rank_group as is or subgroup it based
       * on cycles.
       */

      const tbox::RankGroup *rank_group_in_use = &rank_group;
      int number_of_groups = 1;
      int group_num = 0;

      tbox::RankGroup cycle_rank_group(d_mpi);
      if ( !last_cycle && rank_group.containsAllRanks() ) {
         createBalanceRankGroupBasedOnCycles(
            cycle_rank_group,
            number_of_groups,
            group_num,
            icycle,
            number_of_cycles);
         rank_group_in_use = &cycle_rank_group;
      }


      /*
       * Compute the load for the group.  If this is the last cycle,
       * the group must include all processes, and the group's load
       * is the global sum load.  Else, use all-reduce to get the
       * group load.
       */
      t_compute_tree_load->start();

      double group_sum_load;

      if (icycle == number_of_cycles - 1) {

         group_sum_load = global_sum_load;

      } else {

         t_compute_tree_load_for_cycle[icycle]->start();

         /*
          * Use MPI's vector all-reduce to get individual group loads.
          * This gives more info than the process needs, but because the
          * number of groups << number of procs, it is still faster
          * (probably) than hand coded conmunication.
          */
         std::vector<double> group_loads(number_of_groups, 0.0);
         group_loads[group_num] = local_load;
         if (d_mpi.getSize() > 1) {
            d_mpi.AllReduce(&group_loads[0],
                            static_cast<int>(group_loads.size()),
                            MPI_SUM);
         }
         group_sum_load = group_loads[group_num];

         t_compute_tree_load_for_cycle[icycle]->stop();

      }

      t_compute_tree_load->stop();


      /*
       * Compute the load-balancing map.
       */

      loadBalanceWithinRankGroup(
         balance_box_level,
         balance_to_anchor,
         anchor_to_balance,
         rank_group,
         group_sum_load );

      if (d_barrier_after) {
         t_barrier_after->start();
         d_mpi.Barrier();
         t_barrier_after->stop();
      }

   }


   /*
    * If max_size is given (positive), constrain boxes to the given
    * max_size.  If not given, skip the enforcement step to save some
    * communications.
    */

   if (max_size > hier::IntVector::getZero(d_dim)) {

      t_constrain_size->start();
      constrainMaxBoxSizes(
         balance_box_level,
         anchor_to_balance,
         balance_to_anchor );
      t_constrain_size->stop();

      if (d_print_steps) {
         tbox::plog << " TreeLoadBalancerOld completed constraining box sizes."
                    << "\n";
      }

   }


   /*
    * Finished load balancing.  Clean up and wrap up.
    */

   d_min_size = hier::IntVector(d_dim, -1);
   d_max_size = hier::IntVector(d_dim, -1);
   d_block_domain_boxes.clear();
   d_bad_interval = hier::IntVector(d_dim, -1);
   d_cut_factor = hier::IntVector(d_dim, -1);

   t_load_balance_box_level->stop();

   local_load = computeLocalLoads(balance_box_level);
   d_load_stat.push_back(local_load);
   d_box_count_stat.push_back(
      static_cast<int>(balance_box_level.getBoxes().size()));

   if (d_report_load_balance) {
      t_report_loads->start();
      tbox::plog
      << "TreeLoadBalancerOld::loadBalanceBoxLevel results after "
      << number_of_cycles << " cycles:" << std::endl;
      BalanceUtilities::gatherAndReportLoadBalance(local_load,
         balance_box_level.getMPI());
      t_report_loads->stop();
   }

   if (d_check_connectivity && anchor_to_balance.isFinalized()) {
      hier::OverlapConnectorAlgorithm oca;
      tbox::plog << "TreeLoadBalancerOld checking balance-anchor connectivity."
                 << std::endl;
      int errs = 0;
      if (oca.checkOverlapCorrectness(anchor_to_balance, false, true, true)) {
         ++errs;
         tbox::perr << "Error found in anchor_to_balance!\n";
      }
      if (oca.checkOverlapCorrectness(balance_to_anchor, false, true, true)) {
         ++errs;
         tbox::perr << "Error found in balance_to_anchor!\n";
      }
      if (anchor_to_balance.checkTransposeCorrectness(
             balance_to_anchor)) {
         ++errs;
         tbox::perr << "Error found in balance-anchor transpose!\n";
      }
      if (errs != 0) {
         TBOX_ERROR(
            "Errors in load balance mapping found."
            << "anchor_box_level:\n" << anchor_to_balance.getBase().format("", 2)
            << "balance_box_level:\n" << balance_box_level.format("", 2)
            << "anchor_to_balance:\n" << anchor_to_balance.format("", 2)
            << "balance_to_anchor:\n" << balance_to_anchor.format("", 2));
      }
      tbox::plog << "TreeLoadBalancerOld checked balance-anchor connectivity."
                 << std::endl;
   }

   if (d_barrier_after) {
      t_barrier_after->start();
      d_mpi.Barrier();
      t_barrier_after->stop();
   }

}



/*
 *************************************************************************
 * Constrain maximum box sizes in the given BoxLevel and
 * update given Connectors to the changed BoxLevel.
 *************************************************************************
 */
void
TreeLoadBalancerOld::constrainMaxBoxSizes(
   hier::BoxLevel& box_level,
   hier::Connector &anchor_to_level,
   hier::Connector &level_to_anchor ) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, box_level);

   t_map_big_boxes->start();

   if (d_print_break_steps) {
      tbox::plog << "Mapping oversized boxes starting with "
                 << box_level.getBoxes().size() << " boxes."
                 << std::endl;
   }

   const hier::IntVector& zero_vector(hier::IntVector::getZero(d_dim));

   hier::BoxLevel constrained(box_level.getDim());
   hier::Connector unconstrained_to_constrained;

   constrained.initialize(
      box_level.getRefinementRatio(),
      box_level.getGridGeometry(),
      box_level.getMPI());
   unconstrained_to_constrained.clearNeighborhoods();
   unconstrained_to_constrained.setBase(box_level);
   unconstrained_to_constrained.setHead(constrained);
   unconstrained_to_constrained.setWidth(zero_vector, true);

   const hier::BoxContainer& unconstrained_boxes = box_level.getBoxes();

   hier::LocalId next_available_index = box_level.getLastLocalId() + 1;

   for (hier::BoxContainer::const_iterator ni = unconstrained_boxes.begin();
        ni != unconstrained_boxes.end(); ++ni) {

      const hier::Box& box = *ni;

      const hier::IntVector box_size = box.numberCells();

      /*
       * If box already conform to max size constraint, keep it.
       * Else chop it up and keep the parts.
       */

      if (box_size <= d_max_size) {

         if (d_print_break_steps) {
            tbox::plog << "    Not oversized: " << box
                       << box.numberCells() << "\n";
         }
         constrained.addBox(box);

      } else {

         if (d_print_break_steps) {
            tbox::plog << "    Breaking oversized " << box
                       << box.numberCells() << " ->";
         }
         hier::BoxContainer chopped(box);
         hier::BoxUtilities::chopBoxes(
            chopped,
            d_max_size,
            d_min_size,
            d_cut_factor,
            d_bad_interval,
            d_block_domain_boxes[box.getBlockId().getBlockValue()]);
         TBOX_ASSERT( chopped.size() != 0 );

         if (chopped.size() != 1) {

            hier::Connector::NeighborhoodIterator base_box_itr =
               unconstrained_to_constrained.makeEmptyLocalNeighborhood(
                  box.getId());

            for (hier::BoxContainer::iterator li(chopped);
                 li != chopped.end(); ++li) {

               const hier::Box fragment = *li;

               const hier::Box new_box(fragment,
                                       next_available_index++,
                                       d_mpi.getRank());

               if (d_print_break_steps) {
                  tbox::plog << "  " << new_box
                             << new_box.numberCells();
               }

               constrained.addBox(new_box);

               unconstrained_to_constrained.insertLocalNeighbor(
                  new_box,
                  base_box_itr);

            }

            if (d_print_break_steps) {
               tbox::plog << "\n";
            }

         } else {
            TBOX_ASSERT( box.isSpatiallyEqual( chopped.front() ) );
            if (d_print_break_steps) {
               tbox::plog << " Unbreakable!" << "\n";
            }
            constrained.addBox(box);
         }

      }

   }

   unconstrained_to_constrained.setConnectorType(hier::Connector::MAPPING);

   if (d_print_steps) {
      tbox::plog
      << " TreeLoadBalancerOld::constrainMaxBoxSizes completed building unconstrained_to_constrained"
      << "\n";
   }

   if (anchor_to_level.isFinalized()) {
      // Modify anchor<==>level Connectors and swap box_level with constrained.
      const hier::MappingConnectorAlgorithm mca;
      mca.modify(anchor_to_level,
                 level_to_anchor,
                 unconstrained_to_constrained,
                 &box_level,
                 &constrained);
   } else {
      // Swap box_level and constrained without touching anchor<==>level.
      hier::BoxLevel::swap(box_level, constrained);
   }

   t_map_big_boxes->stop();
}



/*
 *************************************************************************
 * Given an "unbalanced" BoxLevel, compute the BoxLevel that is
 * load-balanced and compute the mapping between the unbalanced and
 * balanced BoxLevels.
 *
 * If given a RankGroup with less than all ranks, we treat it as a
 * specific user request to balance only within the RankGroup and just
 * use the RankGroup as is.  Otherwise, we may generate sub-groups
 * based on the cycle number and balance within the generated
 * sub-group.
 *
 * The objective of balancing over multiple cycles is to avoid
 * unscalable performance in the cases where just a few processes own
 * all the initial load.  By slowly spreading out the load, no process
 * has to set up Connector unbalanced_to_balanced with number of
 * relationships that scales with the machine size.
 *
 * If the local process is not a member of the RankGroup, it does not
 * participate in the work and just sets the output objects to be
 * locally empty.
 *************************************************************************
 */
void
TreeLoadBalancerOld::loadBalanceWithinRankGroup(
   hier::BoxLevel& balance_box_level,
   hier::Connector& balance_to_anchor,
   hier::Connector& anchor_to_balance,
   const tbox::RankGroup& rank_group,
   const double group_sum_load ) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim,
      balance_box_level);

   double group_avg_load = group_sum_load / rank_group.size();


   hier::BoxLevel balanced_box_level(balance_box_level.getDim());
   hier::Connector unbalanced_to_balanced;
   hier::Connector balanced_to_unbalanced;

   /*
    * Initialize empty balanced_box_level and mappings.
    */
   balanced_box_level.initialize(
      balance_box_level.getRefinementRatio(),
      balance_box_level.getGridGeometry(),
      balance_box_level.getMPI());
   balanced_to_unbalanced.setConnectorType(hier::Connector::MAPPING);
   balanced_to_unbalanced.clearNeighborhoods();
   balanced_to_unbalanced.setBase(balanced_box_level);
   balanced_to_unbalanced.setHead(balance_box_level);
   balanced_to_unbalanced.setWidth(hier::IntVector::getZero(d_dim), true);
   unbalanced_to_balanced.setConnectorType(hier::Connector::MAPPING);
   unbalanced_to_balanced.clearNeighborhoods();
   unbalanced_to_balanced.setBase(balance_box_level);
   unbalanced_to_balanced.setHead(balanced_box_level);
   unbalanced_to_balanced.setWidth(hier::IntVector::getZero(d_dim), true);


   if ( !rank_group.isMember(d_mpi.getRank()) ) {
      /*
       * The following assert should be guaranteed by an earlier call
       * to prebalanceBoxLevel.  Having boxes without being in the
       * given rank group leads to undefined results.
       */
      TBOX_ASSERT( balance_box_level.getLocalNumberOfBoxes() == 0 );

      if (anchor_to_balance.isFinalized()) {
         const hier::MappingConnectorAlgorithm mca;
         t_use_map->start();
         mca.modify(
            anchor_to_balance,
            balance_to_anchor,
            unbalanced_to_balanced,
            balanced_to_unbalanced,
            &balance_box_level,
            &balanced_box_level);
         t_use_map->stop();
      } else {
         hier::BoxLevel::swap(balance_box_level, balanced_box_level);
      }
      return;
   }


   t_get_map->start();

   t_load_distribution->start();

   /*
    * Before the last cycle, it is possible for the group average load
    * to be below the global average, if the group just happens to have
    * underloaded processors.  However, there is no point in driving the
    * processor loads down below the global average just to have it
    * brought back up by the last cycle.  It just unnecessarily fragments
    * the boxes and costs more to do.  To prevent this, reset the group
    * average to the global average if it is below.
    */
   group_avg_load =
      tbox::MathUtilities<double>::Max(group_avg_load, d_global_avg_load);


   /*
    * Arrange the group ranks in a BalancedDepthFirstTree in order to get
    * the parent/children in the group.
    *
    * By the way, the BalancedDepthFirstTree currently assumes a binary
    * tree, d_degree = 2.
    */
   TBOX_ASSERT(d_degree == 2);
   // FIXME: BalancedDepthFirstTree could use a constructor that takes a RankGroup.
   tbox::BalancedDepthFirstTree bdfs(0,
                                     rank_group.size()-1,
                                     rank_group.getMapIndex(d_mpi.getRank()),
                                     true);
   const int num_children = bdfs.getNumberOfChildren();


   /*
    * Communication objects for sending to/receiving from
    * parent/children: We could combine all of these AsyncCommStages
    * and most of the AsyncCommPeers, but we intentionally keep them
    * separate to aid performance analysis.
    */

   tbox::AsyncCommStage child_send_stage;
   tbox::AsyncCommPeer<int>* child_sends = NULL;
   tbox::AsyncCommStage parent_send_stage;
   tbox::AsyncCommPeer<int>* parent_send = NULL;

   setupAsyncCommObjects(
      child_send_stage,
      child_sends,
      parent_send_stage,
      parent_send,
      rank_group,
      bdfs );
   child_send_stage.setCommunicationWaitTimer(t_child_send_wait);
   parent_send_stage.setCommunicationWaitTimer(t_parent_send_wait);

   tbox::AsyncCommStage child_recv_stage;
   tbox::AsyncCommPeer<int>* child_recvs = NULL;
   tbox::AsyncCommStage parent_recv_stage;
   tbox::AsyncCommPeer<int>* parent_recv = NULL;

   setupAsyncCommObjects(
      child_recv_stage,
      child_recvs,
      parent_recv_stage,
      parent_recv,
      rank_group,
      bdfs );
   child_recv_stage.setCommunicationWaitTimer(t_child_recv_wait);
   parent_recv_stage.setCommunicationWaitTimer(t_parent_send_wait);



   /*
    * Outline of the tree load balancing algorithm as implemented:
    *
    * 1. For each child of the local process:
    * Receive data from subtree rooted at child (number in
    * subtree, excess work, remaining work in subtree, etc.).
    *
    * 2. Compute data for subtree rooted at self by combining
    * local data with children subtree data.
    *
    * 3. If parent exists:
    * Send subtree info (number in subtree, excess work,
    * remaining work in subtree, etc.) to parent.
    *
    * 4. If parent exists and we need more work:
    * Receive additional work from parent.
    *
    * 5. Partition additional work among children and self.
    *
    * 6. For each child:
    * Send additional work (if any).
    */

   /*
    * Step 1:
    *
    * Post receive for data from subtree rooted at children.
    * We have to do a few local setups, but post the receive
    * now to overlap communication.
    */
   t_get_load_from_children->start();
   for (int c = 0; c < num_children; ++c) {
      child_recvs[c].beginRecv();
      if (child_recvs[c].isDone()) {
         child_recvs[c].pushToCompletionQueue();
      }
   }
   t_get_load_from_children->stop();


   /*
    * Step 2, local part:
    *
    * The local process must generate indices for new and imported
    * boxes.  To do it deterministically, no generated index should
    * depend on message arrival order.  To achieve this, we maintain
    * 2+d_degree values in next_available_index: one for the local
    * process, one for the parent and one for each child.  The first
    * index given to a locally generated Box is some index unused by
    * balance_box_level.  The first index given to a Box
    * from the parent is the same value plus 1.  The first index given
    * to a box from child 0 is the same value plus 2.  And so on.
    * Each time a value from next_available_index is used, we
    * increment it by 2+d_degree so that the 2+d_degree available
    * values can never be the same.  Moreover, boxes from a certain
    * source always take indices from its own set, independent of when
    * boxes from other sources arrive.
    */
   std::vector<hier::LocalId> next_available_index(2 + d_degree);
   next_available_index[0] = balance_box_level.getLastLocalId() + 1;

   /*
    * The next line makes next_available_index[0] divisible by 2+d_degree.
    * It is not strictly necessary but makes debugging much easier because
    * we can quickly associate any value with the source of its Box.
    */
   next_available_index[0] +=
      hier::LocalId(2+d_degree) - (next_available_index[0] % (2 + d_degree));

   for (int c = 1; c < d_degree + 2; ++c) {
      next_available_index[c] = next_available_index[0] + c;
   }


   /*
    * Data for storing and transfering subtree info.
    */
   SubtreeLoadData* child_load_data = new SubtreeLoadData[num_children];
   SubtreeLoadData my_load_data;


   /*
    * Compute local proc's Boxes and loads and store in
    * my_load_data.  This will eventually include data for the subtree.
    * We will add the rest of the subtree's work when we receive that
    * data from the children.
    */
   my_load_data.d_num_procs = 1;
   my_load_data.d_total_work = static_cast<int>(computeLocalLoads(balance_box_level));


   /*
    * unassigned is a container of BoxInTransit that has been released by
    * a process and has not yet been assigned to another.  First, put
    * excess local work (if any) in unassigned.  Imported
    * BoxInTransits are placed here before determining whether to keep
    * them or send them to another part of the tree.
    */
   TransitSet unassigned;


   t_local_balancing->start();

   if (my_load_data.d_total_work <= group_avg_load) {

      /*
       * Local process is underloaded, so put all of balance_box_level into
       * the balanced_box_level (and add more later).
       */
      const hier::BoxContainer& unbalanced_boxes =
         balance_box_level.getBoxes();
      for (hier::BoxContainer::const_iterator ni = unbalanced_boxes.begin();
           ni != unbalanced_boxes.end(); ++ni) {
         balanced_box_level.addBox(*ni);
      }

   } else {
      /*
       * Local process is overloaded, so remove excess loads:
       * - sort BoxInTransit by load
       * - reassignLoads (put excess loads in unassigned container) and
       * - put remainder in balanced_box_level.
       *
       * Note: This algorithm would also work if we put all local
       * Boxes into unassigned (instead of just the excess load).  In
       * the end, all remaining unassigned Boxes get assigned to the
       * local process anyway.  In fact, having more Boxes in
       * unassigned may let reassignLoads do a better job in
       * reassigning loads to children and parents, because it would
       * have more choices.  The reason we place only the excess load
       * into unassigned is to help preserve locality, assuming that
       * current local Boxes may have more local neighbors.  However,
       * practical evidence so far suggest that the lost locality is
       * not that bad, at least for explicit solvers.  Transfering all
       * local Boxes into unassigned at the start may affect the
       * performance of this algorithms though, because it bypasses
       * one reassignLoads call but makes the unassigned Box set
       * bigger, which may make other calls to reassignLoads work
       * harder.  Which one is a bigger effect and whether there are
       * any significant effects at all remains to be seen.
       */

      const hier::BoxContainer& unbalanced_boxes =
         balance_box_level.getBoxes();

      int ideal_transfer = int(0.5 + my_load_data.d_total_work - group_avg_load);

      if (d_print_steps) {
         tbox::plog << "  Reassigning initial overload of " << ideal_transfer
                    << " to unassigned.\n";
      }

      TransitSet
      local_loads(unbalanced_boxes.begin(), unbalanced_boxes.end());

      int actual_transfer = reassignLoads(
         local_loads,
         unassigned,
         next_available_index[d_degree],
         ideal_transfer );

      for (TransitSet::const_iterator
           ni = local_loads.begin(); ni != local_loads.end(); ++ni) {
         const BoxInTransit& box_in_transit = *ni;
         balanced_box_level.addBox(box_in_transit.d_box);
         /*
          * Create edges only for box_in_transit that changed.
          */
         if (box_in_transit.d_box.getLocalId() !=
             box_in_transit.d_orig_box.getLocalId()) {
            balanced_to_unbalanced.insertLocalNeighbor(
               box_in_transit.d_orig_box,
               box_in_transit.d_box.getId());
            unbalanced_to_balanced.insertLocalNeighbor(
               box_in_transit.d_box,
               box_in_transit.d_orig_box.getId());
         }
      }

      if (d_print_steps) {
         tbox::plog << "    Unassigning " << unassigned.size()
                    << " boxes (" << actual_transfer << " / "
                    << ideal_transfer << " units):";
         for (TransitSet::const_iterator
              ni = unassigned.begin(); ni != unassigned.end(); ++ni) {
            tbox::plog << "  " << *ni;
         }
         tbox::plog << std::endl;
      }
   }

   t_local_balancing->stop();



   /*
    * Step 2, remote part:
    *
    * Complete load-receive communications with children.
    * Add imported BoxInTransit to unassigned bin.
    */
   t_get_load_from_children->start();
   while ( child_recv_stage.numberOfCompletedMembers() > 0 ||
           child_recv_stage.advanceSome() ) {

      tbox::AsyncCommPeer<int>* child_recv =
         dynamic_cast<tbox::AsyncCommPeer<int> *>(child_recv_stage.popCompletionQueue());

      TBOX_ASSERT(child_recv != NULL);
      TBOX_ASSERT(child_recv >= child_recvs);
      TBOX_ASSERT(child_recv < child_recvs + num_children);

      const int cindex = static_cast<int>(child_recv - child_recvs);

      TBOX_ASSERT(cindex >= 0 && cindex < num_children);

      /*
       * Extract data from the child cindex, storing it in
       * child_load_data[cindex].  If child sent up any excess Box,
       * put them in unassigned.
       */
      int old_size = static_cast<int>(unassigned.size());
      unpackSubtreeLoadData(
         child_load_data[cindex],
         unassigned,
         next_available_index[cindex],
         child_recv->getRecvData(),
         child_recv->getRecvSize() );

      child_load_data[cindex].d_ideal_work =
         int(group_avg_load * child_load_data[cindex].d_num_procs + 0.5);

      if (d_print_steps) {
         TransitSet::const_iterator
            ni = unassigned.begin();
         for (int ii = 0; ii < old_size; ++ii) { ++ni; }
         if (d_print_steps) {
            tbox::plog << "    Got " << unassigned.size() - old_size
                       << " boxes (" << child_load_data[cindex].d_load_imported
                       << " units) from child "
                       << child_recv->getPeerRank() << ":";
            for ( ; ni != unassigned.end(); ++ni) {
               const BoxInTransit& box_in_transit = *ni;
               tbox::plog << "  " << box_in_transit;
            }
            tbox::plog << std::endl;
         }
      }

      // Sum children load into my_load_data.
      my_load_data.d_num_procs += child_load_data[cindex].d_num_procs;
      my_load_data.d_total_work +=
         child_load_data[cindex].d_total_work
         + child_load_data[cindex].d_load_imported;

   }



   // We should have received everything at this point.
   TBOX_ASSERT(!child_recv_stage.hasPendingRequests());

   my_load_data.d_ideal_work = int(group_avg_load * my_load_data.d_num_procs + 0.5);

   if (d_print_steps) {
      tbox::plog << "Received children subtree data." << std::endl;
      for (int c = 0; c < num_children; ++c) {
         tbox::plog << "Child " << child_recvs[c].getPeerRank()
                    << " subtree data: " << child_load_data[c].d_total_work
                    << "/" << child_load_data[c].d_ideal_work
                    << " for " << child_load_data[c].d_num_procs << " procs averaging "
                    << child_load_data[c].d_total_work / child_load_data[c].d_num_procs
                    << " after sending up " << child_load_data[c].d_load_imported
                    << std::endl;
      }
      tbox::plog << "Initial subtree data: " << my_load_data.d_total_work
                 << " / " << my_load_data.d_ideal_work
                 << " for " << my_load_data.d_num_procs << " procs averaging "
                 << my_load_data.d_total_work / my_load_data.d_num_procs
                 << " before sending up anything."
                 << std::endl;
   }

   t_get_load_from_children->stop();



   /*
    * Step 3:
    *
    * Send subtree info and excess work (if any) up to parent.
    */
   t_send_load_to_parent->start();
   if (parent_send != NULL) {

      /*
       * Compute the excess work we want to send to parent.
       * If it's positive, reassign some boxes to the parent.
       */
      int ideal_transfer = my_load_data.d_total_work > my_load_data.d_ideal_work ?
         my_load_data.d_total_work - my_load_data.d_ideal_work : 0;

      if (ideal_transfer > 0) {

         if (d_print_steps) {
            tbox::plog << "  Attempting to reassign " << ideal_transfer
                       << " of unassigned load to parent.\n";
         }

         int actual_transfer = reassignLoads(
            unassigned,
            my_load_data.d_for_export /* to parent */,
            next_available_index[d_degree],
            ideal_transfer );
         my_load_data.d_load_exported = actual_transfer;
         my_load_data.d_total_work -= actual_transfer;

         if (d_print_steps) {
            tbox::plog << "  Giving " << my_load_data.d_for_export.size()
                       << " boxes (" << actual_transfer << " / "
                       << ideal_transfer << " units) to parent "
                       << parent_send->getPeerRank() << ":";
            for (TransitSet::const_iterator
                 ni = my_load_data.d_for_export.begin();
                 ni != my_load_data.d_for_export.end(); ++ni) {
               tbox::plog << "  " << *ni;
            }
            tbox::plog << std::endl;
         }

      }

      /*
       * Send local process's load info, along with any exported work,
       * up to parent.
       */
      std::vector<int> msg;
      packSubtreeLoadData(msg, my_load_data);
      parent_send->beginSend(&msg[0], static_cast<int>(msg.size()));

   }
   t_send_load_to_parent->stop();



   /*
    * Step 4:
    *
    * Finish the send-up.
    * To preclude sending work in both directions, the parent
    * will *not* send a work message down if we sent work up.
    */
   if (parent_send != NULL && my_load_data.d_load_exported == 0) {
      t_parent_load_comm->start();
      t_get_load_from_parent->start();

      parent_recv->beginRecv();

      t_get_load_from_parent->stop();
      t_parent_load_comm->stop();
   }


   /*
    * May do some things here that do not depend on message from
    * parents.  Steve Smith suggested sending work down to underloaded
    * children at this point, even if it makes the local process
    * underloaded as a result.  The local process can recover the
    * correct amount of work when it comes down from the parent.  This
    * would allow some children subtrees to wait less, but it has
    * other consequences.
    */


   if (parent_recv != NULL && my_load_data.d_load_exported == 0) {

      /*
       * Receive and unpack message from parent.  Since we did not
       * export work to parent, parent may import some to us.  Put
       * imported work in unassigned.
       */
      t_get_load_from_parent->start();

      t_parent_load_comm->start();
      parent_recv->completeCurrentOperation();
      t_parent_load_comm->stop();

      int old_size = static_cast<int>(unassigned.size());
      SubtreeLoadData parent_load_data;
      unpackSubtreeLoadData(
         parent_load_data,
         unassigned,
         next_available_index[1 + d_degree],
         parent_recv->getRecvData(),
         parent_recv->getRecvSize() );
      my_load_data.d_load_imported = parent_load_data.d_load_imported;
      my_load_data.d_total_work += parent_load_data.d_load_imported;

      if (d_print_steps) {
         TransitSet::const_iterator
            ni = unassigned.begin();
         for (int i = 0; i < old_size; ++i) {
            ++ni;
         }
         if (d_print_steps) {
            tbox::plog << "    Got " << unassigned.size() - old_size
                       << " boxes (" << parent_load_data.d_load_imported
                       << " units) from parent "
                       << parent_recv->getPeerRank() << ":";
            for ( ; ni != unassigned.end(); ++ni) {
               const BoxInTransit& box_in_transit = *ni;
               tbox::plog << "  " << box_in_transit;
            }
            tbox::plog << std::endl;
         }
      }

      t_get_load_from_parent->stop();
   }

   if (d_print_steps) {
      tbox::plog << "  After parent/children, my total work is "
                 << my_load_data.d_total_work << " / "
                 << my_load_data.d_ideal_work << std::endl;
   }


   /*
    * Step 5 and 6:
    *
    * Reassign unassigned load to children subtrees as needed.
    */

   t_send_load_to_children->start();

   for (int ichild = 0; ichild < num_children; ++ichild) {

      SubtreeLoadData& recip_data = child_load_data[ichild];

      /*
       * Note: To preclude unneeded messages, we do *not* send a work
       * message down if the child sent work up to us.
       */
      if (recip_data.d_load_imported == 0) {

         int ideal_transfer = recip_data.d_ideal_work - recip_data.d_total_work;
         int actual_transfer = 0;

         if (d_print_steps) {
            tbox::plog << "  Attempting to reassign " << ideal_transfer
                       << " of unassigned load to child "
                       << child_sends[ichild].getPeerRank() << "\n";
         }

         if (ideal_transfer > 0) {
            actual_transfer = reassignLoads(
               unassigned,
               recip_data.d_for_export,
               next_available_index[d_degree],
               ideal_transfer );
            recip_data.d_load_exported += actual_transfer;
            recip_data.d_total_work += actual_transfer;
         }

         if (d_print_steps) {
            tbox::plog << "  Giving " << recip_data.d_for_export.size()
                       << " boxes (" << actual_transfer << " / " << ideal_transfer
                       << " units) to child " << ichild << ':'
                       << child_sends[ichild].getPeerRank() << " for "
                       << recip_data.d_num_procs
                       << " procs:";
            for (TransitSet::const_iterator ni = recip_data.d_for_export.begin();
                 ni != recip_data.d_for_export.end(); ++ni) {
               tbox::plog << "  " << *ni;
            }
            tbox::plog << std::endl;
         }

         std::vector<int> msg;
         packSubtreeLoadData(msg, recip_data);
         child_sends[ichild].beginSend(&msg[0], static_cast<int>(msg.size()));

      }

   }

   t_send_load_to_children->stop();


   /*
    * All unassigned boxes should go into balanced_box_level.
    *
    * Put unassigned boxes into balanced_box_level and generate
    * relationships in balanced<==>unbalanced mapping Connectors where
    * required.
    *
    * We remove boxes from unassigned when we no longer need to keep
    * track of them.  We do leave behind boxes originating remotely so
    * we can later notify the originating processes about where the
    * box finally ended up.
    */
   t_local_balancing->start();

   for (TransitSet::iterator
        ni = unassigned.begin();
        ni != unassigned.end(); /* incremented in loop */) {

      const BoxInTransit& box_in_transit = *ni;
      balanced_box_level.addBox(box_in_transit.d_box);

      if (box_in_transit.d_box.isIdEqual(box_in_transit.d_orig_box)) {
         // Unchanged box requires no edge.  Nothing else need to be done.
         unassigned.erase(ni++);
      } else {

         balanced_to_unbalanced.insertLocalNeighbor(
            box_in_transit.d_orig_box,
            box_in_transit.d_box.getId());

         if (box_in_transit.d_orig_box.getOwnerRank() == d_mpi.getRank()) {
            unbalanced_to_balanced.insertLocalNeighbor(
               box_in_transit.d_box,
               box_in_transit.d_orig_box.getId());
            unassigned.erase(ni++);
         }
         else {
            // Leave this box in unassigned for notifying originating
            // process of where it landed.
            ++ni;
         }
      }

   }

   t_local_balancing->stop();

   t_load_distribution->stop();


   /*
    * Finish messages before starting edge info exchange.
    * We have only sends to complete, so it should not take
    * long to advance them all to completion.
    */
   t_finish_sends->start();
   child_send_stage.advanceAll();
   parent_send_stage.advanceAll();
   t_finish_sends->stop();
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < num_children; ++i) {
      TBOX_ASSERT(child_sends[i].isDone());
      TBOX_ASSERT(child_recvs[i].isDone());
   }
   if (parent_send != NULL) {
      TBOX_ASSERT(parent_send->isDone());
      TBOX_ASSERT(parent_recv->isDone());
   }
#endif


   /*
    * Ranks that we exported work to and imported work from, for use
    * when constructing semilocal unbalanced--->balanced
    * relationships.
    */
   std::vector<int> exported_to;
   std::vector<int> imported_from;

   for ( int ichild=0; ichild<num_children; ++ichild ) {
      if ( child_load_data[ichild].d_load_imported > 0 ) {
         imported_from.push_back(child_recvs[ichild].getPeerRank());
      }
      if ( child_load_data[ichild].d_load_exported > 0 ) {
         exported_to.push_back(child_sends[ichild].getPeerRank());
      }
   }

   if ( my_load_data.d_load_imported > 0 ) {
      imported_from.push_back(parent_send->getPeerRank());
   }
   if ( my_load_data.d_load_exported > 0 ) {
      exported_to.push_back(parent_send->getPeerRank());
   }

   constructSemilocalUnbalancedToBalanced(
      unbalanced_to_balanced,
      exported_to,
      imported_from,
      unassigned );


   if (d_check_connectivity) {
      const hier::OverlapConnectorAlgorithm oca;
      tbox::plog
      << "TreeLoadBalancerOld checking unbalanced-balanced connectivity."
      << std::endl;
      int errs = 0;
      if (oca.checkOverlapCorrectness(unbalanced_to_balanced, true, true)) {
         ++errs;
         tbox::perr << "Error found in unbalanced_to_balanced!\n";
      }
      if (oca.checkOverlapCorrectness(balanced_to_unbalanced, true, true)) {
         ++errs;
         tbox::perr << "Error found in balanced_to_unbalanced!\n";
      }
      if (unbalanced_to_balanced.checkTransposeCorrectness(
             balanced_to_unbalanced)) {
         ++errs;
         tbox::perr << "Error found in balanced-unbalanced transpose!\n";
      }
      if (errs != 0) {
         TBOX_ERROR(
            "Errors in load balance mapping found."
            << "balance_box_level:\n" << balance_box_level.format("", 2)
            << "balanced_box_level:\n" << balanced_box_level.format("", 2)
            << "unbalanced_to_balanced:\n" << unbalanced_to_balanced.format("", 2)
            << "balanced_to_unbalanced:\n" << balanced_to_unbalanced.format("", 2));
      }
   }

   if (d_check_map) {
      const hier::MappingConnectorAlgorithm mca;
      if (mca.findMappingErrors(unbalanced_to_balanced) != 0) {
         TBOX_ERROR(
            "TreeLoadBalancerOld::loadBalanceWithinRankGroup Mapping errors found in unbalanced_to_balanced!");
      }
      if (unbalanced_to_balanced.checkTransposeCorrectness(
             balanced_to_unbalanced)) {
         TBOX_ERROR(
            "TreeLoadBalancerOld::loadBalanceWithinRankGroup Transpose errors found!");
      }
   }


   delete[] child_load_data;
   destroyAsyncCommObjects(child_sends, parent_send);
   destroyAsyncCommObjects(child_recvs, parent_recv);

   t_get_map->stop();

   if (anchor_to_balance.isFinalized()) {
      t_use_map->start();
      const hier::MappingConnectorAlgorithm mca;
      mca.modify(
         anchor_to_balance,
         balance_to_anchor,
         unbalanced_to_balanced,
         balanced_to_unbalanced,
         &balance_box_level,
         &balanced_box_level);
      t_use_map->stop();
   } else {
      hier::BoxLevel::swap(balance_box_level, balanced_box_level);
   }

   return;
}



/*
 *************************************************************************
 *
 * This is the tree load balancer's reassignment method, essentially a
 * two-bin load balancer.  Given two sets of BoxInTransit (the bins)
 * and an amount of work to move from one set to the other, this
 * method makes a best effort to effect the work transfer between the
 * two bins.  It can move BoxInTransit between given sets and, if
 * needed, break some BoxInTransit up to move part of the work.
 *
 * This method is purely local--it reassigns the load but does not
 * communicate the change to any remote process.
 *
 *************************************************************************
 */
int
TreeLoadBalancerOld::reassignLoads(
   TransitSet& src,
   TransitSet& dst,
   hier::LocalId& next_available_index,
   int ideal_transfer ) const
{
   if (d_print_steps) {
      tbox::plog << "  reassignLoads attempting to reassign "
                 << ideal_transfer << " from src to dst."
                 << std::endl;
   }

   int actual_transfer = 0;

   if ((ideal_transfer >= 0 && src.empty()) ||
       (ideal_transfer <= 0 && dst.empty())) {
      return actual_transfer;
   }

   t_reassign_loads->start();

   /*
    * The algorithm cycles through a do-loop.  Each time around, we try
    * to swap some BoxInTransit between src and dst until we cannot improve the
    * actual_transfer any further.  Then, we try breaking up some BoxInTransit to
    * improve the results.  If we break some BoxInTransit, we generate some more
    * swapping options that were not there before, so we loop back to
    * try swapping again.
    *
    * If a break phase does not break any Box (and does not generate more
    * swap options), the loop will stop making changes.  We exit the loop
    * at that point (and whenever we reached the ideal transfer).
    */
   do {

      /*
       * Try to balance load through swapping.
       */
      int swap_transfer = shiftLoadsBySwapping(
         src,
         dst,
         ideal_transfer - actual_transfer );

      actual_transfer += swap_transfer;

      if (d_print_steps) {
         double balance_penalty = computeBalancePenalty(
            src,
            dst,
            actual_transfer - ideal_transfer);
         tbox::plog << "  Balance penalty after shiftLoadsBySwapping = "
                    << balance_penalty
                    << ", needs " << (ideal_transfer - actual_transfer)
                    << " more with " << src.size() << " source and "
                    << dst.size() << " dst Boxes remaining."
                    << std::endl;
      }

      if (actual_transfer == ideal_transfer) break;

      /*
       * Assuming that we did the best we could, swapping
       * some BoxInTransit without breaking any, we now break up a Box
       * in the overloaded side for partial transfer to the
       * underloaded side.
       */
      int brk_transfer = shiftLoadsByBreaking(
         src,
         dst,
         next_available_index,
         ideal_transfer - actual_transfer );
      actual_transfer += brk_transfer;

      if (d_print_steps) {
         double balance_penalty = computeBalancePenalty(
            src,
            dst,
            actual_transfer - ideal_transfer);
         tbox::plog << "    Balance penalty after shiftLoadsByBreaking = "
                    << balance_penalty
                    << ", needs " << (ideal_transfer - actual_transfer)
                    << " more with " << src.size() << " source and "
                    << dst.size() << " dst Boxes remaining."
                    << std::endl;
      }
      if (brk_transfer == 0) {
         /*
          * If no box can be broken to improve the actual_transfer,
          * there is nothing further we can do.  The swap phase, tried
          * before the break phase, also generated no transfer, so
          * there's no point trying again.  Break out now to save
          * retrying the swap phase.
          */
         break;
      }

      /*
       * Now that we have broken up a Box, redo this loop to
       * see if swapping can produce a better result.
       */
   } while (ideal_transfer != actual_transfer);

   t_reassign_loads->stop();

   return actual_transfer;
}



/*
 *************************************************************************
 *************************************************************************
 */
void
TreeLoadBalancerOld::packSubtreeLoadData(
   std::vector<int>& msg,
   const SubtreeLoadData& load_data) const
{
   t_pack_load->start();
   msg.push_back(load_data.d_num_procs);
   msg.push_back(load_data.d_total_work);
   msg.push_back(load_data.d_load_exported);
   const TransitSet& for_export = load_data.d_for_export;
   msg.push_back( static_cast<int>(for_export.size()));
   for (TransitSet::const_iterator
        ni = for_export.begin(); ni != for_export.end(); ++ni) {
      const BoxInTransit& box_in_transit = *ni;
      int i0 = static_cast<int>(msg.size());
      msg.insert(msg.end(), box_in_transit.commBufferSize(), 0);
      box_in_transit.putToIntBuffer(&msg[i0]);
   }
   t_pack_load->stop();
}



/*
 *************************************************************************
 *************************************************************************
 */
void
TreeLoadBalancerOld::unpackSubtreeLoadData(
   SubtreeLoadData& load_data,
   TransitSet& receiving_bin,
   hier::LocalId& next_available_index,
   const int* received_data,
   int received_data_length ) const
{
#ifndef DEBUG_CHECK_ASSERTIONS
   NULL_USE(received_data_length);
#endif

   t_unpack_load->start();
   const int *buffer = received_data;
   load_data.d_num_procs = *(buffer++);
   load_data.d_total_work = *(buffer++);
   load_data.d_load_imported = *(buffer++);
   const int num_boxes = *(buffer++);
   /*
    * As we pull each BoxInTransit out, give it a new id that reflects
    * its new owner.  Place all BoxInTransits in the receiving_bin.
    */
   BoxInTransit received_box(d_dim);
   for (int i = 0; i < num_boxes; ++i) {
      buffer = received_box.getFromIntBuffer(buffer);
      BoxInTransit renamed_box(received_box,
                               received_box.getBox(),
                               d_mpi.getRank(),
                               next_available_index);
      next_available_index += 2 + d_degree;
      receiving_bin.insert(renamed_box);
   }
   t_unpack_load->stop();
   TBOX_ASSERT( buffer == received_data + received_data_length );
}





/*
 *************************************************************************
 * Construct semilocal relationships in unbalanced--->balanced
 * Connector.  Constructing these relationships require communication
 * to determine where exported work ended up.
 *************************************************************************
 */
void
TreeLoadBalancerOld::constructSemilocalUnbalancedToBalanced(
   hier::Connector &unbalanced_to_balanced,
   const std::vector<int> &export_dsts,
   const std::vector<int> &import_srcs,
   const TreeLoadBalancerOld::TransitSet &kept_imports ) const
{
   t_construct_semilocal->start();

   /*
    * Determine edges from unbalanced to balanced BoxLevels by sending
    * balanced Boxes back to the owners of the unbalanced Boxes that
    * originated them.  Use the proc_hist data from each balanced
    * BoxInTransit to send it back along the path it took to get to
    * the process that eventually kept it.
    *
    * Any process we exported work to should send back a message about
    * where that work went.  Any process we imported work from expects
    * us to send a message about what happened to that work.  So we
    * receive messages from procs in export_dsts and send messages to
    * procs in import_srcs.
    *
    * For work that originated locally, construct relationships when
    * an incoming message tells us where that work ended up.  For work
    * that originated elsewhere, forward the information to the
    * process that sent it to us.
    *
    * TODO: Implement alternate edge generation method where terminal
    * owners send edge info directly back to initial owners.  No need
    * for trails.  Initial owners receive from MPI_ANY and stops
    * receiving when it can account for the number of cells it started
    * with.  We should retain both methods for edge generation because
    * we don't know which is faster on any given machine.  BTNG.
    */


   /*
    * Set up objects for asynchronous communication.
    *
    * NOTE: In the code review, it was suggested that the setting up
    * of asynchronous communication objects can be factored out for
    * reuse.
    */
   tbox::AsyncCommStage comm_stage;
   tbox::AsyncCommPeer<int> *importer_comms = import_srcs.empty() ? NULL : new tbox::AsyncCommPeer<int>[import_srcs.size()];
   tbox::AsyncCommPeer<int> *exporter_comms = export_dsts.empty() ? NULL : new tbox::AsyncCommPeer<int>[export_dsts.size()];

   for ( size_t i=0; i<export_dsts.size(); ++i ) {
      exporter_comms[i].initialize(&comm_stage);
      exporter_comms[i].setPeerRank(export_dsts[i]);
      exporter_comms[i].setMPI(d_mpi);
      exporter_comms[i].setMPITag(TreeLoadBalancerOld_EDGETAG0,
                                  TreeLoadBalancerOld_EDGETAG1);
      exporter_comms[i].limitFirstDataLength(TreeLoadBalancerOld_FIRSTDATALEN);
      exporter_comms[i].beginRecv();
      if ( exporter_comms[i].isDone() ) {
         exporter_comms[i].pushToCompletionQueue();
      }
   }

   for ( size_t i=0; i<import_srcs.size(); ++i ) {
      importer_comms[i].initialize(&comm_stage);
      importer_comms[i].setPeerRank(import_srcs[i]);
      importer_comms[i].setMPI(d_mpi);
      importer_comms[i].setMPITag(TreeLoadBalancerOld_EDGETAG0,
                                  TreeLoadBalancerOld_EDGETAG1);
      importer_comms[i].limitFirstDataLength(TreeLoadBalancerOld_FIRSTDATALEN);
   }


   // Buffers for staging data to send.  Indexed by recipient rank.
   std::map<int,std::vector<int> > outgoing_messages;


   /*
    * Prepare messages to tell processes that sent work to us about
    * the work that we kept.
    */
   for (TransitSet::const_iterator ni = kept_imports.begin();
        ni != kept_imports.end(); ++ni) {
      const BoxInTransit &box_in_transit = *ni;
      TBOX_ASSERT(!box_in_transit.d_proc_hist.empty());
      TBOX_ASSERT(box_in_transit.d_orig_box.getOwnerRank() != d_mpi.getRank());
      box_in_transit.packForPreviousOwner(outgoing_messages);
   }

   /*
    * Receive and unpack incoming messages.
    */
   while ( comm_stage.numberOfCompletedMembers() > 0 ||
           comm_stage.advanceSome() ) {

      tbox::AsyncCommPeer<int>* peer_comm =
         dynamic_cast<tbox::AsyncCommPeer<int> *>(comm_stage.popCompletionQueue());
      TBOX_ASSERT(peer_comm != NULL);

#ifdef DEBUG_CHECK_ASSERTIONS
      size_t j;
      for ( j=0; j<export_dsts.size(); ++j ) {
         if ( export_dsts[j] == peer_comm->getPeerRank() ) break;
      }
      TBOX_ASSERT( j < export_dsts.size() );
#endif

      const int* received_data = peer_comm->getRecvData();
      int received_data_length = peer_comm->getRecvSize();
      unpackAndRouteNeighborhoodSets(
         outgoing_messages,
         unbalanced_to_balanced,
         received_data,
         received_data_length );
   }



   /*
    * Send outgoing messages.
    */
   TBOX_ASSERT( outgoing_messages.size() == import_srcs.size() );
   for ( size_t i=0; i<import_srcs.size(); ++i ) {
      tbox::AsyncCommPeer<int> &peer_comm = importer_comms[i];
      std::map<int,std::vector<int> >::const_iterator recip =
         outgoing_messages.find(peer_comm.getPeerRank());
      TBOX_ASSERT( recip != outgoing_messages.end() );
      const std::vector<int> &message = recip->second;
      peer_comm.beginSend( &message[0], static_cast<int>(message.size()) );
   }

   t_construct_semilocal_comm_wait->start();
   t_finish_sends->start();
   comm_stage.advanceAll();
   t_finish_sends->stop();
   t_construct_semilocal_comm_wait->stop();

   delete [] importer_comms;
   delete [] exporter_comms;

   t_construct_semilocal->stop();
   return;
}




/*
 *************************************************************************
 * Unpack BoxInTransit objects from relationship message and either
 * use data to set up local edges in unbalanced--->balanced or route
 * data in the direction of the origin of the BoxInTransit.
 *************************************************************************
 */
void
TreeLoadBalancerOld::unpackAndRouteNeighborhoodSets(
   std::map<int,std::vector<int> > &outgoing_messages,
   hier::Connector& unbalanced_to_balanced,
   const int* received_data,
   int received_data_length ) const
{
   t_unpack_edge->start();

   const int* beg = received_data;
   const int* end = received_data + received_data_length;

   while (beg < end) {

      const BoxInTransit box_in_transit(beg, d_dim);
      if (d_print_edge_steps) {
         tbox::plog << "Unpacked edge " << box_in_transit << std::endl;
      }

      TBOX_ASSERT(beg <= end);

      // Empty proc_hist must mean that local process is the originator.
      TBOX_ASSERT(box_in_transit.d_proc_hist.empty() ==
                  ( box_in_transit.d_orig_box.getOwnerRank() ==
                    d_mpi.getRank() ) );

      /*
       * If the local process originated this box, generate the
       * relationship here.  Else, pack the box to forward to its
       * previous owner (and eventually to its originator).
       */

      if (box_in_transit.d_orig_box.getOwnerRank() == d_mpi.getRank()) {

         if (d_print_edge_steps) {
            tbox::plog << "Keeping edge " << box_in_transit << std::endl;
         }

         unbalanced_to_balanced.insertLocalNeighbor(
            box_in_transit.d_box,
            box_in_transit.d_orig_box.getId());

      } else {
         box_in_transit.packForPreviousOwner(outgoing_messages);
      }
   }
   t_unpack_edge->stop();
}



/*
 *************************************************************************
 * Set the MPI commuicator.  If there's a private communicator, free
 * it first.  It's safe to free the private communicator because no
 * other code have access to it.
 *************************************************************************
 */
void
TreeLoadBalancerOld::setSAMRAI_MPI(
   const tbox::SAMRAI_MPI& samrai_mpi)
{
   if (samrai_mpi.getCommunicator() == tbox::SAMRAI_MPI::commNull) {
      TBOX_ERROR("TreeLoadBalancerOld::setSAMRAI_MPI error: Given\n"
         << "communicator is invalid.");
   }

   if ( d_mpi_is_dupe ) {
      d_mpi.freeCommunicator();
   }

   // Enable private communicator.
   d_mpi.dupCommunicator(samrai_mpi);
   d_mpi_is_dupe = true;
}



/*
 *************************************************************************
 * Set the MPI commuicator.
 *************************************************************************
 */
void
TreeLoadBalancerOld::freeMPICommunicator()
{
   if ( d_mpi_is_dupe ) {
      // Free the private communicator (if MPI has not been finalized).
      int flag;
      tbox::SAMRAI_MPI::Finalized(&flag);
      if (!flag) {
         d_mpi.freeCommunicator();
      }
   }
   d_mpi.setCommunicator(tbox::SAMRAI_MPI::commNull);
   d_mpi_is_dupe = false;
}



/*
 *************************************************************************
 * RankGroups load-balances within their membership and ignore other
 * groups.  When we balance over multiple cycles, the RankGroup for
 * the local process depends on the cycle, as computed by this method.
 *
 * The RankGroup size increases exponentially with the cycle number
 * such that for the last cycle the rank group includes all processes
 * in d_mpi.  It's a heuristic formula, as follows:
 *
 * Partition all ranks into similar sized groups.  With each cycle,
 * the group size grows exponentially while the number of groups
 * shrinks.  The last cycle_number has a single group of
 * d_mpi.getSize() processors.
 *
 * Let m = number of cycles
 * i = cycle number, [0,m)
 * p = communicator size
 *
 * The group size is p^((i+1)/m)
 *************************************************************************
 */
void
TreeLoadBalancerOld::createBalanceRankGroupBasedOnCycles(
   tbox::RankGroup &rank_group,
   int &number_of_groups,
   int &group_num,
   const int cycle_number,
   const int number_of_cycles) const
{

   /*
    * Compute the number of group and, implicitly, the group sizes.
    * Tiny groups tend to leave the members with possibly large
    * overloads.  In order to make all groups similar in size we round
    * down the number of groups (and round up the group size).
    */
   number_of_groups =
      static_cast<int>(pow(static_cast<double>(d_mpi.getSize()),
                           1.0 - double(cycle_number + 1) / number_of_cycles));

   /*
    * All groups will have a base population count of
    * d_mpi.getSize()/number_of_groups.  The remainder from the
    * integer division is distributed to a subset of groups, starting
    * from group 0, so these groups will have one more than the base.
    */
   const int base_group_size = d_mpi.getSize() / number_of_groups;
   const int first_base_sized_group = d_mpi.getSize() % number_of_groups;
   const int first_rank_in_base_sized_group =
      first_base_sized_group * (1 + base_group_size);

   if (d_mpi.getRank() < first_rank_in_base_sized_group) {
      group_num = d_mpi.getRank() / (1 + base_group_size);
      const int group_first_rank = group_num * (1 +base_group_size);
      rank_group.setMinMax( group_first_rank,
                            group_first_rank + base_group_size );
   } else {
      group_num = first_base_sized_group
         + (d_mpi.getRank() - first_rank_in_base_sized_group) / base_group_size;
      const int group_first_rank = first_rank_in_base_sized_group +
         (group_num - first_base_sized_group)*(1+base_group_size);
      rank_group.setMinMax( group_first_rank,
                            group_first_rank + base_group_size - 1 );
   }

   return;
}



/*
 *************************************************************************
 * Set up the asynchronous communication objects for the process tree
 * containing ranks defined by the RankGroup.
 *
 * The process tree lay-out is defined by the BalancedDepthFirstTree
 * class, thus defining parent and children of the local process.
 * This method sets the AsyncCommPeer objects for communication with
 * children and parent.
 *************************************************************************
 */
void
TreeLoadBalancerOld::setupAsyncCommObjects(
   tbox::AsyncCommStage& child_stage,
   tbox::AsyncCommPeer<int> *& child_comms,
   tbox::AsyncCommStage& parent_stage,
   tbox::AsyncCommPeer<int> *& parent_comm,
   const tbox::RankGroup &rank_group,
   const tbox::BalancedDepthFirstTree &bdfs ) const
{

   child_comms = parent_comm = NULL;

   const int num_children = bdfs.getNumberOfChildren();

   if ( num_children > 0 ) {

      child_comms = new tbox::AsyncCommPeer<int>[num_children];

      for (int child_num = 0; child_num < num_children; ++child_num) {

         const int child_rank_in_grp = bdfs.getChildRank(child_num);
         const int child_true_rank = rank_group.getMappedRank(child_rank_in_grp);

         child_comms[child_num].initialize(&child_stage);
         child_comms[child_num].setPeerRank(child_true_rank);
         child_comms[child_num].setMPI(d_mpi);
         child_comms[child_num].setMPITag(TreeLoadBalancerOld_LOADTAG0,
                                          TreeLoadBalancerOld_LOADTAG1);
         child_comms[child_num].limitFirstDataLength(TreeLoadBalancerOld_FIRSTDATALEN);
      }
   }

   if (bdfs.getParentRank() != bdfs.getInvalidRank()) {

      const int parent_rank_in_grp = bdfs.getParentRank();
      int parent_true_rank = rank_group.getMappedRank(parent_rank_in_grp);

      parent_comm = new tbox::AsyncCommPeer<int>;
      parent_comm->initialize(&parent_stage);
      parent_comm->setPeerRank(parent_true_rank);
      parent_comm->setMPI(d_mpi);
      parent_comm->setMPITag(TreeLoadBalancerOld_LOADTAG0,
         TreeLoadBalancerOld_LOADTAG1);
      parent_comm->limitFirstDataLength(TreeLoadBalancerOld_FIRSTDATALEN);

   }

   return;
}



/*
 *************************************************************************
 *************************************************************************
 */
void
TreeLoadBalancerOld::destroyAsyncCommObjects(
   tbox::AsyncCommPeer<int> *& child_comms,
   tbox::AsyncCommPeer<int> *& parent_comm) const
{
   if (d_mpi.getSize() == 1) {
      TBOX_ASSERT(child_comms == NULL);
      TBOX_ASSERT(parent_comm == NULL);
   } else {
      if ( child_comms != NULL ) {
         delete[] child_comms;
      }
      if ( parent_comm != NULL ) {
         delete parent_comm;
      }
      child_comms = parent_comm = NULL;
   }
}



/*
 *************************************************************************
 * Attempt to shift a specified ammount of load from one TransitSet to
 * another by breaking a single box from the overloaded set.  Examine
 * multiple NodeInTransit and breakages to
 *
 * - filter out the ones that worsens the balance (see breakOffLoad) for
 * reasons why this can happen.
 *
 * - find the one that results in the smallest combinedBreakingPenalty.
 *
 * Return whether any changes were made.
 *************************************************************************
 */
int
TreeLoadBalancerOld::shiftLoadsByBreaking(
   TransitSet& src,
   TransitSet& dst,
   hier::LocalId& next_available_index,
   const int ideal_transfer ) const
{
   int actual_transfer = 0;

   if (ideal_transfer < 0) {
      // The logic below does not handle bi-directional transfers, so handle it here.
      actual_transfer = -shiftLoadsByBreaking(
         dst,
         src,
         next_available_index,
         -ideal_transfer );
      return actual_transfer;
   }

   TBOX_ASSERT(src.size() + dst.size() > 0);

   t_shift_loads_by_breaking->start();

   if (d_print_steps) {
      tbox::plog << "    shiftLoadsByBreaking asked to break off "
                 << ideal_transfer
                 << " from one of " << src.size()
                 << " source Boxes to add to set of " << dst.size()
                 << " Boxes."
                 << std::endl;
   }

   /*
    * The best results so far (from not transfering anything).
    */
   TransitSet best_src = src;
   TransitSet best_dst = dst;
   int best_actual_transfer = 0;

   double best_combined_penalty = combinedBreakingPenalty(
      computeBalancePenalty(best_src, best_dst, static_cast<double>(ideal_transfer)
                            - best_actual_transfer),
      computeSurfacePenalty(best_src, best_dst),
      computeSlenderPenalty(best_src, best_dst));

   /*
    * Scale the pre-cut penalty.  Scaling makes this method more
    * agressive about producing a cut.  Sometimes the uncut penalty can
    * be low enough to prevent a cut, but the cut may be required for
    * reasonable balancing.  The down side of agressive cutting is that
    * it tends to produce more slivers (but not terribly many).
    * This exchange of load balance and slivers may not be easily
    * eliminated by the weights in the penalty functions.
    */
   if (d_print_steps) {
      tbox::plog.unsetf(std::ios::fixed | std::ios::scientific);
      tbox::plog.precision(6);
      tbox::plog << "    Uncut penalty: " << best_combined_penalty
                 << " scaled by " << d_precut_penalty_wt << " to "
                 << best_combined_penalty * d_precut_penalty_wt
                 << std::endl;
   }
   best_combined_penalty *= d_precut_penalty_wt;

   bool found_breakage = false;

   std::vector<hier::Box> breakoff;
   std::vector<hier::Box> leftover;
   double breakoff_amt;
   for (TransitSet::iterator si = src.begin();
        si != src.end() &&
        (!found_breakage || si->d_load >= ideal_transfer); ++si) {

      const BoxInTransit& candidate = *si;

      breakOffLoad(
         breakoff,
         leftover,
         breakoff_amt,
         candidate.d_box,
         ideal_transfer );

      if (!breakoff.empty()) {

         const BoxInTransit& brk_box_in_transit = *si;

         const bool improves_balance =
            tbox::MathUtilities<double>::Abs(
               static_cast<double>(ideal_transfer) - breakoff_amt) <
            (ideal_transfer - tbox::MathUtilities<double>::getEpsilon());

         if (d_print_steps) {
            tbox::plog << "    Potential to replace " << brk_box_in_transit << " with "
                       << breakoff.size() << " breakoff Boxes and "
                       << leftover.size() << " leftover Boxes."
                       << "  improves_balance=" << improves_balance
                       << std::endl;
         }

         if (!improves_balance) {
            // Reject.
            continue;
         }

         /*
          * Trial modifications of src and dst, for evaluating
          * combined penalties.
          */
         TransitSet trial_src = src;
         TransitSet trial_dst = dst;
         TBOX_ASSERT(trial_src.size() + trial_dst.size() > 0);
         int trial_actual_transfer = actual_transfer;
         trial_src.erase(candidate);

         /*
          * Put breakoff in trial_dst and leftover back into trial_src.
          */
         for (std::vector<hier::Box>::const_iterator bi = breakoff.begin();
              bi != breakoff.end();
              ++bi) {
            BoxInTransit give_box_in_transit(
               brk_box_in_transit,
               *bi,
               d_mpi.getRank(),
               next_available_index);
            give_box_in_transit.d_load = static_cast<int>(computeLoad(
               give_box_in_transit.d_orig_box,
               give_box_in_transit.getBox()));
            next_available_index += 2 + d_degree;
            trial_dst.insert(give_box_in_transit);
            trial_actual_transfer += give_box_in_transit.d_load;
            if (d_print_steps) {
               tbox::plog << "    Breakoff box " << *bi << bi->numberCells()
                          << '|' << bi->size()
                          << " -> " << give_box_in_transit
                          << std::endl;
            }
         }
         TBOX_ASSERT(trial_src.size() + trial_dst.size() > 0);
         for (std::vector<hier::Box>::const_iterator bi = leftover.begin();
              bi != leftover.end();
              ++bi) {
            BoxInTransit keep_box_in_transit(
               brk_box_in_transit,
               *bi,
               d_mpi.getRank(),
               next_available_index);
            keep_box_in_transit.d_load = static_cast<int>(computeLoad(
                  keep_box_in_transit.d_orig_box,
                  keep_box_in_transit.getBox()));
            next_available_index += 2 + d_degree;
            trial_src.insert(keep_box_in_transit);
            if (d_print_steps) {
               tbox::plog << "    Leftover box " << *bi << " -> " << keep_box_in_transit
                          << std::endl;
            }
         }
         TBOX_ASSERT(trial_src.size() + trial_dst.size() > 0);
         trial_src.erase(brk_box_in_transit);

         /*
          * Compute the new penalty to see if it improves our best result so far.
          */
         TBOX_ASSERT(trial_src.size() + trial_dst.size() > 0);
         double trial_balance_penalty = computeBalancePenalty(trial_src,
                                                              trial_dst,
                                                              static_cast<int>(ideal_transfer) - trial_actual_transfer);
         double trial_surface_penalty = computeSurfacePenalty(trial_src,
               trial_dst);
         double trial_slender_penalty = computeSlenderPenalty(trial_src,
               trial_dst);
         double trial_combined_penalty = combinedBreakingPenalty(
               trial_balance_penalty,
               trial_surface_penalty,
               trial_slender_penalty);
         if (d_print_steps) {
            tbox::plog.unsetf(std::ios::fixed | std::ios::scientific);
            tbox::plog.precision(6);
            tbox::plog << "    Trial's imbalance: "
                       << (ideal_transfer - trial_actual_transfer)
                       << " balance,surface,slender,combined penalty: "
                       << trial_balance_penalty << ' '
                       << trial_surface_penalty << ' '
                       << trial_slender_penalty << ' '
                       << trial_combined_penalty
                       << std::endl;
         }

         if (d_print_steps) {
            if (trial_combined_penalty < best_combined_penalty) {
               tbox::plog << "    Keeping this trial." << std::endl;
            } else {
               tbox::plog << "    Rejecting this trial." << std::endl;
            }
         }

         if (trial_combined_penalty < best_combined_penalty) {
            found_breakage = true;
            best_actual_transfer = static_cast<int>(breakoff_amt);
            best_src = trial_src;
            best_dst = trial_dst;
            best_combined_penalty = trial_combined_penalty;
         }

      } else {
         if (d_print_steps) {
            tbox::plog << "    Break step could not break " << ideal_transfer
                       << std::endl;
         }
      }

   }

   if (found_breakage) {
      src.swap(best_src);
      dst.swap(best_dst);
      actual_transfer = best_actual_transfer;
   }

   t_shift_loads_by_breaking->stop();
   return actual_transfer;
}



/*
 *************************************************************************
 * Attempt to swap some NodesInTransit between 2 sets of
 * NodesInTransit (src and dst) to shift ideal_transfer work units.
 *
 * Transfering a BoxInTransit from one TransitSet to another
 * is considered a degenerate "swap" (a BoxInTransit is
 * swapped for nothing) handled by this function.
 *
 * This method can transfer load both ways.
 * ideal_transfer > 0 means to raise the load of dst
 * ideal_transfer < 0 means to raise the load of src
 * The iterative do loop may overshoot the ideal_transfer
 * and may have to swap to shift some of the load
 * back.
 *
 * Return whether any changes were made.
 *************************************************************************
 */
int
TreeLoadBalancerOld::shiftLoadsBySwapping(
   TransitSet& src,
   TransitSet& dst,
   int ideal_transfer ) const
{
   t_shift_loads_by_swapping->start();

   if (d_print_steps) {
      tbox::plog << "  Attempting to swap " << ideal_transfer << std::endl;
   }

   bool found_swap;

   int actual_transfer = 0;

   do {

      /*
       * Ammount we seek to transfer from hi to lo
       * (the "ideal" for this particular iteration).
       * Unlike ideal_transfer and actual_transfer, this quantity is positive.
       */
      int rem_transfer = ideal_transfer - actual_transfer;
      if (d_print_steps) {
         tbox::plog << "    Swap progress: " << actual_transfer
                    << " / " << ideal_transfer << " remaining transfer = "
                    << rem_transfer << std::endl;
      }

      found_swap = false;

      int swap_transfer;
      found_swap = swapLoadPair(
         src,
         dst,
         swap_transfer,
         rem_transfer );

      if (found_swap) {
         actual_transfer += swap_transfer;
      }

   } while (found_swap && (actual_transfer != ideal_transfer));

   if (d_print_steps) {
      tbox::plog << "  Final imbalance for swap " << ideal_transfer
      - actual_transfer << std::endl;
   }

   t_shift_loads_by_swapping->stop();

   return actual_transfer;
}



/*
 *************************************************************************
 * Find a BoxInTransit in src and a BoxInTransit in dst which when swapped
 * results in shifting close to ideal_shift from src to dst.
 * Return whether a swap pair was found.
 *
 * If isrc is set but idst is not, it means that isrc should
 * be moved to dst, but no dst should be moved back.  This is
 * the degenerate case of swapping isrc for a BoxInTransit with zero
 * load.
 *************************************************************************
 */
bool
TreeLoadBalancerOld::swapLoadPair(
   TransitSet& src,
   TransitSet& dst,
   int& actual_transfer,
   int ideal_transfer ) const
{
   if (ideal_transfer < 0) {
      // The logic below does not handle bi-directional transfers, so handle it here.
      bool rval = swapLoadPair(
         dst,
         src,
         actual_transfer,
         -ideal_transfer );
      actual_transfer = -actual_transfer;
      return rval;
   }

   t_find_swap_pair->start();

   if (d_print_steps) {
      tbox::plog << "  swapLoadPair looking for transfer of "
                 << ideal_transfer
                 << " between " << src.size() << "-BoxInTransit src and "
                 << dst.size() << "-BoxInTransit dst." << std::endl;
   }
   if (d_print_swap_steps) {
      tbox::plog << "    src (" << src.size() << "):" << std::endl;
      for (TransitSet::iterator si = src.begin(); si != src.end(); ++si) {
         tbox::plog << *si << std::endl;
      }
      tbox::plog << "    dst (" << dst.size() << "):" << std::endl;
      for (TransitSet::iterator si = dst.begin(); si != dst.end(); ++si) {
         tbox::plog << *si << std::endl;
      }
   }

   /*
    * Look for two swap options.  The "high side" option would
    * transfer at least ideal_transfer.  The "low side" option would
    * transfer up to ideal_transfer.
    *
    * Each option is defined by a box from src and a box from dst,
    * designated by the iterators src_hiside, dst_hiside, src_loside
    * and dst_loside.  src_hiside points to the box in the src for the
    * high-side transfer, and similarly for dst_hiside.  src_loside
    * points to the box in the src for the low-side transfer, and
    * similarly for dst_loside.
    *
    * Note that in the degenerate case, the dst box does not exist,
    * and the swap degenerates to moving a box from the src to the
    * dst.
    *
    * Compute the balance_penalty if high and low were swapped.  Keep
    * looking until we find the pair giving the lowest balance_penalty
    * on swapping.
    *
    * isrc and idst point to the current best pair to swap.  new_balance_penalty
    * is the balance_penalty if we swap them.
    *
    * src_test and dst_test are trial pairs to check to see if we can improve on
    * new_balance_penalty.
    *
    * We will look for two "best" pairs:
    */

   // Initialization indicating no swap pair found yet.
   TransitSet::iterator src_hiside = src.end();
   TransitSet::iterator dst_hiside = dst.end();
   TransitSet::iterator src_loside = src.end();
   TransitSet::iterator dst_loside = dst.end();

   // A dummy BoxInTransit for set searches.
   hier::Box dummy_box(d_dim);
   BoxInTransit dummy_search_target(d_dim);

   // Difference between swap results and ideal, >= 0
   int imbalance_loside = tbox::MathUtilities<int>::getMax();
   int imbalance_hiside = tbox::MathUtilities<int>::getMax();

   if (dst.empty()) {
      /*
       * There is no dst BoxInTransit, so the swap would
       * degnerate to moving a box from src to dst.  Find
       * the best src BoxInTransit to move.
       */
      dummy_search_target = BoxInTransit(hier::Box(dummy_box, hier::LocalId::getZero(), 0));
      dummy_search_target.d_load = ideal_transfer;
      TransitSet::iterator src_test = src.lower_bound(dummy_search_target);

      if (d_print_swap_steps) {
         tbox::plog << "  swapLoadPair with empty dst: ";
      }

      if (src_test != src.begin()) {
         src_hiside = src_test;
         --src_hiside;
         imbalance_hiside = src_hiside->d_load - ideal_transfer;
         if (d_print_swap_steps) {
            tbox::plog << "  hi src: " << (*src_hiside)
                       << " with transfer " << src_hiside->d_load
                       << ", off by " << imbalance_hiside;
         }
      }
      if (src_test != src.end()) {
         src_loside = src_test;
         imbalance_loside = ideal_transfer - src_loside->d_load;
         if (d_print_swap_steps) {
            tbox::plog << "  lo src: " << (*src_loside)
                       << " with transfer " << src_loside->d_load
                       << ", off by " << imbalance_loside;
         }
      }
      if (d_print_swap_steps) {
         tbox::plog << std::endl;
      }

   } else {

      /*
       * Start search through src beginning with the box whose load
       * exceed the biggest dst box by at least ideal_transfer.
       */
      dummy_search_target = *dst.begin();
      dummy_search_target.d_load += ideal_transfer;
      TransitSet::iterator src_beg = src.lower_bound(dummy_search_target);

      for (TransitSet::iterator src_test = src_beg; src_test != src.end(); ++src_test) {

         /*
          * Set dst_test pointing to where we should start looking in dst.
          * Look for a load less than the load of src_test by
          * ideal_transfer.
          */
         dummy_search_target = BoxInTransit(hier::Box(dummy_box, hier::LocalId::getZero(), 0));
         dummy_search_target.d_load = tbox::MathUtilities<int>::Max(
               src_test->d_load - ideal_transfer,
               0);
         TransitSet::iterator dst_test = dst.lower_bound(dummy_search_target);

         if (dst_test != dst.end()) {

            /*
             * lower_bound returned dst_test that would transfer >=
             * ideal_transfer when swapped with src_test.  Check
             * transfererence between src_test and dst_test for the
             * high-side transfer.  Also check the next smaller box in
             * dst for the low-side transfer.
             */

            // tmp_miss is the difference between the test swap amount and ideal_transfer.
            int tmp_miss = (src_test->d_load - dst_test->d_load) - ideal_transfer;
            TBOX_ASSERT(tmp_miss >= 0);

            if ((tmp_miss < imbalance_hiside)) {
               src_hiside = src_test;
               dst_hiside = dst_test;
               imbalance_hiside = tmp_miss;
               if (d_print_swap_steps) {
                  tbox::plog << "    new hi-swap pair: " << (*src_hiside)
                             << " & " << (*dst_hiside) << " with transfer "
                             << (src_hiside->d_load - dst_hiside->d_load)
                             << " missing by " << imbalance_hiside
                             << std::endl;
               }
            }

            if (dst_test != dst.begin()) {
               --dst_test; // Now, src_test and dst_test transferer by *less* than ideal_transfer.
               tmp_miss = ideal_transfer - (src_test->d_load - dst_test->d_load);
               TBOX_ASSERT(tmp_miss >= 0);
               if (tmp_miss < imbalance_loside) {
                  src_loside = src_test;
                  dst_loside = dst_test;
                  imbalance_loside = tmp_miss;
                  if (d_print_swap_steps) {
                     tbox::plog << "    new lo-swap pair: " << (*src_loside)
                                << " & " << (*dst_loside) << " with transfer "
                                << (src_loside->d_load - dst_loside->d_load)
                                << " missing by " << imbalance_loside
                                << std::endl;
                  }
               }
            }

         } else {

            /*
             * The ideal dst to swap is smaller than the smallest dst
             * box.  So the only choice is swapping src_test for nothing.
             * Chech this against the current high- and low-side choices.
             */
            if (src_test->d_load > ideal_transfer) {
               // Moving src_test to src is moving too much--hiside.
               int tmp_miss = src_test->d_load - ideal_transfer;
               TBOX_ASSERT(tmp_miss >= 0);
               if (tmp_miss < imbalance_hiside) {
                  src_hiside = src_test;
                  dst_hiside = dst.end();
                  imbalance_hiside = tmp_miss;
                  if (d_print_swap_steps) {
                     tbox::plog << "    new hi-swap source: " << (*src_hiside)
                                << " & " << "no dst" << " with transfer "
                                << (src_hiside->d_load)
                                << " missing by " << imbalance_hiside
                                << std::endl;
                  }
               }
            } else {
               // Moving src_test to src is moving (just right or) too little--loside.
               int tmp_miss = ideal_transfer - src_test->d_load;
               TBOX_ASSERT(tmp_miss >= 0);
               if (tmp_miss < imbalance_loside) {
                  src_loside = src_test;
                  dst_loside = dst.end();
                  imbalance_loside = tmp_miss;
                  if (d_print_swap_steps) {
                     tbox::plog << "    new lo-swap source: " << (*src_loside)
                                << " & " << "no dst" << " with transfer "
                                << (src_loside->d_load)
                                << " missing by " << imbalance_loside
                                << std::endl;
                  }
               }
               /*
                * Break out of the loop early because there is no
                * point checking smaller src boxes.
                */
               break;
            }
         }
      }

   }

   /*
    * Swapping does not produce new cuts, so it is ok to omit the penalties
    * arising from cutting.
    */
   double current_balance_penalty = static_cast<double>(ideal_transfer);
   double balance_penalty_loside = static_cast<double>(imbalance_loside);
   double balance_penalty_hiside = static_cast<double>(imbalance_hiside);

   if (d_print_swap_steps) {
      tbox::plog << "    Swap candidates give penalties (unswap,lo,hi): "
                 << current_balance_penalty << " , " << balance_penalty_loside
                 << " , " << balance_penalty_hiside << std::endl;
   }

   bool found_swap = false;
   TransitSet::iterator isrc;
   TransitSet::iterator idst;

   if (balance_penalty_loside < current_balance_penalty &&
       balance_penalty_loside <= balance_penalty_hiside) {
      isrc = src_loside;
      idst = dst_loside;
      found_swap = true;
      if (d_print_swap_steps) {
         tbox::plog << "    Taking loside." << std::endl;
      }
   } else if (balance_penalty_hiside < current_balance_penalty &&
              balance_penalty_hiside <= balance_penalty_loside) {
      isrc = src_hiside;
      idst = dst_hiside;
      found_swap = true;
      if (d_print_swap_steps) {
         tbox::plog << "    Taking hiside." << std::endl;
      }
   } else {
      if (d_print_swap_steps) {
         tbox::plog << "    Keeping original (no swap)." << std::endl;
      }
   }

   if (found_swap) {
      actual_transfer = isrc->d_load;
      if (idst != dst.end()) {
         actual_transfer -= idst->d_load;
      }
   }

   if (found_swap) {

      // We can improve balance_penalty by swapping isrc with idst.
      if (d_print_steps) {
         tbox::plog << "    Swapping " << actual_transfer << " units using ";
         if (isrc != src.end()) tbox::plog << *isrc;
         else tbox::plog << "X";
         tbox::plog << " <--> ";
         if (idst != dst.end()) tbox::plog << *idst;
         else tbox::plog << "X";
         tbox::plog << std::endl;
      }

      if (isrc != src.end()) {
         dst.insert(*isrc);
         src.erase(isrc);
      }
      if (idst != dst.end()) {
         src.insert(*idst);
         dst.erase(idst);
      }


   } else {
      if (d_print_steps) {
         tbox::plog << "    Cannot find swap pair for " << ideal_transfer
                    << " units." << std::endl;
      }
   }

   t_find_swap_pair->stop();
   return found_swap;
}



/*
 *************************************************************************
 * Master method for breaking off a load.
 *
 * Try different heuristics and pick the "best" way to break off a
 * load.  The best is defined as the one with the lowest combined
 * penalty.
 *
 * This method always return a breakage if at all possible, without
 * considering whether the break should be used.  For example,
 * requesting breakage of 1 cell in a 100x100 box might return a
 * breakage of a 100-cells sliver!
 *
 * Return whether a successful break was made.
 *************************************************************************
 */
bool
TreeLoadBalancerOld::breakOffLoad(
   std::vector<hier::Box>& breakoff,
   std::vector<hier::Box>& leftover,
   double& brk_load,
   const hier::Box& box,
   double ideal_load_to_break ) const
{
   TBOX_ASSERT(ideal_load_to_break > 0);

   /*
    * NOTE: We need in this method a way to weigh the
    * value of proximity to the ideal breakoff vs the
    * increased area of the cuts.  However, the weight
    * given to area-optimized cuts should be considered
    * only with real application performance data.
    *
    * NOTE: We can compute the amount of new box boundaries
    * generated by computing the box boundary before and
    * after, and subtracting.  Easier than reconstructing
    * the cuts from the box definitions.
    *
    * NOTE: We should weight negatively the production of
    * high surface-to-volume boxes.
    */

   t_break_off_load->start();

   breakoff.clear();
   leftover.clear();

   if (ideal_load_to_break < d_min_size.getProduct()) {
      /*
       * Assuming uniform load balancing, there is no hope
       * if breaking off a piece of the desired size.
       */
      if (d_print_break_steps) {
         tbox::plog << "      ideal_load_to_break " << ideal_load_to_break
                    << " < " << d_min_size.getProduct() << d_min_size
                    << ":  Cannot break Box " << box << std::endl;
      }
      t_break_off_load->stop();
      return false;
   }

   /*
    * To avoid repeated computations of bad cuts,
    * we precompute bad_cuts here to provide to
    * methods that actually use the information.
    */
   tbox::Array<tbox::Array<bool> > bad_cuts(d_dim.getValue());
   t_find_bad_cuts->start();
   hier::BoxUtilities::findBadCutPoints(bad_cuts,
      box,
      d_block_domain_boxes[box.getBlockId().getBlockValue()],
      d_bad_interval);
   t_find_bad_cuts->stop();

   // Penalty for not transfering ideal load.
   double best_balance_penalty = computeBalancePenalty(box,
         ideal_load_to_break);
   // Penalty for new surfaces generated (none generated yet).
   double best_surface_penalty = computeSurfacePenalty(box);
   // Penalty for slender boxes.
   double best_slender_penalty = computeSlenderPenalty(box);

   double best_combined_penalty = tbox::MathUtilities<double>::getMax();

   if (d_print_break_steps) {
      tbox::plog.unsetf(std::ios::fixed | std::ios::scientific);
      tbox::plog.precision(6);
      tbox::plog << "      pre-break imbalance: " << ideal_load_to_break
                 << " balance,surface,slender,combined penalties: "
                 << best_balance_penalty << ", " << best_surface_penalty
                 << ", "
                 << best_slender_penalty << ", " << best_combined_penalty
                 << std::endl;
   }

   brk_load = 0;
   bool found_any_break = false;

   {
      std::vector<hier::Box> planar_breakoff;
      std::vector<hier::Box> planar_leftover;
      double planar_brk_load;

      bool found_this_break = breakOffLoad_planar(
            planar_breakoff,
            planar_leftover,
            planar_brk_load,
            box,
            ideal_load_to_break,
            bad_cuts );

      if (found_this_break) {
         found_any_break = true;
         double planar_balance_penalty = computeBalancePenalty(planar_breakoff,
               planar_leftover,
               static_cast<double>(planar_brk_load - ideal_load_to_break));
         double planar_surface_penalty = computeSurfacePenalty(planar_breakoff,
               planar_leftover);
         double planar_slender_penalty = computeSlenderPenalty(planar_breakoff,
               planar_leftover);
         double planar_combined_penalty =
            combinedBreakingPenalty(planar_balance_penalty,
               planar_surface_penalty,
               planar_slender_penalty);
         if (d_print_break_steps) {
            tbox::plog.unsetf(std::ios::fixed | std::ios::scientific);
            tbox::plog.precision(6);
            tbox::plog << "      Planar-break broke off "
                       << planar_brk_load << " / " << ideal_load_to_break
                       << " from " << box << '|'
                       << box.numberCells() << '|'
                       << box.size() << " into "
                       << planar_breakoff.size()
                       << " breakoff: ";
            for (std::vector<hier::Box>::const_iterator bi =
                    planar_breakoff.begin();
                 bi != planar_breakoff.end();
                 ++bi) {
               tbox::plog << " " << *bi << '|' << bi->numberCells() << '|'
                          << bi->size();
            }
            tbox::plog << "\n        and " << planar_leftover.size()
                       << " leftover boxes:";
            for (std::vector<hier::Box>::const_iterator bi =
                    planar_leftover.begin();
                 bi != planar_leftover.end();
                 ++bi) {
               tbox::plog << " " << *bi << '|' << bi->numberCells() << '|'
                          << bi->size();
            }
            tbox::plog << "\n        imbalance: "
                       << (planar_brk_load - ideal_load_to_break)
                       << " balance,surface,slender,combined penalties: "
                       << planar_balance_penalty << ", "
                       << planar_surface_penalty << ", "
                       << planar_slender_penalty
                       << ", " << planar_combined_penalty << std::endl;
         }
         if (planar_combined_penalty < best_combined_penalty) {
            if (d_print_break_steps) {
               tbox::plog << "      Keeping planar cut result." << std::endl;
            }
            breakoff.swap(planar_breakoff);
            leftover.swap(planar_leftover);
            brk_load = planar_brk_load;
            best_balance_penalty = planar_balance_penalty;
            best_surface_penalty = planar_surface_penalty;
            best_slender_penalty = planar_slender_penalty;
            best_combined_penalty = planar_combined_penalty;
         } else {
            if (d_print_break_steps) {
               tbox::plog << "      Rejecting planar cut result." << std::endl;
            }
         }
      }
   }

   /*
    * If above cut algorithms fail to break or improve the penalty, try
    * more cutting algorithms.
    */
   {

      std::vector<hier::Box> cubic_breakoff;
      std::vector<hier::Box> cubic_leftover;
      double cubic_brk_load;

      bool found_this_break = breakOffLoad_cubic(
            cubic_breakoff,
            cubic_leftover,
            cubic_brk_load,
            box,
            ideal_load_to_break,
            bad_cuts );

      if (found_this_break) {
         found_any_break = true;
         double cubic_balance_penalty = computeBalancePenalty(
               cubic_breakoff,
               cubic_leftover,
               static_cast<double>(cubic_brk_load - ideal_load_to_break));
         double cubic_surface_penalty = computeSurfacePenalty(cubic_breakoff,
               cubic_leftover);
         double cubic_slender_penalty = computeSlenderPenalty(cubic_breakoff,
               cubic_leftover);
         double cubic_combined_penalty =
            combinedBreakingPenalty(cubic_balance_penalty,
               cubic_surface_penalty,
               cubic_slender_penalty);
         if (d_print_break_steps) {
            tbox::plog.unsetf(std::ios::fixed | std::ios::scientific);
            tbox::plog.precision(6);
            tbox::plog << "      Cubic-break broke off "
                       << cubic_brk_load << " / " << ideal_load_to_break
                       << " from " << box << '|'
                       << box.numberCells() << '|'
                       << box.size() << " into "
                       << cubic_breakoff.size()
                       << " breakoff: ";
            for (std::vector<hier::Box>::const_iterator bi =
                    cubic_breakoff.begin();
                 bi != cubic_breakoff.end();
                 ++bi) {
               tbox::plog << " " << *bi << '|' << bi->numberCells() << '|'
                          << bi->size();
            }
            tbox::plog << "\n        and " << cubic_leftover.size()
                       << " leftover boxes:";
            for (std::vector<hier::Box>::const_iterator bi =
                    cubic_leftover.begin();
                 bi != cubic_leftover.end();
                 ++bi) {
               tbox::plog << " " << *bi << '|' << bi->numberCells() << '|'
                          << bi->size();
            }
            tbox::plog << "\n        imbalance: "
                       << (cubic_brk_load - ideal_load_to_break)
                       << " balance,surface,slender,combined penalties: "
                       << cubic_balance_penalty << ", "
                       << cubic_surface_penalty << ", "
                       << cubic_slender_penalty
                       << ", " << cubic_combined_penalty << std::endl;
         }
         if (cubic_combined_penalty < best_combined_penalty) {
            if (d_print_break_steps) {
               tbox::plog << "      choosing breakOffLoad_cubic result."
                          << std::endl;
            }
            breakoff.swap(cubic_breakoff);
            leftover.swap(cubic_leftover);
            brk_load = cubic_brk_load;
            best_balance_penalty = cubic_balance_penalty;
            best_surface_penalty = cubic_surface_penalty;
            best_slender_penalty = cubic_slender_penalty;
            best_combined_penalty = cubic_combined_penalty;
         } else {
            if (d_print_break_steps) {
               tbox::plog << "      Rejecting cubic cut result." << std::endl;
            }
         }
      } else {
         if (d_print_break_steps) {
            tbox::plog << "      breakOffLoad_cubic could not break "
                       << ideal_load_to_break << " from " << box
                       << '/' << box.numberCells()
                       << '/' << box.numberCells().getProduct()
                       << std::endl;
         }
      }

   }

   t_break_off_load->stop();

   return found_any_break;
}



/*
 *************************************************************************
 * Measuring surface area of boxes is used in penalizing
 * the creation of new surfaces.
 *************************************************************************
 */

double
TreeLoadBalancerOld::computeBoxSurfaceArea(
   const std::vector<hier::Box>& boxes) const
{
   int surface_area = 0;
   for (std::vector<hier::Box>::const_iterator bi = boxes.begin();
        bi != boxes.end();
        ++bi) {
      const hier::Box& box = *bi;
      int volume = box.size();
      for (int d = 0; d < d_dim.getValue(); ++d) {
         surface_area += volume / box.numberCells(d);
      }
   }
   surface_area *= 2;
   return surface_area;
}



/*
 *************************************************************************
 * Measuring surface area of boxes is used in penalizing
 * the creation of new surfaces.
 *************************************************************************
 */

int
TreeLoadBalancerOld::computeBoxSurfaceArea(
   const hier::Box& box) const
{
   int surface_area = 0;
   int volume = box.size();
   for (int d = 0; d < d_dim.getValue(); ++d) {
      surface_area += volume / box.numberCells(d);
   }
   surface_area *= 2;
   return surface_area;
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted surface penalty for a box.  The
 * reference zero penalty is for a box of equal sides having the same
 * volume.
 *************************************************************************
 */

double
TreeLoadBalancerOld::computeSurfacePenalty(
   const std::vector<hier::Box>& a,
   const std::vector<hier::Box>& b) const
{
   double surface_penalty = 0;
   for (std::vector<hier::Box>::const_iterator bi = a.begin();
        bi != a.end();
        ++bi) {
      surface_penalty += computeSurfacePenalty(*bi);
   }
   for (std::vector<hier::Box>::const_iterator bi = b.begin();
        bi != b.end();
        ++bi) {
      surface_penalty += computeSurfacePenalty(*bi);
   }
   return surface_penalty;
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted surface penalty for a box.  The
 * reference zero penalty is for a box of equal sides having the same
 * volume.
 *************************************************************************
 */

double
TreeLoadBalancerOld::computeSurfacePenalty(
   const TransitSet& a,
   const TransitSet& b) const
{
   double surface_penalty = 0;
   for (TransitSet::const_iterator bi = a.begin(); bi != a.end(); ++bi) {
      surface_penalty += computeSurfacePenalty(bi->d_box);
   }
   for (TransitSet::const_iterator bi = b.begin(); bi != b.end(); ++bi) {
      surface_penalty += computeSurfacePenalty(bi->d_box);
   }
   return surface_penalty;
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted surface penalty for a box.  The
 * reference zero penalty is for a box of equal sides having the same
 * volume.
 *************************************************************************
 */

double
TreeLoadBalancerOld::computeSurfacePenalty(
   const hier::Box& a) const
{
   int boxvol = a.size();
   double surface_area = computeBoxSurfaceArea(a);
   double best_surface = 2 * d_dim.getValue() * pow(static_cast<double>(boxvol),
         static_cast<double>(d_dim.getValue() - 1) / d_dim.getValue());
   double surface_penalty = surface_area / best_surface - 1.0;
   surface_penalty = surface_penalty * surface_penalty; // Make it blow up.
   return surface_penalty;
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted slenderness penalty for two
 * containers of boxes.  The reference zero penalty refers to a box with
 * all sides the same length.
 *************************************************************************
 */

double
TreeLoadBalancerOld::computeSlenderPenalty(
   const std::vector<hier::Box>& a,
   const std::vector<hier::Box>& b) const
{
   double slender_penalty = 0;
   for (std::vector<hier::Box>::const_iterator bi = a.begin();
        bi != a.end();
        ++bi) {
      slender_penalty += computeSlenderPenalty(*bi);
   }
   for (std::vector<hier::Box>::const_iterator bi = b.begin();
        bi != b.end();
        ++bi) {
      slender_penalty += computeSlenderPenalty(*bi);
   }
   return slender_penalty;
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted slenderness penalty for two
 * containers of boxes.  The reference zero penalty refers to a box with
 * all sides the same length.
 *************************************************************************
 */

double
TreeLoadBalancerOld::computeSlenderPenalty(
   const TransitSet& a,
   const TransitSet& b) const
{
   double slender_penalty = 0;
   for (TransitSet::const_iterator bi = a.begin(); bi != a.end(); ++bi) {
      slender_penalty += computeSlenderPenalty(bi->d_box);
   }
   for (TransitSet::const_iterator bi = b.begin(); bi != b.end(); ++bi) {
      slender_penalty += computeSlenderPenalty(bi->d_box);
   }
   return slender_penalty;
}



/*
 *************************************************************************
 * Return non-dimensional volume-weighted slenderness penalty for two
 * containers of boxes.  The reference zero penalty refers to a box with
 * all sides the same length.
 *************************************************************************
 */

double
TreeLoadBalancerOld::computeSlenderPenalty(
   const hier::Box& a) const
{
   const hier::IntVector boxdim = a.numberCells();
   double slender_penalty = static_cast<double>(boxdim.max()) / boxdim.min() - 1.0;
   slender_penalty = slender_penalty * slender_penalty; // Make it blow up.
   return slender_penalty;
}



/*
 *************************************************************************
 * Break up box bursty against box solid and adds the pieces to list.
 * This version differs from that in BoxContainer in that it tries to minimize
 * slivers.
 *************************************************************************
 */

void
TreeLoadBalancerOld::burstBox(
   std::vector<hier::Box>& boxes,
   const hier::Box& bursty,
   const hier::Box& solid ) const
{
   /*
    * This method lacks logic to handle the case of solid not being
    * completely inside bursty.  That feature is not currently needed.
    */
   TBOX_ASSERT(bursty.contains(solid));

   const hier::IntVector solid_size = solid.numberCells();

   boxes.clear();
   hier::Box cutme = bursty;
   while (!cutme.isSpatiallyEqual(solid)) {

      int cut_dir = 999999;
      bool cut_above_solid = false; // Whether to slice off the piece above solid (vs below).
      /*
       * Find direction and place to cut.  To minimize slivers, cut
       * from cutme the thickest slab (in direction normal to cut)
       * possible.
       */
      int slab_thickness = 0;
      for (int d = 0; d < d_dim.getValue(); ++d) {
         if (cutme.numberCells(d) > solid_size(d)) {
            const int thickness_from_upper_cut = cutme.upper() (d)
               - solid.upper() (d);
            if (thickness_from_upper_cut > slab_thickness) {
               slab_thickness = thickness_from_upper_cut;
               cut_dir = d;
               cut_above_solid = true;
            }
            const int thickness_from_lower_cut = solid.lower() (d)
               - cutme.lower() (d);
            if (thickness_from_lower_cut > slab_thickness) {
               slab_thickness = thickness_from_lower_cut;
               cut_dir = d;
               cut_above_solid = false;
            }
         }
      }
      TBOX_ASSERT(cut_dir >= 0 && cut_dir < d_dim.getValue());

      hier::Box removeme = cutme;
      if (cut_above_solid) {
         cutme.upper() (cut_dir) = solid.upper() (cut_dir);
         removeme.lower() (cut_dir) = solid.upper() (cut_dir) + 1;
      } else {
         cutme.lower() (cut_dir) = solid.lower() (cut_dir);
         removeme.upper() (cut_dir) = solid.lower() (cut_dir) - 1;
      }

      boxes.push_back(removeme);

   }

   if (d_print_break_steps) {
      tbox::plog << "      burstBox: " << bursty << " = " << solid;
      for (std::vector<hier::Box>::const_iterator bi = boxes.begin();
           bi != boxes.end();
           bi++) {
         tbox::plog << " + " << *bi;
      }
      tbox::plog << std::endl;
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   for (std::vector<hier::Box>::const_iterator bi = boxes.begin();
        bi != boxes.end();
        bi++) {
      for (std::vector<hier::Box>::const_iterator bj = boxes.begin();
           bj != boxes.end();
           bj++) {
         if (bi != bj) {
            TBOX_ASSERT(!bi->intersects(*bj));
         }
      }
   }
   hier::BoxContainer l1(bursty);
   hier::BoxContainer l2(solid);
   for (std::vector<hier::Box>::const_iterator bi = boxes.begin();
        bi != boxes.end();
        bi++) {
      l2.pushFront(*bi);
   }
   l1.removeIntersections(l2);
   TBOX_ASSERT(l1.size() == 0);
   l2.removeIntersections(bursty);
   TBOX_ASSERT(l2.size() == 0);
#endif
}



/*
 *************************************************************************
 *************************************************************************
 */
double
TreeLoadBalancerOld::computeLocalLoads(
   const hier::BoxLevel& box_level) const
{
   // Count up workload.
   double load = 0.0;
   const hier::BoxContainer& boxes = box_level.getBoxes();
   for (hier::BoxContainer::const_iterator ni = boxes.begin();
        ni != boxes.end();
        ++ni) {
      double box_load = computeLoad(*ni);
      load += box_load;
   }
   return static_cast<double>(load);
}



/*
 *************************************************************************
 *
 * Print out all attributes of class instance for debugging.
 *
 *************************************************************************
 */

void
TreeLoadBalancerOld::printClassData(
   std::ostream& os) const
{
   os << "\nTreeLoadBalancerOld::printClassData..." << std::endl;
   os << "\nTreeLoadBalancerOld: this = "
      << (TreeLoadBalancerOld *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;

   int ln;

   os << "d_workload_data_id..." << std::endl;
   for (ln = 0; ln < d_workload_data_id.getSize(); ln++) {
      os << "    d_workload_data_id[" << ln << "] = "
         << d_workload_data_id[ln] << std::endl;
   }
}



/*
 *************************************************************************
 *
 * Read values (described in the class header) from input database.
 *
 *************************************************************************
 */

void
TreeLoadBalancerOld::getFromInput(
   const boost::shared_ptr<tbox::Database>& db)
{

   if (db) {

      d_print_steps =
         db->getBoolWithDefault("print_steps",
            d_print_steps);
      d_print_break_steps =
         db->getBoolWithDefault("print_break_steps",
            d_print_break_steps);
      d_print_swap_steps =
         db->getBoolWithDefault("print_swap_steps",
            d_print_swap_steps);
      d_print_edge_steps =
         db->getBoolWithDefault("print_edge_steps",
            d_print_edge_steps);
      d_check_connectivity =
         db->getBoolWithDefault("check_connectivity",
            d_check_connectivity);
      d_check_map =
         db->getBoolWithDefault("check_map",
            d_check_map);

      d_report_load_balance = db->getBoolWithDefault("report_load_balance",
            d_report_load_balance);
      d_barrier_before = db->getBoolWithDefault("barrier_before",
            d_barrier_before);
      d_barrier_after = db->getBoolWithDefault("barrier_after",
            d_barrier_after);

      d_n_root_cycles = db->getIntegerWithDefault("n_root_cycles",
            d_n_root_cycles);

      d_balance_penalty_wt = db->getDoubleWithDefault("balance_penalty_wt",
            d_balance_penalty_wt);
      d_surface_penalty_wt = db->getDoubleWithDefault("surface_penalty_wt",
            d_surface_penalty_wt);
      d_slender_penalty_wt = db->getDoubleWithDefault("slender_penalty_wt",
            d_slender_penalty_wt);
      d_precut_penalty_wt = db->getDoubleWithDefault("precut_penalty_wt",
            d_precut_penalty_wt);

   }
}



/*
 ***************************************************************************
 *
 ***************************************************************************
 */
void
TreeLoadBalancerOld::assertNoMessageForPrivateCommunicator() const
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
         TBOX_ERROR(
            "Library error!\n"
            << "TreeLoadBalancerOld detected before or\n"
            << "after using a private communicator that there\n"
            << "is a message yet to be received.  This is\n"
            << "an error because all messages using the\n"
            << "private communicator should have been\n"
            << "accounted for.  Message status:\n"
            << "source " << mpi_status.MPI_SOURCE << '\n'
            << "tag " << mpi_status.MPI_TAG << '\n'
            << "count " << count << " (assuming integers)\n"
            << "current tags: "
            << ' ' << TreeLoadBalancerOld_LOADTAG0 << ' '
            << TreeLoadBalancerOld_LOADTAG1
            );
      }
   }
}



/*
 *************************************************************************
 * Break off a load from a box by making a single planar cut across
 * the box's longest dimension.
 *
 * Currently assuming uniform load of one unit per cell.
 *
 * Return whether a successful break was made.
 *************************************************************************
 */
bool
TreeLoadBalancerOld::breakOffLoad_planar(
   std::vector<hier::Box>& breakoff,
   std::vector<hier::Box>& leftover,
   double& brk_load,
   const hier::Box& box,
   double ideal_load_to_break,
   const tbox::Array<tbox::Array<bool> >& bad_cuts ) const
{

   const tbox::Dimension dim(d_dim);

   if (d_print_break_steps) {
      tbox::plog << "      breakOffLoad_planar attempting to break "
                 << ideal_load_to_break << " from Box " << box
                 << " min_size=" << d_min_size << std::endl;
   }

   brk_load = 0;
   breakoff.clear();
   leftover.clear();

   const hier::IntVector& box_dims = box.numberCells();

   const int box_vol = box_dims.getProduct();

   if (box_vol <= ideal_load_to_break) {
      // Easy: break off everything.
      breakoff.push_back(box);
      brk_load = box_vol;
      if (d_print_break_steps) {
         tbox::plog << "      breakOffload_planar broke off entire Box "
                    << box
                    << std::endl;
      }
      return true;
   }

   /*
    * Determine ordering of box_dims from shortest to longest.
    */
   hier::IntVector sorted_dirs(dim);
   sorted_dirs.sortIntVector(box_dims);

   /*
    * best_difference is the difference between the best cut found and
    * ideal_load_to_break.  Initialize it for zero breakoff.
    */
   double best_difference = static_cast<double>(ideal_load_to_break);
   hier::Box best_breakoff_box(dim);
   hier::Box best_leftover_box(dim);

   for (int d = d_dim.getValue() - 1; d >= 0; --d) {

      /*
       * Search directions from longest to shortest
       * because we prefer to break across longest dir.
       */
      const int brk_dir = sorted_dirs(d_dim.getValue() - 1);

      const int brk_area = box_vol / box_dims(brk_dir);

      const tbox::Array<bool>& bad = bad_cuts[brk_dir];

      /*
       * Try rounding the break length down (round==0) and up (round==1),
       * looking for the best place to cut.
       */
      for (int round = 0; round <= 1; ++round) {

         /*
          * Rounding up or down must heed d_cut_factor,
          * so brk_len must be an integer multiple of d_cut_factor.
          */
         const int brk_len =
            (static_cast<int>(ideal_load_to_break / brk_area)
             / d_cut_factor(brk_dir) + round)
            * d_cut_factor(brk_dir);

         if (brk_len < d_min_size(brk_dir)) {
            // Breakoff box violates minimum size.
            continue;
         }

         if (box_dims(brk_dir) - brk_len > 0 &&
             box_dims(brk_dir) - brk_len < d_min_size(brk_dir)) {
            // Leftover box violates minimum size.
            continue;
         }

         const int brk_volume = brk_area * brk_len;

         /*
          * Compute the difference between the current breakage
          * and the ideal.
          */
         const double difference = brk_volume <= ideal_load_to_break ?
            static_cast<double>(ideal_load_to_break - brk_volume) :
            static_cast<double>(brk_volume - ideal_load_to_break);

         if (difference < best_difference) {
            // This cut gives better difference, if it can be done.

            TBOX_ASSERT(brk_len >= 0 && brk_len <= bad.size());

            if (brk_len == box_dims(brk_dir) ||
                bad[brk_len] == false) {
               // Cutting brk_len from low side is ok.
               best_difference = difference;
               best_breakoff_box = box;
               best_breakoff_box.upper() (brk_dir) =
                  best_breakoff_box.lower() (brk_dir) + brk_len - 1;
               best_leftover_box = box;
               best_leftover_box.lower() (brk_dir) =
                  best_breakoff_box.upper() (brk_dir) + 1;
               break;
            } else if (bad[box_dims(brk_dir) - brk_len] == false) {
               // Cutting brk_len from high side is ok.
               best_difference = difference;
               best_breakoff_box = box;
               best_breakoff_box.lower() (brk_dir) =
                  best_breakoff_box.upper() (brk_dir) - brk_len + 1;
               best_leftover_box = box;
               best_leftover_box.upper() (brk_dir) =
                  best_breakoff_box.lower() (brk_dir) - 1;
               break;
            }
         }
      }

   }

   bool successful_break = false;

   if (!best_breakoff_box.empty()) {
      breakoff.push_back(best_breakoff_box);
      brk_load = best_breakoff_box.size();
      successful_break = true;
      if (d_print_break_steps) {
         tbox::plog << "      breakOffload_planar broke off box " << box
                    << " for breakoff box " << best_breakoff_box
                    << " and leftover " << best_leftover_box << std::endl;
      }
   } else {
      if (d_print_break_steps) {
         tbox::plog << "      breakOffload_planar could not break "
                    << ideal_load_to_break << " from Box " << box
                    << std::endl;
      }
   }
   if (!best_leftover_box.empty()) {
      leftover.push_back(best_leftover_box);
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   for (std::vector<hier::Box>::iterator bi = breakoff.begin();
        bi != breakoff.end();
        ++bi) {
      const hier::Box& b = *bi;
      const hier::IntVector s = b.numberCells();
      for (int d = 0; d < d_dim.getValue(); ++d) {
         if (((s(d) < d_min_size(d)) && (s(d) != box_dims(d))) ||
             (s(d) > box_dims(d))) {
            TBOX_ERROR("TreeLoadBalancerOld library error:\n"
               << "breakoff box " << b << ", size " << s
               << "\nis not between the min size " << d_min_size
               << "\nand the original box size " << box_dims << "\n"
               << "break box size " << best_breakoff_box.numberCells() << "\n"
               << "ideal brk load " << ideal_load_to_break);
         }
      }
   }
   for (std::vector<hier::Box>::iterator bi = leftover.begin();
        bi != leftover.end();
        ++bi) {
      const hier::Box& b = *bi;
      const hier::IntVector s = b.numberCells();
      for (int d = 0; d < d_dim.getValue(); ++d) {
         if (((s(d) < d_min_size(d)) && (s(d) != box_dims(d))) ||
             (s(d) > box_dims(d))) {
            TBOX_ERROR("TreeLoadBalancerOld library error:\n"
               << "leftover box " << b << ", size " << s
               << "\nis not between the min size " << d_min_size
               << "\nand the original box size " << box_dims << "\n"
               << "break box size " << best_breakoff_box.numberCells() << "\n"
               << "ideal brk load " << ideal_load_to_break);
         }
      }
   }
#endif

   return successful_break;
}



/*
 *************************************************************************
 * Break off a load from a box by making a box cut that is as close
 * to cubic as feasible.
 *
 * Currently assuming uniform load of one unit per cell.
 *
 * Return whether a successful break was made.
 *
 * Differs from breakOffLoad in that it will always
 * performs a break and if needed, break off more than
 * the ideal.  The calling method should take this into
 * account.
 *************************************************************************
 */
bool
TreeLoadBalancerOld::breakOffLoad_cubic(
   std::vector<hier::Box>& breakoff,
   std::vector<hier::Box>& leftover,
   double& brk_load,
   const hier::Box& box,
   double ideal_load_to_break,
   const tbox::Array<tbox::Array<bool> >& bad_cuts ) const
{

   const hier::IntVector box_dims(box.numberCells());

   const double box_load(box_dims.getProduct());

   if (ideal_load_to_break >= box_load) {
      // Easy: break off everything.
      leftover.clear();
      breakoff.clear();
      breakoff.push_back(box);
      brk_load = box_load;
      if (d_print_break_steps) {
         tbox::plog << "      breakOffload_cubic broke off entire Box "
                    << box
                    << std::endl;
      }
      return true;
   }

   if (ideal_load_to_break > 0.5 * box_load) {
      /*
       * This algorithm is better when breaking off a small portion.
       * Since the ideal is a bigger portion, switch breakoff with leftover.
       */
      if (d_print_break_steps) {
         tbox::plog
         << "      breakOffload_cubic reversing direction to break "
         << (box_dims.getProduct() - ideal_load_to_break)
         << " instead of " << ideal_load_to_break << " / "
         << box_dims.getProduct() << std::endl;
      }
      bool success =
         breakOffLoad_cubic(
            leftover,
            breakoff,
            brk_load,
            box,
            box_dims.getProduct() - ideal_load_to_break,
            bad_cuts );
      if (success) {
         brk_load = box_dims.getProduct() - brk_load;
      }
      return success;
   }

   if (d_print_break_steps) {
      tbox::plog << "      breakOffload_cubic attempting to break "
                 << ideal_load_to_break << " from Box " << box
                 << " min_size=" << d_min_size << std::endl;
   }

   breakoff.clear();
   leftover.clear();

   /*
    * brk_size is the size of the box we want to break off of
    * box.  We start with the smallest allowed brk_size that
    * will not create remainders that violate size constraints.
    *
    * In the do loop below, we increase brk_size to bring brk_load
    * closer to ideal_oad_to_break.
    */
   hier::IntVector brk_size(d_min_size);
   brk_size.max(d_cut_factor);
   brk_size.min(box_dims);
   /*
    * If remainder is too small, zero it out to avoid
    * having non-zero remainder smaller than d_min_size.
    */
   for (int d = 0; d < d_dim.getValue(); ++d) {
      if ((box_dims(d) - brk_size(d) > 0) &&
          (box_dims(d) - brk_size(d) < d_min_size(d))) {
         brk_size(d) = box_dims(d);
      }
   }
   brk_load = brk_size.getProduct();

   if (d_print_break_steps) {
      tbox::plog << "      brk: " << std::flush;
      tbox::plog << "  " << brk_size << "->" << brk_load << std::flush;
   }

   /*
    * stop_growing: whether dimensions of brk_size is already
    * big engough so that it cannot not be grown without breaking
    * off too much.
    */
   hier::IntVector stop_growing(d_dim, 0);
   for (int d = 0; d < d_dim.getValue(); ++d) {
      if (brk_size[d] == box_dims[d]) stop_growing[d] = 1;
   }

   if (brk_load < ideal_load_to_break) {
      /*
       * The do loop gradually increases brk_size to bring brk_load
       * closer to ideal_load_to_break.
       *
       * Select dimension to increase in size, inc_dim.  Use the
       * smallest dimension that is still allowed to grow.
       *
       * Try a new break size that is bigger than the current brk_size
       * by the minimal allowed amount.  If this brings us closer to
       * the ideal break amount, mark it.
       *
       * Exit the loop when we cannot grow anymore or we are already
       * breaking off more than the ideal.
       */
      do {

         int inc_dim = -1;
         for (int d = 0; d < d_dim.getValue(); ++d) {
            if (!stop_growing(d) &&
                (inc_dim == -1 || brk_size(inc_dim) > brk_size(d))) inc_dim = d;
         }
         if (inc_dim == -1) break; // No growable dimension.

         TBOX_ASSERT(brk_size(inc_dim) < box_dims(inc_dim));

         hier::IntVector new_brk_size(brk_size);
         new_brk_size(inc_dim) += d_cut_factor(inc_dim);
         if (d_print_break_steps) {
            tbox::plog << "  " << brk_size << "=>" << brk_load << std::flush;
         }

         // Prevent remainder being smaller than d_min_size.
         if (box_dims(inc_dim) - new_brk_size(inc_dim) < d_min_size(inc_dim)) {
            new_brk_size(inc_dim) = box_dims(inc_dim);
            if (d_print_break_steps) {
               tbox::plog << "  " << brk_size << "==>" << brk_load
                          << std::flush;
            }
         }

         if (new_brk_size(inc_dim) == box_dims(inc_dim)) {
            stop_growing(inc_dim) = 1;
         }

         int new_brk_load = new_brk_size.getProduct();

         if (new_brk_load <= ideal_load_to_break) {
            /*
             * new_brk_load is closer to ideal than current brk_load.
             * Don't break out of the loop yet.  We will grow it again
             * and check.
             */
            brk_size = new_brk_size;
            brk_load = new_brk_load;
            if (d_print_break_steps) {
               tbox::plog << "  " << brk_size << "===>" << brk_load
                          << std::flush;
            }
         } else if ((new_brk_load - ideal_load_to_break) <
                    static_cast<double>(ideal_load_to_break - brk_load)) {
            /*
             * new_brk_load is bigger than ideal but is still an
             * improvement over brk_load.  Accept it, but break out of
             * the loop because further growing will not improve
             * result.
             */
            brk_size = new_brk_size;
            brk_load = new_brk_load;
            if (d_print_break_steps) {
               tbox::plog << "  " << brk_size << "====>" << brk_load
                          << std::flush;
            }
            break;
         } else {
            /*
             * Even though dimension inc_dim has not reached the box
             * dimension, stop growing it because any further growth
             * leads to too big a load.
             */
            stop_growing(inc_dim) = 1;
         }

      } while (true);
   }

   if (d_print_break_steps) {
      tbox::plog << std::endl;
   }

   /*
    * Find a place to put the break-off box so that it does not lie
    * across a bad cut.  If no such placement is found, set
    * placement_impossible to true.
    */
   hier::Box breakoff_box(d_dim);
   breakoff_box.setBlockId(box.getBlockId());
   const hier::IntVector& lower(box.lower());
   const hier::IntVector& upper(box.upper());
   bool placement_impossible = false;
   if (d_print_break_steps) {
      tbox::plog << "      Placing " << brk_size
                 << " to avoid bad cut points" << std::flush;
   }
   for (int d = 0; d < d_dim.getValue(); ++d) {
      /*
       * To minimize the number of boxes remaining after breakoff_box
       * is removed, prefer to place breakoff_box against the upper or
       * lower side of the Box.  First try putting the breakoff_box
       * along the upper side of dimension d.  If it cuts the box at a
       * bad point, try the lower side.  If it still cuts at a bad
       * points, slide the box toward the upper side until it does not
       * cut at any bad points, being careful not to be so close to
       * the box boundaries that we generate remainder boxes violating
       * the minimum size.  If no place can be found to put
       * breakoff_box, set placement_impossible and give up.  (We
       * could go back and reshape the box at this point, but we won't
       * because there is probably another box that would work without
       * reshaping.)
       */
      const tbox::Array<bool>& bad = bad_cuts[d];

      int& gl = breakoff_box.lower()[d];
      int& gu = breakoff_box.upper()[d];

      gu = upper[d];
      gl = gu - (brk_size[d] - 1);
      if (!bad[gl - lower[d]]) {
         if (d_print_break_steps) {
            tbox::plog << "  d=" << d << " upper" << std::flush;
         }
         continue;
      }

      gl = lower[d];
      gu = gl + (brk_size[d] - 1);
      if (gu + 1 - lower[d] < bad.size() && !bad[gu + 1 - lower[d]]) {
         if (d_print_break_steps) {
            tbox::plog << "  d=" << d << " lower" << std::flush;
         }
         continue;
      }

      gl = lower[d] + d_min_size[d];
      gu = gl + (brk_size[d] - 1);
      while (gu <= upper[d] - d_min_size[d] &&
             (bad[gl - lower[d]] || bad[gu + 1 - lower[d]])) {
         ++gl;
         ++gu;
      }
      if (gu <= upper[d] - d_min_size[d]) {
         if (d_print_break_steps) {
            tbox::plog << "  d=" << d << " middle" << std::flush;
         }
         continue;
      }

      if (d_print_break_steps) {
         tbox::plog << "  cannot place dim " << d
                    << " without creating bad cuts."
                    << std::flush;
      }
      placement_impossible = true;
      /*
       * Cannot find place for breakoff_box along dimension d.
       * No point in looking at higher dimensions.
       */
      break;
   }
   if (d_print_break_steps) {
      tbox::plog << std::endl;
   }
   if (placement_impossible) {
      return false;
   }

   breakoff.clear();
   breakoff.push_back(breakoff_box);

   burstBox(
      leftover,
      box,
      breakoff_box );

#ifdef DEBUG_CHECK_ASSERTIONS
   for (std::vector<hier::Box>::iterator bi = breakoff.begin();
        bi != breakoff.end();
        ++bi) {
      const hier::Box& b = *bi;
      const hier::IntVector s = b.numberCells();
      for (int d = 0; d < d_dim.getValue(); ++d) {
         if (((s(d) < d_min_size(d)) && (s(d) != box_dims(d))) ||
             (s(d) > box_dims(d))) {
            TBOX_ERROR("TreeLoadBalancerOld library error:\n"
               << "breakoff box " << b << ", with size " << s
               << "\nis not between the min size " << d_min_size
               << "\nand the original box size " << box_dims << "\n"
               << "orig box " << box << "\n"
               << "break box " << breakoff_box << "\n"
               << "break box size " << brk_size << "\n"
               << "ideal brk load " << ideal_load_to_break);
         }
      }
   }
   for (std::vector<hier::Box>::iterator bi = leftover.begin();
        bi != leftover.end();
        ++bi) {
      const hier::Box& b = *bi;
      const hier::IntVector s = b.numberCells();
      for (int d = 0; d < d_dim.getValue(); ++d) {
         if (((s(d) < d_min_size(d)) && (s(d) != box_dims(d))) ||
             (s(d) > box_dims(d))) {
            TBOX_ERROR("TreeLoadBalancerOld library error:\n"
               << "leftover box " << b << ", with size " << s
               << "\nis not between the min size " << d_min_size
               << "\nand the original box size " << box_dims << "\n"
               << "orig box " << box << "\n"
               << "break box " << breakoff_box << "\n"
               << "break box size " << brk_size << "\n"
               << "ideal brk load " << ideal_load_to_break);
         }
      }
   }
#endif

   return true;
}



/*
**************************************************************************
* Move Boxes in balance_box_level from ranks outside of
* rank_group to ranks inside rank_group.  Modify the given connectors
* to make them correct following this moving of boxes.
**************************************************************************
*/

void
TreeLoadBalancerOld::prebalanceBoxLevel(
   hier::BoxLevel& balance_box_level,
   hier::Connector& balance_to_anchor,
   hier::Connector& anchor_to_balance,
   const tbox::RankGroup& rank_group) const
{

   if (balance_to_anchor.isFinalized()) {
      TBOX_ASSERT(anchor_to_balance.checkTransposeCorrectness(balance_to_anchor) == 0);
      TBOX_ASSERT(balance_to_anchor.checkTransposeCorrectness(anchor_to_balance) == 0);
   }

   /*
    * tmp_box_level will contain the same boxes as
    * balance_box_level, but all will live on the processors
    * specified in rank_group.
    */
   hier::BoxLevel tmp_box_level(d_dim);
   tmp_box_level.initialize(
      balance_box_level.getRefinementRatio(),
      balance_box_level.getGridGeometry(),
      balance_box_level.getMPI());

   /*
    * If a rank is not in rank_group it is called a "sending" rank, as
    * it will send any Boxes it has to a rank in rank_group.
    */
   bool is_sending_rank = rank_group.isMember(d_mpi.getRank()) ? false : true;

   int output_nproc = rank_group.size();

   /*
    * the send and receive comm objects
    */
   tbox::AsyncCommStage comm_stage;
   tbox::AsyncCommPeer<int>* box_send = NULL;
   tbox::AsyncCommPeer<int>* box_recv = NULL;
   tbox::AsyncCommPeer<int>* id_send = NULL;
   tbox::AsyncCommPeer<int>* id_recv = NULL;

   /*
    * A sending rank will send its Boxes to a receiving rank, and
    * that receiving processor will add it to its local set of Boxes.
    * When the box is added on the receiving processor, it will receive
    * a new LocalId.  This LocalId value needs to be sent back to
    * the sending processor, in order to construct the mapping connectors.
    *
    * Therefore the sending ranks construct comm objects for sending boxes
    * and receiving LocalIdes.
    *
    * Sending processors send to ranks in the rank_group determined by
    * a modulo heuristic.
    */
   if (is_sending_rank) {
      box_send = new tbox::AsyncCommPeer<int>;
      box_send->initialize(&comm_stage);
      box_send->setPeerRank(rank_group.getMappedRank(d_mpi.getRank() % output_nproc));
      box_send->setMPI(d_mpi);
      box_send->setMPITag(TreeLoadBalancerOld_PREBALANCE0 + 2 * d_mpi.getRank(),
         TreeLoadBalancerOld_PREBALANCE1 + 2 * d_mpi.getRank());

      id_recv = new tbox::AsyncCommPeer<int>;
      id_recv->initialize(&comm_stage);
      id_recv->setPeerRank(rank_group.getMappedRank(d_mpi.getRank() % output_nproc));
      id_recv->setMPI(d_mpi);
      id_recv->setMPITag(TreeLoadBalancerOld_PREBALANCE0 + 2 * d_mpi.getRank(),
         TreeLoadBalancerOld_PREBALANCE1 + 2 * d_mpi.getRank());
   }

   /*
    * The receiving ranks construct comm objects for receiving boxes
    * and sending LocalIdes.
    */
   int num_recvs = 0;
   if (rank_group.isMember(d_mpi.getRank())) {
      std::list<int> recv_ranks;
      for (int i = 0; i < d_mpi.getSize(); i++) {
         if (!rank_group.isMember(i) &&
             rank_group.getMappedRank(i % output_nproc) == d_mpi.getRank()) {
            recv_ranks.push_back(i);
         }
      }
      num_recvs = static_cast<int>(recv_ranks.size());
      if (num_recvs > 0) {
         box_recv = new tbox::AsyncCommPeer<int>[num_recvs];
         id_send = new tbox::AsyncCommPeer<int>[num_recvs];
         int recv_count = 0;
         for (std::list<int>::const_iterator ri(recv_ranks.begin());
              ri != recv_ranks.end(); ri++) {
            const int rank = *ri;
            box_recv[recv_count].initialize(&comm_stage);
            box_recv[recv_count].setPeerRank(rank);
            box_recv[recv_count].setMPI(d_mpi);
            box_recv[recv_count].setMPITag(TreeLoadBalancerOld_PREBALANCE0 + 2 * rank,
               TreeLoadBalancerOld_PREBALANCE1 + 2 * rank);

            id_send[recv_count].initialize(&comm_stage);
            id_send[recv_count].setPeerRank(rank);
            id_send[recv_count].setMPI(d_mpi);
            id_send[recv_count].setMPITag(TreeLoadBalancerOld_PREBALANCE0 + 2 * rank,
               TreeLoadBalancerOld_PREBALANCE1 + 2 * rank);

            recv_count++;
         }
         TBOX_ASSERT(num_recvs == recv_count);
      }
   }

   /*
    * Construct the mapping Connectors which describe the mapping from the box
    * configuration of the given balance_box_level, to the new
    * configuration stored in tmp_box_level.  These mapping Connectors
    * are necessary to modify the two Connectors given in the argument list,
    * so that on return from this method, they will be correct for the new
    * balance_box_level.
    */
   hier::Connector balance_to_tmp(balance_box_level,
                                  tmp_box_level,
                                  hier::IntVector::getZero(d_dim));

   hier::Connector tmp_to_balance(tmp_box_level,
                                  balance_box_level,
                                  hier::IntVector::getZero(d_dim));

   /*
    * Where Boxes already exist on ranks in rank_group,
    * move them directly to tmp_box_level.
    */
   if (!is_sending_rank) {
      const hier::BoxContainer& unchanged_boxes =
         balance_box_level.getBoxes();

      for (hier::BoxContainer::const_iterator ni =
              unchanged_boxes.begin();
           ni != unchanged_boxes.end(); ++ni) {

         const hier::Box& box = *ni;
         tmp_box_level.addBox(box);
      }
   }

   const int buf_size = hier::Box::commBufferSize(d_dim);

   /*
    * On sending ranks, pack the Boxes into buffers and send.
    */
   if (is_sending_rank) {
      const hier::BoxContainer& sending_boxes =
         balance_box_level.getBoxes();
      const int num_sending_boxes =
         static_cast<int>(sending_boxes.size());

      int* buffer = new int[buf_size * num_sending_boxes];
      int box_count = 0;
      for (hier::BoxContainer::const_iterator ni = sending_boxes.begin();
           ni != sending_boxes.end(); ++ni) {

         const hier::Box& box = *ni;

         box.putToIntBuffer(&buffer[box_count * buf_size]);
         box_count++;
      }
      box_send->beginSend(buffer, buf_size * num_sending_boxes);

      delete[] buffer;
   }

   /*
    * On receiving ranks, complete the receives, add the boxes to local
    * tmp_box_level, insert boxes into tmp_to_balance, and then
    * send the new LocalIdes back to the sending processors.
    */
   if (!is_sending_rank && num_recvs > 0) {
      for (int i = 0; i < num_recvs; i++) {
         box_recv[i].beginRecv();
      }
      int num_completed_recvs = 0;
      tbox::Array<bool> completed(num_recvs, false);
      while (num_completed_recvs < num_recvs) {
         for (int i = 0; i < num_recvs; i++) {
            if (!completed[i] && box_recv[i].checkRecv()) {
               num_completed_recvs++;
               completed[i] = true;
               const int num_boxes = box_recv[i].getRecvSize() / buf_size;
               const int* buffer = box_recv[i].getRecvData();
               int* id_buffer = new int[num_boxes];

               for (int b = 0; b < num_boxes; b++) {
                  hier::Box box(d_dim);

                  box.getFromIntBuffer(&buffer[b * buf_size]);

                  hier::BoxContainer::const_iterator tmp_iter =
                     tmp_box_level.addBox(box,
                        box.getBlockId());

                  hier::BoxId tmp_box_id = tmp_iter->getId();

                  tmp_to_balance.insertLocalNeighbor(box, tmp_box_id);

                  id_buffer[b] = tmp_box_id.getLocalId().getValue();
               }
               id_send[i].beginSend(id_buffer, num_boxes);

               delete[] id_buffer;
            }
         }
      }
      for (int i = 0; i < num_recvs; i++) {
         if (!id_send[i].checkSend()) {
            id_send[i].completeCurrentOperation();
         }
      }
   }

   /*
    * On sending ranks, receive the LocalIds, and add the edges
    * to balance_to_tmp.
    */
   if (is_sending_rank) {
      if (!box_send->checkSend()) {
         box_send->completeCurrentOperation();
      }

      id_recv->beginRecv();

      if (!id_recv->checkRecv()) {
         id_recv->completeCurrentOperation();
      }
      const int* buffer = id_recv->getRecvData();

      const hier::BoxContainer& sending_boxes =
         balance_box_level.getBoxes();
      TBOX_ASSERT(static_cast<int>(id_recv->getRecvSize()) == sending_boxes.size());

      int box_count = 0;
      for (hier::BoxContainer::const_iterator ni = sending_boxes.begin();
           ni != sending_boxes.end(); ++ni) {

         hier::Box new_box(
            *ni,
            (hier::LocalId)buffer[box_count],
            rank_group.getMappedRank(d_mpi.getRank() % output_nproc));

         balance_to_tmp.insertLocalNeighbor(new_box, (*ni).getId());
         box_count++;
      }
   }
   balance_to_tmp.setConnectorType(hier::Connector::MAPPING);
   tmp_to_balance.setConnectorType(hier::Connector::MAPPING);

   if (anchor_to_balance.isFinalized()) {
      /*
       * This modify operation copies tmp_box_level to
       * balance_box_level, and changes anchor_to_balance and
       * balance_to_anchor such that they are correct for the new state
       * of balance_box_level.
       */
      const hier::MappingConnectorAlgorithm mca;
      mca.modify(anchor_to_balance,
         balance_to_anchor,
         balance_to_tmp,
         tmp_to_balance,
         &balance_box_level,
         &tmp_box_level);

      TBOX_ASSERT(anchor_to_balance.checkTransposeCorrectness(balance_to_anchor) == 0);
      TBOX_ASSERT(balance_to_anchor.checkTransposeCorrectness(anchor_to_balance) == 0);
   } else {
      hier::BoxLevel::swap(balance_box_level, tmp_box_level);
   }

   /*
    * Clean up raw pointer allocation.
    */
   if (is_sending_rank) {
      delete box_send;
      delete id_recv;
   }
   if (num_recvs) {
      delete[] box_recv;
      delete[] id_send;
   }
}



/*
**************************************************************************
**************************************************************************
*/

std::ostream&
operator << (
   std::ostream& co,
   const TreeLoadBalancerOld::BoxInTransit& r)
{
   co << r.d_box
   << r.d_box.numberCells() << '|'
   << r.d_box.numberCells().getProduct();
   for (int i = static_cast<int>(r.d_proc_hist.size()) - 1; i >= 0; --i) {
      co << '-' << r.d_proc_hist[i];
   }
   co << '-' << r.d_orig_box
   << r.d_box.numberCells() << '|'
   << r.d_box.numberCells().getProduct();
   return co;
}



/*
 ***********************************************************************
 ***********************************************************************
 */
void
TreeLoadBalancerOld::setTimers()
{
   /*
    * The first constructor gets timers from the TimerManager.
    * and sets up their deallocation.
    */
   if (!t_load_balance_box_level) {
      t_load_balance_box_level = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::loadBalanceBoxLevel()");
      t_get_map = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::get_map");
      t_use_map = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::use_map");
      t_constrain_size = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::constrain_size");
      t_map_big_boxes = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::mapOversizedBoxes()");
      t_load_distribution = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::load_distribution");
      t_compute_local_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::compute_local_load");
      t_compute_global_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::compute_global_load");
      t_compute_tree_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::compute_tree_load");

      const int max_cycles_to_time = 4;
      t_compute_tree_load_for_cycle.resize(
         max_cycles_to_time,
         boost::shared_ptr<tbox::Timer>() );
      for ( int i=0; i<max_cycles_to_time; ++i ) {
         t_compute_tree_load_for_cycle[i] = tbox::TimerManager::getManager()->
            getTimer(d_object_name + "::compute_tree_load_for_cycle["
                     + tbox::Utilities::intToString(i) + "]");
      }

      t_break_off_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::breakOffLoad()");
      t_find_bad_cuts = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::find_bad_cuts");
      t_reassign_loads = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::reassignLoads()");
      t_shift_loads_by_swapping = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::shiftLoadsBySwapping()");
      t_shift_loads_by_breaking = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::shiftLoadsByBreaking()");
      t_find_swap_pair = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::swapLoadPair()");
      t_send_load_to_children = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::send_load_to_children");
      t_send_load_to_parent = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::send_load_to_parent");
      t_get_load_from_children = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::get_load_from_children");
      t_get_load_from_parent = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::get_load_from_parent");
      t_construct_semilocal = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::constructSemilocalUnbalancedToBalanced()");
      t_construct_semilocal_comm_wait = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::constructSemilocalUnbalancedToBalanced()_comm_wait");
      t_send_edge_to_children = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::send_edge_to_children");
      t_send_edge_to_parent = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::send_edge_to_parent");
      t_get_edge_from_children = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::get_edge_from_children");
      t_get_edge_from_parent = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::get_edge_from_parent");
      t_report_loads = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::report_loads");
      t_finish_sends = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::finish_sends");
      t_local_balancing = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::local_balancing");
      t_pack_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::pack_load");
      t_unpack_load = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::unpack_load");
      t_unpack_edge = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::unpack_edge");
      t_parent_load_comm = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::parent_load_comm");
      t_children_load_comm = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::children_load_comm");
      t_parent_edge_comm = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::parent_edge_comm");
      t_children_edge_comm = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::children_edge_comm");
      t_barrier_before = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::barrier_before");
      t_barrier_after = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::barrier_after");
      t_child_send_wait = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::child_send_wait");
      t_child_recv_wait = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::child_recv_wait");
      t_parent_send_wait = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::parent_send_wait");
      t_parent_recv_wait = tbox::TimerManager::getManager()->
         getTimer(d_object_name + "::parent_recv_wait");
   }
}



/*
 *************************************************************************
 *************************************************************************
 */
TreeLoadBalancerOld::BoxInTransit::BoxInTransit(
   const tbox::Dimension &dim) :
   d_box(dim),
   d_orig_box(dim)
{
}



/*
 *************************************************************************
 *************************************************************************
 */
TreeLoadBalancerOld::BoxInTransit::BoxInTransit(
   const hier::Box& origin):
   d_box(origin),
   d_orig_box(origin),
   d_load(origin.size()),
   d_proc_hist()
{
}



/*
 *************************************************************************
 *************************************************************************
 */
TreeLoadBalancerOld::BoxInTransit::BoxInTransit(
   const int *&ptr,
   const tbox::Dimension &dim ):
   d_box(dim),
   d_orig_box(dim),
   d_load(0),
   d_proc_hist(0)
{
   getFromIntBuffer(ptr);
   ptr += commBufferSize();
}



/*
 *************************************************************************
 * Construct a new BoxInTransit with the history of an existing box.
 *************************************************************************
 */
TreeLoadBalancerOld::BoxInTransit::BoxInTransit(
   const BoxInTransit& other,
   const hier::Box& box,
   int rank,
   hier::LocalId local_id):
   d_box(box, local_id, rank),
   d_orig_box(other.d_orig_box),
   d_load(d_box.size()),
   d_proc_hist(other.d_proc_hist)
{
   if (rank != other.getOwnerRank()) {
      d_proc_hist.push_back(other.getOwnerRank());
   }
}

int*
TreeLoadBalancerOld::BoxInTransit::putToIntBuffer(
   int* buffer,
   bool skip_last_owner) const
{
   const int box_comm_buf_size = hier::Box::commBufferSize(d_box.getDim());
   d_box.putToIntBuffer(buffer);
   buffer += box_comm_buf_size;
   d_orig_box.putToIntBuffer(buffer);
   buffer += box_comm_buf_size;
   *(buffer++) = d_load;
   if (skip_last_owner) {
      TBOX_ASSERT( !d_proc_hist.empty() );
   }
   *(buffer++) = static_cast<int>(d_proc_hist.size()-skip_last_owner);
   for (unsigned int i = 0; i < d_proc_hist.size()-skip_last_owner; ++i) {
      *(buffer++) = d_proc_hist[i];
   }
   return buffer;
}

const int*
TreeLoadBalancerOld::BoxInTransit::getFromIntBuffer(
   const int* buffer)
{
   const int box_comm_buf_size = hier::Box::commBufferSize(d_box.getDim());
   d_box.getFromIntBuffer(buffer);
   buffer += box_comm_buf_size;
   d_orig_box.getFromIntBuffer(buffer);
   buffer += box_comm_buf_size;
   d_load = *(buffer++);
   d_proc_hist.clear();
   d_proc_hist.insert(d_proc_hist.end(), *(buffer++), 0);
   for (unsigned int i = 0; i < d_proc_hist.size(); ++i) {
      d_proc_hist[i] = *(buffer++);
   }
   return buffer;
}

void
TreeLoadBalancerOld::BoxInTransit::packForPreviousOwner(
   std::map<int,std::vector<int> > &outgoing_messages ) const
{
   const int prev_owner = d_proc_hist.back();
   std::vector<int> &msg(outgoing_messages[prev_owner]);
   const int cbs1 = commBufferSize()-1;
   msg.insert(msg.end(), cbs1, 0);
   putToIntBuffer(&msg[msg.size() - cbs1], true);
}

TreeLoadBalancerOld::SubtreeLoadData::SubtreeLoadData():
   d_num_procs(0),
   d_total_work(0),
   d_load_exported(0),
   d_load_imported(0)
{
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
