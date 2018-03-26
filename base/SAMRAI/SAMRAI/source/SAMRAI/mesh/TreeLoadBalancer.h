/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_TreeLoadBalancer
#define included_mesh_TreeLoadBalancer

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/mesh/BalanceUtilities.h"
#include "SAMRAI/mesh/LoadBalanceStrategy.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/BalancedDepthFirstTree.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/RankGroup.h"
#include "SAMRAI/tbox/Statistic.h"
#include "SAMRAI/tbox/Statistician.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Utilities.h"

#include <boost/shared_ptr.hpp>
#include <iostream>
#include <vector>
#include <set>

namespace SAMRAI {
namespace mesh {





/*!
 * @brief Provides load balancing routines for AMR hierarchy by
 * implemementing the LoadBalancerStrategy.
 *
 * This class implements a tree-based load balancer.  The MPI
 * processes are arranged in a tree.  Work load is transmitted from
 * process to process along the edges of the tree.
 *
 * Currently, only uniform load balancing is supported.  Eventually,
 * non-uniform load balancing should be supported.  (Non-uniform load
 * balancing is supported by the CutAndPackLoadBalancer class.)
 *
 * Inputs and their default values:
 *
 * No special inputs are required for this class.
 *
 * @verbatim
 * min_load_fraction_per_box = 0.03
 *                           // Additional restriction on box size.
 *                           // Will not generate a box that has less
 *                           // than this fraction of the global average
 *                           // work load in order to move work load.
 *                           // Set to negative to disable.
 * report_load_balance = TRUE // Write out load balance report in log
 * n_root_cycles = -1         // Number of steps over which to smoothly spread
 *                            // out work load.  This helps scalability when
 *                            // initial work load is grossly unbalanced.
 *                            // Usually 1 step is sufficient.
 *                            // Can be set higher (2 or 3) to reduce negative
 *                            // performance effects of extremely poor initial
 *                            // load balance.  Set to -1 (default) to compute
 *                            // the number of cycles by a simple heuristic.
 *                            // Set to zero to effectively bypass load balancing.
 * @endverbatim
 *
 * @see mesh::LoadBalanceStrategy
 */

class TreeLoadBalancer:
   public LoadBalanceStrategy
{
public:
   /*!
    * @brief Initializing constructor sets object state to default or,
    * if database provided, to parameters in database.
    *
    * @param[in] dim
    *
    * @param[in] name User-defined std::string identifier used for error
    * reporting and timer names.  If omitted, "TreeLoadBalancer"
    * is used.
    *
    * @param[in] input_db (optional) database pointer providing
    * parameters from input file.  This pointer may be null indicating
    * no input is used.
    */
   TreeLoadBalancer(
      const tbox::Dimension& dim,
      const std::string& name = std::string("TreeLoadBalancer"),
      const boost::shared_ptr<tbox::Database>& input_db =
         boost::shared_ptr<tbox::Database>());

   /*!
    * @brief Virtual destructor releases all internal storage.
    */
   virtual ~TreeLoadBalancer();

   /*!
    * @brief Set the internal SAMRAI_MPI to a duplicate of the given
    * SAMRAI_MPI.
    *
    * The given SAMRAI_MPI must have a valid communicator.
    *
    * The given SAMRAI_MPI is duplicated for private use.  This
    * requires a global communication, so all processes in the
    * communicator must call it.  The advantage of a duplicate
    * communicator is that it ensures the communications for the
    * object won't accidentally interact with other communications.
    *
    * If the duplicate SAMRAI_MPI it is set, the TreeLoadBalancer will
    * only balance BoxLevels with congruent SAMRAI_MPI objects and
    * will use the duplicate SAMRAI_MPI for communications.
    * Otherwise, the SAMRAI_MPI of the BoxLevel will be used.  The
    * duplicate MPI communicator is freed when the object is
    * destructed, or freeMPICommunicator() is called.
    */
   void
   setSAMRAI_MPI(
      const tbox::SAMRAI_MPI& samrai_mpi);

   /*!
    * @brief Free the internal MPI communicator, if any has been set.
    *
    * This is automatically done by the destructor, if needed.
    *
    * @see setSAMRAI_MPI().
    */
   void
   freeMPICommunicator();

   /*!
    * @brief Configure the load balancer to use the data stored
    * in the hierarchy at the specified descriptor index
    * for estimating the workload on each cell.
    *
    * Note: This method currently does not affect the results because
    * this class does not yet support uniform load balancing.
    *
    * @param data_id
    * Integer value of patch data identifier for workload
    * estimate on each cell.  An invalid value (i.e., < 0)
    * indicates that a spatially-uniform work estimate
    * will be used.  The default value is -1 (undefined)
    * implying the uniform work estimate.
    *
    * @param level_number
    * Optional integer number for level on which data id
    * is used.  If no value is given, the data will be
    * used for all levels.
    */
   void
   setWorkloadPatchDataIndex(
      int data_id,
      int level_number = -1);

   /*!
    * @brief Configure the load balancer to load balance boxes by
    * assuming all cells on the specified level or all hierarchy
    * levels are weighted equally.
    *
    * @param level_number
    * Optional integer number for level on which uniform
    * workload estimate will be used.  If the level
    * number is not specified, a uniform workload
    * estimate will be used on all levels.
    */
   void
   setUniformWorkload(
      int level_number = -1)
   {
      d_workload_data_id[level_number] = -1;
   }

   /*!
    * @brief Return true if load balancing procedure for given level
    * depends on patch data on mesh; otherwise return false.
    *
    * @param[in] level_number  Integer patch level number.
    */
   bool
   getLoadBalanceDependsOnPatchData(
      int level_number) const;

   /*!
    * @copydoc LoadBalanceStrategy::loadBalanceBoxLevel()
    *
    * Note: This implementation does not yet support non-uniform load
    * balancing.
    */
   void
   loadBalanceBoxLevel(
      hier::BoxLevel& balance_mapped_box_level,
      hier::Connector& balance_to_anchor,
      hier::Connector& anchor_to_balance,
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const hier::Connector& unbalanced_to_attractor,
      const hier::Connector& attractor_to_unbalanced,
      const hier::IntVector& min_size,
      const hier::IntVector& max_size,
      const hier::BoxLevel& domain_mapped_box_level,
      const hier::IntVector& bad_interval,
      const hier::IntVector& cut_factor,
      const tbox::RankGroup& rank_group = tbox::RankGroup()) const;

   /*!
    * @brief Print out all members of the class instance to given
    * output stream.
    *
    * @param[in] output_stream
    */
   virtual void
   printClassData(
      std::ostream& output_stream) const;

   /*!
    * @brief Write out statistics recorded for the most recent load
    * balancing result.
    *
    * @param[in] output_stream
    */
   void
   printStatistics(
      std::ostream& output_stream = tbox::plog) const
   {
      BalanceUtilities::gatherAndReportLoadBalance(d_load_stat,
         tbox::SAMRAI_MPI::getSAMRAIWorld(),
         output_stream);
   }



private:

   typedef int LoadType;

   /*!
    * @brief Data to save for each Box that gets passed along the tree
    * edges.
    *
    * The purpose of the BoxInTransit is to associate extra data with
    * a Box as it is broken up and passed from process to process.  A
    * BoxInTransit is a Box going through these changes.  It has a
    * current work load and an orginating Box.
    */
   struct BoxInTransit {

      /*!
       * @brief Constructor
       *
       * @param[in] dim
       */
      BoxInTransit(const tbox::Dimension& dim);

      /*!
       * @brief Construct a new BoxInTransit from an originating box.
       *
       * @param[in] other
       */
      BoxInTransit(const hier::Box& origin);

      /*!
       * @brief Construct new object having the history an existing
       * object but is otherwise different.
       *
       * @param[in] other
       *
       * @param[in] box
       *
       * @param[in] rank
       *
       * @param[in] local_id
       */
      BoxInTransit(
         const BoxInTransit& other,
         const hier::Box& box,
         int rank,
         hier::LocalId local_id);

      /*!
       * @brief Assignment operator
       *
       * @param[in] other
       */
      const BoxInTransit&
      operator = (const BoxInTransit& other)
      {
         d_box = other.d_box;
         d_orig_box = other.d_orig_box;
         d_boxload = other.d_boxload;
         return *this;
      }

      //! @brief Return the owner rank.
      int
      getOwnerRank() const
      {
         return d_box.getOwnerRank();
      }

      //! @brief Return the LocalId.
      hier::LocalId
      getLocalId() const
      {
         return d_box.getLocalId();
      }

      //! @brief Return the Box.
      hier::Box&
      getBox()
      {
         return d_box;
      }

      //! @brief Return the Box.
      const hier::Box&
      getBox() const
      {
         return d_box;
      }

      /*!
       * @brief Put self into a MessageStream.
       *
       * This is the opposite of getFromMessageStream().
       */
      void
      putToMessageStream(
         tbox::MessageStream &msg) const;

      /*!
       * @brief Set attributes according to data in a MessageStream.
       *
       * This is the opposite of putToMessageStream().
       */
      void
      getFromMessageStream(
         tbox::MessageStream &msg);

      //! @brief The Box.
      hier::Box d_box;

      //! @brief Originating Box.
      hier::Box d_orig_box;

      //! @brief Work load in this box.
      LoadType d_boxload;
   };


   /*!
    * @brief Insert BoxInTransit into an output stream.
    */
   friend std::ostream&
   operator << (
      std::ostream& co,
      const BoxInTransit& r);


   /*!
    * @brief Comparison functor for sorting BoxInTransit from more to
    * less loads.
    */
   struct BoxInTransitMoreLoad {
      /*
       * @brief Compares two BoxInTransit for sorting them from more load
       * to less load.
       */
      bool
      operator () (
         const BoxInTransit& a,
         const BoxInTransit& b) const
      {
         if (a.getBox().size() != b.getBox().size()) {
            return a.d_boxload > b.d_boxload;
         }
         return a.d_box.getId() < b.d_box.getId();
      }
   };



   /*
    * Static integer constants.  Tags are for isolating messages
    * from different phases of the algorithm.
    */
   static const int TreeLoadBalancer_LOADTAG0 = 1;
   static const int TreeLoadBalancer_LOADTAG1 = 2;
   static const int TreeLoadBalancer_EDGETAG0 = 3;
   static const int TreeLoadBalancer_EDGETAG1 = 4;
   static const int TreeLoadBalancer_PREBALANCE0 = 5;
   static const int TreeLoadBalancer_PREBALANCE1 = 6;
   static const int TreeLoadBalancer_FIRSTDATALEN = 1000;

   static const int TreeLoadBalancer_MIN_NPROC_FOR_AUTOMATIC_MULTICYCLE = 65;

   // The following are not implemented, but are provided here for
   // dumb compilers.

   TreeLoadBalancer(
      const TreeLoadBalancer&);

   void
   operator = (
      const TreeLoadBalancer&);

   /*!
    * @brief A set of BoxInTransit, sorted from highest load to lowest load.
    */
   typedef std::set<BoxInTransit, BoxInTransitMoreLoad> TransitSet;


   /*!
    * @brief Data to save for each sending/receiving process and the
    * subtree at that process.
    */
   struct SubtreeLoadData {
      // @brief Constructor.
      SubtreeLoadData();

      /*!
       * @brief Number of processes in subtree
       */
      int d_num_procs;

      /*!
       * @brief Current total work amount in the subtree
       */
      LoadType d_total_work;

      /*!
       * @brief Load exported (or to be exported) to nonlocal process.
       *
       * If the object is for the local process, load_exported means
       * the load exported to the process's *parent*.
       */
      LoadType d_load_exported;

      /*!
       * @brief Load imported from nonlocal process.
       *
       * If the object is for the local process, load_imported means
       * the load imported from the process's *parent*.
       */
      LoadType d_load_imported;

      /*!
       * @brief Ideal work amount for the subtree
       */
      LoadType d_ideal_work;

      /*!
       * @brief Work to export.
       *
       * If the object is for the local process, for_export means for
       * exporting to the process's *parent*.
       */
      TransitSet d_for_export;
   };

   /*
    * @brief Check if there is any pending messages for the private
    * communication and throw an error if there is.
    */
   void
   assertNoMessageForPrivateCommunicator() const;

   /*
    * Read parameters from input database.
    */
   void
   getFromInput(
      const boost::shared_ptr<tbox::Database>& db);

   /*!
    * Move Boxes in balance_box_level from ranks outside of
    * rank_group to ranks inside rank_group.  Modify the given connectors
    * to make them correct following this moving of boxes.
    */
   void
   prebalanceBoxLevel(
      hier::BoxLevel& balance_box_level,
      hier::Connector& balance_to_anchor,
      hier::Connector& anchor_to_balance,
      const tbox::RankGroup& rank_group) const;

   /*!
    * @brief Reassign loads from one TransitSet to another.
    *
    * @param[i] ideal_transfer Amount of load to reassign from src to
    * dst.  If negative, reassign the load from dst to src.
    *
    * @param[io] next_available_index Index for guaranteeing new
    * Boxes are uniquely numbered.
    *
    * @return Amount of load transfered.  If positive, work went from
    * src to dst (if negative, from dst to src).
    *
    */
   LoadType
   reassignLoads(
      TransitSet& src,
      TransitSet& dst,
      hier::LocalId& next_available_index,
      LoadType ideal_transfer ) const;

   /*!
    * @brief Shift load from src to dst by swapping BoxInTransit
    * between them.
    *
    * @param[io] src Source of work, for a positive ideal_transfer.
    *
    * @param[io] dst Destination of work, for a positive ideal_transfer.
    *
    * @param[i] ideal_transfer Amount of load to reassign from src to
    * dst.  If negative, reassign the load from dst to src.
    *
    * @return Amount of load transfered.  If positive, load went
    * from src to dst (if negative, from dst to src).
    */
   LoadType
   shiftLoadsBySwapping(
      TransitSet& src,
      TransitSet& dst,
      LoadType ideal_transfer ) const;

   /*!
    * @brief Shift load from src to dst by various box breaking strategies.
    * choosing the break that gives the best overall penalty.
    *
    * @param[io] src Source of work, for a positive ideal_transfer.
    *
    * @param[io] dst Destination of work, for a positive ideal_transfer.
    *
    * @param next_available_index Index for guaranteeing new
    * Boxes are uniquely numbered.
    *
    * @param ideal_transfer Amount of load to reassign from src to
    * dst.  If negative, reassign the load from dst to src.
    *
    * @return Amount of load transfered.  If positive, load went
    * from src to dst (if negative, from dst to src).
    */
   LoadType
   shiftLoadsByBreaking(
      TransitSet& src,
      TransitSet& dst,
      hier::LocalId& next_available_index,
      LoadType ideal_transfer ) const;

   /*!
    * @brief Find a BoxInTransit in each of the source and destination
    * containers that, when swapped, effects a transfer of the given
    * amount of work from the source to the destination.  Swap the boxes.
    */
   bool
   swapLoadPair(
      TransitSet& src,
      TransitSet& dst,
      LoadType& actual_transfer,
      LoadType ideal_transfer ) const;

   /*!
    * @brief Pack load/boxes for sending.
    */
   void
   packSubtreeLoadData(
      std::vector<int>& msg,
      const SubtreeLoadData& load_data) const;

   /*!
    * @brief Unpack load/boxes received.
    */
   void
   unpackSubtreeLoadData(
      SubtreeLoadData& proc_data,
      TransitSet& receiving_bin,
      hier::LocalId& next_available_index,
      const int* received_data,
      int received_data_length ) const;

   /*!
    * @brief Pack load/boxes for sending.
    */
   void
   packSubtreeLoadData(
      tbox::MessageStream &msg,
      const SubtreeLoadData& load_data) const;

   /*!
    * @brief Unpack load/boxes received.
    */
   void
   unpackSubtreeLoadData(
      SubtreeLoadData& proc_data,
      TransitSet& receiving_bin,
      hier::LocalId& next_available_index,
      tbox::MessageStream &msg ) const;

   /*!
    * @brief Construct semilocal relationships in
    * unbalanced--->balanced Connector.
    *
    * Constructing semilocal unbalanced--->balanced relationships
    * require communication to determine where exported work ended up.
    * This methods does the necessary communication and constructs
    * these relationship in the given Connector.
    *
    * @param [o] unbalanced_to_balanced Connector to store
    * relationships in.
    *
    * @param [i] kept_imports Work that was imported and locally kept.
    */
   void
   constructSemilocalUnbalancedToBalanced(
      hier::Connector &unbalanced_to_balanced,
      const TreeLoadBalancer::TransitSet &kept_imports ) const;

   /*!
    * @brief Break off a given load size from a given Box.
    *
    * @param[i] box Box to break.
    *
    * @param[i] ideal_load_to_break Ideal load to break.
    *
    * @param[o] breakoff Boxes broken off (usually just one).
    *
    * @param[o] leftover Remainder of Box after breakoff is gone.
    *
    * @param[o] brk_load The load broken off.
    */
   bool
   breakOffLoad(
      std::vector<hier::Box>& breakoff,
      std::vector<hier::Box>& leftover,
      double& brk_load,
      const hier::Box& box,
      double ideal_load_to_break ) const;

   /*!
    * @brief Computes surface area of a list of boxes.
    */
   double
   computeBoxSurfaceArea(
      const std::vector<hier::Box>& boxes) const;

   /*!
    * @brief Computes the surface area of a box.
    */
   int
   computeBoxSurfaceArea(
      const hier::Box& box) const;

   double
   combinedBreakingPenalty(
      double balance_penalty,
      double surface_penalty,
      double slender_penalty) const
   {
      double combined_penalty =
         d_balance_penalty_wt * balance_penalty * balance_penalty
         + d_surface_penalty_wt * surface_penalty * surface_penalty
         + d_slender_penalty_wt * slender_penalty * slender_penalty;
      return combined_penalty;
   }

   double
   computeBalancePenalty(
      const std::vector<hier::Box>& a,
      const std::vector<hier::Box>& b,
      double imbalance) const
   {
      NULL_USE(a);
      NULL_USE(b);
      return tbox::MathUtilities<double>::Abs(imbalance);
   }

   double
   computeBalancePenalty(
      const TransitSet& a,
      const TransitSet& b,
      double imbalance) const
   {
      NULL_USE(a);
      NULL_USE(b);
      return tbox::MathUtilities<double>::Abs(imbalance);
   }

   double
   computeBalancePenalty(
      const hier::Box& a,
      double imbalance) const
   {
      NULL_USE(a);
      return tbox::MathUtilities<double>::Abs(imbalance);
   }

   double
   computeSurfacePenalty(
      const std::vector<hier::Box>& a,
      const std::vector<hier::Box>& b) const;

   double
   computeSurfacePenalty(
      const TransitSet& a,
      const TransitSet& b) const;

   double
   computeSurfacePenalty(
      const hier::Box& a) const;

   double
   computeSlenderPenalty(
      const std::vector<hier::Box>& a,
      const std::vector<hier::Box>& b) const;

   double
   computeSlenderPenalty(
      const TransitSet& a,
      const TransitSet& b) const;

   double
   computeSlenderPenalty(
      const hier::Box& a) const;

   bool
   breakOffLoad_planar(
      std::vector<hier::Box>& breakoff,
      std::vector<hier::Box>& leftover,
      double& brk_load,
      const hier::Box& box,
      double ideal_load_to_break,
      const tbox::Array<tbox::Array<bool> >& bad_cuts ) const;

   bool
   breakOffLoad_cubic(
      std::vector<hier::Box>& breakoff,
      std::vector<hier::Box>& leftover,
      double& brk_load,
      const hier::Box& box,
      double ideal_load_to_give,
      const tbox::Array<tbox::Array<bool> >& bad_cuts ) const;

   void
   burstBox(
      std::vector<hier::Box>& boxes,
      const hier::Box& bursty,
      const hier::Box& solid ) const;

   /*
    * Utility functions to determine parameter values for level.
    */
   int
   getWorkloadDataId(
      int level_number) const
   {
      TBOX_ASSERT(level_number >= 0);
      return (level_number < d_workload_data_id.getSize() ?
         d_workload_data_id[level_number] :
         d_master_workload_data_id);
   }

   /*!
    * @brief Compute the load for a Box.
    */
   double
   computeLoad(
      const hier::Box& box) const
   {
      /*
       * Currently only for uniform loads, where the load is equal
       * to the number of cells.  For non-uniform loads, this method
       * needs the patch data index for the load.  It would summ up
       * the individual cell loads in the cell.
       */
      return double(box.size());
   }

   /*!
    * @brief Compute the load for the Box, restricted to where it
    * intersects a given box.
    */
   double
   computeLoad(
      const hier::Box& box,
      const hier::Box& restriction) const
   {
      /*
       * Currently only for uniform loads, where the load is equal
       * to the number of cells.  For non-uniform loads, this method
       * needs the patch data index for the load.  It would summ up
       * the individual cell loads in the overlap region.
       */
      return double((box * restriction).size());
   }

   /*
    * Count the local workload.
    */
   double
   computeLocalLoads(
      const hier::BoxLevel& box_level) const;

   /*!
    * @brief Given an "unbalanced" BoxLevel, compute the BoxLevel that
    * is load-balanced within the given rank_group and compute the
    * mapping between the unbalanced and balanced BoxLevels.
    */
   void
   loadBalanceWithinRankGroup(
      hier::BoxLevel& balance_box_level,
      hier::Connector &balance_to_anchor,
      hier::Connector &anchor_to_balance,
      const tbox::RankGroup& rank_group,
      const double group_sum_load ) const;

   /*!
    * @brief Constrain maximum box sizes in the given BoxLevel and
    * update given Connectors to the changed BoxLevel.
    */
   void
   constrainMaxBoxSizes(
      hier::BoxLevel& box_level,
      hier::Connector &anchor_to_level,
      hier::Connector &level_to_anchor ) const;

   /*!
    * @brief Create the cycle-based RankGroups the local process
    * belongs in.
    *
    * The RankGroup size increases exponentially with the cycle
    * number such that for the last cycle the rank group includes
    * all processes in d_mpi.
    *
    * @param [o] rank_group
    * @param [o] num_groups
    * @param [o] group_num
    * @param [i] cycle_number
    * @param [i] number_of_cycles
    */
   void
   createBalanceRankGroupBasedOnCycles(
      tbox::RankGroup &rank_group,
      int &num_groups,
      int &group_num,
      const int cycle_number,
      const int number_of_cycles) const;

   /*!
    * @brief Set up the asynchronous communication objects for the
    * given RankGroup.
    *
    * Based on a conceptual process tree with num_children children,
    * set the AsyncCommPeer objects for communication with children
    * and parent.
    *
    * @param [o] child_stage
    * @param [o] child_comms
    * @param [o] parent_stage
    * @param [o] parent_comm
    * @param [i] rank_group
    * @param [i] bdfs
    */
   void
   setupAsyncCommObjects(
      tbox::AsyncCommStage& child_stage,
      tbox::AsyncCommPeer<char> *& child_comms,
      tbox::AsyncCommStage& parent_stage,
      tbox::AsyncCommPeer<char> *& parent_comm,
      const tbox::RankGroup &rank_group,
      const tbox::BalancedDepthFirstTree &bdfs ) const;

   /*
    * @brief Undo the set-up done by setupAsyncCommObjects.
    */
   void
   destroyAsyncCommObjects(
      tbox::AsyncCommPeer<char> *& child_comms,
      tbox::AsyncCommPeer<char> *& parent_comm) const;

   /*!
    * @brief Sum up the work in a sequence of boxes.
    *
    * @return Sum of work in [first,last)
    */
   LoadType
   sumWorkInBoxes(
      const TransitSet::const_iterator &first,
      const TransitSet::const_iterator &last ) const
   {
      LoadType sum = 0;
      for ( TransitSet::const_iterator itr=first; itr!=last; ++itr ) {
         sum += itr->d_boxload;
      }
      return sum;
   }

   /*!
    * @brief Set up timers for the object.
    */
   void
   setTimers();

   /*
    * Object dimension.
    */
   const tbox::Dimension d_dim;

   /*
    * String identifier for load balancer object.
    */
   std::string d_object_name;

   //! @brief Duplicated communicator object.  See setSAMRAI_MPI().
   mutable tbox::SAMRAI_MPI d_mpi;

   //! @brief Whether d_mpi is an internal duplicate.  See setSAMRAI_MPI().
   bool d_mpi_is_dupe;

   int d_n_root_cycles;

   //! @brief Degree of the tree.  Two means binary tree.
   const int d_degree;

   /*
    * Values for workload estimate data, workload factor, and bin pack method
    * used on individual levels when specified as such.
    */
   tbox::Array<int> d_workload_data_id;

   int d_master_workload_data_id;

   /*!
    * @brief Additional minimum box size restriction.
    *
    * See input parameter "min_load_fraction_per_box".
    */
   double d_min_load_fraction_per_box;

   /*!
    * @brief Weighting factor for penalizing imbalance.
    *
    * @see combinedBreakingPenalty().
    */
   double d_balance_penalty_wt;

   /*!
    * @brief Weighting factor for penalizing new suraces.
    *
    * @see combinedBreakingPenalty().
    */
   double d_surface_penalty_wt;

   /*!
    * @brief Weighting factor for penalizing slenderness.
    *
    * @see combinedBreakingPenalty().
    */
   double d_slender_penalty_wt;

   /*!
    * @brief How high a slenderness ratio we can tolerate before penalizing.
    */
   double d_slender_penalty_threshold;

   /*!
    * @brief Extra penalty weighting applied before cutting.
    *
    * Set to range [1,ininity).
    * Higher value forces more agressive cutting but can produce more slivers.
    */
   double d_precut_penalty_wt;

   //@{
   //! @name Data shared with private methods during balancing.
   mutable hier::IntVector d_min_size;
   mutable hier::IntVector d_max_size;
   mutable std::vector<hier::BoxContainer> d_block_domain_boxes;
   mutable hier::IntVector d_bad_interval;
   mutable hier::IntVector d_cut_factor;
   mutable double d_global_avg_load;
   //@}

   mutable tbox::Array<int> d_output_procs;
   bool d_using_all_procs;

   /*!
    * @brief Whether to immediately report the results of the load balancing cycles
    * in the log files.
    */
   bool d_report_load_balance;

   /*!
    * @brief See "summarize_map" input parameter.
    */
   char d_summarize_map;

   //@{
   //! @name Used for evaluating peformance.
   bool d_barrier_before;
   bool d_barrier_after;
   //@}

   static const int d_default_data_id;

   /*
    * Performance timers.
    */
   boost::shared_ptr<tbox::Timer> t_load_balance_box_level;
   boost::shared_ptr<tbox::Timer> t_get_map;
   boost::shared_ptr<tbox::Timer> t_use_map;
   boost::shared_ptr<tbox::Timer> t_constrain_size;
   boost::shared_ptr<tbox::Timer> t_map_big_boxes;
   boost::shared_ptr<tbox::Timer> t_load_distribution;
   boost::shared_ptr<tbox::Timer> t_compute_local_load;
   boost::shared_ptr<tbox::Timer> t_compute_global_load;
   boost::shared_ptr<tbox::Timer> t_compute_tree_load;
   std::vector<boost::shared_ptr<tbox::Timer> > t_compute_tree_load_for_cycle;
   boost::shared_ptr<tbox::Timer> t_reassign_loads;
   boost::shared_ptr<tbox::Timer> t_shift_loads_by_swapping;
   boost::shared_ptr<tbox::Timer> t_shift_loads_by_breaking;
   boost::shared_ptr<tbox::Timer> t_find_swap_pair;
   boost::shared_ptr<tbox::Timer> t_break_off_load;
   boost::shared_ptr<tbox::Timer> t_find_bad_cuts;
   boost::shared_ptr<tbox::Timer> t_send_load_to_children;
   boost::shared_ptr<tbox::Timer> t_send_load_to_parent;
   boost::shared_ptr<tbox::Timer> t_get_load_from_children;
   boost::shared_ptr<tbox::Timer> t_get_load_from_parent;
   boost::shared_ptr<tbox::Timer> t_construct_semilocal;
   boost::shared_ptr<tbox::Timer> t_construct_semilocal_comm_wait;
   boost::shared_ptr<tbox::Timer> t_send_edge_to_children;
   boost::shared_ptr<tbox::Timer> t_send_edge_to_parent;
   boost::shared_ptr<tbox::Timer> t_get_edge_from_children;
   boost::shared_ptr<tbox::Timer> t_get_edge_from_parent;
   boost::shared_ptr<tbox::Timer> t_report_loads;
   boost::shared_ptr<tbox::Timer> t_local_balancing;
   boost::shared_ptr<tbox::Timer> t_finish_sends;
   boost::shared_ptr<tbox::Timer> t_pack_load;
   boost::shared_ptr<tbox::Timer> t_unpack_load;
   boost::shared_ptr<tbox::Timer> t_pack_edge;
   boost::shared_ptr<tbox::Timer> t_unpack_edge;
   boost::shared_ptr<tbox::Timer> t_children_load_comm;
   boost::shared_ptr<tbox::Timer> t_parent_load_comm;
   boost::shared_ptr<tbox::Timer> t_children_edge_comm;
   boost::shared_ptr<tbox::Timer> t_parent_edge_comm;
   boost::shared_ptr<tbox::Timer> t_barrier_before;
   boost::shared_ptr<tbox::Timer> t_barrier_after;
   boost::shared_ptr<tbox::Timer> t_child_send_wait;
   boost::shared_ptr<tbox::Timer> t_child_recv_wait;
   boost::shared_ptr<tbox::Timer> t_parent_send_wait;
   boost::shared_ptr<tbox::Timer> t_parent_recv_wait;

   /*
    * Statistics on number of cells and patches generated.
    */
   mutable std::vector<double> d_load_stat;
   mutable std::vector<int> d_box_count_stat;

   // Extra checks independent of optimization/debug.
   char d_print_steps;
   char d_print_break_steps;
   char d_print_swap_steps;
   char d_print_edge_steps;
   char d_check_connectivity;
   char d_check_map;

};

}
}

#endif
