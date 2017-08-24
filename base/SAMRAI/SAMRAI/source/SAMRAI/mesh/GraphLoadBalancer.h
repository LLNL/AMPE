/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Scalable load balancer using tree algorithm.
 *
 ************************************************************************/

#ifndef included_mesh_GraphLoadBalancer
#define included_mesh_GraphLoadBalancer

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/mesh/BalanceUtilities.h"
#include "SAMRAI/mesh/LoadBalanceStrategy.h"
#include "SAMRAI/hier/MappingConnector.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/CommGraphWriter.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/RankGroup.h"
#include "SAMRAI/tbox/Statistic.h"
#include "SAMRAI/tbox/Statistician.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Utilities.h"

#include "boost/shared_ptr.hpp"
#include <iostream>
#include <vector>
#include <set>

namespace SAMRAI {
namespace mesh {

/*!
 * @brief A graph-based load balance using PT-Scotch to partition patches
 * on a level.
 *
 * This class load balances by treating input boxes as weighted nodes on a
 * graph. The PT-Scotch graph partitioning library is used to distribute these
 * nodes across all the processors.
 *
 * If SAMRAI is not configured with PT-Scotch, this class will still compile,
 * but the load balancing calls will effectively be no-ops.
 *
 * This class is primarily intended to be used as a load balancing option
 * within the TilePartitioner class after clustering has been executed by
 * the TileClustering class, though it is not required to be used in that
 * context.
 *
 * Currently, only uniform load balancing is supported.
 *
 *   - \b target_box_size
 *   The boxes that are recevied by this load balancer as input may be very
 *   large, and may need to be chopped in order to provide a sufficient number
 *   of nodes on the graph that is used to partition the level. The load
 *   balancer will attempt to chop any boxes larger than this size into
 *   pieces as close to this size as possible. When used in conjunction with
 *   the TileClustering box generator, the recommended size is the tile size
 *   times the refinement ratio. If no input is given for this entry, an
 *   internal heuristic will choose a target box size, but the chosen value
 *   may not be well-suited to any particular problem.
 *
 *   - \b coalesce_boxes
 *   After load balancing, coalesce boxes on each local processor wherever
 *   possible if this boolean is true. If coalescing is turned off, the
 *   load balancer may result in levels containing a large number of small
 *   boxes. This option can be invoked to spatially merge those boxes where
 *   possible.
 *
 *   - \b tile_size
 *   Tile size when using tile mode.  Tile mode restricts box cuts
 *   to tile boundaries.  Default is 1, which is equivalent to no restriction.
 *
 * <b> Details: </b> <br>
 * <table>
 *   <tr>
 *     <th>parameter</th>
 *     <th>type</th>
 *     <th>default</th>
 *     <th>range</th>
 *     <th>opt/req</th>
 *     <th>behavior on restart</th>
 *   </tr>
 *   <tr>
 *     <td>target_box_size</td>
 *     <td>int[]</td>
 *     <td>default computed by internal heuristic</td>
 *     <td>all values > 0</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>coalesce_boxes</td>
 *     <td>bool</td>
 *     <td>true</td>
 *     <td>true/false</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 *   <tr>
 *     <td>tile_size</td>
 *     <td>IntVector</td>
 *     <td>1</td>
 *     <td>1-</td>
 *     <td>opt</td>
 *     <td>Not written to restart. Value in input db used.</td>
 *   </tr>
 * </table>
 *
 * @see LoadBalanceStrategy
 */

class GraphLoadBalancer:
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
    * reporting and timer names.  If omitted, "GraphLoadBalancer"
    * is used.
    *
    * @param[in] input_db (optional) database pointer providing
    * parameters from input file.  This pointer may be null indicating
    * no input is used.
    *
    * @pre !name.empty()
    */
   GraphLoadBalancer(
      const tbox::Dimension& dim,
      const std::string& name,
      const boost::shared_ptr<tbox::Database>& input_db =
         boost::shared_ptr<tbox::Database>());

   /*!
    * @brief Virtual destructor releases all internal storage.
    */
   virtual ~GraphLoadBalancer();

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
    *
    * @pre hier::VariableDatabase::getDatabase()->getPatchDescriptor()->getPatchDataFactory(data_id) is actually a  boost::shared_ptr<pdat::CellDataFactory<double> >
    */
   void
   setWorkloadPatchDataIndex(
      int data_id,
      int level_number = -1) {}

   /*!
    * @brief Return true if load balancing procedure for given level
    * depends on patch data on mesh; otherwise return false.
    *
    * @param[in] level_number  Integer patch level number.
    */
   bool
   getLoadBalanceDependsOnPatchData(
      int level_number) const
   {
      NULL_USE(level_number);
      return false;
   }

   /*!
    * @copydoc LoadBalanceStrategy::loadBalanceBoxLevel()
    *
    * Note: This implementation does not yet support non-uniform load
    * balancing.
    *
    * @pre !balance_to_anchor || balance_to_anchor->hasTranspose()
    * @pre !balance_to_anchor || balance_to_anchor->isTransposeOf(balance_to_anchor->getTranspose())
    * @pre (d_dim == balance_box_level.getDim()) &&
    *      (d_dim == min_size.getDim()) && (d_dim == max_size.getDim()) &&
    *      (d_dim == domain_box_level.getDim()) &&
    *      (d_dim == bad_interval.getDim()) && (d_dim == cut_factor.getDim())
    * @pre !hierarchy || (d_dim == hierarchy->getDim())
    * @pre !d_mpi_is_dupe || (d_mpi.getSize() == balance_box_level.getMPI().getSize())
    * @pre !d_mpi_is_dupe || (d_mpi.getSize() == balance_box_level.getMPI().getRank())
    */
   void
   loadBalanceBoxLevel(
      hier::BoxLevel& balance_box_level,
      hier::Connector* balance_to_anchor,
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const hier::IntVector& min_size,
      const hier::IntVector& max_size,
      const hier::BoxLevel& domain_box_level,
      const hier::IntVector& bad_interval,
      const hier::IntVector& cut_factor,
      const tbox::RankGroup& rank_group = tbox::RankGroup()) const;

   /*!
    * @brief Get the name of this object.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

private:
   /*!
    * @brief Struct to store state of a box before and after communication.
    */

   struct BoxInTransit {

      /*!
       * @brief Constructor
       *
       * BoxInTransit will be constructed with empty Boxes.
       *
       * @param[in] dim
       */
      explicit BoxInTransit(
         const tbox::Dimension& dim);

      /*!
       * @brief Construct new object having the history an existing
       * object but is otherwise different.
       *
       * @param[in] orig_box    Original box before communication
       *
       * @param[in] box         Provides spatial coordinates for new box
       *
       * @param[in] rank        Rank of new box
       *
       * @param[in] local_id    LocalId for new box
       */
      BoxInTransit(
         const hier::Box& orig_box,
         const hier::Box& box,
         int rank,
         hier::LocalId local_id);

      /*!
       * @brief Assignment operator
       *
       * @param[in] other
       */
      BoxInTransit&
      operator = (
         const BoxInTransit& other)
      {
         d_box = other.d_box;
         d_orig_box = other.d_orig_box;
         return *this;
      }

      //! @brief Return the owner rank.
      int
      getOwnerRank() const
      {
         return d_box.getOwnerRank();
      }

      //! @brief Return the LocalId
      hier::LocalId
      getLocalId() const
      {
         return d_box.getLocalId();
      }

      //! @brief Return the Box
      hier::Box&
      getBox()
      {
         return d_box;
      }

      //! @brief Return the Box
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
         tbox::MessageStream& msg) const;

      /*!
       * @brief Set attributes according to data in a MessageStream.
       *
       * This is the opposite of putToMessageStream().
       */
      void
      getFromMessageStream(
         tbox::MessageStream& msg);

      hier::Box d_box;
      hier::Box d_orig_box;
   };

   static const int GraphLoadBalancer_LOADTAG0 = 1;
   static const int GraphLoadBalancer_LOADTAG1 = 2;
   static const int GraphLoadBalancer_FIRSTDATALEN = 500;

   /*!
    * @brief Renumber Boxes in a BoxLevel.
    *
    * This method renumbers Boxes to give them globally
    * sequential LocalIds.  Modifies the given Connector to account for the
    * change in the level's state.
    *
    * @param[in,out] balance_box_level
    *
    * @param[in,out] anchor_to_balance
    *
    * @param[in] mca MappingConnectorAlgorithm with timer set
    * in the calling context.
    */
   void
   renumberBoxes(
      hier::BoxLevel& balance_box_level,
      hier::Connector& anchor_to_balance,
      const hier::MappingConnectorAlgorithm& mca) const;

   /*!
    * @brief Chop Boxes in a box_level.
    *
    * If any boxes are greater than the given max_size, attempt to chop them
    * to be below the max_size.  Modifies the given Connector to account for the
    * change in the level's state.
    *
    * @param[in,out] box_level
    *
    * @param[in,out] anchor_to_level
    *
    * @param[in] mca MappingConnectorAlgorithm with timer set
    * in the calling context.
    */
   void
   chopBoxes(
      hier::BoxLevel& box_level,
      hier::Connector* anchor_to_level,
      const hier::IntVector& max_size) const;

   /*!
    * @brief Coalesce Boxes in a box_level.
    *
    * Attempts to coalesce the local boxes on the given level.  Modifies
    * the given Connector to account for the change int eh level's state.
    *
    * @param[in,out] box_level
    *
    * @param[in,out] anchor_to_level
    *
    * @param[in] mca MappingConnectorAlgorithm with timer set
    * in the calling context.
    */
   void
   coalesceBoxLevel(
      hier::BoxLevel& box_level,
      hier::Connector& anchor_to_level,
      const hier::MappingConnectorAlgorithm& mca) const;

   /*!
    * @brief Set up the asynchronous communication objects for communicating
    * Boxes from pre-balanced to post-balanced processors.
    *
    * Based on a conceptual process tree with num_children children,
    * set the AsyncCommPeer objects for communication with children
    * and parent.
    *
    * @param [out] send_stage Stage for sending operations
    * @param [out] send_comms Maps ranks for sending to AsyncCommPeer objects
    * @param [out] send_procs Ranks of processes sent to
    * @param [out] recv_stage Stage for receive operations
    * @param [out] recv_comms Maps ranks for receiving to AsyncCommPeer objects
    * @param [in] old_partition array describing pre-balance partition
    * @param [in] new_partition array describing post-balance partition
    * @param [in] new_boxes global number of boxes
    * @param [in] my_rank
    * @param [in] mpi
    */
   void
   setupAsyncCommObjects(
      tbox::AsyncCommStage& send_stage,
      std::map<int, tbox::AsyncCommPeer<char> *>& send_comms,
      std::set<int>& send_procs,
      tbox::AsyncCommStage& recv_stage,
      std::map<int, tbox::AsyncCommPeer<char> *>& recv_comms,
      const int* old_partition,
      const int* new_partition,
      const int num_boxes,
      const int my_rank,
      const tbox::SAMRAI_MPI& mpi) const;

   /*
    * Read parameters from input database.
    */
   void
   getFromInput(
      const boost::shared_ptr<tbox::Database>& input_db);

   tbox::Dimension d_dim;

   std::string d_object_name;

   hier::IntVector d_target_box_size;

   bool d_coalesce_boxes;

   /*!
    * @brief Tile size, when restricting cuts to tile boundaries,
    * Set to 1 when not restricting.
    */
   hier::IntVector d_tile_size;

   mutable hier::IntVector d_min_size;
   mutable hier::IntVector d_cut_factor;
   mutable hier::IntVector d_bad_interval;
   mutable std::vector<hier::BoxContainer> d_block_domain_boxes;

};

}
}

#endif
