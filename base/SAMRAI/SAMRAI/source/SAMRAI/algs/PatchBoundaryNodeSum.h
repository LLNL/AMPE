/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Routines for summing node data at patch boundaries
 *
 ************************************************************************/

#ifndef included_algs_PatchBoundaryNodeSum
#define included_algs_PatchBoundaryNodeSum

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/OuternodeVariable.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/RefineTransactionFactory.h"
#include "SAMRAI/tbox/Utilities.h"

#include <boost/shared_ptr.hpp>
#include <string>

namespace SAMRAI {
namespace algs {

/*!
 *  @brief Class PatchBoundaryNodeSum provides operations for summing node data
 *  values at nodes that are shared by multiple patches on a single level or
 *  across multiple hierarchy levels.
 *
 *  Usage of a patch boundry node sum involves the following sequence of steps:
 *
 *  -# Construct a patch boundry node sum object.  For example,
 *     \verbatim
 *         PatchBoundaryNodeSum my_node_sum("My Node Sum");
 *     \endverbatim
 *  -# Register node data quantities to sum.  For example,
 *     \verbatim
 *         my_node_sum.registerSum(node_data_id1);
 *         my_node_sum.registerSum(node_data_id2);
 *         etc...
 *     \endverbatim
 *  -# Setup the sum operations for either single level or a range of levels
 *     in a patch hierarchy.  For example,
 *     \verbatim
 *         my_node_sum.setupSum(level);    // single level
 *         -- or --
 *         my_node_sum.setupSum(hierarchy, coarsest_ln, finest_ln);  // multiple levels
 *     \endverbatim
 *  -# Execute the sum operation.  For example,
 *     \verbatim
 *         my_node_sum.computeSum()
 *     \endverbatim
 *
 *  The result of these operations is that each node patch data value
 *  associated with the registered ids at patch boundaries, on either the
 *  single level or range of hierarchy levels, is replaced by the sum of all
 *  data values at the node.
 *
 *  Note that only one of the setupSum() functions may be called once a
 *  PatchBoundaryNodeSum object is created.
 */

class PatchBoundaryNodeSum
{
public:
   /*!
    *  @brief Static function used to predetermine number of patch data
    *         slots ahared among all PatchBoundaryNodeSum
    *         objects (i.e., static members).  To get a correct count,
    *         this routine should only be called once.
    *
    *  @return integer number of internal patch data slots required
    *          to perform sum.
    *  @param max_variables_to_register integer value indicating
    *          maximum number of patch data ids that will be registered
    *          with node sum objects.TO
    */
   static int
   getNumSharedPatchDataSlots(
      int max_variables_to_register)
   {
      // node boundary sum requires two internal outernode variables
      // (source and destination) for each registered variable.
      return 2 * max_variables_to_register;
   }

   /*!
    *  @brief Static function used to predetermine number of patch data
    *         slots unique to each PatchBoundaryNodeSum
    *         object (i.e., non-static members).  To get a correct count,
    *         this routine should be called exactly once for each object
    *         that will be constructed.
    *
    *  @return integer number of internal patch data slots required
    *          to perform sum.
    *  @param max_variables_to_register integer value indicating
    *          maximum number of patch data ids that will be registered
    *          with node sum objects.
    */
   static int
   getNumUniquePatchDataSlots(
      int max_variables_to_register)
   {
      NULL_USE(max_variables_to_register);
      // all patch data slots used by node boundary sum are static
      // and shared among all objects.
      return 0;
   }

   /*!
    *  @brief Constructor initializes object to default (mostly undefined)
    *  state.
    *
    *  @param object_name const std::string reference for name of object used
    *  in error reporting.  When assertion checking is on, the string
    *  cannot be empty.
    */
   explicit PatchBoundaryNodeSum(
      const std::string& object_name);

   /*!
    *  @brief Destructor for the schedule releases all internal storage.
    */
   ~PatchBoundaryNodeSum();

   /*!
    *  @brief Register node data with given patch data identifier for summing.
    *
    *  @param node_data_id  integer patch data index for node data to sum
    *
    *  The node data id must be a valid patch data id (>=0) and must
    *  correspond to node-centered double data.  If not, an error will result.
    */
   void
   registerSum(
      int node_data_id);

   /*!
    *  @brief Set up summation operations for node data across shared nodes
    *         on a single level.
    *
    *  If the other setupSum() function for a range of hierarchy levels has
    *  been called previously for this object, an error will result.
    *
    *  @param level         pointer to level on which to perform node sum
    *
    *  When assertion checking is active, the level pointer cannot be null.
    */
   void
   setupSum(
      const boost::shared_ptr<hier::PatchLevel>& level);

   /*!
    *  @brief Set up for summation operations for node data at shared nodes
    *         across a range of hierarchy levels.
    *
    *  If the other setupSum() function for a single level has been called
    *  previously for this object, an error will result.
    *
    *  @param hierarchy      pointer to hierarchy on which to perform node sum
    *  @param coarsest_level coarsest level number for node sum
    *  @param finest_level   finest level number for node sum
    *
    *  When assertion checking is active, the hierarchy pointer cannot be null,
    *  and the range of levels must be valid.
    */
   void
   setupSum(
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level);

   /*!
    *  @brief Compute sum of node values at each shared node and replace
    *         each such node value with the corresponding sum.
    *
    *  At the end of this method, all values at shared node locations on
    *  patch boundaries (on levels indicated by the call to one of the
    *  setupSum() routines) will have the same value.
    *
    *  When the setupSum() method taking a range of patch levels in a
    *  hierarchy is called, this method will compute the sum of nodal
    *  quantities at all the specified patch boundaries.  For nodes at a
    *  coarse-fine boundary, nodal sums will only be performed where the
    *  coarse and fine nodes overlap.  A node on a fine level that is not
    *  also a node on the next coarser level (a so-called "hanging node")
    *  will not be summed.
    *
    *  The boolean "fill_hanging_nodes" argument specifies whether the
    *  the hanging nodes should be filled using linearly interpolated values
    *  from neighboring non-hanging nodes (i.e. those that overlap nodes on
    *  a coarse level). The correct steps required to deal with hanging
    *  nodes is algorithm dependent so, if left unspecified, values at the
    *  hanging nodes will not be adjusted.  However, because many algorithms
    *  average hanging nodes we provide the capability to do it here.  Note
    *  that the hanging node interpolation provided does not take into
    *  consideration the spatial location of the nodes.  So the interpolation
    *  may not be correct for coordinate systems other than standard Cartesian
    *  grid geometry.
    *
    *  @param fill_hanging_nodes Optional boolean value specifying whether
    *         hanging node values should be set to values interpolated from
    *         neighboring non-hanging node values.  The default is false.
    */
   void
   computeSum(
      const bool fill_hanging_nodes = false) const;

   /*!
    * @brief Returns the object name.
    *
    * @return The object name.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

private:
   /*
    * Private member function to perform node sum across single level --
    * called from computeSum().
    */
   void
   doLevelSum(
      const boost::shared_ptr<hier::PatchLevel>& level) const;

   /*
    * Private member function to set node node data on a fine level at a
    * coarse-fine boundary to the sum of the node values and the associated
    * outernode values on a coarsened version of the fine level.
    */
   void
   doLocalCoarseFineBoundarySum(
      const boost::shared_ptr<hier::PatchLevel>& fine_level,
      const boost::shared_ptr<hier::PatchLevel>& coarsened_fine_level,
      const tbox::Array<int>& node_data_id,
      const tbox::Array<int>& onode_data_id,
      bool fill_hanging_nodes) const;

   /*
    * Private member function to copy node data to outernode data
    * on all patches on a level.
    */
   void
   copyNodeToOuternodeOnLevel(
      const boost::shared_ptr<hier::PatchLevel>& level,
      const tbox::Array<int>& node_data_id,
      const tbox::Array<int>& onode_data_id) const;

   /*
    * Private member function to copy outernode data to node data
    * on all patches on a level.
    */
   void
   copyOuternodeToNodeOnLevel(
      const boost::shared_ptr<hier::PatchLevel>& level,
      const tbox::Array<int>& onode_data_id,
      const tbox::Array<int>& node_data_id) const;

   /*
    * Static members for managing shared temporary data among multiple
    * PatchBoundaryNodeSum objects.
    */
   static int s_instance_counter;
   // These arrays are indexed [data depth][number of variables with depth]
   static tbox::Array<tbox::Array<int> > s_onode_src_id_array;
   static tbox::Array<tbox::Array<int> > s_onode_dst_id_array;

   enum PATCH_BDRY_NODE_SUM_DATA_ID { ID_UNDEFINED = -1 };

   std::string d_object_name;
   bool d_setup_called;

   int d_num_reg_sum;

   // These arrays are indexed [variable registration sequence number]
   tbox::Array<int> d_user_node_data_id;
   tbox::Array<int> d_user_node_depth;

   // These arrays are indexed [data depth]
   tbox::Array<int> d_num_registered_data_by_depth;

   /*
    * Node-centered variables and patch data indices used as internal work
    * quantities.
    */
   // These arrays are indexed [variable registration sequence number]
   tbox::Array<boost::shared_ptr<hier::Variable> > d_tmp_onode_src_variable;
   tbox::Array<boost::shared_ptr<hier::Variable> > d_tmp_onode_dst_variable;

   // These arrays are indexed [variable registration sequence number]
   tbox::Array<int> d_onode_src_id;
   tbox::Array<int> d_onode_dst_id;

   /*
    * Sets of indices for temporary variables to expedite allocation and
    * deallocation.
    */
   hier::ComponentSelector d_onode_src_data_set;
   hier::ComponentSelector d_onode_dst_data_set;

   boost::shared_ptr<hier::PatchLevel> d_level;

   boost::shared_ptr<hier::PatchHierarchy> d_hierarchy;
   int d_coarsest_level;
   int d_finest_level;

   bool d_level_setup_called;
   bool d_hierarchy_setup_called;

   boost::shared_ptr<xfer::RefineTransactionFactory> d_sum_transaction_factory;

   tbox::Array<boost::shared_ptr<xfer::RefineSchedule> >
      d_single_level_sum_schedule;
   tbox::Array<boost::shared_ptr<xfer::RefineSchedule> >
      d_cfbdry_copy_schedule;
   tbox::Array<boost::shared_ptr<xfer::CoarsenSchedule> >
      d_sync_coarsen_schedule;

   tbox::Array<boost::shared_ptr<hier::PatchLevel> > d_cfbdry_tmp_level;

   tbox::Array<boost::shared_ptr<hier::CoarseFineBoundary> >
      d_coarse_fine_boundary;

};

}
}

#endif
