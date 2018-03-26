/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Routines for summing node data at patch boundaries
 *
 ************************************************************************/

#ifndef included_algs_PatchBoundaryNodeSum_C
#define included_algs_PatchBoundaryNodeSum_C

#include "SAMRAI/algs/PatchBoundaryNodeSum.h"

#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeDataFactory.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/pdat/OuternodeData.h"
#include "SAMRAI/pdat/OuternodeDoubleConstantCoarsen.h"
#include "SAMRAI/algs/OuternodeSumTransactionFactory.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <boost/make_shared.hpp>

/*
 *************************************************************************
 *
 * External declarations for FORTRAN 77 routines used to sum node and
 * outernode data.
 *
 *************************************************************************
 */

extern "C" {

#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

// in algs_nodeouternodeops2d.f:
void F77_FUNC(nodeouternodesum2d, NODEOUTERNODESUM2D) (
   const int&, const int&,     // fine patch lo
   const int&, const int&,     // fine patch hi
   const int&, const int&,     // coarse patch lo
   const int&, const int&,     // coarse patch hi
   const int *,                // ratio vector
   const int&,                 // data depth
   const int&, const int&,     // node data gcw
   double *,                   // node array
   const double *, const double *,   // onode arrays
   const double *, const double *);

void F77_FUNC(nodehangnodeinterp2d, NODEHANGNODEINTERP2D) (
   const int&, const int&,    // fine patch lo
   const int&, const int&,    // fine patch hi
   const int&, const int&,    // coarse patch lo
   const int&, const int&,    // coarse patch hi
   const int&, const int&,    // bbox lo
   const int&, const int&,    // bbox hi
   const int&,                // bbox location
   const int *,               // ratio vector
   const int&,                // data depth
   const int&, const int&,    // node data gcw
   double *);                 // node array

// in algs_nodeouternodeops3d.f:
void F77_FUNC(nodeouternodesum3d, NODEOUTERNODESUM3D) (
   const int&, const int&, const int&,    // fine patch lo
   const int&, const int&, const int&,    // fine patch hi
   const int&, const int&, const int&,    // coarse patch lo
   const int&, const int&, const int&,    // coarse patch hi
   const int *,                           // ratio vector
   const int&,                            // data depth
   const int&, const int&, const int&,    // node data gcw
   double *,                              // node array
   const double *, const double *,        // onode arrays
   const double *, const double *,
   const double *, const double *);

void F77_FUNC(nodehangnodeinterp3d, NODEHANGNODEINTERP3D) (
   const int&, const int&, const int&,    // fine patch lo
   const int&, const int&, const int&,    // fine patch hi
   const int&, const int&, const int&,    // coarse patch lo
   const int&, const int&, const int&,    // coarse patch hi
   const int&, const int&, const int&,    // bbox lo
   const int&, const int&, const int&,    // bbox hi
   const int&,                            // bbox location
   const int *,                           // ratio vector
   const int&,                            // data depth
   const int&, const int&, const int&,    // node data gcw
   double *);                             // node array
}

namespace SAMRAI {
namespace algs {

/*
 *************************************************************************
 *
 * Initialize the static data members.
 *
 *************************************************************************
 */

int PatchBoundaryNodeSum::s_instance_counter = 0;

tbox::Array<tbox::Array<int> >
PatchBoundaryNodeSum::s_onode_src_id_array =
   tbox::Array<tbox::Array<int> >(0);
tbox::Array<tbox::Array<int> >
PatchBoundaryNodeSum::s_onode_dst_id_array =
   tbox::Array<tbox::Array<int> >(0);

/*
 *************************************************************************
 *
 * Constructor patch boundary node sum objects initializes data members
 * to default (undefined) states.
 *
 *************************************************************************
 */

PatchBoundaryNodeSum::PatchBoundaryNodeSum(
   const std::string& object_name) :
   d_setup_called(false),
   d_num_reg_sum(0),
   d_coarsest_level(-1),
   d_finest_level(-1),
   d_level_setup_called(false),
   d_hierarchy_setup_called(false),
   d_sum_transaction_factory(boost::make_shared<OuternodeSumTransactionFactory>())
{
   TBOX_ASSERT(!object_name.empty());

   d_object_name = object_name;

   s_instance_counter++;
}

/*
 *************************************************************************
 *
 * Destructor removes temporary outernode patch data ids from
 * variable database, if defined.
 *
 *************************************************************************
 */

PatchBoundaryNodeSum::~PatchBoundaryNodeSum()
{

   s_instance_counter--;
   if (s_instance_counter == 0) {
      const int arr_length_depth = s_onode_src_id_array.size();

      for (int id = 0; id < arr_length_depth; id++) {
         const int arr_length_nvar = s_onode_src_id_array[id].size();

         for (int iv = 0; iv < arr_length_nvar; iv++) {

            if (s_onode_src_id_array[id][iv] >= 0) {
               hier::VariableDatabase::getDatabase()->
               removeInternalSAMRAIVariablePatchDataIndex(
                  s_onode_src_id_array[id][iv]);
            }
            if (s_onode_dst_id_array[id][iv] >= 0) {
               hier::VariableDatabase::getDatabase()->
               removeInternalSAMRAIVariablePatchDataIndex(
                  s_onode_dst_id_array[id][iv]);
            }

            s_onode_src_id_array[id].resizeArray(0);
            s_onode_dst_id_array[id].resizeArray(0);

         }

      }

      s_onode_src_id_array.resizeArray(0);
      s_onode_dst_id_array.resizeArray(0);

   }

}

/*
 *************************************************************************
 *
 * Register node patch data index for summation.
 *
 *************************************************************************
 */

void
PatchBoundaryNodeSum::registerSum(
   int node_data_id)
{

   if (d_setup_called) {
      TBOX_ERROR("PatchBoundaryNodeSum::register error..."
         << "\nobject named " << d_object_name
         << "\nCannot call registerSum with this PatchBoundaryNodeSum"
         << "\nobject since it has already been used to create communication"
         << "\nschedules; i.e., setupSum() has been called."
         << std::endl);
   }

   if (node_data_id < 0) {
      TBOX_ERROR("PatchBoundaryNodeSum register error..."
         << "\nobject named " << d_object_name
         << "\n node_data_id = " << node_data_id
         << " is an invalid patch data identifier." << std::endl);
   }

   hier::VariableDatabase* var_db = hier::VariableDatabase::getDatabase();

   boost::shared_ptr<pdat::NodeDataFactory<double> > node_factory(
      var_db->getPatchDescriptor()->getPatchDataFactory(node_data_id),
      boost::detail::dynamic_cast_tag());

   if (!node_factory) {

      TBOX_ERROR("PatchBoundaryNodeSum register error..."
         << "\nobject named " << d_object_name
         << "\n node_data_id = " << node_data_id
         << " does not correspond to node data of type double." << std::endl);

   } else {

      const tbox::Dimension& dim(node_factory->getDim());

      static std::string tmp_onode_src_variable_name(
         "PatchBoundaryNodeSum__internal-onode-src");
      static std::string tmp_onode_dst_variable_name(
         "PatchBoundaryNodeSum__internal-onode-dst");

      const int reg_sum_id = d_num_reg_sum;

      d_num_reg_sum++;

      d_user_node_data_id.resizeArray(d_num_reg_sum);
      d_user_node_data_id[reg_sum_id] = ID_UNDEFINED;
      d_user_node_depth.resizeArray(d_num_reg_sum);
      d_user_node_depth[reg_sum_id] = ID_UNDEFINED;
      d_tmp_onode_src_variable.resizeArray(d_num_reg_sum);
      d_tmp_onode_dst_variable.resizeArray(d_num_reg_sum);
      d_onode_src_id.resizeArray(d_num_reg_sum);
      d_onode_src_id[reg_sum_id] = ID_UNDEFINED;
      d_onode_dst_id.resizeArray(d_num_reg_sum);
      d_onode_dst_id[reg_sum_id] = ID_UNDEFINED;

      const int data_depth = node_factory->getDepth();
      const int array_by_depth_size = data_depth + 1;

      if (d_num_registered_data_by_depth.size() < array_by_depth_size) {
         const int old_size = d_num_registered_data_by_depth.size();
         const int new_size = array_by_depth_size;

         d_num_registered_data_by_depth.resizeArray(new_size);
         for (int i = old_size; i < new_size; i++) {
            d_num_registered_data_by_depth[i] = 0;
         }
      }

      const int data_depth_id = d_num_registered_data_by_depth[data_depth];
      const int num_data_at_depth = data_depth_id + 1;

      if (s_onode_src_id_array.size() < array_by_depth_size) {
         s_onode_src_id_array.resizeArray(array_by_depth_size);
         s_onode_dst_id_array.resizeArray(array_by_depth_size);
      }

      if (s_onode_src_id_array[data_depth].size() < num_data_at_depth) {
         const int old_size = s_onode_src_id_array[data_depth].size();
         const int new_size = num_data_at_depth;

         s_onode_src_id_array[data_depth].resizeArray(new_size);
         s_onode_dst_id_array[data_depth].resizeArray(new_size);
         for (int i = old_size; i < new_size; i++) {
            s_onode_src_id_array[data_depth][i] = ID_UNDEFINED;
            s_onode_dst_id_array[data_depth][i] = ID_UNDEFINED;
         }
      }

      std::string var_suffix = tbox::Utilities::intToString(data_depth_id, 4)
         + "__depth=" + tbox::Utilities::intToString(data_depth);

      std::string tonode_src_var_name = tmp_onode_src_variable_name
         + var_suffix;
      d_tmp_onode_src_variable[reg_sum_id] = var_db->getVariable(
            tonode_src_var_name);
      if (!d_tmp_onode_src_variable[reg_sum_id]) {
         d_tmp_onode_src_variable[reg_sum_id].reset(
            new pdat::OuternodeVariable<double>(dim,
                                                tonode_src_var_name,
                                                data_depth));
      }

      std::string tonode_dst_var_name = tmp_onode_dst_variable_name
         + var_suffix;
      d_tmp_onode_dst_variable[reg_sum_id] = var_db->getVariable(
            tonode_dst_var_name);
      if (!d_tmp_onode_dst_variable[reg_sum_id]) {
         d_tmp_onode_dst_variable[reg_sum_id].reset(
            new pdat::OuternodeVariable<double>(dim,
                                                tonode_dst_var_name,
                                                data_depth));
      }

      if (s_onode_src_id_array[data_depth][data_depth_id] < 0) {
         s_onode_src_id_array[data_depth][data_depth_id] =
            var_db->registerInternalSAMRAIVariable(
               d_tmp_onode_src_variable[reg_sum_id],
               hier::IntVector::getZero(dim));
      }
      if (s_onode_dst_id_array[data_depth][data_depth_id] < 0) {
         s_onode_dst_id_array[data_depth][data_depth_id] =
            var_db->registerInternalSAMRAIVariable(
               d_tmp_onode_dst_variable[reg_sum_id],
               hier::IntVector::getZero(dim));
      }

      d_user_node_data_id[reg_sum_id] = node_data_id;
      d_user_node_depth[reg_sum_id] = data_depth;

      d_num_registered_data_by_depth[data_depth] = num_data_at_depth;

      d_onode_src_id[reg_sum_id] =
         s_onode_src_id_array[data_depth][data_depth_id];
      d_onode_dst_id[reg_sum_id] =
         s_onode_dst_id_array[data_depth][data_depth_id];

      d_onode_src_data_set.setFlag(d_onode_src_id[reg_sum_id]);
      d_onode_dst_data_set.setFlag(d_onode_dst_id[reg_sum_id]);

   }

}

/*
 *************************************************************************
 *
 * Set up schedule to sum node data around patch boundaries
 * on a single level.
 *
 *************************************************************************
 */

void
PatchBoundaryNodeSum::setupSum(
   const boost::shared_ptr<hier::PatchLevel>& level)
{
   TBOX_ASSERT(level);

   const tbox::Dimension& dim(level->getDim());

   if (d_hierarchy_setup_called) {

      TBOX_ERROR("PatchBoundaryNodeSum::setupSum error...\n"
         << " object named " << d_object_name
         << " already initialized via setupSum() for hierarchy" << std::endl);

   } else {

      d_setup_called = true;
      d_level_setup_called = true;

      d_level = level;

      d_single_level_sum_schedule.resizeArray(1);

      // Communication algorithm for summing outernode values on a level
      xfer::RefineAlgorithm single_level_sum_algorithm(dim);

      for (int i = 0; i < d_num_reg_sum; i++) {
         single_level_sum_algorithm.registerRefine(d_onode_dst_id[i],  // dst data
            d_onode_src_id[i],                                         // src data
            d_onode_dst_id[i],                                         // scratch data
            boost::shared_ptr<hier::RefineOperator>());
      }

      d_single_level_sum_schedule[0] =
         single_level_sum_algorithm.createSchedule(
            d_level,
            (xfer::RefinePatchStrategy *)NULL,
            d_sum_transaction_factory);

   }

}

/*
 *************************************************************************
 *
 * Set up schedule to sum node data around patch boundaries
 * for set of consecutive hierarchy levels.
 *
 *************************************************************************
 */

void
PatchBoundaryNodeSum::setupSum(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int coarsest_level,
   const int finest_level)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((coarsest_level >= 0)
      && (finest_level >= coarsest_level)
      && (finest_level <= hierarchy->getFinestLevelNumber()));

   const tbox::Dimension& dim(hierarchy->getDim());

   if (d_level_setup_called) {

      TBOX_ERROR("PatchBoundaryNodeSum::setupSum error...\n"
         << " object named " << d_object_name
         << " already initialized via setupSum() for single level" << std::endl);

   } else {

      d_setup_called = true;
      d_hierarchy_setup_called = true;

      d_hierarchy = hierarchy;
      d_coarsest_level = coarsest_level;
      d_finest_level = finest_level;

      const int num_levels = d_finest_level + 1;

      d_single_level_sum_schedule.resizeArray(num_levels);
      d_cfbdry_copy_schedule.resizeArray(num_levels);
      d_sync_coarsen_schedule.resizeArray(num_levels);
      d_cfbdry_tmp_level.resizeArray(num_levels);

      d_coarse_fine_boundary.resizeArray(num_levels);

      // Communication algorithm for summing outernode values on each level
      xfer::RefineAlgorithm single_level_sum_algorithm(dim);

      // Communication algorithm for copying node values on each coarser
      // level to outernode values on coarsened version of patches on
      // next finer level
      xfer::RefineAlgorithm cfbdry_copy_algorithm(dim);

      // Communication algorithm for coarsening outernode values on
      // each finer level to node data on next coarser level
      xfer::CoarsenAlgorithm sync_coarsen_algorithm(dim, false);
      boost::shared_ptr<pdat::OuternodeDoubleConstantCoarsen> coarsen_op(
         boost::make_shared<pdat::OuternodeDoubleConstantCoarsen>(dim));

      for (int i = 0; i < d_num_reg_sum; i++) {
         single_level_sum_algorithm.registerRefine(d_onode_dst_id[i],  // dst data
            d_onode_src_id[i],                                         // src data
            d_onode_dst_id[i],                                         // scratch data
            boost::shared_ptr<hier::RefineOperator>());

         cfbdry_copy_algorithm.registerRefine(d_onode_dst_id[i],      // dst data
            d_user_node_data_id[i],                                   // src data
            d_onode_dst_id[i],                                        // scratch data
            boost::shared_ptr<hier::RefineOperator>());

         sync_coarsen_algorithm.registerCoarsen(d_user_node_data_id[i], // dst data
            d_onode_dst_id[i],                                          // src data
            coarsen_op);
      }

      d_single_level_sum_schedule[d_coarsest_level] =
         single_level_sum_algorithm.createSchedule(
            d_hierarchy->getPatchLevel(d_coarsest_level),
            (xfer::RefinePatchStrategy *)NULL,
            d_sum_transaction_factory);

      for (int ln = d_coarsest_level + 1; ln <= d_finest_level; ln++) {

         const int crse_level_num = ln - 1;
         const int fine_level_num = ln;

         boost::shared_ptr<hier::PatchLevel> crse_level(
            d_hierarchy->getPatchLevel(crse_level_num));
         boost::shared_ptr<hier::PatchLevel> fine_level(
            d_hierarchy->getPatchLevel(fine_level_num));

         d_single_level_sum_schedule[fine_level_num] =
            single_level_sum_algorithm.createSchedule(
               fine_level,
               (xfer::RefinePatchStrategy *)NULL,
               d_sum_transaction_factory);

         d_cfbdry_tmp_level[fine_level_num].reset(new hier::PatchLevel(dim));
         d_cfbdry_tmp_level[fine_level_num]->
            setCoarsenedPatchLevel(fine_level,
               fine_level->getRatioToCoarserLevel());
         hier::IntVector crse_tmp_gcw =
            d_hierarchy->getConnector(crse_level_num,
               fine_level_num).getConnectorWidth();
         // Create persistent overlap Connectors for use in schedule construction.
         // FIXME: There are faster ways to get these edges.  BTNG.
         d_cfbdry_tmp_level[fine_level_num]->getBoxLevel()->
            getPersistentOverlapConnectors().createConnector(
               *crse_level->getBoxLevel(),
               crse_tmp_gcw);
         crse_level->getBoxLevel()->getPersistentOverlapConnectors().
            createConnector(
               *d_cfbdry_tmp_level[fine_level_num]->getBoxLevel(),
               crse_tmp_gcw);
         const hier::Connector& crse_to_domain =
            d_cfbdry_tmp_level[fine_level_num]->getBoxLevel()->
               getPersistentOverlapConnectors().createConnector(
                  d_hierarchy->getDomainBoxLevel(),
                  hier::IntVector::getZero(dim));
         const hier::Connector& crse_to_crse =
            d_cfbdry_tmp_level[fine_level_num]->getBoxLevel()->
               getPersistentOverlapConnectors().createConnector(
                  *d_cfbdry_tmp_level[fine_level_num]->getBoxLevel(),
                  hier::IntVector::getOne(dim));

         d_cfbdry_copy_schedule[fine_level_num] =
            cfbdry_copy_algorithm.createSchedule(
               d_cfbdry_tmp_level[fine_level_num],     // dst level
               crse_level);                            // src level

         d_sync_coarsen_schedule[fine_level_num] =
            sync_coarsen_algorithm.createSchedule(
               crse_level,                             // dst level
               fine_level);                            // src level

         d_coarse_fine_boundary[fine_level_num].reset(
            new hier::CoarseFineBoundary(
               *(d_cfbdry_tmp_level[fine_level_num]),
               crse_to_domain,
               crse_to_crse,
               hier::IntVector::getOne(dim)));
#if 1
         tbox::plog << "PatchBoundaryNodeSum::d_coarse_fine_boundary["
                    << fine_level_num << "]:" << std::endl;
         d_coarse_fine_boundary[fine_level_num]->printClassData(tbox::plog);
#endif

      }

   }

}

/*
 *************************************************************************
 *
 * Perform patch boundary node sum across single level or multiple
 * hierarchy levels depending on how object was initialized.  In the
 * single level case, values at shared nodes are summed.  In the
 * multiple-level case, the algorithm involves the following steps:
 *
 *    1) Sum values at shared nodes on each level.
 *    2) Set node values at coarse-fine boundary on each finer level
 *       to sum of fine level values and coarse level values at all
 *       nodes that are shared between the coarse and fine level.
 *
 *       2a) Copy coarser level node values to finer level (outer)node
 *           values at nodes on boundary of patches on a temporary
 *           level that represents the finer level coarsened to the
 *           index space of the coarser level.
 *       2b) Sum (outer)node values at patch boundaries on finer level
 *           and (outer)node values at patch boundaries on coarsened
 *           finer level and set values on finer level to sum.  Note
 *           that the hanging nodes on the finer level may be treated
 *           at this point if specified to do so by the user.
 *
 *    3) Inject (outer)node values around each finer level patch
 *       boundary to corresponding node values on each coarser level.
 *
 *************************************************************************
 */

void
PatchBoundaryNodeSum::computeSum(
   const bool fill_hanging_nodes) const
{

   if (d_level_setup_called) {

      d_level->allocatePatchData(d_onode_src_data_set);
      d_level->allocatePatchData(d_onode_dst_data_set);

      doLevelSum(d_level);

      d_level->deallocatePatchData(d_onode_src_data_set);
      d_level->deallocatePatchData(d_onode_dst_data_set);

   } else {  // assume d_hierarchy_setup_called

      int ln;

      for (ln = d_coarsest_level; ln <= d_finest_level; ln++) {

         boost::shared_ptr<hier::PatchLevel> level(
            d_hierarchy->getPatchLevel(ln));

         level->allocatePatchData(d_onode_src_data_set);
         level->allocatePatchData(d_onode_dst_data_set);

         doLevelSum(level);

      }

      for (ln = d_coarsest_level + 1; ln <= d_finest_level; ln++) {

         boost::shared_ptr<hier::PatchLevel> level(
            d_hierarchy->getPatchLevel(ln));

         d_cfbdry_tmp_level[ln]->allocatePatchData(d_onode_dst_data_set);

         d_cfbdry_copy_schedule[ln]->fillData(0.0, false);

         doLocalCoarseFineBoundarySum(level,
            d_cfbdry_tmp_level[ln],
            d_user_node_data_id,
            d_onode_dst_id,
            fill_hanging_nodes);

         d_cfbdry_tmp_level[ln]->deallocatePatchData(d_onode_dst_data_set);

      }

      for (ln = d_finest_level; ln > d_coarsest_level; ln--) {

         boost::shared_ptr<hier::PatchLevel> level(
            d_hierarchy->getPatchLevel(ln));

         copyNodeToOuternodeOnLevel(level,
            d_user_node_data_id,
            d_onode_dst_id);

         d_sync_coarsen_schedule[ln]->coarsenData();

         level->deallocatePatchData(d_onode_src_data_set);
         level->deallocatePatchData(d_onode_dst_data_set);

      }

      d_hierarchy->getPatchLevel(d_coarsest_level)->
         deallocatePatchData(d_onode_src_data_set);
      d_hierarchy->getPatchLevel(d_coarsest_level)->
         deallocatePatchData(d_onode_dst_data_set);

   }  // if d_hierarchy_setup_called

}

/*
 *************************************************************************
 *
 * Private member function that performs node sum across single level.
 *
 * 1. Copy node data to local outernode data.
 * 2. Transfer and sum outernode data.
 * 3. Copy outernode data back to node data.
 *
 *************************************************************************
 */

void
PatchBoundaryNodeSum::doLevelSum(
   const boost::shared_ptr<hier::PatchLevel>& level) const
{
   TBOX_ASSERT(level);

   copyNodeToOuternodeOnLevel(level,
      d_user_node_data_id,
      d_onode_src_id);

   int schedule_level_number = 0;
   if (!d_level_setup_called) {
      schedule_level_number =
         tbox::MathUtilities<int>::Max(0, level->getLevelNumber());
   }
   d_single_level_sum_schedule[schedule_level_number]->fillData(0.0, false);

   copyOuternodeToNodeOnLevel(level,
      d_onode_dst_id,
      d_user_node_data_id);
}

/*
 *************************************************************************
 *
 * Private member function to set node node data on a fine level at a
 * coarse-fine boundary to the sum of the node values and the associated
 * outernode values on a coarsened version of the fine level.
 *
 * TODO: This function needs more detailed explanations.  BTNG
 *
 *************************************************************************
 */

void
PatchBoundaryNodeSum::doLocalCoarseFineBoundarySum(
   const boost::shared_ptr<hier::PatchLevel>& fine_level,
   const boost::shared_ptr<hier::PatchLevel>& coarsened_fine_level,
   const tbox::Array<int>& node_data_id,
   const tbox::Array<int>& onode_data_id,
   bool fill_hanging_nodes) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(fine_level);
   TBOX_ASSERT(coarsened_fine_level);
   TBOX_ASSERT(node_data_id.size() == onode_data_id.size());
   for (int i = 0; i < node_data_id.size(); i++) {
      TBOX_ASSERT(fine_level->checkAllocated(node_data_id[i]));
      TBOX_ASSERT(coarsened_fine_level->checkAllocated(onode_data_id[i]));
   }
#endif
   TBOX_DIM_ASSERT_CHECK_ARGS2(*fine_level, *coarsened_fine_level);

   const tbox::Dimension& dim(fine_level->getDim());

   const hier::IntVector& ratio(fine_level->getRatioToCoarserLevel());
   const int level_number(fine_level->getLevelNumber());

   for (hier::PatchLevel::iterator ip(fine_level->begin());
        ip != fine_level->end(); ++ip) {

      const tbox::Array<hier::BoundaryBox>& pboundaries =
         d_coarse_fine_boundary[level_number]->getBoundaries(ip->getGlobalId(), 1);
      const int num_bdry_boxes = pboundaries.getSize();

      if (num_bdry_boxes > 0) {

         const boost::shared_ptr<hier::Patch>& fpatch = *ip;
         boost::shared_ptr<hier::Patch> cfpatch(
            coarsened_fine_level->getPatch(fpatch->getGlobalId()));

         const hier::Index& filo = fpatch->getBox().lower();
         const hier::Index& fihi = fpatch->getBox().upper();
         const hier::Index& cilo = cfpatch->getBox().lower();
         const hier::Index& cihi = cfpatch->getBox().upper();

         for (int i = 0; i < node_data_id.size(); i++) {

            boost::shared_ptr<pdat::NodeData<double> > node_data(
               fpatch->getPatchData(node_data_id[i]),
               boost::detail::dynamic_cast_tag());
            boost::shared_ptr<pdat::OuternodeData<double> > onode_data(
               cfpatch->getPatchData(onode_data_id[i]),
               boost::detail::dynamic_cast_tag());

            const hier::IntVector& node_gcw(node_data->getGhostCellWidth());

            pdat::OuternodeData<double> tmp_onode_data(
               cfpatch->getBox(), onode_data->getDepth());
            tmp_onode_data.fillAll(0.0);

            // Copy "coarse" node values on coarse-fine boundary to
            // temporary outernode data arrays.
            for (int ibb0 = 0; ibb0 < num_bdry_boxes; ibb0++) {

               const hier::BoundaryBox& bbox = pboundaries[ibb0];
               const int bbox_loc = bbox.getLocationIndex();

               hier::Box node_bbox =
                  pdat::NodeGeometry::toNodeBox(bbox.getBox());

               switch (bbox_loc) {

                  case 0: {

                     pdat::ArrayData<double>& tmp_onode_data_side_00 =
                        tmp_onode_data.getArrayData(0, 0);
                     if (tmp_onode_data_side_00.isInitialized()) {
                        tmp_onode_data_side_00.
                        copy(onode_data->getArrayData(0, 0), node_bbox);
                     }
                     pdat::ArrayData<double>& tmp_onode_data_side_10 =
                        tmp_onode_data.getArrayData(1, 0);
                     if (tmp_onode_data_side_10.isInitialized()) {
                        tmp_onode_data_side_10.
                        copy(onode_data->getArrayData(1, 0), node_bbox);
                     }
                     pdat::ArrayData<double>& tmp_onode_data_side_11 =
                        tmp_onode_data.getArrayData(1, 1);
                     if (tmp_onode_data_side_11.isInitialized()) {
                        tmp_onode_data_side_11.
                        copy(onode_data->getArrayData(1, 1), node_bbox);
                     }

                     if ((dim == tbox::Dimension(3))) {
                        pdat::ArrayData<double>& tmp_onode_data_side_20 =
                           tmp_onode_data.getArrayData(2, 0);
                        if (tmp_onode_data_side_20.isInitialized()) {
                           tmp_onode_data_side_20.
                           copy(onode_data->getArrayData(2, 0), node_bbox);
                        }
                        pdat::ArrayData<double>& tmp_onode_data_side_21 =
                           tmp_onode_data.getArrayData(2, 1);
                        if (tmp_onode_data_side_21.isInitialized()) {
                           tmp_onode_data_side_21.
                           copy(onode_data->getArrayData(2, 1), node_bbox);
                        }
                     }
                     break;
                  }  // case 0

                  case 1: {
                     pdat::ArrayData<double>& tmp_onode_data_side_01 =
                        tmp_onode_data.getArrayData(0, 1);
                     if (tmp_onode_data_side_01.isInitialized()) {
                        tmp_onode_data_side_01.
                        copy(onode_data->getArrayData(0, 1), node_bbox);
                     }
                     pdat::ArrayData<double>& tmp_onode_data_side_10 =
                        tmp_onode_data.getArrayData(1, 0);
                     if (tmp_onode_data_side_10.isInitialized()) {
                        tmp_onode_data_side_10.
                        copy(onode_data->getArrayData(1, 0), node_bbox);
                     }
                     pdat::ArrayData<double>& tmp_onode_data_side_11 =
                        tmp_onode_data.getArrayData(1, 1);
                     if (tmp_onode_data_side_11.isInitialized()) {
                        tmp_onode_data_side_11.
                        copy(onode_data->getArrayData(1, 1), node_bbox);
                     }

                     if ((dim == tbox::Dimension(3))) {
                        pdat::ArrayData<double>& tmp_onode_data_side_20 =
                           tmp_onode_data.getArrayData(2, 0);
                        if (tmp_onode_data_side_20.isInitialized()) {
                           tmp_onode_data_side_20.
                           copy(onode_data->getArrayData(2, 0), node_bbox);
                        }
                        pdat::ArrayData<double>& tmp_onode_data_side_21 =
                           tmp_onode_data.getArrayData(2, 1);
                        if (tmp_onode_data_side_21.isInitialized()) {
                           tmp_onode_data_side_21.
                           copy(onode_data->getArrayData(2, 1), node_bbox);
                        }
                     }
                     break;
                  } // case 1

                  case 2: {
                     pdat::ArrayData<double>& tmp_onode_data_side_10 =
                        tmp_onode_data.getArrayData(1, 0);
                     if (tmp_onode_data_side_10.isInitialized()) {
                        tmp_onode_data_side_10.
                        copy(onode_data->getArrayData(1, 0), node_bbox);
                     }

                     if ((dim == tbox::Dimension(3))) {
                        pdat::ArrayData<double>& tmp_onode_data_side_20 =
                           tmp_onode_data.getArrayData(2, 0);
                        if (tmp_onode_data_side_20.isInitialized()) {
                           tmp_onode_data_side_20.
                           copy(onode_data->getArrayData(2, 0), node_bbox);
                        }
                        pdat::ArrayData<double>& tmp_onode_data_side_21 =
                           tmp_onode_data.getArrayData(2, 1);
                        if (tmp_onode_data_side_21.isInitialized()) {
                           tmp_onode_data_side_21.
                           copy(onode_data->getArrayData(2, 1), node_bbox);
                        }
                     }
                     break;
                  }  // case 2

                  case 3: {
                     pdat::ArrayData<double>& tmp_onode_data_side_11 =
                        tmp_onode_data.getArrayData(1, 1);
                     if (tmp_onode_data_side_11.isInitialized()) {
                        tmp_onode_data_side_11.
                        copy(onode_data->getArrayData(1, 1), node_bbox);
                     }

                     if ((dim == tbox::Dimension(3))) {
                        pdat::ArrayData<double>& tmp_onode_data_side_20 =
                           tmp_onode_data.getArrayData(2, 0);
                        if (tmp_onode_data_side_20.isInitialized()) {
                           tmp_onode_data_side_20.
                           copy(onode_data->getArrayData(2, 0), node_bbox);
                        }
                        pdat::ArrayData<double>& tmp_onode_data_side_21 =
                           tmp_onode_data.getArrayData(2, 1);
                        if (tmp_onode_data_side_21.isInitialized()) {
                           tmp_onode_data_side_21.
                           copy(onode_data->getArrayData(2, 1), node_bbox);
                        }
                     }
                     break;
                  }  // case 3

                  case 4: {
                     if ((dim == tbox::Dimension(3))) {
                        pdat::ArrayData<double>& tmp_onode_data_side_20 =
                           tmp_onode_data.getArrayData(2, 0);
                        if (tmp_onode_data_side_20.isInitialized()) {
                           tmp_onode_data_side_20.
                           copy(onode_data->getArrayData(2, 0), node_bbox);
                        }
                     }
                     break;
                  }  // case 4

                  case 5: {
                     if ((dim == tbox::Dimension(3))) {
                        pdat::ArrayData<double>& tmp_onode_data_side_21 =
                           tmp_onode_data.getArrayData(2, 1);
                        if (tmp_onode_data_side_21.isInitialized()) {
                           tmp_onode_data_side_21.
                           copy(onode_data->getArrayData(2, 1), node_bbox);
                        }
                     }
                     break;
                  }  // case 5

                  default: {
                     TBOX_ERROR("PatchBoundaryNodeSum::computeSum error...\n"
                        << " object named " << d_object_name
                        << "\n unrecognized coarse-fine boundary box location "
                        << bbox_loc << std::endl);
                  }

               }  // switch(box_loc)

            } // for (int ibb0 ...  iterate over coarse-fine boundary box regions

            // Sum "coarse" node values on coarse-fine boundary.

            if ((dim == tbox::Dimension(2))) {

               if (tmp_onode_data.getArrayData(0, 0).isInitialized() &&
                   tmp_onode_data.getArrayData(0, 1).isInitialized() &&
                   tmp_onode_data.getArrayData(1, 0).isInitialized() &&
                   tmp_onode_data.getArrayData(1, 1).isInitialized()) {

                  F77_FUNC(nodeouternodesum2d, NODEOUTERNODESUM2D) (
                     filo(0), filo(1),
                     fihi(0), fihi(1),
                     cilo(0), cilo(1),
                     cihi(0), cihi(1),
                     &ratio[0],
                     node_data->getDepth(),
                     node_gcw(0), node_gcw(1),
                     node_data->getPointer(),     // node data dst
                     tmp_onode_data.getPointer(0, 0), // x lower src
                     tmp_onode_data.getPointer(0, 1), // x upper src
                     tmp_onode_data.getPointer(1, 0), // y lower src
                     tmp_onode_data.getPointer(1, 1)); // y upper src
               }

            } // (dim == tbox::Dimension(2))

            if ((dim == tbox::Dimension(3))) {

               if (tmp_onode_data.getArrayData(0, 0).isInitialized() &&
                   tmp_onode_data.getArrayData(0, 1).isInitialized() &&
                   tmp_onode_data.getArrayData(1, 0).isInitialized() &&
                   tmp_onode_data.getArrayData(1, 1).isInitialized() &&
                   tmp_onode_data.getArrayData(2, 0).isInitialized() &&
                   tmp_onode_data.getArrayData(2, 1).isInitialized()) {

                  F77_FUNC(nodeouternodesum3d, NODEOUTERNODESUM3D) (
                     filo(0), filo(1), filo(2),
                     fihi(0), fihi(1), fihi(2),
                     cilo(0), cilo(1), cilo(2),
                     cihi(0), cihi(1), cihi(2),
                     &ratio[0],
                     node_data->getDepth(),
                     node_gcw(0), node_gcw(1), node_gcw(2),
                     node_data->getPointer(),     // node data dst
                     tmp_onode_data.getPointer(0, 0), // x lower src
                     tmp_onode_data.getPointer(0, 1), // x upper src
                     tmp_onode_data.getPointer(1, 0), // y lower src
                     tmp_onode_data.getPointer(1, 1), // y upper src
                     tmp_onode_data.getPointer(2, 0), // z lower src
                     tmp_onode_data.getPointer(2, 1)); // z upper src
               }

            } // (dim == tbox::Dimension(3))

            // If desired, fill "hanging" nodes on fine patch by
            // linear interpolation between "coarse" nodes on
            // coarse-fine boundary.
            if (fill_hanging_nodes) {

               for (int ibb1 = 0; ibb1 < num_bdry_boxes; ibb1++) {

                  const hier::BoundaryBox& bbox = pboundaries[ibb1];
                  const hier::Index& bboxilo = bbox.getBox().lower();
                  const hier::Index& bboxihi = bbox.getBox().upper();

                  const int bbox_loc = bbox.getLocationIndex();

                  if ((dim == tbox::Dimension(2))) {
                     F77_FUNC(nodehangnodeinterp2d, NODEHANGNODEINTERP2D) (
                        filo(0), filo(1),
                        fihi(0), fihi(1),
                        cilo(0), cilo(1),
                        cihi(0), cihi(1),
                        bboxilo(0), bboxilo(1),
                        bboxihi(0), bboxihi(1),
                        bbox_loc,
                        &ratio[0],
                        node_data->getDepth(),
                        node_gcw(0), node_gcw(1),
                        node_data->getPointer());
                  }

                  if ((dim == tbox::Dimension(3))) {
                     F77_FUNC(nodehangnodeinterp3d, NODEHANGNODEINTERP3D) (
                        filo(0), filo(1), filo(2),
                        fihi(0), fihi(1), fihi(2),
                        cilo(0), cilo(1), cilo(2),
                        cihi(0), cihi(1), cihi(2),
                        bboxilo(0), bboxilo(1), bboxilo(2),
                        bboxihi(0), bboxihi(1), bboxihi(2),
                        bbox_loc,
                        &ratio[0],
                        node_data->getDepth(),
                        node_gcw(0), node_gcw(1), node_gcw(2),
                        node_data->getPointer());
                  }

               } // iterate over coarse-fine boundary box regions

            } // if fill hanging nodes

         }  //  for (int i ...  iterate over node data

      }  // if (num_bdry_boxes > 0) .... if patch lies on coarse-fine boundary

   } // iterate over fine level patches

}

/*
 *************************************************************************
 *
 * Private member functions to copy between node and outernode data
 * over an entire level.
 *
 *************************************************************************
 */

void
PatchBoundaryNodeSum::copyNodeToOuternodeOnLevel(
   const boost::shared_ptr<hier::PatchLevel>& level,
   const tbox::Array<int>& node_data_id,
   const tbox::Array<int>& onode_data_id) const
{
   TBOX_ASSERT(level);
   TBOX_ASSERT(node_data_id.size() == onode_data_id.size());

   for (hier::PatchLevel::iterator ip(level->begin());
        ip != level->end(); ++ip) {
      const boost::shared_ptr<hier::Patch>& patch = *ip;

      for (int i = 0; i < node_data_id.size(); i++) {
         boost::shared_ptr<pdat::NodeData<double> > node_data(
            patch->getPatchData(node_data_id[i]),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::OuternodeData<double> > onode_data(
            patch->getPatchData(onode_data_id[i]),
            boost::detail::dynamic_cast_tag());

         onode_data->copy(*node_data);
      }
   }

}

void
PatchBoundaryNodeSum::copyOuternodeToNodeOnLevel(
   const boost::shared_ptr<hier::PatchLevel>& level,
   const tbox::Array<int>& onode_data_id,
   const tbox::Array<int>& node_data_id) const
{
   TBOX_ASSERT(level);
   TBOX_ASSERT(node_data_id.size() == onode_data_id.size());

   for (hier::PatchLevel::iterator ip(level->begin());
        ip != level->end(); ++ip) {
      const boost::shared_ptr<hier::Patch>& patch = *ip;

      for (int i = 0; i < node_data_id.size(); i++) {
         boost::shared_ptr<pdat::OuternodeData<double> > onode_data(
            patch->getPatchData(onode_data_id[i]),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::NodeData<double> > node_data(
            patch->getPatchData(node_data_id[i]),
            boost::detail::dynamic_cast_tag());

         onode_data->copy2(*node_data);
      }
   }

}

}
}

#endif
