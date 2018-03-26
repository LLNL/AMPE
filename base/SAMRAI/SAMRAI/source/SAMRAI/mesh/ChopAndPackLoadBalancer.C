/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Load balance routines for uniform and non-uniform workloads.
 *
 ************************************************************************/

#ifndef included_mesh_ChopAndPackLoadBalancer_C
#define included_mesh_ChopAndPackLoadBalancer_C

#define ChopAndPackLoadBalancer_MARKLOADFORPOSTPROCESSING

#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/mesh/BalanceUtilities.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellDataFactory.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <boost/make_shared.hpp>
#include <cstdlib>
#include <fstream>
#include <list>

namespace SAMRAI {
namespace mesh {

tbox::StartupShutdownManager::Handler
ChopAndPackLoadBalancer::s_initialize_handler(
   ChopAndPackLoadBalancer::initializeCallback,
   0,
   0,
   ChopAndPackLoadBalancer::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

boost::shared_ptr<tbox::Timer> ChopAndPackLoadBalancer::t_load_balance_boxes;
boost::shared_ptr<tbox::Timer> ChopAndPackLoadBalancer::t_load_balance_boxes_remove_intersection;
boost::shared_ptr<tbox::Timer> ChopAndPackLoadBalancer::t_bin_pack_boxes;
boost::shared_ptr<tbox::Timer> ChopAndPackLoadBalancer::t_bin_pack_boxes_sort;
boost::shared_ptr<tbox::Timer> ChopAndPackLoadBalancer::t_bin_pack_boxes_pack;
boost::shared_ptr<tbox::Timer> ChopAndPackLoadBalancer::t_chop_boxes;

/*
 *************************************************************************
 *
 * Constructors and destructor for ChopAndPackLoadBalancer.
 *
 *************************************************************************
 */

ChopAndPackLoadBalancer::ChopAndPackLoadBalancer(
   const tbox::Dimension& dim,
   const std::string& name,
   const boost::shared_ptr<tbox::Database>& input_db):
   d_dim(dim),
   d_object_name(name),
   d_processor_layout_specified(false),
   d_processor_layout(d_dim),
   d_ignore_level_box_union_is_single_box(false),
   d_master_workload_data_id(-1),
   d_master_max_workload_factor(1.0),
   d_master_workload_tolerance(0.0),
   d_master_bin_pack_method("SPATIAL")
{
   TBOX_ASSERT(!name.empty());

   d_workload_data_id.resizeArray(0);
   d_max_workload_factor.resizeArray(0);
   d_workload_tolerance.resizeArray(0);
   d_bin_pack_method.resizeArray(0);

   getFromInput(input_db);
}

ChopAndPackLoadBalancer::ChopAndPackLoadBalancer(
   const tbox::Dimension& dim,
   const boost::shared_ptr<tbox::Database>& input_db):
   d_dim(dim),
   d_object_name("ChopAndPackLoadBalancer"),
   d_processor_layout_specified(false),
   d_processor_layout(hier::IntVector::getZero(d_dim)),
   d_master_workload_data_id(-1),
   d_master_max_workload_factor(1.0),
   d_master_workload_tolerance(0.0),
   d_master_bin_pack_method("SPATIAL")

{

   d_workload_data_id.resizeArray(0);
   d_max_workload_factor.resizeArray(0);
   d_workload_tolerance.resizeArray(0);
   d_bin_pack_method.resizeArray(0);

   d_ignore_level_box_union_is_single_box = false;

   getFromInput(input_db);
}

ChopAndPackLoadBalancer::~ChopAndPackLoadBalancer()
{
}

/*
 *************************************************************************
 *
 * Accessory functions to get/set load balancing parameters.
 *
 *************************************************************************
 */

bool
ChopAndPackLoadBalancer::getLoadBalanceDependsOnPatchData(
   int level_number) const
{
   return getWorkloadDataId(level_number) < 0 ? false : true;
}

void
ChopAndPackLoadBalancer::setMaxWorkloadFactor(
   double factor,
   int level_number)
{
   TBOX_ASSERT(factor > 0.0);
   if (level_number >= 0) {
      int asize = d_max_workload_factor.getSize();
      if (asize < level_number + 1) {
         d_max_workload_factor.resizeArray(level_number + 1);
         for (int i = asize; i < level_number - 1; i++) {
            d_max_workload_factor[i] =
               d_master_max_workload_factor;
         }
         d_max_workload_factor[level_number] = factor;
      }
   } else {
      d_master_max_workload_factor = factor;
      for (int ln = 0; ln < d_max_workload_factor.getSize(); ln++) {
         d_max_workload_factor[ln] = d_master_max_workload_factor;
      }
   }
}

void
ChopAndPackLoadBalancer::setWorkloadTolerance(
   double tolerance,
   int level_number)
{
   TBOX_ASSERT(tolerance > 0.0);
   if (level_number >= 0) {
      int asize = d_workload_tolerance.getSize();
      if (asize < level_number + 1) {
         d_workload_tolerance.resizeArray(level_number + 1);
         for (int i = asize; i < level_number - 1; i++) {
            d_workload_tolerance[i] =
               d_master_workload_tolerance;
         }
         d_workload_tolerance[level_number] = tolerance;
      }
   } else {
      d_master_workload_tolerance = tolerance;
      for (int ln = 0; ln < d_workload_tolerance.getSize(); ln++) {
         d_workload_tolerance[ln] = d_master_workload_tolerance;
      }
   }
}

void
ChopAndPackLoadBalancer::setWorkloadPatchDataIndex(
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

void
ChopAndPackLoadBalancer::setUniformWorkload(
   int level_number)
{
   if (level_number >= 0) {
      int asize = d_workload_data_id.getSize();
      if (asize < level_number + 1) {
         d_workload_data_id.resizeArray(level_number + 1);
         for (int i = asize; i < level_number - 1; i++) {
            d_workload_data_id[i] =
               d_master_workload_data_id;
         }
         d_workload_data_id[level_number] = -1;
      }
   } else {
      d_master_workload_data_id = -1;
      for (int ln = 0; ln < d_workload_data_id.getSize(); ln++) {
         d_workload_data_id[ln] = d_master_workload_data_id;
      }
   }
}

void
ChopAndPackLoadBalancer::setBinPackMethod(
   const std::string& method,
   int level_number)
{

   if (!(method == "GREEDY" ||
         method == "SPATIAL")) {
      TBOX_ERROR(
         d_object_name << " error: "
                       << "\n   String " << method
                       << " passed to setBinPackMethod()"
                       << " is not a valid method string identifier."
                       << std::endl);

   }

   if (level_number >= 0) {
      int asize = d_bin_pack_method.getSize();
      if (asize < level_number + 1) {
         d_bin_pack_method.resizeArray(level_number + 1);
         for (int i = asize; i < level_number - 1; i++) {
            d_bin_pack_method[i] = d_master_bin_pack_method;
         }
         d_bin_pack_method[level_number] = method;
      }
   } else {
      d_master_bin_pack_method = method;
      for (int ln = 0; ln < d_bin_pack_method.getSize(); ln++) {
         d_bin_pack_method[ln] = d_master_bin_pack_method;
      }
   }
}

/*
 *************************************************************************
 *************************************************************************
 */
void
ChopAndPackLoadBalancer::loadBalanceBoxLevel(
   hier::BoxLevel& balance_mapped_box_level,
   hier::Connector& balance_to_anchor,
   hier::Connector& anchor_to_balance,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const hier::Connector& balance_to_attractor,
   const hier::Connector& attractor_to_balance,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::BoxLevel& domain_mapped_box_level,
   const hier::IntVector& bad_interval,
   const hier::IntVector& cut_factor,
   const tbox::RankGroup& rank_group) const
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS6(d_dim,
      balance_mapped_box_level,
      min_size,
      max_size,
      domain_mapped_box_level,
      bad_interval,
      cut_factor);
   NULL_USE(balance_to_attractor);
   NULL_USE(attractor_to_balance);
   NULL_USE(rank_group);

   hier::IntVector actual_max_size = max_size;
   for (int d = 0; d < d_dim.getValue(); ++d) {
      if (actual_max_size(d) < 0) {
         actual_max_size(d) = tbox::MathUtilities<int>::getMax();
      }
   }

   hier::BoxLevel globalized_input_mapped_box_level(
      balance_mapped_box_level);
   globalized_input_mapped_box_level.setParallelState(
      hier::BoxLevel::GLOBALIZED);

   hier::BoxContainer in_boxes;
   const hier::BoxContainer globalized_input_mapped_boxes(
      globalized_input_mapped_box_level.getGlobalBoxes());
   for (hier::RealBoxConstIterator bi(globalized_input_mapped_boxes.realBegin());
        bi != globalized_input_mapped_boxes.realEnd(); ++bi) {
      in_boxes.pushBack(*bi);
   }

   hier::BoxContainer physical_domain;
   domain_mapped_box_level.getGlobalBoxes(physical_domain);

   hier::BoxContainer out_boxes;
   hier::ProcessorMapping mapping;

   loadBalanceBoxes(
      out_boxes,
      mapping,
      in_boxes,
      hierarchy,
      level_number,
      physical_domain,
      balance_mapped_box_level.getRefinementRatio(),
      min_size,
      actual_max_size,
      cut_factor,
      bad_interval);

   // Build up balance_mapped_box_level from old-style data.
   balance_mapped_box_level.initialize(
      balance_mapped_box_level.getRefinementRatio(),
      balance_mapped_box_level.getGridGeometry(),
      balance_mapped_box_level.getMPI(),
      hier::BoxLevel::GLOBALIZED);
   int i = 0;
   for (hier::BoxContainer::iterator itr(out_boxes);
        itr != out_boxes.end(); ++itr, ++i) {
      hier::Box node(*itr, hier::LocalId(i),
                     mapping.getProcessorAssignment(i));
      balance_mapped_box_level.addBox(node);
   }
   // Reinitialize Connectors due to changed balance_mapped_box_level.
   balance_to_anchor.clearNeighborhoods();
   balance_to_anchor.setBase(balance_mapped_box_level, true);
   anchor_to_balance.clearNeighborhoods();
   anchor_to_balance.setHead(balance_mapped_box_level, true);
   hier::OverlapConnectorAlgorithm oca;
   oca.findOverlaps(balance_to_anchor);
   oca.findOverlaps(anchor_to_balance, balance_mapped_box_level);

   balance_mapped_box_level.setParallelState(hier::BoxLevel::DISTRIBUTED);
}

/*
 *************************************************************************
 *
 * This main load balance routine performs either uniform or
 * non-uniform load balancing operations on the given level depending
 * on user specifications.   In either case, the goal is to produce
 * a set of boxes and a mapping of those boxes to processors so that
 * the workload on each processor is close to the average workload.
 * The average workload is the total computational workload divided
 * by the number of processors.  In the uniform load balance case
 * (default), the workload is the number of cells in each box.  In the
 * non-uniform case, the workload is computed using weight data on the
 * grid hierarchy (i.e., a cell-centered double array on each patch.
 *
 * Typically, any box whose load is larger than the average is chopped.
 * A user can prescribe a parameter (the 'max workload factor') to alter
 * average load used in this computation. Chopping is done using the
 * BalanceUtilities::recursiveBisection()) method which is similar
 * to the Berger-Rigoutsos algorithm.
 *
 * Once the boxes are chopped into a collection os smaller boxes, they
 * are assigned to processors by a bin packing algorithm.
 *
 * The algorithm is summarized as follows:
 *
 * 1) Compute the estimated workload associated with each box.  In the
 *    uniform workload case, the load in each box region is the number
 *    of cells in the region.  Otherwise, the workload is computed using
 *    patch data defined by the d_workload_data_id array set by the user.
 *
 * 2) Compute the maximum workload allowed on any box.  This quantity is
 *    by default the total workload divided by the number of processors.
 *    The user may provide a maximum workload factor, either through the
 *    input file or through a member function, which can alter the
 *    average workload used in this computation.
 *
 * 3) Chop each box whose workload is more than the max allowed into a
 *    smaller set of boxes.
 *
 * 4) Check constraints placed on the boxes by the problem - i.e.
 *    verify boxes are within the maximum and minimum box size
 *    constraints and maintain a correct cut factor.
 *
 * 5) Sort boxes largest to smallest and form an array.  Also form an
 *    array of the workloads associated with each box.
 *
 * 6) Use a bin packing procedure to construct a processor mapping for
 *    the set of boxes.
 *
 *************************************************************************
 */

void
ChopAndPackLoadBalancer::loadBalanceBoxes(
   hier::BoxContainer& out_boxes,
   hier::ProcessorMapping& mapping,
   const hier::BoxContainer& in_boxes,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int level_number,
   const hier::BoxContainer& physical_domain,
   const hier::IntVector& ratio_to_hierarchy_level_zero,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::IntVector& cut_factor,
   const hier::IntVector& bad_interval) const
{
   t_load_balance_boxes->start();

   TBOX_DIM_ASSERT_CHECK_ARGS5(ratio_to_hierarchy_level_zero,
      min_size,
      max_size,
      cut_factor,
      bad_interval);

   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(level_number >= 0);
   TBOX_ASSERT(physical_domain.size() > 0);
   TBOX_ASSERT(min_size > hier::IntVector::getZero(d_dim));
   TBOX_ASSERT(max_size >= min_size);
   TBOX_ASSERT(cut_factor > hier::IntVector::getZero(d_dim));
   TBOX_ASSERT(bad_interval >= hier::IntVector::getZero(d_dim));

   /*
    * This method assumes in_boxes is not empty and will fail
    * if it is.  So shortcut it for empty in_boxes.
    */
   if (in_boxes.isEmpty()) {
      out_boxes = hier::BoxContainer();
      return;
   }

   const tbox::SAMRAI_MPI& mpi(hierarchy->getMPI());

   /*
    * If uniform load balancing is used and the level domain can be
    * expressed as a single box, we can construct an optimal box
    * layout across processors without more involved chopping operations.
    *
    * Otherwise, we chop each box individually to construct the new array
    * of boxes and associated array of workloads based on either uniform
    * or nonuniform workload estimates.
    */

   int wrk_indx = getWorkloadDataId(level_number);

   tbox::Array<double> workloads;

   if ((wrk_indx < 0) || (hierarchy->getNumberOfLevels() == 0)) {

      if (!d_ignore_level_box_union_is_single_box) {
         hier::Box bbox = in_boxes.getBoundingBox();
         hier::BoxContainer difference(bbox);
         t_load_balance_boxes_remove_intersection->start();
         difference.removeIntersections(in_boxes);
         t_load_balance_boxes_remove_intersection->stop();

         if (difference.isEmpty()) {

            t_chop_boxes->start();
            chopUniformSingleBox(out_boxes,
               workloads,
               bbox,
               min_size,
               max_size,
               cut_factor,
               bad_interval,
               physical_domain,
               mpi);
            t_chop_boxes->stop();

         } else {

            t_chop_boxes->start();
            chopBoxesWithUniformWorkload(out_boxes,
               workloads,
               in_boxes,
               hierarchy,
               level_number,
               min_size,
               max_size,
               cut_factor,
               bad_interval,
               physical_domain,
               mpi);
            t_chop_boxes->stop();

         }
      } else {
         t_chop_boxes->start();
         chopBoxesWithUniformWorkload(out_boxes,
            workloads,
            in_boxes,
            hierarchy,
            level_number,
            min_size,
            max_size,
            cut_factor,
            bad_interval,
            physical_domain,
            mpi);
         t_chop_boxes->stop();
      }

   } else {

      t_chop_boxes->start();
      chopBoxesWithNonuniformWorkload(out_boxes,
         workloads,
         in_boxes,
         hierarchy,
         level_number,
         ratio_to_hierarchy_level_zero,
         wrk_indx,
         min_size,
         max_size,
         cut_factor,
         bad_interval,
         physical_domain,
         mpi);
      t_chop_boxes->stop();

   }

#if 0
   /*
    * Assertion check for additional debugging - make sure new boxes
    * satisfy min size, cut_factor, bad_interval, and physical domain
    * constraints.
    */
#ifdef DEBUG_CHECK_ASSERTIONS
   const int nboxes = out_boxes.size();
   for (hier::BoxContainer::iterator ib(out_boxes);
        ib != out_boxes.end(); ++ib) {
      hier::BoxUtilities::checkBoxConstraints(*ib,
         min_size,
         cut_factor,
         bad_interval,
         physical_domain);

   }
#endif
#endif

   /*
    * Generate mapping of boxes to processors using workload estimates.
    * Note boxes in array may be reordered during this operation.
    */

   binPackBoxes(out_boxes,
      mapping,
      workloads,
      getBinPackMethod(level_number));

   t_load_balance_boxes->stop();

#if 0
   /*
    * For debugging, output load balance statistics
    * (assuming uniform load).
    */
   tbox::Array<double> procloads(tbox::SAMRAI_MPI::getNodes());
   for (int i = 0; i < procloads.size(); ++i) {
      procloads[i] = 0;
   }
   int itrCt = 0;
   for (hier::BoxContainer::iterator itr(out_boxes);
        itr != out_boxes.end(); ++itr, ++itrCt) {
      int p = mapping.getProcessorAssignment(itrCt);
      procloads[p] += itr->size();
   }
   tbox::plog << "ChopAndPackLoadBalancer results (after):\n";
   reportLoadBalance(procloads);
#endif

#ifdef ChopAndPackLoadBalancer_MARKLOADFORPOSTPROCESSING
   // Performance: Output loads for global postprocessing.
   const tbox::Array<int>& local_indices = mapping.getLocalIndices();
   double local_load = 0;
   int local_indices_idx = 0;
   int idx = 0;
   for (hier::BoxContainer::iterator itr(out_boxes);
        itr != out_boxes.end() && local_indices_idx < local_indices.size();
        ++itr, ++idx) {
      if (local_indices[local_indices_idx] == idx) {
         local_load += itr->size();
         ++local_indices_idx;
      }
   }
   markLoadForPostprocessing(mpi.getSize(),
      local_load,
      local_indices.size());
#endif
}

/*
 *************************************************************************
 *
 * Private function that chops a single box into a set of boxes
 * that will map to the array of processors as uniformly as possible.
 * The routine assumes a spatially-uniform workload distribution.
 * Note that no error checking is done.
 *
 *************************************************************************
 */

void
ChopAndPackLoadBalancer::chopUniformSingleBox(
   hier::BoxContainer& out_boxes,
   tbox::Array<double>& out_workloads,
   const hier::Box& in_box,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::IntVector& cut_factor,
   const hier::IntVector& bad_interval,
   const hier::BoxContainer& physical_domain,
   const tbox::SAMRAI_MPI& mpi) const
{

   TBOX_DIM_ASSERT_CHECK_ARGS4(min_size, max_size, cut_factor, bad_interval);

   /*
    * Determine processor layout that corresponds to box dimensions.
    */

   hier::IntVector processor_distribution(d_dim);
   if (d_processor_layout_specified) {
      processor_distribution = d_processor_layout;
   } else {
      BalanceUtilities::computeDomainDependentProcessorLayout(
         processor_distribution,
         mpi.getSize(),
         in_box);
   }

   /*
    * The ideal box size will be the dimensions of the input box divided
    * by the number of processors in each direction.  Compute this
    * ideal size and then adjust as necessary to fit within min/max size
    * constraints.
    */

   hier::IntVector ideal_box_size(d_dim);
   for (int i = 0; i < d_dim.getValue(); i++) {
      ideal_box_size(i) = (int)ceil((double)in_box.numberCells(
               i) / (double)processor_distribution(i));

      ideal_box_size(i) = (ideal_box_size(i) > max_size(i) ?
                           max_size(i) : ideal_box_size(i));
      ideal_box_size(i) = (ideal_box_size(i) < min_size(i) ?
                           min_size(i) : ideal_box_size(i));
   }

   /*
    * Chop the single input box into a set of smaller boxes using the
    * ideal_box_size as the maximum size of each of the smaller boxes.
    */

   hier::BoxContainer tmp_box_list(in_box);

   hier::BoxUtilities::chopBoxes(tmp_box_list,
      ideal_box_size,
      min_size,
      cut_factor,
      bad_interval,
      physical_domain);

   /*
    * Set output box array to list of chopped boxes and set workload array.
    */

   out_boxes = tmp_box_list;

   const int nboxes = out_boxes.size();
   out_workloads.resizeArray(nboxes);
   int ibCt = 0;
   for (hier::BoxContainer::iterator ib(out_boxes); ib != out_boxes.end();
        ++ib, ++ibCt) {
      out_workloads[ibCt] = (double)(ib->size());
   }

}

/*
 *************************************************************************
 *
 * Private function that chops boxes on a list into another list of
 * boxes all of which have approximately or less than an average
 * workload estimate.  The routine assumes a spatially-uniform
 * workload distribution.   Note that no error checking is done.
 *
 *************************************************************************
 */

void
ChopAndPackLoadBalancer::chopBoxesWithUniformWorkload(
   hier::BoxContainer& out_boxes,
   tbox::Array<double>& out_workloads,
   const hier::BoxContainer& in_boxes,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int level_number,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::IntVector& cut_factor,
   const hier::IntVector& bad_interval,
   const hier::BoxContainer& physical_domain,
   const tbox::SAMRAI_MPI& mpi) const
{
   NULL_USE(hierarchy);
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS5(d_dim,
      *hierarchy,
      min_size,
      max_size,
      cut_factor,
      bad_interval);

   /*
    * Create copy of input box list to prevent changing it.
    */

   hier::BoxContainer tmp_in_boxes_list(in_boxes);

   /*
    * Chop any boxes in input box list that are larger than max box size
    * if possible.
    */

   hier::BoxUtilities::chopBoxes(tmp_in_boxes_list,
      max_size,
      min_size,
      cut_factor,
      bad_interval,
      physical_domain);

   double total_work = 0.0;
   for (hier::BoxContainer::iterator ib0(tmp_in_boxes_list);
        ib0 != tmp_in_boxes_list.end(); ++ib0) {
      total_work += ib0->size();
   }

   double work_factor = getMaxWorkloadFactor(level_number);
   double average_work = work_factor * total_work / mpi.getSize();

   hier::BoxContainer tmp_box_list;
   std::list<double> tmp_work_list;
   BalanceUtilities::recursiveBisectionUniform(tmp_box_list,
      tmp_work_list,
      tmp_in_boxes_list,
      average_work,
      getWorkloadTolerance(level_number),
      min_size,
      cut_factor,
      bad_interval,
      physical_domain);

   if (tmp_box_list.size() != static_cast<int>(tmp_work_list.size())) {
      TBOX_ERROR(
         d_object_name << ": "
                       << "Number of boxes generated != number of workload values generated."
                       << std::endl);
   }

   /*
    * Set output box array to list of chopped boxes and set workload array.
    */

   out_boxes = tmp_box_list;

   out_workloads.resizeArray(out_boxes.size());
   int i = 0;
   for (std::list<double>::const_iterator il(tmp_work_list.begin());
        il != tmp_work_list.end(); il++) {
      out_workloads[i] = *il;
      i++;
   }

}

/*
 *************************************************************************
 *
 * Private function that chops boxes on a list into another list of
 * boxes all of which have approximately or less than an average
 * workload estimate.  The routine assumes a spatially-nonuniform
 * workload distribution.  Note that no error checking is done.
 *
 *************************************************************************
 */

void
ChopAndPackLoadBalancer::chopBoxesWithNonuniformWorkload(
   hier::BoxContainer& out_boxes,
   tbox::Array<double>& out_workloads,
   const hier::BoxContainer& in_boxes,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int level_number,
   const hier::IntVector& ratio_to_hierarchy_level_zero,
   int wrk_indx,
   const hier::IntVector& min_size,
   const hier::IntVector& max_size,
   const hier::IntVector& cut_factor,
   const hier::IntVector& bad_interval,
   const hier::BoxContainer& physical_domain,
   const tbox::SAMRAI_MPI& mpi) const
{

   TBOX_DIM_ASSERT_CHECK_ARGS5(ratio_to_hierarchy_level_zero,
      min_size,
      max_size,
      cut_factor,
      bad_interval);

   /*
    * Create copy of input box list to prevent changing it.
    */

   hier::BoxContainer tmp_in_boxes_list(in_boxes);

   hier::BoxUtilities::chopBoxes(tmp_in_boxes_list,
      max_size,
      min_size,
      cut_factor,
      bad_interval,
      physical_domain);

   /*
    * Create temporary patch level from in_boxes, distributed using simple
    * uniform workload estimate.  Then, fill the patch data on the level
    * corresponding to the non-uniform workload estimate.  Next, accumulate
    * the total work for the set of boxes.
    */

   hier::BoxContainer tmp_level_boxes(tmp_in_boxes_list);

   const int num_tmp_patches = tmp_level_boxes.size();
   tbox::Array<double> tmp_level_workloads(num_tmp_patches);
   int idx = 0;
   for (hier::BoxContainer::iterator i(tmp_level_boxes);
        i != tmp_level_boxes.end();
        ++i, ++idx) {
      tmp_level_workloads[idx] = i->size();
   }

   hier::ProcessorMapping tmp_level_mapping;

   binPackBoxes(tmp_level_boxes,
      tmp_level_mapping,
      tmp_level_workloads,
      "GREEDY");

   boost::shared_ptr<hier::BoxLevel> tmp_mapped_box_level(
      boost::make_shared<hier::BoxLevel>(
         ratio_to_hierarchy_level_zero,
         hierarchy->getGridGeometry(),
         mpi,
         hier::BoxLevel::GLOBALIZED));
   idx = 0;
   for (hier::BoxContainer::iterator i(tmp_level_boxes);
        i != tmp_level_boxes.end(); ++i, ++idx) {
      hier::Box node(*i, hier::LocalId(idx),
                     tmp_level_mapping.getProcessorAssignment(idx));
      tmp_mapped_box_level->addBox(node);
   }

   boost::shared_ptr<hier::PatchLevel> tmp_level(
      boost::make_shared<hier::PatchLevel>(*tmp_mapped_box_level,
         hierarchy->getGridGeometry(),
         hierarchy->getPatchDescriptor()));

   tmp_level->allocatePatchData(wrk_indx);

   xfer::RefineAlgorithm fill_work_algorithm(d_dim);

   boost::shared_ptr<hier::RefineOperator> work_refine_op(
      boost::make_shared<pdat::CellDoubleConstantRefine>(d_dim));

   fill_work_algorithm.registerRefine(wrk_indx,
      wrk_indx,
      wrk_indx,
      work_refine_op);

   boost::shared_ptr<hier::PatchLevel> current_level;
   if (level_number <= hierarchy->getFinestLevelNumber()) {
      current_level = hierarchy->getPatchLevel(level_number);
   }

   fill_work_algorithm.createSchedule(tmp_level,
      current_level,
      level_number - 1,
      hierarchy)->fillData(0.0);

   double local_work = 0;
   for (hier::PatchLevel::iterator ip(tmp_level->begin());
        ip != tmp_level->end(); ++ip) {
      const boost::shared_ptr<hier::Patch>& patch = *ip;

      double patch_work =
         BalanceUtilities::computeNonUniformWorkload(patch,
            wrk_indx,
            patch->getBox());

      local_work += patch_work;
   }

   double total_work = local_work;
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&total_work, 1, MPI_SUM);
   }

   double work_factor = getMaxWorkloadFactor(level_number);
   double average_work = work_factor * total_work / mpi.getSize();

   hier::BoxContainer tmp_box_list;
   std::list<double> tmp_work_list;
   BalanceUtilities::recursiveBisectionNonuniform(tmp_box_list,
      tmp_work_list,
      tmp_level,
      wrk_indx,
      average_work,
      getWorkloadTolerance(level_number),
      min_size,
      cut_factor,
      bad_interval,
      physical_domain);

   tmp_level->deallocatePatchData(wrk_indx);
   tmp_level.reset();

   if (tmp_box_list.size() != static_cast<int>(tmp_work_list.size())) {
      TBOX_ERROR(
         d_object_name << ": "
                       << "Number of boxes generated != number of workload values generated."
                       << std::endl);
   }

   /*
    * Set local box array to list of chopped boxes and set local workload array.
    */

   hier::BoxContainer local_out_boxes(tmp_box_list);

   tbox::Array<double> local_out_workloads(local_out_boxes.size());

   int i = 0;
   for (std::list<double>::const_iterator il(tmp_work_list.begin());
        il != tmp_work_list.end(); il++) {
      local_out_workloads[i] = *il;
      i++;
   }

   /*
    * Gather local box and work arrays so that each processor has a copy.
    */
   exchangeBoxContainersAndWeightArrays(local_out_boxes,
      local_out_workloads,
      out_boxes,
      out_workloads,
      mpi);

}

/*
 *************************************************************************
 *
 * all-to-all exchange of box arrays and associated weights
 *
 *************************************************************************
 */
void
ChopAndPackLoadBalancer::exchangeBoxContainersAndWeightArrays(
   const hier::BoxContainer& box_list_in,
   const tbox::Array<double>& weights_in,
   hier::BoxContainer& box_list_out,
   tbox::Array<double>& weights_out,
   const tbox::SAMRAI_MPI& mpi) const
{
   TBOX_ASSERT(box_list_in.size() == weights_in.getSize());


   /*
    * allocate send and receive buffers, and set array sizes
    * for the output arrays.
    */
   int size_in = box_list_in.size();
   int size_out = size_in;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&size_in, &size_out, 1, MPI_INT, MPI_SUM);
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   if (size_out <= 0) {
      TBOX_ERROR("ChopAndPackLoadBalancer::exchangeBoxContainersAndWeightArrays() error"
         << "\n All box arrays have zero length!" << std::endl);
   }
#endif

   int buf_size_in = size_in * d_dim.getValue() * 2;
   int buf_size_out = size_out * d_dim.getValue() * 2;

   int curr_box_list_out_size = box_list_out.size();
   if (size_out > curr_box_list_out_size) {
      for (int i = curr_box_list_out_size; i < size_out; ++i) {
         box_list_out.pushBack(hier::Box(d_dim));
      }
   } else if (size_out < curr_box_list_out_size) {
      for (int i = size_out; i < curr_box_list_out_size; ++i) {
         box_list_out.popBack();
      }
   }
   weights_out.resizeArray(size_out);

   tbox::Array<int> buf_in(buf_size_in);
   tbox::Array<int> buf_out(buf_size_out);

   int* buf_in_ptr = (int *)NULL;
   int* buf_out_ptr = (int *)NULL;
   const double* wgts_in_ptr = (const double *)NULL;
   double* wgts_out_ptr = (double *)NULL;

   if (size_in > 0) {
      buf_in_ptr = buf_in.getPointer();
      wgts_in_ptr = weights_in.getPointer();
   }
   if (size_out > 0) {
      wgts_out_ptr = weights_out.getPointer();
      buf_out_ptr = buf_out.getPointer();
   }

   /*
    * populate the buffers with data for sending
    */
   int offset = 0;
   for (hier::BoxContainer::const_iterator x(box_list_in);
        x != box_list_in.end(); ++x) {
      for (int i = 0; i < d_dim.getValue(); ++i) {
         buf_in_ptr[offset++] = x->lower(i);
         buf_in_ptr[offset++] = x->upper(i);
      }
   }

   /*
    * exchange the data
    */
   std::vector<int> counts(mpi.getSize());
   mpi.Allgather(&size_in, 1, MPI_INT, &counts[0], 1, MPI_INT);
   std::vector<int> displs(mpi.getSize());
   displs[0] = 0;
   size_t total_count = counts[0];
   for (size_t i = 1; i < counts.size(); ++i) {
      displs[i] = displs[i - 1] + counts[i - 1];
      total_count += counts[i];
   }
   TBOX_ASSERT(static_cast<unsigned int>(weights_out.size()) == total_count);

   mpi.Allgatherv((void *)wgts_in_ptr,
      size_in,
      MPI_INT,
      wgts_out_ptr,
      &counts[0],
      &displs[0],
      MPI_INT);

   const int ints_per_box = d_dim.getValue() * 2;
   for (size_t i = 0; i < displs.size(); ++i) {
      counts[i] *= ints_per_box;
      displs[i] *= ints_per_box;
   }
   mpi.Allgatherv(buf_in_ptr, buf_size_in, MPI_INT, buf_out_ptr, &counts[0], &displs[0], MPI_INT);

   /*
    * assemble the output array of boxes
    */
   offset = 0;
   for (hier::BoxContainer::iterator b(box_list_out);
        b != box_list_out.end(); ++b) {
      for (int j = 0; j < d_dim.getValue(); ++j) {
         b->lower(j) = buf_out_ptr[offset++];
         b->upper(j) = buf_out_ptr[offset++];
      }
   }

}

/*
 *************************************************************************
 *
 * Print out all attributes of class instance for debugging.
 *
 *************************************************************************
 */

void
ChopAndPackLoadBalancer::printClassData(
   std::ostream& os) const
{
   os << "\nChopAndPackLoadBalancer::printClassData..." << std::endl;
   os << "\nChopAndPackLoadBalancer: this = "
      << (ChopAndPackLoadBalancer *)this << std::endl;
   os << "d_object_name = " << d_object_name << std::endl;
   os << "d_processor_layout_specified = "
      << d_processor_layout_specified << std::endl;
   os << "d_processor_layout = "
      << d_processor_layout << std::endl;
   os << "d_master_workload_data_id = "
      << d_master_workload_data_id << std::endl;
   os << "d_master_max_workload_factor = "
      << d_master_max_workload_factor << std::endl;
   os << "d_workload_tolerance = "
      << d_master_workload_tolerance << std::endl;
   os << "d_master_bin_pack_method = "
      << d_master_bin_pack_method << std::endl;

   int ln;

   os << "d_workload_data_id..." << std::endl;
   for (ln = 0; ln < d_workload_data_id.getSize(); ln++) {
      os << "    d_workload_data_id[" << ln << "] = "
         << d_workload_data_id[ln] << std::endl;
   }
   os << "d_max_workload_factor..." << std::endl;
   for (ln = 0; ln < d_max_workload_factor.getSize(); ln++) {
      os << "    d_max_workload_factor[" << ln << "] = "
         << d_max_workload_factor[ln] << std::endl;
   }
   os << "d_workload_tolerance..." << std::endl;
   for (ln = 0; ln < d_workload_tolerance.getSize(); ln++) {
      os << "    d_workload_tolerance[" << ln << "] = "
         << d_workload_tolerance[ln] << std::endl;
   }
   os << "d_bin_pack_method..." << std::endl;
   for (ln = 0; ln < d_bin_pack_method.getSize(); ln++) {
      os << "    d_bin_pack_method[" << ln << "] = "
         << d_bin_pack_method[ln] << std::endl;
   }
   os << "d_ignore_level_box_union_is_single_box = "
      << d_ignore_level_box_union_is_single_box << std::endl;

}

/*
 *************************************************************************
 *
 * Read values (described in the class header) from input database.
 *
 *************************************************************************
 */

void
ChopAndPackLoadBalancer::getFromInput(
   const boost::shared_ptr<tbox::Database>& db)
{

   if (db) {

      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

      d_master_bin_pack_method =
         db->getStringWithDefault("bin_pack_method",
            d_master_bin_pack_method);
      if (!(d_master_bin_pack_method == "GREEDY" ||
            d_master_bin_pack_method == "SPATIAL")) {
         TBOX_WARNING(
            d_object_name << ": "
                          << "Unknown 'bin_pack_method' "
                          << d_master_bin_pack_method
                          << " found in input. \nDefault 'GREEDY' will be used." << std::endl);
         d_master_bin_pack_method = "GREEDY";
      }

      if (db->keyExists("max_workload_factor")) {
         d_max_workload_factor = db->getDoubleArray("max_workload_factor");
         for (int i = 0; i < d_max_workload_factor.getSize(); i++) {
            if (d_max_workload_factor[i] < 0.0) {
               TBOX_ERROR(
                  d_object_name
                  << "Max workload values should be greater than 0");
            }
         }

         // Use last entry in array as value for finer levels
         d_master_max_workload_factor =
            d_max_workload_factor[d_max_workload_factor.getSize() - 1];
      }

      if (db->keyExists("workload_tolerance")) {
         d_workload_tolerance = db->getDoubleArray("workload_tolerance");
         for (int i = 0; i < d_workload_tolerance.getSize(); i++) {
            if (d_workload_tolerance[i] < 0.0 || d_workload_tolerance[i] >=
                1.0) {
               TBOX_ERROR(
                  d_object_name
                  << "Workload tolerance should be >= 0 and < 1.0");
            }
         }

         // Use last entry in array as value for finer levels
         d_master_workload_tolerance =
            d_workload_tolerance[d_workload_tolerance.getSize() - 1];
      }

      d_ignore_level_box_union_is_single_box =
         db->getBoolWithDefault("ignore_level_box_union_is_single_box",
            d_ignore_level_box_union_is_single_box);

      d_processor_layout_specified = false;
      int temp_processor_layout[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
      if (db->keyExists("processor_layout")) {
         db->getIntegerArray("processor_layout", temp_processor_layout, d_dim.getValue());

         /* consistency check */
         int totprocs = 1;
         for (int n = 0; n < d_dim.getValue(); n++) {
            totprocs *= temp_processor_layout[n];
         }

         if (totprocs != mpi.getSize()) {
            TBOX_WARNING(
               d_object_name << ": "
                             << "Input values for 'processor_layout' are inconsistent with"
                             << "\nnumber of processors.  Processor layout information will"
                             << "\nbe generated when needed."
                             << std::endl);
         } else {
            for (int n = 0; n < d_dim.getValue(); n++) {
               d_processor_layout(n) = temp_processor_layout[n];
            }
            d_processor_layout_specified = true;
         }
      }

      d_ignore_level_box_union_is_single_box =
         db->getBoolWithDefault("ignore_level_box_union_is_single_box",
            d_ignore_level_box_union_is_single_box);

   }

}

/*
 *************************************************************************
 *
 * Private utility function to map boxes to processors.
 * Note that no error checking is done.
 *
 *************************************************************************
 */

void
ChopAndPackLoadBalancer::binPackBoxes(
   hier::BoxContainer& boxes,
   hier::ProcessorMapping& mapping,
   tbox::Array<double>& workloads,
   const std::string& bin_pack_method) const
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   t_bin_pack_boxes->start();
   /*
    * Sort boxes in order of highest to lowest workload and assign
    * boxes to processors by bin packing.
    */
   t_bin_pack_boxes_sort->start();
   BalanceUtilities::sortDescendingBoxWorkloads(boxes,
      workloads);
   t_bin_pack_boxes_sort->stop();

   /*
    * Finally, assign boxes to processors by bin packing.
    */

   int num_procs = mpi.getSize();

   t_bin_pack_boxes_pack->start();
   if (bin_pack_method == "SPATIAL") {

      (void)BalanceUtilities::spatialBinPack(mapping,
         workloads,
         boxes,
         num_procs);

   } else if (bin_pack_method == "GREEDY") {

      (void)BalanceUtilities::binPack(mapping,
         workloads,
         num_procs);

   } else {

      TBOX_ERROR(
         d_object_name << ": "
                       << "Unknown bin pack method "
                       << bin_pack_method << " found." << std::endl);

   }
   t_bin_pack_boxes_pack->stop();

   t_bin_pack_boxes->stop();
}

/*
 *************************************************************************
 *************************************************************************
 */
void
ChopAndPackLoadBalancer::initializeCallback()
{
   t_load_balance_boxes = tbox::TimerManager::getManager()->
      getTimer("mesh::ChopAndPackLoadBalancer::loadBalanceBoxes()");
   t_load_balance_boxes_remove_intersection =
      tbox::TimerManager::getManager()->
      getTimer(
         "mesh::ChopAndPackLoadBalancer::loadBalanceBoxes()_remove_intersection");
   t_bin_pack_boxes = tbox::TimerManager::getManager()->
      getTimer("mesh::ChopAndPackLoadBalancer::binPackBoxes()");
   t_bin_pack_boxes_sort = tbox::TimerManager::getManager()->
      getTimer("mesh::ChopAndPackLoadBalancer::binPackBoxes()_sort");
   t_bin_pack_boxes_pack = tbox::TimerManager::getManager()->
      getTimer("mesh::ChopAndPackLoadBalancer::binPackBoxes()_pack");
   t_chop_boxes = tbox::TimerManager::getManager()->
      getTimer("mesh::ChopAndPackLoadBalancer::chop_boxes");
}

/*
 *************************************************************************
 *************************************************************************
 */
void
ChopAndPackLoadBalancer::finalizeCallback()
{
   t_load_balance_boxes.reset();
   t_load_balance_boxes_remove_intersection.reset();
   t_bin_pack_boxes.reset();
   t_bin_pack_boxes_sort.reset();
   t_bin_pack_boxes_pack.reset();
   t_chop_boxes.reset();
}

}
}
#endif
