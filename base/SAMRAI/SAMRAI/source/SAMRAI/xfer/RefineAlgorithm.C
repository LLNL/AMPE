/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Refine algorithm for data transfer between AMR levels
 *
 ************************************************************************/

#ifndef included_xfer_RefineAlgorithm_C
#define included_xfer_RefineAlgorithm_C

#include "SAMRAI/xfer/RefineAlgorithm.h"

#include "SAMRAI/xfer/BoxGeometryVariableFillPattern.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"
#include "SAMRAI/xfer/StandardRefineTransactionFactory.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Utilities.h"

#include <boost/make_shared.hpp>

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Default constructor creates a new RefineClasses object.
 *
 *************************************************************************
 */

RefineAlgorithm::RefineAlgorithm(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_refine_classes(boost::make_shared<RefineClasses>()),
   d_schedule_created(false)
{
}

/*
 *************************************************************************
 *
 * The destructor implicitly deletes the list storage associated with
 * the refine algorithm.
 *
 *************************************************************************
 */

RefineAlgorithm::~RefineAlgorithm()
{
}

/*
 *************************************************************************
 *
 * Register a refine operation that will not require time interpolation.
 *
 *************************************************************************
 */

void
RefineAlgorithm::registerRefine(
   const int dst,
   const int src,
   const int scratch,
   const boost::shared_ptr<hier::RefineOperator>& oprefine,
   const boost::shared_ptr<VariableFillPattern>& var_fill_pattern)
{
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   if (oprefine) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *oprefine);
   }
#endif

   if (d_schedule_created) {
      TBOX_ERROR("RefineAlgorithm::registerRefine error..."
         << "\nCannot call registerRefine with a RefineAlgorithm"
         << "\nobject that has already been used to create a schedule."
         << std::endl);
   }

   RefineClasses::Data data;

   data.d_dst = dst;
   data.d_src = src;
   data.d_src_told = -1;
   data.d_src_tnew = -1;
   data.d_scratch = scratch;
   data.d_fine_bdry_reps_var = hier::VariableDatabase::getDatabase()->
      getPatchDescriptor()->getPatchDataFactory(dst)->
      fineBoundaryRepresentsVariable();
   data.d_time_interpolate = false;
   data.d_oprefine = oprefine;
   data.d_optime.reset();
   data.d_tag = -1;
   if (var_fill_pattern) {
      data.d_var_fill_pattern = var_fill_pattern;
   } else {
      data.d_var_fill_pattern.reset(new BoxGeometryVariableFillPattern());
   }

   d_refine_classes->insertEquivalenceClassItem(data);
}

/*
 *************************************************************************
 *
 * Register a refine operation that will require time interpolation.
 *
 *************************************************************************
 */

void
RefineAlgorithm::registerRefine(
   const int dst,
   const int src,
   const int src_told,
   const int src_tnew,
   const int scratch,
   const boost::shared_ptr<hier::RefineOperator>& oprefine,
   const boost::shared_ptr<hier::TimeInterpolateOperator>& optime,
   const boost::shared_ptr<VariableFillPattern>& var_fill_pattern)
{
   TBOX_ASSERT(optime);

   if (d_schedule_created) {
      TBOX_ERROR("RefineAlgorithm::registerRefine error..."
         << "\nCannot call registerRefine with a RefineAlgorithm object"
         << "\nthat has already been used to create a schedule."
         << std::endl);
   }

   RefineClasses::Data data;

   data.d_dst = dst;
   data.d_src = src;
   data.d_src_told = src_told;
   data.d_src_tnew = src_tnew;
   data.d_scratch = scratch;
   data.d_fine_bdry_reps_var = hier::VariableDatabase::getDatabase()->
      getPatchDescriptor()->getPatchDataFactory(dst)->
      fineBoundaryRepresentsVariable();
   data.d_time_interpolate = true;
   data.d_oprefine = oprefine;
   data.d_optime = optime;
   data.d_tag = -1;
   if (var_fill_pattern) {
      data.d_var_fill_pattern = var_fill_pattern;
   } else {
      data.d_var_fill_pattern.reset(new BoxGeometryVariableFillPattern());
   }

   d_refine_classes->insertEquivalenceClassItem(data);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that will copy data from the
 * interiors of the specified level into the ghost cells and
 * interiors of the same level.
 *
 *************************************************************************
 */

boost::shared_ptr<RefineSchedule>
RefineAlgorithm::createSchedule(
   const boost::shared_ptr<hier::PatchLevel>& level,
   RefinePatchStrategy* patch_strategy,
   const boost::shared_ptr<RefineTransactionFactory>& transaction_factory)
{
   TBOX_ASSERT(level);
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *level);
   if (patch_strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *patch_strategy);
   }
#endif

   d_schedule_created = true;

   boost::shared_ptr<RefineTransactionFactory> trans_factory(
      transaction_factory);

   if (!trans_factory) {
      trans_factory.reset(new StandardRefineTransactionFactory);
   }

   boost::shared_ptr<PatchLevelFullFillPattern> fill_pattern(
      boost::make_shared<PatchLevelFullFillPattern>());

   return boost::make_shared<RefineSchedule>(
      fill_pattern,
      level,
      level,
      d_refine_classes,
      trans_factory,
      patch_strategy);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that will copy data from the
 * interiors of the specified level into the ghost cells and
 * interiors of the same level.
 *
 *************************************************************************
 */

boost::shared_ptr<RefineSchedule>
RefineAlgorithm::createSchedule(
   const boost::shared_ptr<PatchLevelFillPattern>& fill_pattern,
   const boost::shared_ptr<hier::PatchLevel>& level,
   RefinePatchStrategy* patch_strategy,
   const boost::shared_ptr<RefineTransactionFactory>& transaction_factory)
{
   TBOX_ASSERT(level);
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *level);
   if (patch_strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *patch_strategy);
   }
#endif

   d_schedule_created = true;

   boost::shared_ptr<RefineTransactionFactory> trans_factory(
      transaction_factory);

   if (!trans_factory) {
      trans_factory.reset(new StandardRefineTransactionFactory);
   }

   return boost::make_shared<RefineSchedule>(
      fill_pattern,
      level,
      level,
      d_refine_classes,
      trans_factory,
      patch_strategy);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that will copy data from the
 * interiors of the source level into the ghost cell and interiors
 * of the destination level.
 *
 *************************************************************************
 */

boost::shared_ptr<RefineSchedule>
RefineAlgorithm::createSchedule(
   const boost::shared_ptr<hier::PatchLevel>& dst_level,
   const boost::shared_ptr<hier::PatchLevel>& src_level,
   RefinePatchStrategy* patch_strategy,
   bool use_time_refinement,
   const boost::shared_ptr<RefineTransactionFactory>& transaction_factory)
{
   // TBOX_ERROR("Untried method!  I think this method should work, but it's never been excercised.  When code crashes here, remove this line and rerun.  If problem continues, it could well be due to excercising this code.  --BTNG");

   TBOX_ASSERT(dst_level);
   TBOX_ASSERT(src_level);
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, *dst_level, *src_level);
   if (patch_strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *patch_strategy);
   }
#endif

   d_schedule_created = true;

   boost::shared_ptr<RefineTransactionFactory> trans_factory(
      transaction_factory);

   if (!trans_factory) {
      trans_factory.reset(new StandardRefineTransactionFactory);
   }

   boost::shared_ptr<PatchLevelFullFillPattern> fill_pattern(
      boost::make_shared<PatchLevelFullFillPattern>());

   return boost::make_shared<RefineSchedule>(
      fill_pattern,
      dst_level,
      src_level,
      d_refine_classes,
      trans_factory,
      patch_strategy,
      use_time_refinement);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that will copy data from the
 * interiors of the source level into the ghost cell and interiors
 * of the destination level.
 *
 *************************************************************************
 */

boost::shared_ptr<RefineSchedule>
RefineAlgorithm::createSchedule(
   const boost::shared_ptr<PatchLevelFillPattern>& fill_pattern,
   const boost::shared_ptr<hier::PatchLevel>& dst_level,
   const boost::shared_ptr<hier::PatchLevel>& src_level,
   RefinePatchStrategy* patch_strategy,
   bool use_time_refinement,
   const boost::shared_ptr<RefineTransactionFactory>& transaction_factory)
{
   TBOX_ASSERT(dst_level);
   TBOX_ASSERT(src_level);
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS3(*this, *dst_level, *src_level);
   if (patch_strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *patch_strategy);
   }
#endif

   d_schedule_created = true;

   boost::shared_ptr<RefineTransactionFactory> trans_factory(
      transaction_factory);

   if (!trans_factory) {
      trans_factory.reset(new StandardRefineTransactionFactory);
   }

   return boost::make_shared<RefineSchedule>(
      fill_pattern,
      dst_level,
      src_level,
      d_refine_classes,
      trans_factory,
      patch_strategy,
      use_time_refinement);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that copies data from the interiors
 * of the same level and coarser levels into the interior and boundary
 * cells of the given level.
 *
 *************************************************************************
 */

boost::shared_ptr<RefineSchedule>
RefineAlgorithm::createSchedule(
   const boost::shared_ptr<hier::PatchLevel>& level,
   const int next_coarser_level,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   RefinePatchStrategy* patch_strategy,
   bool use_time_refinement,
   const boost::shared_ptr<RefineTransactionFactory>& transaction_factory)
{

   // Do we all agree on the destination mapped_box_level?
   TBOX_ASSERT(level);
   TBOX_ASSERT((next_coarser_level == -1) || hierarchy);
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *level);
   if (hierarchy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *hierarchy);
   }
   if (patch_strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *patch_strategy);
   }
#endif

   d_schedule_created = true;

   boost::shared_ptr<RefineTransactionFactory> trans_factory(
      transaction_factory);

   if (!trans_factory) {
      trans_factory.reset(new StandardRefineTransactionFactory);
   }

   boost::shared_ptr<PatchLevelFullFillPattern> fill_pattern(
      boost::make_shared<PatchLevelFullFillPattern>());

   return boost::make_shared<RefineSchedule>(
      fill_pattern,
      level,
      level,
      next_coarser_level,
      hierarchy,
      d_refine_classes,
      trans_factory,
      patch_strategy,
      use_time_refinement);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that copies data from the interiors
 * of the same level and coarser levels into the interior and boundary
 * cells of the given level.
 *
 *************************************************************************
 */

boost::shared_ptr<RefineSchedule>
RefineAlgorithm::createSchedule(
   const boost::shared_ptr<PatchLevelFillPattern>& fill_pattern,
   const boost::shared_ptr<hier::PatchLevel>& level,
   const int next_coarser_level,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   RefinePatchStrategy* patch_strategy,
   bool use_time_refinement,
   const boost::shared_ptr<RefineTransactionFactory>& transaction_factory)
{

   // Do we all agree on the destination mapped_box_level?
   TBOX_ASSERT(level);
   TBOX_ASSERT((next_coarser_level == -1) || hierarchy);
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *level);
   if (hierarchy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *hierarchy);
   }
   if (patch_strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *patch_strategy);
   }
#endif

   d_schedule_created = true;

   boost::shared_ptr<RefineTransactionFactory> trans_factory(
      transaction_factory);

   if (!trans_factory) {
      trans_factory.reset(new StandardRefineTransactionFactory);
   }

   return boost::make_shared<RefineSchedule>(
      fill_pattern,
      level,
      level,
      next_coarser_level,
      hierarchy,
      d_refine_classes,
      trans_factory,
      patch_strategy,
      use_time_refinement);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that copies data from the interiors
 * of the old level and coarser levels into the ghost cells and interior
 * cells of the given new level.
 *
 *************************************************************************
 */

boost::shared_ptr<RefineSchedule>
RefineAlgorithm::createSchedule(
   const boost::shared_ptr<hier::PatchLevel>& dst_level,
   const boost::shared_ptr<hier::PatchLevel>& src_level,
   const int next_coarser_level,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   RefinePatchStrategy* patch_strategy,
   bool use_time_refinement,
   const boost::shared_ptr<RefineTransactionFactory>& transaction_factory)
{
   NULL_USE(use_time_refinement);

   TBOX_ASSERT(dst_level);
   TBOX_ASSERT((next_coarser_level == -1) || hierarchy);
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *dst_level);
   if (src_level) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *src_level);
   }
   if (hierarchy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *hierarchy);
   }
   if (patch_strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *patch_strategy);
   }
#endif

   // Do we all agree on the destination mapped_box_level?
   if (src_level) {
      if (next_coarser_level >= 0) {
      }
   }

   d_schedule_created = true;

   boost::shared_ptr<RefineTransactionFactory> trans_factory(
      transaction_factory);

   if (!trans_factory) {
      trans_factory.reset(new StandardRefineTransactionFactory);
   }

   boost::shared_ptr<PatchLevelFullFillPattern> fill_pattern(
      boost::make_shared<PatchLevelFullFillPattern>());

   return boost::make_shared<RefineSchedule>(
      fill_pattern,
      dst_level,
      src_level,
      next_coarser_level,
      hierarchy,
      d_refine_classes,
      trans_factory,
      patch_strategy,
      false);
}

/*
 *************************************************************************
 *
 * Create a communication schedule that copies data from the interiors
 * of the old level and coarser levels into the ghost cells and interior
 * cells of the given new level.
 *
 *************************************************************************
 */

boost::shared_ptr<RefineSchedule>
RefineAlgorithm::createSchedule(
   const boost::shared_ptr<PatchLevelFillPattern>& fill_pattern,
   const boost::shared_ptr<hier::PatchLevel>& dst_level,
   const boost::shared_ptr<hier::PatchLevel>& src_level,
   const int next_coarser_level,
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   RefinePatchStrategy* patch_strategy,
   bool use_time_refinement,
   const boost::shared_ptr<RefineTransactionFactory>& transaction_factory)
{
   NULL_USE(use_time_refinement);

   TBOX_ASSERT(dst_level);
   TBOX_ASSERT((next_coarser_level == -1) || hierarchy);
#ifdef DEBUG_CHECK_DIM_ASSERTIONS
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *dst_level);
   if (src_level) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *src_level);
   }
   if (hierarchy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *hierarchy);
   }
   if (patch_strategy) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *patch_strategy);
   }
#endif

   // Do we all agree on the destination mapped_box_level?
   if (src_level) {
      if (next_coarser_level >= 0) {
      }
   }

   d_schedule_created = true;

   boost::shared_ptr<RefineTransactionFactory> trans_factory(
      transaction_factory);

   if (!trans_factory) {
      trans_factory.reset(new StandardRefineTransactionFactory);
   }

   return boost::make_shared<RefineSchedule>(
      fill_pattern,
      dst_level,
      src_level,
      next_coarser_level,
      hierarchy,
      d_refine_classes,
      trans_factory,
      patch_strategy,
      false);
}

/*
 **************************************************************************
 *
 * Reconfigure refine schedule to perform operations in this algorithm.
 *
 **************************************************************************
 */

bool
RefineAlgorithm::checkConsistency(
   const boost::shared_ptr<RefineSchedule>& schedule) const
{
   TBOX_ASSERT(schedule);
   return d_refine_classes->classesMatch(schedule->getEquivalenceClasses());
}

void RefineAlgorithm::resetSchedule(
   const boost::shared_ptr<RefineSchedule>& schedule) const
{
   TBOX_ASSERT(schedule);
   if (d_refine_classes->classesMatch(schedule->getEquivalenceClasses())) {
      schedule->reset(d_refine_classes);
   } else {
      TBOX_ERROR("RefineAlgorithm::resetSchedule error..."
         << "\n Items in xfer::RefineClasses object passed to reset"
         << "\n routine does not match that in existing schedule."
         << std::endl);
   }
}

/*
 *************************************************************************
 *
 * Print refine algorithm data to the specified output stream.
 *
 *************************************************************************
 */

void
RefineAlgorithm::printClassData(
   std::ostream& stream) const
{
   stream << "RefineAlgorithm::printClassData()" << std::endl;
   stream << "----------------------------------------" << std::endl;
   d_refine_classes->printClassData(stream);
}

}
}
#endif
