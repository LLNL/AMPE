/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Concrete factory to create standard copy and time transactions
 *                for refine schedules.
 *
 ************************************************************************/

#ifndef included_xfer_StandardRefineTransactionFactory_C
#define included_xfer_StandardRefineTransactionFactory_C

#include "SAMRAI/xfer/StandardRefineTransactionFactory.h"

#include "SAMRAI/xfer/RefineCopyTransaction.h"
#include "SAMRAI/xfer/RefineTimeTransaction.h"

#include <boost/make_shared.hpp>

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Default constructor and destructor.
 *
 *************************************************************************
 */

StandardRefineTransactionFactory::StandardRefineTransactionFactory()
{
}

StandardRefineTransactionFactory::~StandardRefineTransactionFactory()
{
}

/*
 *************************************************************************
 *
 * Set/unset information for transactions managed by this factory class.
 *
 *************************************************************************
 */

void
StandardRefineTransactionFactory::setRefineItems(
   const RefineClasses::Data** refine_items,
   int num_refine_items)
{
   RefineCopyTransaction::setRefineItems(refine_items, num_refine_items);
   RefineTimeTransaction::setRefineItems(refine_items, num_refine_items);
   d_refine_items = refine_items;
   d_num_refine_items = num_refine_items;
}

void
StandardRefineTransactionFactory::unsetRefineItems()
{
   RefineCopyTransaction::unsetRefineItems();
   RefineTimeTransaction::unsetRefineItems();
   d_refine_items = (const RefineClasses::Data **)NULL;
   d_num_refine_items = 0;
}

void
StandardRefineTransactionFactory::setTransactionTime(
   double fill_time)
{
   RefineTimeTransaction::setTransactionTime(fill_time);
}

/*
 *************************************************************************
 *
 * Allocate appropriate transaction object.
 *
 *************************************************************************
 */

boost::shared_ptr<tbox::Transaction>
StandardRefineTransactionFactory::allocate(
   const boost::shared_ptr<hier::PatchLevel>& dst_level,
   const boost::shared_ptr<hier::PatchLevel>& src_level,
   const boost::shared_ptr<hier::BoxOverlap>& overlap,
   const hier::Box& dst_mapped_box,
   const hier::Box& src_mapped_box,
   int ritem_id,
   const hier::Box& box,
   bool use_time_interpolation) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS5(*dst_level,
      *src_level,
      dst_mapped_box,
      src_mapped_box,
      box);

   if (use_time_interpolation) {

      return boost::make_shared<RefineTimeTransaction>(
         dst_level,
         src_level,
         overlap,
         dst_mapped_box,
         src_mapped_box,
         box,
         ritem_id);

   } else {

      return boost::make_shared<RefineCopyTransaction>(
         dst_level,
         src_level,
         overlap,
         dst_mapped_box,
         src_mapped_box,
         ritem_id);

   }
}

void
StandardRefineTransactionFactory::preprocessScratchSpace(
   const boost::shared_ptr<hier::PatchLevel>& level,
   double fill_time,
   const hier::ComponentSelector& preprocess_vector) const
{
   NULL_USE(level);
   NULL_USE(fill_time);
   NULL_USE(preprocess_vector);
}

}
}
#endif
