/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Concrete factory to create standard copy and time transactions
 *                for refine schedules.
 *
 ************************************************************************/
#include "SAMRAI/xfer/StandardRefineTransactionFactory.h"

#include "SAMRAI/xfer/RefineCopyTransaction.h"
#include "SAMRAI/xfer/RefineTimeTransaction.h"

#include "boost/make_shared.hpp"

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
   const hier::Box& dst_box,
   const hier::Box& src_box,
   const RefineClasses::Data** refine_data,
   int item_id,
   const hier::Box& box,
   bool use_time_interpolation) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY5(*dst_level,
      *src_level,
      dst_box,
      src_box,
      box);

   if (use_time_interpolation) {

      return boost::make_shared<RefineTimeTransaction>(
                dst_level,
                src_level,
                overlap,
                dst_box,
                src_box,
                box,
                refine_data,
                item_id);

   } else {

      return boost::make_shared<RefineCopyTransaction>(
                dst_level,
                src_level,
                overlap,
                dst_box,
                src_box,
                refine_data,
                item_id);

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
