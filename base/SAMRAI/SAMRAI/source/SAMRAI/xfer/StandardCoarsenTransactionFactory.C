/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Concrete factory for create standard copy transactions
 *                for coarsen schedules.
 *
 ************************************************************************/

#ifndef included_xfer_StandardCoarsenTransactionFactory_C
#define included_xfer_StandardCoarsenTransactionFactory_C

#include "SAMRAI/xfer/StandardCoarsenTransactionFactory.h"

#include "SAMRAI/xfer/CoarsenCopyTransaction.h"

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

StandardCoarsenTransactionFactory::StandardCoarsenTransactionFactory()
{
}

StandardCoarsenTransactionFactory::~StandardCoarsenTransactionFactory()
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
StandardCoarsenTransactionFactory::setCoarsenItems(
   const CoarsenClasses::Data** coarsen_items,
   int num_coarsen_items)
{
   CoarsenCopyTransaction::setCoarsenItems(coarsen_items,
      num_coarsen_items);
   d_coarsen_items = coarsen_items;
   d_num_coarsen_items = num_coarsen_items;
}

void
StandardCoarsenTransactionFactory::unsetCoarsenItems()
{
   CoarsenCopyTransaction::unsetCoarsenItems();
   d_coarsen_items = (const CoarsenClasses::Data **)NULL;
   d_num_coarsen_items = 0;
}

/*
 *************************************************************************
 *
 * Allocate appropriate transaction object.
 *
 *************************************************************************
 */

boost::shared_ptr<tbox::Transaction>
StandardCoarsenTransactionFactory::allocate(
   const boost::shared_ptr<hier::PatchLevel>& dst_level,
   const boost::shared_ptr<hier::PatchLevel>& src_level,
   const boost::shared_ptr<hier::BoxOverlap>& overlap,
   const hier::Box& dst_mapped_box,
   const hier::Box& src_mapped_box,
   int citem_id) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst_level,
      *src_level,
      dst_mapped_box,
      src_mapped_box);

   return boost::make_shared<CoarsenCopyTransaction>(
      dst_level,
      src_level,
      overlap,
      dst_mapped_box,
      src_mapped_box,
      citem_id);
}

}
}
#endif
