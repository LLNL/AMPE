/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Concrete factory for create standard copy transactions
 *                for coarsen schedules.
 *
 ************************************************************************/
#include "SAMRAI/xfer/StandardCoarsenTransactionFactory.h"

#include "SAMRAI/xfer/CoarsenCopyTransaction.h"

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

StandardCoarsenTransactionFactory::StandardCoarsenTransactionFactory()
{
}

StandardCoarsenTransactionFactory::~StandardCoarsenTransactionFactory()
{
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
   const hier::Box& dst_box,
   const hier::Box& src_box,
   const CoarsenClasses::Data** coarsen_data,
   int item_id) const
{
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst_level,
      *src_level,
      dst_box,
      src_box);

   return boost::make_shared<CoarsenCopyTransaction>(
             dst_level,
             src_level,
             overlap,
             dst_box,
             src_box,
             coarsen_data,
             item_id);
}

}
}
