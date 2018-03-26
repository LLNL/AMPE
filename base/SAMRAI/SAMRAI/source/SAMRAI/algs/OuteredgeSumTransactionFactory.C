/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Factory for creating outeredge sum transaction objects
 *
 ************************************************************************/

#ifndef included_algs_OuteredgeSumTransactionFactory_C
#define included_algs_OuteredgeSumTransactionFactory_C

#include "SAMRAI/algs/OuteredgeSumTransactionFactory.h"

#include "SAMRAI/pdat/OuteredgeData.h"
#include "SAMRAI/algs/OuteredgeSumTransaction.h"

#include <boost/make_shared.hpp>

namespace SAMRAI {
namespace algs {

/*
 *************************************************************************
 *
 * Default constructor and destructor.
 *
 *************************************************************************
 */

OuteredgeSumTransactionFactory::OuteredgeSumTransactionFactory()
{
}

OuteredgeSumTransactionFactory::~OuteredgeSumTransactionFactory()
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
OuteredgeSumTransactionFactory::setRefineItems(
   const xfer::RefineClasses::Data** refine_items,
   int num_refine_items)
{
   OuteredgeSumTransaction::setRefineItems(refine_items,
      num_refine_items);
   d_refine_items = refine_items;
   d_number_refine_items = num_refine_items;
}

void
OuteredgeSumTransactionFactory::unsetRefineItems()
{
   OuteredgeSumTransaction::unsetRefineItems();
   d_refine_items = (const xfer::RefineClasses::Data **)NULL;
   d_number_refine_items = 0;
}

/*
 *************************************************************************
 *
 * Allocate outeredge sum transaction object.
 *
 *************************************************************************
 */

boost::shared_ptr<tbox::Transaction>
OuteredgeSumTransactionFactory::allocate(
   const boost::shared_ptr<hier::PatchLevel>& dst_level,
   const boost::shared_ptr<hier::PatchLevel>& src_level,
   const boost::shared_ptr<hier::BoxOverlap>& overlap,
   const hier::Box& dst_node,
   const hier::Box& src_node,
   int ritem_id,
   const hier::Box& box,
   bool use_time_interpolation) const
{
   NULL_USE(box);
   NULL_USE(use_time_interpolation);

   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst_level, *src_level, dst_node, src_node);

   return boost::make_shared<OuteredgeSumTransaction>(dst_level,
      src_level,
      overlap,
      dst_node,
      src_node,
      ritem_id);
}

boost::shared_ptr<tbox::Transaction>
OuteredgeSumTransactionFactory::allocate(
   const boost::shared_ptr<hier::PatchLevel>& dst_level,
   const boost::shared_ptr<hier::PatchLevel>& src_level,
   const boost::shared_ptr<hier::BoxOverlap>& overlap,
   const hier::Box& dst_node,
   const hier::Box& src_node,
   int ritem_id) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst_level, *src_level, dst_node, src_node);

   return allocate(dst_level,
      src_level,
      overlap,
      dst_node,
      src_node,
      ritem_id,
      hier::Box::getEmptyBox(dst_level->getDim()),
      false);
}

/*
 *************************************************************************
 *
 * Initialize (to 0.0) scratch storage for sum transactions.
 *
 *************************************************************************
 */

void
OuteredgeSumTransactionFactory::preprocessScratchSpace(
   const boost::shared_ptr<hier::PatchLevel>& level,
   double fill_time,
   const hier::ComponentSelector& preprocess_vector) const
{
   NULL_USE(fill_time);
   TBOX_ASSERT(level);

   for (hier::PatchLevel::iterator ip(level->begin());
        ip != level->end(); ++ip) {
      const boost::shared_ptr<hier::Patch>& patch = *ip;

      const int ncomponents = preprocess_vector.getSize();
      for (int n = 0; n < ncomponents; ++n) {
         if (preprocess_vector.isSet(n)) {
            boost::shared_ptr<pdat::OuteredgeData<double> > oedge_data(
               patch->getPatchData(n),
               boost::detail::dynamic_cast_tag());
            oedge_data->fillAll(0.0);
         }
      }

   }
}

}
}
#endif
