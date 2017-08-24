/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Factory for creating outeredge sum transaction objects
 *
 ************************************************************************/
#include "SAMRAI/algs/OuteredgeSumTransactionFactory.h"

#include "SAMRAI/pdat/OuteredgeData.h"
#include "SAMRAI/algs/OuteredgeSumTransaction.h"

#include "boost/make_shared.hpp"

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
   const xfer::RefineClasses::Data** refine_data,
   int item_id,
   const hier::Box& box,
   bool use_time_interpolation) const
{
   NULL_USE(box);
   NULL_USE(use_time_interpolation);

   TBOX_ASSERT(dst_level);
   TBOX_ASSERT(src_level);
   TBOX_ASSERT(overlap);
   TBOX_ASSERT(dst_node.getLocalId() >= 0);
   TBOX_ASSERT(src_node.getLocalId() >= 0);
   TBOX_ASSERT(item_id >= 0);
   TBOX_ASSERT(refine_data[item_id] != 0);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst_level, *src_level, dst_node, src_node);

   return boost::make_shared<OuteredgeSumTransaction>(dst_level,
                                                      src_level,
                                                      overlap,
                                                      dst_node,
                                                      src_node,
                                                      refine_data,
                                                      item_id);
}

boost::shared_ptr<tbox::Transaction>
OuteredgeSumTransactionFactory::allocate(
   const boost::shared_ptr<hier::PatchLevel>& dst_level,
   const boost::shared_ptr<hier::PatchLevel>& src_level,
   const boost::shared_ptr<hier::BoxOverlap>& overlap,
   const hier::Box& dst_node,
   const hier::Box& src_node,
   const xfer::RefineClasses::Data** refine_data,
   int item_id) const
{
   TBOX_ASSERT(dst_level);
   TBOX_ASSERT(src_level);
   TBOX_ASSERT(overlap);
   TBOX_ASSERT(dst_node.getLocalId() >= 0);
   TBOX_ASSERT(src_node.getLocalId() >= 0);
   TBOX_ASSERT(item_id >= 0);
   TBOX_ASSERT(refine_data[item_id] != 0);
   TBOX_ASSERT_OBJDIM_EQUALITY4(*dst_level, *src_level, dst_node, src_node);

   return allocate(dst_level,
      src_level,
      overlap,
      dst_node,
      src_node,
      refine_data,
      item_id,
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
               BOOST_CAST<pdat::OuteredgeData<double>, hier::PatchData>(
                  patch->getPatchData(n)));
            TBOX_ASSERT(oedge_data);
            oedge_data->fillAll(0.0);
         }
      }

   }
}

}
}
