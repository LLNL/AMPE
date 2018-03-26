/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Communication transaction for summing outeredge data
 *
 ************************************************************************/

#ifndef included_algs_OuteredgeSumTransaction_C
#define included_algs_OuteredgeSumTransaction_C

#include "SAMRAI/algs/OuteredgeSumTransaction.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/math/ArrayDataBasicOps.h"
#include "SAMRAI/pdat/EdgeGeometry.h"
#include "SAMRAI/pdat/OuteredgeData.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace algs {

/*
 *************************************************************************
 *
 * Initialization, set/unset functions for static array of refine items.
 *
 *************************************************************************
 */

const xfer::RefineClasses::Data **
OuteredgeSumTransaction::s_refine_items =
   (const xfer::RefineClasses::Data **)NULL;
int OuteredgeSumTransaction::s_num_refine_items = 0;

/*
 *************************************************************************
 *
 * Constructor sets state of transaction.
 *
 *************************************************************************
 */

OuteredgeSumTransaction::OuteredgeSumTransaction(
   const boost::shared_ptr<hier::PatchLevel>& dst_level,
   const boost::shared_ptr<hier::PatchLevel>& src_level,
   const boost::shared_ptr<hier::BoxOverlap>& overlap,
   const hier::Box& dst_node,
   const hier::Box& src_node,
   int refine_item_id):
   d_dst_level(dst_level),
   d_src_level(src_level),
   d_overlap(overlap),
   d_dst_node(dst_node),
   d_src_node(src_node),
   d_refine_item_id(refine_item_id),
   d_incoming_bytes(0),
   d_outgoing_bytes(0)
{
   TBOX_ASSERT(dst_level);
   TBOX_ASSERT(src_level);
   TBOX_ASSERT(overlap);
   TBOX_ASSERT(dst_node.getLocalId() >= 0);
   TBOX_ASSERT(src_node.getLocalId() >= 0);
   TBOX_ASSERT(refine_item_id >= 0);
   // Note: s_num_refine_items cannot be used at this point!

   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst_level, *src_level, dst_node, src_node);
}

OuteredgeSumTransaction::~OuteredgeSumTransaction()
{
}

/*
 *************************************************************************
 *
 * Functions overridden in tbox::Transaction base class.
 *
 *************************************************************************
 */

bool
OuteredgeSumTransaction::canEstimateIncomingMessageSize()
{
   bool can_estimate = false;
   if (d_src_node.getOwnerRank() == d_src_level->getBoxLevel()->getMPI().getRank()) {
      can_estimate =
         d_src_level->getPatch(d_src_node.getGlobalId())
         ->getPatchData(s_refine_items[d_refine_item_id]->
            d_src)
         ->canEstimateStreamSizeFromBox();
   } else {
      can_estimate =
         d_dst_level->getPatch(d_dst_node.getGlobalId())
         ->getPatchData(s_refine_items[d_refine_item_id]->
            d_scratch)
         ->canEstimateStreamSizeFromBox();
   }
   return can_estimate;
}

size_t
OuteredgeSumTransaction::computeIncomingMessageSize()
{
   d_incoming_bytes =
      d_dst_level->getPatch(d_dst_node.getGlobalId())
      ->getPatchData(s_refine_items[d_refine_item_id]->
         d_scratch)
      ->getDataStreamSize(*d_overlap);
   return d_incoming_bytes;
}

size_t
OuteredgeSumTransaction::computeOutgoingMessageSize()
{
   d_outgoing_bytes =
      d_src_level->getPatch(d_src_node.getGlobalId())
      ->getPatchData(s_refine_items[d_refine_item_id]->
         d_src)
      ->getDataStreamSize(*d_overlap);
   return d_outgoing_bytes;
}

int
OuteredgeSumTransaction::getSourceProcessor()
{
   return d_src_node.getOwnerRank();
}

int
OuteredgeSumTransaction::getDestinationProcessor()
{
   return d_dst_node.getOwnerRank();
}

void
OuteredgeSumTransaction::packStream(
   tbox::MessageStream& stream)
{
   d_src_level->getPatch(d_src_node.getGlobalId())
   ->getPatchData(s_refine_items[d_refine_item_id]->
      d_src)
   ->packStream(stream, *d_overlap);
}

void
OuteredgeSumTransaction::unpackStream(
   tbox::MessageStream& stream)
{
   boost::shared_ptr<pdat::OuteredgeData<double> > oedge_dst_data(
      d_dst_level->getPatch(d_dst_node.getGlobalId())->
      getPatchData(s_refine_items[d_refine_item_id]->d_scratch),
      boost::detail::dynamic_cast_tag());
   TBOX_ASSERT(oedge_dst_data);

   oedge_dst_data->unpackStreamAndSum(stream, *d_overlap);
}

void
OuteredgeSumTransaction::copyLocalData()
{
   boost::shared_ptr<pdat::OuteredgeData<double> > oedge_dst_data(
      d_dst_level->getPatch(d_dst_node.getGlobalId())->
      getPatchData(s_refine_items[d_refine_item_id]->d_scratch),
      boost::detail::dynamic_cast_tag());
   TBOX_ASSERT(oedge_dst_data);

   boost::shared_ptr<pdat::OuteredgeData<double> > oedge_src_data(
      d_src_level->getPatch(d_src_node.getGlobalId())->
      getPatchData(s_refine_items[d_refine_item_id]->d_src),
      boost::detail::dynamic_cast_tag());
   TBOX_ASSERT(oedge_src_data);

   oedge_dst_data->sum(*oedge_src_data, *d_overlap);
}

/*
 *************************************************************************
 *
 * Function to print state of transaction.
 *
 *************************************************************************
 */

void
OuteredgeSumTransaction::printClassData(
   std::ostream& stream) const
{
   stream << "Outeredge Sum Transaction" << std::endl;
   stream << "   refine item array:        "
          << (xfer::RefineClasses::Data **)s_refine_items
          << std::endl;
   stream << "   num refine items:       " << s_num_refine_items << std::endl;
   stream << "   destination node:       " << d_dst_node << std::endl;
   stream << "   source node:            " << d_src_node << std::endl;
   stream << "   refine item id:         " << d_refine_item_id << std::endl;
   stream << "   destination patch data: "
          << s_refine_items[d_refine_item_id]->d_scratch << std::endl;
   stream << "   source patch data:      "
          << s_refine_items[d_refine_item_id]->d_src << std::endl;
   stream << "   incoming bytes:         " << d_incoming_bytes << std::endl;
   stream << "   outgoing bytes:         " << d_outgoing_bytes << std::endl;
   stream << "   destination level:           "
          << d_dst_level.get() << std::endl;
   stream << "   source level:           "
          << d_src_level.get() << std::endl;
   stream << "   overlap:                " << std::endl;
   d_overlap->print(stream);
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
