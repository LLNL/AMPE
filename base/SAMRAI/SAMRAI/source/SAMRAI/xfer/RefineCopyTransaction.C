/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Communication transaction for data copies during data refining
 *
 ************************************************************************/

#ifndef included_xfer_RefineCopyTransaction_C
#define included_xfer_RefineCopyTransaction_C

#include "SAMRAI/xfer/RefineCopyTransaction.h"

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace xfer {

const RefineClasses::Data ** RefineCopyTransaction::s_refine_items =
   (const RefineClasses::Data **)NULL;
int RefineCopyTransaction::s_num_refine_items = 0;

/*
 *************************************************************************
 *
 * Constructor sets state of transaction.
 *
 *************************************************************************
 */

RefineCopyTransaction::RefineCopyTransaction(
   const boost::shared_ptr<hier::PatchLevel>& dst_level,
   const boost::shared_ptr<hier::PatchLevel>& src_level,
   const boost::shared_ptr<hier::BoxOverlap>& overlap,
   const hier::Box& dst_mapped_box,
   const hier::Box& src_mapped_box,
   const int refine_item_id):
   d_dst_patch_rank(dst_mapped_box.getOwnerRank()),
   d_src_patch_rank(src_mapped_box.getOwnerRank()),
   d_overlap(overlap),
   d_refine_item_id(refine_item_id),
   d_incoming_bytes(0),
   d_outgoing_bytes(0)
{
   TBOX_ASSERT(dst_level);
   TBOX_ASSERT(src_level);
   TBOX_ASSERT(overlap);
   TBOX_ASSERT(dst_mapped_box.getLocalId() >= 0);
   TBOX_ASSERT(src_mapped_box.getLocalId() >= 0);
   TBOX_ASSERT(refine_item_id >= 0);
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst_level,
      *src_level,
      dst_mapped_box,
      src_mapped_box);

   // Note: s_num_coarsen_items cannot be used at this point!

   if (d_dst_patch_rank == dst_level->getBoxLevel()->getMPI().getRank()) {
      d_dst_patch = dst_level->getPatch(dst_mapped_box.getGlobalId());
   }
   if (d_src_patch_rank == dst_level->getBoxLevel()->getMPI().getRank()) {
      d_src_patch = src_level->getPatch(src_mapped_box.getGlobalId());
   }
}

RefineCopyTransaction::~RefineCopyTransaction()
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
RefineCopyTransaction::canEstimateIncomingMessageSize()
{
   bool can_estimate = false;
   if (d_src_patch) {
      can_estimate =
         d_src_patch->getPatchData(s_refine_items[d_refine_item_id]->d_src)
         ->canEstimateStreamSizeFromBox();
   } else {
      can_estimate =
         d_dst_patch->getPatchData(s_refine_items[d_refine_item_id]->d_scratch)
         ->canEstimateStreamSizeFromBox();
   }
   return can_estimate;
}

size_t
RefineCopyTransaction::computeIncomingMessageSize()
{
   d_incoming_bytes =
      d_dst_patch->getPatchData(s_refine_items[d_refine_item_id]->d_scratch)
      ->getDataStreamSize(*d_overlap);
   return d_incoming_bytes;
}

size_t
RefineCopyTransaction::computeOutgoingMessageSize()
{
   d_outgoing_bytes =
      d_src_patch->getPatchData(s_refine_items[d_refine_item_id]->d_src)
      ->getDataStreamSize(*d_overlap);
   return d_outgoing_bytes;
}

int
RefineCopyTransaction::getSourceProcessor() {
   return d_src_patch_rank;
}

int
RefineCopyTransaction::getDestinationProcessor() {
   return d_dst_patch_rank;
}

void
RefineCopyTransaction::packStream(
   tbox::MessageStream& stream)
{
   d_src_patch->getPatchData(s_refine_items[d_refine_item_id]->d_src)
   ->packStream(stream, *d_overlap);
}

void
RefineCopyTransaction::unpackStream(
   tbox::MessageStream& stream)
{
   d_dst_patch->getPatchData(s_refine_items[d_refine_item_id]->d_scratch)
   ->unpackStream(stream, *d_overlap);
}

void
RefineCopyTransaction::copyLocalData()
{
   hier::PatchData& dst_data =
      *d_dst_patch->getPatchData(s_refine_items[d_refine_item_id]->d_scratch);

   const hier::PatchData& src_data =
      *d_src_patch->getPatchData(s_refine_items[d_refine_item_id]->d_src);

   dst_data.copy(src_data, *d_overlap);
}

/*
 *************************************************************************
 *
 * Function to print state of transaction.
 *
 *************************************************************************
 */

void
RefineCopyTransaction::printClassData(
   std::ostream& stream) const
{
   stream << "Refine Copy Transaction" << std::endl;
   stream << "   refine item array:        "
          << (RefineClasses::Data **)s_refine_items << std::endl;
   stream << "   num refine items:       " << s_num_refine_items << std::endl;
   stream << "   destination patch rank:       " << d_dst_patch_rank
          << std::endl;
   stream << "   source patch rank:            " << d_src_patch_rank
          << std::endl;
   stream << "   refine item id:         " << d_refine_item_id << std::endl;
   stream << "   destination patch data id: "
          << s_refine_items[d_refine_item_id]->d_scratch << std::endl;
   stream << "   source patch data id:      "
          << s_refine_items[d_refine_item_id]->d_src << std::endl;
   stream << "   incoming bytes:         " << d_incoming_bytes << std::endl;
   stream << "   outgoing bytes:         " << d_outgoing_bytes << std::endl;
   stream << "   destination patch:           "
          << d_dst_patch.get() << std::endl;
   stream << "   source patch:           "
          << d_src_patch.get() << std::endl;
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
