/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Communication transaction for data copies during data coarsening
 *
 ************************************************************************/

#ifndef included_xfer_CoarsenCopyTransaction_C
#define included_xfer_CoarsenCopyTransaction_C

#include "SAMRAI/xfer/CoarsenCopyTransaction.h"

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

const CoarsenClasses::Data ** CoarsenCopyTransaction::s_coarsen_items =
   (const CoarsenClasses::Data **)NULL;
int CoarsenCopyTransaction::s_num_coarsen_items = 0;

/*
 *************************************************************************
 *
 * Constructor sets state of transaction.
 *
 *************************************************************************
 */

CoarsenCopyTransaction::CoarsenCopyTransaction(
   const boost::shared_ptr<hier::PatchLevel>& dst_level,
   const boost::shared_ptr<hier::PatchLevel>& src_level,
   const boost::shared_ptr<hier::BoxOverlap>& overlap,
   const hier::Box& dst_mapped_box,
   const hier::Box& src_mapped_box,
   const int coarsen_item_id):
   d_dst_patch_rank(dst_mapped_box.getOwnerRank()),
   d_src_patch_rank(src_mapped_box.getOwnerRank()),
   d_overlap(overlap),
   d_coarsen_item_id(coarsen_item_id),
   d_incoming_bytes(0),
   d_outgoing_bytes(0)
{
   TBOX_ASSERT(dst_level);
   TBOX_ASSERT(src_level);
   TBOX_ASSERT(overlap);
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst_level,
      *src_level,
      dst_mapped_box,
      src_mapped_box);
   TBOX_ASSERT(dst_mapped_box.getLocalId() >= 0);
   TBOX_ASSERT(src_mapped_box.getLocalId() >= 0);
   TBOX_ASSERT(coarsen_item_id >= 0);

   // Note: s_num_coarsen_items cannot be used at this point!

   if (d_dst_patch_rank == dst_level->getBoxLevel()->getMPI().getRank()) {
      d_dst_patch = dst_level->getPatch(dst_mapped_box.getGlobalId());
   }
   if (d_src_patch_rank == src_level->getBoxLevel()->getMPI().getRank()) {
      d_src_patch = src_level->getPatch(src_mapped_box.getGlobalId());
   }
}

CoarsenCopyTransaction::~CoarsenCopyTransaction()
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
CoarsenCopyTransaction::canEstimateIncomingMessageSize()
{
   bool can_estimate = false;
   if (d_src_patch) {
      can_estimate =
         d_src_patch->getPatchData(s_coarsen_items[d_coarsen_item_id]->d_src)
         ->canEstimateStreamSizeFromBox();
   } else {
      can_estimate =
         d_dst_patch->getPatchData(s_coarsen_items[d_coarsen_item_id]->d_dst)
         ->canEstimateStreamSizeFromBox();
   }
   return can_estimate;
}

size_t
CoarsenCopyTransaction::computeIncomingMessageSize()
{
   d_incoming_bytes =
      d_dst_patch->getPatchData(s_coarsen_items[d_coarsen_item_id]->d_dst)
      ->getDataStreamSize(*d_overlap);
   return d_incoming_bytes;
}

size_t
CoarsenCopyTransaction::computeOutgoingMessageSize()
{
   d_outgoing_bytes =
      d_src_patch->getPatchData(s_coarsen_items[d_coarsen_item_id]->d_src)
      ->getDataStreamSize(*d_overlap);
   return d_outgoing_bytes;
}

int
CoarsenCopyTransaction::getSourceProcessor()
{
   return d_src_patch_rank;
}

int
CoarsenCopyTransaction::getDestinationProcessor()
{
   return d_dst_patch_rank;
}

void
CoarsenCopyTransaction::packStream(
   tbox::MessageStream& stream)
{
   d_src_patch->getPatchData(s_coarsen_items[d_coarsen_item_id]->d_src)
   ->packStream(stream, *d_overlap);
}

void
CoarsenCopyTransaction::unpackStream(
   tbox::MessageStream& stream)
{
   d_dst_patch->getPatchData(s_coarsen_items[d_coarsen_item_id]->d_dst)
   ->unpackStream(stream, *d_overlap);
}

void
CoarsenCopyTransaction::copyLocalData()
{
   hier::PatchData& dst_data =
      *d_dst_patch->getPatchData(s_coarsen_items[d_coarsen_item_id]->d_dst);

   const hier::PatchData& src_data =
      *d_src_patch->getPatchData(s_coarsen_items[d_coarsen_item_id]->d_src);

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
CoarsenCopyTransaction::printClassData(
   std::ostream& stream) const
{
   stream << "Coarsen Copy Transaction" << std::endl;
   stream << "   coarsen item array:        "
          << (CoarsenClasses::Data **)s_coarsen_items << std::endl;
   stream << "   num coarsen items:      " << s_num_coarsen_items << std::endl;
   stream << "   destination patch rank:       " << d_dst_patch_rank
          << std::endl;
   stream << "   source patch_rank:            " << d_src_patch_rank
          << std::endl;
   stream << "   coarsen item id:        " << d_coarsen_item_id << std::endl;
   stream << "   destination patch data id: "
          << s_coarsen_items[d_coarsen_item_id]->d_dst << std::endl;
   stream << "   source patch data id:      "
          << s_coarsen_items[d_coarsen_item_id]->d_src << std::endl;
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
