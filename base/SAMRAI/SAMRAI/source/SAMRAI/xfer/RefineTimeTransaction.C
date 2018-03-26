/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Communication transaction for time interpolation during data refining
 *
 ************************************************************************/

#ifndef included_xfer_RefineTimeTransaction_C
#define included_xfer_RefineTimeTransaction_C

#include "SAMRAI/xfer/RefineTimeTransaction.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <typeinfo>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * Initialization, set/unset functions for static array of refine items
 * and interpolation time.
 *
 *************************************************************************
 */

double RefineTimeTransaction::s_time = 0.0;

const RefineClasses::Data **
RefineTimeTransaction::s_refine_items = (const RefineClasses::Data **)NULL;
int RefineTimeTransaction::s_num_refine_items = 0;

/*
 *************************************************************************
 *
 * Constructor sets state of transaction.
 *
 *************************************************************************
 */

RefineTimeTransaction::RefineTimeTransaction(
   const boost::shared_ptr<hier::PatchLevel>& dst_level,
   const boost::shared_ptr<hier::PatchLevel>& src_level,
   const boost::shared_ptr<hier::BoxOverlap>& overlap,
   const hier::Box& dst_mapped_box,
   const hier::Box& src_mapped_box,
   const hier::Box& box,
   const int refine_item_id):
   d_dst_patch(),
   d_dst_patch_rank(dst_mapped_box.getOwnerRank()),
   d_src_patch(),
   d_src_patch_rank(src_mapped_box.getOwnerRank()),
   d_overlap(overlap),
   d_box(box),
   d_refine_item_id(refine_item_id)
{
   TBOX_ASSERT(dst_level);
   TBOX_ASSERT(src_level);
   TBOX_ASSERT(overlap);
   TBOX_ASSERT(dst_mapped_box.getLocalId() >= 0);
   TBOX_ASSERT(src_mapped_box.getLocalId() >= 0);
   TBOX_ASSERT(refine_item_id >= 0);
   TBOX_DIM_ASSERT_CHECK_ARGS5(*dst_level,
      *src_level,
      dst_mapped_box,
      src_mapped_box,
      box);

   // Note: s_num_coarsen_items cannot be used at this point!

   if (d_dst_patch_rank == dst_level->getBoxLevel()->getMPI().getRank()) {
      d_dst_patch = dst_level->getPatch(dst_mapped_box.getId());
   }
   if (d_src_patch_rank == dst_level->getBoxLevel()->getMPI().getRank()) {
      d_src_patch = src_level->getPatch(src_mapped_box.getId());
   }
}

RefineTimeTransaction::~RefineTimeTransaction()
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
RefineTimeTransaction::canEstimateIncomingMessageSize()
{
   bool can_estimate = false;
   if (d_src_patch) {
      can_estimate =
         d_src_patch->getPatchData(s_refine_items[d_refine_item_id]->
            d_src_told)
         ->canEstimateStreamSizeFromBox();
   } else {
      can_estimate =
         d_dst_patch->getPatchData(s_refine_items[d_refine_item_id]->
            d_scratch)
         ->canEstimateStreamSizeFromBox();
   }
   return can_estimate;
}

size_t
RefineTimeTransaction::computeIncomingMessageSize()
{
   d_incoming_bytes =
      d_dst_patch->getPatchData(s_refine_items[d_refine_item_id]->
         d_scratch)
      ->getDataStreamSize(*d_overlap);
   return d_incoming_bytes;
}

size_t
RefineTimeTransaction::computeOutgoingMessageSize()
{
   d_outgoing_bytes =
      d_src_patch->getPatchData(s_refine_items[d_refine_item_id]->
         d_src_told)
      ->getDataStreamSize(*d_overlap);
   return d_outgoing_bytes;
}

int
RefineTimeTransaction::getSourceProcessor() {
   return d_src_patch_rank;
}

int
RefineTimeTransaction::getDestinationProcessor() {
   return d_dst_patch_rank;
}

void
RefineTimeTransaction::packStream(
   tbox::MessageStream& stream)
{
   hier::Box temporary_mapped_box(d_box.getDim());
   temporary_mapped_box.initialize(d_box, hier::LocalId(-1), tbox::SAMRAI_MPI::getInvalidRank());

   hier::Patch temporary_patch(
      temporary_mapped_box,
      d_src_patch->getPatchDescriptor());

   boost::shared_ptr<hier::PatchData> temporary_patch_data(
      d_src_patch->getPatchDescriptor()
      ->getPatchDataFactory(s_refine_items[d_refine_item_id]->
         d_src_told)
      ->allocate(temporary_patch));
   temporary_patch_data->setTime(s_time);

   timeInterpolate(
      temporary_patch_data,
      d_src_patch->getPatchData(s_refine_items[d_refine_item_id]->d_src_told),
      d_src_patch->getPatchData(s_refine_items[d_refine_item_id]->d_src_tnew));

   temporary_patch_data->packStream(stream, *d_overlap);
}

void
RefineTimeTransaction::unpackStream(
   tbox::MessageStream& stream)
{
   d_dst_patch->getPatchData(s_refine_items[d_refine_item_id]->d_scratch)
   ->unpackStream(stream, *d_overlap);
}

void
RefineTimeTransaction::copyLocalData()
{
   /*
    * If there is no offset between the source and destination, then
    * time interpolate directly to the destination patchdata.  Otherwise,
    * time interpolate into a temporary patchdata and copy the result
    * to the destination patchdata.
    */
   if (d_overlap->getSourceOffset() ==
       hier::IntVector::getZero(d_box.getDim())) {

      timeInterpolate(
         d_dst_patch->getPatchData(s_refine_items[d_refine_item_id]->
            d_scratch),
         d_src_patch->getPatchData(s_refine_items[d_refine_item_id]->
            d_src_told),
         d_src_patch->getPatchData(s_refine_items[d_refine_item_id]->
            d_src_tnew));

   } else {

      hier::Box temporary_mapped_box(d_box.getDim());
      temporary_mapped_box.initialize(d_box, hier::LocalId(-1), tbox::SAMRAI_MPI::getInvalidRank());

      hier::Patch temporary_patch(
         temporary_mapped_box,
         d_src_patch->getPatchDescriptor());

      boost::shared_ptr<hier::PatchData> temp(
         d_src_patch->getPatchDescriptor()
         ->getPatchDataFactory(s_refine_items[d_refine_item_id]->
            d_src_told)
         ->allocate(temporary_patch));

      temp->setTime(s_time);

      timeInterpolate(
         temp,
         d_src_patch->getPatchData(s_refine_items[d_refine_item_id]->
            d_src_told),
         d_src_patch->getPatchData(s_refine_items[d_refine_item_id]->
            d_src_tnew));

      d_dst_patch->getPatchData(s_refine_items[d_refine_item_id]->d_scratch)
      ->copy(*temp, *d_overlap);

   }

}

void
RefineTimeTransaction::timeInterpolate(
   const boost::shared_ptr<hier::PatchData>& pd_dst,
   const boost::shared_ptr<hier::PatchData>& pd_old,
   const boost::shared_ptr<hier::PatchData>& pd_new)
{
   TBOX_ASSERT(pd_old);
   TBOX_ASSERT(pd_dst);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*pd_dst, *pd_old);
   TBOX_ASSERT(tbox::MathUtilities<double>::equalEps(pd_dst->getTime(), s_time));

   if (tbox::MathUtilities<double>::equalEps(pd_old->getTime(), s_time)) {
      s_refine_items[d_refine_item_id]->
      d_optime->timeInterpolate(*pd_dst, d_box, *pd_old, *pd_old);
   } else {

      TBOX_ASSERT(pd_new);
      TBOX_DIM_ASSERT_CHECK_ARGS2(*pd_dst, *pd_new);
      TBOX_ASSERT(pd_old->getTime() < s_time);
      TBOX_ASSERT(pd_new->getTime() >= s_time);

      s_refine_items[d_refine_item_id]->
      d_optime->timeInterpolate(*pd_dst, d_box, *pd_old, *pd_new);
   }
}

/*
 *************************************************************************
 *
 * Function to print state of transaction.
 *
 *************************************************************************
 */

void
RefineTimeTransaction::printClassData(
   std::ostream& stream) const
{
   stream << "Refine Time Transaction" << std::endl;
   stream << "   transaction time:        " << s_time << std::endl;
   stream << "   refine item array:        "
          << (RefineClasses::Data **)s_refine_items << std::endl;
   stream << "   num refine items:        " << s_num_refine_items << std::endl;
   stream << "   destination patch rank:        " << d_dst_patch_rank
          << std::endl;
   stream << "   source patch rank:             " << d_src_patch_rank
          << std::endl;
   stream << "   time interpolation box:  " << d_box << std::endl;
   stream << "   refine item id:          " << d_refine_item_id << std::endl;
   stream << "   destination patch data id:  "
          << s_refine_items[d_refine_item_id]->d_scratch << std::endl;
   stream << "   source (old) patch data id: "
          << s_refine_items[d_refine_item_id]->d_src_told << std::endl;
   stream << "   source (new) patch data id: "
          << s_refine_items[d_refine_item_id]->d_src_tnew << std::endl;
   stream << "   time interpolation name id: "
          << typeid(*s_refine_items[d_refine_item_id]->d_optime).name() << std::endl;
   stream << "   incoming bytes:          " << d_incoming_bytes << std::endl;
   stream << "   outgoing bytes:          " << d_outgoing_bytes << std::endl;
   stream << "   destination patch:           "
          << d_dst_patch.get() << std::endl;
   stream << "   source level:           "
          << d_src_patch.get() << std::endl;
   stream << "   overlap:                 " << std::endl;
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
