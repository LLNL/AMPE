/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Operations for integer cell-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchCellDataOpsInteger_C
#define included_math_PatchCellDataOpsInteger_C

#include "SAMRAI/math/PatchCellDataOpsInteger.h"

namespace SAMRAI {
namespace math {

PatchCellDataOpsInteger::PatchCellDataOpsInteger()
{
}

PatchCellDataOpsInteger::~PatchCellDataOpsInteger()
{
}

/*
 *************************************************************************
 *
 * General operations for integer cell-centered patch data.
 *
 *************************************************************************
 */

void
PatchCellDataOpsInteger::swapData(
   const boost::shared_ptr<hier::Patch>& patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(patch);

   boost::shared_ptr<pdat::CellData<int> > d1(
      patch->getPatchData(data1_id),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<pdat::CellData<int> > d2(
      patch->getPatchData(data2_id),
      boost::detail::dynamic_cast_tag());

   TBOX_ASSERT(d1 && d2);
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox().isSpatiallyEqual(d2->getBox()));
   TBOX_ASSERT(d1->getGhostBox().isSpatiallyEqual(d2->getGhostBox()));

   patch->setPatchData(data1_id, d2);
   patch->setPatchData(data2_id, d1);
}

}
}
#endif
