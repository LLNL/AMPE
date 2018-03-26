/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Operations for complex cell-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchCellDataOpsComplex_C
#define included_math_PatchCellDataOpsComplex_C

#include "SAMRAI/math/PatchCellDataOpsComplex.h"

namespace SAMRAI {
namespace math {

PatchCellDataOpsComplex::PatchCellDataOpsComplex()
{
}

PatchCellDataOpsComplex::~PatchCellDataOpsComplex()
{
}

/*
 *************************************************************************
 *
 * General operations for complex cell-centered patch data.
 *
 *************************************************************************
 */

void
PatchCellDataOpsComplex::swapData(
   const boost::shared_ptr<hier::Patch>& patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(patch);

   boost::shared_ptr<pdat::CellData<dcomplex> > d1(
      patch->getPatchData(data1_id),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<pdat::CellData<dcomplex> > d2(
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
