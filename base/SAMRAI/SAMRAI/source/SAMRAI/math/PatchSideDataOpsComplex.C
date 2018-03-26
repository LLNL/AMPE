/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Operations for complex side-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchSideDataOpsComplex_C
#define included_math_PatchSideDataOpsComplex_C

#include "SAMRAI/math/PatchSideDataOpsComplex.h"
#include "SAMRAI/pdat/SideGeometry.h"

namespace SAMRAI {
namespace math {

PatchSideDataOpsComplex::PatchSideDataOpsComplex()
{
}

PatchSideDataOpsComplex::~PatchSideDataOpsComplex()
{
}

/*
 *************************************************************************
 *
 * General operations for complex side-centered patch data.
 *
 *************************************************************************
 */

void
PatchSideDataOpsComplex::swapData(
   const boost::shared_ptr<hier::Patch>& patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(patch);

   boost::shared_ptr<pdat::SideData<dcomplex> > d1(
      patch->getPatchData(data1_id),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<pdat::SideData<dcomplex> > d2(
      patch->getPatchData(data2_id),
      boost::detail::dynamic_cast_tag());

   TBOX_ASSERT(d1 && d2);
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox().isSpatiallyEqual(d2->getBox()));
   TBOX_ASSERT(d1->getDirectionVector() == d2->getDirectionVector());
   TBOX_ASSERT(d1->getGhostBox().isSpatiallyEqual(d2->getGhostBox()));

   patch->setPatchData(data1_id, d2);
   patch->setPatchData(data2_id, d1);
}

void
PatchSideDataOpsComplex::printData(
   const boost::shared_ptr<pdat::SideData<dcomplex> >& data,
   const hier::Box& box,
   std::ostream& s) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

void
PatchSideDataOpsComplex::copyData(
   const boost::shared_ptr<pdat::SideData<dcomplex> >& dst,
   const boost::shared_ptr<pdat::SideData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   int dimVal = box.getDim().getValue();
   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         dst->getArrayData(d).copy(src->getArrayData(d),
            pdat::SideGeometry::toSideBox(box, d));
      }
   }
}

}
}
#endif
