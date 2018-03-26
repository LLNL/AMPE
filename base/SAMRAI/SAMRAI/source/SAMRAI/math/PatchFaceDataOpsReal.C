/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Templated operations for real face-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchFaceDataOpsReal_C
#define included_math_PatchFaceDataOpsReal_C

#include "SAMRAI/math/PatchFaceDataOpsReal.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/FaceGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchFaceDataOpsReal<TYPE>::PatchFaceDataOpsReal()
{
}

#if 0
/*
 * This was moved into the header due to what looks like bug in the
 * XLC compiler.
 */
template<class TYPE>
PatchFaceDataOpsReal<TYPE>::~PatchFaceDataOpsReal()
{
}
#endif

/*
 *************************************************************************
 *
 * The const constructor and assignment operator are not actually used
 * but are defined here for compilers that require an implementation for
 * every declaration.
 *
 *************************************************************************
 */

template<class TYPE>
PatchFaceDataOpsReal<TYPE>::PatchFaceDataOpsReal(
   const PatchFaceDataOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
void
PatchFaceDataOpsReal<TYPE>::operator = (
   const PatchFaceDataOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

/*
 *************************************************************************
 *
 * General templated operations for real face-centered patch data.
 *
 *************************************************************************
 */

template<class TYPE>
void
PatchFaceDataOpsReal<TYPE>::swapData(
   const boost::shared_ptr<hier::Patch>& patch,
   const int data1_id,
   const int data2_id) const
{
   TBOX_ASSERT(patch);

   boost::shared_ptr<pdat::FaceData<TYPE> > d1(
      patch->getPatchData(data1_id),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<pdat::FaceData<TYPE> > d2(
      patch->getPatchData(data2_id),
      boost::detail::dynamic_cast_tag());

   TBOX_ASSERT(d1 && d2);
   TBOX_ASSERT(d1->getDepth() && d2->getDepth());
   TBOX_ASSERT(d1->getBox().isSpatiallyEqual(d2->getBox()));
   TBOX_ASSERT(d1->getGhostBox().isSpatiallyEqual(d2->getGhostBox()));

   patch->setPatchData(data1_id, d2);
   patch->setPatchData(data2_id, d1);
}

template<class TYPE>
void
PatchFaceDataOpsReal<TYPE>::printData(
   const boost::shared_ptr<pdat::FaceData<TYPE> >& data,
   const hier::Box& box,
   std::ostream& s) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   s << "Data box = " << box << std::endl;
   data->print(box, s);
   s << "\n";
}

template<class TYPE>
void
PatchFaceDataOpsReal<TYPE>::copyData(
   const boost::shared_ptr<pdat::FaceData<TYPE> >& dst,
   const boost::shared_ptr<pdat::FaceData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   int dimVal = dst->getDim().getValue();

   for (int d = 0; d < dimVal; d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      (dst->getArrayData(d)).copy(src->getArrayData(d), face_box);
   }
}

template<class TYPE>
void
PatchFaceDataOpsReal<TYPE>::setToScalar(
   const boost::shared_ptr<pdat::FaceData<TYPE> >& dst,
   const TYPE& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   dst->fillAll(alpha, box);
}

}
}
#endif
