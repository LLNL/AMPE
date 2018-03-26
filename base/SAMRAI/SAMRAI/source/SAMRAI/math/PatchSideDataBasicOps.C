/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Basic templated side-centered patch data operations.
 *
 ************************************************************************/

#ifndef included_math_PatchSideDataBasicOps_C
#define included_math_PatchSideDataBasicOps_C

#include "SAMRAI/math/PatchSideDataBasicOps.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/SideGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchSideDataBasicOps<TYPE>::PatchSideDataBasicOps()
{
}

template<class TYPE>
PatchSideDataBasicOps<TYPE>::~PatchSideDataBasicOps()
{
}

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
PatchSideDataBasicOps<TYPE>::PatchSideDataBasicOps(
   const PatchSideDataBasicOps<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::operator = (
   const PatchSideDataBasicOps<TYPE>& foo)
{
   NULL_USE(foo);
}

/*
 *************************************************************************
 *
 * General basic templated opertions for side data.
 *
 *************************************************************************
 */

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::scale(
   const boost::shared_ptr<pdat::SideData<TYPE> >& dst,
   const TYPE& alpha,
   const boost::shared_ptr<pdat::SideData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.scale(dst->getArrayData(d),
            alpha, src->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::addScalar(
   const boost::shared_ptr<pdat::SideData<TYPE> >& dst,
   const boost::shared_ptr<pdat::SideData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.addScalar(dst->getArrayData(d),
            src->getArrayData(d), alpha,
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::add(
   const boost::shared_ptr<pdat::SideData<TYPE> >& dst,
   const boost::shared_ptr<pdat::SideData<TYPE> >& src1,
   const boost::shared_ptr<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.add(dst->getArrayData(d),
            src1->getArrayData(d), src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::subtract(
   const boost::shared_ptr<pdat::SideData<TYPE> >& dst,
   const boost::shared_ptr<pdat::SideData<TYPE> >& src1,
   const boost::shared_ptr<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.subtract(dst->getArrayData(d),
            src1->getArrayData(d), src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::multiply(
   const boost::shared_ptr<pdat::SideData<TYPE> >& dst,
   const boost::shared_ptr<pdat::SideData<TYPE> >& src1,
   const boost::shared_ptr<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.multiply(dst->getArrayData(d),
            src1->getArrayData(d), src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::divide(
   const boost::shared_ptr<pdat::SideData<TYPE> >& dst,
   const boost::shared_ptr<pdat::SideData<TYPE> >& src1,
   const boost::shared_ptr<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.divide(dst->getArrayData(d),
            src1->getArrayData(d), src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::reciprocal(
   const boost::shared_ptr<pdat::SideData<TYPE> >& dst,
   const boost::shared_ptr<pdat::SideData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_ASSERT(dst->getDirectionVector() == src->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.reciprocal(dst->getArrayData(d),
            src->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::linearSum(
   const boost::shared_ptr<pdat::SideData<TYPE> >& dst,
   const TYPE& alpha,
   const boost::shared_ptr<pdat::SideData<TYPE> >& src1,
   const TYPE& beta,
   const boost::shared_ptr<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.linearSum(dst->getArrayData(d),
            alpha, src1->getArrayData(d),
            beta, src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::axpy(
   const boost::shared_ptr<pdat::SideData<TYPE> >& dst,
   const TYPE& alpha,
   const boost::shared_ptr<pdat::SideData<TYPE> >& src1,
   const boost::shared_ptr<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.axpy(dst->getArrayData(d),
            alpha, src1->getArrayData(d),
            src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::axmy(
   const boost::shared_ptr<pdat::SideData<TYPE> >& dst,
   const TYPE& alpha,
   const boost::shared_ptr<pdat::SideData<TYPE> >& src1,
   const boost::shared_ptr<pdat::SideData<TYPE> >& src2,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src1 && src2);
   TBOX_ASSERT(dst->getDirectionVector() == src1->getDirectionVector());
   TBOX_ASSERT(dst->getDirectionVector() == src2->getDirectionVector());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*dst, *src1, *src2, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.axmy(dst->getArrayData(d),
            alpha, src1->getArrayData(d),
            src2->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
void
PatchSideDataBasicOps<TYPE>::setRandomValues(
   const boost::shared_ptr<pdat::SideData<TYPE> >& dst,
   const TYPE& width,
   const TYPE& low,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*dst, box);

   int dimVal = dst->getDim().getValue();

   const hier::IntVector& directions = dst->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         d_array_ops.setRandomValues(dst->getArrayData(d),
            width, low, side_box);
      }
   }
}

template<class TYPE>
TYPE
PatchSideDataBasicOps<TYPE>::min(
   const boost::shared_ptr<pdat::SideData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   int dimVal = data->getDim().getValue();

   TYPE minval = tbox::MathUtilities<TYPE>::getMax();
   const hier::IntVector& directions = data->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         minval = tbox::MathUtilities<TYPE>::Min(
               minval, d_array_ops.min(data->getArrayData(d), side_box));
      }
   }
   return minval;
}

template<class TYPE>
TYPE
PatchSideDataBasicOps<TYPE>::max(
   const boost::shared_ptr<pdat::SideData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   int dimVal = data->getDim().getValue();

   TYPE maxval = -tbox::MathUtilities<TYPE>::getMax();
   const hier::IntVector& directions = data->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         maxval = tbox::MathUtilities<TYPE>::Max(
               maxval, d_array_ops.max(data->getArrayData(d), side_box));
      }
   }
   return maxval;
}

}
}
#endif
