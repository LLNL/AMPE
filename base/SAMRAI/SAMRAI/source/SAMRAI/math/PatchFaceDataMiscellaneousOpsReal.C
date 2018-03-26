/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Templated miscellaneous operations for real face-centered data.
 *
 ************************************************************************/

#ifndef included_math_PatchFaceDataMiscellaneousOpsReal_C
#define included_math_PatchFaceDataMiscellaneousOpsReal_C

#include "SAMRAI/math/PatchFaceDataMiscellaneousOpsReal.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/FaceGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchFaceDataMiscellaneousOpsReal<TYPE>::PatchFaceDataMiscellaneousOpsReal()
{
}

template<class TYPE>
PatchFaceDataMiscellaneousOpsReal<TYPE>::~PatchFaceDataMiscellaneousOpsReal()
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
PatchFaceDataMiscellaneousOpsReal<TYPE>::PatchFaceDataMiscellaneousOpsReal(
   const PatchFaceDataMiscellaneousOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
void
PatchFaceDataMiscellaneousOpsReal<TYPE>::operator = (
   const PatchFaceDataMiscellaneousOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

/*
 *************************************************************************
 *
 * Templated miscellaneous opertions for real face-centered data.
 *
 *************************************************************************
 */

template<class TYPE>
int
PatchFaceDataMiscellaneousOpsReal<TYPE>::computeConstrProdPos(
   const boost::shared_ptr<pdat::FaceData<TYPE> >& data1,
   const boost::shared_ptr<pdat::FaceData<TYPE> >& data2,
   const hier::Box& box,
   const boost::shared_ptr<pdat::FaceData<double> >& cvol) const
{
   TBOX_ASSERT(data1 && data2);

   int dimVal = data1->getDim().getValue();

   int retval = 1;
   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box face_box =
            pdat::FaceGeometry::toFaceBox(box, d);
         retval = tbox::MathUtilities<int>::Min(retval,
               d_array_ops.computeConstrProdPos(
                  data1->getArrayData(d),
                  data2->getArrayData(d),
                  face_box));
      }
   } else {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box face_box =
            pdat::FaceGeometry::toFaceBox(box, d);
         retval = tbox::MathUtilities<int>::Min(retval,
               d_array_ops.computeConstrProdPosWithControlVolume(
                  data1->getArrayData(d),
                  data2->getArrayData(d),
                  cvol->getArrayData(d),
                  face_box));
      }
   }
   return retval;
}

template<class TYPE>
void
PatchFaceDataMiscellaneousOpsReal<TYPE>::compareToScalar(
   const boost::shared_ptr<pdat::FaceData<TYPE> >& dst,
   const boost::shared_ptr<pdat::FaceData<TYPE> >& src,
   const TYPE& alpha,
   const hier::Box& box,
   const boost::shared_ptr<pdat::FaceData<double> >& cvol) const
{
   TBOX_ASSERT(dst && src);

   int dimVal = dst->getDim().getValue();

   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         d_array_ops.compareToScalar(dst->getArrayData(d),
            src->getArrayData(d),
            alpha,
            face_box);
      }
   } else {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
         d_array_ops.compareToScalarWithControlVolume(dst->getArrayData(d),
            src->getArrayData(d),
            alpha,
            cvol->getArrayData(d),
            face_box);
      }
   }
}

template<class TYPE>
int
PatchFaceDataMiscellaneousOpsReal<TYPE>::testReciprocal(
   const boost::shared_ptr<pdat::FaceData<TYPE> >& dst,
   const boost::shared_ptr<pdat::FaceData<TYPE> >& src,
   const hier::Box& box,
   const boost::shared_ptr<pdat::FaceData<double> >& cvol) const
{
   TBOX_ASSERT(dst && src);

   int dimVal = dst->getDim().getValue();

   int retval = 1;
   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box face_box =
            pdat::FaceGeometry::toFaceBox(box, d);
         retval = tbox::MathUtilities<int>::Min(retval,
               d_array_ops.testReciprocal(
                  dst->getArrayData(d),
                  src->getArrayData(d),
                  face_box));
      }
   } else {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box face_box =
            pdat::FaceGeometry::toFaceBox(box, d);
         retval = tbox::MathUtilities<int>::Min(retval,
               d_array_ops.testReciprocalWithControlVolume(
                  dst->getArrayData(d),
                  src->getArrayData(d),
                  cvol->getArrayData(d),
                  face_box));
      }
   }
   return retval;
}

template<class TYPE>
TYPE
PatchFaceDataMiscellaneousOpsReal<TYPE>::maxPointwiseDivide(
   const boost::shared_ptr<pdat::FaceData<TYPE> >& numer,
   const boost::shared_ptr<pdat::FaceData<TYPE> >& denom,
   const hier::Box& box) const
{
   TBOX_ASSERT(numer && denom);

   int dimVal = numer->getDim().getValue();

   TYPE retval = 0.0;
   for (int d = 0; d < dimVal; d++) {
      const hier::Box face_box =
         pdat::FaceGeometry::toFaceBox(box, d);
      TYPE dirval = d_array_ops.maxPointwiseDivide(numer->getArrayData(d),
            denom->getArrayData(d),
            face_box);
      retval = tbox::MathUtilities<TYPE>::Max(retval, dirval);
   }
   return retval;
}

template<class TYPE>
TYPE
PatchFaceDataMiscellaneousOpsReal<TYPE>::minPointwiseDivide(
   const boost::shared_ptr<pdat::FaceData<TYPE> >& numer,
   const boost::shared_ptr<pdat::FaceData<TYPE> >& denom,
   const hier::Box& box) const
{
   TBOX_ASSERT(numer && denom);

   int dimVal = numer->getDim().getValue();

   TYPE retval = 0.0;
   for (int d = 0; d < dimVal; d++) {
      const hier::Box face_box = pdat::FaceGeometry::toFaceBox(box, d);
      TYPE dirval = d_array_ops.minPointwiseDivide(numer->getArrayData(d),
            denom->getArrayData(d),
            face_box);
      retval = tbox::MathUtilities<TYPE>::Min(retval, dirval);
   }
   return retval;
}

}
}
#endif
