/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Templated norm operations for real side-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchSideDataNormOpsReal_C
#define included_math_PatchSideDataNormOpsReal_C

#include "SAMRAI/math/PatchSideDataNormOpsReal.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/pdat/SideGeometry.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchSideDataNormOpsReal<TYPE>::PatchSideDataNormOpsReal()
{
}

template<class TYPE>
PatchSideDataNormOpsReal<TYPE>::~PatchSideDataNormOpsReal()
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
PatchSideDataNormOpsReal<TYPE>::PatchSideDataNormOpsReal(
   const PatchSideDataNormOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
void
PatchSideDataNormOpsReal<TYPE>::operator = (
   const PatchSideDataNormOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

/*
 *************************************************************************
 *
 * Compute the number of data entries on a patch in the given box.
 *
 *************************************************************************
 */

template<class TYPE>
int
PatchSideDataNormOpsReal<TYPE>::numberOfEntries(
   const boost::shared_ptr<pdat::SideData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   int dimVal = box.getDim().getValue();

   int retval = 0;
   const hier::Box ibox = box * data->getGhostBox();
   const hier::IntVector& directions = data->getDirectionVector();
   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         const hier::Box dbox = pdat::SideGeometry::toSideBox(ibox, d);
         retval += (dbox.size() * data->getDepth());
      }
   }
   return retval;
}

/*
 *************************************************************************
 *
 * Templated norm operations for real side-centered data.
 *
 *************************************************************************
 */

template<class TYPE>
double
PatchSideDataNormOpsReal<TYPE>::sumControlVolumes(
   const boost::shared_ptr<pdat::SideData<TYPE> >& data,
   const boost::shared_ptr<pdat::SideData<double> >& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT(data && cvol);

   double retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();

   TBOX_ASSERT(directions ==
      hier::IntVector::min(directions, cvol->getDirectionVector()));

   int dimVal = data->getDim().getValue();

   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         retval += d_array_ops.sumControlVolumes(data->getArrayData(d),
               cvol->getArrayData(d),
               side_box);
      }
   }
   return retval;
}

template<class TYPE>
void
PatchSideDataNormOpsReal<TYPE>::abs(
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
         d_array_ops.abs(dst->getArrayData(d),
            src->getArrayData(d),
            side_box);
      }
   }
}

template<class TYPE>
double
PatchSideDataNormOpsReal<TYPE>::L1Norm(
   const boost::shared_ptr<pdat::SideData<TYPE> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   int dimVal = data->getDim().getValue();

   double retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();
   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            retval += d_array_ops.L1Norm(data->getArrayData(d), side_box);
         }
      }
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));

      for (int d = 0; d < dimVal; d++) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            retval += d_array_ops.L1NormWithControlVolume(data->getArrayData(d),
                  cvol->getArrayData(d),
                  side_box);
         }
      }
   }
   return retval;
}

template<class TYPE>
double
PatchSideDataNormOpsReal<TYPE>::L2Norm(
   const boost::shared_ptr<pdat::SideData<TYPE> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   int dimVal = data->getDim().getValue();

   double retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();
   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            double aval = d_array_ops.L2Norm(data->getArrayData(d), side_box);
            retval += aval * aval;
         }
      }
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));

      for (int d = 0; d < dimVal; d++) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            double aval = d_array_ops.L2NormWithControlVolume(
                  data->getArrayData(d),
                  cvol->getArrayData(d),
                  side_box);
            retval += aval * aval;
         }
      }
   }
   return sqrt(retval);
}

template<class TYPE>
double
PatchSideDataNormOpsReal<TYPE>::weightedL2Norm(
   const boost::shared_ptr<pdat::SideData<TYPE> >& data,
   const boost::shared_ptr<pdat::SideData<TYPE> >& weight,
   const hier::Box& box,
   const boost::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data && weight);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*data, *weight, box);

   int dimVal = data->getDim().getValue();

   double retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();

   TBOX_ASSERT(directions ==
      hier::IntVector::min(directions, weight->getDirectionVector()));

   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            double aval = d_array_ops.weightedL2Norm(data->getArrayData(d),
                  weight->getArrayData(d),
                  side_box);
            retval += aval * aval;
         }
      }
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));

      for (int d = 0; d < dimVal; d++) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            double aval = d_array_ops.weightedL2NormWithControlVolume(
                  data->getArrayData(d),
                  weight->getArrayData(d),
                  cvol->getArrayData(d),
                  side_box);
            retval += aval * aval;
         }
      }
   }
   return sqrt(retval);
}

template<class TYPE>
double
PatchSideDataNormOpsReal<TYPE>::RMSNorm(
   const boost::shared_ptr<pdat::SideData<TYPE> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::SideData<double> >& cvol) const
{
// SGS

   TBOX_ASSERT(data);

   double retval = L2Norm(data, box, cvol);
   if (!cvol) {
      retval /= sqrt((double)numberOfEntries(data, box));
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);
      retval /= sqrt(sumControlVolumes(data, cvol, box));
   }
   return retval;
}

template<class TYPE>
double
PatchSideDataNormOpsReal<TYPE>::weightedRMSNorm(
   const boost::shared_ptr<pdat::SideData<TYPE> >& data,
   const boost::shared_ptr<pdat::SideData<TYPE> >& weight,
   const hier::Box& box,
   const boost::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data && weight);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   double retval = weightedL2Norm(data, weight, box, cvol);
   if (!cvol) {
      retval /= sqrt((double)numberOfEntries(data, box));
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);
      retval /= sqrt(sumControlVolumes(data, cvol, box));
   }
   return retval;
}

template<class TYPE>
double
PatchSideDataNormOpsReal<TYPE>::maxNorm(
   const boost::shared_ptr<pdat::SideData<TYPE> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   int dimVal = data->getDim().getValue();

   double retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();
   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         if (directions(d)) {
            const hier::Box side_box =
               pdat::SideGeometry::toSideBox(box, d);
            retval = tbox::MathUtilities<double>::Max(retval,
                  d_array_ops.maxNorm(data->getArrayData(d), side_box));
         }
      }
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));

      for (int d = 0; d < dimVal; d++) {
         if (directions(d)) {
            const hier::Box side_box =
               pdat::SideGeometry::toSideBox(box, d);
            retval = tbox::MathUtilities<double>::Max(retval,
                  d_array_ops.maxNormWithControlVolume(
                     data->getArrayData(d),
                     cvol->getArrayData(d),
                     side_box));
         }
      }
   }
   return retval;
}

template<class TYPE>
TYPE
PatchSideDataNormOpsReal<TYPE>::dot(
   const boost::shared_ptr<pdat::SideData<TYPE> >& data1,
   const boost::shared_ptr<pdat::SideData<TYPE> >& data2,
   const hier::Box& box,
   const boost::shared_ptr<pdat::SideData<double> >& cvol) const
{
   TBOX_ASSERT(data1 && data2);
   TBOX_ASSERT(data1->getDirectionVector() == data2->getDirectionVector());

   int dimVal = data1->getDim().getValue();

   TYPE retval = 0.0;
   const hier::IntVector& directions = data1->getDirectionVector();
   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            retval += d_array_ops.dot(data1->getArrayData(d),
                  data2->getArrayData(d),
                  side_box);
         }
      }
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data1, *cvol);
      TBOX_ASSERT(directions ==
         hier::IntVector::min(directions, cvol->getDirectionVector()));

      for (int d = 0; d < dimVal; d++) {
         if (directions(d)) {
            const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
            retval += d_array_ops.dotWithControlVolume(
                  data1->getArrayData(d),
                  data2->getArrayData(d),
                  cvol->getArrayData(d),
                  side_box);
         }
      }
   }
   return retval;
}

template<class TYPE>
TYPE
PatchSideDataNormOpsReal<TYPE>::integral(
   const boost::shared_ptr<pdat::SideData<TYPE> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::SideData<double> >& vol) const
{
   TBOX_ASSERT(data);
   TBOX_ASSERT(vol);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*data, box, *vol);

   int dimVal = data->getDim().getValue();

   TYPE retval = 0.0;
   const hier::IntVector& directions = data->getDirectionVector();

   TBOX_ASSERT(directions ==
      hier::IntVector::min(directions, vol->getDirectionVector()));

   for (int d = 0; d < dimVal; d++) {
      if (directions(d)) {
         const hier::Box side_box = pdat::SideGeometry::toSideBox(box, d);
         retval += d_array_ops.integral(
               data->getArrayData(d),
               vol->getArrayData(d),
               side_box);
      }
   }

   return retval;
}

}
}
#endif
