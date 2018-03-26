/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Templated norm operations for real cell-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchCellDataNormOpsReal_C
#define included_math_PatchCellDataNormOpsReal_C

#include "SAMRAI/math/PatchCellDataNormOpsReal.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
PatchCellDataNormOpsReal<TYPE>::PatchCellDataNormOpsReal()
{
}

template<class TYPE>
PatchCellDataNormOpsReal<TYPE>::~PatchCellDataNormOpsReal()
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
PatchCellDataNormOpsReal<TYPE>::PatchCellDataNormOpsReal(
   const PatchCellDataNormOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
void
PatchCellDataNormOpsReal<TYPE>::operator = (
   const PatchCellDataNormOpsReal<TYPE>& foo)
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
PatchCellDataNormOpsReal<TYPE>::numberOfEntries(
   const boost::shared_ptr<pdat::CellData<TYPE> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   const hier::Box ibox = box * data->getGhostBox();
   int retval = (ibox.size() * data->getDepth());
   return retval;
}

/*
 *************************************************************************
 *
 * Templated norm operations for real cell-centered data.
 *
 *************************************************************************
 */

template<class TYPE>
double
PatchCellDataNormOpsReal<TYPE>::sumControlVolumes(
   const boost::shared_ptr<pdat::CellData<TYPE> >& data,
   const boost::shared_ptr<pdat::CellData<double> >& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT(data && cvol);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   return d_array_ops.sumControlVolumes(data->getArrayData(),
      cvol->getArrayData(),
      box);
}

template<class TYPE>
void
PatchCellDataNormOpsReal<TYPE>::abs(
   const boost::shared_ptr<pdat::CellData<TYPE> >& dst,
   const boost::shared_ptr<pdat::CellData<TYPE> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   d_array_ops.abs(dst->getArrayData(),
      src->getArrayData(),
      box);
}

template<class TYPE>
double
PatchCellDataNormOpsReal<TYPE>::L1Norm(
   const boost::shared_ptr<pdat::CellData<TYPE> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::CellData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   double retval;
   if (!cvol) {
      retval = d_array_ops.L1Norm(data->getArrayData(), box);
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      retval = d_array_ops.L1NormWithControlVolume(data->getArrayData(),
            cvol->getArrayData(),
            box);
   }
   return retval;
}

template<class TYPE>
double
PatchCellDataNormOpsReal<TYPE>::L2Norm(
   const boost::shared_ptr<pdat::CellData<TYPE> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::CellData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   double retval;
   if (!cvol) {
      retval = d_array_ops.L2Norm(data->getArrayData(), box);
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      retval = d_array_ops.L2NormWithControlVolume(data->getArrayData(),
            cvol->getArrayData(),
            box);
   }
   return retval;
}

template<class TYPE>
double
PatchCellDataNormOpsReal<TYPE>::weightedL2Norm(
   const boost::shared_ptr<pdat::CellData<TYPE> >& data,
   const boost::shared_ptr<pdat::CellData<TYPE> >& weight,
   const hier::Box& box,
   const boost::shared_ptr<pdat::CellData<double> >& cvol) const
{
   TBOX_ASSERT(data && weight);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*data, *weight, box);

   double retval;
   if (!cvol) {
      retval = d_array_ops.weightedL2Norm(data->getArrayData(),
            weight->getArrayData(),
            box);
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      retval = d_array_ops.weightedL2NormWithControlVolume(
            data->getArrayData(),
            weight->getArrayData(),
            cvol->getArrayData(),
            box);
   }
   return retval;
}

template<class TYPE>
double
PatchCellDataNormOpsReal<TYPE>::RMSNorm(
   const boost::shared_ptr<pdat::CellData<TYPE> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::CellData<double> >& cvol) const
{
   TBOX_ASSERT(data);

   double retval = L2Norm(data, box, cvol);
   if (!cvol) {
      retval /= sqrt((double)numberOfEntries(data, box));
   } else {
      retval /= sqrt(sumControlVolumes(data, cvol, box));
   }
   return retval;
}

template<class TYPE>
double
PatchCellDataNormOpsReal<TYPE>::weightedRMSNorm(
   const boost::shared_ptr<pdat::CellData<TYPE> >& data,
   const boost::shared_ptr<pdat::CellData<TYPE> >& weight,
   const hier::Box& box,
   const boost::shared_ptr<pdat::CellData<double> >& cvol) const
{
   TBOX_ASSERT(data && weight);

   double retval = weightedL2Norm(data, weight, box, cvol);
   if (!cvol) {
      retval /= sqrt((double)numberOfEntries(data, box));
   } else {
      retval /= sqrt(sumControlVolumes(data, cvol, box));
   }
   return retval;
}

template<class TYPE>
double
PatchCellDataNormOpsReal<TYPE>::maxNorm(
   const boost::shared_ptr<pdat::CellData<TYPE> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::CellData<double> >& cvol) const
{
   TBOX_ASSERT(data);

   double retval;
   if (!cvol) {
      retval = d_array_ops.maxNorm(data->getArrayData(), box);
   } else {
      retval = d_array_ops.maxNormWithControlVolume(data->getArrayData(),
            cvol->getArrayData(),
            box);
   }
   return retval;
}

template<class TYPE>
TYPE
PatchCellDataNormOpsReal<TYPE>::dot(
   const boost::shared_ptr<pdat::CellData<TYPE> >& data1,
   const boost::shared_ptr<pdat::CellData<TYPE> >& data2,
   const hier::Box& box,
   const boost::shared_ptr<pdat::CellData<double> >& cvol) const
{
   TBOX_ASSERT(data1 && data2);

   TYPE retval;
   if (!cvol) {
      retval = d_array_ops.dot(data1->getArrayData(),
            data2->getArrayData(),
            box);
   } else {
      retval = d_array_ops.dotWithControlVolume(
            data1->getArrayData(),
            data2->getArrayData(),
            cvol->getArrayData(),
            box);
   }
   return retval;
}

template<class TYPE>
TYPE
PatchCellDataNormOpsReal<TYPE>::integral(
   const boost::shared_ptr<pdat::CellData<TYPE> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::CellData<double> >& vol) const
{
   TBOX_ASSERT(data);

   TYPE retval = 0.0;

   retval = d_array_ops.integral(data->getArrayData(),
         vol->getArrayData(),
         box);

   return retval;
}

}
}
#endif
