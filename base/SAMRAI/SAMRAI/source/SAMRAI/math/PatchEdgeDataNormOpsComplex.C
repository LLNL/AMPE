/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Norm operations for complex edge-centered patch data.
 *
 ************************************************************************/

#ifndef included_math_PatchEdgeDataNormOpsComplex_C
#define included_math_PatchEdgeDataNormOpsComplex_C

#include "SAMRAI/math/PatchEdgeDataNormOpsComplex.h"
#include "SAMRAI/pdat/EdgeGeometry.h"
#include "SAMRAI/tbox/MathUtilities.h"

namespace SAMRAI {
namespace math {

PatchEdgeDataNormOpsComplex::PatchEdgeDataNormOpsComplex()
{
}

PatchEdgeDataNormOpsComplex::~PatchEdgeDataNormOpsComplex()
{
}
/*
 *************************************************************************
 *
 * Compute the number of data entries on a patch in the given box.
 *
 *************************************************************************
 */

int
PatchEdgeDataNormOpsComplex::numberOfEntries(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const hier::Box& box) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   int dimVal = box.getDim().getValue();
   int retval = 0;
   const hier::Box ibox = box * data->getGhostBox();
   const int data_depth = data->getDepth();
   for (int d = 0; d < dimVal; d++) {
      retval += (pdat::EdgeGeometry::toEdgeBox(ibox, d).size() * data_depth);
   }
   return retval;
}

/*
 *************************************************************************
 *
 * Norm operations for complex edge-centered data.
 *
 *************************************************************************
 */

double
PatchEdgeDataNormOpsComplex::sumControlVolumes(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const boost::shared_ptr<pdat::EdgeData<double> >& cvol,
   const hier::Box& box) const
{
   TBOX_ASSERT(data && cvol);

   int dimVal = box.getDim().getValue();
   double retval = 0.0;
   for (int d = 0; d < dimVal; d++) {
      retval += d_array_ops.sumControlVolumes(data->getArrayData(d),
         cvol->getArrayData(d),
         pdat::EdgeGeometry::toEdgeBox(box, d));
   }
   return retval;
}

void
PatchEdgeDataNormOpsComplex::abs(
   const boost::shared_ptr<pdat::EdgeData<double> >& dst,
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& src,
   const hier::Box& box) const
{
   TBOX_ASSERT(dst && src);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*dst, *src, box);

   int dimVal = box.getDim().getValue();
   for (int d = 0; d < dimVal; d++) {
      d_array_ops.abs(dst->getArrayData(d),
         src->getArrayData(d),
         pdat::EdgeGeometry::toEdgeBox(box, d));
   }
}

double
PatchEdgeDataNormOpsComplex::L1Norm(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::EdgeData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   int dimVal = box.getDim().getValue();

   double retval = 0.0;
   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         retval += d_array_ops.L1Norm(data->getArrayData(d), edge_box);
      }
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         retval += d_array_ops.L1NormWithControlVolume(data->getArrayData(d),
               cvol->getArrayData(d),
               edge_box);
      }
   }
   return retval;
}

double
PatchEdgeDataNormOpsComplex::L2Norm(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::EdgeData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   TBOX_DIM_ASSERT_CHECK_ARGS2(*data, box);

   int dimVal = box.getDim().getValue();

   double retval = 0.0;
   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         double aval = d_array_ops.L2Norm(data->getArrayData(d), edge_box);
         retval += aval * aval;
      }
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         double aval = d_array_ops.L2NormWithControlVolume(
               data->getArrayData(d),
               cvol->getArrayData(d),
               edge_box);
         retval += aval * aval;
      }
   }
   return sqrt(retval);
}

double
PatchEdgeDataNormOpsComplex::weightedL2Norm(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& weight,
   const hier::Box& box,
   const boost::shared_ptr<pdat::EdgeData<double> >& cvol) const
{
   TBOX_ASSERT(data && weight);
   TBOX_DIM_ASSERT_CHECK_ARGS3(*data, *weight, box);

   int dimVal = box.getDim().getValue();

   double retval = 0.0;
   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         double aval = d_array_ops.weightedL2Norm(data->getArrayData(d),
               weight->getArrayData(d),
               edge_box);
         retval += aval * aval;
      }
   } else {
      TBOX_DIM_ASSERT_CHECK_ARGS2(*data, *cvol);

      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         double aval = d_array_ops.weightedL2NormWithControlVolume(
               data->getArrayData(d),
               weight->getArrayData(d),
               cvol->getArrayData(d),
               edge_box);
         retval += aval * aval;
      }
   }
   return sqrt(retval);
}

double
PatchEdgeDataNormOpsComplex::RMSNorm(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::EdgeData<double> >& cvol) const
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

double
PatchEdgeDataNormOpsComplex::weightedRMSNorm(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& weight,
   const hier::Box& box,
   const boost::shared_ptr<pdat::EdgeData<double> >& cvol) const
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

double
PatchEdgeDataNormOpsComplex::maxNorm(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::EdgeData<double> >& cvol) const
{
   TBOX_ASSERT(data);
   int dimVal = box.getDim().getValue();

   double retval = 0.0;
   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box =
            pdat::EdgeGeometry::toEdgeBox(box, d);
         retval = tbox::MathUtilities<double>::Max(retval,
               d_array_ops.maxNorm(data->getArrayData(d), edge_box));
      }
   } else {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box =
            pdat::EdgeGeometry::toEdgeBox(box, d);
         retval = tbox::MathUtilities<double>::Max(retval,
               d_array_ops.maxNormWithControlVolume(
                  data->getArrayData(d), cvol->getArrayData(d), edge_box));
      }
   }
   return retval;
}

dcomplex
PatchEdgeDataNormOpsComplex::dot(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data1,
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data2,
   const hier::Box& box,
   const boost::shared_ptr<pdat::EdgeData<double> >& cvol) const
{
   TBOX_ASSERT(data1 && data2);
   int dimVal = box.getDim().getValue();

   dcomplex retval = dcomplex(0.0, 0.0);
   if (!cvol) {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         retval += d_array_ops.dot(data1->getArrayData(d),
               data2->getArrayData(d),
               edge_box);
      }
   } else {
      for (int d = 0; d < dimVal; d++) {
         const hier::Box edge_box = pdat::EdgeGeometry::toEdgeBox(box, d);
         retval += d_array_ops.dotWithControlVolume(
               data1->getArrayData(d),
               data2->getArrayData(d),
               cvol->getArrayData(d),
               edge_box);
      }
   }
   return retval;
}

dcomplex
PatchEdgeDataNormOpsComplex::integral(
   const boost::shared_ptr<pdat::EdgeData<dcomplex> >& data,
   const hier::Box& box,
   const boost::shared_ptr<pdat::EdgeData<double> >& vol) const
{
   TBOX_ASSERT(data);

   int dimVal = box.getDim().getValue();
   dcomplex retval = dcomplex(0.0, 0.0);
   for (int d = 0; d < dimVal; d++) {
      retval += d_array_ops.integral(
            data->getArrayData(d),
            vol->getArrayData(d),
            pdat::EdgeGeometry::toEdgeBox(box, d));
   }
   return retval;
}

}
}
#endif
