/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Norm operations for integer data arrays.
 *
 ************************************************************************/

#ifndef included_math_ArrayDataNormOpsInteger_C
#define included_math_ArrayDataNormOpsInteger_C

#include "SAMRAI/math/ArrayDataNormOpsInteger.h"

#include "SAMRAI/tbox/MathUtilities.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace math {

ArrayDataNormOpsInteger::ArrayDataNormOpsInteger()
{
}

ArrayDataNormOpsInteger::~ArrayDataNormOpsInteger()
{
}

/*
 *************************************************************************
 *
 * Norm operations for integer array data.
 *
 *************************************************************************
 */

void
ArrayDataNormOpsInteger::abs(
   pdat::ArrayData<int>& dst,
   const pdat::ArrayData<int>& src,
   const hier::Box& box) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(dst, src, box);
   TBOX_ASSERT(dst.getDepth() == src.getDepth());

   int dimVal = dst.getDim().getValue();

   const hier::Box dst_box = dst.getBox();
   const hier::Box src_box = src.getBox();
   const hier::Box ibox = box * dst_box * src_box;

   if (!ibox.empty()) {

      int box_w[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
      int dst_w[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
      int src_w[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
      int dim_counter[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
      for (int i = 0; i < dimVal; i++) {
         box_w[i] = ibox.numberCells(i);
         dst_w[i] = dst_box.numberCells(i);
         src_w[i] = src_box.numberCells(i);
         dim_counter[i] = 0;
      }

      const int dst_offset = dst.getOffset();
      const int src_offset = src.getOffset();

      const int num_d0_blocks = ibox.size() / box_w[0];

      int dst_begin = dst_box.offset(ibox.lower());
      int src_begin = src_box.offset(ibox.lower());

      int* dd = dst.getPointer();
      const int* sd = src.getPointer();

      const int ddepth = dst.getDepth();
      for (int d = 0; d < ddepth; d++) {

         int dst_counter = dst_begin;
         int src_counter = src_begin;

         int dst_b[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
         int src_b[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
         for (int nd = 0; nd < dimVal; nd++) {
            dst_b[nd] = dst_counter;
            src_b[nd] = src_counter;
         }

         for (int nb = 0; nb < num_d0_blocks; nb++) {

            for (int i0 = 0; i0 < box_w[0]; i0++) {
               dd[dst_counter + i0] =
                  tbox::MathUtilities<int>::Abs(sd[src_counter + i0]);
            }

            int dim_jump = 0;

            for (int j = 1; j < dimVal; j++) {
               if (dim_counter[j] < box_w[j] - 1) {
                  ++dim_counter[j];
                  dim_jump = j;
                  break;
               } else {
                  dim_counter[j] = 0;
               }
            }

            if (dim_jump > 0) {
               int dst_step = 1;
               int src_step = 1;
               for (int k = 0; k < dim_jump; k++) {
                  dst_step *= dst_w[k];
                  src_step *= src_w[k];
               }
               dst_counter = dst_b[dim_jump - 1] + dst_step;
               src_counter = src_b[dim_jump - 1] + src_step;

               for (int m = 0; m < dim_jump; m++) {
                  dst_b[m] = dst_counter;
                  src_b[m] = src_counter;
               }
            }
         }

         dst_begin += dst_offset;
         src_begin += src_offset;

      }
   }
}

}
}
#endif
