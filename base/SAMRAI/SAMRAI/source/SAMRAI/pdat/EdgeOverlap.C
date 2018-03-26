/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_EdgeOverlap_C
#define included_pdat_EdgeOverlap_C

#include "SAMRAI/pdat/EdgeOverlap.h"

namespace SAMRAI {
namespace pdat {

EdgeOverlap::EdgeOverlap(
   const tbox::Array<hier::BoxContainer>& boxes,
   const hier::Transformation& transformation):
   d_is_overlap_empty(true),
   d_transformation(transformation)
{
   const tbox::Dimension dim(transformation.getOffset().getDim());
   d_dst_boxes.resizeArray(boxes.getSize());

   for (int d = 0; d < boxes.getSize(); d++) {
      d_dst_boxes[d] = boxes[d];
      if (!d_dst_boxes[d].isEmpty()) d_is_overlap_empty = false;
   }
}

EdgeOverlap::~EdgeOverlap()
{
}

bool
EdgeOverlap::isOverlapEmpty() const
{
   return d_is_overlap_empty;
}

const hier::BoxContainer&
EdgeOverlap::getDestinationBoxContainer(
   const int axis) const
{
   TBOX_ASSERT((axis >= 0) && (axis < d_dst_boxes.getSize()));

   return d_dst_boxes[axis];
}

const hier::IntVector&
EdgeOverlap::getSourceOffset() const
{
   return d_transformation.getOffset();
}

const hier::Transformation&
EdgeOverlap::getTransformation() const
{
   return d_transformation;
}

}
}
#endif
