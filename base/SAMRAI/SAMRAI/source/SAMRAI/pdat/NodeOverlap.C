/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_NodeOverlap_C
#define included_pdat_NodeOverlap_C

#include "SAMRAI/pdat/NodeOverlap.h"

namespace SAMRAI {
namespace pdat {

NodeOverlap::NodeOverlap(
   const hier::BoxContainer& boxes,
   const hier::Transformation& transformation):
   d_is_overlap_empty(boxes.isEmpty()),
   d_transformation(transformation),
   d_dst_boxes(boxes)
{
}

NodeOverlap::~NodeOverlap()
{
}

bool
NodeOverlap::isOverlapEmpty() const
{
   return d_is_overlap_empty;
}

const hier::BoxContainer&
NodeOverlap::getDestinationBoxContainer() const
{
   return d_dst_boxes;
}

const hier::IntVector&
NodeOverlap::getSourceOffset() const
{
   return d_transformation.getOffset();
}

const hier::Transformation&
NodeOverlap::getTransformation() const
{
   return d_transformation;
}

}
}
#endif
