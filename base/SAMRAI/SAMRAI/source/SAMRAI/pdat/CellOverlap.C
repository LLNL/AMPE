/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   hier
 *
 ************************************************************************/

#ifndef included_pdat_CellOverlap_C
#define included_pdat_CellOverlap_C

#include "SAMRAI/pdat/CellOverlap.h"

#include "SAMRAI/hier/BoxContainer.h"

namespace SAMRAI {
namespace pdat {

CellOverlap::CellOverlap(
   const hier::BoxContainer& boxes,
   const hier::Transformation& transformation):
   d_is_overlap_empty(boxes.isEmpty()),
   d_transformation(transformation),
   d_dst_boxes(boxes)
{
}

CellOverlap::~CellOverlap()
{
}

bool
CellOverlap::isOverlapEmpty() const
{
   return d_is_overlap_empty;
}

const hier::BoxContainer&
CellOverlap::getDestinationBoxContainer() const
{
   return d_dst_boxes;
}

const hier::IntVector&
CellOverlap::getSourceOffset() const
{
   return d_transformation.getOffset();
}

const hier::Transformation&
CellOverlap::getTransformation() const
{
   return d_transformation;
}

void
CellOverlap::print(
   std::ostream& os) const
{
   os << "CellOverlap boxes:";
   for (hier::BoxContainer::const_iterator b(d_dst_boxes);
        b != d_dst_boxes.end(); ++b) {
      const hier::Box& box = *b;
      os << "  " << box << std::endl;
   }
}

}
}
#endif
