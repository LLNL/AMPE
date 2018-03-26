/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Describes boundaries for a patch
 *
 ************************************************************************/

#ifndef included_hier_PatchBoundaries_C
#define included_hier_PatchBoundaries_C

#include "SAMRAI/hier/PatchBoundaries.h"

namespace SAMRAI {
namespace hier {

/*
 *************************************************************************
 *
 * Default constructor
 *
 *************************************************************************
 */
PatchBoundaries::PatchBoundaries():
   d_dim(tbox::Dimension::getInvalidDimension()),
   d_array_of_bboxes(0)
{
}

/*
 *************************************************************************
 *
 * Constructor leaves the arrays empty.
 *
 *************************************************************************
 */
PatchBoundaries::PatchBoundaries(
   const tbox::Dimension& dim):
   d_dim(dim),
   d_array_of_bboxes(dim.getValue())
{
}

/*
 *************************************************************************
 *
 * Copy constructor
 *
 *************************************************************************
 */
PatchBoundaries::PatchBoundaries(
   const PatchBoundaries& r):
   d_dim(r.d_dim),
   d_array_of_bboxes(r.d_dim.getValue())
{
   for (unsigned int d = 0; d < d_dim.getValue(); ++d) {
      d_array_of_bboxes[d] = r.d_array_of_bboxes[d];
   }
}

} // SAMRAI namespace
} // hier namespace

#endif
