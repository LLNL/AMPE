/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for coarsening AMR data.
 *
 ************************************************************************/

#ifndef included_xfer_CoarsenPatchStrategy_C
#define included_xfer_CoarsenPatchStrategy_C

#include "SAMRAI/xfer/CoarsenPatchStrategy.h"

namespace SAMRAI {
namespace xfer {

CoarsenPatchStrategy::CoarsenPatchStrategy(
   const tbox::Dimension& dim):
   d_dim(dim)
{
   registerObject();
}

CoarsenPatchStrategy::~CoarsenPatchStrategy()
{
}

/*
 *************************************************************************
 * Compute the max coarsen stencil width from all constructed
 * coarsen patch strategies.
 *************************************************************************
 */
hier::IntVector
CoarsenPatchStrategy::getMaxCoarsenOpStencilWidth(
   const tbox::Dimension& dim)
{
   hier::IntVector max_width(dim, 0);

   std::set<CoarsenPatchStrategy *>& current_objects =
      CoarsenPatchStrategy::getCurrentObjects();
   for (std::set<CoarsenPatchStrategy *>::const_iterator
        si = current_objects.begin(); si != current_objects.end(); ++si) {
      const CoarsenPatchStrategy* op = *si;
      if (op->getDim() == dim) {
         max_width.max(op->getCoarsenOpStencilWidth());
      }
   }

   return max_width;
}

}
}
#endif
