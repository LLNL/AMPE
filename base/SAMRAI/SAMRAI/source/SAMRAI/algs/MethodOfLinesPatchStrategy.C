/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Interface to application-specific patch functions to support
 *                MethodOfLines integration algorithm
 *
 ************************************************************************/

#ifndef included_algs_MethodOfLinesPatchStrategy_C
#define included_algs_MethodOfLinesPatchStrategy_C

#include "SAMRAI/algs/MethodOfLinesPatchStrategy.h"
#include "SAMRAI/hier/VariableDatabase.h"

namespace SAMRAI {
namespace algs {

/*
 *************************************************************************
 *
 * Note: hier::Variable contexts should be consistent with those in
 *       MethodOfLinesIntegrator class.
 *
 *************************************************************************
 */

MethodOfLinesPatchStrategy::MethodOfLinesPatchStrategy(
   const tbox::Dimension& dim):
   xfer::RefinePatchStrategy(dim),
   xfer::CoarsenPatchStrategy(dim),
   d_dim(dim)
{
   d_interior_with_ghosts.reset();
   d_interior.reset();
}

MethodOfLinesPatchStrategy::~MethodOfLinesPatchStrategy()
{
}

hier::IntVector
MethodOfLinesPatchStrategy::getRefineOpStencilWidth() const
{
   return hier::IntVector::getZero(d_dim);
}

hier::IntVector
MethodOfLinesPatchStrategy::getCoarsenOpStencilWidth() const
{
   return hier::IntVector::getZero(d_dim);
}

}
}
#endif
