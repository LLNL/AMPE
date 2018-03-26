/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for refining AMR data.
 *
 ************************************************************************/

#ifndef included_xfer_RefinePatchStrategy_C
#define included_xfer_RefinePatchStrategy_C

#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/hier/BoxContainer.h"

namespace SAMRAI {
namespace xfer {

/*
 *************************************************************************
 *
 * The default constructor and virtual destructor do nothing
 * particularly interesting.
 *
 *************************************************************************
 */

RefinePatchStrategy::RefinePatchStrategy(
   const tbox::Dimension& dim):
   d_dim(dim)
{
   registerObject();
}

RefinePatchStrategy::~RefinePatchStrategy()
{
   unregisterObject();
}

void
RefinePatchStrategy::fillSingularityBoundaryConditions(
   hier::Patch& patch,
   const hier::PatchLevel& encon_level,
   const hier::Connector& dst_to_encon,
   const double fill_time,
   const hier::Box& fill_box,
   const hier::BoundaryBox& boundary_box,
   const boost::shared_ptr<hier::BaseGridGeometry>& grid_geometry)
{
   NULL_USE(patch);
   NULL_USE(encon_level);
   NULL_USE(dst_to_encon);
   NULL_USE(fill_time);
   NULL_USE(fill_box);
   NULL_USE(boundary_box);
   NULL_USE(grid_geometry);
   TBOX_ERROR(
      "The abstract RefinePatchLevelStragey::fillSingularityBoudaryConditions:\n"
      << "must be implemented whenever the concrete derived\n"
      << "class supports multiblock and singularities.");
}

/*
 *************************************************************************
 * Compute the max refine stencil width from all constructed
 * refine patch strategies.
 *************************************************************************
 */
hier::IntVector
RefinePatchStrategy::getMaxRefineOpStencilWidth(
   const tbox::Dimension& dim)
{
   hier::IntVector max_width(dim, 0);

   std::set<RefinePatchStrategy *>& current_objects =
      RefinePatchStrategy::getCurrentObjects();
   for (std::set<RefinePatchStrategy *>::const_iterator
        si = current_objects.begin(); si != current_objects.end(); ++si) {
      const RefinePatchStrategy* op = *si;
      if (op->getDim() == dim) {
         max_width.max(op->getRefineOpStencilWidth());
      }
   }

   return max_width;
}

}
}
#endif
