/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Fill pattern class to provide interface for stencils
 *
 ************************************************************************/

#ifndef included_pdat_SecondLayerNodeVariableFillPattern_C
#define included_pdat_SecondLayerNodeVariableFillPattern_C

#include "SAMRAI/pdat/SecondLayerNodeVariableFillPattern.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/tbox/Utilities.h"

#include <boost/make_shared.hpp>

namespace SAMRAI {
namespace pdat {

const std::string SecondLayerNodeVariableFillPattern::s_name_id =
   "SECOND_LAYER_NODE_FILL_PATTERN";

/*
 *************************************************************************
 *
 * Constructor
 *
 *************************************************************************
 */

SecondLayerNodeVariableFillPattern::SecondLayerNodeVariableFillPattern(
   const tbox::Dimension& dim):
   d_dim(dim)
{
}

/*
 *************************************************************************
 *
 * Destructor
 *
 *************************************************************************
 */

SecondLayerNodeVariableFillPattern::~SecondLayerNodeVariableFillPattern()
{
}

/*
 *************************************************************************
 *
 * Calculate the overlap according to the desired pattern
 *
 *************************************************************************
 */

boost::shared_ptr<hier::BoxOverlap>
SecondLayerNodeVariableFillPattern::calculateOverlap(
   const hier::BoxGeometry& dst_geometry,
   const hier::BoxGeometry& src_geometry,
   const hier::Box& dst_patch_box,
   const hier::Box& src_mask,
   const hier::Box& fill_box,
   const bool overwrite_interior,
   const hier::Transformation& transformation) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(dst_patch_box, src_mask);
   NULL_USE(overwrite_interior);

   const tbox::Dimension dim(dst_patch_box.getDim());

   hier::Box dst_node_box(NodeGeometry::toNodeBox(dst_patch_box));
   hier::Box src_node_mask(NodeGeometry::toNodeBox(src_mask));
   bool corner_overlap = ((dst_node_box * src_node_mask).size() == 1)
      ? true : false;

   hier::BoxContainer stencil_boxes;
   computeStencilBoxes(stencil_boxes, dst_patch_box);
   if (corner_overlap) {
      hier::IntVector grow_vec(dim, 0);
      hier::Box grow_box(dim);
      hier::BoxContainer remove_list;
      for (unsigned int i = 0; i < dim.getValue(); i++) {
         grow_box = dst_node_box;
         grow_vec(i) = 1;
         grow_box.grow(grow_vec);
         remove_list.pushFront(grow_box);
         grow_vec(i) = 0;
      }
      stencil_boxes.removeIntersections(remove_list);
   }

   hier::BoxContainer dst_boxes;

   const NodeGeometry* t_dst =
      dynamic_cast<const NodeGeometry *>(&dst_geometry);
   const NodeGeometry* t_src =
      dynamic_cast<const NodeGeometry *>(&src_geometry);

   TBOX_ASSERT(t_dst);
   TBOX_ASSERT(t_src);

   t_dst->computeDestinationBoxes(dst_boxes, *t_src, src_mask, fill_box,
      false, transformation);

   dst_boxes.intersectBoxes(stencil_boxes);

   return boost::make_shared<NodeOverlap>(dst_boxes, transformation);

}

/*
 *************************************************************************
 *
 * Return the stencil width (1)
 *
 *************************************************************************
 */

const hier::IntVector&
SecondLayerNodeVariableFillPattern::getStencilWidth()
{
   return hier::IntVector::getOne(d_dim);
}

/*
 *************************************************************************
 *
 * Return the string name identifier
 *
 *************************************************************************
 */

const std::string&
SecondLayerNodeVariableFillPattern::getPatternName() const
{
   return s_name_id;
}

/*
 *************************************************************************
 *
 * Compute the boxes for the stencil around a given patch box
 *
 *************************************************************************
 */

void
SecondLayerNodeVariableFillPattern::computeStencilBoxes(
   hier::BoxContainer& stencil_boxes,
   const hier::Box& dst_box) const
{
   TBOX_ASSERT(stencil_boxes.size() == 0);

   hier::Box dst_node_box(NodeGeometry::toNodeBox(dst_box));

   hier::Box ghost_box(dst_node_box);
   ghost_box.grow(hier::IntVector::getOne(dst_box.getDim()));
   stencil_boxes.removeIntersections(ghost_box, dst_node_box);
}

/*
 *************************************************************************
 *
 * Compute BoxOverlap that specifies data to be filled by refinement
 * operator.
 *
 *************************************************************************
 */

boost::shared_ptr<hier::BoxOverlap>
SecondLayerNodeVariableFillPattern::computeFillBoxesOverlap(
   const hier::BoxContainer& fill_boxes,
   const hier::Box& patch_box,
   const hier::Box& data_box,
   const hier::PatchDataFactory& pdf) const
{
   NULL_USE(pdf);

   const tbox::Dimension& dim = patch_box.getDim();

   hier::BoxContainer stencil_boxes;
   computeStencilBoxes(stencil_boxes, patch_box);

   hier::BoxContainer overlap_boxes(fill_boxes);

   /*
    * This is the equivalent of converting every box in overlap_boxes
    * to a node centering, which must be done before intersecting with
    * stencil_boxes, which is node-centered.
    */
   for (hier::BoxContainer::iterator b(overlap_boxes);
        b != overlap_boxes.end(); ++b) {
      b->growUpper(hier::IntVector::getOne(patch_box.getDim()));
   }

   overlap_boxes.intersectBoxes(NodeGeometry::toNodeBox(data_box));

   overlap_boxes.intersectBoxes(stencil_boxes);

   overlap_boxes.coalesce();

   return boost::make_shared<NodeOverlap>(
      overlap_boxes,
      hier::Transformation(hier::IntVector::getZero(dim)));
}

}
}
#endif
