/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Base class for geometry management on patches
 *
 ************************************************************************/

#ifndef included_hier_PatchGeometry_C
#define included_hier_PatchGeometry_C

#include "SAMRAI/hier/PatchGeometry.h"

#include "SAMRAI/hier/BoundaryLookupTable.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

PatchGeometry::PatchGeometry(
   const IntVector& ratio_to_level_zero,
   const TwoDimBool& touches_regular_bdry,
   const TwoDimBool& touches_periodic_bdry):
   d_dim(ratio_to_level_zero.getDim()),
   d_ratio_to_level_zero(ratio_to_level_zero),
   d_patch_boundaries(ratio_to_level_zero.getDim()),
   d_touches_regular_bdry(ratio_to_level_zero.getDim())

{
#ifndef DEBUG_CHECK_ASSERTIONS
   NULL_USE(touches_periodic_bdry);
#endif

   TBOX_DIM_ASSERT_CHECK_ARGS3(ratio_to_level_zero,
      touches_regular_bdry,
      touches_periodic_bdry);

#ifdef DEBUG_CHECK_ASSERTIONS

   /*
    * All components of ratio must me nonzero.  Additionally, all components
    * of ratio not equal to 1 must have the same sign.
    */
   int i;
   for (i = 0; i < d_dim.getValue(); i++) {
      TBOX_ASSERT(ratio_to_level_zero(i) != 0);
   }
   if (d_dim.getValue() > 1) {
      for (i = 0; i < d_dim.getValue(); i++) {
         TBOX_ASSERT((ratio_to_level_zero(i)
                      * ratio_to_level_zero((i + 1) % d_dim.getValue()) > 0)
            || (ratio_to_level_zero(i) == 1)
            || (ratio_to_level_zero((i + 1) % d_dim.getValue()) == 1));
      }
   }
#endif

   d_has_regular_boundary = false;
   d_has_periodic_boundary = false;

   for (int axis = 0; axis < d_dim.getValue(); axis++) {
      for (int dir = 0; dir < 2; dir++) {
         d_touches_regular_bdry(axis, dir) = touches_regular_bdry(axis, dir);

         if (d_touches_regular_bdry(axis, dir)) {
            d_has_regular_boundary = true;
         }
      }
   }
}

PatchGeometry::~PatchGeometry()
{
}

Box
PatchGeometry::getBoundaryFillBox(
   const BoundaryBox& bbox,
   const Box& patch_box,
   const IntVector& gcw) const
{

   TBOX_DIM_ASSERT_CHECK_ARGS3(bbox, patch_box, gcw);

#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < d_dim.getValue(); i++) {
      TBOX_ASSERT(gcw(i) >= 0);
   }
#endif
   Box tmp_box(patch_box);
   tmp_box.grow(gcw);
   Box fill_box(bbox.getBox() * tmp_box);

   int bdry_type = bbox.getBoundaryType();
   int location_index = bbox.getLocationIndex();

   // Get the singleton class lookup table
   const BoundaryLookupTable* blut;
   blut = BoundaryLookupTable::getLookupTable(d_dim);

#ifdef DEBUG_CHECK_ASSERTIONS
   const tbox::Array<int>& location_index_max = blut->getMaxLocationIndices();
   TBOX_ASSERT(bdry_type > 0);
   TBOX_ASSERT(bdry_type <= d_dim.getValue());
   TBOX_ASSERT(location_index >= 0);
#endif

   if (!fill_box.empty()) {

      // Loop over codimension (a.k.a. boundary type)
      for (int codim = 1; codim <= d_dim.getValue(); codim++) {

         // When we get a match on the boundary type
         if (bdry_type == codim) {

            TBOX_ASSERT(location_index < location_index_max[codim - 1]);

            // Get the directions involved in this boundary type from the
            // lookup table.
            const tbox::Array<int>& dir =
               blut->getDirections(location_index, codim);

            // For each direction, identify this as an upper or lower boundary.
            for (int i = 0; i < codim; i++) {
               if (blut->isUpper(location_index, codim, i)) {
                  fill_box.growUpper(dir[i], gcw(dir[i]) - 1);
               } else {
                  fill_box.growLower(dir[i], gcw(dir[i]) - 1);
               }
            }

            // We've found boundary type, so break out of the loop.
            break;
         }
      }
   }

   return fill_box;
}

void
PatchGeometry::setCodimensionBoundaries(
   const tbox::Array<BoundaryBox>& bdry_boxes,
   int codim)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int i = 0; i < bdry_boxes.size(); i++) {
      TBOX_ASSERT(bdry_boxes[i].getBoundaryType() == codim);
   }
   TBOX_ASSERT(codim <= d_dim.getValue());
   TBOX_ASSERT(codim > 0);
#endif

   d_patch_boundaries[codim - 1].resizeArray(bdry_boxes.size(),
      BoundaryBox(d_dim));

   for (int b = 0; b < bdry_boxes.size(); b++) {
      d_patch_boundaries[codim - 1][b] = bdry_boxes[b];
   }
}

void
PatchGeometry::setBoundaryBoxesOnPatch(
   const tbox::Array<tbox::Array<BoundaryBox> > bdry)
{
   for (int i = 0; i < d_dim.getValue(); i++) {
      setCodimensionBoundaries(bdry[i], i + 1);
   }
}

void
PatchGeometry::printClassData(
   std::ostream& stream) const
{
   stream << "\nPatchGeometry::printClassData..." << std::endl;
   stream << "Ratio to level zero = " << d_ratio_to_level_zero << std::endl;
   stream << "d_has_regular_boundary = "
          << d_has_regular_boundary << std::endl;
   stream << "Boundary boxes for patch..." << std::endl;
   for (int d = 0; d < d_dim.getValue(); d++) {
      const int n = d_patch_boundaries[d].getSize();
      stream << "Boundary box array " << d << " has " << n << " boxes"
             << std::endl;
      for (int i = 0; i < n; i++) {
         stream << "box " << i << " = "
                << d_patch_boundaries[d][i].getBox() << std::endl;
      }
   }
}

PatchGeometry::TwoDimBool::TwoDimBool():
   d_dim(tbox::Dimension::getInvalidDimension())
{
}

PatchGeometry::TwoDimBool::TwoDimBool(
   const tbox::Dimension& dim):
   d_dim(dim)
{
   TBOX_DIM_ASSERT_CHECK_DIM(dim);
   setAll(false);
}

PatchGeometry::TwoDimBool::TwoDimBool(
   const tbox::Dimension& dim,
   bool v):
   d_dim(dim)
{
   for (int i = 0; i < 2 * d_dim.getValue(); ++i) {
      d_data[i] = v;
   }
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif

#endif
