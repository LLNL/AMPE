/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface to user routines for refining AMR data.
 *
 ************************************************************************/

#ifndef included_mesh_MultiblockGriddingTagger_C
#define included_mesh_MultiblockGriddingTagger_C

#include "SAMRAI/mesh/MultiblockGriddingTagger.h"

#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace mesh {

/*
 *************************************************************************
 *
 * The default constructor and virtual destructor do nothing
 * particularly interesting.
 *
 *************************************************************************
 */

MultiblockGriddingTagger::MultiblockGriddingTagger(
   const tbox::Dimension& dim):
   xfer::RefinePatchStrategy(dim),
   d_dim(dim)
{
}

MultiblockGriddingTagger::~MultiblockGriddingTagger()
{
}

hier::IntVector
MultiblockGriddingTagger::getRefineOpStencilWidth() const
{
   return hier::IntVector::getOne(d_dim);
}

void
MultiblockGriddingTagger::setScratchTagPatchDataIndex(
   int buf_tag_indx)
{

   boost::shared_ptr<hier::Variable> check_var;
   bool indx_maps_to_variable =
      hier::VariableDatabase::getDatabase()->mapIndexToVariable(buf_tag_indx,
         check_var);
   if (!indx_maps_to_variable || !check_var) {
      TBOX_ERROR(
         "MultiblockGriddingTagger::setScratchTagPatchDataIndex error...\n"
         << "Given patch data index = " << buf_tag_indx
         << " is not in VariableDatabase."
         << std::endl);
   } else {
      boost::shared_ptr<pdat::CellVariable<int> > t_check_var(
         check_var,
         boost::detail::dynamic_cast_tag());
      if (!t_check_var) {
         TBOX_ERROR(
            "MultiblockGriddingTagger::setScratchTagPatchDataIndex error...\n"
            << "Given patch data index = " << buf_tag_indx
            << " does not map to cell-centered"
            << "\ninteger data in VariableDatabase." << std::endl);
      }
   }

   d_buf_tag_indx = buf_tag_indx;
}

void
MultiblockGriddingTagger::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double fill_time,
   const hier::IntVector& ghost_width_to_fill)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(d_dim, patch);

   NULL_USE(fill_time);

   const boost::shared_ptr<pdat::CellData<int> > tag_data(
      patch.getPatchData(d_buf_tag_indx),
      boost::detail::dynamic_cast_tag());

   hier::IntVector gcw =
      hier::IntVector::min(ghost_width_to_fill,
         tag_data->getGhostCellWidth());

   boost::shared_ptr<hier::PatchGeometry> pgeom(patch.getPatchGeometry());

   for (int d = 0; d < d_dim.getValue(); d++) {

      tbox::Array<hier::BoundaryBox> bbox =
         pgeom->getCodimensionBoundaries(d + 1);

      for (int b = 0; b < bbox.size(); b++) {
         if (!bbox[b].getIsMultiblockSingularity()) {
            hier::Box fill_box = pgeom->getBoundaryFillBox(bbox[b],
                  patch.getBox(),
                  gcw);

            tag_data->fillAll(0, fill_box);
         }
      }
   }
}

void
MultiblockGriddingTagger::fillSingularityBoundaryConditions(
   hier::Patch& patch,
   const hier::PatchLevel& encon_level,
   const hier::Connector& dst_to_encon,
   const double fill_time,
   const hier::Box& fill_box,
   const hier::BoundaryBox& boundary_box,
   const boost::shared_ptr<hier::BaseGridGeometry>& grid_geometry)
{
   NULL_USE(fill_time);
   NULL_USE(boundary_box);
   NULL_USE(grid_geometry);

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS3(d_dim, patch, fill_box, boundary_box);

   const tbox::Dimension& dim = fill_box.getDim();

   const hier::BoxId& dst_mb_id = patch.getBox().getId();

   const hier::BlockId& patch_blk_id = patch.getBox().getBlockId();

   const boost::shared_ptr<pdat::CellData<int> > tag_data(
      patch.getPatchData(d_buf_tag_indx),
      boost::detail::dynamic_cast_tag());

   hier::Box sing_fill_box(tag_data->getGhostBox() * fill_box);
   tag_data->fillAll(0, sing_fill_box);

   if (grid_geometry->hasEnhancedConnectivity()) {

      const std::list<hier::BaseGridGeometry::Neighbor>& neighbors =
         grid_geometry->getNeighbors(patch_blk_id);

      hier::Connector::ConstNeighborhoodIterator ni =
         dst_to_encon.findLocal(dst_mb_id);

      if (ni != dst_to_encon.end()) {

         for (hier::Connector::ConstNeighborIterator ei = dst_to_encon.begin(ni);
              ei != dst_to_encon.end(ni); ++ei) {

            boost::shared_ptr<hier::Patch> encon_patch(
               encon_level.getPatch(ei->getId()));

            const hier::BlockId& encon_blk_id = ei->getBlockId();

            hier::Transformation::RotationIdentifier rotation =
               hier::Transformation::NO_ROTATE;
            hier::IntVector offset(dim);

            for (std::list<hier::BaseGridGeometry::Neighbor>::const_iterator
                 nbri = neighbors.begin(); nbri != neighbors.end(); nbri++) {

               if (nbri->getBlockId() == encon_blk_id) {
                  rotation = nbri->getRotationIdentifier();
                  offset = nbri->getShift();
                  break;
               }
            }

            offset *= patch.getPatchGeometry()->getRatio();

            hier::Transformation transformation(
               rotation, offset, encon_blk_id, patch_blk_id);
            hier::Box encon_patch_box(encon_patch->getBox());
            transformation.transform(encon_patch_box);

            hier::Box encon_fill_box(encon_patch_box * sing_fill_box);
            if (!encon_fill_box.empty()) {

               const hier::Transformation::RotationIdentifier back_rotate =
                  hier::Transformation::getReverseRotationIdentifier(
                     rotation, dim);

               hier::IntVector back_shift(dim);

               hier::Transformation::calculateReverseShift(
                  back_shift, offset, rotation);

               hier::Transformation back_trans(back_rotate, back_shift,
                                               encon_fill_box.getBlockId(),
                                               encon_patch->getBox().getBlockId()); 
                                               

               boost::shared_ptr<pdat::CellData<int> > sing_data(
                  encon_patch->getPatchData(d_buf_tag_indx),
                  boost::detail::dynamic_cast_tag());

               pdat::CellIterator ciend(encon_fill_box, false);
               for (pdat::CellIterator ci(encon_fill_box, true);
                    ci != ciend; ++ci) {
                  pdat::CellIndex src_index(*ci);
                  pdat::CellGeometry::transform(src_index, back_trans);

                  int sing_val = (*sing_data)(src_index);
                  if (sing_val != 0 && (*tag_data)(*ci) == 0) {
                     (*tag_data)(*ci) = sing_val;
                  }
               }
            }
         }
      }
   }
}

void
MultiblockGriddingTagger::preprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio)
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}

void
MultiblockGriddingTagger::postprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(d_dim, fine, coarse, fine_box, ratio);

   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
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
