/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Base class for patch data test operations.
 *
 ************************************************************************/

#include "PatchMultiblockTestStrategy.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Utilities.h"

using namespace SAMRAI;

// These are used in the cell tagging routine.
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

/*
 *************************************************************************
 *
 * The constructor and destructor.
 *
 *************************************************************************
 */

PatchMultiblockTestStrategy::PatchMultiblockTestStrategy(
   const tbox::Dimension& dim):
   d_dim(dim)
{
   d_variable_src_name.resizeArray(0);
   d_variable_dst_name.resizeArray(0);
   d_variable_depth.resizeArray(0);
   d_variable_src_ghosts.resizeArray(0, hier::IntVector(d_dim));
   d_variable_dst_ghosts.resizeArray(0, hier::IntVector(d_dim));
   d_variable_coarsen_op.resizeArray(0);
   d_variable_refine_op.resizeArray(0);
}

PatchMultiblockTestStrategy::~PatchMultiblockTestStrategy()
{
}

/*
 *************************************************************************
 *
 * Routines for reading variable and refinement data from input.
 *
 *************************************************************************
 */

void PatchMultiblockTestStrategy::readVariableInput(
   boost::shared_ptr<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(db);
#endif

   tbox::Array<string> var_keys = db->getAllKeys();
   int nkeys = var_keys.getSize();

   d_variable_src_name.resizeArray(nkeys);
   d_variable_dst_name.resizeArray(nkeys);
   d_variable_depth.resizeArray(nkeys);
   d_variable_src_ghosts.resizeArray(nkeys, hier::IntVector(d_dim, 0));
   d_variable_dst_ghosts.resizeArray(nkeys, hier::IntVector(d_dim, 0));
   d_variable_coarsen_op.resizeArray(nkeys);
   d_variable_refine_op.resizeArray(nkeys);

   for (int i = 0; i < nkeys; i++) {

      boost::shared_ptr<tbox::Database> var_db(db->getDatabase(var_keys[i]));

      if (var_db->keyExists("src_name")) {
         d_variable_src_name[i] = var_db->getString("src_name");
      } else {
         TBOX_ERROR("Variable input error: No `src_name' string found for "
            << "key = " << var_keys[i] << endl);
      }

      if (var_db->keyExists("dst_name")) {
         d_variable_dst_name[i] = var_db->getString("dst_name");
      } else {
         TBOX_ERROR("Variable input error: No `dst_name' string found for "
            << "key = " << var_keys[i] << endl);
      }

      if (var_db->keyExists("depth")) {
         d_variable_depth[i] = var_db->getInteger("depth");
      } else {
         d_variable_depth[i] = 1;
      }

      if (var_db->keyExists("src_ghosts")) {
         int* tmp_ghosts = &d_variable_src_ghosts[i][0];
         var_db->getIntegerArray("src_ghosts", tmp_ghosts, d_dim.getValue());
      }

      if (var_db->keyExists("dst_ghosts")) {
         int* tmp_ghosts = &d_variable_dst_ghosts[i][0];
         var_db->getIntegerArray("dst_ghosts", tmp_ghosts, d_dim.getValue());
      }

      if (var_db->keyExists("coarsen_operator")) {
         d_variable_coarsen_op[i] = var_db->getString("coarsen_operator");
      } else {
         d_variable_coarsen_op[i] = "NO_COARSEN";
      }

      if (var_db->keyExists("refine_operator")) {
         d_variable_refine_op[i] = var_db->getString("refine_operator");
      } else {
         d_variable_refine_op[i] = "NO_REFINE";
      }

   }

}

void PatchMultiblockTestStrategy::readRefinementInput(
   boost::shared_ptr<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(db);
#endif

   tbox::Array<string> box_keys = db->getAllKeys();
   int nkeys = box_keys.getSize();

   d_refine_level_boxes.resizeArray(nkeys);
   for (int i = 0; i < nkeys; i++) {
      d_refine_level_boxes[i] = db->getDatabaseBoxArray(box_keys[i]);
   }

}

/*
 *************************************************************************
 *
 * Tag cells on level specified in input box array for refinement.
 *
 *************************************************************************
 */

void PatchMultiblockTestStrategy::tagCellsInInputBoxes(
   hier::Patch& patch,
   int level_number,
   int tag_index)
{

   if (level_number < d_refine_level_boxes.getSize()) {

      boost::shared_ptr<pdat::CellData<int> > tags(
         patch.getPatchData(tag_index),
         boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(tags);
#endif
      tags->fillAll(0);

      const hier::Box pbox = patch.getBox();

      for (hier::BoxContainer::iterator k(d_refine_level_boxes[level_number]);
           k != d_refine_level_boxes[level_number].end(); ++k) {
         tags->fill(1, *k * pbox, 0);
      }

   }

}

/*
 *************************************************************************
 *
 * Blank physical boundary and pre/postprocess coarsen/refine operations
 * so tester isn't required to implement them when not needed.
 *
 *************************************************************************
 */

void PatchMultiblockTestStrategy::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double time,
   const hier::IntVector& gcw) const
{
   NULL_USE(patch);
   NULL_USE(time);
   NULL_USE(gcw);
}

void PatchMultiblockTestStrategy::preprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const boost::shared_ptr<hier::VariableContext>& context,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(context);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}

void PatchMultiblockTestStrategy::postprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const boost::shared_ptr<hier::VariableContext>& context,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(context);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}

void PatchMultiblockTestStrategy::preprocessCoarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const boost::shared_ptr<hier::VariableContext>& context,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   NULL_USE(coarse);
   NULL_USE(fine);
   NULL_USE(context);
   NULL_USE(coarse_box);
   NULL_USE(ratio);
}

void PatchMultiblockTestStrategy::postprocessCoarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const boost::shared_ptr<hier::VariableContext>& context,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   NULL_USE(coarse);
   NULL_USE(fine);
   NULL_USE(context);
   NULL_USE(coarse_box);
   NULL_USE(ratio);
}
