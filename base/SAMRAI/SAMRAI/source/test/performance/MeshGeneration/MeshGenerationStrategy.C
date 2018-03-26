/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   MeshGeneration class implementation
 *
 ************************************************************************/
#include "MeshGenerationStrategy.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxLevelConnectorUtils.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"

#include <iomanip>

using namespace SAMRAI;



void MeshGenerationStrategy::setTagsByShrinkingLevel(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   int tag_ln,
   int tag_data_id,
   const hier::IntVector &shrink_cells,
   const double *shrink_distance )
{

   const tbox::Dimension dim(hierarchy->getDim());

   boost::shared_ptr<geom::CartesianGridGeometry> grid_geometry(
      hierarchy->getGridGeometry(),
      boost::detail::dynamic_cast_tag());

   const int tag_val = 1;

   const boost::shared_ptr<hier::PatchLevel> &tag_level(
      hierarchy->getPatchLevel(tag_ln));

   const hier::BoxLevel &Ltag = *tag_level->getBoxLevel();


   /*
    * Compute shrinkage in terms of coarse cell count.  It should be
    * the largest of properly converted values for shrink_cells,
    * shrink_distance and nesting width.
    */
   hier::IntVector shrink_width(dim, hierarchy->getProperNestingBuffer(tag_ln));
   shrink_width.max(shrink_cells);

   const double *ref_dx = grid_geometry->getDx();
   for ( int i=0; i<dim.getValue(); ++i ) {
      double h = ref_dx[i]/Ltag.getRefinementRatio()(i);
      shrink_width(i) = tbox::MathUtilities<int>::Max(
         static_cast<int>(0.5 + shrink_distance[i]/h),
         shrink_width(i) );
   }


   hier::BoxLevel tagfootprint(dim);
   hier::Connector Ltag_to_tagfootprint;
   const hier::Connector &Ltag_to_Ltag =
      Ltag.getPersistentOverlapConnectors().findOrCreateConnector(
         Ltag,
         shrink_width );

   hier::BoxLevelConnectorUtils blcu;
   blcu.computeInternalParts( tagfootprint,
                              Ltag_to_tagfootprint,
                              Ltag_to_Ltag,
                              -shrink_width,
                              grid_geometry->getDomainSearchTree() );
   tbox::plog << "Ltag_to_tagfootprint:\n" << Ltag_to_tagfootprint.format("Ltag->tagfootprint: ", 2);


   for (hier::PatchLevel::iterator pi(tag_level->begin());
        pi != tag_level->end(); ++pi) {

      boost::shared_ptr<hier::Patch> patch = *pi;
      boost::shared_ptr<pdat::CellData<int> > tag_data(
         patch->getPatchData(tag_data_id),
         boost::detail::dynamic_cast_tag());

      tag_data->getArrayData().fillAll(0);

      if ( !Ltag_to_tagfootprint.hasNeighborSet(patch->getBox().getId()) ) {
         tag_data->getArrayData().fillAll(1);
      }
      else {
         hier::Connector::ConstNeighborhoodIterator ni =
            Ltag_to_tagfootprint.find(patch->getBox().getId());

         for ( hier::Connector::ConstNeighborIterator na = Ltag_to_tagfootprint.begin(ni);
               na != Ltag_to_tagfootprint.end(ni); ++na ) {

            const hier::Box &tag_box = *na;
            tag_data->getArrayData().fillAll(tag_val, tag_box);

         }
      }

   }

   return;
}
