/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   A flattened representation of a hierarchy
 *
 ************************************************************************/
#include "SAMRAI/hier/FlattenedHierarchy.h"

#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/hier/HierarchyNeighbors.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"

#include "boost/make_shared.hpp"

namespace SAMRAI {
namespace hier {

/*
 ***************************************************************************
 * Constructor evaluates the hierarchy and fills the visible boxes container
 ***************************************************************************
 */

FlattenedHierarchy::FlattenedHierarchy(
   const PatchHierarchy& hierarchy,
   int coarsest_level,
   int finest_level)
: d_coarsest_level(coarsest_level),
  d_finest_level(finest_level),
  d_patch_hierarchy(&hierarchy)
{
   int num_levels = hierarchy.getNumberOfLevels();
   TBOX_ASSERT(coarsest_level >= 0);
   TBOX_ASSERT(coarsest_level <= finest_level);
   TBOX_ASSERT(finest_level < num_levels);

   d_visible_boxes.resize(num_levels);

   LocalId local_id(0);

   for (int ln = finest_level; ln >= coarsest_level; --ln) {
      const boost::shared_ptr<PatchLevel>& current_level =
         hierarchy.getPatchLevel(ln);

      if (ln != finest_level) {

         const Connector& coarse_to_fine =
            current_level->findConnector(
               *(hierarchy.getPatchLevel(ln+1)),
               IntVector::getOne(hierarchy.getDim()),
               CONNECTOR_IMPLICIT_CREATION_RULE,
               true);

         const IntVector& connector_ratio = coarse_to_fine.getRatio();

         for (PatchLevel::iterator ip(current_level->begin());
              ip != current_level->end(); ++ip) {

            const boost::shared_ptr<Patch>& patch = *ip;
            const Box& box = patch->getBox();
            const BlockId& block_id = box.getBlockId();
            const BoxId& box_id = box.getBoxId();
            BoxContainer& visible_boxes = d_visible_boxes[ln][box_id];

            BoxContainer coarse_boxes(box);

            BoxContainer fine_nbr_boxes;
            if (coarse_to_fine.hasNeighborSet(box_id)) {
               coarse_to_fine.getNeighborBoxes(box_id, fine_nbr_boxes);
            }
            if (!fine_nbr_boxes.empty()) {
               BoxContainer fine_boxes;
               for (SAMRAI::hier::RealBoxConstIterator
                    nbr_itr = fine_nbr_boxes.realBegin();
                    nbr_itr != fine_nbr_boxes.realEnd(); ++nbr_itr) {
                  if (nbr_itr->getBlockId() == block_id) {
                     fine_boxes.pushBack(*nbr_itr);
                  }
               }

               fine_boxes.coarsen(connector_ratio);

               coarse_boxes.removeIntersections(fine_boxes);
               coarse_boxes.coalesce();
            }

            for (BoxContainer::iterator itr =
                 coarse_boxes.begin(); itr != coarse_boxes.end(); ++itr) {

               Box new_box(*itr, local_id, box_id.getOwnerRank());
               ++local_id;
               visible_boxes.insert(visible_boxes.end(), new_box);
            }
         }
      } else {
         for (PatchLevel::iterator ip(current_level->begin());
              ip != current_level->end(); ++ip) {
            const boost::shared_ptr<Patch>& patch = *ip;
            const Box& box = patch->getBox();
            const BoxId& box_id = box.getBoxId();
            BoxContainer& visible_boxes = d_visible_boxes[ln][box_id];

            Box new_box(box, local_id, box.getOwnerRank());
            ++local_id;
            visible_boxes.insert(visible_boxes.end(), new_box);
         }
      }
   }
}

/*
 **************************************************************************
 * Destructor
 **************************************************************************
 */

FlattenedHierarchy::~FlattenedHierarchy()
{
}

}
}
