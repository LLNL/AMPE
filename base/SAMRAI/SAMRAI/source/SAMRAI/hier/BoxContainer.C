/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   A container of boxes with basic domain calculus operations
 *
 ************************************************************************/

#ifndef included_hier_BoxContainer_C
#define included_hier_BoxContainer_C

#include "SAMRAI/hier/BoxContainer.h"

#include "SAMRAI/hier/BoxContainerSingleBlockIterator.h"
#include "SAMRAI/hier/BoxContainerSingleOwnerIterator.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

const int BoxContainer::HIER_BOX_CONTAINER_VERSION = 0;


/*
 *************************************************************************
 * Constructors and destructor
 *************************************************************************
 */

BoxContainer::BoxContainer() :
   d_ordered(false)
{
}

BoxContainer::BoxContainer(
   const bool ordered) :
   d_ordered(ordered)
{
}

BoxContainer::BoxContainer(
   const Box& box,
   const bool ordered):
   d_list(1, box),
   d_ordered(false)
{
   if (ordered) {
      order();
   }
}

BoxContainer::BoxContainer(
   const BoxContainer& other):
   d_list(other.d_list),
   d_ordered(false) 
{
   if (other.d_ordered) {
      order();
   }
   if (other.d_tree) {
      makeTree();
   }
}

/*
 ***********************************************************************
 * Construct a BoxContainer consisting of the Boxes in other BoxContainer
 * with the given BlockId.
 ***********************************************************************
 */
BoxContainer::BoxContainer(
   const BoxContainer& other,
   const BlockId& block_id):
   d_ordered(false)
{
   BoxContainerSingleBlockIterator itr(other.begin(block_id));
   while (itr != other.end(block_id)) {
      const Box& mapped_box = *itr;
      pushBack(mapped_box);
      ++itr;
   }
   if (other.d_ordered) {
      order();
   }
}

/*
 *************************************************************************
 *
 * Construct BoxContainer from a range.
 *
 *************************************************************************
 */

BoxContainer::BoxContainer(
   const_iterator first,
   const_iterator last,
   const bool ordered):
   d_ordered(false)
{
   while (first != last) {
      pushBack(*first);
      ++first;
   }
   if (ordered) {
      order();
   }
}

/*
 *************************************************************************
 *
 * Construct from DatabaseBox array.
 * 
 *************************************************************************
 */

BoxContainer::BoxContainer(
   const tbox::Array<tbox::DatabaseBox>& other):
   d_ordered(false)
{
   const int n = other.size();
   for (int j = 0; j < n; j++) {
      pushBack(Box(other[j]));
   }
}

BoxContainer::~BoxContainer()
{
}

/*
 *************************************************************************
 *
 * Assignment.
 * 
 *************************************************************************
 */

BoxContainer&
BoxContainer::operator = (
   const BoxContainer& rhs)
{
   if (this != &rhs) {
      clear();
      d_list = rhs.d_list;
      if (rhs.d_ordered) { 
         order();
      } else {
         d_ordered = false;
      }
   }
   return *this;
}

BoxContainer&
BoxContainer::operator = (
   const tbox::Array<tbox::DatabaseBox>& rhs)
{
   clear();

   const int n = rhs.size();
   for (int j = 0; j < n; j++) {
      pushBack(Box(rhs[j]));
   }
   d_ordered = false;

   return *this;
}

/*
 *************************************************************************
 * Equality/inequality
 *************************************************************************
 */

bool
BoxContainer::operator == (
   const BoxContainer& rhs) const
{
   bool isEqual = (d_ordered && rhs.d_ordered); 
   if (isEqual) {
      isEqual = size() == rhs.size();
   }
   if (isEqual) {
      if (!d_ordered) { 
         isEqual = std::equal(d_list.begin(), d_list.end(),
                      rhs.d_list.begin(), Box::box_equality());
      } else {
         isEqual = std::equal(d_set.begin(), d_set.end(),
                      rhs.d_set.begin(), Box::id_equal());
      }
   }

   return isEqual;
}

/*
 *************************************************************************
 *
 * Insert Box.
 *
 *************************************************************************
 */

BoxContainer::iterator
BoxContainer::insert(
   iterator position,
   const Box& box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(box.getId().isValid());
   TBOX_ASSERT(box.getBlockId() != BlockId::invalidId());
   if (size() > 0) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(front(), box);
   }
#endif

   if (!d_ordered && size() == 0) {
      order();
   }

   if (!d_ordered) {
      TBOX_ERROR("insert attempted on unordered container.");
   }

   if (d_tree) {
      d_tree.reset();
   }

   const std::list<Box>::iterator& list_iter =
      d_list.insert(d_list.end(), box);

   iterator insert_iter;
   insert_iter.d_ordered = true;

   std::set<int>::size_type old_size = d_set.size();
   insert_iter.d_set_iter = d_set.insert(position.d_set_iter, &(*list_iter));
   if (d_set.size() == old_size) {
      d_list.erase(list_iter);
   } else {
      list_iter->lockId();
   }
   return insert_iter;
}

bool
BoxContainer::insert(
   const Box& box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(box.getId().isValid());
   if (size() > 0) {
      TBOX_DIM_ASSERT_CHECK_ARGS2(front(), box);
   }
#endif

   if (!d_ordered && size() == 0) {
      order();
   }

   if (!d_ordered) {
      TBOX_ERROR("BoxContainer::insert called on unordered BoxContainer.");
   }

   if (d_tree) {
      d_tree.reset();
   }

   const std::list<Box>::iterator& iter = d_list.insert(d_list.end(), box);
   Box* box_ptr(&(*iter));
   if (d_set.insert(box_ptr).second) {
      box_ptr->lockId();
      return true;
   } else {
      d_list.erase(iter);
      return false;
   }

}

/*
 *************************************************************************
 *
 * Insert a range.
 *
 *************************************************************************
 */

void
BoxContainer::insert (
   const_iterator first,
   const_iterator last )
{

   if (!d_ordered && size() == 0) {
      order();
   }

   if (!d_ordered) {
      TBOX_ERROR("insert attempted on unordered container.");
   }

   if (d_tree) {
      d_tree.reset();
   }

   for (std::set<Box*>::const_iterator set_iter = first.d_set_iter;
        set_iter != last.d_set_iter; ++set_iter) {

#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT((**set_iter).getId().isValid());
      if (size() > 0) {
         TBOX_DIM_ASSERT_CHECK_ARGS2(front(), **set_iter);
      }
#endif

      const std::list<Box>::iterator& list_iter =
         d_list.insert(d_list.end(), **set_iter);

      if (!d_set.insert(&(*list_iter)).second) {
         d_list.erase(list_iter);
      } else {
         list_iter->lockId();
      }
   }

}

/*
 ************************************************************************
 *
 * Function simplify() takes the complicated container of boxes and
 * coalesces regions together where possible.
 *
 * The canonical ordering for boxes is defined such that boxes which
 * lie next to each other in higher dimensions are coalesced together
 * before boxes which lie next to each other in lower dimensions.
 * Thus, we try to coalesce two boxes together on the higher
 * dimensions first.
 *
 * Assuming that two boxes a and b of dimension DIM are in canonical
 * order for dimensions d+1, ..., D, we can coalesce them together on
 * dimension d if:
 *
 *      (1) the lower and upper bounds for a and b agree for all
 *          dimensions greater than d
 *      (2) boxes a and b overlap or are next to each other in
 *          dimension d
 *      (3) boxes a and b overlap for all dimensions less than d
 *
 * If these conditions hold, then we break up the two boxes and put
 * them into the container of non-canonical boxes.
 *
 *************************************************************************
 */
void
BoxContainer::simplify()
{
   if (d_ordered) {
      TBOX_ERROR("simplify called on ordered BoxContainer.");
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   if (size() > 0) {
      const BlockId& front_block_id = front().getBlockId();
      for (const_iterator itr = begin(); itr != end(); ++itr) {
         TBOX_ASSERT(itr->getBlockId() == front_block_id);
      }
   }
#endif

   if (d_tree) {
      d_tree.reset();
   }

   // Start coalescing on the highest dimension of the containers and work down
   // While there are non-canonical boxes, pick somebody out of the container.

   if (!isEmpty()) {
      const tbox::Dimension dim(d_list.front().getDim());

      BoxContainer notCanonical;
      for (int d = dim.getValue() - 1; d >= 0; d--) {
         notCanonical.spliceBack(*this);
         while (!notCanonical.isEmpty()) {
            Box tryMe = notCanonical.front();
            notCanonical.popFront();

            // Pick somebody off of the canonical container and compare
            // against tryMe.

            if (!tryMe.empty()) {
               bool combineDaPuppies = false;
               iterator l(*this);
               for ( ; l != end(); ++l) {
                  const Box andMe = *l;

                  const Index& al = andMe.lower();
                  const Index& ah = andMe.upper();
                  const Index& bl = tryMe.lower();
                  const Index& bh = tryMe.upper();

                  combineDaPuppies = true;
                  for (int du = d + 1; du < dim.getValue(); du++) {
                     if ((al(du) != bl(du)) || (ah(du) != bh(du))) {
                        combineDaPuppies = false;
                        break;
                     }
                  }
                  if (combineDaPuppies) {
                     if ((bl(d) > ah(d) + 1) || (bh(d) < al(d) - 1)) {
                        combineDaPuppies = false;
                     } else {
                        for (int dl = 0; dl < d; dl++) {
                           if ((bl(dl) > ah(dl)) || (bh(dl) < al(dl))) {
                              combineDaPuppies = false;
                              break;
                           }
                        }
                     }
                  }
                  if (combineDaPuppies) {
                     break;
                  }
               }

               // If we are at the end of the canonical container, then just
               // add.  Otherwise, burst tryMe and andMe and put on
               // noncanonical.

               if (!combineDaPuppies) {
                  pushBack(tryMe);
               } else {
                  Box andMe = *l;
                  erase(l);
                  const Index& bl = tryMe.lower();
                  const Index& bh = tryMe.upper();
                  Index il = andMe.lower();
                  Index ih = andMe.upper();
                  for (int dl = 0; dl < d; dl++) {
                     if (il(dl) < bl(dl)) {
                        il(dl) = bl(dl);
                     }
                     if (ih(dl) > bh(dl)) {
                        ih(dl) = bh(dl);
                     }
                  }
                  if (bl(d) < il(d)) {
                     il(d) = bl(d);
                  }
                  if (bh(d) > ih(d)) {
                     ih(d) = bh(d);
                  }
                  Box intersection(il, ih, tryMe.getBlockId());
                  notCanonical.pushFront(intersection);
                  if (d > 0) {
                     notCanonical.burstBoxes(tryMe, intersection, d);
                     notCanonical.burstBoxes(andMe, intersection, d);
                  }
               }
            }
         }
      }
   }
}

/*
 *************************************************************************
 *
 * Coalesce boxes in the container where possible.  The resulting box
 * container will contain a non-overlapping set of boxes covering the
 * identical region of index space covered by the original container.
 * Two boxes may be coalesced if their union is a box (recall that union
 * is not closed over boxes), and they have a non-empty intersection or
 * they are adjacent to each other in index space.  Empty boxes in the
 * container are removed during this process.  Also, the boxes are
 * coalesced in the order in which they appear in the container.  No
 * attempt is made to coalesce boxes in any particular way (e.g., to
 * achieve the smallest number of boxes).
 *
 *************************************************************************
 */
void
BoxContainer::coalesce()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (size() > 0) {
      const BlockId& front_block_id = front().getBlockId();
      for (const_iterator itr = begin(); itr != end(); ++itr) {
         TBOX_ASSERT(itr->getBlockId() == front_block_id);
      }
   }
#endif

   if (d_ordered) {
      TBOX_ERROR("coalesce called on ordered BoxContainer.");
   }

   if (d_tree) {
      d_tree.reset();
   }

   iterator tb(*this);
   while (tb != end()) {

      bool found_match = false;

      iterator tb2 = tb;
      ++tb2;

      while (!found_match && tb2 != end()) {

         if (tb2->coalesceWith(*tb)) {
            found_match = true;
            erase(tb);
         }

         ++tb2;
      }

      if (found_match) {
         tb = begin();
      } else {
         ++tb;
      }
   }
}

/*
 *************************************************************************
 * Remove periodic images from container.
 *************************************************************************
 */

void
BoxContainer::removePeriodicImageBoxes()
{
   if (!d_ordered && size() > 0) {
      TBOX_ERROR("removePeriodicImages attempted on unordered container.");
   }

   if (d_tree) {
      d_tree.reset();
   }

   for (iterator na = begin(); na != end(); ) {
      if (na->isPeriodicImage()) {
         erase(na++);
      }
      else {
         ++na;
      }
   }
}

/*
 *************************************************************************
 * Separate periodic images from real boxes
 *************************************************************************
 */
void
BoxContainer::separatePeriodicImages(
   std::vector<Box>& real_mapped_box_vector,
   std::vector<Box>& periodic_image_mapped_box_vector) const
{
   if (!d_ordered) {
      TBOX_ERROR("separatePeriodicImages called on unordered BoxContainer.");
   }

   if (!isEmpty()) {

      const Box& first_element(*begin());

      const PeriodicId zero_shift_number(PeriodicShiftCatalog::getCatalog(
                                            first_element.getDim())->
                                         getZeroShiftNumber());

      real_mapped_box_vector.reserve(real_mapped_box_vector.size() + size());
      for (const_iterator ni = begin(); ni != end(); ++ni) {
         const Box& mapped_box = *ni;
         if (mapped_box.getPeriodicId() == zero_shift_number) {
            real_mapped_box_vector.push_back(mapped_box);
         } else {
            periodic_image_mapped_box_vector.push_back(mapped_box);
         }
      }
   }
}

/*
 *************************************************************************
 * Rotate all boxes
 *************************************************************************
 */

void
BoxContainer::rotate(
   const Transformation::RotationIdentifier rotation_ident)
{
   if (!isEmpty()) {

   if (d_tree) {
      d_tree.reset();
   }

      const tbox::Dimension& dim = d_list.front().getDim();
      const BlockId& block_id = d_list.front().getBlockId();
      if (dim.getValue() == 2 || dim.getValue() == 3) {
         for (iterator i(*this); i != end(); ++i) {
            if (i->getBlockId() != block_id) {
               TBOX_ERROR("BoxContainer::rotate() error ..."
                 << "\n  Attempted to rotate BoxContainer having Boxes with"
                 << "\n  differing BlockIds " << std::endl);
            } else {
               i->rotate(rotation_ident);
            }
         }
      } else {
         NULL_USE(rotation_ident);

         TBOX_ERROR("BoxContainer::rotate() error ..."
            << "\n   Rotation only implemented for 2D and 3D " << std::endl);
      }
   }
}


/*
 *************************************************************************
 *
 * Return the bounding box for all boxes in the BoxContainer.
 *
 *************************************************************************
 */
Box
BoxContainer::getBoundingBox() const
{
   if (isEmpty()) {
      TBOX_WARNING("Bounding box container is empty");
      const tbox::Dimension dim(tbox::Dimension::getInvalidDimension());
      Box empty(dim);
      return empty;
   } else {
      const_iterator i = begin();
      Box bbox(*i);
      const BlockId& block_id = bbox.getBlockId();
      ++i;
      for ( ; i != end(); ++i) {
         if (i->getBlockId() == block_id) {
            bbox += *i;
         } else {
            TBOX_ERROR("Attempted to find bounding box for BoxContainer with boxes from different blocks");
         }
      }
      return bbox;
   }
}


Box
BoxContainer::getBoundingBox(
   const BlockId& block_id) const
{
   if (isEmpty()) {
      TBOX_WARNING("Bounding box container is empty");
      const tbox::Dimension dim(tbox::Dimension::getInvalidDimension());
      Box empty(dim);
      return empty;
   } else {
      const tbox::Dimension& dim = d_list.front().getDim();
      Box bbox(dim);

      /*
       * First find the first box with the given BlockId
       */ 
      const_iterator i = begin();
      for ( ; i != end(); ++i) {
         if (i->getBlockId() == block_id) {
            bbox = *i;
            break;
         }
      }

      /*
       * If no boxes were found with the desired BlockId, then the returned
       * box will be empty.
       */
      if (i == end()) {
         TBOX_WARNING("Container has no boxes with the given BlockId");
      } else {

         ++i;

         for ( ; i != end(); ++i) {
            if (i->getBlockId() == block_id) {
               bbox += *i;
            }
         }
      }

      return bbox;
   }
}

/*
 *************************************************************************
 *
 * Test the box container for intersections among its boxes.
 *
 *************************************************************************
 */
bool
BoxContainer::boxesIntersect() const
{
   bool intersections = false;

   const_iterator tryMe(*this);
   const_iterator whatAboutMe(*this);
   ++whatAboutMe;
   while (!intersections && tryMe != end()) {
      while (!intersections && whatAboutMe != end()) {
         if (tryMe->getBlockId() == whatAboutMe->getBlockId()) {
            if (!(((*tryMe) * (*whatAboutMe)).size() == 0)) {
               intersections = true;
            }
         }
         ++whatAboutMe;
      }
      ++tryMe;
      whatAboutMe = tryMe;
      ++whatAboutMe;
   }
   return intersections;
}

/*
 *************************************************************************
 *
 * Remove intersections with a single box
 *
 *************************************************************************
 */
void
BoxContainer::removeIntersections(
   const Box& takeaway)
{
   if (isEmpty()) {
      return;
   }

   if (d_ordered) {
      TBOX_ERROR("removeIntersections attempted on ordered container.");
   }

   if (d_tree) {
      d_tree.reset();
   }

   const unsigned short dim = takeaway.getDim().getValue();
   iterator insertion_pt(*this);
   while (insertion_pt != end()) {
      Box& tryme = *insertion_pt;
      if (!tryme.intersects(takeaway)) {
         ++insertion_pt;
      } else {
         iterator tmp = insertion_pt;
         burstBoxes(tryme, takeaway, dim, insertion_pt);
         ++insertion_pt;
         erase(tmp);
      }
   }
}

/*
 *************************************************************************
 *
 * Remove intersections with another BoxContainer
 *
 *************************************************************************
 */
void
BoxContainer::removeIntersections(
   const BoxContainer& takeaway)
{
   if (isEmpty()) {
      return;
   }

   if (d_ordered) {
      TBOX_ERROR("removeIntersections attempted on ordered container.");
   }

   if (d_tree) {
      d_tree.reset();
   }

   if (takeaway.d_tree) {
      removeIntersections(*(takeaway.d_tree));
   } else {
      for (const_iterator remove(takeaway); remove != takeaway.end(); ++remove) {
         const Box& byebye = *remove;
         removeIntersections(byebye);
      }
   }
}

/*
 *************************************************************************
 *
 * Remove intersections with another BoxContainer using multiblock tree
 * representation.
 *
 *************************************************************************
 */
void
BoxContainer::removeIntersections(
   const IntVector& refinement_ratio,
   const BoxContainer& takeaway,
   const bool include_singularity_block_neighbors)
{
   if (isEmpty()) {
      return;
   }

   if (d_ordered) {
      TBOX_ERROR("removeIntersections attempted on ordered container.");
   }

   if (d_tree) {
      d_tree.reset();
   }

   if (takeaway.size() == 0) {
      return;
   }

   if (!takeaway.d_tree) {
      TBOX_ERROR("Multiblock removeIntersections called without makeTree.");
   }

   const BaseGridGeometry& grid_geometry(*takeaway.d_tree->getGridGeometry());

   std::vector<const Box *> overlap_mapped_boxes;
   iterator itr(*this);
   while (itr != end()) {
      const Box& tryme = *itr;
      takeaway.d_tree->findOverlapBoxes(overlap_mapped_boxes,
         tryme,
         refinement_ratio,
         include_singularity_block_neighbors);
      if (overlap_mapped_boxes.empty()) {
         ++itr;
      } else {
         iterator sublist_start = itr;
         iterator sublist_end = sublist_start;
         ++sublist_end;
         for (size_t i = 0;
              i < overlap_mapped_boxes.size() && sublist_start != sublist_end;
              ++i) {
            iterator insertion_pt = sublist_start;
            const BlockId& overlap_box_block_id =
               overlap_mapped_boxes[i]->getBlockId();
            if (overlap_box_block_id != sublist_start->getBlockId()) {
               Box overlap_box = *overlap_mapped_boxes[i];
               grid_geometry.transformBox(overlap_box,
                  refinement_ratio,
                  sublist_start->getBlockId(),
                  overlap_box_block_id);
               removeIntersectionsFromSublist(
                  overlap_box,
                  sublist_start,
                  sublist_end,
                  insertion_pt);
            } else {
               removeIntersectionsFromSublist(
                  *overlap_mapped_boxes[i],
                  sublist_start,
                  sublist_end,
                  insertion_pt);
            }
         }
         overlap_mapped_boxes.clear();
         itr = sublist_end;
      }
   }
}

void
BoxContainer::removeIntersections(
   const Box& box,
   const Box& takeaway)
{
   if (d_ordered) {
      TBOX_ERROR("removeIntersections attempted on ordered container.");
   }

   /*
    * The box container MUST be empty to use this function (see comments
    * in header file for discussion of why). If the two boxes intersect,
    * form a BoxContainer that contains the boxes resulting from removing
    * the intersection of box with takeaway.  If the two boxes do not
    * intersect, simply add box to the box container (no intersection removed).
    */
   TBOX_ASSERT(isEmpty());
   TBOX_ASSERT(box.getBlockId() == takeaway.getBlockId());

   if (d_tree) {
      d_tree.reset();
   }

   if (box.intersects(takeaway)) {
      burstBoxes(box, takeaway, box.getDim().getValue());
   } else {
      pushBack(box);
   }

}

void
BoxContainer::removeIntersections(
   const MultiblockBoxTree& takeaway)
{
   if (isEmpty()) {
      return;
   }

   if (d_ordered) {
      TBOX_ERROR("removeIntersections attempted on ordered container.");
   }

   if (d_tree) {
      d_tree.reset();
   }

   std::vector<const Box *> overlap_mapped_boxes;
   iterator itr(*this);
   while (itr != end()) {
      const Box& tryme = *itr;
      takeaway.findOverlapBoxes(overlap_mapped_boxes, tryme);
      if (overlap_mapped_boxes.empty()) {
         ++itr;
      } else {
         iterator sublist_start = itr;
         iterator sublist_end = sublist_start;
         ++sublist_end;
         for (size_t i = 0;
              i < overlap_mapped_boxes.size() && sublist_start != sublist_end;
              ++i) {
            iterator insertion_pt = sublist_start;
            removeIntersectionsFromSublist(
               *overlap_mapped_boxes[i],
               sublist_start,
               sublist_end,
               insertion_pt);
         }
         overlap_mapped_boxes.clear();
         itr = sublist_end;
      }
   }
}


void
BoxContainer::removeIntersectionsFromSublist(
   const Box& takeaway,
   iterator& sublist_start,
   iterator& sublist_end,
   iterator& insertion_pt)
{
   if (d_ordered) {
      TBOX_ERROR("removeIntersections attempted on ordered container.");
   }

   if (d_tree) {
      d_tree.reset();
   }

   const unsigned short dim = takeaway.getDim().getValue();
   iterator itr = sublist_start;
   while (itr != sublist_end) {
      Box& tryme = *itr;
      if (!tryme.intersects(takeaway)) {
         ++itr;
      } else {
         burstBoxes(tryme, takeaway, dim, insertion_pt);
         iterator tmp = itr;
         ++itr;
         if (tmp == sublist_start) {
            ++sublist_start;
         }
         erase(tmp);
      }
      insertion_pt = itr;
   }
}

/*
 *************************************************************************
 *
 * Return the boxes in the current container that intersect the index
 * space of the argument.
 *
 *************************************************************************
 */
void
BoxContainer::intersectBoxes(
   const Box& keep)
{
   if (isEmpty()) {
      return;
   }

   if (d_ordered) {
      TBOX_ERROR("intersectBoxes attempted on ordered container.");
   }

   if (d_tree) {
      d_tree.reset();
   }

   iterator i(*this);
   Box overlap(i->getDim());
   while (i != end()) {
      Box& tryMe = *i;
      tryMe.intersect(keep, overlap);
      if (!overlap.empty()) {
         tryMe = overlap;
         ++i;
      } else {
         iterator tmp = i;
         ++i;
         erase(tmp);
      }
   }
}

void
BoxContainer::intersectBoxes(
   const BoxContainer& keep)
{
   if (isEmpty()) {
      return;
   }

   if (d_ordered) {
      TBOX_ERROR("intersectBoxes attempted on ordered container.");
   }

   if (d_tree) {
      d_tree.reset();
   }

   if (keep.d_tree) {
      intersectBoxes(*(keep.d_tree));
   } else {
      iterator insertion_pt(*this);
      Box overlap(insertion_pt->getDim());
      while (insertion_pt != end()) {
         iterator tmp = insertion_pt;
         const Box& tryme = *insertion_pt;
         for (const_iterator i(keep); i != keep.end(); ++i) {
            tryme.intersect(*i, overlap);
            if (!overlap.empty()) {
               insertAfter(insertion_pt, overlap);
               ++insertion_pt;
            }
         }
         ++insertion_pt;
         erase(tmp);
      }
   }
}

void
BoxContainer::intersectBoxes(
   const MultiblockBoxTree& keep)
{
   if (isEmpty()) {
      return;
   }

   if (d_ordered) {
      TBOX_ERROR("intersectBoxes attempted on ordered container.");
   }

   if (d_tree) {
      d_tree.reset();
   }

   std::vector<const Box *> overlap_mapped_boxes;
   Box overlap(front().getDim());
   iterator itr(*this);
   iterator insertion_pt = itr;
   while (itr != end()) {
      const Box& tryme = *itr;
      keep.findOverlapBoxes(overlap_mapped_boxes, tryme);
      for (size_t i = 0; i < overlap_mapped_boxes.size(); ++i) {
         tryme.intersect(*overlap_mapped_boxes[i], overlap);
         if (!overlap.empty()) {
            insertAfter(insertion_pt, overlap);
            ++insertion_pt;
         }
      }
      overlap_mapped_boxes.clear();
      iterator tmp = itr;
      ++insertion_pt;
      itr = insertion_pt;
      erase(tmp);
   }
}

void
BoxContainer::intersectBoxes(
   const IntVector& refinement_ratio,
   const BoxContainer& keep,
   const bool include_singularity_block_neighbors)
{
   if (isEmpty()) {
      return;
   }

   if (d_ordered) {
      TBOX_ERROR("intersectBoxes attempted on ordered container.");
   }

   if (d_tree) {
      d_tree.reset();
   }

   if (keep.size() == 0) {
      clear();
      return; 
   }

   if (!keep.d_tree) {
      TBOX_ERROR("Multiblock intersectBoxes called without makeTree.");
   }

   const BaseGridGeometry& grid_geometry(*keep.d_tree->getGridGeometry());

   std::vector<const Box *> overlap_mapped_boxes;
   Box overlap(front().getDim());
   iterator itr(*this);
   iterator insertion_pt = itr;
   while (itr != end()) {
      const Box& tryme = *itr;
      keep.d_tree->findOverlapBoxes(overlap_mapped_boxes,
         tryme,
         refinement_ratio,
         include_singularity_block_neighbors);
      for (size_t i = 0; i < overlap_mapped_boxes.size(); ++i) {
         const BlockId& overlap_box_block_id =
            overlap_mapped_boxes[i]->getBlockId();
         if (overlap_box_block_id != tryme.getBlockId()) {
            Box overlap_box = *overlap_mapped_boxes[i];
            grid_geometry.transformBox(overlap_box,
               refinement_ratio,
               tryme.getBlockId(),
               overlap_box_block_id);
            tryme.intersect(overlap_box, overlap);
            if (!overlap.empty()) {
               insertAfter(insertion_pt, overlap);
               ++insertion_pt;
            }
         } else {
            tryme.intersect(*overlap_mapped_boxes[i], overlap);
            if (!overlap.empty()) {
               insertAfter(insertion_pt, overlap);
               ++insertion_pt;
            }
         }
      }
      overlap_mapped_boxes.clear();
      iterator tmp = itr;
      ++insertion_pt;
      itr = insertion_pt;
      erase(tmp);
   }
}

/*
 *************************************************************************
 *
 * Type conversion from a BoxContainer to an Array of
 * tbox::DatabaseBoxes.
 *
 *************************************************************************
 */
BoxContainer::operator tbox::Array<tbox::DatabaseBox>() const
{
   tbox::Array<tbox::DatabaseBox> new_Array(size());

   int j = 0;
   for (const_iterator i = begin(); i != end(); ++i) {
      new_Array[j++] = (tbox::DatabaseBox)(*i);
   }

   return new_Array;
}

/*
 *************************************************************************
 *
 * Break up box bursty against box solid and adds the pieces to
 * container.  The bursting is done on dimensions 0 through dimension-1,
 * starting with lowest dimensions first to try to maintain the
 * canonical representation for the bursted domains.
 *
 *************************************************************************
 */
void
BoxContainer::burstBoxes(
   const Box& bursty,
   const Box& solid,
   const int dimension)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(bursty, solid);
   TBOX_ASSERT(dimension <= bursty.getDim().getValue());

   // Set up the lower and upper bounds of the regions for ease of access

   Index burstl = bursty.lower();
   Index bursth = bursty.upper();
   const Index& solidl = solid.lower();
   const Index& solidh = solid.upper();
   const BlockId& block_id = bursty.getBlockId();

   // Break bursty region against solid region along low dimensions first

   for (int d = 0; d < dimension; d++) {
      if (bursth(d) > solidh(d)) {
         Index newl = burstl;
         newl(d) = solidh(d) + 1;
         pushBack(Box(newl, bursth, block_id));
         bursth(d) = solidh(d);
      }
      if (burstl(d) < solidl(d)) {
         Index newh = bursth;
         newh(d) = solidl(d) - 1;
         pushBack(Box(burstl, newh, block_id));
         burstl(d) = solidl(d);
      }
   }
}

/*
 *************************************************************************
 *
 * Break up box bursty against box solid and adds the pieces to
 * container starting at the location pointed to by the supplied
 * iterator.  The bursting is done on dimensions 0 through dimension-1,
 * starting with lowest dimensions first to try to maintain the
 * canonical representation for the bursted domains.
 *
 *************************************************************************
 */
void
BoxContainer::burstBoxes(
   const Box& bursty,
   const Box& solid,
   const int dimension,
   iterator& insertion_pt)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(bursty, solid);
   TBOX_ASSERT(dimension <= bursty.getDim().getValue());

   // Set up the lower and upper bounds of the regions for ease of access

   Index burstl = bursty.lower();
   Index bursth = bursty.upper();
   const Index& solidl = solid.lower();
   const Index& solidh = solid.upper();
   const BlockId& block_id = bursty.getBlockId();

   // Break bursty region against solid region along low dimensions first

   for (int d = 0; d < dimension; d++) {
      if (bursth(d) > solidh(d)) {
         Index newl = burstl;
         newl(d) = solidh(d) + 1;
         insertAfter(insertion_pt, Box(newl, bursth, block_id));
         bursth(d) = solidh(d);
         ++insertion_pt;
      }
      if (burstl(d) < solidl(d)) {
         Index newh = bursth;
         newh(d) = solidl(d) - 1;
         insertAfter(insertion_pt, Box(burstl, newh, block_id));
         burstl(d) = solidl(d);
         ++insertion_pt;
      }
   }
}

/*
 ***********************************************************************
 * Insert Box owners into a single set container.
 ***********************************************************************
 */
void
BoxContainer::getOwners(
   std::set<int>& owners) const
{
   if (d_ordered) {
      for (const_iterator i_nabr = begin();
           i_nabr != end(); ++i_nabr) {
         const int owner = i_nabr->getOwnerRank();
         owners.insert(owner);
      }
   } else {
      for (const_iterator i_nabr = begin();
           i_nabr != end(); ++i_nabr) {
         if (i_nabr->getId().isValid()) {
            const int owner = i_nabr->getOwnerRank();
            owners.insert(owner);
         } else {
            TBOX_ERROR("Attempted to get owner of Box with invalid BoxId.");
         }
      }
   }
}

/*
 ***********************************************************************
 * Unshift periodic image Boxes from a BoxContainer.
 ***********************************************************************
 */
void
BoxContainer::unshiftPeriodicImageBoxes(
   BoxContainer& output_mapped_boxes,
   const IntVector& refinement_ratio) const
{
   if (!d_ordered) {
      TBOX_ERROR("unshiftPeriodicImageBoxes called on unordered container.");
   }

   iterator hint = output_mapped_boxes.begin();

   if (!isEmpty()) {
      const Box& first_element(*begin());

      const PeriodicId zero_shift_number(PeriodicShiftCatalog::getCatalog(
                                            first_element.getDim())->
                                         getZeroShiftNumber());

      for (const_iterator na = begin(); na != end(); ++na) {
         if (na->isPeriodicImage()) {
            const Box unshifted_mapped_box(
               *na, zero_shift_number, refinement_ratio);
            hint = output_mapped_boxes.insert(hint, unshifted_mapped_box);
         } else {
            hint = output_mapped_boxes.insert(hint, *na);
         }
      }
   }
}

/*
 ***********************************************************************
 * Switch to ordered state.
 ***********************************************************************
 */
void
BoxContainer::order()
{
   if (!d_ordered) {
      d_set.clear();
      for (iterator i(*this); i != end(); ++i) {
         if (!i->getId().isValid()) {
            TBOX_ERROR("Attempted to order a BoxContainer that has a member with an invalid BoxId.");
         }
         if (d_set.insert(&(*i)).second == false) {
            TBOX_ERROR("Attempted to order a BoxContainer with duplicate BoxIds.");
         }
         i->lockId();
      }
      d_ordered = true;
   }
}

/*
 ***********************************************************************
 * Switch to unordered state.
 ***********************************************************************
 */
void
BoxContainer::unorder()
{
   if (d_ordered) {
      d_set.clear();
      d_ordered = false;
   }
}

/*
 *************************************************************************
 * Erase methods
 *************************************************************************
 */

void
BoxContainer::erase(
   iterator iter)
{
   if (!d_ordered) {
      d_list.erase(iter.d_list_iter);
   } else {
      const Box& box = **(iter.d_set_iter);
      d_set.erase(iter.d_set_iter);

      for (std::list<Box>::iterator bi = d_list.begin(); bi != d_list.end();
           ++bi) {
         if (bi->getId() == box.getId()) {
            d_list.erase(bi);
            break;
         }
      }
   }
   if (d_tree) {
      d_tree.reset();
   }
}

void
BoxContainer::erase(
   iterator first,
   iterator last)
{
   if (!d_ordered) {
      d_list.erase(first.d_list_iter, last.d_list_iter);
   } else {
      for (iterator iter = first; iter != last; ++iter) {
         erase(iter);
      }
   }
   if (d_tree) {
      d_tree.reset();
   }
}


int
BoxContainer::erase(
   const Box& box)
{
   if (!d_ordered) {
      TBOX_ERROR("erase with Box argument attempted on unordered BoxContainer.");
   }

   int ret = static_cast<int>(d_set.erase(const_cast<Box*>(&box)));
   for (std::list<Box>::iterator bi = d_list.begin(); bi != d_list.end();
        ++bi) {
      if (bi->getId() == box.getId()) {
         d_list.erase(bi++);
         break;
      }
   }

   if (d_tree) {
      d_tree.reset();
   }

   return ret;
}

/*
 *************************************************************************
 * Box-based queries
 *************************************************************************
 */

int
BoxContainer::getTotalSizeOfBoxes() const
{
   int size = 0;
   if (!d_ordered) {
      for (const_iterator i(*this); i != end(); ++i) {
         size += i->size();
      }
   } else {
      for (const_iterator i(*this); i != end(); ++i) {
         size += i->size();
      }
   }
   return size;
}

bool
BoxContainer::contains(
   const Index& idx,
   const BlockId& block_id) const
{
   for (const_iterator i(*this); i != end(); ++i) {
      //TODO: Change this when BoxContainer can no longer accept invalid BlockId
      if (i->getBlockId() == block_id ||
          i->getBlockId() == BlockId::invalidId()) {
         if (i->contains(idx)) {
            return true;
         }
      }
   }
   return false;
}

/*
 *************************************************************************
 * Spatial manipulation of Boxes
 ************************************************************************
 */

void
BoxContainer::grow(
   const IntVector& ghosts)
{
   for (iterator i(*this); i != end(); ++i) {
      i->grow(ghosts);
   }

   if (d_tree) {
      d_tree.reset();
   }
}


void
BoxContainer::shift(
   const IntVector& offset)
{
   for (iterator i(*this); i != end(); ++i) {
      i->shift(offset);
   }
   if (d_tree) {
      d_tree.reset();
   }
}

void
BoxContainer::refine(
   const IntVector& ratio)
{
   for (iterator i(*this); i != end(); ++i) {
      i->refine(ratio);
   }
   if (d_tree) {
      d_tree.reset();
   }
}

void
BoxContainer::coarsen(
   const IntVector& ratio)
{
   for (iterator i(*this); i != end(); ++i) {
      i->coarsen(ratio);
   }
   if (d_tree) {
      d_tree.reset();
   }
}

void
BoxContainer::makeTree(
   const BaseGridGeometry* grid_geometry,
   const int min_number) const
{
   TBOX_ASSERT(min_number > 0);

   if (!d_tree && size() > 0) {
      d_tree.reset(new MultiblockBoxTree(*this, grid_geometry, min_number));
   }
}

bool
BoxContainer::hasBoxInBlock(
   const BlockId& block_id) const
{
   if (d_tree) {
      return d_tree->hasBoxInBlock(block_id);
   } else {
      for (const_iterator bi = begin(); bi != end(); ++bi) {
         if (bi->getBlockId() == block_id) {
            return true;
         }
      }
      return false;
   }
}


/*
 ***********************************************************************
 * Write the BoxContainer to a database.
 ***********************************************************************
 */
void
BoxContainer::putUnregisteredToDatabase(
   const boost::shared_ptr<tbox::Database>& database) const
{
   database->putInteger(
      "HIER_BOX_CONTAINER_VERSION", HIER_BOX_CONTAINER_VERSION);

   const int mbs_size = size();
   database->putInteger("mapped_box_set_size", mbs_size);
   if (mbs_size > 0) {

      std::vector<int> local_ids;
      std::vector<int> ranks;
      std::vector<int> block_ids;
      std::vector<int> periodic_ids;
      local_ids.reserve(mbs_size);
      ranks.reserve(mbs_size);
      block_ids.reserve(mbs_size);
      periodic_ids.reserve(mbs_size);

      tbox::Array<tbox::DatabaseBox> db_box_array(mbs_size);

      int counter = -1;
      for (BoxContainer::const_iterator ni = begin(); ni != end(); ++ni) {
         local_ids.push_back(ni->getLocalId().getValue());
         ranks.push_back(ni->getOwnerRank());
         block_ids.push_back(ni->getBlockId().getBlockValue());
         periodic_ids.push_back(ni->getPeriodicId().getPeriodicValue());
         db_box_array[++counter] = *ni;
      }

      database->putIntegerArray("local_indices", &local_ids[0], mbs_size);
      database->putIntegerArray("ranks", &ranks[0], mbs_size);
      database->putIntegerArray("block_ids", &block_ids[0], mbs_size);
      database->putIntegerArray("periodic_ids", &periodic_ids[0], mbs_size);
      database->putDatabaseBoxArray("boxes", &db_box_array[0], mbs_size);
   }
}


/*
 ***********************************************************************
 * Read the BoxContainer from a database.
 ***********************************************************************
 */
void
BoxContainer::getFromDatabase(
   tbox::Database& database)
{
   const unsigned int mbs_size = database.getInteger("mapped_box_set_size");
   if (mbs_size > 0) {
      std::vector<int> local_ids(mbs_size);
      std::vector<int> ranks(mbs_size);
      std::vector<int> block_ids(mbs_size);
      std::vector<int> periodic_ids(mbs_size);
      tbox::Array<tbox::DatabaseBox> db_box_array(mbs_size);

      database.getIntegerArray(
         "local_indices", &local_ids[0], mbs_size);
      database.getIntegerArray(
         "ranks", &ranks[0], mbs_size);
      database.getIntegerArray(
         "block_ids", &block_ids[0], mbs_size);
      database.getIntegerArray(
         "periodic_ids", &periodic_ids[0], mbs_size);
      database.getDatabaseBoxArray(
         "boxes", &db_box_array[0], mbs_size);

      for (unsigned int i = 0; i < mbs_size; ++i) {
         Box box(db_box_array[i]);
         box.setBlockId(BlockId(block_ids[i]));
         Box mapped_box(
            box,
            LocalId(local_ids[i]),
            ranks[i],
            PeriodicId(periodic_ids[i]));
         insert(end(), mapped_box);
      }
   }
}


/*
 ***********************************************************************
 * Prind contents of the BoxContainer
 ***********************************************************************
 */
void
BoxContainer::print(
   std::ostream& co) const
{
   co << size() << " boxes, " << (d_ordered?"ordered":"unordered") << '\n';
   for (const_iterator bi = begin(); bi != end(); ++bi) {
      Box box(*bi);
      co << "    "
         << box << "   "
         << box.numberCells() << '\n';
   }
}

/*
 ***********************************************************************
 * Construct a BoxContainer Outputter with formatting parameters.
 ***********************************************************************
 */

BoxContainer::Outputter::Outputter(
   const BoxContainer& mapped_box_set,
   const std::string& border,
   int detail_depth):
   d_set(mapped_box_set),
   d_border(border),
   d_detail_depth(detail_depth)
{
}

/*
 ***********************************************************************
 * Print out a BoxContainer according to settings in the Outputter.
 ***********************************************************************
 */

std::ostream&
operator << (
   std::ostream& s,
   const BoxContainer::Outputter& format)
{
   format.d_set.print(s);
   return s;
}


/*
 ***********************************************************************
 * Return a Outputter that can dump the BoxContainer to a stream.
 ***********************************************************************
 */

BoxContainer::Outputter
BoxContainer::format(
   const std::string& border,
   int detail_depth) const
{
   return Outputter(*this, border, detail_depth);
}

void
BoxContainer::findOverlapBoxes(
   BoxContainer& container,
   const Box& box) const
{
   if (d_tree) {
      d_tree->findOverlapBoxes(container, box);
   } else {
      if (container.isOrdered()) {
         for (const_iterator ni = begin(); ni != end(); ++ni) {
            const Box& my_box = *ni;
            if (box.intersects(my_box)) {
               container.insert(container.end(), my_box);
            }
         }
      } else {
         for (const_iterator ni = begin(); ni != end(); ++ni) {
            const Box& my_box = *ni;
            if (box.intersects(my_box)) {
               container.pushBack(my_box);
            }
         }
      }
   }
}

void
BoxContainer::findOverlapBoxes(
   std::vector<const Box*>& box_vector,
   const Box& box) const
{
   if (d_tree) {
      TBOX_ASSERT(d_tree->getNumberBlocksInTree() == 1);
      d_tree->findOverlapBoxes(box_vector, box);
   } else {
      for (const_iterator ni = begin(); ni != end(); ++ni) {
         const Box& my_box = *ni;
         if (box.intersects(my_box)) {
            box_vector.push_back(&my_box);
         }
      }
   }
}

void
BoxContainer::findOverlapBoxes(
   BoxContainer& overlap_boxes,
   const Box& box,
   const IntVector& refinement_ratio,
   bool include_singularity_block_neighbors) const
{
   if (isEmpty()) {
      return;
   }

   if (!d_tree) {
      TBOX_ERROR("Must call makeTree before calling findOverlapBoxes with refinement ratio argument.");
   }

   d_tree->findOverlapBoxes(overlap_boxes,
                            box,
                            refinement_ratio,
                            include_singularity_block_neighbors); 
}


bool
BoxContainer::hasOverlap(
   const Box& box) const
{
   if (d_tree) {
      return (d_tree->hasOverlap(box));
   } else {
      bool ret_val = false;      
      for (const_iterator ni = begin(); ni != end(); ++ni) {
         const Box& my_box = *ni;
         if (box.intersects(my_box)) {
            ret_val = true;
            break;
         } 
      }
      return ret_val;
   } 
}

BoxContainerSingleBlockIterator
BoxContainer::begin(
   const BlockId& block_id) const
{
   return BoxContainerSingleBlockIterator(*this, block_id, true);
}

BoxContainerSingleBlockIterator
BoxContainer::end(
   const BlockId& block_id) const
{
   return BoxContainerSingleBlockIterator(*this, block_id, false);
}

BoxContainerSingleOwnerIterator
BoxContainer::begin(
   const int& owner_rank) const
{
   return BoxContainerSingleOwnerIterator(*this, owner_rank, true);
}

BoxContainerSingleOwnerIterator
BoxContainer::end(
   const int& owner_rank) const
{
   return BoxContainerSingleOwnerIterator(*this, owner_rank, false);
}

RealBoxConstIterator
BoxContainer::realBegin() const
{
   return RealBoxConstIterator(*this, true);
}

RealBoxConstIterator
BoxContainer::realEnd() const
{
   return RealBoxConstIterator(*this, false);
}

BoxContainer::BoxContainerIterator::BoxContainerIterator(
   BoxContainer& container,
   bool from_start):
   d_list_iter(from_start ? container.d_list.begin() :
               container.d_list.end()),
   d_set_iter(from_start ? container.d_set.begin() :
               container.d_set.end()),
   d_ordered(container.d_ordered)
{
}

BoxContainer::BoxContainerIterator::BoxContainerIterator(
   const BoxContainerIterator& other):
   d_list_iter(other.d_list_iter),
   d_set_iter(other.d_set_iter),
   d_ordered(other.d_ordered)
{
}

BoxContainer::BoxContainerIterator::BoxContainerIterator():
   d_ordered(false)
{
}

BoxContainer::BoxContainerIterator::~BoxContainerIterator()
{
}

BoxContainer::BoxContainerConstIterator::BoxContainerConstIterator(
   const BoxContainer& container,
   bool from_start):
   d_list_iter(from_start ? container.d_list.begin() :
               container.d_list.end()),
   d_set_iter(from_start ? container.d_set.begin() :
              container.d_set.end()),
   d_ordered(container.d_ordered)
{
}

BoxContainer::BoxContainerConstIterator::BoxContainerConstIterator(
   const BoxContainerConstIterator& other):
   d_list_iter(other.d_list_iter),
   d_set_iter(other.d_set_iter),
   d_ordered(other.d_ordered)
{
}

BoxContainer::BoxContainerConstIterator::BoxContainerConstIterator(
   const BoxContainerIterator& other):
   d_list_iter(other.d_list_iter),
   d_set_iter(other.d_set_iter),
   d_ordered(other.d_ordered)
{
}

BoxContainer::BoxContainerConstIterator::~BoxContainerConstIterator()
{
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
