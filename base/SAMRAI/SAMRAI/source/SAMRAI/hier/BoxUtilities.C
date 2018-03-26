/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Routines for processing boxes within a domain of index space.
 *
 ************************************************************************/

#ifndef included_hier_BoxUtilities_C
#define included_hier_BoxUtilities_C

#include "SAMRAI/hier/BoxUtilities.h"

#include "SAMRAI/hier/BoxContainer.h" 
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"


namespace SAMRAI {
namespace hier {

/*
 *************************************************************************
 *
 * This static private member function is called by findBadCutPoints(),
 * and the findBadCutPointsForDirection() member functions.  It sets bad
 * cut points near the lower and upper ends of the border box in the
 * given coordinate direction.
 *
 *************************************************************************
 */

void
BoxUtilities::findBadCutPointsForBorderAndDirection(
   const int id,
   tbox::Array<bool>& bad_cuts,
   const Box& box,
   const Box& border,
   const int bad_interval)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, border);

   TBOX_ASSERT((0 <= id) && (id < box.getDim().getValue()));
   TBOX_ASSERT(bad_cuts.getSize() == box.numberCells(id));
   TBOX_ASSERT(bad_interval >= 0);

   if (bad_interval > 0) {

      const int ilo = box.lower(id);
      const int ihi = box.upper(id);

      int iclo, ichi, ic;

      /*
       * Set bad cut points near lower end of border box.
       */
      int mark = border.lower(id);
      if (mark > (ilo - bad_interval)) {

         iclo =
            tbox::MathUtilities<int>::Max(ilo, (mark - bad_interval + 1)) - ilo;
         ichi =
            tbox::MathUtilities<int>::Min(ihi, (mark - 1)) - ilo + 1;
         for (ic = iclo; ic < ichi; ic++) bad_cuts[ic] = true;

         iclo =
            tbox::MathUtilities<int>::Max(ilo, (mark + 1)) - ilo;
         ichi =
            tbox::MathUtilities<int>::Min(ihi,
               (mark + bad_interval - 1)) - ilo + 1;
         for (ic = iclo; ic < ichi; ic++) bad_cuts[ic] = true;

      }

      /*
       * Set bad cut points near upper end of border box.
       */
      mark = border.upper(id) + 1;
      if (mark < (ihi + bad_interval + 1)) {

         iclo =
            tbox::MathUtilities<int>::Max(ilo, (mark - bad_interval + 1)) - ilo;
         ichi =
            tbox::MathUtilities<int>::Min(ihi, (mark - 1)) - ilo + 1;
         for (ic = iclo; ic < ichi; ic++) bad_cuts[ic] = true;

         iclo =
            tbox::MathUtilities<int>::Max(ilo, (mark + 1)) - ilo;
         ichi =
            tbox::MathUtilities<int>::Min(ihi,
               (mark + bad_interval - 1)) - ilo + 1;
         for (ic = iclo; ic < ichi; ic++) bad_cuts[ic] = true;

      }

   }
}

/*
 *************************************************************************
 *
 * Check min size, cut factor, and physical domain constraints for
 * given box.  If a patch is generated from a box that violates any
 * of these constraints, then some other routine (e.g., ghost cell
 * filling, or inter-patch communication) may fail.  Thus, an error
 * message will be generated describing the violation and the program
 * will abort.
 *
 *************************************************************************
 */

void
BoxUtilities::checkBoxConstraints(
   const Box& box,
   const IntVector& min_size,
   const IntVector& cut_factor,
   const IntVector& bad_interval,
   const BoxContainer& physical_boxes)
{

   TBOX_DIM_ASSERT_CHECK_ARGS3(min_size, cut_factor, bad_interval);

   TBOX_ASSERT(min_size > IntVector::getZero(min_size.getDim()));
   TBOX_ASSERT(cut_factor > IntVector::getZero(min_size.getDim()));
   TBOX_ASSERT(bad_interval >= IntVector::getZero(min_size.getDim()));

   const tbox::Dimension& dim(box.getDim());

   int id;

   /*
    * Test box against minimum size constraint.
    */
   tbox::Array<bool> min_is_bad(dim.getValue());
   bool min_violation = false;
   for (id = 0; id < dim.getValue(); id++) {
      if (box.numberCells(id) < min_size(id)) {
         min_is_bad[id] = true;
         min_violation = true;
      } else {
         min_is_bad[id] = false;
      }
   }

   if (min_violation) {
      tbox::perr << "\nBox = " << box << " -- minimum size = " << min_size
                 << std::endl;
      for (id = 0; id < dim.getValue(); id++) {
         if (min_is_bad[id]) {
            tbox::perr << "min size violated in direction " << id << std::endl;
         }
      }
      TBOX_ERROR("BoxUtilities::checkBoxConstraints() error:\n"
         << "  Box violates minimum size restriction" << std::endl);
   }

   /*
    * Test box against cut factor constraint.
    */
   tbox::Array<bool> factor_is_bad(dim.getValue());
   bool factor_violation = false;
   for (id = 0; id < dim.getValue(); id++) {
      if ((box.numberCells(id) % cut_factor(id)) != 0) {
         factor_is_bad[id] = true;
         factor_violation = true;
      } else {
         factor_is_bad[id] = false;
      }
   }

   if (factor_violation) {
      tbox::perr << "\nBox = " << box << " -- cut factor = " << cut_factor
                 << std::endl;
      for (id = 0; id < dim.getValue(); id++) {
         if (factor_is_bad[id]) {
            tbox::perr << "factor bad in direction " << id << std::endl;
         }
      }
      TBOX_ERROR("BoxUtilities::checkBoxConstraints() error:\n"
         << "  Box violates cut factor restriction" << std::endl);
   }

   if (physical_boxes.size() > 0) {

      tbox::Array<bool> cut_is_bad(dim.getValue());
      for (id = 0; id < dim.getValue(); id++) {
         cut_is_bad[id] = false;
      }

      bool bad_cut_violation = false;

      /*
       * Test box for bad cut point violation.
       */

      Box test_border = box;
      test_border.grow(bad_interval);

      BoxContainer border_boxes(test_border);
      border_boxes.removeIntersections(physical_boxes);

      if (!border_boxes.isEmpty()) {

         /*
          * Test individual box faces in each direction for bad cuts.
          */

         id = 0;
         while ((id < dim.getValue()) && !bad_cut_violation) {

            int blo = box.lower(id);
            int bhi = box.upper(id);
            int bad = bad_interval(id);

            /*
             * Test lower box face in single direction.
             */

            Box test_box = box;
            test_box.grow(bad_interval);

            test_box.upper(id) = box.lower(id) - 1;

            BoxContainer test_boxes(test_box);
            test_boxes.intersectBoxes(border_boxes);
            test_boxes.simplify();

            BoxContainer::iterator tb = test_boxes.begin();
            while (!bad_cut_violation && tb != test_boxes.end()) {
               if ((tb->lower(id) > (blo - bad))
                   || (tb->upper(id) < (blo - 1))) {
                  bad_cut_violation = true;
                  cut_is_bad[id] = true;
               }
               ++tb;
            }

            if (!bad_cut_violation) {

               /*
                * Test upper box face in single direction.
                */

               test_box = box;
               test_box.grow(bad_interval);

               test_box.lower(id) = box.upper(id) + 1;

               test_boxes = BoxContainer(test_box); 
               test_boxes.intersectBoxes(border_boxes);
               test_boxes.simplify();

               tb = test_boxes.begin();
               while (!bad_cut_violation && tb != test_boxes.end()) {
                  if ((tb->lower(id) > (bhi + 1))
                      || (tb->upper(id) < (bhi + bad))) {
                     bad_cut_violation = true;
                     cut_is_bad[id] = true;
                  }
                  ++tb;
               }

            }

            id++;
         }

      }

      if (bad_cut_violation) {

         tbox::perr << "Box violates bad cut restriction in directions...";
         for (id = 0; id < dim.getValue(); id++) {
            if (cut_is_bad[id]) tbox::perr << "\n" << id;
         }
         tbox::perr << "\nBox = " << box << " -- bad cut interval = "
                    << bad_interval << std::endl;
         tbox::perr << "Physical domain boxes ... " << std::endl;
         int ib = 0;
         for (BoxContainer::const_iterator itr(physical_boxes);
              itr != physical_boxes.end(); ++itr, ++ib) {
            tbox::perr << "Box # " << ib << " -- " << *itr << std::endl;
         }
         TBOX_ERROR("BoxUtilities::checkBoxConstraints() error:\n"
            << "  Box violates bad cut restriction" << std::endl);
      }

   }
}

/*
 *************************************************************************
 *
 * Replace each box in the list that is too large with a list of
 * nonoverlapping smaller boxes whose union covers the same region of
 * index space as the original box.  The resulting boxes will obey the
 * minimum size, and cut factor restrictions if the original box does.
 * However, the maximum size restriction may be sacrified if the box
 * cannot be chopped at appropriate points.
 *
 * For each box in the list, we perform the following operations
 *
 *    (1) Determine a set of cut points for each coordinate direction.
 *        The ideal cuts satisfy all min, max, and factor restrictions
 *        assuming the box does too.
 *
 *    (2) If step (1) finds that the box may be chopped, we determine
 *        the bad cut points for the box and adjust the original cut
 *        points if necessary.  Note that this operation uses the
 *        physical domain and the bad interval information.
 *
 *    (3) The box is chopped if this is still possible after (1) and (2).
 *
 *    (4) If the box is chopped, set the box list to the resulting
 *        boxes.  Otherwise, put the original box on the list.
 *
 *************************************************************************
 */

void
BoxUtilities::chopBoxes(
   BoxContainer& boxes,
   const IntVector& max_size,
   const IntVector& min_size,
   const IntVector& cut_factor,
   const IntVector& bad_interval,
   const BoxContainer& physical_boxes)
{
   TBOX_DIM_ASSERT_CHECK_ARGS4(max_size, min_size, cut_factor, bad_interval);

   TBOX_ASSERT(min_size > IntVector::getZero(min_size.getDim()));
   TBOX_ASSERT(max_size >= min_size);
   TBOX_ASSERT(cut_factor > IntVector::getZero(min_size.getDim()));
   TBOX_ASSERT(bad_interval >= IntVector::getZero(min_size.getDim()));
   TBOX_ASSERT(physical_boxes.size() > 0);
   TBOX_ASSERT(!boxes.isOrdered());

   const tbox::Dimension& dim(max_size.getDim());

   BoxContainer in_boxes(boxes);
   boxes.clear();

   while (!in_boxes.isEmpty()) {

      Box box = in_boxes.front();
      in_boxes.popFront();

      BoxContainer tmp_boxes;

      tbox::Array<std::list<int> > cut_points(dim.getValue());
      bool chop_box = findBestCutPointsGivenMax(cut_points,
            box,
            max_size,
            min_size,
            cut_factor);

      if (chop_box) {
         TBOX_ASSERT(box.getBlockId().isValid());
         BoxContainer phys_block_boxes(physical_boxes, box.getBlockId());

         for (int id = 0; id < dim.getValue(); id++) {

            if (!cut_points[id].empty()) {

               tbox::Array<bool> bad_cut_points;

               findBadCutPointsForDirection(id,
                  bad_cut_points,
                  box,
                  phys_block_boxes,
                  bad_interval);
               fixBadCutPointsForDirection(id,
                  cut_points[id],
                  bad_cut_points,
                  box,
                  min_size(id),
                  cut_factor(id));

            }

         }

         chopBox(tmp_boxes,
            box,
            cut_points);

         boxes.spliceBack(tmp_boxes);

      } else {

         boxes.pushBack(box);

      }

   }

}

/*
 *************************************************************************
 *
 * Chop given box into a collection of boxes according to the collection
 * of cut points specified along each coordinate direction.   This box
 * list is formed from the resulting boxes.
 *
 *************************************************************************
 */

void
BoxUtilities::chopBox(
   BoxContainer& boxes,
   const Box& box,
   const tbox::Array<std::list<int> > cut_points)
{
   const tbox::Dimension& dim(box.getDim());

   TBOX_ASSERT(cut_points.getSize() == dim.getValue());

   if (!box.empty()) {

      boxes.clear();
      boxes.pushBack(box);

      BoxContainer tmp_boxes;
      for (int id = 0; id < dim.getValue(); id++) {

         tmp_boxes.clear();

         while (!boxes.isEmpty()) {

            Box chop_box = boxes.front();
            boxes.popFront();

            TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(dim, chop_box);

            if (!cut_points[id].empty()) {

               Index ilo = chop_box.lower();
               Index ihi = chop_box.upper();
               Index boxhi = chop_box.upper();

               const std::list<int>& cut_points_list = cut_points[id];
               std::list<int>::const_iterator cut = cut_points_list.begin();
#ifdef DEBUG_CHECK_ASSERTIONS
               int last_cut = tbox::MathUtilities<int>::getMin();
#endif
               while (cut != cut_points_list.end()) {
                  int cut_val = *cut;
#ifdef DEBUG_CHECK_ASSERTIONS
                  TBOX_ASSERT(last_cut <= cut_val);
                  last_cut = cut_val;
#endif
                  ihi(id) = cut_val - 1;
                  if ((ilo(id) < cut_val) && (ihi(id) <= boxhi(id))) {
                     Box new_box(ilo, ihi, box.getBlockId());
                     tmp_boxes.pushBack(new_box);
                     ilo(id) = cut_val;
                  }
                  cut++;
               }

               ihi(id) = chop_box.upper(id);
               Box last_box(ilo, ihi, box.getBlockId());
               tmp_boxes.pushBack(last_box);

            } else {
               tmp_boxes.pushBack(chop_box);
            }

         }

         boxes = tmp_boxes;

      }

   }
}

/*
 *************************************************************************
 *
 * Test each box in this box list for its intersection with the physical
 * domain boundary when it is grown by the given ghost width.  If the
 * ghost box lies entirely within the domain, or if all of its ghost
 * cells intersect the domain boundary appropriately, then the box will
 * not be changed.  Otherwise, the box is removed from the list and is
 * replaced by a new box formed by growing the original box to boundary.
 * This process eliminates domain boundary intersections which are
 * deemed unacceptable.  Intersections that are disallowed are those in
 * which a portion of the domain boundary is parallel to a box face and
 * lies strictly in the interior of the ghost cell mapped_box_level adjacent to
 * that face.  In other words, we eliminate ghost cell regions residing
 * outside of the domain and which are narrower than the ghost width.
 *
 *************************************************************************
 */

bool
BoxUtilities::extendBoxesToDomainBoundary(
   BoxContainer& boxes,
   const BoxContainer& domain,
   const IntVector& ext_ghosts)
{
   TBOX_ASSERT(!domain.isEmpty());
   TBOX_ASSERT(ext_ghosts >= IntVector::getZero(ext_ghosts.getDim()));

   bool out_val = false;

   BoxContainer out_boxes;

   while (!boxes.isEmpty()) {

      Box try_box = boxes.front();
      boxes.popFront();

      out_val = extendBoxToDomainBoundary(try_box, domain, ext_ghosts) ||
         out_val;

      out_boxes.pushBack(try_box);

   }

   boxes = out_boxes;

   return out_val;
}

bool
BoxUtilities::extendBoxToDomainBoundary(
   Box& box,
   const BoxContainer& domain,
   const IntVector& ext_ghosts)
{

   TBOX_ASSERT(!domain.isEmpty());
   TBOX_ASSERT(ext_ghosts >= IntVector::getZero(ext_ghosts.getDim()));

   const tbox::Dimension& dim(box.getDim());

   int id;
   bool out_val = false;

   if (!box.empty()) {

      Box test_ghost_box = box;
      test_ghost_box.grow(ext_ghosts);

      BoxContainer outside_domain(test_ghost_box);
      outside_domain.removeIntersections(domain);

      if (!outside_domain.isEmpty()) {

         for (id = 0; id < dim.getValue(); id++) {
            BoxContainer outside_boxes;

            // Test whether lower end of ghost box extends outside domain
            Box test_region = test_ghost_box;
            test_region.upper(id) = box.lower(id) - 1;

            outside_boxes = outside_domain;
            outside_boxes.intersectBoxes(test_region);

            int box_lo = box.lower(id);
            BoxContainer::iterator lb = outside_boxes.begin();
            for ( ; lb != outside_boxes.end(); ++lb) {
               box_lo = tbox::MathUtilities<int>::Min(box_lo, lb->upper(
                        id) + 1);
            }

            // Test whether upper end of ghost box extends outside domain
            test_region = test_ghost_box;
            test_region.lower(id) = box.upper(id) + 1;

            outside_boxes = outside_domain;
            outside_boxes.intersectBoxes(test_region);

            int box_hi = box.upper(id);
            for (lb = outside_boxes.begin(); lb != outside_boxes.end(); ++lb) {
               box_hi = tbox::MathUtilities<int>::Max(box_hi, lb->lower(
                        id) - 1);
            }

            if (!out_val) {
               out_val = ((box.lower(id) != box_lo) ||
                          (box.upper(id) != box_hi));
            }

            // Adjust box dimensions as necessary
            box.lower(id) = box_lo;
            box.upper(id) = box_hi;

         }

      }

   }

   return out_val;
}

/*
 *************************************************************************
 *
 * Grow each box in the list that is smaller than the specified minimum
 * size.  Each box that is grown must remain within the union of the
 * boxes of the given domain.  If the specified domain is an empty box
 * list, then each box will be grown to be as large as the minimum size
 * with no particular restrictions applied.  Note that this operation
 * may produce overlap regions among boxes on the list in either case.
 *
 *************************************************************************
 */

void
BoxUtilities::growBoxesWithinDomain(
   BoxContainer& boxes,
   const BoxContainer& domain,
   const IntVector& min_size)
{
   const tbox::Dimension& dim(min_size.getDim());

   int id;

   TBOX_ASSERT(min_size > IntVector::getZero(dim));

   if (!boxes.isEmpty()) {

      BoxContainer out_boxes;

      BoxContainer outside_domain;
      if (domain.isEmpty()) {
         Box big_box(boxes.getBoundingBox());
         big_box.grow(min_size);
         outside_domain.pushBack(big_box);
         outside_domain.grow(IntVector::getOne(dim));
         outside_domain.removeIntersections(big_box);
      } else {
         outside_domain = domain;
         outside_domain.unorder();
         outside_domain.grow(IntVector::getOne(dim));
         outside_domain.removeIntersections(domain);
      }

      while (!boxes.isEmpty()) {

         Box try_box = boxes.front();
         boxes.popFront();

         for (id = 0; id < dim.getValue(); id++) {

            int grow = min_size(id) - try_box.numberCells(id);

            if (grow > 0) {

               BoxContainer outside_boxes;
               Box test_region(dim);

               // How far may box be grown within domain in lower direction?
               test_region = try_box;
               test_region.lower(id) -= grow;
               test_region.upper(id) = try_box.lower(id) - 1;

               outside_boxes = outside_domain;
               outside_boxes.intersectBoxes(test_region);

               int grow_lo = try_box.lower(id) - grow;
               BoxContainer::iterator lb = outside_boxes.begin();
               for ( ; lb != outside_boxes.end(); ++lb) {
                  grow_lo =
                     tbox::MathUtilities<int>::Max(grow_lo, lb->upper(id) + 1);
               }

               // How far may box be grown within domain in upper direction?
               test_region = try_box;
               test_region.upper(id) += grow;
               test_region.lower(id) = try_box.upper(id) + 1;

               outside_boxes = outside_domain;
               outside_boxes.intersectBoxes(test_region);

               int grow_up = try_box.upper(id) + grow;
               for (lb = outside_boxes.begin(); lb != outside_boxes.end(); ++lb) {
                  grow_up =
                     tbox::MathUtilities<int>::Min(grow_up, lb->lower(id) - 1);
               }

               // Adjust box dimensions as necessary
               if ((grow_up - grow_lo + 1) < min_size(id)) {
                  try_box.lower(id) = grow_lo;
                  try_box.upper(id) = grow_up;
               } else {
                  int left = try_box.lower(id) - grow_lo;
                  int right = grow_up - try_box.upper(id);
                  int grow_half = grow / 2;

                  if (left < right) {
                     try_box.lower(id) -= ((left < grow_half) ? left
                                           : grow_half);
                     try_box.upper(id) = try_box.lower(id) + min_size(id) - 1;
                  } else {
                     try_box.upper(id) += ((right < grow_half) ? right
                                           : grow_half);
                     try_box.lower(id) = try_box.upper(id) - min_size(id) + 1;
                  }
               }

            }

         }

         out_boxes.pushBack(try_box);

      }

      boxes = out_boxes;

   }
}

/*
 *************************************************************************
 *
 * Grow each box in the list that is smaller than the specified minimum
 * size.  Each box that is grown must remain within the union of the
 * boxes of the given domain.  The domain is specified by the complement
 * of the local portion of the domain.
 *
 *************************************************************************
 */

void
BoxUtilities::growBoxWithinDomain(
   Box& box,
   const BoxContainer& local_domain_complement,
   const IntVector& min_size)
{
   const tbox::Dimension& dim(min_size.getDim());
   int id;

   TBOX_ASSERT(min_size > IntVector::getZero(dim));

   Box try_box = box;

   for (id = 0; id < dim.getValue(); id++) {

      int grow = min_size(id) - try_box.numberCells(id);

      if (grow > 0) {

         BoxContainer outside_boxes;
         Box test_region(dim);

         // How far may box be grown within domain in lower direction?
         test_region = try_box;
         test_region.lower(id) -= grow;
         test_region.upper(id) = try_box.lower(id) - 1;

         outside_boxes = local_domain_complement;
         outside_boxes.unorder();
         outside_boxes.intersectBoxes(test_region);

         BoxContainer::iterator lb = outside_boxes.begin(); 
         int grow_lo = try_box.lower(id) - grow;
         for ( ; lb != outside_boxes.end(); ++lb) {
            grow_lo =
               tbox::MathUtilities<int>::Max(grow_lo, lb->upper(id) + 1);
         }

         // How far may box be grown within domain in upper direction?
         test_region = try_box;
         test_region.upper(id) += grow;
         test_region.lower(id) = try_box.upper(id) + 1;

         outside_boxes = local_domain_complement;
         outside_boxes.unorder();
         outside_boxes.intersectBoxes(test_region);

         int grow_up = try_box.upper(id) + grow;
         for (lb = outside_boxes.begin(); lb != outside_boxes.end(); ++lb) {
            grow_up =
               tbox::MathUtilities<int>::Min(grow_up, lb->lower(id) - 1);
         }

         // Adjust box dimensions as necessary
         if ((grow_up - grow_lo + 1) < min_size(id)) {
            try_box.lower(id) = grow_lo;
            try_box.upper(id) = grow_up;
         } else {
            int left = try_box.lower(id) - grow_lo;
            int right = grow_up - try_box.upper(id);
            int grow_half = grow / 2;

            if (left < right) {
               try_box.lower(id) -= ((left < grow_half) ? left
                                     : grow_half);
               try_box.upper(id) = try_box.lower(id) + min_size(id) - 1;
            } else {
               try_box.upper(id) += ((right < grow_half) ? right
                                     : grow_half);
               try_box.lower(id) = try_box.upper(id) - min_size(id) + 1;
            }
         }

      }

   }

   box = try_box;
}

/*
 *************************************************************************
 *
 * Determine whether this box can be chopped according to specified
 * max, min, and factor constraints.  If the box may be chopped along
 * any face, true is returned.  Otherwise, false is returned.  For those
 * directions along which the box may be chopped, the cut points are
 * computed.  The procedure is as follows:
 *
 *    (1) Determine which directions chopping is allowed.
 *    (2) For each direction to chop, determine list of cut points.
 *
 * Important note: By convention, each integer cut point that is
 * computed corresponds to the cell index to the right of cut point.
 *
 *************************************************************************
 */

bool
BoxUtilities::findBestCutPointsGivenMax(
   tbox::Array<std::list<int> >& cut_points,
   const Box& box,
   const IntVector& max_size,
   const IntVector& min_size,
   const IntVector& cut_factor)
{
   const tbox::Dimension& dim(max_size.getDim());

   TBOX_DIM_ASSERT_CHECK_ARGS3(max_size, min_size, cut_factor);

   TBOX_ASSERT(min_size > IntVector::getZero(dim));
   TBOX_ASSERT(min_size <= max_size);
   TBOX_ASSERT(cut_factor > IntVector::getZero(dim));

   int id;
   bool chop_ok = false;

   cut_points.resizeArray(dim.getValue());

   for (id = 0; id < dim.getValue(); id++) {
      if (findBestCutPointsForDirectionGivenMax(id,
             cut_points[id],
             box,
             max_size(id),
             min_size(id),
             cut_factor(id))) {
         chop_ok = true;
      }
   }

   return chop_ok;

}

/*
 *************************************************************************
 *
 * Determine whether this box can be chopped according to specified
 * max, min, and factor constraints along given coordinate direction.
 * If the box may be chopped, true is returned; otherwise, false is
 * returned.  The procedure for determining the cuts is as follows:
 *
 *    (1) Adjust min and max values so that they are integer
 *        multiples of the cut factor.
 *    (2) Determine number of boxes, min and max box widths.
 *    (3) Determine list of cut points.
 *
 * Important note: By convention, each integer cut point that is
 * computed corresponds to the cell index to the right of cut point.
 *
 *************************************************************************
 */

bool
BoxUtilities::findBestCutPointsForDirectionGivenMax(
   const int idir,
   std::list<int>& cut_points,
   const Box& box,
   const int max_size,
   const int min_size,
   const int cut_factor)
{
   TBOX_ASSERT(!box.empty());
   TBOX_ASSERT(min_size > 0);
   TBOX_ASSERT(max_size >= min_size);
   TBOX_ASSERT(cut_factor > 0);

   cut_points.clear();

   bool chop_ok = (((box.numberCells(idir) % cut_factor)
                    || (box.numberCells(idir) <= max_size)
                    || (box.numberCells(idir) < 2 * min_size))
                   ? false : true);

   if (chop_ok) {

      int min = min_size;
      int max = max_size;

      chop_ok = false;

      int len = box.numberCells(idir);

      if (min % cut_factor) min = (min / cut_factor + 1) * cut_factor;
      if (max % cut_factor) max = (max / cut_factor) * cut_factor;

      /* make sure that max >= min.  In the case that
       * max equals min, max is increased only if the len is
       * not divisible by max.  This choice ensures that we
       * choose cut points that satisfy the min constraint
       * but possibly at the expense of breaking the max constraint.
       */
      if ((max < min) || ((max == min) && ((len % max) != 0))) {
         max = tbox::MathUtilities<int>::Min(2 * min, len / 2);
      }

      int num_boxes = 1;
      int max_width = min;
      int num_wide_boxes = num_boxes;
      int min_width = min;

      num_boxes = (len - 1) / max + 1;
      int len_remaining = len - num_boxes * min;

      if (len_remaining > 0) {
         int len_mult = len_remaining / cut_factor;
         num_wide_boxes = len_mult % num_boxes;
         if (num_wide_boxes != 0) {
            max_width += (len_mult / num_boxes + 1) * cut_factor;
            min_width = max_width - cut_factor;
         } else {
            max_width += (len_mult / num_boxes) * cut_factor;
            num_wide_boxes = num_boxes;
            min_width = 0;
         }
      }

      if (num_boxes > 1) {
         int mark = box.lower(idir);
         int wide_count = 0;
         for (int ic = 0; ic < num_boxes - 1; ic++) {
            int width = ((wide_count < num_wide_boxes)
                         ? max_width : min_width);
            mark += width;
            cut_points.push_back(mark);
            wide_count++;
         }

         chop_ok = true;
      }

   }

   return chop_ok;

}

/*
 *************************************************************************
 *
 * Determine whether this box may be chopped according to requested
 * number of cuts along each side.  If the box may be chopped along any
 * coordinate direction, true is returned.  Otherwise, false is
 * returned.  For those directions along which the box may be chopped,
 * the cut points are computed.  The procedure is as follows:
 *
 *    (1) Determine for which directions shopping is allowed.
 *    (2) For each direction to chop, determine list of cut points.
 *
 * Important note: By convention, each integer cut point that is
 * computed corresponds to the cell index to the right of cut point.
 *
 *************************************************************************
 */

bool
BoxUtilities::findBestCutPointsGivenNumber(
   tbox::Array<std::list<int> >& cut_points,
   const Box& box,
   const IntVector& number_boxes,
   const IntVector& min_size,
   const IntVector& cut_factor)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(number_boxes, min_size, cut_factor);

   const tbox::Dimension& dim(number_boxes.getDim());

   TBOX_ASSERT(!box.empty());
   TBOX_ASSERT(min_size > IntVector::getZero(dim));
   TBOX_ASSERT(number_boxes > IntVector::getZero(dim));
   TBOX_ASSERT(cut_factor > IntVector::getZero(dim));

   int id;

   cut_points.resizeArray(dim.getValue());

   tbox::Array<bool> chop_dir(dim.getValue());
   for (id = 0; id < dim.getValue(); id++) {
      cut_points[id].clear();
      chop_dir[id] = (((number_boxes(id) <= 1)
                       || (box.numberCells(id) % cut_factor(id))
                       || (box.numberCells(id) < 2 * min_size(id))
                       || (box.numberCells(id) <
                           (number_boxes(id) * min_size(id))))
                      ? false : true);
   }

   bool chop_ok = false;

   for (id = 0; id < dim.getValue(); id++) {

      if (chop_dir[id]) {

         if (findBestCutPointsForDirectionGivenNumber(id,
                cut_points[id],
                box,
                number_boxes(id),
                min_size(id),
                cut_factor(id))) {
            chop_ok = true;
         }

      }

   }

   return chop_ok;

}

/*
 *************************************************************************
 *
 * Determine whether this box may be chopped according to requested
 * number of cuts along given direction.  If the box may be chopped,
 * true is returned; otherwise, false is returned.  The procedure for
 * determining the cuts is as follows:
 *
 *    (1) Adjust min value so that it is an integer multiple of
 *        the cut factor.
 *    (2) Determine number of boxes, min and max box widths.
 *    (3) Determine list of cut points.
 *
 * Important note: By convention, each integer cut point that is
 * computed corresponds to the cell index to the right of cut point.
 *
 *************************************************************************
 */

bool
BoxUtilities::findBestCutPointsForDirectionGivenNumber(
   const int idir,
   std::list<int>& cut_points,
   const Box& box,
   const int num_boxes,
   const int min_size,
   const int cut_factor)
{
   TBOX_ASSERT(min_size > 0);
   TBOX_ASSERT(num_boxes > 0);
   TBOX_ASSERT(cut_factor > 0);

   cut_points.clear();

   bool chop_ok = (((num_boxes <= 1)
                    || (box.numberCells(idir) % cut_factor)
                    || (box.numberCells(idir) < 2 * min_size)
                    || (box.numberCells(idir) < num_boxes * min_size))
                   ? false : true);

   if (chop_ok) {

      chop_ok = false;

      int len = box.numberCells(idir);
      int min = min_size;

      if (min % cut_factor) min = (min / cut_factor + 1) * cut_factor;

      int max_width = min;
      int num_wide_boxes = num_boxes;
      int min_width = min;

      int len_remaining = len - num_boxes * min;

      if (len_remaining > 0) {
         int len_mult = len_remaining / cut_factor;
         num_wide_boxes = len_mult % num_boxes;
         if (num_wide_boxes != 0) {
            max_width += (len_mult / num_boxes + 1) * cut_factor;
            min_width = max_width - cut_factor;
         } else {
            max_width += (len_mult / num_boxes) * cut_factor;
            num_wide_boxes = num_boxes;
            min_width = 0;
         }
      }

      if (num_boxes > 1) {
         int mark = box.lower(idir);
         int wide_count = 0;
         for (int ic = 0; ic < num_boxes - 1; ic++) {
            int width = ((wide_count < num_wide_boxes)
                         ? max_width : min_width);
            mark += width;
            cut_points.push_back(mark);
            wide_count++;
         }

         chop_ok = true;
      }

   }

   return chop_ok;

}

/*
 *************************************************************************
 *
 * Return true if the box may have bad cut points, potentially.
 * Otherwise, return false.  Information about which directions may
 * have bad cut points is returned in the integer vector.  An entry of
 * zero indicates that there are no bad cut points for the box along
 * that coordinate direction.  An entry of one indicates that there
 * may be a bad cut point along that direction.
 *
 *************************************************************************
 */

bool
BoxUtilities::checkBoxForBadCutPoints(
   IntVector& bad_cut_information,
   const Box& box,
   const BoxContainer& physical_boxes,
   const IntVector& bad_interval)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(bad_cut_information, box,
      bad_interval);

   const tbox::Dimension& dim(box.getDim());

   bool found_bad = false;

   int id;

   bad_cut_information = IntVector::getZero(dim);
   for (id = 0; id < dim.getValue(); id++) {
      if (checkBoxForBadCutPointsInDirection(id,
             box,
             physical_boxes,
             bad_interval)) {
         bad_cut_information(id) = 1;
         found_bad = true;
      }
   }

   return found_bad;
}

/*
 *************************************************************************
 *
 * Return true if the box may have bad cut points along the given
 * coordinate direction, potentially.  Otherwise, return false.
 *
 *************************************************************************
 */

bool
BoxUtilities::checkBoxForBadCutPointsInDirection(
   const int id,
   const Box& box,
   const BoxContainer& physical_boxes,
   const IntVector& bad_interval)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, bad_interval);

   const tbox::Dimension& dim(box.getDim());

   TBOX_ASSERT(!box.empty());
   TBOX_ASSERT(bad_interval >= IntVector::getZero(dim));

   bool found_bad = false;

   if (physical_boxes.size() > 0) {

      int bad = bad_interval(id);

      int id2 = 0;
      while ((id2 < dim.getValue()) && !found_bad) {
         if (id2 != id) {

            int blo = box.lower(id);
            int bhi = box.upper(id);

            /*
             * Test lower box face in direction id2.
             */

            Box border = box;
            border.grow(bad_interval);
            border.upper(id2) = box.lower(id2) - 1;

            BoxContainer border_boxes(border);
            border_boxes.removeIntersections(physical_boxes);
            border_boxes.simplify();

            BoxContainer::iterator bb = border_boxes.begin();
            while (!found_bad && bb != border_boxes.end()) {
               found_bad = ((bb->lower(id) > (blo - bad))
                            || (bb->upper(id) < (bhi + bad)));
               ++bb;
            }

            if (!found_bad) {

               /*
                * Test upper box face in direction id2.
                */

               border = box;
               border.grow(bad_interval);
               border.lower(id2) = box.upper(id2) + 1;

               border_boxes.clear();
               border_boxes.pushBack(border);
               border_boxes.removeIntersections(physical_boxes);
               border_boxes.simplify();

               bb = border_boxes.begin();
               while (!found_bad && bb != border_boxes.end()) {
                  found_bad = ((bb->lower(id) > (blo - bad))
                               || (bb->upper(id) < (bhi + bad)));
                  ++bb;
               }

            }

         }
         id2++;
      }

   }
   return found_bad;
}

/*
 *************************************************************************
 *
 * Determine bad cut points for box based on the specified physical
 * domain and bad interval.   The cut information is returned as an
 * array (size = dim) of arrays (size = number of cells along edge
 * of the box) of boolean values.  A value of false indicates a
 * good cut point, a true value indicates that the cut is bad.
 *
 * Important notes: By convention, each integer cut point that is
 * computed corresponds to the cell index to the right of cut point.
 *
 *************************************************************************
 */

void
BoxUtilities::findBadCutPoints(
   tbox::Array<tbox::Array<bool> >& bad_cuts,
   const Box& box,
   const BoxContainer& physical_boxes,
   const IntVector& bad_interval)
{
   const tbox::Dimension& dim(box.getDim());

   TBOX_ASSERT(!box.empty());
   TBOX_ASSERT(bad_cuts.getSize() == dim.getValue());

   for (int id = 0; id < dim.getValue(); id++) {
      findBadCutPointsForDirection(id,
         bad_cuts[id],
         box,
         physical_boxes,
         bad_interval);
   }
}

/*
 *************************************************************************
 *
 * Determine bad cut points for box for given coordinate direction
 * based on the specified physical domain and bad interval.  The cut
 * information is returned as an array of integer values (size = number
 * of cells along edge of box.  A value of zero (0) indicates a good
 * cut point, a non-zero value indicates that the cut is bad.  The
 * process works as follows:
 *
 *    (1) Initialize all cut points to zero (zero = good).
 *    (2) Determine bad cut points based on domain configuration.
 *
 * Important notes: By convention, each integer cut point that is
 * computed corresponds to the cell index to the right of cut point.
 *
 *************************************************************************
 */

void
BoxUtilities::findBadCutPointsForDirection(
   const int id,
   tbox::Array<bool>& bad_cuts,
   const Box& box,
   const BoxContainer& physical_boxes,
   const IntVector& bad_interval)
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(box, bad_interval);

   const tbox::Dimension& dim(box.getDim());

   TBOX_ASSERT(!box.empty());
   TBOX_ASSERT(bad_interval >= IntVector::getZero(dim));

   int ic;

   /*
    * Initialize all bad cut points to false; i.e., all are good.
    */
   const int ncells = box.numberCells(id);
   bad_cuts.resizeArray(ncells);
   for (ic = 0; ic < ncells; ic++) {
      bad_cuts[ic] = false;
   }

   if (physical_boxes.size() == 0) {
      return; // Avoid the code below, which may crash for zero boxes.
   }

   /*
    * Determine whether box intersects physical boundary in such a way
    * that a bad cut point may result when the box is grown by the bad
    * interval.  To determine bad cut points for direction i, the box
    * must be intersected against the domain exterior in each directions
    * j not equal to i.  First, we check the lower end of box in direction j.
    * Then, we check the upper end of the box in direction j.  In each case,
    * the bad cut points are generated by considering each region that
    * intersects the domain exterior.
    */

   Box level_bounding_box = physical_boxes.getBoundingBox(box.getBlockId());

   for (int id2 = 0; id2 < dim.getValue(); id2++) {

      if (((dim.getValue() == 1) && id2 == id) ||
          ((dim.getValue() != 1) && (id2 != id))) {
         /*
          * Test lower box face in direction id2.
          */

         Box border = box;
         border.grow(bad_interval);
         border.upper(id2) = box.lower(id2) - 1;

         /*
          * limit the width of the border box to the width of the
          * domain to ensure that bad cut points near the boundary
          * of the box are not missed.
          */
         border.upper(id) = level_bounding_box.upper(id);
         border.lower(id) = level_bounding_box.lower(id);

         BoxContainer border_boxes(border);

         if (dim.getValue() > 1) {
            /*
             * only remove the level interior if the dimensionality of
             * the problem is greater than 1.
             */
            border_boxes.removeIntersections(physical_boxes);
         }

         if (!border_boxes.isEmpty()) {
            border_boxes.simplify();

            for (BoxContainer::iterator bbox(border_boxes);
                 bbox != border_boxes.end(); ++bbox) {
               findBadCutPointsForBorderAndDirection(id,
                  bad_cuts,
                  box,
                  *bbox,
                  bad_interval(id));
            }
         }

         /*
          * Test upper box face in direction id2.
          */

         border = box;
         border.grow(bad_interval);
         border.lower(id2) = box.upper(id2) + 1;

         /*
          * limit the width of the border box to the width of the
          * domain to ensure that bad cut points near the boundary
          * of the box are not missed.
          */
         border.upper(id) = level_bounding_box.upper(id);
         border.lower(id) = level_bounding_box.lower(id);

         border_boxes.clear();
         border_boxes.pushBack(border);

         if (dim.getValue() > 1) {
            /*
             * only remove the level interior if the dimensionality of
             * the problem is greater than 1.
             */
            border_boxes.removeIntersections(physical_boxes);
         }

         if (!border_boxes.isEmpty()) {
            border_boxes.simplify();
            for (BoxContainer::iterator bbox(border_boxes);
                 bbox != border_boxes.end(); ++bbox) {
               findBadCutPointsForBorderAndDirection(id,
                  bad_cuts,
                  box,
                  *bbox,
                  bad_interval(id));
            }
         }

      }
   }
}

/*
 *************************************************************************
 *
 * Adjust cut points if they coincide with bad cut points.
 *
 *************************************************************************
 */

void
BoxUtilities::fixBadCutPoints(
   tbox::Array<std::list<int> >& cuts,
   const tbox::Array<tbox::Array<bool> >& bad_cuts,
   const Box& box,
   const IntVector& min_size,
   const IntVector& cut_factor)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(box, min_size, cut_factor);

   const tbox::Dimension& dim(box.getDim());

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(cuts.getSize() == dim.getValue());
   TBOX_ASSERT(bad_cuts.getSize() == dim.getValue());
   bool bad_cuts_ok = true;
   for (int id = 0; id < dim.getValue(); id++) {
      bad_cuts_ok = bad_cuts_ok &&
         (bad_cuts[id].getSize() == box.numberCells(id));
   }
   TBOX_ASSERT(bad_cuts_ok);
   TBOX_ASSERT(!box.empty());
   TBOX_ASSERT(min_size > IntVector::getZero(dim));
   TBOX_ASSERT(cut_factor > IntVector::getZero(dim));
#endif

   for (int id = 0; id < dim.getValue(); id++) {
      fixBadCutPointsForDirection(id,
         cuts[id],
         bad_cuts[id],
         box,
         min_size(id),
         cut_factor(id));
   }
}

/*
 *************************************************************************
 *
 * For specified coordinate direction, adjust cut points if they
 * coincide with bad cut points.  This routine processes cut points
 * from the beginning of the list and the end of the list simultaneously.
 * When a bad cut is found when processing from the list beginning,
 * a good cut point is searched for by moving toward the lower end of
 * the box.  The opposite holds when processing from list end.  The
 * In either case, a new cut point will be inserted in the list if one
 * is found.  Otherwise, there will be one less cut point along the box
 * side.  This routine may be made more robust in the future.
 *
 *************************************************************************
 */

void
BoxUtilities::fixBadCutPointsForDirection(
   const int id,
   std::list<int>& cuts,
   const tbox::Array<bool>& bad_cuts,
   const Box& box,
   const int min_in,
   const int fact)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   std::list<int>::iterator cut = cuts.begin();
   TBOX_ASSERT(bad_cuts.getSize() == box.numberCells(id));
   bool cuts_strictly_increase = true;
   if (cut != cuts.end()) {
      int prev = *cut;
      cut++;
      while (cut != cuts.end() && cuts_strictly_increase) {
         if (*cut <= prev) {
            cuts_strictly_increase = false;
         }
         prev = *cut;
         cut++;
      }
   }
   TBOX_ASSERT(cuts_strictly_increase);
   TBOX_ASSERT(!box.empty());
   TBOX_ASSERT(min_in > 0);
   TBOX_ASSERT(fact > 0);
   bool cuts_satisfy_factor = true;
   cut = cuts.begin();
   while (cut != cuts.end() && cuts_satisfy_factor) {
      if ((((*cut) - box.lower(id)) % fact) != 0) {
         cuts_satisfy_factor = false;
      }
      cut++;
   }
   TBOX_ASSERT(cuts_satisfy_factor);
#endif

   /*
    * Do a quick check to see whether there are any bad cut points for
    * the box in the specified coordinate direction.  If not, we are done.
    */
   bool bad_point_exists = false;
   const int ncells = box.numberCells(id);
   for (int ic = 0; ic < ncells; ic++) {
      if (bad_cuts[ic]) {
         bad_point_exists = true;
      }
   }

   if (bad_point_exists) {

      std::list<int>::iterator cutlo = cuts.begin();

      if (cutlo != cuts.end()) {

         int min = min_in;

         if (min % fact) {
            min = (min / fact + 1) * fact;
         }

         const int offset = box.lower(id);
         const int ilo = box.lower(id);
         const int ihi = box.upper(id) + 1;

         int foo = 0;
         std::list<int>::iterator cuthi = cuts.insert(cuts.end(), foo);
         cuthi--;
         cuts.pop_back();

         while (cutlo != cuts.end() && cuthi != cuts.end() &&
                (*cutlo <= *cuthi)) {

            int bad_cut_val, below, above, try_cut;

            if (cutlo == cuthi) {

               if (bad_cuts[*cutlo - offset]) {

                  bool found_good_cut = false;

                  bad_cut_val = *cutlo;
                  std::list<int>::iterator tmplo = cutlo;
                  std::list<int>::iterator tmphi = cutlo;
                  tmplo--;
                  tmphi++;
                  cuts.erase(cutlo);

                  below = (tmplo != cuts.end() ? *tmplo : ilo);

                  try_cut = bad_cut_val - fact;
                  while ((try_cut >= (below + min))
                         && bad_cuts[try_cut - offset]) {
                     try_cut -= fact;
                  }

                  if (try_cut >= (below + min)) {
                     found_good_cut = true;
                     if (tmplo != cuts.end()) {
                        std::list<int>::iterator tmp = tmplo;
                        tmp++;
                        cuts.insert(tmp, try_cut);
                        cutlo = tmplo;
                        cutlo++;
                     } else {
                        cuts.push_front(try_cut);
                        cutlo = cuts.begin();
                     }
                     cutlo++;
                  } else {
                     cutlo = tmphi;
                  }

                  if (!found_good_cut) {
                     above = (tmphi != cuts.end() ? *tmphi : ihi);

                     try_cut = bad_cut_val + fact;
                     while ((try_cut <= (above - min))
                            && bad_cuts[try_cut - offset]) {
                        try_cut += fact;
                     }

                     if (try_cut <= (above - min)) {
                        if (tmphi != cuts.end()) {
                           cuts.insert(tmphi, try_cut);
                           cuthi = tmphi;
                           cuthi--;
                        } else {
                           cuthi = cuts.insert(cuts.end(), try_cut);
                        }
                        cuthi--;
                     } else {
                        cuthi = tmplo;
                     }
                  }

               } else {
                  cutlo++;
                  cuthi--;
               }

            } else {

               if (bad_cuts[*cutlo - offset]) {

                  bad_cut_val = *cutlo;
                  std::list<int>::iterator tmplo = cutlo;
                  tmplo--;
                  cuts.erase(cutlo);

                  below = (tmplo != cuts.end() ? *tmplo : ilo);

                  try_cut = bad_cut_val - fact;
                  while ((try_cut >= (below + min))
                         && bad_cuts[try_cut - offset]) {
                     try_cut -= fact;
                  }

                  if (try_cut >= (below + min)) {
                     if (tmplo != cuts.end()) {
                        std::list<int>::iterator tmp = tmplo;
                        tmp++;
                        cuts.insert(tmplo, try_cut);
                        cutlo = tmplo;
                        cutlo++;
                     } else {
                        cuts.push_front(try_cut);
                        cutlo = cuts.begin();
                     }
                     cutlo++;
                  } else {
                     if (tmplo != cuts.end()) {
                        cutlo = tmplo;
                        cutlo++;
                     } else {
                        cutlo = cuts.begin();
                     }
                  }

               } else {
                  cutlo++;
               }

               if (bad_cuts[*cuthi - offset]) {

                  bad_cut_val = *cuthi;
                  std::list<int>::iterator tmphi = cuthi;
                  tmphi++;
                  cuts.erase(cuthi);

                  above = (tmphi != cuts.end() ? *tmphi : ihi);

                  try_cut = bad_cut_val + fact;
                  while ((try_cut <= (above - min))
                         && bad_cuts[try_cut - offset]) {
                     try_cut += fact;
                  }

                  if (try_cut <= (above - min)) {
                     if (tmphi != cuts.end()) {
                        cuts.insert(tmphi, try_cut);
                        cuthi = tmphi;
                        cuthi--;
                     } else {
                        cuthi = cuts.insert(cuts.end(), try_cut);
                     }
                     cuthi--;
                  } else {
                     if (tmphi != cuts.end()) {
                        cuthi = tmphi;
                        cuthi--;
                     } else {
                        cuthi = cuts.insert(cuts.end(), foo);
                        cuthi--;
                        cuts.pop_back();
                     }
                  }

               } else {
                  cuthi--;
               }

            }

         }

      }

   }
}

/*
 *************************************************************************
 *
 * Decompose each box in this box array into a list of non overlapping
 * boxes.  Moreover, the regions of index space formed by composing the
 * union of boxes on each box list are mutually disjoint.
 *
 *************************************************************************
 */

void
BoxUtilities::makeNonOverlappingBoxContainers(
   tbox::Array<BoxContainer>& box_list_array,
   const BoxContainer& boxes)
{
   const int nb = boxes.size();

   for (int i = 0; i < box_list_array.getSize(); i++) {
      box_list_array[i].clear();
   }

   box_list_array.resizeArray(nb);

   // Copy boxes into a list to preserve the original box array.
   BoxContainer box_list(boxes);

   // Remove portion of index space represented by array box from list.
   // Keep unique pieces on box list.
   BoxContainer::const_iterator itr(boxes);
   for (int ib = 0; ib < nb; ++ib, ++itr) {
      Box remove = *itr;

      for (BoxContainer::iterator l(box_list); l != box_list.end(); ++l) {
         Box intersection = remove * (*l);
         if (intersection.isSpatiallyEqual(*l)) {
            box_list_array[ib].pushBack(*l);
         }
      }
      box_list_array[ib].coalesce();

      box_list.removeIntersections(remove);
   }
}

}
}

#endif
