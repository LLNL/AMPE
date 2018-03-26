/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Utilities for working on DLBG edges.
 *
 ************************************************************************/
#ifndef included_hier_BoxLevelConnectorUtils_C
#define included_hier_BoxLevelConnectorUtils_C

#include "SAMRAI/hier/BoxLevelConnectorUtils.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/MappingConnectorAlgorithm.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/TimerManager.h"

#include <boost/make_shared.hpp>
#include <limits>
#include <cstdlib>
#include <list>

namespace SAMRAI {
namespace hier {

boost::shared_ptr<tbox::Timer> BoxLevelConnectorUtils::t_make_sorting_map;
boost::shared_ptr<tbox::Timer> BoxLevelConnectorUtils::t_compute_external_parts;
boost::shared_ptr<tbox::Timer> BoxLevelConnectorUtils::t_compute_external_parts_intersection;
boost::shared_ptr<tbox::Timer> BoxLevelConnectorUtils::t_compute_internal_parts;
boost::shared_ptr<tbox::Timer> BoxLevelConnectorUtils::t_compute_internal_parts_intersection;

tbox::StartupShutdownManager::Handler
BoxLevelConnectorUtils::s_initialize_finalize_handler(
   BoxLevelConnectorUtils::initializeCallback,
   0,
   0,
   BoxLevelConnectorUtils::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 ***********************************************************************
 ***********************************************************************
 */
BoxLevelConnectorUtils::BoxLevelConnectorUtils():
   d_sanity_check_precond(false),
   d_sanity_check_postcond(false)
{
}

/*
 ***********************************************************************
 * Given the base and head levels, determine whether the base nests
 * nests in the head.
 ***********************************************************************
 */
bool
BoxLevelConnectorUtils::baseNestsInHead(
   bool* locally_nests,
   const BoxLevel& base,
   const BoxLevel& head,
   const IntVector& base_swell,
   const IntVector& head_swell,
   const IntVector& head_nesting_margin,
   const BoxContainer* domain) const
{

   tbox::Dimension dim(head.getDim());

   TBOX_DIM_ASSERT_CHECK_ARGS2(head_nesting_margin, base_swell);

   TBOX_ASSERT(head.getMPI() == base.getMPI());

#ifdef DEBUG_CHECK_ASSERTIONS
   const IntVector& zero_vector(IntVector::getZero(dim));
   TBOX_ASSERT(base_swell >= zero_vector);
   TBOX_ASSERT(head_swell >= zero_vector);
   TBOX_ASSERT(head_nesting_margin >= zero_vector);
#endif

   IntVector required_gcw = base_swell;
   if (head.getRefinementRatio() <= base.getRefinementRatio()) {
      const IntVector ratio = base.getRefinementRatio()
         / head.getRefinementRatio();
      required_gcw += (head_swell + head_nesting_margin) * ratio;
   } else if (head.getRefinementRatio() >= base.getRefinementRatio()) {
      const IntVector ratio = head.getRefinementRatio()
         / base.getRefinementRatio();
      required_gcw += IntVector::ceilingDivide(
         (head_swell + head_nesting_margin),
         ratio);
   } else {
      TBOX_ERROR("BoxLevelConnectorUtils::baseNestsInHead: head index space\n"
         << "must be either a refinement or a coarsening of\n"
         << "base, but not both.");
   }

   Connector base_to_head(
      base,
      head,
      required_gcw);

   OverlapConnectorAlgorithm oca;
   oca.findOverlaps(base_to_head);

   bool rval = baseNestsInHead(
         locally_nests,
         base_to_head,
         base_swell,
         head_swell,
         head_nesting_margin,
         domain);

   return rval;
}

/*
 ***********************************************************************
 * Given a Connector, determine the extent to which the base nests in the
 * head.  The Connector is assumed, without verification, to be complete.
 *
 * This method returns true if the base, grown by (non-negative)
 * base_swell nests inside the head by a margin of head_nesting_margin.
 * base_swell should be in the base index space and head_nesting_margin
 * should be in the head index space.  The Connector GCW must be at least
 * the sum of the base_swell and the appropriately converted
 * head_nesting_margin.  base_swell and head_nesting_margin can be
 * interchangable if you do the index space conversion yourself.
 *
 * If the domain is given, disregard non-nesting parts that are outside
 * the domain.
 ***********************************************************************
 */
bool
BoxLevelConnectorUtils::baseNestsInHead(
   bool* locally_nests,
   const Connector& connector,
   const IntVector& base_swell,
   const IntVector& head_swell,
   const IntVector& head_nesting_margin,
   const BoxContainer* domain) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(
      connector.getBase(), base_swell, head_nesting_margin);
   TBOX_ASSERT(connector.isFinalized());
   const tbox::Dimension& dim(connector.getBase().getDim());
   TBOX_ASSERT(base_swell >= IntVector::getZero(dim));
   TBOX_ASSERT(head_nesting_margin >= IntVector::getZero(dim));

   /*
    * To ensure correct results, connector must be sufficiently wide.
    * It should be at least as wide as the combination of base_swell,
    * head_swell and head_nesting_margin.
    */
   const IntVector required_gcw =
      base_swell
      + (connector.getHeadCoarserFlag() ?
         head_swell * connector.getRatio() :
         IntVector::ceilingDivide(head_swell, connector.getRatio()))
      + (connector.getHeadCoarserFlag() ?
         head_nesting_margin * connector.getRatio() :
         IntVector::ceilingDivide(head_nesting_margin, connector.getRatio()))
   ;
   if (!(connector.getConnectorWidth() >= required_gcw)) {
      TBOX_ERROR("BoxLevelConnectorUtils::baseNestsInHead: connector lacks\n"
         << "sufficient ghost cell width for determining whether its base\n"
         << "nests inside its head.");
   }

   OverlapConnectorAlgorithm oca;

   const BoxLevel& base = connector.getBase();
   const BoxLevel& head = connector.getHead();
   const boost::shared_ptr<const BaseGridGeometry>& grid_geom(
      base.getGridGeometry());

   /*
    * We swell the base then check for the parts outside the head if
    * the head domain is grown by head_swell then shrunken by
    * head_nesting_margin.
    *
    * TODO: We can probably remove the base swelling step by converting
    * the base_swell into head index space and add it to
    * head_nesting_margin.
    */

   BoxLevel swelledbase(dim);
   if (base_swell == IntVector::getZero(dim)) {
      swelledbase = base;
   } else {
      const BoxContainer& base_mapped_boxes = base.getBoxes();
      swelledbase.initialize(
         base.getRefinementRatio(),
         grid_geom,
         base.getMPI());
      for (BoxContainer::const_iterator ni = base_mapped_boxes.begin();
           ni != base_mapped_boxes.end(); ++ni) {
         Box swelledbase_mapped_box(*ni);
         swelledbase_mapped_box.grow(base_swell);
         swelledbase.addBoxWithoutUpdate(swelledbase_mapped_box);
      }
      swelledbase.finalize();
   }

   BoxLevel swelledhead(dim);
   if (head_swell == IntVector::getZero(dim)) {
      swelledhead = head;
   } else {
      const BoxContainer& head_mapped_boxes = head.getBoxes();

      swelledhead.initialize(
         head.getRefinementRatio(),
         grid_geom,
         head.getMPI());

      for (BoxContainer::const_iterator ni = head_mapped_boxes.begin();
           ni != head_mapped_boxes.end(); ++ni) {
         Box swelledhead_mapped_box(*ni);

         swelledhead_mapped_box.grow(head_swell);
         swelledhead.addBoxWithoutUpdate(swelledhead_mapped_box);
      }
      swelledhead.finalize();
   }

   Connector swelledbase_to_swelledhead(connector);
   swelledbase_to_swelledhead.setBase(swelledbase);
   swelledbase_to_swelledhead.setHead(swelledhead);
   swelledbase_to_swelledhead.setWidth(
      connector.getConnectorWidth() - base_swell,
      true);

   swelledbase_to_swelledhead.growLocalNeighbors(head_swell);

   if (d_sanity_check_precond &&
       head_swell == IntVector::getZero(dim)) {
      /*
       * If head was swelled, it may generate undetected overlaps that
       * cannot be compensated for by shrinking the connector width.
       * The additional overlaps do not matter to the nesting check,
       * so it does not affect our result.  Nevertheless, because they
       * are not detected, don't make this check if head was swelled.
       */
      oca.assertOverlapCorrectness(swelledbase_to_swelledhead);
   }

   BoxLevel external(dim);
   Connector swelledbase_to_external;
   if (domain) {
      computeExternalParts(
         external,
         swelledbase_to_external,
         swelledbase_to_swelledhead,
         -head_nesting_margin,
         *domain);
   } else {
      computeExternalParts(
         external,
         swelledbase_to_external,
         swelledbase_to_swelledhead,
         -head_nesting_margin,
         BoxContainer());
   }
   if (domain) {
      /*
       * If domain is given, do not count external parts that are
       * outside the domain.  In many usages, part of base is outside
       * the domain and we want to ignore those parts.
       */
      MappingConnectorAlgorithm mca;
      BoxLevel domain_mapped_box_level(
         IntVector::getOne(dim),
         grid_geom,
         connector.getMPI(),
         BoxLevel::GLOBALIZED);
      for (BoxContainer::const_iterator bi = domain->begin();
           bi != domain->end(); ++bi) {
         domain_mapped_box_level.addBox(*bi);
      }
      Connector external_to_domain(
         external,
         domain_mapped_box_level,
         base_swell);
      oca.findOverlaps(external_to_domain);
      BoxLevel finalexternal(dim);
      Connector external_to_finalexternal;
      computeInternalParts(
         finalexternal,
         external_to_finalexternal,
         external_to_domain,
         IntVector::getZero(dim),
         *domain);
      mca.modify(swelledbase_to_external,
         external_to_finalexternal,
         &external,
         &finalexternal);
   }

   if (locally_nests) {
      *locally_nests = external.getLocalNumberOfBoxes() == 0;
   }
   bool globally_nests = external.getGlobalNumberOfBoxes() == 0;

   return globally_nests;
}

/*
 ***********************************************************************
 * Make a Connector object for changing the Box indices of a BoxLevel.
 *
 * If sequentialize_global_indices is true, the indices are changed
 * such that they become globally sequential, with processor n
 * starting where processor n-1 ended.  In order to determine what the
 * global indices should be, an allgather communication is used to
 * determine how many mapped_boxes each processor has.  This is a
 * utility function for resetting Box indices to correspond to
 * patch indices while we try to be backward compatible with non-DLBG
 * parts of SAMRAI.
 *
 * If sort_mapped_boxes_by_corner is true, the local Boxes are
 * sorted by their corner indices.  This helps to de-randomize
 * Boxes that may be randomly ordered by non-deterministic
 * algorithms.
 ***********************************************************************
 */
void
BoxLevelConnectorUtils::makeSortingMap(
   BoxLevel& sorted_mapped_box_level,
   Connector& output_map,
   const BoxLevel& unsorted_mapped_box_level,
   bool sort_mapped_boxes_by_corner,
   bool sequentialize_global_indices,
   LocalId initial_sequential_index) const
{
   const tbox::Dimension& dim(unsorted_mapped_box_level.getDim());

   if (!sort_mapped_boxes_by_corner && !sequentialize_global_indices) {
      // Make a blank map.
      sorted_mapped_box_level = unsorted_mapped_box_level;
      output_map.clearNeighborhoods();
      output_map.setConnectorType(Connector::MAPPING);
      output_map.setBase(unsorted_mapped_box_level);
      output_map.setHead(sorted_mapped_box_level);
      output_map.setWidth(IntVector::getZero(dim), true);
      return;
   }

   t_make_sorting_map->start();

   const BoxContainer& cur_mapped_boxes =
      unsorted_mapped_box_level.getBoxes();

   LocalId last_index = initial_sequential_index - 1;

   if (sequentialize_global_indices) {
      // Increase last_index by the box count of all lower MPI ranks.

      int box_count =
         static_cast<int>(unsorted_mapped_box_level.getLocalNumberOfBoxes());
      unsorted_mapped_box_level.getMPI().parallelPrefixSum(&box_count, 1, 0);
      box_count -= static_cast<int>(unsorted_mapped_box_level.getLocalNumberOfBoxes());

      last_index += box_count;
   }

   std::vector<Box> real_mapped_box_vector;
   std::vector<Box> periodic_image_mapped_box_vector;
   if (!cur_mapped_boxes.isEmpty()) {
      /*
       * Bypass qsort if we have no mapped_boxes (else there is a memory
       * warning).
       */
      cur_mapped_boxes.separatePeriodicImages(
         real_mapped_box_vector,
         periodic_image_mapped_box_vector);
      if (sort_mapped_boxes_by_corner) {
         qsort((void *)&real_mapped_box_vector[0],
            real_mapped_box_vector.size(),
            sizeof(Box),
            qsortBoxCompare);
      }
   }

   sorted_mapped_box_level.initialize(
      unsorted_mapped_box_level.getRefinementRatio(),
      unsorted_mapped_box_level.getGridGeometry(),
      unsorted_mapped_box_level.getMPI(),
      BoxLevel::DISTRIBUTED);
   output_map.clearNeighborhoods();
   output_map.setConnectorType(Connector::MAPPING);
   output_map.setBase(unsorted_mapped_box_level);
   output_map.setHead(sorted_mapped_box_level);
   output_map.setWidth(IntVector::getZero(dim), true);

   for (std::vector<Box>::const_iterator ni = real_mapped_box_vector.begin();
        ni != real_mapped_box_vector.end(); ++ni) {

      const Box& cur_mapped_box = *ni;
      const Box new_mapped_box(cur_mapped_box,
                               ++last_index,
                               cur_mapped_box.getOwnerRank(),
                               cur_mapped_box.getPeriodicId());
      sorted_mapped_box_level.addBoxWithoutUpdate(new_mapped_box);

      /*
       * Now, add cur_mapped_box's periodic images, but give them
       * cur_mapped_box's new LocalId.  In finding the image mapped_boxes, we
       * use the fact that a real mapped_box's image follows the real
       * mapped_box in a BoxContainer.
       */
      BoxContainer::const_iterator ini = cur_mapped_boxes.find(cur_mapped_box);
      TBOX_ASSERT(ini != cur_mapped_boxes.end());
      ++ini; // Skip the real mapped_box to look for its image mapped_boxes.
      while (ini != cur_mapped_boxes.end() &&
             ini->getGlobalId() == cur_mapped_box.getGlobalId()) {
         const Box& image_mapped_box = *ini;
         const Box new_image_mapped_box(image_mapped_box,
                                        new_mapped_box.getLocalId(),
                                        new_mapped_box.getOwnerRank(),
                                        image_mapped_box.
                                        getPeriodicId());
         TBOX_ASSERT(new_image_mapped_box.getBlockId() ==
                     cur_mapped_box.getBlockId());
         sorted_mapped_box_level.addBoxWithoutUpdate(new_image_mapped_box);
         ++ini;
      }

      /*
       * Edge for the mapping.  By convention, image mapped_boxes are
       * not explicitly mapped.  Also by convention, we don't create
       * edges unless there is a change.
       */
      if (cur_mapped_box.getLocalId() != new_mapped_box.getLocalId()) {
         output_map.insertLocalNeighbor(new_mapped_box,
            cur_mapped_box.getId());
      }
      sorted_mapped_box_level.finalize();
   }

   t_make_sorting_map->stop();
}

/*
 *************************************************************************
 * for use when sorting integers using the C-library qsort
 *************************************************************************
 */
int
BoxLevelConnectorUtils::qsortBoxCompare(
   const void* v,
   const void* w)
{
   const Box& mapped_box_v(*(const Box *)v);
   const Box& mapped_box_w(*(const Box *)w);

   if (mapped_box_v.getBlockId() > mapped_box_w.getBlockId()) return 1;

   if (mapped_box_v.getBlockId() < mapped_box_w.getBlockId()) return -1;

   const tbox::Dimension& dim(mapped_box_v.getDim());

   const IntVector& lowv = mapped_box_v.lower();
   const IntVector& loww = mapped_box_w.lower();
   for (int i = 0; i < dim.getValue(); ++i) {
      if (lowv[i] > loww[i]) return 1;

      if (lowv[i] < loww[i]) return -1;
   }

   const IntVector& upv = mapped_box_v.upper();
   const IntVector& upw = mapped_box_w.upper();
   for (int i = 0; i < dim.getValue(); ++i) {
      if (upv[i] > upw[i]) return 1;

      if (upv[i] < upw[i]) return -1;
   }

   return 0;
}

/*
 *************************************************************************
 * Methods computeInternalParts and
 * computeExternalParts delegates to this method.
 *
 * Compare an input BoxLevel to a "reference" BoxLevel.
 * Identify parts of the input that are internal or external (depending
 * on the value of internal_or_external) to the reference
 * BoxLevel, and store the in/external parts in a BoxLevel.
 * Create the input_to_parts Connector between the input and these
 * parts.
 *
 * For generality, the reference BoxLevel can be grown a
 * specified amount (nesting_width) before comparing.  nesting_width
 * must be in the index space of the input BoxLevel (not the
 * reference BoxLevel, despite the name).  A negative growth
 * indicates shrinking the reference layer at its boundary.
 *
 * As a practical consideration of how this method is used, we do not
 * shrink the reference layer where it touches the domain boundary.
 * This feature can be disabled by specifying an uninitialized domain
 * object.
 *
 * On return, input_to_parts is set to an appropriate mapping for use
 * in MappingConnectorAlgorithm::modify().
 *
 * This method does not require any communication.
 *
 * Formula for computing external parts:
 *
 * Definitions:
 * L = input BoxLevel
 * R = reference BoxLevel
 * g = nesting width (non-negative or non-positive, but not mixed)
 * E = parts of L external to R^g (R^g is R grown by g)
 * I = parts of L internal to R^g (R^g is R grown by g)
 * O = domain (without periodic images).  Universe, if not specified.
 * \ = set theory notation.  x \ y means set x with y removed from it.
 *
 * For non-negative g:
 *
 * E := L \ { (R^g) <intersection> O }
 * I := L <intersection> { (R^g) <intersection> O }
 *
 * For non-positive g:
 *
 * E := L <intersection> { ( ( (R^1) \ R ) <intersection> O )^(-g) }
 * I := L \ { ( ( (R^1) \ R ) <intersection> O )^(-g) }
 *
 * A requirement of the computation for negative g is that input must
 * nest in R^(1-g).  In other words: L \ (R^(1-g)} = <empty>.  If not
 * satisfied, this method may classify some external parts as
 * internal.
 *
 *************************************************************************
 */
void
BoxLevelConnectorUtils::computeInternalOrExternalParts(
   BoxLevel& parts,
   Connector& input_to_parts,
   char internal_or_external,
   const Connector& input_to_reference,
   const IntVector& nesting_width,
   const BoxContainer& domain) const
{
   const BoxLevel& input = input_to_reference.getBase();

   const boost::shared_ptr<const BaseGridGeometry>& grid_geometry(
      input.getGridGeometry());

   const tbox::Dimension& dim(input.getDim());
   const IntVector& zero_vec(IntVector::getZero(input.getDim()));
   const IntVector& one_vec(IntVector::getOne(dim));

   const bool nonnegative_nesting_width = nesting_width >= zero_vec;

   const char* caller = internal_or_external == 'i' ?
      "computInternalParts" : "computeExternalparts";

   // Sanity check inputs.

   if (!(nesting_width >= zero_vec) && !(nesting_width <= zero_vec)) {
      TBOX_ERROR(
         "BoxLevelConnectorUtils::computeInternalOrExternalParts:" << caller
                                                                   <<
         ": error:\n"
                                                                   <<
         "nesting_width may not have mix of positive\n"
                                                                   <<
         "and negative values.");
   }

   if ( nesting_width != zero_vec &&
           input_to_reference.getConnectorWidth() < one_vec ) {
      TBOX_ERROR(
         "BoxLevelConnectorUtils::computeInternalOrExternalParts:" << caller
                                                                   <<
         ": error:\n"
         << "If nesting width (" << nesting_width << " is non-zero,\n"
         << "width of input_to_reference, " << input_to_reference.getConnectorWidth() << ",\n"
         << "must be at least 1.  Otherwise, correct results cannot be guaranteed."
         );
   }

   if ( !( input_to_reference.getConnectorWidth() >=
          (nonnegative_nesting_width ? nesting_width : -nesting_width)) ) {
      TBOX_ERROR(
         "BoxLevelConnectorUtils::computeInternalOrExternalParts:"
         << caller << ": error:\n"
         << "input_to_reference width, " << input_to_reference.getConnectorWidth()
         << ",\nmust be greater than the absolute value of nesting_width, "
         << nesting_width << ",\nto avoid erroneous results.");
   }

   parts.initialize(input.getRefinementRatio(),
      input.getGridGeometry(), input.getMPI());

   /*
    * Get the set of neighboring boxes on the reference BoxLevel.  We first
    * store these boxes in a NeighborSet in order to remove duplicate entries.
    * Then we move them into BoxContainer for each block for box manipulation.
    */
   BoxContainer reference_box_list;
   reference_box_list.order();
   input_to_reference.getLocalNeighbors(reference_box_list);

   /*
    * Bring reference_box_list into refinement ratio of input
    * (for intersection checks).
    */
   if (input_to_reference.getRatio() != IntVector::getOne(dim)) {
      if (input_to_reference.getHeadCoarserFlag()) {
         reference_box_list.refine(input_to_reference.getRatio());
      } else {
         reference_box_list.coarsen(input_to_reference.getRatio());
      }
   }

   /*
    * Build a search tree containing either the internal or external
    * parts of the reference BoxLevel, depending on the sign of
    * nesting_width.  If the nesting_width is non-negative, the
    * internal parts are the same as the reference, possibly after
    * growing.  If it is negative, shrinking the reference boxes does
    * not work.  We take its complement and grow the complement by
    * -nesting_width.  The result represents the external parts of the
    * reference.
    */

   const bool search_tree_represents_internal = nonnegative_nesting_width;

   if (search_tree_represents_internal) {

      if (!(nesting_width == zero_vec)) {
         reference_box_list.grow(nesting_width);
      }

   } else {

      /*
       * nesting_width is non-positive.  The external parts are given
       * by the grown boundary boxes.
       */

      computeBoxesAroundBoundary(
         reference_box_list,
         input.getRefinementRatio(),
         grid_geometry);
      // ... reference_boundary is now ( (R^1) \ R )

      if (domain.size() > 0) {

         if (input.getRefinementRatio() == one_vec) {
            reference_box_list.intersectBoxes(
               input.getRefinementRatio(),
               domain);
         } else {
            BoxContainer refined_domain(domain);
            refined_domain.refine(input.getRefinementRatio());
            refined_domain.makeTree(grid_geometry.get());
            reference_box_list.intersectBoxes(input.getRefinementRatio(),
                                              refined_domain);
         }

      }
      // ... reference_boundary is now ( ( (R^1) \ R ) <intersection> O )

      reference_box_list.grow(-nesting_width);
      // ... reference_boundary is now ( ( (R^1) \ R ) <intersection> O )^(-g)
   } // search_tree_represents_internal == false

   reference_box_list.makeTree(grid_geometry.get());

   /*
    * Keep track of last index so we don't give parts an index
    * used by input.  This is required because if we allow parts
    * to select an index it is not using, it may select one that is
    * used by input, creating an invalid mapping that prevents
    * MappingConnectorAlgorithm::modify() from working correctly.
    */
   LocalId last_used_index = input.getLastLocalId();

   /*
    * The output mapping has zero width.  For a mapping, zero width
    * means that no Box is mapped to something outside its
    * extent.
    */
   input_to_parts.clearNeighborhoods();
   input_to_parts.setConnectorType(Connector::UNKNOWN);
   input_to_parts.setBase(input);
   input_to_parts.setHead(parts);
   input_to_parts.setWidth(zero_vec, true);

   const bool compute_overlaps =
      search_tree_represents_internal == (internal_or_external == 'i');

   /*
    * For each Box in input, compare it to the search tree to
    * compute its overlapping (or non-overlapping) parts.
    */

   const BoxContainer& input_mapped_boxes = input.getBoxes();

   for (RealBoxConstIterator ni(input_mapped_boxes.realBegin());
        ni != input_mapped_boxes.realEnd(); ++ni) {

      const Box& input_mapped_box = *ni;
      const BoxId& input_mapped_box_id = input_mapped_box.getId();

      if (!input_to_reference.hasNeighborSet(input_mapped_box_id)) {
         /*
          * Absence of a reference neighbor set in the overlap
          * Connector means the input Box does not overlap the
          * reference BoxLevel.
          */
         if (compute_overlaps) {
            /*
             * Trying to get the overlapping parts.  Create empty
             * neighbor list to indicate there are no such parts.
             */
            input_to_parts.makeEmptyLocalNeighborhood(input_mapped_box_id);
         } else {
            /*
             * Trying to get the non-overlapping parts.
             * Non-overlapping parts is the whole box.
             */
            parts.addBox(input_mapped_box);
         }

      } else {

         BoxContainer parts_list(input_mapped_box);
         /*
          * Compute parts of input_mapped_box either overlapping
          * or nor overlapping the reference_box_list.
          *
          * Note about intersections in singularity neighbor blocks:
          * Cells from multiple singularity neighbor blocks can
          * coincide when transformed into input_mapped_box's block.
          * There is no way to specify that a cell in input_mapped_box
          * intersects in some singularity block neighbors but not
          * others.  By comparing to singularity neighbor blocks, we
          * take the convention that intersection in one singularity
          * block neighbor is considered intersection in all at the
          * same singularity.  When compute_overlaps == true,
          * this can lead to overspecifying parts and
          * underspecifying external parts, and vice versa.
          */

         t_compute_internal_parts_intersection->start();
         if (compute_overlaps) {
            parts_list.intersectBoxes(
               input.getRefinementRatio(),
               reference_box_list,
               true /* Count singularity neighbors */);
         } else {
            parts_list.removeIntersections(
               input.getRefinementRatio(),
               reference_box_list,
               true /* Count singularity neighbors */);
         }
         t_compute_internal_parts_intersection->stop();

         /*
          * Make Boxes from parts_list and create
          * Connector from input.
          */
         parts_list.simplify();
         if (parts_list.size() == 1 &&
             parts_list.front().isSpatiallyEqual(input_mapped_box)) {

            /*
             * The entire input_mapped_box is the part we want.
             * The input_mapped_box should be mapped to itself.
             * We can create such a map, but a missing map
             * means the same thing, so we omit the map
             */
            parts.addBox(input_mapped_box);

         } else {

            Connector::NeighborhoodIterator base_box_itr =
               input_to_parts.makeEmptyLocalNeighborhood(input_mapped_box_id);
            for (BoxContainer::iterator bi(parts_list);
                 bi != parts_list.end(); ++bi) {
               const Box
               parts_mapped_box((*bi),
                                ++last_used_index,
                                input_mapped_box.getOwnerRank());
               TBOX_ASSERT(parts_mapped_box.getBlockId() ==
                           input_mapped_box.getBlockId());
               parts.addBox(parts_mapped_box);

               // Set connectivities between input and internal.
               input_to_parts.insertLocalNeighbor(
                  parts_mapped_box,
                  base_box_itr);
            }

         } // parts_list

      } // !input_to_reference.hasNeighborSet(ni->getId())

   } // Loop through input_mapped_boxes

#ifdef DEBUG_CHECK_ASSERTIONS
   if (parts.getBoxes().isEmpty()) {
      /*
       * If there are no parts, then all in input
       * should be mapped to empty neighbor containers according
       * to the definition of a map in MappingConnectorAlgorithm::modify().
       */
      int a = input_to_parts.getLocalNumberOfNeighborSets();
      int b = static_cast<int>(input.getLocalNumberOfBoxes());
      if (a != b) {
         tbox::perr << "BoxLevelConnectorUtils::" << caller
                    << ": library error:\n"
                    <<
         "There are no parts, so input BoxLevel should be completely mapped away.\n"
                    << "However, not all input Boxes have been mapped.\n"
                    << "input BoxLevel:\n" << input.format("", 2)
                    << "input_to_parts:\n";
         input_to_parts.writeNeighborhoodsToErrorStream("");
         TBOX_ERROR("Library error\n");
      }
      TBOX_ASSERT(a == input_to_parts.numLocalEmptyNeighborhoods());
   }
#endif

   TBOX_ASSERT(input_to_parts.isLocal());
}

/*
 *************************************************************************
 * Given a BoxContainer, compute its boundary as a set of boxes located
 * just outside it.
 *
 * Given a set of boxes R, the boundary is computed as (R^1)\R.
 * R^1 means grown boxes in R by a width of one.
 *************************************************************************
 */
void
BoxLevelConnectorUtils::computeBoxesAroundBoundary(
   BoxContainer& boundary,
   const IntVector& refinement_ratio,
   const boost::shared_ptr<const BaseGridGeometry>& grid_geometry,
   const bool simplify_boundary_boxes) const
{

   const tbox::Dimension& dim(grid_geometry->getDim());
   const IntVector& one_vec(IntVector::getOne(dim));

   BoxContainer reference_mapped_boxes_tree(boundary);
   reference_mapped_boxes_tree.makeTree(grid_geometry.get());

   std::map<BlockId, BoxContainer> single_block_reference;
   if (grid_geometry->getNumberOfBlockSingularities() > 0) {
      for (BoxContainer::const_iterator bi = boundary.begin();
           bi != boundary.end(); ++bi) {
         single_block_reference[bi->getBlockId()].pushBack(*bi); 
      }
   }

   // Boundary starts as R

   boundary.grow(one_vec);
   // ... boundary is now (R^1)

   /*
    * Remove R from R^1, leaving alone boundary boxes in singularity
    * neighbor blocks.  These are specially handled in the following
    * if-block.
    */
   boundary.unorder();
   boundary.removeIntersections(
      refinement_ratio,
      reference_mapped_boxes_tree,
      false /* excludes singularity neighbors */);
   // ... boundary is now ( (R^1) \ R )


   /*
    * Separate boundary into containers for individual blocks required
    * by the block_id loop.  At the end of each block_id loop, stuff
    * the results for block_id back into boundary.
    */
   std::map<BlockId, BoxContainer> boundary_by_blocks;
   for (BoxContainer::const_iterator bi = boundary.begin();
        bi != boundary.end(); ++bi) {
      boundary_by_blocks[bi->getBlockId()].pushBack(*bi);
   }
   boundary.clear();


   if (grid_geometry->getNumberOfBlockSingularities() > 0) {
      /*
       * The boundary obtained by the formula (R^1)\R can have
       * errors at multiblock singularities.  Fix it here.
       *
       * Boundaries with codimension > 1 (node boundaries in 2D; node
       * and edge boundaries in 3D) passing through a singularity
       * cannot be computed correctly by the formula.
       *
       * What we want to determine is whether the boundary touches the
       * singularity point.
       *
       * - Problem 1: If R touches the singularity in all blocks that
       * touch the singularity, then the boundary does not pass
       * through the singularity.  However, at a reduced connectivity,
       * the formula computes that the boundary does pass through the
       * singularity.
       *
       * - Problem 2: If R touches the singularity in some, but not
       * all, blocks that touch the singularity, then the boundary
       * passes through the singularity point.  However, the formula
       * (R^1)\R may compute that the boundary does not.
       *
       * In both cases, the problem is that the cells representing the
       * codimension > 1 boundaries are removed or left behind when
       * they should not be.  At reduced connectivity singularities,
       * the fix is to always remove them because they do not live in
       * the index space of any block.  At enhanced connectivity, we determine
       * through some box calculus portions of the boundary that are touched
       * by Boxes in all blocks.  The boundary cannot touch these portions
       * so we remove these parts of the boundary.
       */

      for (std::map<BlockId, BoxContainer>::iterator bi = boundary_by_blocks.begin();
           bi != boundary_by_blocks.end(); ++bi) {

         const BlockId& block_id(bi->first);

         /*
          * Compute a version of singularity boxes for reduced
          * connectivity by removing enhanced connectivity boxes from
          * the singularity box list.  We will remove reduced
          * connectivity singularity boxes from the boundary
          * description, because they do not live in a valid index
          * space.
          */
         BoxContainer reduced_connectivity_singularity_boxes(
            grid_geometry->getSingularityBoxContainer(block_id));
         const std::list<BaseGridGeometry::Neighbor>& neighbors(
            grid_geometry->getNeighbors(block_id));

         for (std::list<BaseGridGeometry::Neighbor>::const_iterator ni =
              neighbors.begin();
              ni != neighbors.end(); ni++) {
            const BaseGridGeometry::Neighbor& neighbor(*ni);
            if (neighbor.isSingularity()) {
               reduced_connectivity_singularity_boxes.removeIntersections(
                  neighbor.getTransformedDomain());
            }
         }

         if (!reduced_connectivity_singularity_boxes.isEmpty()) {
            if (refinement_ratio != one_vec) {
               reduced_connectivity_singularity_boxes.refine(refinement_ratio);
            }
            bi->second.removeIntersections(
               reduced_connectivity_singularity_boxes);
         }

         /*
          * Intersect singularity_boxes with Boxes from each
          * singularity neighbor.  What remains is where all
          * singularity neighbors have Boxes touching the
          * singularity.  The remains tell us where the boundary does
          * not touch the singularity, overriding what the (R^1)\R
          * formula says.
          */
         BoxContainer singularity_boxes(
            grid_geometry->getSingularityBoxContainer(block_id));
         if (refinement_ratio != one_vec) {
            singularity_boxes.refine(refinement_ratio);
         }

         for (std::list<BaseGridGeometry::Neighbor>::const_iterator ni =
              neighbors.begin();
              ni != neighbors.end(); ni++) {
            const BaseGridGeometry::Neighbor& neighbor(*ni);
            const BlockId neighbor_block_id(neighbor.getBlockId());
            if (neighbor.isSingularity() &&
                reference_mapped_boxes_tree.hasBoxInBlock(neighbor_block_id)) {

               grid_geometry->transformBoxContainer(singularity_boxes,
                  refinement_ratio,
                  neighbor_block_id,
                  block_id);

               if (!single_block_reference[neighbor_block_id].hasTree()) {
                  single_block_reference[neighbor_block_id].makeTree(
                     grid_geometry.get());
               }
               singularity_boxes.intersectBoxes(
                  single_block_reference[neighbor_block_id]);

               grid_geometry->transformBoxContainer(singularity_boxes,
                  refinement_ratio,
                  block_id,
                  neighbor_block_id);
            }
         }

         bi->second.removeIntersections(singularity_boxes);

      } // for std::map<BlockId, ...

   } // grid_geometry->getNumberOfBlockSingularities() > 0

   if (simplify_boundary_boxes) {
      for (std::map<BlockId, BoxContainer>::iterator mi = boundary_by_blocks.begin();
           mi != boundary_by_blocks.end(); ++mi) {
         mi->second.simplify();
      }
   }

   // Set correct box ids.
   for (std::map<BlockId, BoxContainer>::iterator bi = boundary_by_blocks.begin();
        bi != boundary_by_blocks.end(); ++bi) {
      BoxContainer &boxes(bi->second);
      for ( BoxContainer::iterator bj=boxes.begin(); bj!=boxes.end(); ++bj ) {
         bj->setId(BoxId( bj->getLocalId(), bj->getOwnerRank(),
                          bj->getPeriodicId() ));
      }
   }

   for (std::map<BlockId, BoxContainer>::iterator bi = boundary_by_blocks.begin();
        bi != boundary_by_blocks.end(); ++bi) {
      boundary.spliceBack(bi->second);
   }

}

/*
 *************************************************************************
 * Given a mapping from an original BoxLevel to parts to be
 * removed, construct the remainder BoxLevel and the mapping from
 * the original to a remainder.
 *************************************************************************
 */

void
BoxLevelConnectorUtils::makeRemainderMap(
   BoxLevel& remainder,
   Connector& orig_to_remainder,
   const Connector& orig_to_rejection) const
{
   TBOX_ASSERT(orig_to_rejection.isLocal());

   const tbox::Dimension& dim(remainder.getDim());

   /*
    * remainder_nodes starts as a copy of orig codes.
    * It will be modified to become the remainder version.
    *
    * orig_to_remainder is the mapping between orig and
    * its properly remainder version.
    */

   const BoxLevel& orig = orig_to_rejection.getBase();
   const BoxContainer& orig_nodes = orig.getBoxes();
   const int rank = orig.getMPI().getRank();

   remainder = orig;

   orig_to_remainder.clearNeighborhoods();
   orig_to_remainder.setConnectorType(Connector::MAPPING);
   orig_to_remainder.setBase(orig);
   orig_to_remainder.setHead(remainder);
   orig_to_remainder.setWidth(IntVector::getZero(dim), true);

   /*
    * Track last used index to ensure we use unique indices for new
    * nodes, so that MappingConnectorAlgorithm::modify() works
    * properly.
    */
   LocalId last_used_index = orig.getLastLocalId();

   for (BoxContainer::const_iterator ni = orig_nodes.begin();
        ni != orig_nodes.end(); ++ni) {

      const Box& orig_node = *ni;
      const BoxId mapped_box_id = orig_node.getId();

      if (!orig_to_rejection.hasNeighborSet(mapped_box_id)) {
         /*
          * By the definition of a mapping Connector, no mapping means
          * the entire orig_node is rejected.
          *
          * - Erase rejected node from remainder
          * - Build connectivities in orig_to_remainder (empty neighbor list).
          */
         remainder.eraseBoxWithoutUpdate(orig_node);

         TBOX_ASSERT(!orig_to_remainder.hasNeighborSet(mapped_box_id));

         orig_to_remainder.makeEmptyLocalNeighborhood(mapped_box_id);
      } else if (orig_to_rejection.numLocalNeighbors(mapped_box_id) == 0) {
         /*
          * By the definition of a mapping Connector, empty mapping
          * means entire orig_node remains.
          *
          * No orig<==>remainder mapping is required.
          */
      } else {
         /*
          * The orig_node is partially rejected.
          *
          * - Erase rejected node from remainder
          * - Remove rejected parts to obtain remainder parts
          * - Add remaining parts to remainder
          * - Build connectivities in orig_to_remainder
          */

         Connector::ConstNeighborhoodIterator ci =
            orig_to_rejection.findLocal(mapped_box_id);

         remainder.eraseBoxWithoutUpdate(orig_node);

         BoxContainer remaining_parts_list(orig_node);

         for (Connector::ConstNeighborIterator vi = orig_to_rejection.begin(ci);
              vi != orig_to_rejection.end(ci); ++vi) {
            remaining_parts_list.removeIntersections((*vi));
         }
         /*
          * Coalesce the remaining_parts_list, because it may have unneeded
          * cuts.  The coalesce algorithm is O(N^2) or O(N^3), but we expect
          * the length of remaining_parts_list to be very small.
          */
         if (remaining_parts_list.size() > 1) {
            remaining_parts_list.coalesce();
         }

         /*
          * Create neighborhood of mapped_box_id in orig_to_remainder even if
          * remaining_parts_list is empty because its existence defines the
          * required mapping from the orig node to a (possibly empty)
          * container of nesting parts.
          */
         Connector::NeighborhoodIterator base_box_itr =
            orig_to_remainder.makeEmptyLocalNeighborhood(mapped_box_id);
         for (BoxContainer::iterator bi(remaining_parts_list);
              bi != remaining_parts_list.end(); ++bi) {
            Box new_box = (*bi);
            Box new_node(new_box,
                         ++last_used_index,
                         rank);
            TBOX_ASSERT(new_node.getBlockId() == orig_node.getBlockId());
            remainder.addBoxWithoutUpdate(new_node);
            orig_to_remainder.insertLocalNeighbor(new_node, base_box_itr);
         }
      }

   }
   remainder.finalize();
}

/*
 *************************************************************************
 * Add periodic images to a BoxLevel.
 *
 * We add the periodic images by examining real Mapped_boxes in the
 * BoxLevel.  For each real mapped_box, consider all of its
 * possible periodic images and add those that are within the
 * given width of the domain.
 *************************************************************************
 */

void
BoxLevelConnectorUtils::addPeriodicImages(
   BoxLevel& mapped_box_level,
   const BoxContainer& domain_search_tree,
   const IntVector& threshold_distance) const
{
   const PeriodicShiftCatalog* shift_catalog =
      PeriodicShiftCatalog::getCatalog(mapped_box_level.getDim());

   if (!shift_catalog->isPeriodic()) {
      return; // No-op.
   }

   boost::shared_ptr<BoxContainer> domain_tree_for_mapped_box_level(
      boost::make_shared<BoxContainer>(domain_search_tree));
   domain_tree_for_mapped_box_level->refine(
      mapped_box_level.getRefinementRatio());
   domain_tree_for_mapped_box_level->makeTree(NULL);

   const BoxContainer& domain_tree =
      *domain_tree_for_mapped_box_level;

   const IntVector& mapped_box_level_growth = threshold_distance;

   const BoxContainer& level_boxes(mapped_box_level.getBoxes());
   for (RealBoxConstIterator ni(level_boxes.realBegin());
        ni != level_boxes.realEnd(); ++ni) {

      const Box& mapped_box = *ni;
      for (int s = 1; s < shift_catalog->getNumberOfShifts(); ++s) {
         PeriodicId id(s);
         const IntVector& try_shift =
            shift_catalog->shiftNumberToShiftDistance(id);
         Box box = mapped_box;
         box.shift(try_shift);
         box.grow(mapped_box_level_growth);
         if (domain_tree.hasOverlap(box)) {
            mapped_box_level.addPeriodicBox(mapped_box, id);
         }
      }
   }

}

/*
 *************************************************************************
 * Add periodic images to a BoxLevel, and update Connectors that
 * require new edges incident on the additions to the BoxLevel.
 *
 * We add the periodic images by examining real Mapped_boxes in the
 * BoxLevel.  For each real mapped_box, consider all of its
 * possible periodic images and add those that are within the
 * Connector width distance of the domain.  (We are not interested in
 * periodic images so far from the domain that they are never used.)
 *
 * After adding periodic images, we bridge through
 * mapped_box_level<==>anchor<==>anchor so bridge can find the periodic
 * edges.
 *************************************************************************
 */

void
BoxLevelConnectorUtils::addPeriodicImagesAndRelationships(
   BoxLevel& mapped_box_level,
   Connector& mapped_box_level_to_anchor,
   Connector& anchor_to_mapped_box_level,
   const BoxContainer& domain_search_tree,
   const Connector& anchor_to_anchor) const
{
   OverlapConnectorAlgorithm oca;

   if (d_sanity_check_precond) {
      if (!mapped_box_level_to_anchor.isTransposeOf(anchor_to_mapped_box_level)) {
         TBOX_ERROR(
            "BoxLevelConnectorUtils::addPeriodicImages: non-transposed connector inputs.\n"
            << "mapped_box_level_to_anchor and anchor_to_mapped_box_level\n"
            << "must be mutual transposes.");
      }
      mapped_box_level_to_anchor.assertTransposeCorrectness(
         anchor_to_mapped_box_level);
      if (oca.checkOverlapCorrectness(anchor_to_anchor)) {
         TBOX_ERROR(
            "BoxLevelConnectorUtils::addPeriodicImages: input anchor_to_anchor\n"
            << "Connector failed edge correctness check.");
      }
      if (oca.checkOverlapCorrectness(anchor_to_mapped_box_level, false, true,
             true)) {
         TBOX_ERROR(
            "BoxLevelConnectorUtils::addPeriodicImages: input anchor_to_mapped_box_level\n"
            << "Connector failed edge correctness check.");
      }
      if (oca.checkOverlapCorrectness(mapped_box_level_to_anchor, false, true,
             true)) {
         TBOX_ERROR(
            "BoxLevelConnectorUtils::addPeriodicImages: input mapped_box_level_to_anchor\n"
            << "Connector failed edge correctness check.");
      }
   }
   if (!(anchor_to_anchor.getConnectorWidth() >=
         anchor_to_mapped_box_level.getConnectorWidth())) {
      TBOX_ERROR("BoxLevelConnectorUtils::addPeriodicImages: anchor_to_anchor width\n"
         << anchor_to_anchor.getConnectorWidth() << " is insufficient for\n"
         << "generating periodic edges for anchor_to_mapped_box_level's width of "
         << anchor_to_mapped_box_level.getConnectorWidth() << ".\n");
   }

   const PeriodicShiftCatalog* shift_catalog =
      PeriodicShiftCatalog::getCatalog(anchor_to_anchor.getConnectorWidth(
            ).getDim());

   if (!shift_catalog->isPeriodic()) {
      return; // No-op.
   }

   const BoxLevel& anchor = anchor_to_mapped_box_level.getBase();

   BoxContainer domain_tree_for_mapped_box_level(domain_search_tree);
   domain_tree_for_mapped_box_level.refine(
      mapped_box_level.getRefinementRatio());
   domain_tree_for_mapped_box_level.makeTree(NULL);

   {
      /*
       * Add the periodic image mapped_boxes for mapped_box_level.
       *
       * Adding images to a Box without neighbors in the anchor
       * means that this method will not find any edges to the added
       * images.
       */
      mapped_box_level.clearForBoxChanges(false);
      if (0) {
         tbox::perr << "mapped_box_level:\n"
                    << mapped_box_level.format("BEFORE-> ", 3);
      }
      const BoxContainer& domain_tree =
         domain_tree_for_mapped_box_level;

      const IntVector& mapped_box_level_growth =
         mapped_box_level_to_anchor.getConnectorWidth();

      const BoxContainer& level_boxes(mapped_box_level.getBoxes());
      for (RealBoxConstIterator ni(level_boxes.realBegin());
           ni != level_boxes.realEnd(); ++ni) {

         const Box& mapped_box = *ni;
         Box grown_box = mapped_box;
         grown_box.grow(mapped_box_level_growth);
         bool images_added(false);
         for (int s = 1; s < shift_catalog->getNumberOfShifts(); ++s) {
            PeriodicId id(s);
            const IntVector& try_shift =
               shift_catalog->shiftNumberToShiftDistance(id);
            Box box = grown_box;
            box.shift(try_shift);
            if (domain_tree.hasOverlap(box)) {
               mapped_box_level.addPeriodicBox(mapped_box, id);
               images_added = true;
            }
         }
         if (d_sanity_check_precond) {
            if (images_added &&
                (!mapped_box_level_to_anchor.hasNeighborSet(mapped_box.getId()) ||
                 mapped_box_level_to_anchor.isEmptyNeighborhood(mapped_box.getId()))) {
               TBOX_WARNING(
                  "BoxLevelConnectorUtils::addPeriodicImages: Box " << mapped_box
                                                                    <<
                  "\nhas periodic images in or close to the domain\n"
                                                                    <<
                  "but it does not have any neighbors in the anchor BoxLevel.\n"
                                                                    <<
                  "This will lead to missing neighbors in the output.\n"
                                                                    <<
                  "If post-condition checking is enabled, this will\n"
                                                                    <<
                  "result in an error.\n");
            }
         }

      }
      if (0) {
         tbox::perr << "mapped_box_level:\n" << mapped_box_level.format("AFTER-> ", 3);
      }
   }

   if (0) {
      tbox::plog << "Before bridging for periodic edges:\n"
                 << "anchor_to_anchor:\n" << anchor_to_anchor.format("DBG-> ", 3)
                 << "anchor_to_anchor:\n" << anchor_to_anchor.format("DBG-> ", 3)
                 << "mapped_box_level_to_anchor:\n" << mapped_box_level_to_anchor.format("DBG-> ",
         3)
                 << "anchor_to_mapped_box_level:\n" << anchor_to_mapped_box_level.format("DBG-> ",
         3);
   }

   IntVector width_limit =
      anchor_to_mapped_box_level.getHeadCoarserFlag() ?
      mapped_box_level_to_anchor.getConnectorWidth() :
      anchor_to_mapped_box_level.getConnectorWidth();

   oca.setSanityCheckMethodPreconditions(d_sanity_check_precond);
   oca.bridge(mapped_box_level_to_anchor,
      anchor_to_anchor,
      anchor_to_anchor,
      anchor_to_mapped_box_level,
      width_limit);
   anchor_to_mapped_box_level.eraseEmptyNeighborSets();
   mapped_box_level_to_anchor.eraseEmptyNeighborSets();

   if (d_sanity_check_postcond) {
      // Expensive sanity check for consistency.
      size_t err1 = anchor_to_mapped_box_level.checkConsistencyWithBase();
      if (err1) {
         tbox::perr << "OverlapConnectorAlgorithm found " << err1
                    << " edge-base consistency errors in\n"
                    << "anchor_to_mapped_box_level after computing periodic images.\n";
      }
      size_t err2 = mapped_box_level_to_anchor.checkConsistencyWithBase();
      if (err2) {
         tbox::perr << "OverlapConnectorAlgorithm found " << err2
                    << " edge-base consistency errors in\n"
                    << "mapped_box_level_to_anchor after computing periodic images.\n";
      }
      size_t err3 = anchor_to_mapped_box_level.checkConsistencyWithHead();
      if (err3) {
         tbox::perr << "OverlapConnectorAlgorithm found " << err3
                    << " edge-mapped_box consistency errors in\n"
                    << "anchor_to_mapped_box_level after computing periodic images.\n";
      }
      size_t err4 = mapped_box_level_to_anchor.checkConsistencyWithHead();
      if (err4) {
         tbox::perr << "OverlapConnectorAlgorithm found " << err4
                    << " edge-mapped_box consistency errors in\n"
                    << "mapped_box_level_to_anchor after computing periodic images.\n";
      }
      if (err1 + err2 + err3 + err4) {
         TBOX_ERROR(
            "OverlapConnectorAlgorithm found consistency errors in\n"
            << "addPeriodicImages\n"
            << "anchor:\n" << anchor.format("ERR-> ", 3)
            << "mapped_box_level:\n" << mapped_box_level.format("ERR-> ", 3)
            << "anchor_to_anchor:\n" << anchor_to_anchor.format("ERR-> ", 3)
            << "anchor_to_mapped_box_level:\n"
            << anchor_to_mapped_box_level.format("ERR-> ", 3)
            << "mapped_box_level_to_anchor:\n"
            << mapped_box_level_to_anchor.format("ERR-> ", 3));
      }
   }
   if (d_sanity_check_postcond) {
      // Expensive sanity check for correctness.
      int err1 = oca.checkOverlapCorrectness(anchor_to_mapped_box_level);
      if (err1) {
         tbox::perr << "BoxLevelConnectorUtils::addPeriodicImages found " << err1
                    << " errors\n"
                    << "in anchor_to_mapped_box_level after\n"
                    << "computing periodic images.  If you enabled\n"
                    << "precondition checking, this is probably a\n"
                    << "library error.\n";
      }
      int err2 = oca.checkOverlapCorrectness(mapped_box_level_to_anchor);
      if (err2) {
         tbox::perr << "BoxLevelConnectorUtils::addPeriodicImages found " << err2
                    << " errors\n"
                    << "in mapped_box_level_to_anchor after\n"
                    << "computing periodic images.  If you enabled\n"
                    << "precondition checking, this is probably a\n"
                    << "library error.\n";
      }
      if (err1 + err2) {
         TBOX_ERROR(
            "BoxLevelConnectorUtils::addPeriodicImages found edge errors\n"
            << "in output data\n"
            << "anchor:\n" << anchor.format("ERR-> ", 3)
            << "mapped_box_level:\n" << mapped_box_level.format("ERR-> ", 3)
            << "anchor_to_anchor:\n" << anchor_to_anchor.format("ERR-> ", 3)
            << "anchor_to_mapped_box_level:\n" << anchor_to_mapped_box_level.format("ERR-> ", 3)
            << "mapped_box_level_to_anchor:\n" << mapped_box_level_to_anchor.format("ERR-> ", 3)
            << "anchor_to_mapped_box_level:\n"
            << anchor_to_mapped_box_level.format("ERR-> ", 3)
            << "mapped_box_level_to_anchor:\n"
            << mapped_box_level_to_anchor.format("ERR-> ", 3));
      }
   }
}

}
}
#endif
