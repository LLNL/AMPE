/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Algorithms for working with overlap Connectors.
 *
 ************************************************************************/
#ifndef included_hier_OverlapConnectorAlgorithm_C
#define included_hier_OverlapConnectorAlgorithm_C

#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxContainerUtils.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/hier/RealBoxConstIterator.h"
#include "SAMRAI/tbox/AsyncCommStage.h"
#include "SAMRAI/tbox/AsyncCommPeer.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/TimerManager.h"

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace hier {

boost::shared_ptr<tbox::Timer> OverlapConnectorAlgorithm::t_find_overlaps_rbbt;
boost::shared_ptr<tbox::Timer> OverlapConnectorAlgorithm::t_bridge;
boost::shared_ptr<tbox::Timer> OverlapConnectorAlgorithm::t_bridge_setup_comm;
boost::shared_ptr<tbox::Timer> OverlapConnectorAlgorithm::t_bridge_remove_and_cache;
boost::shared_ptr<tbox::Timer> OverlapConnectorAlgorithm::t_bridge_discover;
boost::shared_ptr<tbox::Timer> OverlapConnectorAlgorithm::
t_bridge_discover_get_neighbors;
boost::shared_ptr<tbox::Timer> OverlapConnectorAlgorithm::t_bridge_discover_form_rbbt;
boost::shared_ptr<tbox::Timer> OverlapConnectorAlgorithm::
t_bridge_discover_find_overlaps;
boost::shared_ptr<tbox::Timer> OverlapConnectorAlgorithm::t_bridge_share;
boost::shared_ptr<tbox::Timer> OverlapConnectorAlgorithm::t_bridge_receive_and_unpack;
boost::shared_ptr<tbox::Timer> OverlapConnectorAlgorithm::t_bridge_MPI_wait;

int OverlapConnectorAlgorithm::s_operation_mpi_tag = 0;
/*
 * Do we even need to use different tags each time we bridge???
 * Unique tags were used to help debug, but the methods may work
 * with reused tags anyway.
 */

tbox::SAMRAI_MPI OverlapConnectorAlgorithm::s_class_mpi(
   tbox::SAMRAI_MPI::commNull);

tbox::StartupShutdownManager::Handler
OverlapConnectorAlgorithm::s_initialize_finalize_handler(
   OverlapConnectorAlgorithm::initializeCallback,
   0,
   0,
   OverlapConnectorAlgorithm::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 ***********************************************************************
 ***********************************************************************
 */

OverlapConnectorAlgorithm::OverlapConnectorAlgorithm():
   d_sanity_check_method_preconditions(false),
   d_sanity_check_method_postconditions(false)
{
   /*
    * While we figure out how to use multiple communicators in SAMRAI,
    * we are still assuming that all communications use congruent
    * communicators.  This class just makes a duplicate communicator
    * to protect itself from unrelated communications in shared
    * communicators.
    */
   if (s_class_mpi.getCommunicator() == tbox::SAMRAI_MPI::commNull) {
      if (tbox::SAMRAI_MPI::usingMPI()) {
         const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
         s_class_mpi.dupCommunicator(mpi);
      }
   }

}

/*
 ***********************************************************************
 ***********************************************************************
 */

OverlapConnectorAlgorithm::~OverlapConnectorAlgorithm()
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::setSanityCheckMethodPreconditions(
   bool do_check)
{
   d_sanity_check_method_preconditions = do_check;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::setSanityCheckMethodPostconditions(
   bool do_check)
{
   d_sanity_check_method_postconditions = do_check;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::extractNeighbors(
   NeighborSet& neighbors,
   const Connector& connector,
   const BoxId& mapped_box_id,
   const IntVector& gcw) const
{
   const tbox::Dimension& dim(gcw.getDim());

#ifdef DEBUG_CHECK_ASSERTIONS
   if (!(gcw <= connector.getConnectorWidth())) {
      TBOX_ERROR("OverlapConnectorAlgorithm::extractNeighbors cannot provide neighbors for\n"
         << "a wider ghost cell width that used to initialize it.\n");
   }
   if (connector.getParallelState() != BoxLevel::GLOBALIZED &&
       mapped_box_id.getOwnerRank() != connector.getMPI().getRank()) {
      TBOX_ERROR("OverlapConnectorAlgorithm::extractNeighbors cannot get neighbor data\n"
         << "for a remote mapped_box unless in GLOBALIZED mode.\n");
   }
   if (!connector.getBase().hasBox(mapped_box_id)) {
      std::string dbgbord;
      TBOX_ERROR(
         "\nOverlapConnectorAlgorithm::extractNeighbors: mapped_box_id " << mapped_box_id
                                                         <<
         " is not in the base of the mapped_box_level.\n"
                                                         << "base:\n"
                                                         << connector.getBase().format(dbgbord, 2)
                                                         << "head:\n"
                                                         << connector.getHead().format(dbgbord, 2)
                                                         << "connector:\n"
                                                         << connector.format(dbgbord, 2));
   }
#endif

   /*
    * Temporarily disable extracting neighbors for remote boxes.  This
    * method functionality is not much used and prrobably should be
    * removed.
    */
   TBOX_ASSERT(mapped_box_id.getOwnerRank() == connector.getMPI().getRank());

   const boost::shared_ptr<const BaseGridGeometry>& grid_geom(
      connector.getBase().getGridGeometry());

   const Box& mapped_box(*connector.getBase().getBox(Box(dim, mapped_box_id)));
   Connector::ConstNeighborhoodIterator ins =
      connector.findLocal(mapped_box_id);
   neighbors.clear();
   if (ins != connector.end()) {
      if (gcw == connector.getConnectorWidth()) {
         for (Connector::ConstNeighborIterator ni = connector.begin(ins);
              ni != connector.end(ins); ++ni) {
            neighbors.insert(neighbors.end(), *ni);
         }
      }
      else {
         Box grown_mapped_box = mapped_box;
         grown_mapped_box.grow(gcw);
         if (connector.getHeadCoarserFlag() == false) {
            grown_mapped_box.refine(connector.getRatio());
         }
         for (Connector::ConstNeighborIterator ni = connector.begin(ins);
              ni != connector.end(ins); ++ni) {
            const Box& neighbor(*ni);
            Box nabr_box(neighbor);
            if (neighbor.getBlockId() != mapped_box.getBlockId()) {
               grid_geom->transformBox(nabr_box,
                  connector.getHead().getRefinementRatio(),
                  mapped_box.getBlockId(),
                  neighbor.getBlockId());
            }
            if (connector.getHeadCoarserFlag() == true) {
               nabr_box.refine(connector.getRatio());
            }
            if (grown_mapped_box.intersects(nabr_box)) {
               neighbors.insert(neighbors.end(), neighbor);
            }
         }
      }
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::extractNeighbors(
   Connector& other,
   const Connector& connector,
   const IntVector& gcw) const
{
   other.clearNeighborhoods();
   for (Connector::ConstNeighborhoodIterator ni = connector.begin();
        ni != connector.end(); ++ni) {

      const BoxId& mapped_box_id = *ni;
      Connector::NeighborhoodIterator base_box_itr =
         other.makeEmptyLocalNeighborhood(mapped_box_id);

      const tbox::Dimension& dim(gcw.getDim());

#ifdef DEBUG_CHECK_ASSERTIONS
      if (!(gcw <= connector.getConnectorWidth())) {
         TBOX_ERROR("OverlapConnectorAlgorithm::extractNeighbors cannot provide neighbors for\n"
            << "a wider ghost cell width that used to initialize it.\n");
      }
      if (connector.getParallelState() != BoxLevel::GLOBALIZED &&
          mapped_box_id.getOwnerRank() != connector.getMPI().getRank()) {
         TBOX_ERROR("OverlapConnectorAlgorithm::extractNeighbors cannot get neighbor data\n"
            << "for a remote mapped_box unless in GLOBALIZED mode.\n");
      }
      if (!connector.getBase().hasBox(mapped_box_id)) {
         std::string dbgbord;
         TBOX_ERROR(
            "\nOverlapConnectorAlgorithm::extractNeighbors: mapped_box_id " << mapped_box_id
                                                            <<
            " is not in the base of the mapped_box_level.\n"
                                                            << "base:\n"
                                                            << connector.getBase().format(dbgbord, 2)
                                                            << "head:\n"
                                                            << connector.getHead().format(dbgbord, 2)
                                                            << "connector:\n"
                                                            << connector.format(dbgbord, 2));
      }
#endif

      /*
       * Temporarily disable extracting neighbors for remote boxes.  This
       * method functionality is not much used and prrobably should be
       * removed.
       */
      TBOX_ASSERT(mapped_box_id.getOwnerRank() == connector.getMPI().getRank());

      const boost::shared_ptr<const BaseGridGeometry>& grid_geom (
         connector.getBase().getGridGeometry());

      const Box& mapped_box =
         *connector.getBase().getBox(Box(dim, mapped_box_id));

      if (gcw == connector.getConnectorWidth()) {
         for (Connector::ConstNeighborIterator si = connector.begin(ni);
              si != connector.end(ni); ++si) {
            other.insertLocalNeighbor(*si, base_box_itr);
         }
      } else {
         Box grown_mapped_box = mapped_box;
         grown_mapped_box.grow(gcw);
         if (connector.getHeadCoarserFlag() == false) {
            grown_mapped_box.refine(connector.getRatio());
         }
         for (Connector::ConstNeighborIterator si = connector.begin(ni);
              si != connector.end(ni); ++si) {
            const Box& neighbor = *si;
            Box nabr_box(neighbor);
            if (neighbor.getBlockId() != mapped_box.getBlockId()) {
               grid_geom->transformBox(nabr_box,
                  connector.getHead().getRefinementRatio(),
                  mapped_box.getBlockId(),
                  neighbor.getBlockId());
            }
            if (connector.getHeadCoarserFlag() == true) {
               nabr_box.refine(connector.getRatio());
            }
            if (grown_mapped_box.intersects(nabr_box)) {
               other.insertLocalNeighbor(neighbor, base_box_itr);
            }
         }
      }
   }
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::findOverlaps(
   Connector& connector,
   const bool ignore_self_overlap) const
{
   findOverlaps(connector,
      connector.getHead().getGlobalizedVersion(),
      ignore_self_overlap);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::findOverlaps(
   Connector& connector,
   const BoxLevel& globalized_head,
   const bool ignore_self_overlap) const
{
   findOverlaps_rbbt(connector, globalized_head, ignore_self_overlap);
}

/*
 ***********************************************************************
 * ignore_self_overlap should be set to true only if
 * - the base and head mapped_box_levels represent the same mapped_box_level.
 *   Two different mapped_box_level objects may represent the same
 *   mapped_box_level if they are of the same refinement ratio.
 * - you want to ignore overlaps between a mapped_box and itself.
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::findOverlaps_rbbt(
   Connector& connector,
   const BoxLevel& head,
   const bool ignore_self_overlap) const
{
   const tbox::Dimension dim(head.getDim());

   t_find_overlaps_rbbt->start();

   /*
    * Finding overlaps for this object, using
    * an externally provided head BoxLevel
    * meant to represent d_head.  We allow the
    * substitution of an external head because
    * we require the head is GLOBALIZED.  The
    * user may have a GLOBALIZED version already,
    * in which case we want to avoid the expense
    * of creating a temporary GLOBALIZED version.
    *
    * Global mapped_boxes provided by head are sorted in a BoxContainer
    * so they can be quickly searched to see which intersects the
    * boxes in this object.
    */
   if (head.getParallelState() != BoxLevel::GLOBALIZED) {
      TBOX_ERROR("OverlapConnectorAlgorithm::findOverlaps_rbbt() requires given head\n"
         << "to be GLOBALIZED.\n");
   }

   /*
    * The nomenclature "base" refers to the *this mapped_box_level
    * and "head" refer to the mapped_box_level in the argument.
    */
   const BoxLevel& base(connector.getBase());

   /*
    * Determine relationship between base and head index spaces.
    */
   const bool head_is_finer =
      head.getRefinementRatio() >= base.getRefinementRatio() &&
      head.getRefinementRatio() != base.getRefinementRatio();
   const bool base_is_finer =
      base.getRefinementRatio() >= head.getRefinementRatio() &&
      base.getRefinementRatio() != head.getRefinementRatio();

   /*
    * Create single container of visible head mapped_boxes
    * to generate the search tree.
    */
   const BoxContainer& rbbt = head.getGlobalBoxes();
   rbbt.makeTree(head.getGridGeometry().get());

   /*
    * A neighbor of a Box would be discarded if
    * - ignore_self_overlap is true,
    * - the two are equal by comparison, and
    * - they are from mapped_box_levels with the same refinement ratio
    *   (we cannot compare mapped_box mapped_box_level pointers because that
    *   does not work when a mapped_box mapped_box_level is a temporary globalized object)
    */
   const bool discard_self_overlap =
      ignore_self_overlap &&
      (connector.getBase().getRefinementRatio() == head.getRefinementRatio());

   /*
    * Discard current overlaps in connector (if any).
    */
   connector.clearNeighborhoods();

   /*
    * Use BoxTree to find local base Boxes intersecting head Boxes.
    */
   NeighborSet nabrs_for_box;
   const BoxContainer& base_mapped_boxes = base.getBoxes();
   for (RealBoxConstIterator ni(base_mapped_boxes.realBegin());
        ni != base_mapped_boxes.realEnd(); ++ni) {

      const Box& base_mapped_box = *ni;

      // Grow the base_mapped_box and put it in the head refinement ratio.
      Box box = base_mapped_box;
      box.grow(connector.getConnectorWidth());
      if (head_is_finer) {
         box.refine(connector.getRatio());
      } else if (base_is_finer) {
         box.coarsen(connector.getRatio());
      }

      // Add found overlaps to neighbor set for mapped_box.
      rbbt.findOverlapBoxes(nabrs_for_box,
         box,
                            // base_mapped_box.getBlockId(),
         head.getRefinementRatio(),
         true);
      if (discard_self_overlap) {
         nabrs_for_box.order();
         nabrs_for_box.erase(base_mapped_box);
      }
      if (!nabrs_for_box.isEmpty()) {
         connector.insertNeighbors(nabrs_for_box, base_mapped_box.getId());
         nabrs_for_box.clear();
      }

   }
   connector.setConnectorType(Connector::COMPLETE_OVERLAP);

   if (d_sanity_check_method_postconditions) {
      connector.assertConsistencyWithBase();
      connector.assertConsistencyWithHead();
   }

   t_find_overlaps_rbbt->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::bridge(
   Connector& west_to_east,
   Connector& east_to_west,
   const Connector& west_to_cent,
   const Connector& cent_to_east,
   const Connector& east_to_cent,
   const Connector& cent_to_west) const
{
   const IntVector& zero_vector(
      IntVector::getZero(cent_to_east.getConnectorWidth().getDim()));
   const IntVector connector_width_limit(
      cent_to_east.getConnectorWidth().getDim(), -1); // No user-imposed limit.
   privateBridge(
      west_to_east,
      &east_to_west == &west_to_east ? NULL : &east_to_west,
      west_to_cent,
      cent_to_east,
      east_to_cent,
      cent_to_west,
      false,
      zero_vector,
      false,
      zero_vector,
      connector_width_limit);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::bridgeWithNesting(
   Connector& west_to_east,
   Connector& east_to_west,
   const Connector& west_to_cent,
   const Connector& cent_to_east,
   const Connector& east_to_cent,
   const Connector& cent_to_west,
   const IntVector& cent_growth_to_nest_west,
   const IntVector& cent_growth_to_nest_east,
   const IntVector& connector_width_limit) const
{
   privateBridge(
      west_to_east,
      &east_to_west == &west_to_east ? NULL : &east_to_west,
      west_to_cent,
      cent_to_east,
      east_to_cent,
      cent_to_west,
      (cent_growth_to_nest_west(0) >= 0),
      cent_growth_to_nest_west,
      (cent_growth_to_nest_east(0) >= 0),
      cent_growth_to_nest_east,
      connector_width_limit);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::bridge(
   Connector& west_to_east,
   Connector& east_to_west,
   const Connector& west_to_cent,
   const Connector& cent_to_east,
   const Connector& east_to_cent,
   const Connector& cent_to_west,
   const IntVector& connector_width_limit) const
{
   const IntVector& zero_vector(
      IntVector::getZero(cent_to_east.getConnectorWidth().getDim()));
   privateBridge(
      west_to_east,
      &east_to_west == &west_to_east ? NULL : &east_to_west,
      west_to_cent,
      cent_to_east,
      east_to_cent,
      cent_to_west,
      false,
      zero_vector,
      false,
      zero_vector,
      connector_width_limit);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::bridge(
   Connector& west_to_cent,
   const Connector& cent_to_east,
   const Connector& east_to_cent,
   Connector& cent_to_west,
   const IntVector& connector_width_limit) const
{
   const IntVector& zero_vector(
      IntVector::getZero(cent_to_east.getConnectorWidth().getDim()));
   privateBridge(
      west_to_cent,
      cent_to_east,
      east_to_cent,
      cent_to_west,
      false,
      zero_vector,
      false,
      zero_vector,
      connector_width_limit);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::bridge(
   Connector& west_to_east,
   const Connector& west_to_cent,
   const Connector& cent_to_east,
   const Connector& east_to_cent,
   const Connector& cent_to_west,
   const IntVector& connector_width_limit) const
{
   const IntVector& zero_vector(
      IntVector::getZero(cent_to_east.getConnectorWidth().getDim()));
   privateBridge(
      west_to_east,
      NULL,
      west_to_cent,
      cent_to_east,
      east_to_cent,
      cent_to_west,
      false,
      zero_vector,
      false,
      zero_vector,
      connector_width_limit);
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::bridge(
   Connector& west_to_east,
   const Connector& west_to_cent,
   const Connector& cent_to_east,
   const Connector& east_to_cent,
   const Connector& cent_to_west) const
{
   const IntVector& zero_vector(
      IntVector::getZero(cent_to_east.getConnectorWidth().getDim()));
   const IntVector connector_width_limit(
      cent_to_east.getConnectorWidth().getDim(), -1); // No user-imposed limit.
   privateBridge(
      west_to_east,
      NULL,
      west_to_cent,
      cent_to_east,
      east_to_cent,
      cent_to_west,
      false,
      zero_vector,
      false,
      zero_vector,
      connector_width_limit);
}

/*
 ***********************************************************************
 *
 *                           west to east
 * (west mapped_box_level)  ------------------> (east mapped_box_level)
 *                       ^  <------------------ ^
 *                        \    east to west    /
 *                         \                  /
 *           center to west \                / center to east
 *                           \              /
 *                            \            /
 *                      (center mapped_box_level)
 *
 * Bridge operation is in two phases, discovery and
 * sharing.  The discovery phase loops through local
 * Boxes in the center and comparing the west and east neighbors
 * for overlaps.  Local overlaps are stored immediately.
 * Remote overlaps are placed in messages to be sent to appropriate
 * processors by the sharing phase.
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::privateBridge(
   Connector& west_to_east,
   Connector* east_to_west,
   const Connector& west_to_cent,
   const Connector& cent_to_east,
   const Connector& east_to_cent,
   const Connector& cent_to_west,
   bool west_nesting_is_known,
   const IntVector& cent_growth_to_nest_west,
   bool east_nesting_is_known,
   const IntVector& cent_growth_to_nest_east,
   const IntVector& connector_width_limit) const
{
   t_bridge->barrierAndStart();

   if (s_print_steps == 'y') {
      std::string dbgbord("bridge->  ");
      tbox::plog
      << "bridge west:\n" << west_to_cent.getBase().format(dbgbord, 3)
      << "bridge east:\n" << east_to_cent.getBase().format(dbgbord, 3)
      << "bridge center:\n" << cent_to_west.getBase().format(dbgbord, 3)
      << "bridge west_to_cent:\n" << west_to_cent.format(dbgbord, 3)
      << "bridge cent_to_west:\n" << cent_to_west.format(dbgbord, 3)
      << "bridge cent_to_east:\n" << cent_to_east.format(dbgbord, 3)
      << "bridge east_to_cent:\n" << cent_to_east.format(dbgbord, 3);
   }

   privateBridge_checkParameters(
      west_to_cent,
      cent_to_east,
      east_to_cent,
      cent_to_west);

   const BoxLevel& cent = cent_to_west.getBase();
   const BoxLevel& west = cent_to_west.getHead();
   const BoxLevel& east = cent_to_east.getHead();
   const IntVector& cent_refinement_ratio = cent.getRefinementRatio();
   const IntVector& west_refinement_ratio = west.getRefinementRatio();
   const IntVector& east_refinement_ratio = east.getRefinementRatio();

   const tbox::Dimension& dim(connector_width_limit.getDim());
   const tbox::SAMRAI_MPI& mpi(cent.getMPI());

   const IntVector& zero_vector(IntVector::getZero(dim));

   const IntVector finest_refinement_ratio =
      IntVector::max(
         cent_refinement_ratio,
         IntVector::max(west_refinement_ratio, east_refinement_ratio));

   /*
    * Using the bridge theorem, compute the largest bridge width for
    * which we can guarantee discovering all the overlaps (when
    * nesting is satisfied).  If either the east or west
    * BoxLevel's nesting in the center is known, compute the
    * output width by the bridge theorem, and use the bigger one.  If
    * neither is known, we assume that both east and west nest in
    * center, and just to do something reasonable.
    */
   IntVector output_width1(dim, 0), output_width2(dim, 0);
   if (west_nesting_is_known || east_nesting_is_known) {
      if (west_nesting_is_known) {
         output_width1 =
            cent_to_east.getConnectorWidth() - cent_growth_to_nest_west;
      }
      if (east_nesting_is_known) {
         output_width2 =
            cent_to_west.getConnectorWidth() - cent_growth_to_nest_east;
      }
      if (!(output_width1 >= zero_vector || output_width2 >= zero_vector)) {
         TBOX_ERROR("OverlapConnectorAlgorithm::privateBridge:\n"
            << "Useless nesting specifications!\n"
            << "Neither west nor east BoxLevel nest with enough\n"
            << "margin to guarantee finding all overlaps.\n"
            << "To ensure you understand completness is not guaranteed,\n"
            << "this is considered an error.  To proceed anyway and live\n"
            << "with potential incompleteness, use a bridge interface\n"
            << "that does not claim any nesting.  Or, you can specify\n"
            << "a different nesting claim (but don't enable sanity\n"
            << "checking, which will catch your fib).\n");
      }
   } else {
      output_width1 = cent_to_east.getConnectorWidth();
      output_width2 = cent_to_west.getConnectorWidth();
   }
   IntVector output_width_in_finest_refinement_ratio =
      IntVector::max(output_width1, output_width2) *
      finest_refinement_ratio / cent_refinement_ratio;

   /*
    * Reduce the output width to the user-specified width limit.  Note
    * that the width limit is specified in the coarser of the east and
    * west refinement ratios.
    */
   if (connector_width_limit >= IntVector::getZero(dim)) {
      const IntVector coarser_refinement_ratio =
         IntVector::min(west_refinement_ratio, east_refinement_ratio);
      const IntVector width_limit_in_finest_refinement_ratio(
         connector_width_limit * finest_refinement_ratio / coarser_refinement_ratio);
      if (width_limit_in_finest_refinement_ratio > output_width_in_finest_refinement_ratio) {
         /*
          * If user specifies a width limit, he is probably assuming
          * that the bridge's allowable width is bigger.  If that is
          * not the case, this method won't crash, but it will give
          * bad results that result in elusive bugs.  Therefore, we
          * catch it immediately.
          */
         TBOX_ERROR("OverlapConnectorAlgorithm::privateBridge found input error:\n"
            << "The given connector width limit, " << connector_width_limit
            << " (" << width_limit_in_finest_refinement_ratio << " in finest index space)\n"
            << "is smaller than the width of the bridge, "
            << output_width_in_finest_refinement_ratio << " (in finest index space).");
      }
      output_width_in_finest_refinement_ratio.min(width_limit_in_finest_refinement_ratio);
   }

   const IntVector west_to_east_width = IntVector::ceilingDivide(
      output_width_in_finest_refinement_ratio,
      finest_refinement_ratio / west_refinement_ratio);
   const IntVector east_to_west_width = IntVector::ceilingDivide(
      output_width_in_finest_refinement_ratio,
      finest_refinement_ratio / east_refinement_ratio);


   const int rank = mpi.getRank();

   /*
    * Compute the reverse bridge (east_to_west) if it is given and is
    * distinct.
    */
   const bool compute_reverse =
      (east_to_west != NULL && east_to_west != &west_to_east);

   /*
    * Owners we have to exchange information with are the ones
    * owning east/west Boxes visible to the local process.
    */
   std::set<int> incoming_ranks, outgoing_ranks;
   cent_to_west.getLocalOwners(outgoing_ranks);
   west_to_cent.getLocalOwners(incoming_ranks);
   if (compute_reverse && &cent_to_west != &cent_to_east) {
      cent_to_east.getLocalOwners(outgoing_ranks);
      east_to_cent.getLocalOwners(incoming_ranks);
   }
   outgoing_ranks.erase(rank);
   incoming_ranks.erase(rank);

   /*
    * Create BoxContainers which will later be used to initialize the search trees
    * for visible east and west neighbors:
    * visible_west_nabrs and visible_east_nabrs.
    */
   t_bridge_discover_get_neighbors->start();
   bool ordered = true;
   NeighborSet visible_west_nabrs(ordered), visible_east_nabrs(ordered);
   cent_to_west.getLocalNeighbors(visible_west_nabrs);
   cent_to_east.getLocalNeighbors(visible_east_nabrs);
   t_bridge_discover_get_neighbors->stop();

   /*
    * Setup the output Connectors' head, base, width and type.  Add
    * overlaps below as they are discovered or received from other procs.
    */
   west_to_east.clearNeighborhoods();
   west_to_east.setConnectorType(Connector::UNKNOWN);
   west_to_east.setBase(west);
   west_to_east.setHead(east);
   west_to_east.setWidth(west_to_east_width, true);
   if (compute_reverse) {
      east_to_west->clearNeighborhoods();
      east_to_west->setConnectorType(Connector::UNKNOWN);
      east_to_west->setBase(east);
      east_to_west->setHead(west);
      east_to_west->setWidth(east_to_west_width, true);
#ifdef DEBUG_CHECK_ASSERTIONS
      if (west_refinement_ratio / east_refinement_ratio * east_refinement_ratio == west_refinement_ratio ||
          east_refinement_ratio / west_refinement_ratio * west_refinement_ratio == east_refinement_ratio) {
         /*
          * If it's possible to make west<==>east transposes, it
          * should happen.  The requirement is that one refinement ratio is
          * an IntVector times the other.
          */
         west_to_east.isTransposeOf(*east_to_west);
         east_to_west->isTransposeOf(west_to_east);
      }
#endif
   }


   /*
    * Set up communication mechanism and post receives.
    * Note that in all_comms, all the incoming_comm come
    * first, the outgoing_comm later.
    */

   tbox::AsyncCommStage comm_stage;
   tbox::AsyncCommPeer<int>* all_comms(NULL);

   t_bridge_share->start();
   t_bridge_setup_comm->start();

   setupCommunication(
      all_comms,
      comm_stage,
      mpi,
      incoming_ranks,
      outgoing_ranks,
      t_bridge_MPI_wait,
      s_operation_mpi_tag);

   t_bridge_setup_comm->stop();
   t_bridge_share->stop();

   /*
    * Data for caching relationship removal messages.  This
    * temporary object holds data computed when removing neighbors,
    * to be used when sending out the relationship removal
    * information.
    */
   std::map<int, std::vector<int> > neighbor_removal_mesg;

   /*
    * First step: Remove neighbor data for Boxes that are
    * going away and cache information to be sent out.
    */
   privateModify_removeAndCache(
      neighbor_removal_mesg,
      west_to_east,
      east_to_west,
      cent_to_east);

   privateBridge_discoverAndSend(
      neighbor_removal_mesg,
      west_to_east,
      east_to_west,
      incoming_ranks,
      outgoing_ranks,
      all_comms,
      visible_west_nabrs,
      visible_east_nabrs);

   t_bridge_share->start();

   receiveAndUnpack(
      west_to_east,
      east_to_west,
      incoming_ranks,
      all_comms,
      comm_stage,
      t_bridge_receive_and_unpack);

   t_bridge_share->stop();

   delete[] all_comms;

   if (d_sanity_check_method_postconditions) {
      west_to_east.assertConsistencyWithBase();
      west_to_east.assertConsistencyWithHead();
      if (compute_reverse) {
         east_to_west->assertConsistencyWithBase();
         east_to_west->assertConsistencyWithHead();
         east_to_west->assertTransposeCorrectness(west_to_east, true);
      }
   }

   t_bridge->stop();
}

/*
 ***********************************************************************
 *
 *                           west to east
 * (west mapped_box_level)  ------------------> (east mapped_box_level)
 *                       ^  <------------------ ^
 *                        \    east to west    /
 *                         \                  /
 *           center to west \                / center to east
 *                           \              /
 *                            \            /
 *                      (center mapped_box_level)
 *
 * Bridge operation is in two phases, discovery and
 * sharing.  The discovery phase loops through local
 * Boxes in the center and comparing the west and east neighbors
 * for overlaps.  Local overlaps are stored immediately.
 * Remote overlaps are placed in messages to be sent to appropriate
 * processors by the sharing phase.
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::privateBridge(
   Connector& west_to_cent,
   const Connector& cent_to_east,
   const Connector& east_to_cent,
   Connector& cent_to_west,
   bool west_nesting_is_known,
   const IntVector& cent_growth_to_nest_west,
   bool east_nesting_is_known,
   const IntVector& cent_growth_to_nest_east,
   const IntVector& connector_width_limit) const
{
   t_bridge->barrierAndStart();

   if (s_print_steps == 'y') {
      std::string dbgbord("bridge->  ");
      tbox::plog
      << "bridge west:\n" << west_to_cent.getBase().format(dbgbord, 3)
      << "bridge east:\n" << east_to_cent.getBase().format(dbgbord, 3)
      << "bridge center:\n" << cent_to_west.getBase().format(dbgbord, 3)
      << "bridge west_to_cent:\n" << west_to_cent.format(dbgbord, 3)
      << "bridge cent_to_west:\n" << cent_to_west.format(dbgbord, 3)
      << "bridge cent_to_east:\n" << cent_to_east.format(dbgbord, 3)
      << "bridge east_to_cent:\n" << east_to_cent.format(dbgbord, 3);
   }

   privateBridge_checkParameters(
      west_to_cent,
      cent_to_east,
      east_to_cent,
      cent_to_west);

   const BoxLevel& cent = cent_to_west.getBase();
   const BoxLevel& west = cent_to_west.getHead();
   const BoxLevel& east = cent_to_east.getHead();
   const IntVector& cent_refinement_ratio = cent.getRefinementRatio();
   const IntVector& west_refinement_ratio = west.getRefinementRatio();
   const IntVector& east_refinement_ratio = east.getRefinementRatio();

   const tbox::Dimension& dim(connector_width_limit.getDim());
   const tbox::SAMRAI_MPI& mpi(cent.getMPI());

   const IntVector& zero_vector(IntVector::getZero(dim));

   const IntVector finest_refinement_ratio =
      IntVector::max(
         cent_refinement_ratio,
         IntVector::max(west_refinement_ratio, east_refinement_ratio));

   /*
    * Using the bridge theorem, compute the largest bridge width for
    * which we can guarantee discovering all the overlaps (when
    * nesting is satisfied).  If either the east or west
    * BoxLevel's nesting in the center is known, compute the
    * output width by the bridge theorem, and use the bigger one.  If
    * neither is known, we assume that both east and west nest in
    * center, and just to do something reasonable.
    */
   IntVector output_width1(dim, 0), output_width2(dim, 0);
   if (west_nesting_is_known || east_nesting_is_known) {
      if (west_nesting_is_known) {
         output_width1 =
            cent_to_east.getConnectorWidth() - cent_growth_to_nest_west;
      }
      if (east_nesting_is_known) {
         output_width2 =
            cent_to_west.getConnectorWidth() - cent_growth_to_nest_east;
      }
      if (!(output_width1 >= zero_vector || output_width2 >= zero_vector)) {
         TBOX_ERROR("OverlapConnectorAlgorithm::privateBridge:\n"
            << "Useless nesting specifications!\n"
            << "Neither west nor east BoxLevel nest with enough\n"
            << "margin to guarantee finding all overlaps.\n"
            << "To ensure you understand completness is not guaranteed,\n"
            << "this is considered an error.  To proceed anyway and live\n"
            << "with potential incompleteness, use a bridge interface\n"
            << "that does not claim any nesting.  Or, you can specify\n"
            << "a different nesting claim (but don't enable sanity\n"
            << "checking, which will catch your fib).\n");
      }
   } else {
      output_width1 = cent_to_east.getConnectorWidth();
      output_width2 = cent_to_west.getConnectorWidth();
   }
   IntVector output_width_in_finest_refinement_ratio =
      IntVector::max(output_width1, output_width2) *
      finest_refinement_ratio / cent_refinement_ratio;

   /*
    * Reduce the output width to the user-specified width limit.  Note
    * that the width limit is specified in the coarser of the east and
    * west refinement ratios.
    */
   if (connector_width_limit >= IntVector::getZero(dim)) {
      const IntVector coarser_refinement_ratio =
         IntVector::min(west_refinement_ratio, east_refinement_ratio);
      const IntVector width_limit_in_finest_refinement_ratio(
         connector_width_limit * finest_refinement_ratio / coarser_refinement_ratio);
      if (width_limit_in_finest_refinement_ratio > output_width_in_finest_refinement_ratio) {
         /*
          * If user specifies a width limit, he is probably assuming
          * that the bridge's allowable width is bigger.  If that is
          * not the case, this method won't crash, but it will give
          * bad results that result in elusive bugs.  Therefore, we
          * catch it immediately.
          */
         TBOX_ERROR("OverlapConnectorAlgorithm::privateBridge found input error:\n"
            << "The given connector width limit, " << connector_width_limit
            << " (" << width_limit_in_finest_refinement_ratio << " in finest index space)\n"
            << "is smaller than the width of the bridge, "
            << output_width_in_finest_refinement_ratio << " (in finest index space).");
      }
      output_width_in_finest_refinement_ratio.min(width_limit_in_finest_refinement_ratio);
   }

   const IntVector west_to_east_width = IntVector::ceilingDivide(
      output_width_in_finest_refinement_ratio,
      finest_refinement_ratio / west_refinement_ratio);
   const IntVector east_to_west_width = IntVector::ceilingDivide(
      output_width_in_finest_refinement_ratio,
      finest_refinement_ratio / east_refinement_ratio);


   const int rank = mpi.getRank();


   /*
    * Compute the reverse bridge (east_to_west) if it is given and is
    * distinct.
    */
   const bool compute_reverse = (&cent_to_west != &west_to_cent);

   /*
    * Owners we have to exchange information with are the ones
    * owning east/west Boxes visible to the local process.
    */
   std::set<int> incoming_ranks, outgoing_ranks;
   cent_to_west.getLocalOwners(outgoing_ranks);
   west_to_cent.getLocalOwners(incoming_ranks);
   if (compute_reverse && &cent_to_west != &cent_to_east) {
      cent_to_east.getLocalOwners(outgoing_ranks);
      east_to_cent.getLocalOwners(incoming_ranks);
   }
   outgoing_ranks.erase(rank);
   incoming_ranks.erase(rank);

   /*
    * Create BoxContainers which will later be used to initialize the search trees
    * for visible east and west neighbors:
    * visible_west_nabrs and visible_east_nabrs.
    */
   t_bridge_discover_get_neighbors->start();
   bool ordered = true;
   NeighborSet visible_west_nabrs(ordered), visible_east_nabrs(ordered);
   cent_to_west.getLocalNeighbors(visible_west_nabrs);
   cent_to_east.getLocalNeighbors(visible_east_nabrs);
   t_bridge_discover_get_neighbors->stop();

   /*
    * Initialize the output Connectors without overlaps.  Add overlaps
    * below as they are discovered or received from other procs.
    */
   west_to_cent.clearNeighborhoods();
   west_to_cent.setConnectorType(Connector::UNKNOWN);
   west_to_cent.setBase(west);
   west_to_cent.setHead(east);
   west_to_cent.setWidth(west_to_east_width, true);
   if (compute_reverse) {
      cent_to_west.clearNeighborhoods();
      cent_to_west.setConnectorType(Connector::UNKNOWN);
      cent_to_west.setBase(east);
      cent_to_west.setHead(west);
      cent_to_west.setWidth(east_to_west_width, true);
#ifdef DEBUG_CHECK_ASSERTIONS
      if (west_refinement_ratio / east_refinement_ratio * east_refinement_ratio == west_refinement_ratio ||
          east_refinement_ratio / west_refinement_ratio * west_refinement_ratio == east_refinement_ratio) {
         /*
          * If it's possible to make west<==>east transposes, it
          * should happen.  The requirement is that one refinement ratio is
          * an IntVector times the other.
          */
         TBOX_ASSERT(west_to_cent.isTransposeOf(cent_to_west));
         TBOX_ASSERT(cent_to_west.isTransposeOf(west_to_cent));
      }
#endif
   }


   /*
    * Set up communication mechanism and post receives.
    * Note that in all_comms, all the incoming_comm come
    * first, the outgoing_comm later.
    */

   tbox::AsyncCommStage comm_stage;
   tbox::AsyncCommPeer<int>* all_comms(NULL);

   t_bridge_share->start();
   t_bridge_setup_comm->start();

   setupCommunication(
      all_comms,
      comm_stage,
      mpi,
      incoming_ranks,
      outgoing_ranks,
      t_bridge_MPI_wait,
      s_operation_mpi_tag);

   t_bridge_setup_comm->stop();
   t_bridge_share->stop();

   /*
    * Data for caching relationship removal messages.  This
    * temporary object holds data computed when removing neighbors,
    * to be used when sending out the relationship removal
    * information.
    */
   std::map<int, std::vector<int> > neighbor_removal_mesg;

   /*
    * First step: Remove neighbor data for Boxes that are
    * going away and cache information to be sent out.
    */
   privateModify_removeAndCache(
      neighbor_removal_mesg,
      west_to_cent,
      &cent_to_west,
      cent_to_east);

   privateBridge_discoverAndSend(
      neighbor_removal_mesg,
      west_to_cent,
      &cent_to_west,
      incoming_ranks,
      outgoing_ranks,
      all_comms,
      visible_west_nabrs,
      visible_east_nabrs);

   t_bridge_share->start();

   receiveAndUnpack(
      west_to_cent,
      &cent_to_west,
      incoming_ranks,
      all_comms,
      comm_stage,
      t_bridge_receive_and_unpack);

   t_bridge_share->stop();

   delete[] all_comms;

   if (d_sanity_check_method_postconditions) {
      west_to_cent.assertConsistencyWithBase();
      west_to_cent.assertConsistencyWithHead();
      if (compute_reverse) {
         cent_to_west.assertConsistencyWithBase();
         cent_to_west.assertConsistencyWithHead();
         cent_to_west.assertTransposeCorrectness(west_to_cent, true);
      }
   }

   t_bridge->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::privateBridge_checkParameters(
   const Connector& west_to_cent,
   const Connector& cent_to_east,
   const Connector& east_to_cent,
   const Connector& cent_to_west) const
{
   const BoxLevel& cent = cent_to_west.getBase();

   /*
    * Ensure that Connectors incident to and from the center agree on
    * what the center is.  This can be an expensive (though still
    * scalable) check, so we only check in debug mode, unless debugging,
    * in which it can be independently enabled in any mode.
    */
   if (cent != cent_to_east.getBase() ||
       cent != east_to_cent.getHead() ||
       cent != west_to_cent.getHead()) {
      TBOX_ERROR("Bad input for OverlapConnectorAlgorithm::bridge:\n"
         << "Given Connectors to base and head of bridge are not incident\n"
         << "from the same center in OverlapConnectorAlgorithm::bridge:\n"
         << "west_to_cent is  TO  " << &west_to_cent.getHead() << "\n"
         << "cent_to_east is FROM " << &cent_to_east.getBase() << "\n"
         << "east_to_cent is  TO  " << &east_to_cent.getHead() << "\n"
         << "cent_to_west is FROM " << &cent_to_west.getBase() << "\n"
         );
   }
   /*
    * Ensure that head and base mapped_box_levels in argument agree with
    * head and base in the object.
    */
   if (cent_to_west.getHead() != west_to_cent.getBase()) {
      TBOX_ERROR("Bad input for OverlapConnectorAlgorithm::bridge:\n"
         << "Given Connectors to and from base of bridge do not refer\n"
         << "to the base of the bridge in OverlapConnectorAlgorithm::bridge:\n"
         << "west_to_cent is FROM " << &west_to_cent.getBase() << "\n"
         << "cent_to_west is  TO  " << &cent_to_west.getHead() << "\n"
         );
   }
   if (cent_to_east.getHead() != east_to_cent.getBase()) {
      TBOX_ERROR("Bad input for OverlapConnectorAlgorithm::bridge:\n"
         << "Given Connectors to and from head of bridge do not refer\n"
         << "to the head of the bridge in OverlapConnectorAlgorithm::bridge:\n"
         << "east_to_cent is FROM " << &east_to_cent.getBase() << "\n"
         << "cent_to_east is  TO  " << &cent_to_east.getHead() << "\n"
         );
   }
   if (!west_to_cent.isTransposeOf(cent_to_west)) {
      TBOX_ERROR("Bad input for OverlapConnectorAlgorithm::bridge:\n"
         << "Given Connectors between base and center of bridge\n"
         << "are not transposes of each other.\n"
         << "See OverlapConnectorAlgorithm::isTransposeOf().\n"
         );
   }
   if (!east_to_cent.isTransposeOf(cent_to_east)) {
      TBOX_ERROR("Bad input for OverlapConnectorAlgorithm::bridge:\n"
         << "Given Connectors between head and center of bridge\n"
         << "are not transposes of each other.\n"
         << "See OverlapConnectorAlgorithm::isTransposeOf().\n"
         );
   }

   // Expensive sanity checks:
   if (d_sanity_check_method_preconditions) {
      west_to_cent.assertTransposeCorrectness(cent_to_west);
      cent_to_west.assertTransposeCorrectness(west_to_cent);
      east_to_cent.assertTransposeCorrectness(cent_to_east);
      cent_to_east.assertTransposeCorrectness(east_to_cent);
   }
}

/*
 ***********************************************************************
 * Remove relationships from resulting overlap.  Cache outgoing
 * information in message buffers.
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::privateModify_removeAndCache(
   std::map<int, std::vector<int> >& neighbor_removal_mesg,
   Connector& overlap_connector,
   Connector* overlap_connector_transpose,
   const Connector& misc_connector) const
{
   t_bridge_remove_and_cache->start();

   NULL_USE(neighbor_removal_mesg); 
   NULL_USE(overlap_connector); 
   NULL_USE(overlap_connector_transpose); 
   NULL_USE(misc_connector); 
  /*
   * As the overlap relationships are empty to start there are never any
   * that need to be deleted.
   */
   t_bridge_remove_and_cache->stop();
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::privateBridge_discoverAndSend(
   std::map<int, std::vector<int> > neighbor_removal_mesg,
   Connector& west_to_east,
   Connector* east_to_west,
   std::set<int>& incoming_ranks,
   std::set<int>& outgoing_ranks,
   tbox::AsyncCommPeer<int> all_comms[],
   NeighborSet& visible_west_nabrs,
   NeighborSet& visible_east_nabrs) const
{
   if (!visible_west_nabrs.isEmpty() || !visible_east_nabrs.isEmpty()) {

      /*
       * Discover overlaps.  Overlaps are either locally stored or
       * packed into a message for sending.
       */

      t_bridge_discover->start();

      if (s_print_steps == 'y') {
         tbox::plog << "Before building RBBTs:\n"
                    << "visible_west_nabrs:" << visible_west_nabrs.format("\n  ")
                    << "visible_east_nabrs:" << visible_east_nabrs.format("\n  ")
                    << std::endl;
      }

      bool compute_reverse =
         (east_to_west != NULL && east_to_west != &west_to_east);

      const BoxLevel& east(west_to_east.getBase());
      const boost::shared_ptr<const BaseGridGeometry> &grid_geometry(
         east.getGridGeometry());

      const tbox::Dimension& dim(east.getDim());
      const int rank = east.getMPI().getRank();
      const int nproc = east.getMPI().getSize();

      t_bridge_discover_form_rbbt->start();
      const BoxContainer east_rbbt(visible_east_nabrs);
      east_rbbt.makeTree(grid_geometry.get());
      // Note: west_rbbt only needed when compute_reverse is true.
      BoxContainer empty_nabrs(true);
      const BoxContainer west_rbbt(compute_reverse ? visible_west_nabrs : empty_nabrs);
      west_rbbt.makeTree(grid_geometry.get());
      t_bridge_discover_form_rbbt->stop();

      /*
       * Iterators west_ni and east_ni point to the west/east
       * Box whose neighbors are being sought.  If we are not
       * interested in the east-->west connector, then east_ni will
       * be unused.
       */
      NeighborSet::iterator west_ni(visible_west_nabrs);
      NeighborSet::iterator east_ni(visible_east_nabrs);
      /*
       * Local process can find some neighbors for the (local and
       * remote) Boxes in visible_west_nabrs and
       * visible_east_nabrs.  We loop through the visible_west_nabrs
       * and compare each to visible_east_nabrs, looking for overlaps.
       * Then vice versa.
       *
       * Looping through the NeighborSets is like looping through
       * their owners, since they are ordered by owners first.  As an
       * optimization measure, start loop on the first owner with
       * higher rank than the local rank.  This avoid the higher-end
       * ranks from having to wait for messages at the beginning and
       * the lower-end ranks from having to wait for messages at the
       * end.  After the highest rank owner has been handled, continue
       * at the beginning and do the remaining.  (If local rank is
       * highest of all owners of the visible Boxes, start at
       * the beginning.)
       */
      const Box start_loop_here(dim, GlobalId(LocalId::getZero(), rank + 1));
      west_ni = visible_west_nabrs.lowerBound(start_loop_here);
      if (compute_reverse) {
         east_ni = visible_east_nabrs.lowerBound(start_loop_here);
      }

      if (west_ni == visible_west_nabrs.end() &&
          (!compute_reverse || east_ni == visible_east_nabrs.end())) {
         /*
          * There are no visible Boxes owned by rank higher than
          * local process.  So loop from the beginning.
          */
         west_ni = visible_west_nabrs.begin();
         east_ni = visible_east_nabrs.begin();
      }

      /*
       * Set send_comm_idx to reference the first outgoing rank in all_comms.
       * It will be incremented to correpond to the rank whose overlaps
       * are being searched for.
       */
      int send_comm_idx = static_cast<int>(incoming_ranks.size());

#ifdef DEBUG_CHECK_ASSERTIONS
      // Owners that we have sent messages to.  Used for debugging.
      std::set<int> owners_sent_to;
#endif

      /*
       * Loop until all visible neighbors have their neighbors
       * searched for.  But only do this for the east Boxes if
       * we are actively seeking neighbor data for them.
       */
      while ((west_ni != visible_west_nabrs.end()) ||
             (compute_reverse && east_ni != visible_east_nabrs.end())) {

         /*
          * curr_owner is the owner whose neighbors is currently
          * being searched for.  It should be the owner of the
          * next west or east Box in our cyclic-type looping.
          */
         int curr_owner = nproc; // Start with invalid value.
         if (west_ni != visible_west_nabrs.end() &&
             curr_owner > west_ni->getOwnerRank()) {
            curr_owner = west_ni->getOwnerRank();
         }
         if (compute_reverse) {
            if (east_ni != visible_east_nabrs.end() &&
                curr_owner > east_ni->getOwnerRank()) {
               curr_owner = east_ni->getOwnerRank();
            }
         }
         if (s_print_steps == 'y') {
            tbox::plog << "cur_owner set to " << curr_owner << std::endl;
         }

         TBOX_ASSERT(curr_owner < nproc);
         TBOX_ASSERT(curr_owner > -1);
         TBOX_ASSERT(owners_sent_to.find(curr_owner) == owners_sent_to.end());

         /*
          * Set up send_message to contain info discovered
          * locally but needed by curr_owner.
          *
          * Content of send_mesg:
          * - neighbor-removal section cached in neighbor_removal_mesg.
          * - offset to the reference section (see below)
          * - number of west mapped_boxes for which neighbors are found
          * - number of east mapped_boxes for which neighbors are found
          *   - id of west/east mapped_box
          *   - number of neighbors found for west/east mapped_box.
          *     - BoxId of neighbors found.
          *       Boxes of these found neighbors are given in the
          *       reference section of the message.
          * - reference section: all the Boxes referenced as
          *   neighbors (accumulated in referenced_west_nabrs
          *   and referenced_east_nabrs).
          *   - number of referenced west neighbors
          *   - number of referenced east neighbors
          *   - referenced west neighbors
          *   - referenced east neighbors
          *
          * The purpose of factoring out info on the neighbors referenced
          * is to reduce redundant data that can eat up lots of memory
          * and message passing bandwidth when there are lots of Boxes
          * with the same neighbors.
          */
         std::vector<int> send_mesg;

         /*
          * The first section of the send_mesg is the remote neighbor-removal
          * section (computed above).
          */
         if (curr_owner != rank) {
            // swap( send_mesg, neighbor_removal_mesg[curr_owner] );
            // We could use swap, but assign instead to leave data for debugging.
            send_mesg = neighbor_removal_mesg[curr_owner];
            if (send_mesg.empty()) {
               // No neighbor-removal data found for curr_owner.
               send_mesg.insert(send_mesg.end(), 0);
            }
         }

         // Indices of certain positions in send_mesg.
         const int idx_offset_to_ref = static_cast<int>(send_mesg.size());
         const int idx_num_west_mapped_boxes = idx_offset_to_ref + 1;
         const int idx_num_east_mapped_boxes = idx_offset_to_ref + 2;
         send_mesg.insert(send_mesg.end(), 3, 0);

         // Mapped_boxes referenced in the message, used when adding ref section.
         BoxContainer referenced_west_nabrs;
         BoxContainer referenced_east_nabrs;

         t_bridge_discover_find_overlaps->start();

         if (s_print_steps == 'y') {
            tbox::plog << "Finding west --> east overlaps for owner "
                       << curr_owner << std::endl;
         }

         // Find neighbors for all west mapped_boxes owned by curr_owner.
         privateBridge_findOverlapsForOneProcess(
            curr_owner,
            visible_west_nabrs,
            west_ni,
            send_mesg,
            idx_num_west_mapped_boxes,
            west_to_east,
            referenced_east_nabrs,
            east_rbbt);

         // Find neighbors for all east mapped_boxes owned by curr_owner.
         if (compute_reverse) {
            if (s_print_steps == 'y') {
               tbox::plog << "Finding west <-- east overlaps for owner "
                          << curr_owner << std::endl;
            }
            privateBridge_findOverlapsForOneProcess(
               curr_owner,
               visible_east_nabrs,
               east_ni,
               send_mesg,
               idx_num_east_mapped_boxes,
               *east_to_west,
               referenced_west_nabrs,
               west_rbbt);
         }

         t_bridge_discover_find_overlaps->stop();

         if (curr_owner != rank) {
            /*
             * Send discoveries to the curr_owner.
             */

            t_bridge_discover->stop();
            t_bridge_share->start();

            /*
             * Find the communication object by increasing send_comm_idx
             * (cyclically) until it corresponds to curr_owner.
             */
            while (all_comms[send_comm_idx].getPeerRank() != curr_owner) {
               ++send_comm_idx;
               if (send_comm_idx == static_cast<int>(incoming_ranks.size() +
                                                     outgoing_ranks.size())) {
                  send_comm_idx -= static_cast<int>(outgoing_ranks.size());
               }
            }
            tbox::AsyncCommPeer<int>& outgoing_comm = all_comms[send_comm_idx];
            TBOX_ASSERT(outgoing_comm.getPeerRank() == curr_owner);

            sendDiscoveryToOneProcess(
               send_mesg,
               idx_offset_to_ref,
               referenced_east_nabrs,
               referenced_west_nabrs,
               outgoing_comm,
               dim);

            TBOX_ASSERT(owners_sent_to.find(curr_owner) == owners_sent_to.end());
#ifdef DEBUG_CHECK_ASSERTIONS
            owners_sent_to.insert(curr_owner);
#endif

            t_bridge_share->stop();
            t_bridge_discover->start();

         } // Block to send discoveries to curr_owner.

         /*
          * If we come to the end of visible Boxes, go back and
          * work on the Boxes owned by processors with lower rank
          * than the local rank.  (This is part of the optimization
          * to reduce communication time.)
          */
         if (west_ni == visible_west_nabrs.end() &&
             (!compute_reverse || east_ni == visible_east_nabrs.end())) {
            /*
             * There are no Boxes that owned by rank higher than
             * local process and that we still need to find neighbors
             * for.  So loop from the beginning.
             */
            west_ni = visible_west_nabrs.begin();
            east_ni = visible_east_nabrs.begin();
         }

      } // Loop through visible neighbors.

      t_bridge_discover->stop();

   }
}

/*
 ***********************************************************************
 *
 * Find overlaps from visible_base_nabrs to head_rbbt.  Find only
 * overlaps for Boxes owned by owner_rank.
 *
 * On input, base_ni points to the first Box in visible_base_nabrs
 * owned by owner_rank.  Increment base_ni past those Boxes
 * processed and remove them from visible_base_nabrs.
 *
 * Save local and semilocal overlaps in bridging_connector.  For
 * remote overlaps, pack in send_mesg, add head Box to
 * referenced_head_nabrs and increment
 * send_mesg[remote_mapped_box_counter_index].
 *
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::privateBridge_findOverlapsForOneProcess(
   const int owner_rank,
   NeighborSet& visible_base_nabrs,
   NeighborSet::iterator& base_ni,
   std::vector<int>& send_mesg,
   const int remote_mapped_box_counter_index,
   Connector& bridging_connector,
   NeighborSet& referenced_head_nabrs,
   const BoxContainer& head_rbbt) const
{
   const IntVector &head_refinement_ratio(bridging_connector.getHead().getRefinementRatio());

   bool refine_base = false;
   bool coarsen_base = false;
   if (bridging_connector.getHead().getRefinementRatio() ==
       bridging_connector.getBase().getRefinementRatio()) {
      /*
       * Don't do any coarsen/refine because head and base have same
       * refinement ratio.
       */
   } else if (bridging_connector.getHead().getRefinementRatio() <=
              bridging_connector.getBase().getRefinementRatio()) {
      coarsen_base = true;
   } else if (bridging_connector.getHead().getRefinementRatio() >=
              bridging_connector.getBase().getRefinementRatio()) {
      refine_base = true;
   }
#ifdef DEBUG_CHECK_ASSERTIONS
   else {
      TBOX_ERROR("Can't coarsen in one dimension and refine in another");
   }
#endif

   //std::vector<Box> found_nabrs, scratch_found_nabrs; // Should be made a member to avoid repetitive alloc/dealloc.  Reserve in privateBridge and used here.
   BoxContainer found_nabrs, scratch_found_nabrs; // Should be made a member to avoid repetitive alloc/dealloc.  Reserve in privateBridge and used here.

   while (base_ni != visible_base_nabrs.end() &&
          base_ni->getOwnerRank() == owner_rank) {
      const Box& base_mapped_box = *base_ni;
      if (s_print_steps == 'y') {
         tbox::plog << "Finding neighbors for non-periodic base_mapped_box "
                    << base_mapped_box << std::endl;
      }
      Box base_box = base_mapped_box;
      base_box.grow(bridging_connector.getConnectorWidth());
      if (refine_base) {
         base_box.refine(bridging_connector.getRatio());
      }
      else if (coarsen_base) {
         base_box.coarsen(bridging_connector.getRatio());
      }
      found_nabrs.clear();
      head_rbbt.findOverlapBoxes(found_nabrs, base_box, // base_box.getBlockId(),
                                 head_refinement_ratio,
                                 true /* include singularity block neighbors */ );
      if (s_print_steps == 'y') {
         tbox::plog << "Found " << found_nabrs.size() << " neighbors:";
         found_nabrs.print(tbox::plog);
         //BoxContainerUtils::recursivePrintBoxVector(found_nabrs, tbox::plog, "\n ");
         tbox::plog << std::endl;
      }
      if (!found_nabrs.isEmpty()) {
         if (base_mapped_box.isPeriodicImage()) {
            privateBridge_unshiftOverlappingNeighbors(
               base_mapped_box,
               found_nabrs,
               scratch_found_nabrs,
               bridging_connector.getHead().getRefinementRatio());
         }
         if (owner_rank != bridging_connector.getMPI().getRank()) {
            // Pack up info for sending.
            ++send_mesg[remote_mapped_box_counter_index];
            const int subsize = 3 +
               BoxId::commBufferSize() * static_cast<int>(found_nabrs.size());
            send_mesg.insert(send_mesg.end(), subsize, -1);
            int* submesg = &send_mesg[send_mesg.size() - subsize];
            *(submesg++) = base_mapped_box.getLocalId().getValue();
            *(submesg++) = base_mapped_box.getBlockId().getBlockValue();
            *(submesg++) = static_cast<int>(found_nabrs.size());
            for (BoxContainer::const_iterator na = found_nabrs.begin();
                 na != found_nabrs.end(); ++na) {
               const Box& head_nabr = *na;
               referenced_head_nabrs.insert(head_nabr);
               head_nabr.getId().putToIntBuffer(submesg);
               submesg += BoxId::commBufferSize();
            }
         }
         else {
            // Save neighbor info locally.
            BoxId unshifted_base_mapped_box_id;
            if (!base_mapped_box.isPeriodicImage()) {
               unshifted_base_mapped_box_id = base_mapped_box.getId();
            }
            else {
               unshifted_base_mapped_box_id.initialize(
                  base_mapped_box.getLocalId(),
                  base_mapped_box.getOwnerRank(),
                  PeriodicId::zero());
            }
            // Add found neighbors for base_mapped_box.
            if (!found_nabrs.isEmpty()) {
               Connector::NeighborhoodIterator base_box_itr =
                  bridging_connector.makeEmptyLocalNeighborhood(
                     unshifted_base_mapped_box_id);
               for (BoxContainer::const_iterator na = found_nabrs.begin();
                    na != found_nabrs.end(); ++na) {
                  bridging_connector.insertLocalNeighbor(*na, base_box_itr);
               }
            }
         }
      }
      if (s_print_steps == 'y') {
         tbox::plog << "Erasing visible base nabr " << (*base_ni) << std::endl;
      }
      visible_base_nabrs.erase(base_ni++);
      if (s_print_steps == 'y') {
         if (base_ni == visible_base_nabrs.end()) {
            tbox::plog << "Next base nabr: end" << std::endl;
         }
         else {
            tbox::plog << "Next base nabr: " << *base_ni << std::endl;
         }
      }

   }
}

/*
 ***********************************************************************
 * Shift neighbors by amount equal and opposite of a Box's shift so that
 * they become neighbors of the unshifed mapped_box.  If this results in a
 * neighbor shift that is not in the shift catalog, discard the neighbor.
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::privateBridge_unshiftOverlappingNeighbors(
   const Box& mapped_box,
   BoxContainer& neighbors,
   BoxContainer& scratch_space,
   const IntVector& neighbor_refinement_ratio) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(mapped_box, neighbor_refinement_ratio);

   const PeriodicShiftCatalog* shift_catalog =
      PeriodicShiftCatalog::getCatalog(mapped_box.getDim());

   scratch_space.clear();
//   scratch_space.reserve(neighbors.size());
   for (BoxContainer::iterator na = neighbors.begin();
        na != neighbors.end(); ++na) {
      Box& nabr = *na;
      IntVector sum_shift =
         shift_catalog->shiftNumberToShiftDistance(nabr.getPeriodicId())
         - shift_catalog->shiftNumberToShiftDistance(mapped_box.getPeriodicId());
      const PeriodicId new_shift_number = shift_catalog->shiftDistanceToShiftNumber(sum_shift);
      if (new_shift_number.getPeriodicValue() != shift_catalog->getInvalidShiftNumber()) {
         nabr.initialize(nabr, new_shift_number, neighbor_refinement_ratio);
         scratch_space.pushBack(nabr);
      }
   }
   if (scratch_space.size() != neighbors.size()) {
      // We have discarded some neighbors due to invalid shift.
      neighbors.swap(scratch_space);
   }
}

/*
 ***********************************************************************
 * Checking is done as follows:
 *   - Rebuild the overlap containers using findOverlaps().
 *     Note that the rebuilt overlap set is complete.
 *   - Check the current overlap set against the rebuilt overlap set
 *     to find missing overlaps and extra overlaps.
 *
 * Currently, the rebuilt overlaps are rebuilt using findOverlaps().
 * Thus, it may be pointless to use this method as a check for that
 * method.
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::findOverlapErrors(
   const Connector& connector,
   Connector& missing,
   Connector& extra,
   bool ignore_self_overlap) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (!connector.getBase().isInitialized() ||
       !connector.getHead().isInitialized()) {
      TBOX_ERROR(
         "OverlapConnectorAlgorithm::findOverlapErrors: Cannot check overlaps\n"
         << "when base or head mapped_box_level is uninitialized.");
   }
#endif

   /*
    * Obtain a globalized version of the head for checking.
    */
   const BoxLevel& head = connector.getHead().getGlobalizedVersion();

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * Before checking on overlap errors, make sure the user gave a
    * valid Connector.
    *
    * Each neighbor set should correspond to a base mapped box.
    *
    * Each referenced neighbor should exist in the head.
    */
   size_t num_base_consistency_errors = connector.checkConsistencyWithBase();
   size_t num_head_consistency_errors = connector.checkConsistencyWithHead();
   if (num_base_consistency_errors > 0) {
      tbox::perr
      << "OverlapConnectorAlgorithm::findOverlapErrors: cannot check overlap errors\n"
      << "for inconsistent base data.\n";
   }
   if (num_head_consistency_errors > 0) {
      tbox::perr
      << "OverlapConnectorAlgorithm::findOverlapErrors: cannot check overlap errors\n"
      << "for inconsistent head data.\n";
   }
   if (num_base_consistency_errors || num_head_consistency_errors) {
      TBOX_ERROR(
         "OverlapConnectorAlgorithm::findOverlapErrors exiting due to\n"
         << "inconsistent data.\n"
         << "Base:\n" << connector.getBase().format("B->", 2)
         << "Head:\n" << connector.getHead().format("H->", 2)
         << "Connector:\n" << connector.format("C->", 3));
   }
#endif

   /*
    * Rebuild the overlap Connector for checking.
    */
   Connector rebuilt(
      connector.getBase(),
      connector.getHead(),
      connector.getConnectorWidth());
   findOverlaps_rbbt(rebuilt, head, ignore_self_overlap);

   /*
    * Check that the rebuilt overlaps match the existing overlaps.
    *
    * Currently, we use findOverlaps to rebuild the overlaps.
    * Thus, it may be pointless to use this method
    * as a check for that method.
    */
   Connector::computeNeighborhoodDifferences(
      extra,
      connector,
      rebuilt);
   TBOX_ASSERT(&extra.getBase() == &connector.getBase());
   TBOX_ASSERT(&extra.getHead() == &connector.getHead());
   Connector::computeNeighborhoodDifferences(
      missing,
      rebuilt,
      connector);
   TBOX_ASSERT(&missing.getBase() == &connector.getBase());
   TBOX_ASSERT(&missing.getHead() == &connector.getHead());
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
OverlapConnectorAlgorithm::assertOverlapCorrectness(
   const Connector& connector,
   bool ignore_self_overlap,
   bool assert_completeness,
   bool ignore_periodic_images) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (!connector.getBase().isInitialized() ||
       !connector.getHead().isInitialized()) {
      TBOX_ERROR(
         "OverlapConnectorAlgorithm::findOverlapErrors: Cannot check overlaps\n"
         << "when base or head mapped_box_level is uninitialized.");

   }
#endif

   int local_error_count =
      checkOverlapCorrectness(connector,
         ignore_self_overlap,
         assert_completeness,
         ignore_periodic_images);

   const tbox::SAMRAI_MPI& mpi(connector.getMPI());
   int max_error_count = local_error_count;
   int rank_of_max = mpi.getRank();
   if (mpi.getSize() > 1) {
      IntIntStruct send, recv;
      send.rank = recv.rank = mpi.getRank();
      send.i = local_error_count;
      mpi.Allreduce(&send, &recv, 1, MPI_2INT, MPI_MAXLOC);
      max_error_count = recv.i;
      rank_of_max = recv.rank;
   }
   if (max_error_count > 0) {
      std::string dbgbord;
      TBOX_ERROR(
         "OverlapConnectorAlgorithm::assertOverlapCorrectness found missing and/or extra overlaps."
         << "Error in connector, " << local_error_count
         << " local errors, "
         << max_error_count << " max errors on proc " << rank_of_max
         << ":\n"
         << connector.format(dbgbord, 2)
         << "base mapped_box_level:\n" << connector.getBase().format(dbgbord, 2)
         << "head mapped_box_level:\n" << connector.getHead().format(dbgbord, 2));
   }
}

/*
 ***********************************************************************
 * Return number of missing and number of extra overlaps.
 ***********************************************************************
 */

int
OverlapConnectorAlgorithm::checkOverlapCorrectness(
   const Connector& connector,
   bool ignore_self_overlap,
   bool assert_completeness,
   bool ignore_periodic_images) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (!connector.getBase().isInitialized() ||
       !connector.getHead().isInitialized()) {
      TBOX_ERROR(
         "OverlapConnectorAlgorithm::findOverlapErrors: Cannot check overlaps when\n"
         << "base or head mapped_box_level is uninitialized.");

   }
#endif
   TBOX_ASSERT(!connector.hasPeriodicLocalNeighborhoodBaseBoxes());

   Connector missing, extra;
   findOverlapErrors(connector, missing, extra, ignore_self_overlap);

   if (!assert_completeness) {
      // Disregard missing overlaps by resetting missing to empty.
      missing.clearNeighborhoods();
   } else if (ignore_periodic_images) {
      // Disregard missing overlaps if they are incident on a periodic box.
      missing.removePeriodicLocalNeighbors();
      missing.eraseEmptyNeighborSets();
   }

   const BoxId dummy_mapped_box_id;

   /*
    * Report the errors found, ordered by the Box where the
    * error appears.  In order to do this, we have to loop through
    * the neighborhoods of missing and extra at the same time.
    */

   Connector::ConstNeighborhoodIterator im = missing.begin();
   Connector::ConstNeighborhoodIterator ie = extra.begin();
   for (; im != missing.end() || ie != extra.end();
        /* incremented in loop */) {

      const BoxId& global_id_missing =
         im == missing.end() ? dummy_mapped_box_id : *im;
      const BoxId& global_id_extra =
         ie == extra.end() ? dummy_mapped_box_id : *ie;

      if (im != missing.end() && ie != extra.end() &&
          *im == *ie) {

         /*
          * im and ie are pointing at the same Box.  Report the
          * error for this Box.
          */

         const Box& mapped_box = *connector.getBase().getBoxStrict(
               global_id_missing);
         tbox::perr << "Found " << missing.numLocalNeighbors(*im)
                    << " missing and "
                    << extra.numLocalNeighbors(*ie)
                    << " extra overlaps for "
                    << mapped_box << std::endl;
         Connector::ConstNeighborhoodIterator it =
            connector.findLocal(global_id_missing);
         if (it == connector.end()) {
            tbox::perr << "  Current Neighbors (no neighbor set)." << std::endl;
         } else {
            tbox::perr << "  Current Neighbors ("
                       << connector.numLocalNeighbors(*it) << "):"
                       << std::endl;
            Box ghost_box = mapped_box;
            ghost_box.grow(connector.getConnectorWidth());
            for (Connector::ConstNeighborIterator na = connector.begin(it);
                 na != connector.end(it); ++na) {
               const Box& nabr = *na;
               Box nabr_box = nabr;
               if (connector.getHeadCoarserFlag()) {
                  nabr_box.refine(connector.getRatio());
               } else if (connector.getRatio() != 1) {
                  nabr_box.coarsen(connector.getRatio());
               }
               Box ovlap = ghost_box * nabr_box;
               tbox::perr << "    " << nabr << '_' << nabr.numberCells()
                          << "\tov" << ovlap << '_' << ovlap.numberCells()
                          << std::endl;
            }
         }
         {
            tbox::perr << "  Missing Neighbors ("
                       << missing.numLocalNeighbors(*im) << "):"
                       << std::endl;
            Box ghost_box = mapped_box;
            ghost_box.grow(connector.getConnectorWidth());
            for (Connector::ConstNeighborIterator na = missing.begin(im);
                 na != missing.end(im); ++na) {
               const Box& nabr = *na;
               Box nabr_box = nabr;
               if (connector.getHeadCoarserFlag()) {
                  nabr_box.refine(connector.getRatio());
               } else if (connector.getRatio() != 1) {
                  nabr_box.coarsen(connector.getRatio());
               }
               Box ovlap = ghost_box * nabr_box;
               tbox::perr << "    " << nabr << '_' << nabr.numberCells()
                          << "\tov" << ovlap << '_' << ovlap.numberCells()
                          << std::endl;
            }
         }
         {
            tbox::perr << "  Extra Neighbors ("
                       << extra.numLocalNeighbors(*ie) << "):"
                       << std::endl;
            Box ghost_box = mapped_box;
            ghost_box.grow(connector.getConnectorWidth());
            for (Connector::ConstNeighborIterator na = extra.begin(ie);
                 na != extra.end(ie); ++na) {
               const Box& nabr = *na;
               Box nabr_box = nabr;
               if (connector.getHeadCoarserFlag()) {
                  nabr_box.refine(connector.getRatio());
               } else if (connector.getRatio() != 1) {
                  nabr_box.coarsen(connector.getRatio());
               }
               Box ovlap = ghost_box * nabr_box;
               tbox::perr << "    " << nabr << '_' << nabr.numberCells()
                          << "\tov" << ovlap << '_' << ovlap.numberCells()
                          << std::endl;
            }
         }
         ++im;
         ++ie;

      } else if ((ie == extra.end()) ||
                 (im != missing.end() && *im < *ie)) {

         /*
          * im goes before ie (or ie has reached the end).  Report the
          * errors for the Box at im.
          */

         const Box& mapped_box = *connector.getBase().getBoxStrict(
               global_id_missing);
         tbox::perr << "Found " << missing.numLocalNeighbors(*im)
                    << " missing overlaps for " << mapped_box << std::endl;
         Connector::ConstNeighborhoodIterator it =
            connector.findLocal(global_id_missing);
         if (it == connector.end()) {
            tbox::perr << "    Current Neighbors (no neighbor set)."
                       << std::endl;
         } else {
            tbox::perr << "  Current Neighbors ("
                       << connector.numLocalNeighbors(*it) << "):"
                       << std::endl;
            Box ghost_box = mapped_box;
            ghost_box.grow(connector.getConnectorWidth());
            for (Connector::ConstNeighborIterator na = connector.begin(it);
                 na != connector.end(it); ++na) {
               const Box& nabr = *na;
               Box nabr_box = nabr;
               if (connector.getHeadCoarserFlag()) {
                  nabr_box.refine(connector.getRatio());
               } else if (connector.getRatio() != 1) {
                  nabr_box.coarsen(connector.getRatio());
               }
               Box ovlap = ghost_box * nabr_box;
               tbox::perr << "    " << nabr << '_' << nabr.numberCells()
                          << "\tov" << ovlap << '_' << ovlap.numberCells()
                          << std::endl;
            }
         }
         {
            tbox::perr << "  Missing Neighbors ("
                       << missing.numLocalNeighbors(*im) << "):"
                       << std::endl;
            Box ghost_box = mapped_box;
            ghost_box.grow(connector.getConnectorWidth());
            for (Connector::ConstNeighborIterator na = missing.begin(im);
                 na != missing.end(im); ++na) {
               const Box& nabr = *na;
               Box nabr_box = nabr;
               if (connector.getHeadCoarserFlag()) {
                  nabr_box.refine(connector.getRatio());
               } else if (connector.getRatio() != 1) {
                  nabr_box.coarsen(connector.getRatio());
               }
               Box ovlap = ghost_box * nabr_box;
               tbox::perr << "    " << nabr << '_' << nabr.numberCells()
                          << "\tov" << ovlap << '_' << ovlap.numberCells()
                          << std::endl;
            }
         }
         ++im;
      } else if ((im == missing.end()) ||
                 (ie != extra.end() && *ie < *im)) {

         /*
          * ie goes before im (or im has reached the end).  Report the
          * errors for the Box at ie.
          */

         const Box& mapped_box = *connector.getBase().getBoxStrict(
               global_id_extra);
         tbox::perr << "Found " << extra.numLocalNeighbors(*ie)
                    << " extra overlaps for " << mapped_box << std::endl;
         Connector::ConstNeighborhoodIterator it =
            connector.findLocal(global_id_extra);
         if (it == connector.end()) {
            tbox::perr << "  Current Neighbors (no neighbor set)." << std::endl;
         } else {
            tbox::perr << "  Current Neighbors ("
                       << connector.numLocalNeighbors(*it) << "):"
                       << std::endl;
            Box ghost_box = mapped_box;
            ghost_box.grow(connector.getConnectorWidth());
            for (Connector::ConstNeighborIterator na = connector.begin(it);
                 na != connector.end(it); ++na) {
               const Box& nabr = *na;
               Box nabr_box = nabr;
               if (connector.getHeadCoarserFlag()) {
                  nabr_box.refine(connector.getRatio());
               } else if (connector.getRatio() != 1) {
                  nabr_box.coarsen(connector.getRatio());
               }
               Box ovlap = ghost_box * nabr_box;
               tbox::perr << "    " << nabr << '_' << nabr.numberCells()
                          << "\tov" << ovlap << '_' << ovlap.numberCells()
                          << std::endl;
            }
         }
         {
            tbox::perr << "  Extra Neighbors ("
                       << extra.numLocalNeighbors(*ie) << "):"
                       << std::endl;
            Box ghost_box = mapped_box;
            ghost_box.grow(connector.getConnectorWidth());
            for (Connector::ConstNeighborIterator na = extra.begin(ie);
                 na != extra.end(ie); ++na) {
               const Box& nabr = *na;
               Box nabr_box = nabr;
               if (connector.getHeadCoarserFlag()) {
                  nabr_box.refine(connector.getRatio());
               } else if (connector.getRatio() != 1) {
                  nabr_box.coarsen(connector.getRatio());
               }
               Box ovlap = ghost_box * nabr_box;
               tbox::perr << "    " << nabr << '_' << nabr.numberCells()
                          << "\tov" << ovlap << '_' << ovlap.numberCells()
                          << std::endl;
            }
         }
         ++ie;
      }

   }

   return missing.getLocalNumberOfNeighborSets() +
          extra.getLocalNumberOfNeighborSets();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
OverlapConnectorAlgorithm::initializeCallback()
{

   if (s_print_steps == '\0') {
      if (tbox::InputManager::inputDatabaseExists()) {
         s_print_steps = 'n';
         boost::shared_ptr<tbox::Database> idb(
            tbox::InputManager::getInputDatabase());
         if (idb->isDatabase("OverlapConnectorAlgorithm")) {
            boost::shared_ptr<tbox::Database> ocu_db(
               idb->getDatabase("OverlapConnectorAlgorithm"));
            s_print_steps = ocu_db->getCharWithDefault("print_bridge_steps",
                  s_print_steps);
         }
      }
   }

   t_find_overlaps_rbbt = tbox::TimerManager::getManager()->
      getTimer("hier::OverlapConnectorAlgorithm::findOverlaps_rbbt()");
   t_bridge = tbox::TimerManager::getManager()->
      getTimer("hier::OverlapConnectorAlgorithm::privateBridge()");
   t_bridge_setup_comm = tbox::TimerManager::getManager()->
      getTimer("hier::OverlapConnectorAlgorithm::setupCommunication()");
   t_bridge_remove_and_cache = tbox::TimerManager::getManager()->
      getTimer("hier::OverlapConnectorAlgorithm::privateBridge_removeAndCache()");
   t_bridge_discover = tbox::TimerManager::getManager()->
      getTimer("hier::OverlapConnectorAlgorithm::privateBridge()_discover");
   t_bridge_discover_get_neighbors = tbox::TimerManager::getManager()->
      getTimer(
         "hier::OverlapConnectorAlgorithm::privateBridge()_discover_get_neighbors");
   t_bridge_discover_form_rbbt = tbox::TimerManager::getManager()->
      getTimer(
         "hier::OverlapConnectorAlgorithm::privateBridge()_discover_form_rbbt");
   t_bridge_discover_find_overlaps = tbox::TimerManager::getManager()->
      getTimer(
         "hier::OverlapConnectorAlgorithm::privateBridge()_discover_find_overlaps");
   t_bridge_share = tbox::TimerManager::getManager()->
      getTimer("hier::OverlapConnectorAlgorithm::privateBridge()_share");
   t_bridge_receive_and_unpack = tbox::TimerManager::getManager()->
      getTimer("hier::OverlapConnectorAlgorithm::receiveAndUnpack");
   t_bridge_MPI_wait = tbox::TimerManager::getManager()->
      getTimer("hier::OverlapConnectorAlgorithm::privateBridge()_MPI_wait");
}

/*
 ***************************************************************************
 * Release static timers.  To be called by shutdown registry to make sure
 * memory for timers does not leak.
 ***************************************************************************
 */

void
OverlapConnectorAlgorithm::finalizeCallback()
{
   t_find_overlaps_rbbt.reset();
   t_bridge.reset();
   t_bridge_setup_comm.reset();
   t_bridge_remove_and_cache.reset();
   t_bridge_discover.reset();
   t_bridge_discover_get_neighbors.reset();
   t_bridge_discover_form_rbbt.reset();
   t_bridge_discover_find_overlaps.reset();
   t_bridge_share.reset();
   t_bridge_receive_and_unpack.reset();
   t_bridge_MPI_wait.reset();

   if (s_class_mpi.getCommunicator() != tbox::SAMRAI_MPI::commNull) {
      s_class_mpi.freeCommunicator();
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
