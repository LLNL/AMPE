/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Set of edges incident from a mapped_box_level of a distributed box graph.
 *
 ************************************************************************/
#ifndef included_hier_Connector_C
#define included_hier_Connector_C

#include "SAMRAI/hier/Connector.h"

#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/ConnectorStatistics.h"
#include "SAMRAI/hier/PeriodicShiftCatalog.h"
#include "SAMRAI/tbox/StartupShutdownManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <algorithm>
//#include <iomanip>

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

static const std::string dbgbord;

namespace SAMRAI {
namespace hier {

const int Connector::HIER_CONNECTOR_VERSION = 0;

boost::shared_ptr<tbox::Timer> Connector::t_acquire_remote_relationships;
boost::shared_ptr<tbox::Timer> Connector::t_cache_global_reduced_data;

tbox::StartupShutdownManager::Handler
Connector::s_initialize_finalize_handler(
   Connector::initializeCallback,
   0,
   0,
   Connector::finalizeCallback,
   tbox::StartupShutdownManager::priorityTimers);

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector::Connector():
   d_base_handle(),
   d_head_handle(),
   d_base_width(tbox::Dimension::getInvalidDimension()),
   d_ratio(tbox::Dimension::getInvalidDimension()),
   d_ratio_is_exact(false),
   d_head_coarser(false),
   d_relationships(),
   d_global_relationships(),
   d_mpi(tbox::SAMRAI_MPI::commNull),
   d_parallel_state(BoxLevel::DISTRIBUTED),
   d_finalized(false),
   d_global_number_of_neighbor_sets(0),
   d_global_number_of_relationships(0),
   d_global_data_up_to_date(false),
   d_connector_type(UNKNOWN)
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector::Connector(
   const Connector& other):
   d_base_handle(other.d_base_handle),
   d_head_handle(other.d_head_handle),
   d_base_width(other.d_base_width),
   d_ratio(other.d_ratio),
   d_ratio_is_exact(other.d_ratio_is_exact),
   d_head_coarser(other.d_head_coarser),
   d_relationships(other.d_relationships),
   d_global_relationships(other.d_global_relationships),
   d_mpi(other.d_mpi),
   d_parallel_state(other.d_parallel_state),
   d_finalized(other.d_finalized),
   d_global_number_of_neighbor_sets(other.d_global_number_of_neighbor_sets),
   d_global_number_of_relationships(other.d_global_number_of_relationships),
   d_global_data_up_to_date(other.d_global_data_up_to_date),
   d_connector_type(other.d_connector_type)
{
}

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector::Connector(
   const BoxLevel& base_mapped_box_level,
   const BoxLevel& head_mapped_box_level,
   const IntVector& base_width,
   const BoxLevel::ParallelState parallel_state):
   d_base_width(base_width.getDim(), 0),
   d_ratio(base_width.getDim(), 0),
   d_head_coarser(false),
   d_relationships(),
   d_global_relationships(),
   d_mpi(base_mapped_box_level.getMPI()),
   d_parallel_state(parallel_state),
   d_finalized(false),
   d_global_number_of_neighbor_sets(0),
   d_global_number_of_relationships(0),
   d_global_data_up_to_date(true),
   d_connector_type(UNKNOWN)
{
   TBOX_DIM_ASSERT_CHECK_ARGS3(base_mapped_box_level,
      head_mapped_box_level,
      base_width);

   setBase(base_mapped_box_level);
   setHead(head_mapped_box_level);
   setWidth(base_width, true);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector::~Connector()
{
   clear();
}

/*
 ***********************************************************************
 ***********************************************************************
 */
const Connector&
Connector::operator = (
   const Connector& rhs)
{
   if (this != &rhs) {
      d_base_handle = rhs.d_base_handle;
      d_global_data_up_to_date = rhs.d_global_data_up_to_date;
      d_global_number_of_neighbor_sets = rhs.d_global_number_of_neighbor_sets;
      d_head_handle = rhs.d_head_handle;
      d_global_number_of_relationships = rhs.d_global_number_of_relationships;
      d_relationships = rhs.d_relationships;
      d_global_relationships = rhs.d_global_relationships;
      d_mpi = rhs.d_mpi;
      d_base_width = rhs.d_base_width;
      d_ratio = rhs.d_ratio;
      d_head_coarser = rhs.d_head_coarser;
      d_parallel_state = rhs.d_parallel_state;
      d_finalized = rhs.d_finalized;
      d_connector_type = rhs.d_connector_type;
   }
   return *this;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
bool
Connector::operator == (
   const Connector& rhs) const
{
   if (this == &rhs) {
      return true;
   }
   // Note: two unfinalized Connectors always compare equal.
   if (!isFinalized() && !rhs.isFinalized()) {
      return true;
   }
   if (!isFinalized() && rhs.isFinalized()) {
      return false;
   }
   if (isFinalized() && !rhs.isFinalized()) {
      return false;
   }

   // Compare only independent attributes.
   if (d_base_width != rhs.d_base_width) {
      return false;
   }
   if (d_base_handle->getBoxLevel() !=
       rhs.d_base_handle->getBoxLevel()) {
      return false;
   }
   if (d_head_handle->getBoxLevel() !=
       rhs.d_head_handle->getBoxLevel()) {
      return false;
   }
   if (d_relationships != rhs.d_relationships) {
      return false;
   }

   return true;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
bool
Connector::operator != (
   const Connector& rhs) const
{
   if (this == &rhs) {
      return false;
   }
   // Note: two unfinalized Connectors always compare equal.
   if (!isFinalized() && !rhs.isFinalized()) {
      return false;
   }
   if (!isFinalized() && rhs.isFinalized()) {
      return true;
   }
   if (isFinalized() && !rhs.isFinalized()) {
      return true;
   }

   // Compare only independent attributes.
   if (d_base_width != rhs.d_base_width) {
      return true;
   }
   if (d_base_handle->getBoxLevel() !=
       rhs.d_base_handle->getBoxLevel()) {
      return true;
   }
   if (d_head_handle->getBoxLevel() !=
       rhs.d_head_handle->getBoxLevel()) {
      return true;
   }
   if (d_relationships != rhs.d_relationships) {
      return true;
   }

   return false;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::insertNeighbors(
   const BoxContainer& neighbors,
   const BoxId& base_box)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state == BoxLevel::DISTRIBUTED &&
       base_box.getOwnerRank() != getBase().getMPI().getRank()) {
      TBOX_ERROR("Connector::insertNeighbors error: Cannot work on remote\n"
         << "data in DISTRIBUTED mode.");
   }
   if (!getBase().hasBox(base_box)) {
      TBOX_ERROR(
         "Exiting due to above reported error."
         << "Connector::insertNeighbors: Cannot access neighbors for\n"
         << "id " << base_box << " because it does not "
         << "exist in the base.\n"
         << "base:\n" << getBase().format("", 2));
   }
#endif
   if (d_parallel_state == BoxLevel::GLOBALIZED) {
      d_global_relationships.insert(base_box, neighbors);
   }
   if (base_box.getOwnerRank() == getBase().getMPI().getRank()) {
      d_relationships.insert(base_box, neighbors);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::eraseNeighbor(
   const Box& neighbor,
   const BoxId& mapped_box_id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_parallel_state == BoxLevel::DISTRIBUTED &&
       mapped_box_id.getOwnerRank() != getBase().getMPI().getRank()) {
      TBOX_ERROR("Connector::eraseNeighbor error: Cannot work on remote\n"
         << "data in DISTRIBUTED mode.");
   }
   if (!getBase().hasBox(mapped_box_id)) {
      TBOX_ERROR(
         "Connector::eraseNeighbors: Cannot access neighbors for\n"
         << "id " << mapped_box_id << " because it does not "
         << "exist in the base.\n"
         << "base:\n" << getBase().format("", 2));
   }
#endif
   if (d_parallel_state == BoxLevel::GLOBALIZED) {
      d_global_relationships.erase(mapped_box_id, neighbor);
   }
   if (mapped_box_id.getOwnerRank() == getBase().getMPI().getRank()) {
      d_relationships.erase(mapped_box_id, neighbor);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::shrinkWidth(
const IntVector& new_width)
{
   if (!(new_width <= getConnectorWidth())) {
      TBOX_ERROR("Connector::shrinkWidth: new ghost cell\n"
         << "width " << new_width << " involves an\n"
         << "enlargement of the current cell width "
         << getConnectorWidth());
   }
   else if (new_width == getConnectorWidth()) {
      // This is a no-op.
      return;
   }

   // Have not yet written this for GLOBALIZED mode.
   TBOX_ASSERT(getParallelState() == BoxLevel::DISTRIBUTED);

   /*
    * Remove overlaps that disappeared given the new GCW.
    * Swap out the overlaps, modify them then swap them back in.
    */

   const bool head_coarser = getHeadCoarserFlag();
   const bool base_coarser = !getHeadCoarserFlag() &&
      getBase().getRefinementRatio() != getHead().getRefinementRatio();

   const boost::shared_ptr<const BaseGridGeometry>& grid_geom(
      getBase().getGridGeometry());

   for (NeighborhoodIterator ei = begin(); ei != end(); ++ei) {
      const BoxId& mapped_box_id = *ei;
      const Box& mapped_box = *getBase().getBoxStrict(
            mapped_box_id);
      Box mapped_box_box = mapped_box;
      mapped_box_box.grow(new_width);
      if (base_coarser) {
         mapped_box_box.refine(getRatio());
      }
      for (NeighborIterator na = begin(ei);
           na != end(ei); /* incremented in loop */) {
         const Box& nabr = *na;
         Box nabr_box = nabr;
         if (nabr.getBlockId() != mapped_box.getBlockId()) {
            grid_geom->transformBox(nabr_box,
               getHead().getRefinementRatio(),
               mapped_box.getBlockId(),
               nabr.getBlockId());
         }
         if (head_coarser) {
            nabr_box.refine(getRatio());
         }
         ++na;
         if (!mapped_box_box.intersects(nabr_box)) {
            d_relationships.erase(ei, nabr);
         }
      }
   }

   d_base_width = new_width;
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::acquireRemoteNeighborhoods()
{
   tbox::SAMRAI_MPI mpi(getMPI());
   if (mpi.getSize() == 1) {
      // In single-proc mode, we already have all the relationships already.
      d_global_relationships = d_relationships;
      return;
   }

   t_acquire_remote_relationships->start();

   std::vector<int> send_mesg;
   std::vector<int> recv_mesg;
   /*
    * Pack relationships from all mapped_box_level relationship sets into a single message.
    * Note that each mapped_box_level relationship set object packs the size of its
    * sub-message into send_mesg.
    */
   acquireRemoteNeighborhoods_pack(send_mesg);
   int send_mesg_size = static_cast<int>(send_mesg.size());

   /*
    * Send and receive the data.
    */

   std::vector<int> recv_mesg_size(getMPI().getSize());
   mpi.Allgather(&send_mesg_size,
      1,
      MPI_INT,
      &recv_mesg_size[0],
      1,
      MPI_INT);

   std::vector<int> proc_offset(getMPI().getSize());
   int totl_size = 0;
   for (int n = 0; n < getMPI().getSize(); ++n) {
      proc_offset[n] = totl_size;
      totl_size += recv_mesg_size[n];
   }
   recv_mesg.resize(totl_size, BAD_INT);
   mpi.Allgatherv(&send_mesg[0],
      send_mesg_size,
      MPI_INT,
      &recv_mesg[0],
      &recv_mesg_size[0],
      &proc_offset[0],
      MPI_INT);

   /*
    * Extract relationship info received from other processors.
    */
   acquireRemoteNeighborhoods_unpack(recv_mesg, proc_offset);

   t_acquire_remote_relationships->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::acquireRemoteNeighborhoods_pack(
   std::vector<int>& send_mesg) const
{
   const tbox::Dimension& dim = getBase().getDim();
   d_relationships.putToIntBuffer(send_mesg, dim, BAD_INT);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::acquireRemoteNeighborhoods_unpack(
   const std::vector<int>& recv_mesg,
   const std::vector<int>& proc_offset)
{
   const tbox::Dimension& dim = getBase().getDim();
   const int num_procs = getMPI().getSize();
   const int rank = getMPI().getRank();
   d_global_relationships = d_relationships;
   d_global_relationships.getFromIntBuffer(
      recv_mesg,
      proc_offset,
      dim,
      num_procs,
      rank);
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::setParallelState(
   const BoxLevel::ParallelState parallel_state)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (!isFinalized()) {
      TBOX_ERROR(
         "Connector::setParallelState: Cannot change the parallel state of\n"
         << "an unfinalized Connector.  See Connector::finalizeContext()");
   }
#endif
   if (parallel_state != BoxLevel::DISTRIBUTED && parallel_state !=
       BoxLevel::GLOBALIZED) {
      TBOX_ERROR("Connector::setParallelState: Invalid distribution state: "
         << parallel_state << "\n");
   }

   if (d_parallel_state == BoxLevel::DISTRIBUTED && parallel_state ==
       BoxLevel::GLOBALIZED) {
      acquireRemoteNeighborhoods();
   } else if (d_parallel_state == BoxLevel::GLOBALIZED && parallel_state ==
              BoxLevel::DISTRIBUTED) {
      d_global_relationships.clear();
   }
   d_parallel_state = parallel_state;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::finalizeContext()
{
   TBOX_ASSERT(d_base_handle);
   TBOX_ASSERT(d_head_handle);
   TBOX_ASSERT(d_base_width.getDim().isValid());

   const BoxLevel& base = d_base_handle->getBoxLevel();
   const BoxLevel& head = d_head_handle->getBoxLevel();
   const IntVector& baseRefinementRatio = base.getRefinementRatio();
   const IntVector& headRefinementRatio = head.getRefinementRatio();

   if (base.getGridGeometry() != head.getGridGeometry()) {
      TBOX_ERROR("Connector::finalizeContext():\n"
         << "Connector must be finalized with\n"
         << "BoxLevels using the same grid geometry.");
   }
   if (!(baseRefinementRatio >= headRefinementRatio ||
         baseRefinementRatio <= headRefinementRatio)) {
      TBOX_ERROR("Connector::finalizeContext():\n"
         << "Refinement ratio between base and head box_levels\n"
         << "cannot be mixed (bigger in some dimension and\n"
         << "smaller in others).\n"
         << "Input base ratio = " << baseRefinementRatio
         << "\n"
         << "Input head ratio = " << headRefinementRatio
         << "\n");
   }
   if (d_parallel_state == BoxLevel::GLOBALIZED &&
       base.getParallelState() != BoxLevel::GLOBALIZED) {
      TBOX_ERROR(
         "Connector::finalizeContext: base BoxLevel must be in\n"
         << "GLOBALIZED state for the Connector to be in\n"
         << "GLOBALIZED state.");
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   bool errf = false;
   for (NeighborhoodIterator ci = begin(); ci != end(); ++ci) {
      if (!base.hasBox(*ci)) {
         tbox::perr << "\nConnector::finalizeContext: NeighborhoodSet "
                    << "provided for non-existent box " << *ci
                    << "\n" << "Neighbors ("
                    << d_relationships.numNeighbors(*ci) << "):\n";
         for (NeighborIterator na = begin(ci); na != end(ci); ++na) {
            tbox::perr << (*na) << "\n";
         }
         errf = true;
      }
   }
   if (errf) {
      TBOX_ERROR(
         "Exiting due to errors."
         << "\nConnector::finalizeContext base box_level:\n"
         << base.format(dbgbord, 2));
   }
#endif
   computeRatioInfo(
      baseRefinementRatio,
      headRefinementRatio,
      d_ratio,
      d_head_coarser,
      d_ratio_is_exact);

   if (d_parallel_state == BoxLevel::DISTRIBUTED) {
      d_global_relationships.clear();
   }
   else {
      if (&d_relationships != &d_global_relationships) {
         d_global_relationships = d_relationships;
      }
   }

   // Erase remote relationships, if any, from d_relationships.
   // Note that we don't call getMPI here to get the rank as d_finalized isn't
   // set until after this step is complete.  It may be picky to insist that
   // d_finalized be set at the very end of the method but it's more correct.
   d_relationships.eraseNonLocalNeighborhoods(
      d_base_handle->getBoxLevel().getMPI().getRank());

   d_finalized = true;
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::setBase(
   const BoxLevel& new_base,
   bool finalize_context)
{
   if (!new_base.isInitialized()) {
      TBOX_ERROR("Connector::setBase():\n"
         << "Connector may not be finalized with\n"
         << "an uninitialized BoxLevel.");
   }
   d_finalized = false;
   d_base_handle = new_base.getBoxLevelHandle();
   d_mpi = new_base.getMPI();
   if (finalize_context) {
      finalizeContext();
   }
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::setHead(
   const BoxLevel& new_head,
   bool finalize_context)
{
   if (!new_head.isInitialized()) {
      TBOX_ERROR("Connector::setHead():\n"
         << "Connector may not be finalized with\n"
         << "an uninitialized BoxLevel.");
   }
   d_finalized = false;
   d_head_handle = new_head.getBoxLevelHandle();
   if (finalize_context) {
      finalizeContext();
   }
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::setWidth(
   const IntVector& new_width,
   bool finalize_context)
{
   if (!(new_width >= IntVector::getZero(new_width.getDim()))) {
      TBOX_ERROR("Connector::setWidth():\n"
         << "Invalid ghost cell width: "
         << new_width << "\n");
   }
   d_finalized = false;
   d_base_width = new_width;
   if (finalize_context) {
      finalizeContext();
   }
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::computeRatioInfo(
   const IntVector& baseRefinementRatio,
   const IntVector& headRefinementRatio,
   IntVector& ratio,
   bool& head_coarser,
   bool& ratio_is_exact)
{
   if (baseRefinementRatio <= headRefinementRatio) {
      ratio = headRefinementRatio / baseRefinementRatio;
      head_coarser = false;
      ratio_is_exact = (ratio * baseRefinementRatio) == headRefinementRatio;
   }
   else {
      ratio = baseRefinementRatio / headRefinementRatio;
      head_coarser = true;
      ratio_is_exact = (ratio * headRefinementRatio) == baseRefinementRatio;
   }
   if (baseRefinementRatio * headRefinementRatio <
       IntVector::getZero(baseRefinementRatio.getDim())) {
      // Note that negative ratios like -N really mean 1/N (negative reciprocal).
      ratio = -headRefinementRatio * baseRefinementRatio;
      ratio_is_exact = true;
   }
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::writeNeighborhoodsToErrorStream(
   const std::string& border) const
{
   const std::string indented_border = border + "  ";
   tbox::perr << border << "  " << d_relationships.numBoxNeighborhoods()
              << " neigborhoods:\n";
   for (ConstNeighborhoodIterator ei = begin(); ei != end(); ++ei) {
      tbox::perr << border << "  " << *ei << "\n";
      tbox::perr << border << "    Neighbors ("
                 << d_relationships.numNeighbors(ei) << "):\n";
      for (ConstNeighborIterator bi = begin(ei); bi != end(ei); ++bi) {
         const Box& box = *bi;
         tbox::perr << border << "    "
                    << box << "   "
                    << box.numberCells() << '\n';
      }
   }
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::writeNeighborhoodToErrorStream(
   const BoxId& box_id) const
{
   const BoxNeighborhoodCollection& relationships = getRelations(box_id);
   BoxId non_per_id(box_id.getGlobalId(),
                    PeriodicId::zero());
   ConstNeighborhoodIterator ei = relationships.find(non_per_id);
   if (ei == relationships.end()) {
      TBOX_ERROR("Connector::find: No neighbor set exists for\n"
         << "box " << box_id << ".\n");
   }
   for (ConstNeighborIterator bi = relationships.begin(ei);
        bi != relationships.end(ei); ++bi) {
      const Box& box = *bi;
      tbox::perr << "    "
                 << box << "   "
                 << box.numberCells() << '\n';
   }
   return;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::initializeToLocalTranspose(
   const Connector& connector)
{
   const IntVector my_gcw = convertHeadWidthToBase(
         connector.getHead().getRefinementRatio(),
         connector.getBase().getRefinementRatio(),
         connector.getConnectorWidth());

   clearNeighborhoods();
   setBase(connector.d_head_handle->getBoxLevel());
   setHead(connector.d_base_handle->getBoxLevel());
   setWidth(my_gcw, true);
   TBOX_ASSERT(isTransposeOf(connector));

   const tbox::Dimension dim(my_gcw.getDim());
   const PeriodicShiftCatalog* shift_catalog =
      PeriodicShiftCatalog::getCatalog(dim);

   for (ConstNeighborhoodIterator ci = connector.begin();
        ci != connector.end(); ++ci) {

      const BoxId& mapped_box_id = *ci;
      const BoxContainer::const_iterator ni = getHead().getBox(mapped_box_id);
      if (ni == getHead().getBoxes().end()) {
         TBOX_ERROR(
            "Connector::initializeToLocalTranspose: mapped_box index\n"
            << mapped_box_id
            << " not found in local part of head mapped_box_level.\n"
            << "This means that the incoming Connector data was not a\n"
            << "self-consistent local mapping.\n");
      }
      const Box& my_head_mapped_box = *ni;

      for (ConstNeighborIterator na = connector.begin(ci);
           na != connector.end(ci); ++na) {
         const Box& my_base_mapped_box = *na;
         if (my_base_mapped_box.getOwnerRank() != getMPI().getRank()) {
            TBOX_ERROR(
               "Connector::initializeToLocalTranspose: base mapped_box "
               << my_head_mapped_box << "\n"
               << "has remote neighbor " << my_base_mapped_box
               << " which is disallowed.\n"
               << "Mapped_boxes must have only local neighbors in this method.");
         }
         if (my_base_mapped_box.isPeriodicImage()) {
            Box my_shifted_head_mapped_box(
               my_head_mapped_box,
               shift_catalog->getOppositeShiftNumber(
                  my_base_mapped_box.getPeriodicId()),
               getHead().getRefinementRatio());
            if (getHead().hasBox(my_shifted_head_mapped_box)) {
               BoxId base_non_per_id(
                  my_base_mapped_box.getGlobalId(),
                  PeriodicId::zero());
               d_relationships.insert(
                  base_non_per_id,
                  my_shifted_head_mapped_box);
            }
         } else {
            d_relationships.insert(
               my_base_mapped_box.getId(),
               my_head_mapped_box);
         }
      }

   }

   if (0) {
      tbox::perr << "end of initializeToLocalTranspose:\n"
                 << "base:\n" << getBase().format("BASE->", 3)
                 << "head:\n" << getHead().format("HEAD->", 3)
                 << "this:\n" << format("THIS->", 3)
                 << "r:\n" << connector.format("RRRR->", 3)
                 << "Checking this transpose correctness:" << std::endl;
      assertTransposeCorrectness(connector, false);
      tbox::perr << "Checking r's transpose correctness:" << std::endl;
      connector.assertTransposeCorrectness(*this, false);
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

bool
Connector::isTransposeOf(
   const Connector& other) const
{
   bool rval = false;
   if (d_base_handle == other.d_head_handle &&
       d_head_handle == other.d_base_handle) {
      if (d_head_coarser) {
         IntVector transpose_base_width = convertHeadWidthToBase(
               getHead().getRefinementRatio(),
               getBase().getRefinementRatio(),
               d_base_width);
         rval = other.d_base_width == transpose_base_width;
      } else {
         IntVector transpose_base_width = convertHeadWidthToBase(
               other.getHead().getRefinementRatio(),
               other.getBase().getRefinementRatio(),
               other.d_base_width);
         rval = d_base_width == transpose_base_width;
      }
   }
   return rval;
}

/*
 ***********************************************************************
 ***********************************************************************
 */
void
Connector::cacheGlobalReducedData() const
{
   TBOX_ASSERT(isFinalized());

   if (d_global_data_up_to_date) {
      return;
   }

   t_cache_global_reduced_data->barrierAndStart();

   tbox::SAMRAI_MPI mpi(getMPI());

   if (d_parallel_state == BoxLevel::GLOBALIZED) {
      d_global_number_of_relationships =
         d_global_relationships.sumNumNeighbors();
      d_global_number_of_neighbor_sets =
         d_global_relationships.numBoxNeighborhoods();
   } else {
      if (mpi.getSize() > 1) {
         int tmpa[2], tmpb[2];
         tmpa[0] = getLocalNumberOfNeighborSets();
         tmpa[1] = getLocalNumberOfRelationships();

         TBOX_ASSERT(tmpa[0] >= 0);
         TBOX_ASSERT(tmpa[0] >= 0);

         mpi.Allreduce(tmpa,
            tmpb,                        // Better to use MPI_IN_PLACE, but not some MPI's do not support.
            2,
            MPI_INT,
            MPI_SUM);
         d_global_number_of_neighbor_sets = tmpb[0];
         d_global_number_of_relationships = tmpb[1];
      } else {
         d_global_number_of_neighbor_sets = getLocalNumberOfNeighborSets();
         d_global_number_of_relationships = getLocalNumberOfRelationships();
      }

      TBOX_ASSERT(d_global_number_of_neighbor_sets >= 0);
      TBOX_ASSERT(d_global_number_of_relationships >= 0);
   }

   d_global_data_up_to_date = true;

   t_cache_global_reduced_data->stop();
}

/*
 ***********************************************************************
 ***********************************************************************
 */

IntVector
Connector::convertHeadWidthToBase(
   const IntVector& base_refinement_ratio,
   const IntVector& head_refinement_ratio,
   const IntVector& head_gcw)
{
   if (!(base_refinement_ratio >= head_refinement_ratio ||
         base_refinement_ratio <= head_refinement_ratio)) {
      TBOX_ERROR("Connector::convertHeadWidthToBase:\n"
         << "head mapped_box_level must be either\n"
         << "finer or coarser than base.\n"
         << "Combined refinement and coarsening not allowed.");
   }

   tbox::Dimension dim(head_refinement_ratio.getDim());

   IntVector ratio(dim); // Ratio between head and base.

   if (head_refinement_ratio * base_refinement_ratio >
       IntVector::getZero(dim)) {
      // Same signs for both ratios -> simple to compute head-base ratio.
      if (base_refinement_ratio >= head_refinement_ratio) {
         ratio = base_refinement_ratio / head_refinement_ratio;
      } else {
         ratio = head_refinement_ratio / base_refinement_ratio;
      }
   } else {
      // Note that negative ratios like -N really mean 1/N (negative reciprocal).
      ratio = -base_refinement_ratio * head_refinement_ratio;
   }
   TBOX_ASSERT(ratio >= IntVector::getOne(dim));

   const IntVector base_width =
      (base_refinement_ratio >= head_refinement_ratio) ?
      (head_gcw * ratio) : IntVector::ceilingDivide(head_gcw, ratio);

   return base_width;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::recursivePrint(
   std::ostream& os,
   const std::string& border,
   int detail_depth) const
{
   if (detail_depth < 0) {
      return;
   }

   if (!isFinalized()) {
      os << border << "Unfinalized.\n";
      return;
   }
   bool head_coarser = d_head_coarser;
   const IntVector head_gcw =
      convertHeadWidthToBase(
         getHead().getRefinementRatio(),
         getBase().getRefinementRatio(),
         d_base_width);
   os << border << "Parallel state     : "
      << (getParallelState() == BoxLevel::DISTRIBUTED ? "DIST" : "GLOB")
      << '\n'
      << border << "Rank,nproc         : " << getMPI().getRank() << ", " << getMPI().getSize() << '\n'
      << border << "Base,head objects  :"
      << " ("
      << (d_base_handle == d_head_handle ? "same" : "different") << ") "
      << (void *)&d_base_handle->getBoxLevel() << ", "
      << (void *)&d_head_handle->getBoxLevel() << "\n"
      << border << "Base,head,/ ratios : "
      << getBase().getRefinementRatio() << ", "
      << getHead().getRefinementRatio() << ", "
      << d_ratio << (d_head_coarser ? " (head coarser)" : "") << '\n'
      << border << "Base,head widths   : " << d_base_width << ", "
      << head_gcw << '\n'
      << border << "Box count    : " << getBase().getLocalNumberOfBoxes()
      << " (" << getLocalNumberOfNeighborSets() << " with neighbor lists)\n"
   ;
   if (detail_depth > 0) {
      os << border << "Mapped_boxes with neighbors:\n";
      for (ConstNeighborhoodIterator ei = begin(); ei != end(); ++ei) {
         const BoxId& box_id = *ei;
         BoxContainer::const_iterator ni = getBase().getBox(box_id);
         if (ni != getBase().getBoxes().end()) {
            os << border << "  "
               << (*ni) << "_"
               << (*ni).numberCells() << '\n';
         } else {
            os << border << "  #"
               << box_id
               << ": INVALID DATA WARNING: No base box with this index!\n";
         }
         os << border << "    Neighbors (" << numLocalNeighbors(box_id) << "):"
            << ((detail_depth > 1) ? "\n" : " ...\n");
         if (detail_depth > 1) {
            for (ConstNeighborIterator i_nabr = begin(ei);
                 i_nabr != end(ei); ++i_nabr) {
               os << border << "      "
                  << (*i_nabr) << "_"
                  << i_nabr->numberCells();
               if (ni != getBase().getBoxes().end()) {
                  Box ovlap = *i_nabr;
                  if (ni->getBlockId() != i_nabr->getBlockId()) {
                     d_base_handle->getBoxLevel().getGridGeometry()->
                        transformBox(
                           ovlap,
                           d_head_handle->getBoxLevel().getRefinementRatio(),
                           ni->getBlockId(),
                           i_nabr->getBlockId());
                  }
                  if (head_coarser) {
                     ovlap.refine(d_ratio);
                  }
                  else if (d_ratio != 1) {
                     ovlap.coarsen(d_ratio);
                  }
                  Box ghost_box = (*ni);
                  ghost_box.grow(d_base_width);
                  ovlap = ovlap * ghost_box;
                  os << "\tov" << ovlap << "_" << ovlap.numberCells();
               }
               os << '\n';
            }
         }
      }
   }
}

/*
 ***********************************************************************
 * Construct a Connector Outputter with formatting parameters.
 ***********************************************************************
 */

Connector::Outputter::Outputter(
   const Connector& connector,
   const std::string& border,
   int detail_depth,
   bool output_statistics):
   d_conn(connector),
   d_border(border),
   d_detail_depth(detail_depth),
   d_output_statistics(output_statistics)
{
}

/*
 ***********************************************************************
 * Print out a Connector according to settings in the Outputter.
 ***********************************************************************
 */

std::ostream&
operator << (
   std::ostream& os,
   const Connector::Outputter& format)
{
   if ( format.d_output_statistics ) {
      ConnectorStatistics cs(format.d_conn);
      cs.printNeighborStats(os, format.d_border);
   }
   else {
      format.d_conn.recursivePrint(os, format.d_border, format.d_detail_depth);
   }
   return os;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

Connector*
Connector::makeGlobalizedCopy(
   const Connector& other) const
{
   // Prevent wasteful accidental use when this method is not needed.
   TBOX_ASSERT(other.getParallelState() != BoxLevel::GLOBALIZED);

   Connector* copy = new Connector(other);
   copy->setParallelState(BoxLevel::GLOBALIZED);
   return copy;
}

/*
 ***********************************************************************
 * Run checkTransposeCorrectness and assert that no errors are found.
 ***********************************************************************
 */

void
Connector::assertTransposeCorrectness(
   const Connector& input_transpose,
   const bool ignore_periodic_relationships) const
{
   size_t err_count =
      checkTransposeCorrectness(input_transpose, ignore_periodic_relationships);
   if (err_count) {
      TBOX_ERROR(
         "Connector::assertTransposeCorrectness:\n"
         << "Aborting with " << err_count << " transpose errors found:\n"
         << "this base:\n" << getBase().format("B:", 3)
         << "this head:\n" << getHead().format("H:", 3)
         << "this Connector:\n" << format("B->H:", 3)
         << "\ntranspose Connector:\n" << input_transpose.format("H->B:", 3));
   }
}

/*
 ***********************************************************************
 *
 * For every relationship in this, there should be reverse relationship in
 * transpose.
 *
 * This method does not check whether the Connectors are defined to
 * form logical transposes (based on their widths and their base and
 * head mapped_box_levels).  For that, see isTransposeOf().
 *
 ***********************************************************************
 */

size_t
Connector::checkTransposeCorrectness(
   const Connector& input_transpose,
   const bool ignore_periodic_relationships) const
{
   const tbox::Dimension dim(getBase().getDim());

   const Connector* transpose =
      (input_transpose.d_parallel_state == BoxLevel::GLOBALIZED) ?
      &input_transpose : makeGlobalizedCopy(input_transpose);

   const BoxLevel& head = getHead().getGlobalizedVersion();

   const PeriodicShiftCatalog* shift_catalog =
      PeriodicShiftCatalog::getCatalog(dim);

   /*
    * Check for extraneous relationships.
    * For every relationship in this, there should be reverse relationship in transpose.
    */
   Box shifted_mapped_box(dim);   // Shifted version of an unshifted Box.
   Box unshifted_mapped_box(dim); // Unhifted version of a shifted Box.

   size_t err_count = 0;


   const BoxNeighborhoodCollection& tran_relationships =
      transpose->getGlobalNeighborhoodSets();
   for (ConstNeighborhoodIterator ci = begin(); ci != end(); ++ci) {

      const BoxId& mapped_box_id = *ci;
      const Box& mapped_box = *getBase().getBox(mapped_box_id);

      size_t err_count_for_current_index = 0;

      for (ConstNeighborIterator ni = begin(ci); ni != end(ci); ++ni) {

         if (ignore_periodic_relationships && ni->isPeriodicImage()) {
            continue;
         }

         const Box& nabr = *ni;

         /*
          * Key for find in NeighborhoodSet must be non-periodic.
          */
         BoxId non_per_nabr_id(nabr.getGlobalId(),
                               PeriodicId::zero());

         ConstNeighborhoodIterator cn =
            tran_relationships.find(non_per_nabr_id);

         if (cn == tran_relationships.end()) {
            tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
            << "Local mapped_box " << mapped_box
            << " has relationship to " << nabr
            << " but " << nabr << " has no relationship container.\n";
            ++err_count_for_current_index;
            continue;
         }

         TBOX_ASSERT(*cn == non_per_nabr_id);

         bool nabr_has_box;
         if (nabr.isPeriodicImage()) {
            shifted_mapped_box.initialize(
               mapped_box,
               shift_catalog->getOppositeShiftNumber(nabr.getPeriodicId()),
               getBase().getRefinementRatio());
            nabr_has_box =
               tran_relationships.hasNeighbor(cn, shifted_mapped_box);
         } else {
            nabr_has_box = tran_relationships.hasNeighbor(cn, mapped_box);
         }

         if (!nabr_has_box) {
            tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
            << "Local mapped_box " << mapped_box;
            if (nabr.isPeriodicImage()) {
               tbox::perr << " (shifted version " << shifted_mapped_box << ")";
            }
            tbox::perr << " has relationship to " << nabr << " but "
            << nabr << " does not have the reverse relationship.\n"
            ;
            tbox::perr << "Neighbors of " << nabr << " are:\n";
            for (ConstNeighborIterator nj = tran_relationships.begin(cn);
                 nj != tran_relationships.end(cn); ++nj) {
               tbox::perr << "   " << *nj << std::endl;
            }
            ++err_count_for_current_index;
            continue;
         }

      }

      if (err_count_for_current_index > 0) {
         tbox::perr << "Mapped_box " << mapped_box << " had "
         << err_count_for_current_index
         << " errors.  Neighbors are:\n";
         for (ConstNeighborIterator nj = begin(ci); nj != end(ci); ++nj) {
            tbox::perr << "  " << *nj << std::endl;
         }
         err_count += err_count_for_current_index;
      }

   }

   /*
    * Check for missing relationships:
    * Transpose should not contain any relationship that does not correspond to
    * one in this.
    */

   for (ConstNeighborhoodIterator ci = tran_relationships.begin();
        ci != tran_relationships.end(); ++ci) {

      const BoxId& mapped_box_id = *ci;

      size_t err_count_for_current_index = 0;

      if (!head.hasBox(mapped_box_id)) {
         TBOX_ASSERT(head.hasBox(mapped_box_id));
      }
      const Box& head_mapped_box = *head.getBoxStrict(mapped_box_id);

      for (ConstNeighborIterator na = tran_relationships.begin(ci);
           na != tran_relationships.end(ci); ++na) {

         const Box nabr = *na;

         if (nabr.getOwnerRank() == getMPI().getRank()) {

            if (ignore_periodic_relationships && nabr.isPeriodicImage()) {
               continue;
            }

            if (!getBase().hasBox(nabr)) {
               tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
               << "Head mapped_box " << head_mapped_box
               << " has neighbor " << nabr << "\n"
               << " but the neighbor does not exist "
               << "in the base mapped_box_level.\n";
               tbox::perr << "Neighbors of head mapped_box "
               << mapped_box_id << " are:\n";
               for (ConstNeighborIterator nj = tran_relationships.begin(ci);
                    nj != tran_relationships.end(ci); ++nj) {
                  tbox::perr << "   " << *nj << std::endl;
               }
               ++err_count_for_current_index;
               continue;
            }

            const Box& base_mapped_box = *getBase().getBoxStrict(nabr);

            /*
             * Non-periodic BoxId needed for NeighborhoodSet::find()
             */
            BoxId base_non_per_id(base_mapped_box.getGlobalId(),
                                  PeriodicId::zero());

            if (!d_relationships.isBaseBox(base_non_per_id)) {
               tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
               << "Head mapped_box " << head_mapped_box << "\n"
               << " has base mapped_box "
               << base_mapped_box << " as a neighbor.\n"
               << "But " << base_mapped_box
               << " has no neighbor container.\n";
               tbox::perr << "Neighbors of head mapped_box " << BoxId(
                  mapped_box_id)
               << ":" << std::endl;
               for (ConstNeighborIterator nj = tran_relationships.begin(ci);
                    nj != tran_relationships.end(ci); ++nj) {
                  tbox::perr << "   " << *nj << std::endl;
               }
               ++err_count_for_current_index;
               continue;
            }

            const Box nabr_nabr(dim, mapped_box_id.getGlobalId(),
                                shift_catalog->getOppositeShiftNumber(
                                   base_mapped_box.getPeriodicId()));

            if (!d_relationships.hasNeighbor(base_non_per_id, nabr_nabr)) {
               tbox::perr << "\nConnector::checkTransposeCorrectness:\n"
               << "Head mapped_box " << head_mapped_box << "\n"
               << " has base mapped_box " << base_mapped_box
               << " as a neighbor.\n"
               << "But base mapped_box " << base_mapped_box
               << " does not have a mapped_box indexed "
               << nabr_nabr.getId()
               << " in its neighbor list." << std::endl;
               tbox::perr << "Neighbors of head mapped_box " << nabr_nabr.getId()
               << ":" << std::endl;
               for (ConstNeighborIterator nj = tran_relationships.begin(ci);
                    nj != tran_relationships.end(ci); ++nj) {
                  tbox::perr << "   " << *nj << std::endl;
               }
               tbox::perr << "Neighbors of base mapped_box ";
               if (nabr.isPeriodicImage()) {
                  unshifted_mapped_box.initialize(
                     nabr,
                     shift_catalog->getZeroShiftNumber(),
                     getBase().getRefinementRatio());
                  tbox::perr << unshifted_mapped_box;
               }
               tbox::perr << ":" << std::endl;
               ConstNeighborhoodIterator nabr_nabrs_ =
                  d_relationships.find(base_non_per_id);
               for (ConstNeighborIterator nj = begin(nabr_nabrs_);
                    nj != end(nabr_nabrs_); ++nj) {
                  tbox::perr << "   " << *nj << std::endl;
               }
               ++err_count_for_current_index;
               continue;
            }

         }

      }

      if (err_count_for_current_index > 0) {
         err_count += err_count_for_current_index;
      }

   }

   if (transpose != &input_transpose) {
      delete transpose;
   }

   int global_err_count = static_cast<int>(err_count);
   getBase().getMPI().AllReduce( &global_err_count, 1, MPI_SUM );

   return static_cast<size_t>(global_err_count);
}

/*
 ***********************************************************************
 ***********************************************************************
 */

size_t
Connector::checkConsistencyWithBase() const
{
   size_t num_errors = 0;
   for (ConstNeighborhoodIterator i_relationships = begin();
        i_relationships != end(); ++i_relationships) {
      const BoxId& mapped_box_id = *i_relationships;
      if (!getBase().hasBox(mapped_box_id)) {
         ++num_errors;
         tbox::plog << "ERROR->"
         << "Connector::assertConsistencyWithBase: Neighbor data given "
         << "\nfor mapped_box " << mapped_box_id
         << " but the mapped_box does not exist.\n";
      }
   }
   return num_errors;
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::assertConsistencyWithBase() const
{
   if (checkConsistencyWithBase() > 0) {
      TBOX_ERROR(
         "Connector::assertConsistencyWithBase() found inconsistencies.\n"
         << "Base mapped box level:\n" << getBase().format("ERROR->", 2));
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::computeNeighborhoodDifferences(
   Connector& left_minus_right,
   const Connector& left,
   const Connector& right)
{
   if (0) {
      tbox::plog << "Computing relationship differences, a:\n" << left.format(dbgbord, 3)
      << "Computing relationship differences, b:\n" << right.format(dbgbord, 3);
   }
   left_minus_right.clearNeighborhoods();
   left_minus_right.d_parallel_state = left.getParallelState();
   left_minus_right.setBase(left.d_base_handle->getBoxLevel());
   left_minus_right.setHead(left.d_head_handle->getBoxLevel());
   left_minus_right.setWidth(left.d_base_width, true);

   for (ConstNeighborhoodIterator ai = left.begin(); ai != left.end(); ++ai) {

      const BoxId& mapped_box_id = *ai;

      ConstNeighborhoodIterator bi = right.findLocal(mapped_box_id);
      if (bi != right.end()) {
         // Remove bi from ai.  Put results in a_minus_b.
         NeighborhoodIterator base_box_itr =
            left_minus_right.makeEmptyLocalNeighborhood(mapped_box_id);

         //equivalent of stl set_difference
         ConstNeighborIterator na = left.begin(ai);
         ConstNeighborIterator nb = right.begin(bi);
         while (na != left.end(ai) && nb != right.end(bi)) {
            if (na->getId() < nb->getId()) {
               left_minus_right.insertLocalNeighbor(*na, base_box_itr);
               ++na;
            } else if (nb->getId() < na->getId()) {
               ++nb;
            } else {
               ++na;
               ++nb;
            } 
         }
             
         if (left_minus_right.isEmptyNeighborhood(mapped_box_id)) {
            left_minus_right.eraseLocalNeighborhood(mapped_box_id);
         }
      } else if (left.numLocalNeighbors(mapped_box_id) != 0) {
         NeighborhoodIterator base_box_itr =
            left_minus_right.makeEmptyLocalNeighborhood(mapped_box_id);
         for (ConstNeighborIterator na = left.begin(ai);
              na != left.end(ai); ++na) {
            left_minus_right.insertLocalNeighbor(*na, base_box_itr);
         }
      }

   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

void
Connector::assertConsistencyWithHead() const
{
   const int number_of_inconsistencies = checkConsistencyWithHead();
   if (number_of_inconsistencies > 0) {
      TBOX_ERROR(
         "Connector::assertConsistencyWithHead() found inconsistencies.\n"
         << getBase().format("base-> ", 3)
         << getHead().format("head-> ", 3)
         << format("E-> ", 3));
   }
}

/*
 ***********************************************************************
 ***********************************************************************
 */

size_t
Connector::checkConsistencyWithHead() const
{
   const BoxLevel& head_mapped_box_level = getHead().getGlobalizedVersion();

   TBOX_ASSERT(head_mapped_box_level.getParallelState() ==
      BoxLevel::GLOBALIZED);

   const BoxContainer& head_mapped_boxes =
      head_mapped_box_level.getGlobalBoxes();

   size_t number_of_inconsistencies = 0;

   /*
    * For each neighbor in each neighbor list,
    * check that the neighbor is in the head_mapped_box_level.
    */

   for (ConstNeighborhoodIterator ei = begin(); ei != end(); ++ei) {

      const BoxId& mapped_box_id = *ei;

      for (ConstNeighborIterator na = begin(ei); na != end(ei); ++na) {

         const Box& nabr = *na;
         const Box unshifted_nabr(
            nabr, PeriodicId::zero(), head_mapped_box_level.getRefinementRatio());

         BoxContainer::const_iterator na_in_head =
            head_mapped_boxes.find(unshifted_nabr);

         if (na_in_head == head_mapped_boxes.end()) {
            tbox::perr << "\nConnector::checkConsistencyWithHead:\n"
            << "Neighbor list for mapped_box " << mapped_box_id << "\n"
            << "referenced nonexistent neighbor "
            << nabr << "\n";
            tbox::perr << "Neighbors of mapped_box " << mapped_box_id << ":\n";
            for (ConstNeighborIterator nb = begin(ei); nb != end(ei); ++nb) {
               tbox::perr << "   " << *nb << '\n';
            }
            ++number_of_inconsistencies;
            continue;
         }

         const Box& nabr_in_head = *na_in_head;
         if (!unshifted_nabr.isIdEqual(nabr_in_head) ||
             !unshifted_nabr.isSpatiallyEqual(nabr_in_head)) {
            tbox::perr << "\nConnector::checkConsistencyWithHead:\n"
            << "Inconsistent mapped_box data at mapped_box "
            << mapped_box_id << "\n"
            << "Neighbor " << nabr << "(unshifted to "
            << unshifted_nabr << ") does not match "
            << "head mapped_box " << nabr_in_head
            << "\n";
            ++number_of_inconsistencies;
         }

      }
   }

   return number_of_inconsistencies;
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
