/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Operations for complex edge data on multiple levels.
 *
 ************************************************************************/

#ifndef included_math_HierarchyEdgeDataOpsComplex_C
#define included_math_HierarchyEdgeDataOpsComplex_C

#include "SAMRAI/math/HierarchyEdgeDataOpsComplex.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/BoxUtilities.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/pdat/EdgeDataFactory.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include <typeinfo>
#include <stdlib.h>
#include <float.h>
#include <math.h>

namespace SAMRAI {
namespace math {

HierarchyEdgeDataOpsComplex::HierarchyEdgeDataOpsComplex(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int coarsest_level,
   const int finest_level):
   HierarchyDataOpsComplex(),
   d_hierarchy(hierarchy)
{
   TBOX_ASSERT(hierarchy);

   if ((coarsest_level < 0) || (finest_level < 0)) {
      if (d_hierarchy->getNumberOfLevels() == 0) {
         d_coarsest_level = coarsest_level;
         d_finest_level = finest_level;
      } else {
         resetLevels(0, d_hierarchy->getFinestLevelNumber());
      }
   } else {
      resetLevels(coarsest_level, finest_level);
   }
}

HierarchyEdgeDataOpsComplex::~HierarchyEdgeDataOpsComplex()
{
}

/*
 *************************************************************************
 *
 * Routines to set the hierarchy and level information.
 *
 *************************************************************************
 */

void
HierarchyEdgeDataOpsComplex::setPatchHierarchy(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
   TBOX_ASSERT(hierarchy);

   d_hierarchy = hierarchy;
}

void
HierarchyEdgeDataOpsComplex::resetLevels(
   const int coarsest_level,
   const int finest_level)
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((coarsest_level >= 0)
      && (finest_level >= coarsest_level)
      && (finest_level <= d_hierarchy->getFinestLevelNumber()));

   int dimVal = d_hierarchy->getDim().getValue();

   d_coarsest_level = coarsest_level;
   d_finest_level = finest_level;

   for (int d = 0; d < dimVal; d++) {
      d_nonoverlapping_edge_boxes[d].resizeArray(d_finest_level + 1);
   }

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      hier::BoxContainer edge_boxes;

      for (int nd = 0; nd < dimVal; nd++) {
         edge_boxes = level->getBoxes();
         for (hier::BoxContainer::iterator i(edge_boxes);
              i != edge_boxes.end(); ++i) {
            *i = pdat::EdgeGeometry::toEdgeBox(*i, nd);
         }
         hier::BoxUtilities::makeNonOverlappingBoxContainers(
            d_nonoverlapping_edge_boxes[nd][ln],
            edge_boxes);
      }
   }
}

const boost::shared_ptr<hier::PatchHierarchy>
HierarchyEdgeDataOpsComplex::getPatchHierarchy() const
{
   return d_hierarchy;
}

/*
 *************************************************************************
 *
 * Basic generic operations.
 *
 *************************************************************************
 */

void
HierarchyEdgeDataOpsComplex::copyData(
   const int dst_id,
   const int src_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > s(
            p->getPatchData(src_id),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.copyData(d, s, box);
      }
   }
}

void
HierarchyEdgeDataOpsComplex::swapData(
   const int data1_id,
   const int data2_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   boost::shared_ptr<pdat::EdgeDataFactory<dcomplex> > d1fact(
      d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data1_id),
      boost::detail::dynamic_cast_tag());
   TBOX_ASSERT(d1fact);
   boost::shared_ptr<pdat::EdgeDataFactory<dcomplex> > d2fact(
      d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data2_id),
      boost::detail::dynamic_cast_tag());
   TBOX_ASSERT(d2fact);
   TBOX_ASSERT(d1fact->getDepth() == d2fact->getDepth());
   TBOX_ASSERT(d1fact->getGhostCellWidth() == d2fact->getGhostCellWidth());
#endif

   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         d_patch_ops.swapData(p, data1_id, data2_id);
      }
   }
}

void
HierarchyEdgeDataOpsComplex::printData(
   const int data_id,
   std::ostream& s,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   s << "Patch descriptor id = " << data_id << std::endl;
   s << "Factory = " << typeid(*d_hierarchy->getPatchDescriptor()->
                               getPatchDataFactory(data_id)).name()
     << std::endl;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      s << "Level number = " << ln << std::endl;
       boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > d(
            p->getPatchData(data_id),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.printData(d, box, s);
      }
   }
}

void
HierarchyEdgeDataOpsComplex::setToScalar(
   const int data_id,
   const dcomplex& alpha,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > d(
            p->getPatchData(data_id),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.setToScalar(d, alpha, box);
      }
   }
}

/*
 *************************************************************************
 *
 * Basic generic arithmetic operations.
 *
 *************************************************************************
 */

void
HierarchyEdgeDataOpsComplex::scale(
   const int dst_id,
   const dcomplex& alpha,
   const int src_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > dst(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > src(
            p->getPatchData(src_id),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(dst);

         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.scale(dst, alpha, src, box);
      }
   }
}

void
HierarchyEdgeDataOpsComplex::addScalar(
   const int dst_id,
   const int src_id,
   const dcomplex& alpha,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > dst(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > src(
            p->getPatchData(src_id),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(dst);

         hier::Box box = (interior_only ? p->getBox() : dst->getGhostBox());

         d_patch_ops.addScalar(dst, src, alpha, box);
      }
   }
}

void
HierarchyEdgeDataOpsComplex::add(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > s1(
            p->getPatchData(src1_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > s2(
            p->getPatchData(src2_id),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.add(d, s1, s2, box);
      }
   }
}

void
HierarchyEdgeDataOpsComplex::subtract(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > s1(
            p->getPatchData(src1_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > s2(
            p->getPatchData(src2_id),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.subtract(d, s1, s2, box);
      }
   }
}

void
HierarchyEdgeDataOpsComplex::multiply(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > s1(
            p->getPatchData(src1_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > s2(
            p->getPatchData(src2_id),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.multiply(d, s1, s2, box);
      }
   }
}

void
HierarchyEdgeDataOpsComplex::divide(
   const int dst_id,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > s1(
            p->getPatchData(src1_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > s2(
            p->getPatchData(src2_id),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.divide(d, s1, s2, box);
      }
   }
}

void
HierarchyEdgeDataOpsComplex::reciprocal(
   const int dst_id,
   const int src_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > src(
            p->getPatchData(src_id),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.reciprocal(d, src, box);
      }
   }
}

void
HierarchyEdgeDataOpsComplex::linearSum(
   const int dst_id,
   const dcomplex& alpha,
   const int src1_id,
   const dcomplex& beta,
   const int src2_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > s1(
            p->getPatchData(src1_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > s2(
            p->getPatchData(src2_id),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.linearSum(d, alpha, s1, beta, s2, box);
      }
   }
}

void
HierarchyEdgeDataOpsComplex::axpy(
   const int dst_id,
   const dcomplex& alpha,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > s1(
            p->getPatchData(src1_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > s2(
            p->getPatchData(src2_id),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.axpy(d, alpha, s1, s2, box);
      }
   }
}

void
HierarchyEdgeDataOpsComplex::axmy(
   const int dst_id,
   const dcomplex& alpha,
   const int src1_id,
   const int src2_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > s1(
            p->getPatchData(src1_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > s2(
            p->getPatchData(src2_id),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.axmy(d, alpha, s1, s2, box);
      }
   }
}

void
HierarchyEdgeDataOpsComplex::abs(
   const int dst_id,
   const int src_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<double> > d(
            p->getPatchData(dst_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > src(
            p->getPatchData(src_id),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.abs(d, src, box);
      }
   }
}

void
HierarchyEdgeDataOpsComplex::setRandomValues(
   const int data_id,
   const dcomplex& width,
   const dcomplex& low,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > d(
            p->getPatchData(data_id),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(d);

         hier::Box box = (interior_only ? p->getBox() : d->getGhostBox());

         d_patch_ops.setRandomValues(d, width, low, box);
      }
   }
}

/*
 *************************************************************************
 *
 * Generic norm and order operations.
 *
 *************************************************************************
 */

int
HierarchyEdgeDataOpsComplex::numberOfEntries(
   const int data_id,
   const bool interior_only) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   int dimVal = d_hierarchy->getDim().getValue();
   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   int entries = 0;

   if (interior_only) {

     boost::shared_ptr<pdat::EdgeDataFactory<dcomplex> > dfact(
        d_hierarchy->getPatchDescriptor()->getPatchDataFactory(data_id),
        boost::detail::dynamic_cast_tag());

      TBOX_ASSERT(dfact);

      for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
         boost::shared_ptr<hier::PatchLevel> level(
            d_hierarchy->getPatchLevel(ln));
         const int npatches = level->getNumberOfPatches();
#ifdef DEBUG_CHECK_ASSERTIONS
         for (int nd = 0; nd < dimVal; nd++) {
            TBOX_ASSERT(npatches == d_nonoverlapping_edge_boxes[nd][ln].getSize());
         }
#endif
         for (int il = 0; il < npatches; il++) {

            for (int eb = 0; eb < dimVal; eb++) {
               hier::BoxContainer::const_iterator lb =
                  ((d_nonoverlapping_edge_boxes[eb][ln])[il]).begin();
               for ( ; lb != ((d_nonoverlapping_edge_boxes[eb][ln])[il]).end();
                    ++lb) {
                  entries += lb->size();
               }
            }
         }
      }

      entries *= dfact->getDepth();

   } else {

      for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
         boost::shared_ptr<hier::PatchLevel> level(
            d_hierarchy->getPatchLevel(ln));
         for (hier::PatchLevel::iterator ip(level->begin());
              ip != level->end(); ++ip) {
            boost::shared_ptr<pdat::EdgeData<dcomplex> > d(
               (*ip)->getPatchData(data_id),
               boost::detail::dynamic_cast_tag());

            TBOX_ASSERT(d);

            entries += d_patch_ops.numberOfEntries(d, d->getGhostBox());
         }
      }

      int global_entries = entries;
      if (mpi.getSize() > 1) {
         mpi.Allreduce(&entries, &global_entries, 1, MPI_INT, MPI_SUM);
      }
      entries = global_entries;

   }

   return entries;
}

double
HierarchyEdgeDataOpsComplex::sumControlVolumes(
   const int data_id,
   const int vol_id) const
{
   TBOX_ASSERT(vol_id >= 0);
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   double sum = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > d(
            p->getPatchData(data_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<double> > cv(
            p->getPatchData(vol_id),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(cv);

         hier::Box box = cv->getGhostBox();

         sum += d_patch_ops.sumControlVolumes(d, cv, box);
      }
   }

   double global_sum = sum;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM);
   }
   return global_sum;
}

double
HierarchyEdgeDataOpsComplex::L1Norm(
   const int data_id,
   const int vol_id) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   double norm = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > d(
            p->getPatchData(data_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<hier::PatchData> pd;

         hier::Box box = p->getBox();
         if (vol_id >= 0) {

            TBOX_ASSERT(d);

            box = d->getGhostBox();
            pd = p->getPatchData(vol_id);
         }

         boost::shared_ptr<pdat::EdgeData<double> > cv(
            pd,
            boost::detail::dynamic_cast_tag());
         norm += d_patch_ops.L1Norm(d, box, cv);
      }
   }

   double global_norm = norm;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM);
   }
   return global_norm;
}

double
HierarchyEdgeDataOpsComplex::L2Norm(
   const int data_id,
   const int vol_id) const
{
   dcomplex dotprod = HierarchyEdgeDataOpsComplex::dot(data_id,
         data_id,
         vol_id);

   return sqrt(real(dotprod));
}

double
HierarchyEdgeDataOpsComplex::weightedL2Norm(
   const int data_id,
   const int wgt_id,
   const int vol_id) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   double norm_squared = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > d(
            p->getPatchData(data_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > w(
            p->getPatchData(wgt_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<hier::PatchData> pd;

         hier::Box box = p->getBox();
         if (vol_id >= 0) {

            TBOX_ASSERT(d);

            box = d->getGhostBox();
            pd = p->getPatchData(vol_id);
         }

         boost::shared_ptr<pdat::EdgeData<double> > cv(
            pd,
            boost::detail::dynamic_cast_tag());
         double pnorm = d_patch_ops.weightedL2Norm(d, w, box, cv);

         norm_squared += pnorm * pnorm;
      }
   }

   double global_norm_squared = norm_squared;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&norm_squared, &global_norm_squared, 1, MPI_DOUBLE, MPI_SUM);
   }
   return sqrt(global_norm_squared);
}

double
HierarchyEdgeDataOpsComplex::RMSNorm(
   const int data_id,
   const int vol_id) const
{
   double l2_norm = L2Norm(data_id, vol_id);

   double volume = ((vol_id < 0) ? (double)numberOfEntries(data_id, true)
                    : sumControlVolumes(data_id, vol_id));

   double rms_norm = l2_norm / sqrt(volume);
   return rms_norm;
}

double
HierarchyEdgeDataOpsComplex::weightedRMSNorm(
   const int data_id,
   const int wgt_id,
   const int vol_id) const
{

   double l2_norm = weightedL2Norm(data_id, wgt_id, vol_id);

   double volume = ((vol_id < 0) ? (double)numberOfEntries(data_id, true)
                    : sumControlVolumes(data_id, vol_id));

   double rms_norm = l2_norm / sqrt(volume);
   return rms_norm;
}

double
HierarchyEdgeDataOpsComplex::maxNorm(
   const int data_id,
   const int vol_id) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   double norm = 0.0;

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > d(
            p->getPatchData(data_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<hier::PatchData> pd;

         hier::Box box = p->getBox();
         if (vol_id >= 0) {

            TBOX_ASSERT(d);

            box = d->getGhostBox();
            pd = p->getPatchData(vol_id);
         }

         boost::shared_ptr<pdat::EdgeData<double> > cv(
            pd,
            boost::detail::dynamic_cast_tag());
         norm = tbox::MathUtilities<double>::Max(norm,
               d_patch_ops.maxNorm(d, box, cv));
      }
   }

   double global_norm = norm;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&norm, &global_norm, 1, MPI_DOUBLE, MPI_MAX);
   }
   return global_norm;
}

dcomplex
HierarchyEdgeDataOpsComplex::dot(
   const int data1_id,
   const int data2_id,
   const int vol_id) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   dcomplex dprod = dcomplex(0.0, 0.0);

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > d1(
            p->getPatchData(data1_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<dcomplex> > d2(
            p->getPatchData(data2_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<hier::PatchData> pd;

         hier::Box box = p->getBox();
         if (vol_id >= 0) {

            TBOX_ASSERT(d1);

            box = d1->getGhostBox();
            pd = p->getPatchData(vol_id);
         }

         boost::shared_ptr<pdat::EdgeData<double> > cv(
            pd,
            boost::detail::dynamic_cast_tag());
         dprod += d_patch_ops.dot(d1, d2, box, cv);
      }
   }

   dcomplex global_dot = dprod;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&dprod, &global_dot, 1, MPI_DOUBLE_COMPLEX, MPI_SUM);
   }
   return global_dot;
}

dcomplex
HierarchyEdgeDataOpsComplex::integral(
   const int data_id,
   const int vol_id) const
{
   TBOX_ASSERT(d_hierarchy);
   TBOX_ASSERT((d_coarsest_level >= 0)
      && (d_finest_level >= d_coarsest_level)
      && (d_finest_level <= d_hierarchy->getFinestLevelNumber()));

   const tbox::SAMRAI_MPI& mpi(d_hierarchy->getMPI());

   dcomplex local_integral = dcomplex(0.0, 0.0);

   for (int ln = d_coarsest_level; ln <= d_finest_level; ln++) {
      boost::shared_ptr<hier::PatchLevel> level(
         d_hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& p = *ip;

         boost::shared_ptr<pdat::EdgeData<dcomplex> > data(
            p->getPatchData(data_id),
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::EdgeData<double> > vol(
            p->getPatchData(vol_id),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(data);
         TBOX_ASSERT(vol);

         hier::Box box = data->getGhostBox();

         local_integral += d_patch_ops.integral(data, box, vol);
      }
   }

   dcomplex global_integral = local_integral;
   if (mpi.getSize() > 1) {
      mpi.Allreduce(&local_integral, &global_integral, 1, MPI_DOUBLE_COMPLEX, MPI_SUM);
   }
   return global_integral;
}

}
}
#endif
