/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Constant refine operator for cell-centered double data on
 *                a  mesh.
 *
 ************************************************************************/

#ifndef included_pdat_CellDoubleConstantRefine_C
#define included_pdat_CellDoubleConstantRefine_C

#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/Utilities.h"

#include <float.h>
#include <math.h>

/*
 *************************************************************************
 *
 * External declarations for FORTRAN  routines.
 *
 *************************************************************************
 */
extern "C" {

#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif

// in conrefine1d.f:
void F77_FUNC(conrefcelldoub1d, CONREFCELLDOUB1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const double *, double *);
// in conrefine2d.f:
void F77_FUNC(conrefcelldoub2d, CONREFCELLDOUB2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const double *, double *);
// in conrefine3d.f:
void F77_FUNC(conrefcelldoub3d, CONREFCELLDOUB3D) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const double *, double *);
}

namespace SAMRAI {
namespace pdat {

CellDoubleConstantRefine::CellDoubleConstantRefine(
   const tbox::Dimension& dim):
   hier::RefineOperator(dim, "CONSTANT_REFINE")
{
}

CellDoubleConstantRefine::~CellDoubleConstantRefine()
{
}

int
CellDoubleConstantRefine::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CellDoubleConstantRefine::getStencilWidth() const
{
   return hier::IntVector::getZero(getDim());
}

void
CellDoubleConstantRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::BoxOverlap& fine_overlap,
   const hier::IntVector& ratio) const
{
   const CellOverlap* t_overlap =
      dynamic_cast<const CellOverlap *>(&fine_overlap);

   TBOX_ASSERT(t_overlap != NULL);

   const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer();
   for (hier::BoxContainer::const_iterator b(boxes); b != boxes.end(); ++b) {
      refine(fine,
         coarse,
         dst_component,
         src_component,
         *b,
         ratio);
   }
}

void
CellDoubleConstantRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   boost::shared_ptr<CellData<double> > cdata(
      coarse.getPatchData(src_component),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<CellData<double> > fdata(
      fine.getPatchData(dst_component),
      boost::detail::dynamic_cast_tag());

   TBOX_ASSERT(cdata);
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
   TBOX_DIM_ASSERT_CHECK_ARGS5(*this, fine, coarse, fine_box, ratio);

   const hier::Box cgbox(cdata->getGhostBox());

   const hier::Index cilo = cgbox.lower();
   const hier::Index cihi = cgbox.upper();
   const hier::Index filo = fdata->getGhostBox().lower();
   const hier::Index fihi = fdata->getGhostBox().upper();

   const hier::Box coarse_box = hier::Box::coarsen(fine_box, ratio);
   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();
   const hier::Index ifirstf = fine_box.lower();
   const hier::Index ilastf = fine_box.upper();

   for (int d = 0; d < fdata->getDepth(); d++) {
      if (getDim() == tbox::Dimension(1)) {
         F77_FUNC(conrefcelldoub1d, CONREFCELLDOUB1D) (ifirstc(0), ilastc(0),
            ifirstf(0), ilastf(0),
            cilo(0), cihi(0),
            filo(0), fihi(0),
            &ratio[0],
            cdata->getPointer(d),
            fdata->getPointer(d));
      } else if (getDim() == tbox::Dimension(2)) {
         F77_FUNC(conrefcelldoub2d, CONREFCELLDOUB2D) (ifirstc(0), ifirstc(1),
            ilastc(0), ilastc(1),
            ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            filo(0), filo(1), fihi(0), fihi(1),
            &ratio[0],
            cdata->getPointer(d),
            fdata->getPointer(d));
      } else if (getDim() == tbox::Dimension(3)) {
         F77_FUNC(conrefcelldoub3d, CONREFCELLDOUB3D) (ifirstc(0), ifirstc(1),
            ifirstc(2),
            ilastc(0), ilastc(1), ilastc(2),
            ifirstf(0), ifirstf(1), ifirstf(2),
            ilastf(0), ilastf(1), ilastf(2),
            cilo(0), cilo(1), cilo(2),
            cihi(0), cihi(1), cihi(2),
            filo(0), filo(1), filo(2),
            fihi(0), fihi(1), fihi(2),
            &ratio[0],
            cdata->getPointer(d),
            fdata->getPointer(d));
      } else {
         TBOX_ERROR(
            "CellDoubleConstantRefine::refine dimension > 3 not supported"
            << std::endl);
      }
   }
}

}
}
#endif
