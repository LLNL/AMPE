/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Linear refine operator for cell-centered complex data on
 *                a Cartesian mesh.
 *
 ************************************************************************/

#ifndef included_geom_CartesianCellComplexLinearRefine_C
#define included_geom_CartesianCellComplexLinearRefine_C

#include "SAMRAI/geom/CartesianCellComplexLinearRefine.h"
#include "SAMRAI/tbox/Complex.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"

#include "SAMRAI/tbox/Utilities.h"

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

// in cartrefine1d.f:
void F77_FUNC(cartlinrefcellcplx1d, CARTLINREFCELLCPLX1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *);
// in cartrefine2d.f:
void F77_FUNC(cartlinrefcellcplx2d, CARTLINREFCELLCPLX2D) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *);
// in cartrefine3d.f:
void F77_FUNC(cartlinrefcellcplx3d, CARTLINREFCELLCPLX3D) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *);
}

namespace SAMRAI {
namespace geom {

// using namespace std;

CartesianCellComplexLinearRefine::CartesianCellComplexLinearRefine(
   const tbox::Dimension& dim):
   hier::RefineOperator(dim, "LINEAR_REFINE")
{
}

CartesianCellComplexLinearRefine::~CartesianCellComplexLinearRefine()
{
}

int
CartesianCellComplexLinearRefine::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianCellComplexLinearRefine::getStencilWidth() const
{
   return hier::IntVector::getOne(getDim());
}

void
CartesianCellComplexLinearRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::BoxOverlap& fine_overlap,
   const hier::IntVector& ratio) const
{
   const pdat::CellOverlap* t_overlap =
      dynamic_cast<const pdat::CellOverlap *>(&fine_overlap);

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
CartesianCellComplexLinearRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::Box& fine_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(getDim());
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim, fine, coarse, fine_box, ratio);

   boost::shared_ptr<pdat::CellData<dcomplex> > cdata(
      coarse.getPatchData(src_component),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<pdat::CellData<dcomplex> > fdata(
      fine.getPatchData(dst_component),
      boost::detail::dynamic_cast_tag());

   TBOX_ASSERT(cdata);
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());

   const hier::Box cgbox(cdata->getGhostBox());

   const hier::Index cilo = cgbox.lower();
   const hier::Index cihi = cgbox.upper();
   const hier::Index filo = fdata->getGhostBox().lower();
   const hier::Index fihi = fdata->getGhostBox().upper();

   const boost::shared_ptr<CartesianPatchGeometry> cgeom(
      coarse.getPatchGeometry(),
      boost::detail::dynamic_cast_tag());
   const boost::shared_ptr<CartesianPatchGeometry> fgeom(
      fine.getPatchGeometry(),
      boost::detail::dynamic_cast_tag());

   const hier::Box coarse_box = hier::Box::coarsen(fine_box, ratio);
   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();
   const hier::Index ifirstf = fine_box.lower();
   const hier::Index ilastf = fine_box.upper();

   for (int d = 0; d < fdata->getDepth(); d++) {
      if ((dim == tbox::Dimension(1))) {
         F77_FUNC(cartlinrefcellcplx1d, CARTLINREFCELLCPLX1D) (ifirstc(0),
            ilastc(0),
            ifirstf(0), ilastf(0),
            cilo(0), cihi(0),
            filo(0), fihi(0),
            &ratio[0],
            cgeom->getDx(),
            fgeom->getDx(),
            cdata->getPointer(d),
            fdata->getPointer(d));
      } else if ((dim == tbox::Dimension(2))) {
         F77_FUNC(cartlinrefcellcplx2d, CARTLINREFCELLCPLX2D) (ifirstc(0),
            ifirstc(1), ilastc(0), ilastc(1),
            ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            filo(0), filo(1), fihi(0), fihi(1),
            &ratio[0],
            cgeom->getDx(),
            fgeom->getDx(),
            cdata->getPointer(d),
            fdata->getPointer(d));
      } else if ((dim == tbox::Dimension(3))) {
         F77_FUNC(cartlinrefcellcplx3d, CARTLINREFCELLCPLX3D) (ifirstc(0),
            ifirstc(1), ifirstc(2),
            ilastc(0), ilastc(1), ilastc(2),
            ifirstf(0), ifirstf(1), ifirstf(2),
            ilastf(0), ilastf(1), ilastf(2),
            cilo(0), cilo(1), cilo(2),
            cihi(0), cihi(1), cihi(2),
            filo(0), filo(1), filo(2),
            fihi(0), fihi(1), fihi(2),
            &ratio[0],
            cgeom->getDx(),
            fgeom->getDx(),
            cdata->getPointer(d),
            fdata->getPointer(d));
      } else {
         TBOX_ERROR("CartesianCellComplexLinearLinearRefine error...\n"
            << "dim > 3 not supported." << std::endl);
      }

   }
}

}
}
#endif
