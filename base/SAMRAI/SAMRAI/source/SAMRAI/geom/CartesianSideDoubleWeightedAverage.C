/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for side-centered double data on
 *                a Cartesian mesh.
 *
 ************************************************************************/

#ifndef included_geom_CartesianSideDoubleWeightedAverage_C
#define included_geom_CartesianSideDoubleWeightedAverage_C

#include "SAMRAI/geom/CartesianSideDoubleWeightedAverage.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
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

// in cartcoarsen1d.f:
void F77_FUNC(cartwgtavgsidedoub1d, CARTWGTAVGSIDEDOUB1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
// in cartcoarsen2d.f:
void F77_FUNC(cartwgtavgsidedoub2d0, CARTWGTAVGSIDEDOUB2D0) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);

void F77_FUNC(cartwgtavgsidedoub2d1, CARTWGTAVGSIDEDOUB2D1) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
// in cartcoarsen3d.f:
void F77_FUNC(cartwgtavgsidedoub3d0, CARTWGTAVGSIDEDOUB3D0) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
void F77_FUNC(cartwgtavgsidedoub3d1, CARTWGTAVGSIDEDOUB3D1) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
void F77_FUNC(cartwgtavgsidedoub3d2, CARTWGTAVGSIDEDOUB3D2) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
}

namespace SAMRAI {
namespace geom {

// using namespace std;

CartesianSideDoubleWeightedAverage::CartesianSideDoubleWeightedAverage(
   const tbox::Dimension& dim):
   hier::CoarsenOperator(dim, "CONSERVATIVE_COARSEN")
{
}

CartesianSideDoubleWeightedAverage::~CartesianSideDoubleWeightedAverage()
{
}

int
CartesianSideDoubleWeightedAverage::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianSideDoubleWeightedAverage::getStencilWidth() const
{
   return hier::IntVector::getZero(getDim());
}

void
CartesianSideDoubleWeightedAverage::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(getDim());

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim, coarse, fine, coarse_box, ratio);

   boost::shared_ptr<pdat::SideData<double> > fdata(
      fine.getPatchData(src_component),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<pdat::SideData<double> > cdata(
      coarse.getPatchData(dst_component),
      boost::detail::dynamic_cast_tag());

   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());

   const hier::IntVector& directions(cdata->getDirectionVector());

   TBOX_ASSERT(directions ==
      hier::IntVector::min(directions, fdata->getDirectionVector()));

   const hier::Index filo = fdata->getGhostBox().lower();
   const hier::Index fihi = fdata->getGhostBox().upper();
   const hier::Index cilo = cdata->getGhostBox().lower();
   const hier::Index cihi = cdata->getGhostBox().upper();

   const boost::shared_ptr<CartesianPatchGeometry> fgeom(
      fine.getPatchGeometry(),
      boost::detail::dynamic_cast_tag());
   const boost::shared_ptr<CartesianPatchGeometry> cgeom(
      coarse.getPatchGeometry(),
      boost::detail::dynamic_cast_tag());

   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();

   for (int d = 0; d < cdata->getDepth(); d++) {
      if ((dim == tbox::Dimension(1))) {
         if (directions(0)) {
            F77_FUNC(cartwgtavgsidedoub1d, CARTWGTAVGSIDEDOUB1D) (ifirstc(0),
               ilastc(0),
               filo(0), fihi(0),
               cilo(0), cihi(0),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(0, d),
               cdata->getPointer(0, d));
         }
      } else if ((dim == tbox::Dimension(2))) {
         if (directions(0)) {
            F77_FUNC(cartwgtavgsidedoub2d0, CARTWGTAVGSIDEDOUB2D0) (ifirstc(0),
               ifirstc(1), ilastc(0), ilastc(1),
               filo(0), filo(1), fihi(0), fihi(1),
               cilo(0), cilo(1), cihi(0), cihi(1),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(0, d),
               cdata->getPointer(0, d));
         }
         if (directions(1)) {
            F77_FUNC(cartwgtavgsidedoub2d1, CARTWGTAVGSIDEDOUB2D1) (ifirstc(0),
               ifirstc(1), ilastc(0), ilastc(1),
               filo(0), filo(1), fihi(0), fihi(1),
               cilo(0), cilo(1), cihi(0), cihi(1),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(1, d),
               cdata->getPointer(1, d));
         }
      } else if ((dim == tbox::Dimension(3))) {
         if (directions(0)) {
            F77_FUNC(cartwgtavgsidedoub3d0, CARTWGTAVGSIDEDOUB3D0) (ifirstc(0),
               ifirstc(1), ifirstc(2),
               ilastc(0), ilastc(1), ilastc(2),
               filo(0), filo(1), filo(2),
               fihi(0), fihi(1), fihi(2),
               cilo(0), cilo(1), cilo(2),
               cihi(0), cihi(1), cihi(2),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(0, d),
               cdata->getPointer(0, d));
         }
         if (directions(1)) {
            F77_FUNC(cartwgtavgsidedoub3d1, CARTWGTAVGSIDEDOUB3D1) (ifirstc(0),
               ifirstc(1), ifirstc(2),
               ilastc(0), ilastc(1), ilastc(2),
               filo(0), filo(1), filo(2),
               fihi(0), fihi(1), fihi(2),
               cilo(0), cilo(1), cilo(2),
               cihi(0), cihi(1), cihi(2),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(1, d),
               cdata->getPointer(1, d));
         }
         if (directions(2)) {
            F77_FUNC(cartwgtavgsidedoub3d2, CARTWGTAVGSIDEDOUB3D2) (ifirstc(0),
               ifirstc(1), ifirstc(2),
               ilastc(0), ilastc(1), ilastc(2),
               filo(0), filo(1), filo(2),
               fihi(0), fihi(1), fihi(2),
               cilo(0), cilo(1), cilo(2),
               cihi(0), cihi(1), cihi(2),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(2, d),
               cdata->getPointer(2, d));
         }
      } else {
         TBOX_ERROR("CartesianCellSideDoubleWeightedAverage error...\n"
            << "dim > 3 not supported." << std::endl);
      }
   }
}

}
}
#endif
