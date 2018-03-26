/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for outerface double data on
 *                a Cartesian mesh.
 *
 ************************************************************************/

#ifndef included_geom_CartesianOuterfaceDoubleWeightedAverage_C
#define included_geom_CartesianOuterfaceDoubleWeightedAverage_C

#include "SAMRAI/geom/CartesianOuterfaceDoubleWeightedAverage.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/OuterfaceData.h"
#include "SAMRAI/pdat/OuterfaceVariable.h"
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
void F77_FUNC(cartwgtavgoutfacedoub1d, CARTWGTAVGOUTFACEDOUB1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
// in cartcoarsen2d.f:
void F77_FUNC(cartwgtavgoutfacedoub2d0, CARTWGTAVGOUTFACEDOUB2D0) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);

void F77_FUNC(cartwgtavgoutfacedoub2d1, CARTWGTAVGOUTFACEDOUB2D1) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
// in cartcoarsen3d.f:
void F77_FUNC(cartwgtavgoutfacedoub3d0, CARTWGTAVGOUTFACEDOUB3D0) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
void F77_FUNC(cartwgtavgoutfacedoub3d1, CARTWGTAVGOUTFACEDOUB3D1) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
void F77_FUNC(cartwgtavgoutfacedoub3d2, CARTWGTAVGOUTFACEDOUB3D2) (const int&,
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

CartesianOuterfaceDoubleWeightedAverage::
CartesianOuterfaceDoubleWeightedAverage(
   const tbox::Dimension& dim):
   hier::CoarsenOperator(dim, "CONSERVATIVE_COARSEN")
{
}

CartesianOuterfaceDoubleWeightedAverage::~
CartesianOuterfaceDoubleWeightedAverage()
{
}

int
CartesianOuterfaceDoubleWeightedAverage::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianOuterfaceDoubleWeightedAverage::getStencilWidth() const
{
   return hier::IntVector::getZero(getDim());
}

void
CartesianOuterfaceDoubleWeightedAverage::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(getDim());

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim, coarse, fine, coarse_box, ratio);

   boost::shared_ptr<pdat::OuterfaceData<double> > fdata(
      fine.getPatchData(src_component),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<pdat::OuterfaceData<double> > cdata(
      coarse.getPatchData(dst_component),
      boost::detail::dynamic_cast_tag());
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());

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
      // loop over lower and upper outerface arrays
      for (int i = 0; i < 2; i++) {
         if ((dim == tbox::Dimension(1))) {
            F77_FUNC(cartwgtavgoutfacedoub1d,
               CARTWGTAVGOUTFACEDOUB1D) (ifirstc(0), ilastc(0),
               filo(0), fihi(0),
               cilo(0), cihi(0),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(0, i, d),
               cdata->getPointer(0, i, d));
         } else if ((dim == tbox::Dimension(2))) {
            F77_FUNC(cartwgtavgoutfacedoub2d0,
               CARTWGTAVGOUTFACEDOUB2D0) (ifirstc(0), ifirstc(1), ilastc(0),
               ilastc(1),
               filo(0), filo(1), fihi(0), fihi(1),
               cilo(0), cilo(1), cihi(0), cihi(1),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(0, i, d),
               cdata->getPointer(0, i, d));
            F77_FUNC(cartwgtavgoutfacedoub2d1,
               CARTWGTAVGOUTFACEDOUB2D1) (ifirstc(0), ifirstc(1), ilastc(0),
               ilastc(1),
               filo(0), filo(1), fihi(0), fihi(1),
               cilo(0), cilo(1), cihi(0), cihi(1),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(1, i, d),
               cdata->getPointer(1, i, d));
         } else if ((dim == tbox::Dimension(3))) {
            F77_FUNC(cartwgtavgoutfacedoub3d0,
               CARTWGTAVGOUTFACEDOUB3D0) (ifirstc(0), ifirstc(1), ifirstc(2),
               ilastc(0), ilastc(1), ilastc(2),
               filo(0), filo(1), filo(2),
               fihi(0), fihi(1), fihi(2),
               cilo(0), cilo(1), cilo(2),
               cihi(0), cihi(1), cihi(2),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(0, i, d),
               cdata->getPointer(0, i, d));
            F77_FUNC(cartwgtavgoutfacedoub3d1,
               CARTWGTAVGOUTFACEDOUB3D1) (ifirstc(0), ifirstc(1), ifirstc(2),
               ilastc(0), ilastc(1), ilastc(2),
               filo(0), filo(1), filo(2),
               fihi(0), fihi(1), fihi(2),
               cilo(0), cilo(1), cilo(2),
               cihi(0), cihi(1), cihi(2),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(1, i, d),
               cdata->getPointer(1, i, d));
            F77_FUNC(cartwgtavgoutfacedoub3d2,
               CARTWGTAVGOUTFACEDOUB3D2) (ifirstc(0), ifirstc(1), ifirstc(2),
               ilastc(0), ilastc(1), ilastc(2),
               filo(0), filo(1), filo(2),
               fihi(0), fihi(1), fihi(2),
               cilo(0), cilo(1), cilo(2),
               cihi(0), cihi(1), cihi(2),
               &ratio[0],
               fgeom->getDx(),
               cgeom->getDx(),
               fdata->getPointer(2, i, d),
               cdata->getPointer(2, i, d));
         } else {
            TBOX_ERROR("CartesianOuterfaceDoubleWeightedAverage error...\n"
               << "dim > 3 not supported." << std::endl);
         }
      }
   }
}

}
}
#endif
