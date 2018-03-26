/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for cell-centered complex data on
 *                a Cartesian mesh.
 *
 ************************************************************************/

#ifndef included_geom_CartesianCellComplexWeightedAverage_C
#define included_geom_CartesianCellComplexWeightedAverage_C

#include "SAMRAI/geom/CartesianCellComplexWeightedAverage.h"
#include "SAMRAI/tbox/Complex.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
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

// in cartcoarsen1d.f:
void F77_FUNC(cartwgtavgcellcplx1d, CARTWGTAVGCELLCPLX1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *);
// in cartcoarsen2d.f:
void F77_FUNC(cartwgtavgcellcplx2d, CARTWGTAVGCELLCPLX2D) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const dcomplex *, dcomplex *);
// in cartcoarsen3d.f:
void F77_FUNC(cartwgtavgcellcplx3d, CARTWGTAVGCELLCPLX3D) (const int&,
   const int&, const int&,
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

CartesianCellComplexWeightedAverage::CartesianCellComplexWeightedAverage(
   const tbox::Dimension& dim):
   hier::CoarsenOperator(dim, "CONSERVATIVE_COARSEN")
{
}

CartesianCellComplexWeightedAverage::~CartesianCellComplexWeightedAverage()
{
}

int
CartesianCellComplexWeightedAverage::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianCellComplexWeightedAverage::getStencilWidth() const
{
   return hier::IntVector::getZero(getDim());
}

void
CartesianCellComplexWeightedAverage::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(getDim());
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS4(dim, coarse, fine, coarse_box, ratio);

   boost::shared_ptr<pdat::CellData<dcomplex> > fdata(
      fine.getPatchData(src_component),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<pdat::CellData<dcomplex> > cdata(
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
      if ((dim == tbox::Dimension(1))) {
         F77_FUNC(cartwgtavgcellcplx1d, CARTWGTAVGCELLCPLX1D) (ifirstc(0),
            ilastc(0),
            filo(0), fihi(0),
            cilo(0), cihi(0),
            &ratio[0],
            fgeom->getDx(),
            cgeom->getDx(),
            fdata->getPointer(d),
            cdata->getPointer(d));
      }
      if ((dim == tbox::Dimension(2))) {
         F77_FUNC(cartwgtavgcellcplx2d, CARTWGTAVGCELLCPLX2D) (ifirstc(0),
            ifirstc(1), ilastc(0), ilastc(1),
            filo(0), filo(1), fihi(0), fihi(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            &ratio[0],
            fgeom->getDx(),
            cgeom->getDx(),
            fdata->getPointer(d),
            cdata->getPointer(d));
      }
      if ((dim == tbox::Dimension(3))) {
         F77_FUNC(cartwgtavgcellcplx3d, CARTWGTAVGCELLCPLX3D) (ifirstc(0),
            ifirstc(1), ifirstc(2),
            ilastc(0), ilastc(1), ilastc(2),
            filo(0), filo(1), filo(2),
            fihi(0), fihi(1), fihi(2),
            cilo(0), cilo(1), cilo(2),
            cihi(0), cihi(1), cihi(2),
            &ratio[0],
            fgeom->getDx(),
            cgeom->getDx(),
            fdata->getPointer(d),
            cdata->getPointer(d));
      } else {
         TBOX_ERROR("CartesianEdgeComplexWeightedAverage error...\n"
            << "dim > 3 not supported." << std::endl);
      }
   }
}

}
}
#endif
