/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Weighted averaging operator for cell-centered double data on
 *                a Skeleton mesh.
 *
 ************************************************************************/

#include "SkeletonCellDoubleWeightedAverage.h"

#include <float.h>
#include <math.h>
#include "SAMRAI/hier/PatchGeometry.h"
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
// in cartcoarsen1d.f:
void F77_FUNC(cartwgtavgcelldoub1d, CARTWGTAVGCELLDOUB1D) (
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
// in cartcoarsen2d.f:
void F77_FUNC(cartwgtavgcelldoub2d, CARTWGTAVGCELLDOUB2D) (
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
// in cartcoarsen3d.f:
void F77_FUNC(cartwgtavgcelldoub3d, CARTWGTAVGCELLDOUB3D) (
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *);
}

using namespace SAMRAI;

SkeletonCellDoubleWeightedAverage::SkeletonCellDoubleWeightedAverage(
   const tbox::Dimension& dim):
   hier::CoarsenOperator(dim, "SKELETON_CONSERVATIVE_COARSEN")
{
}

SkeletonCellDoubleWeightedAverage::~SkeletonCellDoubleWeightedAverage()
{
}

int SkeletonCellDoubleWeightedAverage::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
SkeletonCellDoubleWeightedAverage::getStencilWidth() const {
   return hier::IntVector(getDim(), 0);
}

void SkeletonCellDoubleWeightedAverage::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   boost::shared_ptr<pdat::CellData<double> > fdata(
      fine.getPatchData(src_component),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<pdat::CellData<double> > cdata(
      coarse.getPatchData(dst_component),
      boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
#endif

   const hier::Index filo = fdata->getGhostBox().lower();
   const hier::Index fihi = fdata->getGhostBox().upper();
   const hier::Index cilo = cdata->getGhostBox().lower();
   const hier::Index cihi = cdata->getGhostBox().upper();

   const boost::shared_ptr<hier::PatchGeometry> fgeom(
      fine.getPatchGeometry(),
      boost::detail::dynamic_cast_tag());
   const boost::shared_ptr<hier::PatchGeometry> cgeom(
      coarse.getPatchGeometry(),
      boost::detail::dynamic_cast_tag());

   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();

   int flev_num = fine.getPatchLevelNumber();
   int clev_num = coarse.getPatchLevelNumber();

   // deal with levels not in hierarchy
   if (flev_num < 0) flev_num = clev_num + 1;
   if (clev_num < 0) clev_num = flev_num - 1;

   double cdx[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   double fdx[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   getDx(clev_num, cdx);
   getDx(flev_num, fdx);

   for (int d = 0; d < cdata->getDepth(); d++) {
      if (getDim() == tbox::Dimension(1)) {
         F77_FUNC(cartwgtavgcelldoub1d, CARTWGTAVGCELLDOUB1D) (
            ifirstc(0), ilastc(0),
            filo(0), fihi(0),
            cilo(0), cihi(0),
            &ratio[0],
            fdx,
            cdx,
            fdata->getPointer(d),
            cdata->getPointer(d));
      } else if (getDim() == tbox::Dimension(2)) {
         F77_FUNC(cartwgtavgcelldoub2d, CARTWGTAVGCELLDOUB2D) (
            ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
            filo(0), filo(1), fihi(0), fihi(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            &ratio[0],
            fdx,
            cdx,
            fdata->getPointer(d),
            cdata->getPointer(d));
      } else if (getDim() == tbox::Dimension(3)) {
         F77_FUNC(cartwgtavgcelldoub3d, CARTWGTAVGCELLDOUB3D) (
            ifirstc(0), ifirstc(1), ifirstc(2),
            ilastc(0), ilastc(1), ilastc(2),
            filo(0), filo(1), filo(2),
            fihi(0), fihi(1), fihi(2),
            cilo(0), cilo(1), cilo(2),
            cihi(0), cihi(1), cihi(2),
            &ratio[0],
            fdx,
            cdx,
            fdata->getPointer(d),
            cdata->getPointer(d));
      } else {
         TBOX_ERROR("SkeletonCellDoubleWeightedAverage error...\n"
            << "getDim() > 3 not supported." << endl);

      }
   }
}

void SkeletonCellDoubleWeightedAverage::setDx(
   const int level_number,
   const double* dx)
{
   if (level_number >= d_dx.getSize()) {
      d_dx.resizeArray(level_number + 1);
      d_dx[level_number].resizeArray(getDim().getValue());
      for (int i = 0; i < getDim().getValue(); i++) {
         d_dx[level_number][i] = dx[i];
      }
   }
}

void SkeletonCellDoubleWeightedAverage::getDx(
   const int level_number,
   double* dx) const
{
   for (int i = 0; i < getDim().getValue(); i++) {
      dx[i] = d_dx[level_number][i];
   }
}
