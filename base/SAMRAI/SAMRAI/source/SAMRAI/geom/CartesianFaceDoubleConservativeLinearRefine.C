/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Conservative linear refine operator for face-centered
 *                double data on a Cartesian mesh.
 *
 ************************************************************************/

#ifndef included_geom_CartesianFaceDoubleConservativeLinearRefine_C
#define included_geom_CartesianFaceDoubleConservativeLinearRefine_C

#include "SAMRAI/geom/CartesianFaceDoubleConservativeLinearRefine.h"
#include <float.h>
#include <math.h>
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/tbox/Array.h"
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
void F77_FUNC(cartclinreffacedoub1d, CARTCLINREFFACEDOUB1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *,
   double *, double *);
// in cartrefine2d.f:
void F77_FUNC(cartclinreffacedoub2d0, CARTCLINREFFACEDOUB2D0) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *,
   double *, double *, double *, double *);
void F77_FUNC(cartclinreffacedoub2d1, CARTCLINREFFACEDOUB2D1) (const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *,
   double *, double *, double *, double *);
// in cartrefine3d.f:
void F77_FUNC(cartclinreffacedoub3d0, CARTCLINREFFACEDOUB3D0) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *,
   double *, double *, double *,
   double *, double *, double *);
void F77_FUNC(cartclinreffacedoub3d1, CARTCLINREFFACEDOUB3D1) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *,
   double *, double *, double *,
   double *, double *, double *);
void F77_FUNC(cartclinreffacedoub3d2, CARTCLINREFFACEDOUB3D2) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *, const double *, const double *,
   const double *, double *,
   double *, double *, double *,
   double *, double *, double *);
}

namespace SAMRAI {
namespace geom {

// using namespace std;

CartesianFaceDoubleConservativeLinearRefine::
CartesianFaceDoubleConservativeLinearRefine(
   const tbox::Dimension& dim):
   hier::RefineOperator(dim, "CONSERVATIVE_LINEAR_REFINE")
{
}

CartesianFaceDoubleConservativeLinearRefine::~
CartesianFaceDoubleConservativeLinearRefine()
{
}

int
CartesianFaceDoubleConservativeLinearRefine::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
CartesianFaceDoubleConservativeLinearRefine::getStencilWidth() const
{
   return hier::IntVector::getOne(getDim());
}

void
CartesianFaceDoubleConservativeLinearRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::BoxOverlap& fine_overlap,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(getDim());
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS3(dim, fine, coarse, ratio);

   boost::shared_ptr<pdat::FaceData<double> > cdata(
      coarse.getPatchData(src_component),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<pdat::FaceData<double> > fdata(
      fine.getPatchData(dst_component),
      boost::detail::dynamic_cast_tag());

   const pdat::FaceOverlap* t_overlap =
      dynamic_cast<const pdat::FaceOverlap *>(&fine_overlap);

   TBOX_ASSERT(t_overlap != NULL);

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

   for (int axis = 0; axis < dim.getValue(); axis++) {
      const hier::BoxContainer& boxes = t_overlap->getDestinationBoxContainer(axis);

      for (hier::BoxContainer::const_iterator b(boxes);
           b != boxes.end(); ++b) {

         const hier::Box& face_box = *b;
         TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(dim, face_box);

         hier::Box fine_box(dim);
         for (int i = 0; i < dim.getValue(); i++) {
            fine_box.lower((axis + i) % dim.getValue()) = face_box.lower(i);
            fine_box.upper((axis + i) % dim.getValue()) = face_box.upper(i);
         }

         fine_box.upper(axis) -= 1;

         const hier::Box coarse_box = hier::Box::coarsen(fine_box, ratio);
         const hier::Index ifirstc = coarse_box.lower();
         const hier::Index ilastc = coarse_box.upper();
         const hier::Index ifirstf = fine_box.lower();
         const hier::Index ilastf = fine_box.upper();

         const hier::IntVector tmp_ghosts(dim, 0);
         tbox::Array<double> diff0(cgbox.numberCells(0) + 2);
         pdat::FaceData<double> slope0(cgbox, 1, tmp_ghosts);

         for (int d = 0; d < fdata->getDepth(); d++) {
            if ((dim == tbox::Dimension(1))) {
               F77_FUNC(cartclinreffacedoub1d, CARTCLINREFFACEDOUB1D) (
                  ifirstc(0), ilastc(0),
                  ifirstf(0), ilastf(0),
                  cilo(0), cihi(0),
                  filo(0), fihi(0),
                  &ratio[0],
                  cgeom->getDx(),
                  fgeom->getDx(),
                  cdata->getPointer(0, d),
                  fdata->getPointer(0, d),
                  diff0.getPointer(), slope0.getPointer(0));
            } else if ((dim == tbox::Dimension(2))) {
               tbox::Array<double> diff1(cgbox.numberCells(1) + 2);
               pdat::FaceData<double> slope1(cgbox, 1, tmp_ghosts);

               if (axis == 0) {
                  F77_FUNC(cartclinreffacedoub2d0, CARTCLINREFFACEDOUB2D0) (
                     ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
                     ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
                     cilo(0), cilo(1), cihi(0), cihi(1),
                     filo(0), filo(1), fihi(0), fihi(1),
                     &ratio[0],
                     cgeom->getDx(),
                     fgeom->getDx(),
                     cdata->getPointer(0, d),
                     fdata->getPointer(0, d),
                     diff0.getPointer(), slope0.getPointer(0),
                     diff1.getPointer(), slope1.getPointer(0));
               } else if (axis == 1) {
                  F77_FUNC(cartclinreffacedoub2d1, CARTCLINREFFACEDOUB2D1) (
                     ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
                     ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
                     cilo(0), cilo(1), cihi(0), cihi(1),
                     filo(0), filo(1), fihi(0), fihi(1),
                     &ratio[0],
                     cgeom->getDx(),
                     fgeom->getDx(),
                     cdata->getPointer(1, d),
                     fdata->getPointer(1, d),
                     diff1.getPointer(), slope1.getPointer(1),
                     diff0.getPointer(), slope0.getPointer(1));
               }
            } else if ((dim == tbox::Dimension(3))) {
               tbox::Array<double> diff1(cgbox.numberCells(1) + 2);
               pdat::FaceData<double> slope1(cgbox, 1, tmp_ghosts);

               tbox::Array<double> diff2(cgbox.numberCells(2) + 2);
               pdat::FaceData<double> slope2(cgbox, 1, tmp_ghosts);

               if (axis == 0) {
                  F77_FUNC(cartclinreffacedoub3d0, CARTCLINREFFACEDOUB3D0) (
                     ifirstc(0), ifirstc(1), ifirstc(2),
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
                     cdata->getPointer(0, d),
                     fdata->getPointer(0, d),
                     diff0.getPointer(), slope0.getPointer(0),
                     diff1.getPointer(), slope1.getPointer(0),
                     diff2.getPointer(), slope2.getPointer(0));
               } else if (axis == 1) {
                  F77_FUNC(cartclinreffacedoub3d1, CARTCLINREFFACEDOUB3D1) (
                     ifirstc(0), ifirstc(1), ifirstc(2),
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
                     cdata->getPointer(1, d),
                     fdata->getPointer(1, d),
                     diff1.getPointer(), slope1.getPointer(1),
                     diff2.getPointer(), slope2.getPointer(1),
                     diff0.getPointer(), slope0.getPointer(1));
               } else if (axis == 2) {
                  F77_FUNC(cartclinreffacedoub3d2, CARTCLINREFFACEDOUB3D2) (
                     ifirstc(0), ifirstc(1), ifirstc(2),
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
                     cdata->getPointer(2, d),
                     fdata->getPointer(2, d),
                     diff2.getPointer(), slope2.getPointer(2),
                     diff0.getPointer(), slope0.getPointer(2),
                     diff1.getPointer(), slope1.getPointer(2));
               }
            } else {
               TBOX_ERROR(
                  "CartesianFaceDoubleConservativeLinearRefine error...\n"
                  << "dim > 3 not supported." << std::endl);
            }
         }
      }
   }
}

}
}
#endif
