/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Constant refine operator for face-centered complex data on
 *                a  mesh.
 *
 ************************************************************************/

#ifndef included_pdat_FaceComplexConstantRefine_C
#define included_pdat_FaceComplexConstantRefine_C

#include "SAMRAI/pdat/FaceComplexConstantRefine.h"

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/FaceData.h"
#include "SAMRAI/pdat/FaceVariable.h"

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

// in refine1d.f:
void F77_FUNC(conreffacecplx1d, CONREFFACECPLX1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
// in refine2d.f:
void F77_FUNC(conreffacecplx2d0, CONREFFACECPLX2D0) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
void F77_FUNC(conreffacecplx2d1, CONREFFACECPLX2D1) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
// in refine3d.f:
void F77_FUNC(conreffacecplx3d0, CONREFFACECPLX3D0) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
void F77_FUNC(conreffacecplx3d1, CONREFFACECPLX3D1) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
void F77_FUNC(conreffacecplx3d2, CONREFFACECPLX3D2) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const dcomplex *, dcomplex *);
}

namespace SAMRAI {
namespace pdat {

FaceComplexConstantRefine::FaceComplexConstantRefine(
   const tbox::Dimension& dim):
   hier::RefineOperator(dim, "CONSTANT_REFINE")
{
}

FaceComplexConstantRefine::~FaceComplexConstantRefine()
{
}

int
FaceComplexConstantRefine::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
FaceComplexConstantRefine::getStencilWidth() const
{
   return hier::IntVector::getZero(getDim());
}

void
FaceComplexConstantRefine::refine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const int dst_component,
   const int src_component,
   const hier::BoxOverlap& fine_overlap,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(getDim());

   boost::shared_ptr<FaceData<dcomplex> > cdata(
      coarse.getPatchData(src_component),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<FaceData<dcomplex> > fdata(
      fine.getPatchData(dst_component),
      boost::detail::dynamic_cast_tag());

   const FaceOverlap* t_overlap =
      dynamic_cast<const FaceOverlap *>(&fine_overlap);

   TBOX_ASSERT(t_overlap != NULL);

   TBOX_ASSERT(cdata);
   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
   TBOX_DIM_ASSERT_CHECK_ARGS4(*this, fine, coarse, ratio);

   const hier::Box cgbox(cdata->getGhostBox());

   const hier::Index cilo = cgbox.lower();
   const hier::Index cihi = cgbox.upper();
   const hier::Index filo = fdata->getGhostBox().lower();
   const hier::Index fihi = fdata->getGhostBox().upper();

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

         for (int d = 0; d < fdata->getDepth(); d++) {
            if (dim == tbox::Dimension(1)) {
               F77_FUNC(conreffacecplx1d, CONREFFACECPLX1D) (
                  ifirstc(0), ilastc(0),
                  ifirstf(0), ilastf(0),
                  cilo(0), cihi(0),
                  filo(0), fihi(0),
                  &ratio[0],
                  cdata->getPointer(0, d),
                  fdata->getPointer(0, d));
            } else if (dim == tbox::Dimension(2)) {
               if (axis == 0) {
                  F77_FUNC(conreffacecplx2d0, CONREFFACECPLX2D0) (
                     ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
                     ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
                     cilo(0), cilo(1), cihi(0), cihi(1),
                     filo(0), filo(1), fihi(0), fihi(1),
                     &ratio[0],
                     cdata->getPointer(0, d),
                     fdata->getPointer(0, d));
               } else if (axis == 1) {
                  F77_FUNC(conreffacecplx2d1, CONREFFACECPLX2D1) (
                     ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
                     ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
                     cilo(0), cilo(1), cihi(0), cihi(1),
                     filo(0), filo(1), fihi(0), fihi(1),
                     &ratio[0],
                     cdata->getPointer(1, d),
                     fdata->getPointer(1, d));
               }
            } else if (dim == tbox::Dimension(3)) {
               if (axis == 0) {
                  F77_FUNC(conreffacecplx3d0, CONREFFACECPLX3D0) (
                     ifirstc(0), ifirstc(1), ifirstc(2),
                     ilastc(0), ilastc(1), ilastc(2),
                     ifirstf(0), ifirstf(1), ifirstf(2),
                     ilastf(0), ilastf(1), ilastf(2),
                     cilo(0), cilo(1), cilo(2),
                     cihi(0), cihi(1), cihi(2),
                     filo(0), filo(1), filo(2),
                     fihi(0), fihi(1), fihi(2),
                     &ratio[0],
                     cdata->getPointer(0, d),
                     fdata->getPointer(0, d));
               } else if (axis == 1) {
                  F77_FUNC(conreffacecplx3d1, CONREFFACECPLX3D1) (
                     ifirstc(0), ifirstc(1), ifirstc(2),
                     ilastc(0), ilastc(1), ilastc(2),
                     ifirstf(0), ifirstf(1), ifirstf(2),
                     ilastf(0), ilastf(1), ilastf(2),
                     cilo(0), cilo(1), cilo(2),
                     cihi(0), cihi(1), cihi(2),
                     filo(0), filo(1), filo(2),
                     fihi(0), fihi(1), fihi(2),
                     &ratio[0],
                     cdata->getPointer(1, d),
                     fdata->getPointer(1, d));
               } else if (axis == 2) {
                  F77_FUNC(conreffacecplx3d2, CONREFFACECPLX3D2) (
                     ifirstc(0), ifirstc(1), ifirstc(2),
                     ilastc(0), ilastc(1), ilastc(2),
                     ifirstf(0), ifirstf(1), ifirstf(2),
                     ilastf(0), ilastf(1), ilastf(2),
                     cilo(0), cilo(1), cilo(2),
                     cihi(0), cihi(1), cihi(2),
                     filo(0), filo(1), filo(2),
                     fihi(0), fihi(1), fihi(2),
                     &ratio[0],
                     cdata->getPointer(2, d),
                     fdata->getPointer(2, d));
               }
            } else {
               TBOX_ERROR(
                  "FaceComplexConstantRefine::refine() dimension > 3 not supported"
                  << std::endl);
            }
         }
      }
   }
}

}
}
#endif
