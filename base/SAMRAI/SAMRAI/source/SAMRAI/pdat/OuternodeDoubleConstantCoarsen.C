/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   ConstantCoarsen averaging operator for outernode-centered
 *                double data on a mesh.
 *
 ************************************************************************/

#ifndef included_pdat_OuternodeDoubleConstantCoarsen_C
#define included_pdat_OuternodeDoubleConstantCoarsen_C

#include "SAMRAI/pdat/OuternodeDoubleConstantCoarsen.h"

#include "SAMRAI/pdat/OuternodeData.h"
#include "SAMRAI/pdat/OuternodeVariable.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/Utilities.h"

#include <cfloat>
#include <cmath>

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

// in pdat_concoarsen1d.f:
void F77_FUNC(conavgouternodedoub1d, CONAVGOUTERNODEDOUB1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const double *, double *);
// in pdat_concoarsen2d.f:
void F77_FUNC(conavgouternodedoub2d0, CONAVGOUTERNODEDOUB2D0) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const double *, double *);
void F77_FUNC(conavgouternodedoub2d1, CONAVGOUTERNODEDOUB2D1) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const double *, double *);

// in pdat_concoarsen3d.f:
void F77_FUNC(conavgouternodedoub3d0, CONAVGOUTERNODEDOUB3D0) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const double *, double *);
void F77_FUNC(conavgouternodedoub3d1, CONAVGOUTERNODEDOUB3D1) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const double *, double *);
void F77_FUNC(conavgouternodedoub3d2, CONAVGOUTERNODEDOUB3D2) (const int&,
   const int&, const int&,
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

OuternodeDoubleConstantCoarsen::OuternodeDoubleConstantCoarsen(
   const tbox::Dimension& dim):
   hier::CoarsenOperator(dim, "CONSTANT_COARSEN")
{
}

OuternodeDoubleConstantCoarsen::~OuternodeDoubleConstantCoarsen()
{
}

int
OuternodeDoubleConstantCoarsen::getOperatorPriority() const
{
   return 0;
}

hier::IntVector
OuternodeDoubleConstantCoarsen::getStencilWidth() const
{
   return hier::IntVector::getZero(getDim());
}

void
OuternodeDoubleConstantCoarsen::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(getDim());

   boost::shared_ptr<OuternodeData<double> > fdata(
      fine.getPatchData(src_component),
      boost::detail::dynamic_cast_tag());
   boost::shared_ptr<OuternodeData<double> > cdata(
      coarse.getPatchData(dst_component),
      boost::detail::dynamic_cast_tag());

   TBOX_ASSERT(fdata);
   TBOX_ASSERT(cdata);
   TBOX_ASSERT(cdata->getDepth() == fdata->getDepth());
   TBOX_DIM_ASSERT_CHECK_ARGS5(*this, coarse, fine, coarse_box, ratio);

   const hier::Index filo = fine.getBox().lower();
   const hier::Index fihi = fine.getBox().upper();
   const hier::Index cilo = coarse.getBox().lower();
   const hier::Index cihi = coarse.getBox().upper();

   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();

   for (int i = 0; i < 2; i++) {

      for (int axis = 0; axis < dim.getValue(); axis++) {

         if (cdata->dataExists(axis)) {

            for (int d = 0; d < cdata->getDepth(); d++) {

               if (dim == tbox::Dimension(1)) {
                  F77_FUNC(conavgouternodedoub1d,
                     CONAVGOUTERNODEDOUB1D) (ifirstc(0), ilastc(0),
                     filo(0), fihi(0),
                     cilo(0), cihi(0),
                     &ratio[0],
                     fdata->getPointer(axis, i, d),
                     cdata->getPointer(axis, i, d));
               } else if (dim == tbox::Dimension(2)) {
                  if (axis == 0) {
                     F77_FUNC(conavgouternodedoub2d0,
                        CONAVGOUTERNODEDOUB2D0) (ifirstc(0), ifirstc(1),
                        ilastc(0), ilastc(1),
                        filo(0), filo(1),
                        fihi(0), fihi(1),
                        cilo(0), cilo(1),
                        cihi(0), cihi(1),
                        &ratio[0],
                        fdata->getPointer(axis, i, d),
                        cdata->getPointer(axis, i, d));
                  }

                  if (axis == 1) {
                     F77_FUNC(conavgouternodedoub2d1,
                        CONAVGOUTERNODEDOUB2D1) (ifirstc(0), ifirstc(1),
                        ilastc(0), ilastc(1),
                        filo(0), filo(1),
                        fihi(0), fihi(1),
                        cilo(0), cilo(1),
                        cihi(0), cihi(1),
                        &ratio[0],
                        fdata->getPointer(axis, i, d),
                        cdata->getPointer(axis, i, d));
                  }
               } else if (dim == tbox::Dimension(3)) {
                  if (axis == 0) {
                     F77_FUNC(conavgouternodedoub3d0,
                        CONAVGOUTERNODEDOUB3D0) (ifirstc(0), ifirstc(1),
                        ifirstc(2),
                        ilastc(0), ilastc(1), ilastc(2),
                        filo(0), filo(1), filo(2),
                        fihi(0), fihi(1), fihi(2),
                        cilo(0), cilo(1), cilo(2),
                        cihi(0), cihi(1), cihi(2),
                        &ratio[0],
                        fdata->getPointer(axis, i, d),
                        cdata->getPointer(axis, i, d));
                  }
                  if (axis == 1) {
                     F77_FUNC(conavgouternodedoub3d1,
                        CONAVGOUTERNODEDOUB3D1) (ifirstc(0), ifirstc(1),
                        ifirstc(2),
                        ilastc(0), ilastc(1), ilastc(2),
                        filo(0), filo(1), filo(2),
                        fihi(0), fihi(1), fihi(2),
                        cilo(0), cilo(1), cilo(2),
                        cihi(0), cihi(1), cihi(2),
                        &ratio[0],
                        fdata->getPointer(axis, i, d),
                        cdata->getPointer(axis, i, d));
                  }
                  if (axis == 2) {
                     F77_FUNC(conavgouternodedoub3d2,
                        CONAVGOUTERNODEDOUB3D2) (ifirstc(0), ifirstc(1),
                        ifirstc(2),
                        ilastc(0), ilastc(1), ilastc(2),
                        filo(0), filo(1), filo(2),
                        fihi(0), fihi(1), fihi(2),
                        cilo(0), cilo(1), cilo(2),
                        cihi(0), cihi(1), cihi(2),
                        &ratio[0],
                        fdata->getPointer(axis, i, d),
                        cdata->getPointer(axis, i, d));
                  }
               }

            }
         }
      }
   }
}

}
}

#endif
