/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Linear time interp operator for double outerface patch data.
 *
 ************************************************************************/

#ifndef included_pdat_OuterfaceDoubleLinearTimeInterpolateOp_C
#define included_pdat_OuterfaceDoubleLinearTimeInterpolateOp_C

#include "SAMRAI/pdat/OuterfaceDoubleLinearTimeInterpolateOp.h"

#include "SAMRAI/pdat/OuterfaceData.h"
#include "SAMRAI/pdat/OuterfaceVariable.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <boost/shared_ptr.hpp>

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

// in lintimint1d.f:
void F77_FUNC(lintimeintoutfacedoub1d, LINTIMEINTOUTFACEDOUB1D) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double&,
   const double *, const double *,
   double *);
// in lintimint2d.f:
void F77_FUNC(lintimeintoutfacedoub2d0, LINTIMEINTOUTFACEDOUB2D0) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double&,
   const double *, const double *,
   double *);
void F77_FUNC(lintimeintoutfacedoub2d1, LINTIMEINTOUTFACEDOUB2D1) (const int&,
   const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const double&,
   const double *, const double *,
   double *);
// in lintimint3d.f:
void F77_FUNC(lintimeintoutfacedoub3d0, LINTIMEINTOUTFACEDOUB3D0) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double&,
   const double *, const double *,
   double *);
void F77_FUNC(lintimeintoutfacedoub3d1, LINTIMEINTOUTFACEDOUB3D1) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double&,
   const double *, const double *,
   double *);
void F77_FUNC(lintimeintoutfacedoub3d2, LINTIMEINTOUTFACEDOUB3D2) (const int&,
   const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const double&,
   const double *, const double *,
   double *);
}

namespace SAMRAI {
namespace pdat {

OuterfaceDoubleLinearTimeInterpolateOp::OuterfaceDoubleLinearTimeInterpolateOp():
   hier::TimeInterpolateOperator()
{
}

OuterfaceDoubleLinearTimeInterpolateOp::~OuterfaceDoubleLinearTimeInterpolateOp()
{
}

void
OuterfaceDoubleLinearTimeInterpolateOp::timeInterpolate(
   hier::PatchData& dst_data,
   const hier::Box& where,
   const hier::PatchData& src_data_old,
   const hier::PatchData& src_data_new) const
{
   const tbox::Dimension& dim(where.getDim());

   const OuterfaceData<double>* old_dat =
      dynamic_cast<const OuterfaceData<double> *>(&src_data_old);
   const OuterfaceData<double>* new_dat =
      dynamic_cast<const OuterfaceData<double> *>(&src_data_new);
   OuterfaceData<double>* dst_dat =
      dynamic_cast<OuterfaceData<double> *>(&dst_data);

   TBOX_ASSERT(old_dat != NULL);
   TBOX_ASSERT(new_dat != NULL);
   TBOX_ASSERT(dst_dat != NULL);
   TBOX_ASSERT((where * old_dat->getGhostBox()).isSpatiallyEqual(where));
   TBOX_ASSERT((where * new_dat->getGhostBox()).isSpatiallyEqual(where));
   TBOX_ASSERT((where * dst_dat->getGhostBox()).isSpatiallyEqual(where));
   TBOX_DIM_ASSERT_CHECK_ARGS4(dst_data, where, src_data_old, src_data_new);

   const hier::Index old_ilo = old_dat->getGhostBox().lower();
   const hier::Index old_ihi = old_dat->getGhostBox().upper();
   const hier::Index new_ilo = new_dat->getGhostBox().lower();
   const hier::Index new_ihi = new_dat->getGhostBox().upper();

   const hier::Index dst_ilo = dst_dat->getGhostBox().lower();
   const hier::Index dst_ihi = dst_dat->getGhostBox().upper();

   const hier::Index ifirst = where.lower();
   const hier::Index ilast = where.upper();

   const double old_time = old_dat->getTime();
   const double new_time = new_dat->getTime();
   const double dst_time = dst_dat->getTime();

   TBOX_ASSERT((old_time < dst_time ||
                tbox::MathUtilities<double>::equalEps(old_time, dst_time)) &&
      (dst_time < new_time ||
       tbox::MathUtilities<double>::equalEps(dst_time, new_time)));

   double tfrac = dst_time - old_time;
   double denom = new_time - old_time;
   if (denom > tbox::MathUtilities<double>::getMin()) {
      tfrac /= denom;
   } else {
      tfrac = 0.0;
   }

   for (int d = 0; d < dst_dat->getDepth(); d++) {
      // loop over lower and upper outerface arrays
      for (int i = 0; i < 2; i++) {
         if (dim == tbox::Dimension(1)) {
            F77_FUNC(lintimeintoutfacedoub1d,
               LINTIMEINTOUTFACEDOUB1D) (ifirst(0), ilast(0),
               old_ilo(0), old_ihi(0),
               new_ilo(0), new_ihi(0),
               dst_ilo(0), dst_ihi(0),
               tfrac,
               old_dat->getPointer(0, i, d),
               new_dat->getPointer(0, i, d),
               dst_dat->getPointer(0, i, d));
         } else if (dim == tbox::Dimension(2)) {
            F77_FUNC(lintimeintoutfacedoub2d0,
               LINTIMEINTOUTFACEDOUB2D0) (ifirst(0), ifirst(1), ilast(0),
               ilast(1),
               old_ilo(0), old_ilo(1), old_ihi(0), old_ihi(1),
               new_ilo(0), new_ilo(1), new_ihi(0), new_ihi(1),
               dst_ilo(0), dst_ilo(1), dst_ihi(0), dst_ihi(1),
               tfrac,
               old_dat->getPointer(0, i, d),
               new_dat->getPointer(0, i, d),
               dst_dat->getPointer(0, i, d));
            F77_FUNC(lintimeintoutfacedoub2d1,
               LINTIMEINTOUTFACEDOUB2D1) (ifirst(0), ifirst(1), ilast(0),
               ilast(1),
               old_ilo(0), old_ilo(1), old_ihi(0), old_ihi(1),
               new_ilo(0), new_ilo(1), new_ihi(0), new_ihi(1),
               dst_ilo(0), dst_ilo(1), dst_ihi(0), dst_ihi(1),
               tfrac,
               old_dat->getPointer(1, i, d),
               new_dat->getPointer(1, i, d),
               dst_dat->getPointer(1, i, d));
         } else if (dim == tbox::Dimension(3)) {
            F77_FUNC(lintimeintoutfacedoub3d0,
               LINTIMEINTOUTFACEDOUB3D0) (ifirst(0), ifirst(1), ifirst(2),
               ilast(0), ilast(1), ilast(2),
               old_ilo(0), old_ilo(1), old_ilo(2),
               old_ihi(0), old_ihi(1), old_ihi(2),
               new_ilo(0), new_ilo(1), new_ilo(2),
               new_ihi(0), new_ihi(1), new_ihi(2),
               dst_ilo(0), dst_ilo(1), dst_ilo(2),
               dst_ihi(0), dst_ihi(1), dst_ihi(2),
               tfrac,
               old_dat->getPointer(0, i, d),
               new_dat->getPointer(0, i, d),
               dst_dat->getPointer(0, i, d));
            F77_FUNC(lintimeintoutfacedoub3d1,
               LINTIMEINTOUTFACEDOUB3D1) (ifirst(0), ifirst(1), ifirst(2),
               ilast(0), ilast(1), ilast(2),
               old_ilo(0), old_ilo(1), old_ilo(2),
               old_ihi(0), old_ihi(1), old_ihi(2),
               new_ilo(0), new_ilo(1), new_ilo(2),
               new_ihi(0), new_ihi(1), new_ihi(2),
               dst_ilo(0), dst_ilo(1), dst_ilo(2),
               dst_ihi(0), dst_ihi(1), dst_ihi(2),
               tfrac,
               old_dat->getPointer(1, i, d),
               new_dat->getPointer(1, i, d),
               dst_dat->getPointer(1, i, d));
            F77_FUNC(lintimeintoutfacedoub3d2,
               LINTIMEINTOUTFACEDOUB3D2) (ifirst(0), ifirst(1), ifirst(2),
               ilast(0), ilast(1), ilast(2),
               old_ilo(0), old_ilo(1), old_ilo(2),
               old_ihi(0), old_ihi(1), old_ihi(2),
               new_ilo(0), new_ilo(1), new_ilo(2),
               new_ihi(0), new_ihi(1), new_ihi(2),
               dst_ilo(0), dst_ilo(1), dst_ilo(2),
               dst_ihi(0), dst_ihi(1), dst_ihi(2),
               tfrac,
               old_dat->getPointer(2, i, d),
               new_dat->getPointer(2, i, d),
               dst_dat->getPointer(2, i, d));
         } else {
            TBOX_ERROR(
               "OuterfaceDoubleLinearTimeInterpolateOp::TimeInterpolate dim > 3 not supported"
               << std::endl);
         }
      }
   }
}

}
}
#endif
