/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Simple Cartesian grid geometry for an AMR hierarchy.
 *
 ************************************************************************/

#ifndef included_geom_CartesianPatchGeometry_C
#define included_geom_CartesianPatchGeometry_C

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace geom {

// using namespace std;

/*
 *************************************************************************
 *
 * Constructor for CartesianPatchGeometry allocates and sets
 * patch coordinate system information.
 *
 *************************************************************************
 */
CartesianPatchGeometry::CartesianPatchGeometry(
   const hier::IntVector& ratio_to_level_zero,
   const TwoDimBool& touches_regular_bdry,
   const TwoDimBool& touches_periodic_bdry,
   const double* dx,
   const double* x_lo,
   const double* x_up):
   hier::PatchGeometry(ratio_to_level_zero,
                       touches_regular_bdry,
                       touches_periodic_bdry)
{
   TBOX_ASSERT(!(dx == (double *)NULL));
   TBOX_ASSERT(!(x_lo == (double *)NULL));
   TBOX_ASSERT(!(x_up == (double *)NULL));

   const tbox::Dimension& dim(ratio_to_level_zero.getDim());

   for (int id = 0; id < dim.getValue(); id++) {
      d_dx[id] = dx[id];
      d_x_lo[id] = x_lo[id];
      d_x_up[id] = x_up[id];
   }
}

/*
 *************************************************************************
 *
 * Destructor for CartesianPatchGeometry deallocates dx array.
 *
 *************************************************************************
 */
CartesianPatchGeometry::~CartesianPatchGeometry()
{
}

/*
 *************************************************************************
 *
 * Print CartesianPatchGeometry class data.
 *
 *************************************************************************
 */
void
CartesianPatchGeometry::printClassData(
   std::ostream& os) const
{
   const tbox::Dimension& dim(getRatio().getDim());

   os << "Printing CartesianPatchGeometry data: this = "
      << (CartesianPatchGeometry *)this << std::endl;
   os << "x_lo = ";
   for (int id1 = 0; id1 < dim.getValue(); id1++) {
      os << d_x_lo[id1] << "   ";
   }
   os << std::endl;
   os << "x_up = ";
   for (int id2 = 0; id2 < dim.getValue(); id2++) {
      os << d_x_up[id2] << "   ";
   }
   os << std::endl;
   os << "dx = ";
   for (int id3 = 0; id3 < dim.getValue(); id3++) {
      os << d_dx[id3] << "   ";
   }
   os << std::endl;

   hier::PatchGeometry::printClassData(os);
}

}
}
#endif
