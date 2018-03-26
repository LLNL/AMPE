/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Misc printing functions in FAC solver test.
 *
 ************************************************************************/
#include "SAMRAI/SAMRAI_config.h"

#include "printObject.h"

#include "SAMRAI/pdat/MDA_Access.h"

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/pdat/ArrayData.h"

#include <boost/shared_ptr.hpp>

using namespace SAMRAI;

/*!
 * \brief Print a box
 */
int printObject(
   std::ostream& os,
   const hier::Box& box,
   const std::string& border,
   unsigned short depth)
{
   NULL_USE(depth);
   const tbox::Dimension& dim = box.getDim();
   if ((dim == tbox::Dimension(1))) {
      os << border << "( " << box.lower(0)
         << " ) to ( "
         << box.upper(0)
         << " )\n";
   }
   if ((dim == tbox::Dimension(2))) {
      os << border << "( " << box.lower(0) << ' ' << box.lower(1)
         << " ) to ( "
         << box.upper(0) << ' ' << box.upper(1)
         << " )\n";
   }
   if ((dim == tbox::Dimension(3))) {
      os << border << "( "
         << box.lower(0) << ' ' << box.lower(1) << ' ' << box.lower(2)
         << " ) to ( "
         << box.upper(0) << ' ' << box.upper(1) << ' ' << box.upper(2)
         << " )\n";
   }
   return 0;
}

/*!
 * \brief Print a patch data object
 */
int printObject(
   std::ostream& os,
   const hier::PatchData& pdat,
   const std::string& border,
   unsigned short depth)
{
   NULL_USE(depth);
   const hier::Box& rbox = pdat.getBox();
   const hier::Box& gbox = pdat.getBox();
   os << border
      << "Box = ( "
      << rbox.lower(0) << ' ' << rbox.lower(1)
      << " ) to ( "
      << rbox.upper(0) << ' ' << rbox.upper(1)
      << " )\n";
   os << border
      << "Ghost hier::Box = ( "
      << gbox.lower(0) << ' ' << gbox.lower(1)
      << " ) to ( "
      << gbox.upper(0) << ' ' << gbox.upper(1)
      << " )\n";
   return 0;
}

/*!
 * \brief Print an array data object
 */
template<class T>
int printObject(
   std::ostream& os,
   const pdat::ArrayData<T>& adat,
   const int depth,
   const std::string& border)
{
   os << border << "printObject(pdat::ArayData)... " << std::endl;
   const hier::Box& rbox = adat.getBox();
   os << border
      << "Box = ( "
      << rbox.lower(0) << ' ' << rbox.lower(1)
      << " ) to ( "
      << rbox.upper(0) << ' ' << rbox.upper(1)
      << " )\n";
   os << border << "depth " << depth << "\n";
   printObject(os,
      adat.getDim(),
      (const T *)adat.getPointer(depth),
      &adat.getBox().lower()[0],
      &adat.getBox().upper()[0]
      );
   return 0;
}
template int printObject<double
                         >(
   std::ostream& os,
   const pdat::ArrayData<double>& adat,
   const int depth,
   const std::string& border);

/*!
 * \brief Print an array
 */
int printObject(
   std::ostream& os
   ,
   const tbox::Dimension& dim
   ,
   const double* a_ptr,
   const int* a_lower,
   const int* a_upper
   ,
   const int* lower
   ,
   const int* upper) {
   switch (dim.getValue()) {
      case 1:
      {
         MDA_AccessConst<double, 1, MDA_OrderColMajor<1> > a(a_ptr,
                                                             a_lower,
                                                             a_upper);
         if (!lower) lower = a_lower;
         if (!upper) upper = a_upper;
         os << "\narray data... " << std::endl;
         for (int j = lower[1]; j <= upper[1]; ++j) {
            for (int i = lower[0]; i <= upper[0]; ++i) {
               os << i << ' ' << j << ' ' << a(i, j) << std::endl;
            }
         }
      }
      break;
      case 2:
      {
         MDA_AccessConst<double, 2, MDA_OrderColMajor<2> > a(a_ptr,
                                                             a_lower,
                                                             a_upper);
         if (!lower) lower = a_lower;
         if (!upper) upper = a_upper;
         os << "\narray data... " << std::endl;
         for (int j = lower[1]; j <= upper[1]; ++j) {
            for (int i = lower[0]; i <= upper[0]; ++i) {
               os << i << ' ' << j << ' ' << a(i, j) << std::endl;
            }
         }
      }
      break;
      case 3:
      {
         MDA_AccessConst<double, 3, MDA_OrderColMajor<3> > a(a_ptr,
                                                             a_lower,
                                                             a_upper);
         if (!lower) lower = a_lower;
         if (!upper) upper = a_upper;
         os << "\narray data... " << std::endl;
         for (int j = lower[1]; j <= upper[1]; ++j) {
            for (int i = lower[0]; i <= upper[0]; ++i) {
               os << i << ' ' << j << ' ' << a(i, j) << std::endl;
            }
         }
      }
      break;
   }
   return 0;
}
