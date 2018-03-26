/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Iterator for side centered patch data types
 *
 ************************************************************************/

#ifndef included_pdat_SideIterator_C
#define included_pdat_SideIterator_C

#include "SAMRAI/pdat/SideIterator.h"

namespace SAMRAI {
namespace pdat {

SideIterator::SideIterator(
   const hier::Box& box,
   const int axis,
   bool begin):
   d_index(box.lower(), axis, SideIndex::Lower),
   d_box(SideGeometry::toSideBox(box, axis))
{
   if (!d_box.empty() && !begin) {
      d_index(d_box.getDim().getValue()-1) =
         d_box.upper(d_box.getDim().getValue()-1) + 1;
   }
}

SideIterator::SideIterator(
   const SideIterator& iter):
   d_index(iter.d_index),
   d_box(iter.d_box)
{
}

SideIterator::~SideIterator()
{
}

SideIterator&
SideIterator::operator ++ ()
{
   d_index(0)++;
   for (int i = 0; i < d_box.getDim().getValue() - 1; i++) {
      if (d_index(i) > d_box.upper(i)) {
         d_index(i) = d_box.lower(i);
         d_index(i + 1)++;
      } else {
         break;
      }
   }
   return *this;
}

SideIterator
SideIterator::operator ++ (
   int)
{
   SideIterator tmp = *this;
   d_index(0)++;
   for (int i = 0; i < d_box.getDim().getValue() - 1; i++) {
      if (d_index(i) > d_box.upper(i)) {
         d_index(i) = d_box.lower(i);
         d_index(i + 1)++;
      } else {
         break;
      }
   }
   return tmp;
}

}
}
#endif
