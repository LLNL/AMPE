/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Iterator for array patch data types
 *
 ************************************************************************/

#ifndef included_pdat_ArrayDataIterator_C
#define included_pdat_ArrayDataIterator_C

#include "SAMRAI/pdat/ArrayDataIterator.h"

namespace SAMRAI {
namespace pdat {

ArrayDataIterator::ArrayDataIterator(
   const hier::Box& box,
   bool begin):
   d_index(box.lower()),
   d_box(box)
{
   if (!d_box.empty() && !begin) {
      d_index(d_box.getDim().getValue()-1) =
         d_box.upper(d_box.getDim().getValue()-1) + 1;
   }
}

ArrayDataIterator::ArrayDataIterator(
   const ArrayDataIterator& iter):
   d_index(iter.d_index),
   d_box(iter.d_box)
{
}

ArrayDataIterator::~ArrayDataIterator()
{
}

ArrayDataIterator&
ArrayDataIterator::operator ++ ()
{
   const tbox::Dimension& dim(d_box.getDim());
   d_index(0)++;
   for (int i = 0; i < dim.getValue() - 1; i++) {
      if (d_index(i) > d_box.upper(i)) {
         d_index(i) = d_box.lower(i);
         d_index(i + 1)++;
      } else {
         break;
      }
   }
   return *this;
}

ArrayDataIterator
ArrayDataIterator::operator ++ (
   int)
{
   ArrayDataIterator tmp = *this;
   const tbox::Dimension& dim(d_box.getDim());
   d_index(0)++;
   for (int i = 0; i < dim.getValue() - 1; i++) {
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
