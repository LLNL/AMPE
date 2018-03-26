/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   A box structure representing a portion of the AMR index space
 *
 ************************************************************************/

#include "SAMRAI/tbox/DatabaseBox.h"

namespace SAMRAI {
namespace tbox {

DatabaseBox::DatabaseBox():
   d_dim(Dimension::getInvalidDimension())
{
   d_data.d_dimension = 0;
   d_data.d_lo[0] = d_data.d_hi[0] = 0;
   d_data.d_lo[1] = d_data.d_hi[1] = 0;
   d_data.d_lo[2] = d_data.d_hi[2] = 0;
}

DatabaseBox::DatabaseBox(
   const DatabaseBox& box):
   d_dim(box.d_dim)
{
   d_data.d_dimension = box.d_data.d_dimension;
   d_data.d_lo[0] = box.d_data.d_lo[0];
   d_data.d_lo[1] = box.d_data.d_lo[1];
   d_data.d_lo[2] = box.d_data.d_lo[2];
   d_data.d_hi[0] = box.d_data.d_hi[0];
   d_data.d_hi[1] = box.d_data.d_hi[1];
   d_data.d_hi[2] = box.d_data.d_hi[2];
}

DatabaseBox::DatabaseBox(
   const Dimension& dim,
   const int* lower,
   const int* upper):
   d_dim(dim)
{
   TBOX_ASSERT(dim.getValue() <= DatabaseBox_MAX_DIM);

   d_data.d_dimension = dim.getValue();
   for (int i = 0; i < d_data.d_dimension; i++) {
      d_data.d_lo[i] = lower[i];
      d_data.d_hi[i] = upper[i];
   }
   for (int j = d_data.d_dimension; j < DatabaseBox_MAX_DIM; j++) {
      d_data.d_lo[j] = 0;
      d_data.d_hi[j] = 0;
   }
}

DatabaseBox::~DatabaseBox()
{
}

DatabaseBox&
DatabaseBox::operator = (
   const DatabaseBox& box)
{
   d_dim = box.d_dim;
   d_data.d_dimension = box.d_data.d_dimension;
   d_data.d_lo[0] = box.d_data.d_lo[0];
   d_data.d_lo[1] = box.d_data.d_lo[1];
   d_data.d_lo[2] = box.d_data.d_lo[2];
   d_data.d_hi[0] = box.d_data.d_hi[0];
   d_data.d_hi[1] = box.d_data.d_hi[1];
   d_data.d_hi[2] = box.d_data.d_hi[2];
   return *this;
}

bool
DatabaseBox::empty() const
{
   bool is_empty = (d_data.d_dimension == 0 ? true : false);
   for (int i = 0; i < d_data.d_dimension; i++) {
      if (d_data.d_hi[i] < d_data.d_lo[i]) is_empty = true;
   }
   return is_empty;
}

int
DatabaseBox::operator == (
   const DatabaseBox& box) const
{
   bool equals = (d_data.d_dimension == box.d_data.d_dimension);
   for (int i = 0; i < d_data.d_dimension; i++) {
      if (d_data.d_lo[i] != box.d_data.d_lo[i]) {
         equals = false;
      }
      if (d_data.d_hi[i] != box.d_data.d_hi[i]) {
         equals = false;
      }
   }
   return equals;
}

}
}
