/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   A smart pointer template class with RTTI
 *
 ************************************************************************/

#ifndef included_tbox_Dimension_C
#define included_tbox_Dimension_C

#include "SAMRAI/tbox/Dimension.h"

namespace SAMRAI {
namespace tbox {

Dimension::Dimension():
   d_dim(Dimension::getInvalidDimValue())
{
}

Dimension::Dimension(
   const unsigned short& dim):d_dim(dim) {
   TBOX_DIM_ASSERT((!isValid()) ||
      (d_dim > 0 && d_dim <= Dimension::getMaxDimValue()));
}

Dimension::Dimension(
   const Dimension& dimension):d_dim(dimension.d_dim) {
   TBOX_DIM_ASSERT((!isValid()) ||
      (d_dim > 0 && d_dim <= Dimension::getMaxDimValue()));
}

std::ostream&
operator << (
   std::ostream& s,
   const Dimension& dim)
{
   s << dim.getValue() << 'D';
   return s;
}

}
}

#endif
