/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Dimension class for abstracting dimension
 *
 ************************************************************************/
#include "SAMRAI/tbox/Dimension.h"

namespace SAMRAI {
namespace tbox {

Dimension::Dimension(
   const unsigned short& dim):d_dim(dim)
{
   TBOX_DIM_ASSERT(dim > 0 && dim <= SAMRAI::MAX_DIM_VAL);
}

Dimension::Dimension(
   const Dimension& dimension):d_dim(dimension.d_dim)
{
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
