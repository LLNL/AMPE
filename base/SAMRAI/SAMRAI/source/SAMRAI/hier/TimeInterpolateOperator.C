/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Abstract base class for time interpolation operators.
 *
 ************************************************************************/

#ifndef included_hier_TimeInterpolateOperator_C
#define included_hier_TimeInterpolateOperator_C

#include "SAMRAI/hier/TimeInterpolateOperator.h"

namespace SAMRAI {
namespace hier {

TimeInterpolateOperator::TimeInterpolateOperator(
   const std::string& name) :
   d_name(name)
{
}

TimeInterpolateOperator::~TimeInterpolateOperator()
{
}

}
}
#endif
