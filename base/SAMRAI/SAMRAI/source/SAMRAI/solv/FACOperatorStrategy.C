/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Interface to user-defined operations used in FAC solve.
 *
 ************************************************************************/

#ifndef included_solv_FACOperatorStrategy_C
#define included_solv_FACOperatorStrategy_C

#include "SAMRAI/solv/FACOperatorStrategy.h"

namespace SAMRAI {
namespace solv {

FACOperatorStrategy::FACOperatorStrategy()
{
}

FACOperatorStrategy::~FACOperatorStrategy()
{
}

void
FACOperatorStrategy::deallocateOperatorState()
{
}

}
}
#endif
