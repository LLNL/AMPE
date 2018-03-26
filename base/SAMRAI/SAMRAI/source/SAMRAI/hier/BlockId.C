/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Block identifier in multiblock domain.
 *
 ************************************************************************/

#ifndef included_hier_BlockId_C
#define included_hier_BlockId_C

#include "SAMRAI/hier/BlockId.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

const BlockId BlockId::s_invalid_id(tbox::MathUtilities<int>::getMax());
const BlockId BlockId::s_zero_id(0);

/*
 *******************************************************************************
 *******************************************************************************
 */
BlockId::BlockId():
   d_value(invalidId().d_value)
{
}

/*
 *******************************************************************************
 *******************************************************************************
 */
BlockId::BlockId(
   const BlockId& other):
   d_value(other.d_value)
{
}

/*
 *******************************************************************************
 *******************************************************************************
 */
BlockId::BlockId(
   const int& value):
   d_value(value)
{
}

/*
 *******************************************************************************
 *******************************************************************************
 */
BlockId::~BlockId()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   d_value = s_invalid_id.d_value;
#endif
}

}
}
#endif
