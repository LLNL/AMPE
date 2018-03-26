/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Interface to templated operations for real data on hierarchy.
 *
 ************************************************************************/

#ifndef included_math_HierarchyDataOpsReal_C
#define included_math_HierarchyDataOpsReal_C

#include "SAMRAI/math/HierarchyDataOpsReal.h"

namespace SAMRAI {
namespace math {

template<class TYPE>
HierarchyDataOpsReal<TYPE>::HierarchyDataOpsReal()
{
}

template<class TYPE>
HierarchyDataOpsReal<TYPE>::~HierarchyDataOpsReal()
{
}

/*
 *************************************************************************
 *
 * The following are private and cannot be used, but they are defined
 * here for compilers that require that every template declaration have
 * a definition (a stupid requirement, if you ask me).
 *
 *************************************************************************
 */

template<class TYPE>
HierarchyDataOpsReal<TYPE>::HierarchyDataOpsReal(
   const HierarchyDataOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

template<class TYPE>
void
HierarchyDataOpsReal<TYPE>::operator = (
   const HierarchyDataOpsReal<TYPE>& foo)
{
   NULL_USE(foo);
}

}
}
#endif
