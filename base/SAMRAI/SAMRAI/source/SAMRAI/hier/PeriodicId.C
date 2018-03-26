/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Periodic shift identifier in periodic domain.
 *
 ************************************************************************/

#ifndef included_hier_PeriodicId_C
#define included_hier_PeriodicId_C

#include "SAMRAI/hier/PeriodicId.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

const PeriodicId PeriodicId::s_invalid_id(-1);
const PeriodicId PeriodicId::s_zero_id(0);

/*
 ******************************************************************************
 ******************************************************************************
 */
PeriodicId::PeriodicId():
   d_value(invalidId().d_value) {
}

/*
 ******************************************************************************
 ******************************************************************************
 */
PeriodicId::PeriodicId(
   const PeriodicId& other):
   d_value(other.d_value) {
}

/*
 ******************************************************************************
 ******************************************************************************
 */
PeriodicId::PeriodicId(
   const int& value):
   d_value(value) {
}

/*
 ******************************************************************************
 ******************************************************************************
 */
PeriodicId::~PeriodicId()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   d_value = s_invalid_id.d_value;
#endif
}

/*
 ******************************************************************************
 ******************************************************************************
 */
std::ostream&
operator << (
   std::ostream& co,
   const PeriodicId& r)
{
   co << r.d_value;
   return co;
}

}
}
#endif
