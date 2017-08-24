/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Generic identifier used on a single process.
 *
 ************************************************************************/
#include "SAMRAI/hier/LocalId.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <iostream>

namespace SAMRAI {
namespace hier {

const LocalId
LocalId::s_invalid_id(
   tbox::MathUtilities<int>::getMax());
const LocalId LocalId::s_zero_id(0);

/*
 *******************************************************************************
 *******************************************************************************
 */
LocalId::LocalId():
   d_value(getInvalidId().d_value) {
}

/*
 *******************************************************************************
 *******************************************************************************
 */
LocalId::LocalId(
   const LocalId& other):
   d_value(other.d_value) {
}

/*
 *******************************************************************************
 *******************************************************************************
 */
LocalId::LocalId(
   const int& value):
   d_value(value) {
}

/*
 *******************************************************************************
 *******************************************************************************
 */
LocalId::~LocalId()
{
#ifdef DEBUG_CHECK_ASSERTIONS
   d_value = s_invalid_id.d_value;
#endif
}

/*
 *******************************************************************************
 *******************************************************************************
 */
std::ostream&
operator << (
   std::ostream& co,
   const LocalId& r)
{
   if (r.isValid()) {
      co << r.d_value;
   } else {
      co << 'X';
   }
   return co;
}

}
}
