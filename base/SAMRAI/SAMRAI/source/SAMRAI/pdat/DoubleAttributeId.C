/**********************************************************************
*
* This file is part of the SAMRAI distribution.  For full copyright
* information, see COPYRIGHT and COPYING.LESSER.
*
* Copyright:     (c) 1997 - 2011 Lawrence Livermore National Security, LLC
* Description:   pdat
**********************************************************************/
#ifndef included_pdat_DoubleAttributeId_C
#define included_pdat_DoubleAttributeId_C

#include "SAMRAI/pdat/DoubleAttributeId.h"

namespace SAMRAI {
namespace pdat {

/**********************************************************************
 * c'tor
 *********************************************************************/
DoubleAttributeId::DoubleAttributeId(
   int value):
   d_val(value)
{
}

/**********************************************************************
 * DoubleAttributeId c'tor (private)
 *********************************************************************/
DoubleAttributeId::DoubleAttributeId():
   d_val(-1)
{
}

/**********************************************************************
 * copy c'tor
 *********************************************************************/
DoubleAttributeId::DoubleAttributeId(
   const DoubleAttributeId& other):
   d_val(other.d_val)
{
}

/**********************************************************************
 * d'tor
 *********************************************************************/
DoubleAttributeId::~DoubleAttributeId()
{
}

}
}
#endif
