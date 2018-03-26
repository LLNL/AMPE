/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Factory abstract base class for creating patch data objects
 *
 ************************************************************************/

#ifndef included_hier_PatchDataFactory_C
#define included_hier_PatchDataFactory_C

#include "SAMRAI/hier/PatchDataFactory.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace hier {

PatchDataFactory::PatchDataFactory(
   const IntVector& ghosts):
   d_ghosts(ghosts)
{
   TBOX_ASSERT(ghosts.min() >= 0);
}

PatchDataFactory::~PatchDataFactory()
{
}

/**********************************************************************
* Default implementation
**********************************************************************/

MultiblockDataTranslator *
PatchDataFactory::getMultiblockDataTranslator()
{
   return (MultiblockDataTranslator *)NULL;
}

}
}
#endif
