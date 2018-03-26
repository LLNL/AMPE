/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Abstract factory class for creating patch classes
 *
 ************************************************************************/

#ifndef included_hier_PatchFactory_C
#define included_hier_PatchFactory_C

#include "SAMRAI/hier/PatchFactory.h"

#include <boost/make_shared.hpp>

namespace SAMRAI {
namespace hier {

PatchFactory::PatchFactory()
{
}

PatchFactory::~PatchFactory()
{
}

boost::shared_ptr<Patch>
PatchFactory::allocate(
   const Box& mapped_box_level_mapped_box,
   const boost::shared_ptr<PatchDescriptor>& descriptor) const
{
   return boost::make_shared<Patch>(mapped_box_level_mapped_box, descriptor);
}

}
}
#endif
