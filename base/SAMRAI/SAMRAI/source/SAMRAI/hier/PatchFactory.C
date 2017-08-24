/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Abstract factory class for creating patch classes
 *
 ************************************************************************/
#include "SAMRAI/hier/PatchFactory.h"

#include "boost/make_shared.hpp"

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
   const Box& box_level_box,
   const boost::shared_ptr<PatchDescriptor>& descriptor) const
{
   return boost::make_shared<Patch>(box_level_box, descriptor);
}

}
}
