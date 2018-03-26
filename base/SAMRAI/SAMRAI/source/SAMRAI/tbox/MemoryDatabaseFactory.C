/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   A factory for building MemoryDatabases
 *
 ************************************************************************/

#include "SAMRAI/tbox/MemoryDatabaseFactory.h"
#include "SAMRAI/tbox/MemoryDatabase.h"

#include <boost/make_shared.hpp>

namespace SAMRAI {
namespace tbox {

/**
 * Build a new MemoryDatabase object.
 */
boost::shared_ptr<Database>
MemoryDatabaseFactory::allocate(
   const std::string& name) {
   return boost::make_shared<MemoryDatabase>(name);
}

}
}
