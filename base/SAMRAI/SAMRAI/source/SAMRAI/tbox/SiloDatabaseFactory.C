/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   A factory for building SiloDatabases
 *
 ************************************************************************/

#include "SAMRAI/tbox/SiloDatabaseFactory.h"
#include "SAMRAI/tbox/SiloDatabase.h"

#include <boost/make_shared.hpp>

namespace SAMRAI {
namespace tbox {

#ifdef HAVE_SILO
/**
 * Build a new SiloDatabase object.
 */
boost::shared_ptr<Database>
SiloDatabaseFactory::allocate(
   const std::string& name) {
#ifdef HAVE_SILO
   return boost::make_shared<SiloDatabase>(name);

#else
   return NULL;

#endif
}
#endif

}
}
