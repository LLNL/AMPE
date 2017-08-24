/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   A factory for building SiloDatabases
 *
 ************************************************************************/

#include "SAMRAI/tbox/SiloDatabaseFactory.h"
#include "SAMRAI/tbox/SiloDatabase.h"
#include "SAMRAI/tbox/Utilities.h"

#include "boost/make_shared.hpp"

namespace SAMRAI {
namespace tbox {

SiloDatabaseFactory::SiloDatabaseFactory()
{
}

SiloDatabaseFactory::~SiloDatabaseFactory()
{
}

SiloDatabaseFactory::SiloDatabaseFactory(
   const SiloDatabaseFactory& other):
   DatabaseFactory()
{
   NULL_USE(other);
}

SiloDatabaseFactory&
SiloDatabaseFactory::operator = (
   const SiloDatabaseFactory& rhs)
{
   NULL_USE(rhs);
   return *this;
}

/**
 * Build a new SiloDatabase object.
 */
boost::shared_ptr<Database>
SiloDatabaseFactory::allocate(
   const std::string& name) {
#ifdef HAVE_SILO
   return boost::make_shared<SiloDatabase>(name);

#else
   NULL_USE(name);
   TBOX_WARNING("Silo5DatabaseFactory: Cannot allocate a SiloDatabase.\n"
      << "SAMRAI was not configured with Silo." << std::endl);
   return boost::shared_ptr<Database>();

#endif
}

}
}
