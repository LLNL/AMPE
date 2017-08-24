/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   A factory for building MemoryDatabases
 *
 ************************************************************************/

#include "SAMRAI/tbox/MemoryDatabaseFactory.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/Utilities.h"

#include "boost/make_shared.hpp"

namespace SAMRAI {
namespace tbox {

MemoryDatabaseFactory::MemoryDatabaseFactory()
{
}

MemoryDatabaseFactory::~MemoryDatabaseFactory()
{
}

MemoryDatabaseFactory::MemoryDatabaseFactory(
   const MemoryDatabaseFactory& other):
   DatabaseFactory()
{
   NULL_USE(other);
}

MemoryDatabaseFactory&
MemoryDatabaseFactory::operator = (
   const MemoryDatabaseFactory& rhs)
{
   NULL_USE(rhs);
   return *this;
}

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
