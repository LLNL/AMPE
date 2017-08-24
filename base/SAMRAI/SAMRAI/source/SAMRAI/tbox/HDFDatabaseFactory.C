/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   An abstract base class for a HDFDatabaseFactory
 *
 ************************************************************************/

#include "SAMRAI/tbox/HDFDatabaseFactory.h"
#include "SAMRAI/tbox/HDFDatabase.h"
#include "SAMRAI/tbox/Utilities.h"

#include "boost/make_shared.hpp"

namespace SAMRAI {
namespace tbox {

HDFDatabaseFactory::HDFDatabaseFactory()
{
}

HDFDatabaseFactory::~HDFDatabaseFactory()
{
}

HDFDatabaseFactory::HDFDatabaseFactory(
   const HDFDatabaseFactory& other):
   DatabaseFactory()
{
   NULL_USE(other);
}

HDFDatabaseFactory&
HDFDatabaseFactory::operator = (
   const HDFDatabaseFactory& rhs)
{
   NULL_USE(rhs);
   return *this;
}

/**
 * Build a new Database object.
 */
boost::shared_ptr<Database>
HDFDatabaseFactory::allocate(
   const std::string& name) {
#ifdef HAVE_HDF5
   boost::shared_ptr<HDFDatabase> database(
      boost::make_shared<HDFDatabase>(name));
   return database;

#else
   NULL_USE(name);
   TBOX_WARNING("HDF5DatabaseFactory: Cannot allocate an HDFDatabase.\n"
      << "SAMRAI was not configured with HDF." << std::endl);
   return boost::shared_ptr<Database>();

#endif
}

}
}
