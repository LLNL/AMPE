/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   A factory for building MemoryDatabases
 *
 ************************************************************************/

#ifndef included_tbox_MemoryDatabaseFactory
#define included_tbox_MemoryDatabaseFactory

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/DatabaseFactory.h"

namespace SAMRAI {
namespace tbox {

/**
 * @brief MemoryDatabase factory.
 *
 * Builds a new MemoryDatabase.
 */
class MemoryDatabaseFactory:public DatabaseFactory
{
   /**
    * Build a new Database object.
    */
   virtual boost::shared_ptr<Database>
   allocate(
      const std::string& name);
};

}
}

#endif
