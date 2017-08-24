/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   An abstract base class for a HDFDatabaseFactory
 *
 ************************************************************************/

#ifndef included_tbox_HDFDatabaseFactory
#define included_tbox_HDFDatabaseFactory

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/DatabaseFactory.h"

namespace SAMRAI {
namespace tbox {

/**
 * @brief HDFDatabase factory.
 *
 * Builds a new HDFDatabase.
 */
class HDFDatabaseFactory:public DatabaseFactory
{
public:
   /**
    * Default constructor.
    */
   HDFDatabaseFactory();

   /**
    * Copy constructor.
    */
   HDFDatabaseFactory(
      const HDFDatabaseFactory& other);

   /**
    * Assignment operator.
    */
   HDFDatabaseFactory&
   operator = (
      const HDFDatabaseFactory& rhs);

   /**
    * Destructor.
    */
   ~HDFDatabaseFactory();

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
