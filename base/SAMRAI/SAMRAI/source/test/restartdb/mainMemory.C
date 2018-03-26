/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Tests Memory database in SAMRAI
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/DatabaseBox.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"

#include <boost/shared_ptr.hpp>
#include <string>

using namespace std;
using namespace SAMRAI;

#include "database_tests.h"

class RestartTester:public tbox::Serializable
{
public:
   RestartTester()
   {
      tbox::RestartManager::getManager()->registerRestartItem("RestartTester",
         this);
   }

   virtual ~RestartTester() {
   }

   void putToDatabase(
      const boost::shared_ptr<tbox::Database>& db) const
   {
      writeTestData(db);
   }

   void getFromDatabase()
   {
      boost::shared_ptr<tbox::Database> root_db(
         tbox::RestartManager::getManager()->getRootDatabase());

      boost::shared_ptr<tbox::Database> db;
      if (root_db->isDatabase("RestartTester")) {
         db = root_db->getDatabase("RestartTester");
      }

      readTestData(db);
   }

};

int main(
   int argc,
   char* argv[])
{
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      tbox::PIO::logAllNodes("Memorytest.log");

      tbox::plog << "\n--- Memory database tests BEGIN ---" << endl;

      tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();

      RestartTester memory_tester;

      tbox::plog << "\n--- Memory write database tests BEGIN ---" << endl;

      setupTestData();

      boost::shared_ptr<tbox::MemoryDatabase> database(
         new tbox::MemoryDatabase("SAMRAI Restart"));

      restart_manager->setRootDatabase(database);

      restart_manager->writeRestartToDatabase();

      tbox::plog << "\n--- Memory write database tests END ---" << endl;

      tbox::plog << "\n--- Memory read database tests BEGIN ---" << endl;

      // In this test just read the database stored in memory that
      // was just created.
      memory_tester.getFromDatabase();

      database->close();

      tbox::plog << "\n--- Memory read database tests END ---" << endl;

      tbox::plog << "\n--- Memory database tests END ---" << endl;

      if (number_of_failures == 0) {
         tbox::pout << "\nPASSED:  Memory" << endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return number_of_failures;

}
