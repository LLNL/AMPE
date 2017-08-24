/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   Tests Silo database in SAMRAI
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/DatabaseBox.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/SiloDatabase.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"

#include "boost/shared_ptr.hpp"
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

   void putToRestart(
      const boost::shared_ptr<tbox::Database>& db) const
   {
      writeTestData(db);
   }

   void getFromRestart()
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
      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

      tbox::PIO::logAllNodes("Silotest.log");

#ifdef HAVE_SILO

      tbox::plog << "\n--- Silo database tests BEGIN ---" << endl;

      tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();

      RestartTester silo_tester;

      tbox::plog << "\n--- Silo write database tests BEGIN ---" << endl;

      setupTestData();

      boost::shared_ptr<tbox::SiloDatabase> database(
         new tbox::SiloDatabase("SAMRAI Restart"));
      std::string name = "./restart." + tbox::Utilities::processorToString(
            mpi.getRank()) + ".silo";
      DBfile* silo_file = DBCreate(
            name.c_str(), DB_CLOBBER, DB_LOCAL, 0, DB_PDB);
      database->attachToFile(silo_file, "/");

      restart_manager->setRootDatabase(database);

      restart_manager->writeRestartToDatabase();

      database->close();
      restart_manager->setRootDatabase(boost::shared_ptr<tbox::Database>());

      DBClose(silo_file);

      tbox::plog << "\n--- Silo write database tests END ---" << endl;

      tbox::plog << "\n--- Silo read database tests BEGIN ---" << endl;

      database.reset(new tbox::SiloDatabase("SAMRAI Restart"));
      silo_file = DBOpen(name.c_str(), DB_UNKNOWN, DB_READ);
      database->attachToFile(silo_file, "/");

      restart_manager->setRootDatabase(database);

      silo_tester.getFromRestart();

      database->close();
      restart_manager->setRootDatabase(boost::shared_ptr<tbox::Database>());

      DBClose(silo_file);

      tbox::plog << "\n--- Silo read database tests END ---" << endl;

      tbox::plog << "\n--- Silo database tests END ---" << endl;

#endif

      if (number_of_failures == 0) {
         tbox::pout << "\nPASSED:  Silo" << endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return number_of_failures;

}
