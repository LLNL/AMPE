/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Tests HDF database in SAMRAI
 *
 ************************************************************************/

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/DatabaseBox.h"
#include "SAMRAI/tbox/Complex.h"
#include "SAMRAI/tbox/HDFDatabase.h"
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
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   /*
    * Create block to force pointer deallocation.  If this is not done
    * then there will be memory leaks reported.
    */
   {

      tbox::PIO::logAllNodes("HDF5test.log");

#ifdef HAVE_HDF5

      tbox::plog << "\n--- HDF5 database tests BEGIN ---" << endl;

      tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();

      RestartTester hdf_tester;

      tbox::plog << "\n--- HDF5 write database tests BEGIN ---" << endl;

      setupTestData();

      restart_manager->writeRestartFile("test_dir", 0);

      tbox::plog << "\n--- HDF5 write database tests END ---" << endl;

      tbox::plog << "\n--- HDF5 read database tests BEGIN ---" << endl;

      restart_manager->closeRestartFile();

      restart_manager->openRestartFile("test_dir",
         0,
         mpi.getSize());

      hdf_tester.getFromDatabase();

      restart_manager->closeRestartFile();

      tbox::plog << "\n--- HDF5 read database tests END ---" << endl;

      tbox::plog << "\n--- HDF5 database tests END ---" << endl;

#endif

      if (number_of_failures == 0) {
         tbox::pout << "\nPASSED:  HDF5" << endl;
      }
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return number_of_failures;

}
