#include "Database2JSON.h"

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MemoryDatabase.h"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <string>

using namespace SAMRAI;
namespace pt = boost::property_tree;


int main(int argc, char* argv[])
{
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   {
      std::string calphad_filename = argv[1];

      std::shared_ptr<tbox::MemoryDatabase> database(
          new tbox::MemoryDatabase("database"));
      tbox::InputManager::getManager()->parseInputFile(calphad_filename,
                                                       database);

      //      std::shared_ptr<tbox::MemoryDatabase> new_database(
      //          new tbox::MemoryDatabase("new_database"));
      //      copyDatabase(database, new_database);
      //      new_database->printClassData(std::cout);

      // copy MemoryDatabase into a boost property_tree
      pt::ptree troot;
      copyDatabase(database, troot);
      pt::write_json(std::cout, troot);
   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return 0;
}
