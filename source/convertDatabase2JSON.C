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

void copyDatabase(std::shared_ptr<tbox::Database> database,
                  pt::ptree& new_database)
{
   std::vector<std::string> keys(database->getAllKeys());

   for (std::vector<std::string>::const_iterator k_itr = keys.begin();
        k_itr != keys.end(); ++k_itr) {

      const std::string& key = *k_itr;
      tbox::Database::DataType my_type = database->getArrayType(key);
      size_t size = database->getArraySize(key);

      if (my_type == tbox::Database::SAMRAI_DATABASE) {
         std::shared_ptr<tbox::Database> child_db = database->getDatabase(key);
         pt::ptree new_db;
         copyDatabase(child_db, new_db);
         new_database.add_child(key, new_db);
      } else if (my_type == tbox::Database::SAMRAI_CHAR) {
         if (size == 1) {
            new_database.put(key, database->getChar(key));
         } else {
            std::vector<char> bvec(database->getCharVector(key));
            pt::ptree nodes;
            for (auto& vec : bvec) {
               pt::ptree node;
               node.put("", vec);
               nodes.push_back(std::make_pair("", node));
            }
            new_database.add_child(key, nodes);
         }
      } else if (my_type == tbox::Database::SAMRAI_INT) {
         if (size == 1) {
            new_database.put(key, database->getInteger(key));
         } else {
            std::vector<int> bvec(database->getIntegerVector(key));
            pt::ptree nodes;
            for (auto& vec : bvec) {
               pt::ptree node;
               node.put("", vec);
               nodes.push_back(std::make_pair("", node));
            }
            new_database.add_child(key, nodes);
         }
      } else if (my_type == tbox::Database::SAMRAI_DOUBLE) {
         if (size == 1) {
            new_database.put(key, database->getDouble(key));
         } else {
            std::vector<double> bvec(database->getDoubleVector(key));
            pt::ptree nodes;
            for (auto& vec : bvec) {
               pt::ptree node;
               node.put("", vec);
               nodes.push_back(std::make_pair("", node));
            }
            new_database.add_child(key, nodes);
         }
      } else if (my_type == tbox::Database::SAMRAI_FLOAT) {
         if (size == 1) {
            new_database.put(key, database->getFloat(key));
         } else {
            std::vector<float> bvec(database->getFloatVector(key));
            pt::ptree nodes;
            for (auto& vec : bvec) {
               pt::ptree node;
               node.put("", vec);
               nodes.push_back(std::make_pair("", node));
            }
            new_database.add_child(key, nodes);
         }
      } else if (my_type == tbox::Database::SAMRAI_STRING) {
         if (size == 1) {
            new_database.put(key, database->getString(key));
         } else {
            std::vector<std::string> bvec(database->getStringVector(key));
            pt::ptree nodes;
            for (auto& vec : bvec) {
               pt::ptree node;
               node.put("", vec);
               nodes.push_back(std::make_pair("", node));
            }
            new_database.add_child(key, nodes);
         }
      }
   }
}

void copyDatabase(std::shared_ptr<tbox::Database> database,
                  std::shared_ptr<tbox::Database> new_database)
{
   std::vector<std::string> keys(database->getAllKeys());

   for (std::vector<std::string>::const_iterator k_itr = keys.begin();
        k_itr != keys.end(); ++k_itr) {

      const std::string& key = *k_itr;
      tbox::Database::DataType my_type = database->getArrayType(key);
      size_t size = database->getArraySize(key);

      if (my_type == tbox::Database::SAMRAI_DATABASE) {
         std::shared_ptr<tbox::Database> child_db = database->getDatabase(key);
         std::shared_ptr<tbox::Database> new_db =
             new_database->putDatabase(key);
         copyDatabase(child_db, new_db);
      } else if (my_type == tbox::Database::SAMRAI_BOOL) {
         if (size == 1) {
            new_database->putBool(key, database->getBool(key));
         } else {
            std::vector<bool> bvec(database->getBoolVector(key));
            new_database->putBoolVector(key, bvec);
         }
      } else if (my_type == tbox::Database::SAMRAI_CHAR) {
         if (size == 1) {
            new_database->putChar(key, database->getChar(key));
         } else {
            std::vector<char> bvec(database->getCharVector(key));
            new_database->putCharVector(key, bvec);
         }
      } else if (my_type == tbox::Database::SAMRAI_INT) {
         if (size == 1) {
            new_database->putInteger(key, database->getInteger(key));
         } else {
            std::vector<int> bvec(database->getIntegerVector(key));
            new_database->putIntegerVector(key, bvec);
         }
      } else if (my_type == tbox::Database::SAMRAI_DOUBLE) {
         if (size == 1) {
            new_database->putDouble(key, database->getDouble(key));
         } else {
            std::vector<double> bvec(database->getDoubleVector(key));
            new_database->putDoubleVector(key, bvec);
         }
      } else if (my_type == tbox::Database::SAMRAI_FLOAT) {
         if (size == 1) {
            new_database->putFloat(key, database->getFloat(key));
         } else {
            std::vector<float> bvec(database->getFloatVector(key));
            new_database->putFloatVector(key, bvec);
         }
      } else if (my_type == tbox::Database::SAMRAI_STRING) {
         if (size == 1) {
            new_database->putString(key, database->getString(key));
         } else {
            std::vector<std::string> bvec(database->getStringVector(key));
            new_database->putStringVector(key, bvec);
         }
      }
   }
}


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
