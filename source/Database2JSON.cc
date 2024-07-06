// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/MemoryDatabase.h"

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <string>
#include <sstream>

using namespace SAMRAI;
namespace pt = boost::property_tree;

void copyDatabase(std::shared_ptr<tbox::Database> database,
                  pt::ptree& new_database)
{
   assert(database);

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
            double value = database->getDouble(key);
            std::stringstream ss;
            ss << std::setprecision(10);
            ss << value;
            new_database.put(key, ss.str());
         } else {
            std::vector<double> bvec(database->getDoubleVector(key));
            std::vector<std::string> sbvec;
            for (auto& value : bvec) {
               std::stringstream ss;
               ss << std::setprecision(10);
               ss << value;
               sbvec.push_back(ss.str());
            }
            pt::ptree nodes;
            for (auto& vec : sbvec) {
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
