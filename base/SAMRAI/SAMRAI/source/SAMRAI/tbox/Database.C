/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   An abstract base class for the SAMRAI database objects
 *
 ************************************************************************/

#include "SAMRAI/tbox/Database.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace tbox {

Database::Database()
{
}

Database::~Database()
{
}

/*
 ************************************************************************
 *
 * Get database entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a database type.
 *
 ************************************************************************
 */

boost::shared_ptr<Database>
Database::getDatabaseWithDefault(
   const std::string& key,
   const boost::shared_ptr<Database>& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   if (keyExists(key)) {
      return getDatabase(key);
   } else {
      return defaultvalue;
   }

}

/*
 * Boolean
 */

/*
 * Create a boolean scalar entry in the database with the specified
 * key name.  A scalar entry is an array of one.
 *
 */

void
Database::putBool(
   const std::string& key,
   const bool& data)
{
   TBOX_ASSERT(!key.empty());

   putBoolArray(key, &data, 1);
}

/*
 * Create a boolean array entry in the database with the specified
 * key
 */

void
Database::putBoolVector(
   const std::string& key,
   const std::vector<bool>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.size() > 0) {
      int nbools = static_cast<int>(data.size());
      bool* bool_array = new bool[nbools];
      for (int i = 0; i < nbools; ++i) {
         bool_array[i] = data[i];
      }
      putBoolArray(key, bool_array, nbools);
      delete[] bool_array;
   } else {
      TBOX_ERROR("Database::putBoolVector() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 * Get boolean scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a boolean type.
 */

bool
Database::getBool(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   bool ret_val;
   getBoolArray(key, &ret_val, 1);

   return ret_val;
}

/*
 * Get boolean scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a boolean type.
 */

bool
Database::getBoolWithDefault(
   const std::string& key,
   const bool& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   if (keyExists(key)) {
      std::vector<bool> local_bool = getBoolVector(key);
      return local_bool.empty() ? defaultvalue : local_bool[0];
   } else {
      return defaultvalue;
   }
}

void
Database::getBoolArray(
   const std::string& key,
   bool* data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());

   std::vector<bool> tmp = getBoolVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getBoolVector() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }

}

/*
 *************************************************************************
 *
 * Create a box entry in the database with the specified
 * key name.  A box entry is an array of one.
 *
 *************************************************************************
 */

void
Database::putDatabaseBox(
   const std::string& key,
   const DatabaseBox& data)
{
   TBOX_ASSERT(!key.empty());

   putDatabaseBoxArray(key, &data, 1);
}

/*
 *************************************************************************
 *
 * Create a box vector entry in the database with the specified key name.
 *
 *************************************************************************
 */

void
Database::putDatabaseBoxVector(
   const std::string& key,
   const std::vector<DatabaseBox>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.size() > 0) {
      putDatabaseBoxArray(key, &data[0], static_cast<int>(data.size()));
   } else {
      TBOX_ERROR("Database::putDatabaseBoxVector() error in database "
         << getName()
         << "\n    Attempt to put zero-length vector with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Get box scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a box type.
 *
 ************************************************************************
 */

DatabaseBox
Database::getDatabaseBox(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   DatabaseBox ret_val;
   getDatabaseBoxArray(key, &ret_val, 1);

   return ret_val;
}

/*
 ************************************************************************
 *
 * Get box scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a box type.
 *
 ************************************************************************
 */

DatabaseBox
Database::getDatabaseBoxWithDefault(
   const std::string& key,
   const DatabaseBox& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   if (keyExists(key)) {
      std::vector<DatabaseBox> local_box = getDatabaseBoxVector(key);
      return local_box.empty() ? defaultvalue : local_box[0];
   } else {
      return defaultvalue;
   }

}

void
Database::getDatabaseBoxArray(
   const std::string& key,
   DatabaseBox* data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());

   std::vector<DatabaseBox> tmp = getDatabaseBoxVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getDatabaseBoxArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 * Char
 */

/*
 *************************************************************************
 *
 * Create a char scalar entry in the database with the specified
 * key name.  A scalar entry is an array of one.
 *
 *************************************************************************
 */

void
Database::putChar(
   const std::string& key,
   const char& data)
{
   TBOX_ASSERT(!key.empty());

   putCharArray(key, &data, 1);

}

/*
 *************************************************************************
 *
 * Create a char vector entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
Database::putCharVector(
   const std::string& key,
   const std::vector<char>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.size() > 0) {
      putCharArray(key, &data[0], static_cast<int>(data.size()));
   } else {
      TBOX_ERROR("Database::putCharVector() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Get char scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a char type.
 *
 ************************************************************************
 */

char
Database::getChar(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   char ret_val;
   getCharArray(key, &ret_val, 1);

   return ret_val;
}

/*
 ************************************************************************
 *
 * Get char scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a char type.
 *
 ************************************************************************
 */

char
Database::getCharWithDefault(
   const std::string& key,
   const char& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   if (keyExists(key)) {
      std::vector<char> local_char = getCharVector(key);
      return local_char.empty() ? defaultvalue : local_char[0];
   } else {
      return defaultvalue;
   }

}

void
Database::getCharArray(
   const std::string& key,
   char* data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());

   std::vector<char> tmp = getCharVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getCharArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 * Complex
 */

/*
 *************************************************************************
 *
 * Create a complex scalar entry in the database with the specified
 * key name.  A scalar entry is an array of one.
 *
 *************************************************************************
 */

void
Database::putComplex(
   const std::string& key,
   const dcomplex& data)
{
   TBOX_ASSERT(!key.empty());

   putComplexArray(key, &data, 1);
}

/*
 *************************************************************************
 *
 * Create a complex vector entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
Database::putComplexVector(
   const std::string& key,
   const std::vector<dcomplex>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.size() > 0) {
      putComplexArray(key, &data[0], static_cast<int>(data.size()));
   } else {
      TBOX_ERROR("Database::putComplexVector() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Get complex scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a complex type.
 *
 ************************************************************************
 */

dcomplex
Database::getComplex(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   dcomplex ret_val;
   getComplexArray(key, &ret_val, 1);

   return ret_val;
}

/*
 ************************************************************************
 *
 * Get complex scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a complex type.
 *
 ************************************************************************
 */

dcomplex
Database::getComplexWithDefault(
   const std::string& key,
   const dcomplex& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   dcomplex retval = defaultvalue;

   if (keyExists(key)) {
      std::vector<dcomplex> local_dcomplex = getComplexVector(key);
      if (!local_dcomplex.empty()) {
         retval = local_dcomplex[0];
      }
   }
   return retval;
}

void
Database::getComplexArray(
   const std::string& key,
   dcomplex* data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());

   std::vector<dcomplex> tmp = getComplexVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getComplexArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 * Float
 */

/*
 *************************************************************************
 *
 * Create a float scalar entry in the database with the specified
 * key name.  A scalar entry is an array of one.
 *
 *************************************************************************
 */

void
Database::putFloat(
   const std::string& key,
   const float& data)
{
   TBOX_ASSERT(!key.empty());

   putFloatArray(key, &data, 1);
}

/*
 *************************************************************************
 *
 * Create a float vector entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
Database::putFloatVector(
   const std::string& key,
   const std::vector<float>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.size() > 0) {
      putFloatArray(key, &data[0], static_cast<int>(data.size()));
   } else {
      TBOX_ERROR("Database::putFloatVector() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }

}

/*
 ************************************************************************
 *
 * Get float scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a float type.
 *
 ************************************************************************
 */

float
Database::getFloat(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   float ret_val;
   getFloatArray(key, &ret_val, 1);

   return ret_val;
}

/*
 ************************************************************************
 *
 * Get float scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a float type.
 *
 ************************************************************************
 */

float
Database::getFloatWithDefault(
   const std::string& key,
   const float& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   if (keyExists(key)) {
      std::vector<float> local_float = getFloatVector(key);
      return local_float.empty() ? defaultvalue : local_float[0];
   } else {
      return defaultvalue;
   }

}

void
Database::getFloatArray(
   const std::string& key,
   float* data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());

   std::vector<float> tmp = getFloatVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getFloatArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 * Double
 */

/*
 *************************************************************************
 *
 * Create a double scalar entry in the database with the specified
 * key name.  A scalar entry is an array of one.
 *
 *************************************************************************
 */

void
Database::putDouble(
   const std::string& key,
   const double& data)
{
   TBOX_ASSERT(!key.empty());

   putDoubleArray(key, &data, 1);
}

/*
 *************************************************************************
 *
 * Create a double vector entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
Database::putDoubleVector(
   const std::string& key,
   const std::vector<double>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.size() > 0) {
      putDoubleArray(key, &data[0], static_cast<int>(data.size()));
   } else {
      TBOX_ERROR("Database::putDoubleVector() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Get double scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a double type.
 *
 ************************************************************************
 */

double
Database::getDouble(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   double ret_val;
   getDoubleArray(key, &ret_val, 1);

   return ret_val;
}

/*
 ************************************************************************
 *
 * Get double scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a double type.
 *
 ************************************************************************
 */

double
Database::getDoubleWithDefault(
   const std::string& key,
   const double& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   if (keyExists(key)) {
      std::vector<double> local_double = getDoubleVector(key);
      return local_double.empty() ? defaultvalue : local_double[0];
   } else {
      return defaultvalue;
   }
}

void
Database::getDoubleArray(
   const std::string& key,
   double* data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());

   std::vector<double> tmp = getDoubleVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getDoubleArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 * Integer
 */

/*
 *************************************************************************
 *
 * Create a integer scalar entry in the database with the specified
 * key name.  A scalar entry is an array of one.
 *
 *************************************************************************
 */

void
Database::putInteger(
   const std::string& key,
   const int& data)
{
   TBOX_ASSERT(!key.empty());

   putIntegerArray(key, &data, 1);
}

/*
 *************************************************************************
 *
 * Create an integer array entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
Database::putIntegerVector(
   const std::string& key,
   const std::vector<int>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.size() > 0) {
      putIntegerArray(key, &data[0], static_cast<int>(data.size()));
   } else {
      TBOX_ERROR("Database::putIntegerVector() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Get integer scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a integer type.
 *
 ************************************************************************
 */

int
Database::getInteger(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   int ret_val;
   getIntegerArray(key, &ret_val, 1);

   return ret_val;
}

/*
 ************************************************************************
 *
 * Get integer scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a integer type.
 *
 ************************************************************************
 */

int
Database::getIntegerWithDefault(
   const std::string& key,
   const int& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   if (keyExists(key)) {
      std::vector<int> local_int = getIntegerVector(key);
      return local_int.empty() ? defaultvalue : local_int[0];
   } else {
      return defaultvalue;
   }

}

void
Database::getIntegerArray(
   const std::string& key,
   int* data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());

   std::vector<int> tmp = getIntegerVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getIntegerArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }
}

/*
 * String
 */

/*
 *************************************************************************
 *
 * Create a string scalar entry in the database with the specified
 * key name.  A scalar entry is an array of one.
 *
 *************************************************************************
 */

void
Database::putString(
   const std::string& key,
   const std::string& data)
{
   TBOX_ASSERT(!key.empty());

   putStringArray(key, &data, 1);
}

/*
 *************************************************************************
 *
 * Create a string vector entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
Database::putStringVector(
   const std::string& key,
   const std::vector<std::string>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.size() > 0) {
      putStringArray(key, &data[0], static_cast<int>(data.size()));
   } else {
      TBOX_ERROR("Database::putStringVector() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
         << key << std::endl);
   }
}

/*
 ************************************************************************
 *
 * Get string scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a string type.
 *
 ************************************************************************
 */

std::string
Database::getString(
   const std::string& key)
{
   TBOX_ASSERT(!key.empty());

   std::string ret_val;
   getStringArray(key, &ret_val, 1);

   return ret_val;
}

/*
 ************************************************************************
 *
 * Get string scalar entry from the database with the specified key
 * name. An error message is printed and the program exits if the
 * specified key does not exist in the database or is not associated
 * with a string type.
 *
 ************************************************************************
 */

std::string
Database::getStringWithDefault(
   const std::string& key,
   const std::string& defaultvalue)
{
   TBOX_ASSERT(!key.empty());

   if (keyExists(key)) {
      std::vector<std::string> local_string = getStringVector(key);
      return local_string.empty() ? defaultvalue : local_string[0];
   } else {
      return defaultvalue;
   }

}

void
Database::getStringArray(
   const std::string& key,
   std::string* data,
   const size_t nelements)
{
   TBOX_ASSERT(!key.empty());

   std::vector<std::string> tmp = getStringVector(key);
   const size_t tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getStringArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (size_t i = 0; i < tsize; ++i) {
      data[i] = tmp[i];
   }

}

bool
Database::isVector(
   const std::string& key)
{
   return isInteger(key + "_size");
}

}
}
