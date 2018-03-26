/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
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
Database::putBoolArray(
   const std::string& key,
   const Array<bool>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.getSize() > 0) {
      putBoolArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("Database::putBoolArray() error in database "
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
      Array<bool> local_bool = getBoolArray(key);
      bool* locptr = local_bool.getPointer();
      return locptr == NULL ? defaultvalue : *locptr;
   } else {
      return defaultvalue;
   }
}

void
Database::getBoolArray(
   const std::string& key,
   bool* data,
   const int nelements)
{
   TBOX_ASSERT(!key.empty());

   Array<bool> tmp = getBoolArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getBoolArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
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
 * Create a box array entry in the database with the specified key name.
 *
 *************************************************************************
 */

void
Database::putDatabaseBoxArray(
   const std::string& key,
   const Array<DatabaseBox>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.getSize() > 0) {
      putDatabaseBoxArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("Database::putDatabaseBoxArray() error in database "
         << getName()
         << "\n    Attempt to put zero-length array with key = "
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
      Array<DatabaseBox> local_box = getDatabaseBoxArray(key);
      DatabaseBox* locptr = local_box.getPointer();
      return locptr == NULL ? defaultvalue : *locptr;
   } else {
      return defaultvalue;
   }

}

void
Database::getDatabaseBoxArray(
   const std::string& key,
   DatabaseBox* data,
   const int nelements)
{
   TBOX_ASSERT(!key.empty());

   Array<DatabaseBox> tmp = getDatabaseBoxArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getDatabaseBoxArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
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
 * Create a char array entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
Database::putCharArray(
   const std::string& key,
   const Array<char>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.getSize() > 0) {
      putCharArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("Database::putCharArray() error in database "
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
      Array<char> local_char = getCharArray(key);
      char* locptr = local_char.getPointer();
      return locptr == NULL ? defaultvalue : *locptr;
   } else {
      return defaultvalue;
   }

}

void
Database::getCharArray(
   const std::string& key,
   char* data,
   const int nelements)
{
   TBOX_ASSERT(!key.empty());

   Array<char> tmp = getCharArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getCharArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
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
 * Create a complex array entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
Database::putComplexArray(
   const std::string& key,
   const Array<dcomplex>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.getSize() > 0) {
      putComplexArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("Database::putComplexArray() error in database "
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
      Array<dcomplex> local_dcomplex = getComplexArray(key);
      dcomplex* locptr = local_dcomplex.getPointer();
      if (locptr != NULL) {
         retval = *locptr;
      }
   }
   return retval;
}

void
Database::getComplexArray(
   const std::string& key,
   dcomplex* data,
   const int nelements)
{
   TBOX_ASSERT(!key.empty());

   Array<dcomplex> tmp = getComplexArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getComplexArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
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
 * Create a float array entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
Database::putFloatArray(
   const std::string& key,
   const Array<float>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.getSize() > 0) {
      putFloatArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("Database::putFloatArray() error in database "
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
      Array<float> local_float = getFloatArray(key);
      float* locptr = local_float.getPointer();
      return locptr == NULL ? defaultvalue : *locptr;
   } else {
      return defaultvalue;
   }

}

void
Database::getFloatArray(
   const std::string& key,
   float* data,
   const int nelements)
{
   TBOX_ASSERT(!key.empty());

   Array<float> tmp = getFloatArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getFloatArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
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
 * Create a double array entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
Database::putDoubleArray(
   const std::string& key,
   const Array<double>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.getSize() > 0) {
      putDoubleArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("Database::putDoubleArray() error in database "
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
      Array<double> local_double = getDoubleArray(key);
      double* locptr = local_double.getPointer();
      return locptr == NULL ? defaultvalue : *locptr;
   } else {
      return defaultvalue;
   }
}

void
Database::getDoubleArray(
   const std::string& key,
   double* data,
   const int nelements)
{
   TBOX_ASSERT(!key.empty());

   Array<double> tmp = getDoubleArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getDoubleArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
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
Database::putIntegerArray(
   const std::string& key,
   const Array<int>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.getSize() > 0) {
      putIntegerArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("Database::putIntegerArray() error in database "
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
      Array<int> local_int = getIntegerArray(key);
      int* locptr = local_int.getPointer();
      return locptr == NULL ? defaultvalue : *locptr;
   } else {
      return defaultvalue;
   }

}

void
Database::getIntegerArray(
   const std::string& key,
   int* data,
   const int nelements)
{
   TBOX_ASSERT(!key.empty());

   Array<int> tmp = getIntegerArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getIntegerArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
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
 * Create a string array entry in the database with the specified
 * key name.
 *
 *************************************************************************
 */

void
Database::putStringArray(
   const std::string& key,
   const Array<std::string>& data)
{
   TBOX_ASSERT(!key.empty());

   if (data.getSize() > 0) {
      putStringArray(key, data.getPointer(), data.getSize());
   } else {
      TBOX_ERROR("Database::putStringArray() error in database "
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
      Array<std::string> local_string = getStringArray(key);
      std::string* locptr = local_string.getPointer();
      return locptr == NULL ? defaultvalue : *locptr;
   } else {
      return defaultvalue;
   }

}

void
Database::getStringArray(
   const std::string& key,
   std::string* data,
   const int nelements)
{
   TBOX_ASSERT(!key.empty());

   Array<std::string> tmp = getStringArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tsize) {
      TBOX_ERROR("Database::getStringArray() error in database "
         << getName()
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << std::endl);
   }

   for (int i = 0; i < tsize; i++) {
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
