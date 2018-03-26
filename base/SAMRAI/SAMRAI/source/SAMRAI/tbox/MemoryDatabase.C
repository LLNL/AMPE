/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   An memory database structure that stores (key,value) pairs in memory
 *
 ************************************************************************/

#include "SAMRAI/tbox/MemoryDatabase.h"

#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/IOStream.h"

#include <stdlib.h>

#include "SAMRAI/tbox/SAMRAI_MPI.h"

#define MEMORY_DB_ERROR(X) \
   do {                                         \
      pout << "MemoryDatabase: " << X << std::endl << std::flush;       \
      printClassData(pout);                                             \
      pout << "Program abort called..." << std::endl << std::flush;     \
      SAMRAI_MPI::abort();                                              \
   } while (0)

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Suppress XLC warnings
 * o */
#pragma report(disable, CPPC5334)
#pragma report(disable, CPPC5328)
#endif

namespace SAMRAI {
namespace tbox {

const int MemoryDatabase::PRINT_DEFAULT = 1;
const int MemoryDatabase::PRINT_INPUT = 2;
const int MemoryDatabase::PRINT_UNUSED = 4;
const int MemoryDatabase::SSTREAM_BUFFER = 4096;

MemoryDatabase::MemoryDatabase(
   const std::string& name) :
   d_database_name(name)
{
}

/*
 *************************************************************************
 *
 * The virtual destructor deallocates database data.
 *
 *************************************************************************
 */

MemoryDatabase::~MemoryDatabase()
{
}

/*
 *************************************************************************
 *
 * Create memory data file specified by name.
 *
 *************************************************************************
 */

bool
MemoryDatabase::create(
   const std::string& name)
{
   d_database_name = name;
   d_keyvalues.clear();

   return true;
}

/*
 *************************************************************************
 *
 * Open memory data file specified by name
 *
 *************************************************************************
 */

bool
MemoryDatabase::open(
   const std::string& name,
   const bool read_write_mode)
{
   if (read_write_mode == false) {
      TBOX_ERROR("MemoryDatabase::open: MemoryDatabase only supports\n"
         << "read-write mode.  The read_write_mode flag must be true.");

   }
   d_database_name = name;
   d_keyvalues.clear();

   return true;
}

/*
 *************************************************************************
 *
 * Close the open data file.
 *
 *************************************************************************
 */

bool
MemoryDatabase::close()
{
   d_database_name = "";
   d_keyvalues.clear();

   return true;
}

/*
 *************************************************************************
 *
 * Return whether the key exists in the database.
 *
 *************************************************************************
 */

bool
MemoryDatabase::keyExists(
   const std::string& key)
{
   return findKeyData(key) ? true : false;
}

/*
 *************************************************************************
 *
 * Return all of the keys in the database.
 *
 *************************************************************************
 */

Array<std::string>
MemoryDatabase::getAllKeys()
{
   const int n = static_cast<int>(d_keyvalues.size());
   Array<std::string> keys(n);

   int k = 0;
   for (std::list<KeyData>::iterator i = d_keyvalues.begin();
         i != d_keyvalues.end(); i++) {
      keys[k++] = i->d_key;
   }

   return keys;
}

/*
 *************************************************************************
 *
 * Get the type of the array entry associated with the specified key
 *
 *************************************************************************
 */
enum Database::DataType
MemoryDatabase::getArrayType(
   const std::string& key) {
   KeyData* keydata = findKeyData(key);

   if (keydata) {
      return keydata->d_type;
   } else {
      return Database::SAMRAI_INVALID;
   }
}

/*
 *************************************************************************
 *
 * Get the size of the array entry associated with the specified key;
 * return 0 if the key does not exist.
 *
 *************************************************************************
 */

int
MemoryDatabase::getArraySize(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   if (keydata && keydata->d_type != Database::SAMRAI_DATABASE) {
      return keydata->d_array_size;
   } else {
      return 0;
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage the database values within the database.
 *
 *************************************************************************
 */

bool
MemoryDatabase::isDatabase(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return keydata ? keydata->d_type == Database::SAMRAI_DATABASE : false;
}

boost::shared_ptr<Database>
MemoryDatabase::putDatabase(
   const std::string& key)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key = key;
   keydata.d_type = Database::SAMRAI_DATABASE;
   keydata.d_array_size = 1;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_database.reset(new MemoryDatabase(key));
   d_keyvalues.push_back(keydata);
   return keydata.d_database;
}

boost::shared_ptr<Database>
MemoryDatabase::getDatabase(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if (keydata->d_type != Database::SAMRAI_DATABASE) {
      MEMORY_DB_ERROR("Key=" << key << " is not a database...");
   }
   keydata->d_accessed = true;
   return keydata->d_database;
}

/*
 *************************************************************************
 *
 * Member functions that manage boolean values within the database.
 *
 *************************************************************************
 */

bool
MemoryDatabase::isBool(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return keydata ? keydata->d_type == Database::SAMRAI_BOOL : false;
}

void
MemoryDatabase::putBool(
   const std::string& key,
   const bool& data)
{
   putBoolArray(key, &data, 1);
}

void
MemoryDatabase::putBoolArray(
   const std::string& key,
   const Array<bool>& data)
{
   putBoolArray(key, data.getPointer(), data.getSize());
}

void
MemoryDatabase::putBoolArray(
   const std::string& key,
   const bool * const data,
   const int nelements)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key = key;
   keydata.d_type = Database::SAMRAI_BOOL;
   keydata.d_array_size = nelements;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_boolean = Array<bool>(nelements);

   for (int i = 0; i < nelements; i++) {
      keydata.d_boolean[i] = data[i];
   }

   d_keyvalues.push_back(keydata);
}

bool
MemoryDatabase::getBool(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if ((keydata->d_type != Database::SAMRAI_BOOL) ||
       (keydata->d_array_size != 1)) {
      MEMORY_DB_ERROR("Key=" << key << " is not a boolean scalar...");
   }
   keydata->d_accessed = true;
   return keydata->d_boolean[0];
}

bool
MemoryDatabase::getBoolWithDefault(
   const std::string& key,
   const bool& defaultvalue)
{
   KeyData* keydata = findKeyData(key);
   if (keydata) return getBool(key);

   putBool(key, defaultvalue);
   d_keyvalues.back().d_from_default = true;
   return defaultvalue;
}

Array<bool>
MemoryDatabase::getBoolArray(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if (keydata->d_type != Database::SAMRAI_BOOL) {
      MEMORY_DB_ERROR("Key=" << key << " is not a boolean...");
   }
   keydata->d_accessed = true;
   return keydata->d_boolean;
}

void
MemoryDatabase::getBoolArray(
   const std::string& key,
   bool* data,
   const int nelements)
{
   Array<bool> tmp = getBoolArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tmp.getSize()) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage box values within the database.
 *
 *************************************************************************
 */

bool
MemoryDatabase::isDatabaseBox(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return keydata ? keydata->d_type == Database::SAMRAI_BOX : false;
}

void
MemoryDatabase::putDatabaseBox(
   const std::string& key,
   const DatabaseBox& data)
{
   putDatabaseBoxArray(key, &data, 1);
}

void
MemoryDatabase::putDatabaseBoxArray(
   const std::string& key,
   const Array<DatabaseBox>& data)
{
   putDatabaseBoxArray(key, data.getPointer(), data.getSize());
}

void
MemoryDatabase::putDatabaseBoxArray(
   const std::string& key,
   const DatabaseBox * const data,
   const int nelements)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key = key;
   keydata.d_type = Database::SAMRAI_BOX;
   keydata.d_array_size = nelements;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_box = Array<DatabaseBox>(nelements);

   for (int i = 0; i < nelements; i++) {
      keydata.d_box[i] = data[i];
   }

   d_keyvalues.push_back(keydata);
}

DatabaseBox
MemoryDatabase::getDatabaseBox(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if ((keydata->d_type != Database::SAMRAI_BOX) ||
       (keydata->d_array_size != 1)) {
      MEMORY_DB_ERROR("Key=" << key << " is not a single box...");
   }
   keydata->d_accessed = true;
   return keydata->d_box[0];
}

DatabaseBox
MemoryDatabase::getDatabaseBoxWithDefault(
   const std::string& key,
   const DatabaseBox& defaultvalue)
{
   KeyData* keydata = findKeyData(key);
   if (keydata) return getDatabaseBox(key);

   putDatabaseBox(key, defaultvalue);
   d_keyvalues.back().d_from_default = true;
   return defaultvalue;
}

Array<DatabaseBox>
MemoryDatabase::getDatabaseBoxArray(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if (keydata->d_type != Database::SAMRAI_BOX) {
      MEMORY_DB_ERROR("Key=" << key << " is not a box...");
   }
   keydata->d_accessed = true;
   return keydata->d_box;
}

void
MemoryDatabase::getDatabaseBoxArray(
   const std::string& key,
   DatabaseBox* data,
   const int nelements)
{
   Array<DatabaseBox> tmp = getDatabaseBoxArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tmp.getSize()) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage character values within the database.
 *
 *************************************************************************
 */

bool
MemoryDatabase::isChar(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return keydata ? keydata->d_type == Database::SAMRAI_CHAR : false;
}

void
MemoryDatabase::putChar(
   const std::string& key,
   const char& data)
{
   putCharArray(key, &data, 1);
}

void
MemoryDatabase::putCharArray(
   const std::string& key,
   const Array<char>& data)
{
   putCharArray(key, data.getPointer(), data.getSize());
}

void
MemoryDatabase::putCharArray(
   const std::string& key,
   const char * const data,
   const int nelements)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key = key;
   keydata.d_type = Database::SAMRAI_CHAR;
   keydata.d_array_size = nelements;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_char = Array<char>(nelements);

   for (int i = 0; i < nelements; i++) {
      keydata.d_char[i] = data[i];
   }

   d_keyvalues.push_back(keydata);
}

char
MemoryDatabase::getChar(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if ((keydata->d_type != Database::SAMRAI_CHAR) ||
       (keydata->d_array_size != 1)) {
      MEMORY_DB_ERROR("Key=" << key << " is not a single character...");
   }
   keydata->d_accessed = true;
   return keydata->d_char[0];
}

char
MemoryDatabase::getCharWithDefault(
   const std::string& key,
   const char& defaultvalue)
{
   KeyData* keydata = findKeyData(key);
   if (keydata) return getChar(key);

   putChar(key, defaultvalue);
   d_keyvalues.back().d_from_default = true;
   return defaultvalue;
}

Array<char>
MemoryDatabase::getCharArray(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if (keydata->d_type != Database::SAMRAI_CHAR) {
      MEMORY_DB_ERROR("Key=" << key << " is not a character...");
   }
   keydata->d_accessed = true;
   return keydata->d_char;
}

void
MemoryDatabase::getCharArray(
   const std::string& key,
   char* data,
   const int nelements)
{
   Array<char> tmp = getCharArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tmp.getSize()) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage complex values within the database.
 * Note that complex numbers may be promoted from integers, floats,
 * and doubles.
 *
 *************************************************************************
 */

bool
MemoryDatabase::isComplex(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return !keydata ? false : (keydata->d_type == Database::SAMRAI_COMPLEX
                              || keydata->d_type == Database::SAMRAI_INT
                              || keydata->d_type == Database::SAMRAI_FLOAT
                              || keydata->d_type == Database::SAMRAI_DOUBLE);
}

void
MemoryDatabase::putComplex(
   const std::string& key,
   const dcomplex& data)
{
   putComplexArray(key, &data, 1);
}

void
MemoryDatabase::putComplexArray(
   const std::string& key,
   const Array<dcomplex>& data)
{
   putComplexArray(key, data.getPointer(), data.getSize());
}

void
MemoryDatabase::putComplexArray(
   const std::string& key,
   const dcomplex * const data,
   const int nelements)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key = key;
   keydata.d_type = Database::SAMRAI_COMPLEX;
   keydata.d_array_size = nelements;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_complex = Array<dcomplex>(nelements);

   for (int i = 0; i < nelements; i++) {
      keydata.d_complex[i] = data[i];
   }

   d_keyvalues.push_back(keydata);
}

dcomplex
MemoryDatabase::getComplex(
   const std::string& key)
{
   dcomplex value(0.0, 0.0);
   KeyData* keydata = findKeyDataOrExit(key);

   if (keydata->d_array_size != 1) {
      MEMORY_DB_ERROR("Key=" << key << " is not a single complex...");
   }

   switch (keydata->d_type) {
      case Database::SAMRAI_INT:
         value = dcomplex((double)keydata->d_integer[0], 0.0);
         break;
      case Database::SAMRAI_FLOAT:
         value = dcomplex((double)keydata->d_float[0], 0.0);
         break;
      case Database::SAMRAI_DOUBLE:
         value = dcomplex(keydata->d_double[0], 0.0);
         break;
      case Database::SAMRAI_COMPLEX:
         value = keydata->d_complex[0];
         break;
      default:
         MEMORY_DB_ERROR("Key=" << key << " is not a single complex...");
   }

   keydata->d_accessed = true;
   return value;
}

dcomplex
MemoryDatabase::getComplexWithDefault(
   const std::string& key,
   const dcomplex& defaultvalue)
{
   KeyData* keydata = findKeyData(key);
   if (keydata) return getComplex(key);

   putComplex(key, defaultvalue);
   d_keyvalues.back().d_from_default = true;
   return defaultvalue;
}

Array<dcomplex>
MemoryDatabase::getComplexArray(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   Array<dcomplex> array;
   switch (keydata->d_type) {
      case Database::SAMRAI_INT: {
         array = Array<dcomplex>(keydata->d_integer.getSize());
         for (int i = 0; i < keydata->d_integer.getSize(); i++) {
            array[i] = dcomplex((double)keydata->d_integer[i], 0.0);
         }
         break;
      }
      case Database::SAMRAI_FLOAT: {
         array = Array<dcomplex>(keydata->d_float.getSize());
         for (int i = 0; i < keydata->d_float.getSize(); i++) {
            array[i] = dcomplex((double)keydata->d_float[i], 0.0);
         }
         break;
      }
      case Database::SAMRAI_DOUBLE: {
         array = Array<dcomplex>(keydata->d_double.getSize());
         for (int i = 0; i < keydata->d_float.getSize(); i++) {
            array[i] = dcomplex(keydata->d_double[i], 0.0);
         }
         break;
      }
      case Database::SAMRAI_COMPLEX:
         array = keydata->d_complex;
         break;
      default:
         MEMORY_DB_ERROR("Key=" << key << " is not a complex...");
   }
   keydata->d_accessed = true;
   return array;
}

void
MemoryDatabase::getComplexArray(
   const std::string& key,
   dcomplex* data,
   const int nelements)
{
   Array<dcomplex> tmp = getComplexArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tmp.getSize()) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage double values within the database.
 * Note that doubles may be promoted from integers or floats.
 *
 *************************************************************************
 */

bool
MemoryDatabase::isDouble(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return !keydata ? false : (keydata->d_type == Database::SAMRAI_DOUBLE
                              || keydata->d_type == Database::SAMRAI_INT
                              || keydata->d_type == Database::SAMRAI_FLOAT);
}

void
MemoryDatabase::putDouble(
   const std::string& key,
   const double& data)
{
   putDoubleArray(key, &data, 1);
}

void
MemoryDatabase::putDoubleArray(
   const std::string& key,
   const Array<double>& data)
{
   putDoubleArray(key, data.getPointer(), data.getSize());
}

void
MemoryDatabase::putDoubleArray(
   const std::string& key,
   const double * const data,
   const int nelements)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key = key;
   keydata.d_type = Database::SAMRAI_DOUBLE;
   keydata.d_array_size = nelements;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_double = Array<double>(nelements);

   for (int i = 0; i < nelements; i++) {
      keydata.d_double[i] = data[i];
   }

   d_keyvalues.push_back(keydata);
}

double
MemoryDatabase::getDouble(
   const std::string& key)
{
   double value = 0.0;
   KeyData* keydata = findKeyDataOrExit(key);

   if (keydata->d_array_size != 1) {
      MEMORY_DB_ERROR("Key=" << key << " is not a single double...");
   }

   switch (keydata->d_type) {
      case Database::SAMRAI_INT:
         value = (double)keydata->d_integer[0];
         break;
      case Database::SAMRAI_FLOAT:
         value = (double)keydata->d_float[0];
         break;
      case Database::SAMRAI_DOUBLE:
         value = keydata->d_double[0];
         break;
      default:
         MEMORY_DB_ERROR("Key=" << key << " is not a single double...");
   }

   keydata->d_accessed = true;
   return value;
}

double
MemoryDatabase::getDoubleWithDefault(
   const std::string& key,
   const double& defaultvalue)
{
   KeyData* keydata = findKeyData(key);
   if (keydata) return getDouble(key);

   putDouble(key, defaultvalue);
   d_keyvalues.back().d_from_default = true;
   return defaultvalue;
}

Array<double>
MemoryDatabase::getDoubleArray(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   Array<double> array;
   switch (keydata->d_type) {
      case Database::SAMRAI_INT: {
         array = Array<double>(keydata->d_integer.getSize());
         for (int i = 0; i < keydata->d_integer.getSize(); i++) {
            array[i] = (double)keydata->d_integer[i];
         }
         break;
      }
      case Database::SAMRAI_FLOAT: {
         array = Array<double>(keydata->d_float.getSize());
         for (int i = 0; i < keydata->d_float.getSize(); i++) {
            array[i] = (double)keydata->d_float[i];
         }
         break;
      }
      case Database::SAMRAI_DOUBLE: {
         array = keydata->d_double;
         break;
      }
      default:
         MEMORY_DB_ERROR("Key=" << key << " is not a double...");
   }
   keydata->d_accessed = true;
   return array;
}

void
MemoryDatabase::getDoubleArray(
   const std::string& key,
   double* data,
   const int nelements)
{
   Array<double> tmp = getDoubleArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tmp.getSize()) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage float values within the database.
 * Note that floats may be promoted from integers or truncated from
 * doubles (without a warning).
 *
 *************************************************************************
 */

bool
MemoryDatabase::isFloat(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return !keydata ? false : (keydata->d_type == Database::SAMRAI_DOUBLE
                              || keydata->d_type == Database::SAMRAI_INT
                              || keydata->d_type == Database::SAMRAI_FLOAT);
}

void
MemoryDatabase::putFloat(
   const std::string& key,
   const float& data)
{
   putFloatArray(key, &data, 1);
}

void
MemoryDatabase::putFloatArray(
   const std::string& key,
   const Array<float>& data)
{
   putFloatArray(key, data.getPointer(), data.getSize());
}

void
MemoryDatabase::putFloatArray(
   const std::string& key,
   const float * const data,
   const int nelements)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key = key;
   keydata.d_type = Database::SAMRAI_FLOAT;
   keydata.d_array_size = nelements;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_float = Array<float>(nelements);

   for (int i = 0; i < nelements; i++) {
      keydata.d_float[i] = data[i];
   }

   d_keyvalues.push_back(keydata);
}

float
MemoryDatabase::getFloat(
   const std::string& key)
{

// Disable Intel warning about conversions
#ifdef __INTEL_COMPILER
#pragma warning (disable:810)
#endif

   float value = 0.0;
   KeyData* keydata = findKeyDataOrExit(key);

   if (keydata->d_array_size != 1) {
      MEMORY_DB_ERROR("Key=" << key << " is not a single float...");
   }

   switch (keydata->d_type) {
      case Database::SAMRAI_INT:
         value = static_cast<float>(keydata->d_integer[0]);
         break;
      case Database::SAMRAI_FLOAT:
         value = keydata->d_float[0];
         break;
      case Database::SAMRAI_DOUBLE:
         value = static_cast<float>(keydata->d_double[0]);
         break;
      default:
         MEMORY_DB_ERROR("Key=" << key << " is not a single float...");
   }

   keydata->d_accessed = true;
   return value;
}

float
MemoryDatabase::getFloatWithDefault(
   const std::string& key,
   const float& defaultvalue)
{
   KeyData* keydata = findKeyData(key);
   if (keydata) return getFloat(key);

   putFloat(key, defaultvalue);
   d_keyvalues.back().d_from_default = true;
   return defaultvalue;
}

Array<float>
MemoryDatabase::getFloatArray(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   Array<float> array;
   switch (keydata->d_type) {
      case Database::SAMRAI_INT: {
         array = Array<float>(keydata->d_integer.getSize());
         for (int i = 0; i < keydata->d_integer.getSize(); i++) {
            array[i] = static_cast<float>(keydata->d_integer[i]);
         }
         break;
      }
      case Database::SAMRAI_FLOAT:
         array = keydata->d_float;
         break;
      case Database::SAMRAI_DOUBLE: {
         array = Array<float>(keydata->d_double.getSize());
         for (int i = 0; i < keydata->d_double.getSize(); i++) {
            array[i] = static_cast<float>(keydata->d_double[i]);
         }
         break;
      }
      default:
         MEMORY_DB_ERROR("Key=" << key << " is not a float...");
   }
   keydata->d_accessed = true;
   return array;
}

void
MemoryDatabase::getFloatArray(
   const std::string& key,
   float* data,
   const int nelements)
{
   Array<float> tmp = getFloatArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tmp.getSize()) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage integer values within the database.
 *
 *************************************************************************
 */

bool
MemoryDatabase::isInteger(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return !keydata ? false : keydata->d_type == Database::SAMRAI_INT;
}

void
MemoryDatabase::putInteger(
   const std::string& key,
   const int& data)
{
   putIntegerArray(key, &data, 1);
}

void
MemoryDatabase::putIntegerArray(
   const std::string& key,
   const Array<int>& data)
{
   putIntegerArray(key, data.getPointer(), data.getSize());
}

void
MemoryDatabase::putIntegerArray(
   const std::string& key,
   const int * const data,
   const int nelements)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key = key;
   keydata.d_type = Database::SAMRAI_INT;
   keydata.d_array_size = nelements;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_integer = Array<int>(nelements);

   for (int i = 0; i < nelements; i++) {
      keydata.d_integer[i] = data[i];
   }

   d_keyvalues.push_back(keydata);
}

int
MemoryDatabase::getInteger(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if ((keydata->d_type != Database::SAMRAI_INT) ||
       (keydata->d_array_size != 1)) {
      MEMORY_DB_ERROR("Key=" << key << " is not an integer scalar...");
   }
   keydata->d_accessed = true;
   return keydata->d_integer[0];
}

int
MemoryDatabase::getIntegerWithDefault(
   const std::string& key,
   const int& defaultvalue)
{
   KeyData* keydata = findKeyData(key);
   if (keydata) return getInteger(key);

   putInteger(key, defaultvalue);
   d_keyvalues.back().d_from_default = true;
   return defaultvalue;
}

Array<int>
MemoryDatabase::getIntegerArray(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if (keydata->d_type != Database::SAMRAI_INT) {
      MEMORY_DB_ERROR("Key=" << key << " is not an integer...");
   }
   keydata->d_accessed = true;
   return keydata->d_integer;
}

void
MemoryDatabase::getIntegerArray(
   const std::string& key,
   int* data,
   const int nelements)
{
   Array<int> tmp = getIntegerArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tmp.getSize()) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
 *************************************************************************
 *
 * Member functions that manage string values within the database.
 *
 *************************************************************************
 */

bool
MemoryDatabase::isString(
   const std::string& key)
{
   KeyData* keydata = findKeyData(key);
   return !keydata ? false : keydata->d_type == Database::SAMRAI_STRING;
}

void
MemoryDatabase::putString(
   const std::string& key,
   const std::string& data)
{
   putStringArray(key, &data, 1);
}

void
MemoryDatabase::putStringArray(
   const std::string& key,
   const Array<std::string>& data)
{
   putStringArray(key, data.getPointer(), data.getSize());
}

void
MemoryDatabase::putStringArray(
   const std::string& key,
   const std::string * const data,
   const int nelements)
{
   deleteKeyIfFound(key);
   KeyData keydata;
   keydata.d_key = key;
   keydata.d_type = Database::SAMRAI_STRING;
   keydata.d_array_size = nelements;
   keydata.d_accessed = false;
   keydata.d_from_default = false;
   keydata.d_string = Array<std::string>(nelements);

   for (int i = 0; i < nelements; i++) {
      keydata.d_string[i] = data[i];
   }

   d_keyvalues.push_back(keydata);
}

std::string
MemoryDatabase::getString(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if ((keydata->d_type != Database::SAMRAI_STRING) ||
       (keydata->d_array_size != 1)) {
      MEMORY_DB_ERROR("Key=" << key << " is not a single string...");
   }
   keydata->d_accessed = true;
   return keydata->d_string[0];
}

std::string
MemoryDatabase::getStringWithDefault(
   const std::string& key,
   const std::string& defaultvalue)
{
   KeyData* keydata = findKeyData(key);
   if (keydata) return getString(key);

   putString(key, defaultvalue);
   d_keyvalues.back().d_from_default = true;
   return defaultvalue;
}

Array<std::string>
MemoryDatabase::getStringArray(
   const std::string& key)
{
   KeyData* keydata = findKeyDataOrExit(key);
   if (keydata->d_type != Database::SAMRAI_STRING) {
      MEMORY_DB_ERROR("Key=" << key << " is not a string...");
   }
   keydata->d_accessed = true;
   return keydata->d_string;
}

void
MemoryDatabase::getStringArray(
   const std::string& key,
   std::string* data,
   const int nelements)
{
   Array<std::string> tmp = getStringArray(key);
   const int tsize = tmp.getSize();

   if (nelements != tmp.getSize()) {
      MEMORY_DB_ERROR(
         "Incorrect array size=" << nelements << " specified for key="
                                 << key << " with array size="
                                 << tsize << "...");
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

std::string
MemoryDatabase::getName()
{
   return d_database_name;
}

std::string
MemoryDatabase::getName() const
{
   return d_database_name;
}

/*
 *************************************************************************
 *
 * Search the current database for a matching key.  If found, delete
 * that key and return true.  If the key does not exist, then return
 * false.
 *
 *************************************************************************
 */

bool MemoryDatabase::deleteKeyIfFound(
   const std::string& key)
{
   for (std::list<KeyData>::iterator i = d_keyvalues.begin();
        i != d_keyvalues.end(); i++) {
      if (i->d_key == key) {
         d_keyvalues.erase(i);
         return true;
      }
   }
   return false;
}

/*
 *************************************************************************
 *
 * Find the key data associated with the specified key and return a
 * pointer to the record.  If no such key data exists, then return NULL.
 *
 *************************************************************************
 */

MemoryDatabase::KeyData *
MemoryDatabase::findKeyData(
   const std::string& key)
{
   for (std::list<KeyData>::iterator i = d_keyvalues.begin();
        i != d_keyvalues.end(); i++) {
      if (key == i->d_key) {
         return &(*i);
      }
   }
   return NULL;
}

/*
 *************************************************************************
 *
 * Find the key data associated with the specified key and return a
 * pointer to the record.  If no such key data exists, then exit with
 * an error message.
 *
 *************************************************************************
 */

MemoryDatabase::KeyData *
MemoryDatabase::findKeyDataOrExit(
   const std::string& key)
{
   for (std::list<KeyData>::iterator i = d_keyvalues.begin();
        i != d_keyvalues.end(); i++) {
      if (key == i->d_key) {
         return &(*i);
      }
   }
   MEMORY_DB_ERROR("Key ``" << key << "'' does not exist in the database...");
   return NULL;
}

/*
 *************************************************************************
 *
 * Print the entire database to the specified output stream.
 *
 *************************************************************************
 */

void
MemoryDatabase::printClassData(
   std::ostream& os)
{
   printDatabase(os, 0, PRINT_DEFAULT | PRINT_INPUT | PRINT_UNUSED);
}

/*
 *************************************************************************
 *
 * Print database data to the specified output stream.
 *
 *************************************************************************
 */

void
MemoryDatabase::printDatabase(
   std::ostream& os,
   const int indent,
   const int toprint) const
{
   /*
    * Get the maximum key width in the output (excluding databases)
    */

   long width = 0;
   for (std::list<KeyData>::const_iterator k = d_keyvalues.begin();
        k != d_keyvalues.end(); k++) {
      if (((k->d_from_default) && (toprint & PRINT_DEFAULT))
          || ((k->d_accessed) && (toprint & PRINT_INPUT))
          || (!(k->d_accessed) && (toprint & PRINT_UNUSED))) {
         if (k->d_type != Database::SAMRAI_DATABASE) {
            const long keywidth = k->d_key.length();
            if (keywidth > width) {
               width = keywidth;
            }
         }
      }
   }

   /*
    * Iterate over all non-database keys in the database and output key values
    */

   indentStream(os, indent);
   os << d_database_name << " {\n";
   for (std::list<KeyData>::const_iterator i = d_keyvalues.begin();
        i != d_keyvalues.end(); i++) {

      if (((i->d_from_default) && (toprint & PRINT_DEFAULT))
          || ((i->d_accessed) && (toprint & PRINT_INPUT))
          || (!(i->d_accessed) && (toprint & PRINT_UNUSED))) {

#ifndef LACKS_SSTREAM
         std::ostringstream sstream;
#else
         char sstream_buffer[SSTREAM_BUFFER];
         std::ostrstream sstream(sstream_buffer, SSTREAM_BUFFER);
#endif

         switch (i->d_type) {

            case Database::SAMRAI_INVALID: {
               break;
            }

            case Database::SAMRAI_DATABASE: {
               break;
            }

            case Database::SAMRAI_BOOL: {
               indentStream(sstream, indent + 3);
               sstream << i->d_key;
               indentStream(sstream, width - i->d_key.length());
               sstream << " = ";
               const int n = i->d_boolean.getSize();
               for (int j = 0; j < n; j++) {
                  sstream << (i->d_boolean[j] ? "TRUE" : "FALSE");
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_BOX: {
               indentStream(sstream, indent + 3);
               sstream << i->d_key;
               indentStream(sstream, width - i->d_key.length());
               sstream << " = ";
               const int n = i->d_box.getSize();
               for (int j = 0; j < n; j++) {
                  const int m = i->d_box[j].getDim().getValue();
                  sstream << "[(";
                  for (int k = 0; k < m; k++) {
                     sstream << i->d_box[j].lower(k);
                     if (k < m - 1) {
                        sstream << ",";
                     }
                  }
                  sstream << "),(";
                  for (int l = 0; l < m; l++) {
                     sstream << i->d_box[j].upper(l);
                     if (l < m - 1) {
                        sstream << ",";
                     }
                  }
                  sstream << ")]";
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_CHAR: {
               indentStream(sstream, indent + 3);
               sstream << i->d_key;
               indentStream(sstream, width - i->d_key.length());
               sstream << " = ";
               const int n = i->d_char.getSize();
               for (int j = 0; j < n; j++) {
                  sstream << "'" << i->d_char[j] << "'";
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_COMPLEX: {
               indentStream(sstream, indent + 3);
               sstream << i->d_key;
               indentStream(sstream, width - i->d_key.length());
               sstream << " = ";
               const int n = i->d_complex.getSize();
               for (int j = 0; j < n; j++) {
                  sstream << i->d_complex[j];
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_DOUBLE: {
               indentStream(sstream, indent + 3);
               sstream << i->d_key;
               indentStream(sstream, width - i->d_key.length());
               sstream << " = ";
               const int n = i->d_double.getSize();
               for (int j = 0; j < n; j++) {
                  sstream << i->d_double[j];
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_FLOAT: {
               indentStream(sstream, indent + 3);
               sstream << i->d_key;
               indentStream(sstream, width - i->d_key.length());
               sstream << " = ";
               const int n = i->d_float.getSize();
               for (int j = 0; j < n; j++) {
                  sstream << i->d_float[j];
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_INT: {
               indentStream(sstream, indent + 3);
               sstream << i->d_key;
               indentStream(sstream, width - i->d_key.length());
               sstream << " = ";
               const int n = i->d_integer.getSize();
               for (int j = 0; j < n; j++) {
                  sstream << i->d_integer[j];
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }

            case Database::SAMRAI_STRING: {
               indentStream(sstream, indent + 3);
               sstream << i->d_key;
               indentStream(sstream, width - i->d_key.length());
               sstream << " = ";
               const int n = i->d_string.getSize();
               for (int j = 0; j < n; j++) {
                  sstream << "\"" << i->d_string[j] << "\"";
                  if (j < n - 1) {
                     sstream << ", ";
                  }
               }
               break;
            }
         }

         /*
          * Output whether the key was used or default in column 60
          */

         if (i->d_type != Database::SAMRAI_DATABASE) {
#ifndef LACKS_SSTREAM
            const long tab = 59 - sstream.str().length();
#else
            const long tab = 59 - sstream.pcount();
#endif
            if (tab > 0) {
               indentStream(sstream, tab);
            }
            if (i->d_from_default) {
               sstream << " // from default";
            } else if (i->d_accessed) {
               sstream << " // input used";
            } else {
               sstream << " // input not used";
            }

//            sstream << std::endl << ends;
            sstream << std::endl;
            os << sstream.str();
         }
      }
   }

   /*
    * Finally, output all databases in the current key list
    */

   for (std::list<KeyData>::const_iterator j = d_keyvalues.begin();
        j != d_keyvalues.end(); j++) {
      if (j->d_type == Database::SAMRAI_DATABASE) {
         boost::shared_ptr<MemoryDatabase> db(
            j->d_database,
            boost::detail::dynamic_cast_tag());
         db->printDatabase(os, indent + 3, toprint);
      }
   }

   indentStream(os, indent);
   os << "}\n";
}

}
}

#if !defined(__BGL_FAMILY__) && defined(__xlC__)
/*
 * Unsuppress XLC warnings
 */
#pragma report(enable, CPPC5334)
#pragma report(enable, CPPC5328)
#endif
