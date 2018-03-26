/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   An input database structure that stores (key,value) pairs
 *
 ************************************************************************/

#include "SAMRAI/tbox/NullDatabase.h"

namespace SAMRAI {
namespace tbox {

NullDatabase::NullDatabase()
{
}

/*
 *************************************************************************
 *
 * The virtual destructor deallocates database data.
 *
 *************************************************************************
 */

NullDatabase::~NullDatabase()
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
NullDatabase::create(
   const std::string& name)
{
   NULL_USE(name);
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
NullDatabase::open(
   const std::string& name,
   const bool read_write_mode)
{
   NULL_USE(name);
   NULL_USE(read_write_mode);
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
NullDatabase::close()
{
   return true;
}

/*
 *************************************************************************
 *
 * Always returns true.
 *
 *************************************************************************
 */

bool
NullDatabase::keyExists(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

/*
 *************************************************************************
 *
 * Return an empty Array<std::string>.
 *
 *************************************************************************
 */

Array<std::string>
NullDatabase::getAllKeys()
{
   Array<std::string> keys(0);
   return keys;
}

/*
 *************************************************************************
 *
 * Always returns INVALID.
 *
 *************************************************************************
 */

Database::DataType
NullDatabase::getArrayType(
   const std::string& key)
{
   NULL_USE(key);
   return Database::SAMRAI_INVALID;
}

/*
 *************************************************************************
 *
 * Always returns 0.
 *
 *************************************************************************
 */

int
NullDatabase::getArraySize(
   const std::string& key)
{
   NULL_USE(key);
   return 0;
}

/*
 *************************************************************************
 *
 * Member functions that manage the database values within the database.
 *
 *************************************************************************
 */

bool
NullDatabase::isDatabase(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

boost::shared_ptr<Database>
NullDatabase::putDatabase(
   const std::string& key)
{
   NULL_USE(key);
   return boost::shared_ptr<Database>(this);
}

boost::shared_ptr<Database>
NullDatabase::getDatabase(
   const std::string& key)
{
   NULL_USE(key);
   return boost::make_shared<NullDatabase>();
}

/*
 *************************************************************************
 *
 * Member functions that manage boolean values within the database.
 *
 *************************************************************************
 */

bool
NullDatabase::isBool(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

void
NullDatabase::putBoolArray(
   const std::string& key,
   const bool * const data,
   const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

Array<bool>
NullDatabase::getBoolArray(
   const std::string& key)
{
   NULL_USE(key);
   Array<bool> empty(0);
   return empty;
}

/*
 *************************************************************************
 *
 * Member functions that manage box values within the database.
 *
 *************************************************************************
 */

bool
NullDatabase::isDatabaseBox(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

void
NullDatabase::putDatabaseBoxArray(
   const std::string& key,
   const DatabaseBox * const data,
   const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

Array<DatabaseBox>
NullDatabase::getDatabaseBoxArray(
   const std::string& key)
{
   NULL_USE(key);

   Array<DatabaseBox> empty(0);
   return empty;
}

/*
 *************************************************************************
 *
 * Member functions that manage character values within the database.
 *
 *************************************************************************
 */

bool
NullDatabase::isChar(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

void
NullDatabase::putCharArray(
   const std::string& key,
   const char * const data,
   const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

Array<char>
NullDatabase::getCharArray(
   const std::string& key)
{
   NULL_USE(key);

   Array<char> empty(0);
   return empty;
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
NullDatabase::isComplex(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

void
NullDatabase::putComplexArray(
   const std::string& key,
   const dcomplex * const data,
   const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

Array<dcomplex>
NullDatabase::getComplexArray(
   const std::string& key)
{
   NULL_USE(key);

   Array<dcomplex> empty(0);
   return empty;
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
NullDatabase::isDouble(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

void
NullDatabase::putDoubleArray(
   const std::string& key,
   const double * const data,
   const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

Array<double>
NullDatabase::getDoubleArray(
   const std::string& key)
{
   NULL_USE(key);
   Array<double> empty(0);
   return empty;
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
NullDatabase::isFloat(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

void
NullDatabase::putFloatArray(
   const std::string& key,
   const float * const data,
   const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

Array<float>
NullDatabase::getFloatArray(
   const std::string& key)
{
   NULL_USE(key);

   Array<float> empty(0);
   return empty;
}

/*
 *************************************************************************
 *
 * Member functions that manage integer values within the database.
 *
 *************************************************************************
 */

bool
NullDatabase::isInteger(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

void
NullDatabase::putIntegerArray(
   const std::string& key,
   const int * const data,
   const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

Array<int>
NullDatabase::getIntegerArray(
   const std::string& key)
{
   NULL_USE(key);

   Array<int> empty(0);
   return empty;
}

/*
 *************************************************************************
 *
 * Member functions that manage string values within the database.
 *
 *************************************************************************
 */

bool
NullDatabase::isString(
   const std::string& key)
{
   NULL_USE(key);
   return true;
}

void
NullDatabase::putStringArray(
   const std::string& key,
   const std::string * const data,
   const int nelements)
{
   NULL_USE(key);
   NULL_USE(data);
   NULL_USE(nelements);
}

Array<std::string>
NullDatabase::getStringArray(
   const std::string& key)
{
   NULL_USE(key);
   Array<std::string> empty(0);
   return empty;
}

std::string
NullDatabase::getName(
   void)
{
   return std::string();
}

/*
 *************************************************************************
 *
 * Does nothing.
 *
 *************************************************************************
 */

void
NullDatabase::printClassData(
   std::ostream& os)
{
   NULL_USE(os);
}

}
}
