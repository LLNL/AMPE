/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   A memory database structure that stores (key,value) pairs in memory
 *
 ************************************************************************/

#ifndef included_tbox_MemoryDatabase
#define included_tbox_MemoryDatabase

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Database.h"

#include <list>

namespace SAMRAI {
namespace tbox {

/**
 * Class MemoryDatabase stores (key,value) pairs in a hierarchical
 * database.  Each value may be another database, boolean, box, character,
 * complex, double, float, integer, or string.  Note that boxes are stored
 * using the toolbox box class that can store boxes of any dimension in the
 * same data structure.
 *
 * See the Database class documentation for a description of the
 * generic database interface.
 *
 * The memory database is used for storing input (a "deck") in memory; the
 * derived class InputDatabase use for this purpose; InputDatabase is
 * kept around for historical reasons and the name is a bit more clear
 * than "memory" database but this class is very generic.  The Parser
 * is really the class that process the input; it parses the input file
 * into a database.
 *
 * Note that the memory database will attempt to promote numerical types
 * where appropriate.  The promotion chain is int -> float -> double ->
 * complex.  For example, an integer key will be promoted to a complex
 * value if isComplex() or getComplex() is called.  Double values will also
 * be truncated to floats (with loss of information) if a float call is
 * made on a double value.
 *
 * It is assumed that all processors will access the database in the same
 * manner.  Thus, all error messages are output to pout instead of perr.
 */

class MemoryDatabase:public Database
{
public:
   /**
    * The memory database constructor creates an empty database with the
    * specified name.
    */
   explicit MemoryDatabase(
      const std::string& name);

   /**
    * The memory database destructor deallocates the data in the database.
    */
   virtual ~MemoryDatabase();

   /**
    * Create a new database file.
    *
    * Returns true if successful.
    *
    * @param name name of database. Normally a filename.
    */
   virtual bool
   create(
      const std::string& name);

   /**
    * Open an existing database file.
    *
    * Returns true if successful.
    *
    * @param name name of database. Normally a filename.
    *
    * @param read_write_mode Open the database in read-write
    * mode instead of read-only mode.  NOTE: This class currently
    * does not support read-only mode, so this flag must be true.
    */
   virtual bool
   open(
      const std::string& name,
      const bool read_write_mode = false);

   /**
    * Close the database.
    *
    * Returns true if successful.
    *
    * If the database is currently open then close it.  This should
    * flush all data to the file (if the database is on disk).
    */
   virtual bool
   close();

   /**
    * Return string name of memory database object.
    */
   virtual std::string
   getName() const;

   /**
    * Return true if the specified key exists in the database and false
    * otherwise.
    */
   virtual bool
   keyExists(
      const std::string& key);

   /**
    * Return all keys in the database.
    */
   virtual Array<std::string>
   getAllKeys();

   /**
    * @brief Return the type of data associated with the key.
    *
    * If the key does not exist, then INVALID is returned
    *
    * @param key Key name in database.
    */
   virtual enum DataType
   getArrayType(
      const std::string& key);

   /**
    * Return the size of the array associated with the key.  If the key
    * does not exist, then zero is returned.
    */
   virtual int
   getArraySize(
      const std::string& key);

   /**
    * Return whether the specified key represents a database entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool
   isDatabase(
      const std::string& key);

   /**
    * Create a new database with the specified key name.  If the key already
    * exists in the database, then the old key record is deleted and the new
    * one is silently created in its place.
    */
   virtual boost::shared_ptr<Database>
   putDatabase(
      const std::string& key);

   /**
    * Get the database with the specified key name.  If the specified
    * key does not exist in the database or it is not a database, then
    * an error message is printed and the program exits.
    */
   virtual boost::shared_ptr<Database>
   getDatabase(
      const std::string& key);

   /**
    * Return whether the specified key represents a boolean entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool
   isBool(
      const std::string& key);

   /**
    * Create a boolean scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putBool(
      const std::string& key,
      const bool& data);

   /**
    * Create a boolean array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putBoolArray(
      const std::string& key,
      const Array<bool>& data);

   /**
    * Create a boolean array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putBoolArray(
      const std::string& key,
      const bool * const data,
      const int nelements);

   /**
    * Get a boolean entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * boolean scalar, then an error message is printed and the program
    * exits.
    */
   virtual bool
   getBool(
      const std::string& key);

   /**
    * Get a boolean entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a boolean scalar,
    * then an error message is printed and the program exits.
    */
   virtual bool
   getBoolWithDefault(
      const std::string& key,
      const bool& defaultvalue);

   /**
    * Get a boolean entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a boolean array, then an error message is printed and
    * the program exits.
    */
   virtual Array<bool>
   getBoolArray(
      const std::string& key);

   /**
    * Get a boolean entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a boolean array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    */
   virtual void
   getBoolArray(
      const std::string& key,
      bool * data,
      const int nelements);

   /**
    * Return whether the specified key represents a box entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool
   isDatabaseBox(
      const std::string& key);

   /**
    * Create a box scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putDatabaseBox(
      const std::string& key,
      const DatabaseBox& data);

   /**
    * Create a box array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putDatabaseBoxArray(
      const std::string& key,
      const Array<DatabaseBox>& data);

   /**
    * Create a box array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putDatabaseBoxArray(
      const std::string& key,
      const DatabaseBox * const data,
      const int nelements);

   /**
    * Get a box entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * box scalar, then an error message is printed and the program
    * exits.
    */
   virtual DatabaseBox
   getDatabaseBox(
      const std::string& key);

   /**
    * Get a box entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a box scalar,
    * then an error message is printed and the program exits.
    */
   virtual DatabaseBox
   getDatabaseBoxWithDefault(
      const std::string& key,
      const DatabaseBox& defaultvalue);

   /**
    * Get a box entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a box array, then an error message is printed and
    * the program exits.
    */
   virtual Array<DatabaseBox>
   getDatabaseBoxArray(
      const std::string& key);

   /**
    * Get a box entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a box array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    */
   virtual void
   getDatabaseBoxArray(
      const std::string& key,
      DatabaseBox * data,
      const int nelements);

   /**
    * Return whether the specified key represents a character entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool
   isChar(
      const std::string& key);

   /**
    * Create a character scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putChar(
      const std::string& key,
      const char& data);

   /**
    * Create a character array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putCharArray(
      const std::string& key,
      const Array<char>& data);

   /**
    * Create a character array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putCharArray(
      const std::string& key,
      const char * const data,
      const int nelements);

   /**
    * Get a character entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * character scalar, then an error message is printed and the program
    * exits.
    */
   virtual char
   getChar(
      const std::string& key);

   /**
    * Get a character entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a character scalar,
    * then an error message is printed and the program exits.
    */
   virtual char
   getCharWithDefault(
      const std::string& key,
      const char& defaultvalue);

   /**
    * Get a character entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a character array, then an error message is printed and
    * the program exits.
    */
   virtual Array<char>
   getCharArray(
      const std::string& key);

   /**
    * Get a character entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a character array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    */
   virtual void
   getCharArray(
      const std::string& key,
      char * data,
      const int nelements);

   /**
    * Return whether the specified key represents a complex entry.  If
    * the key does not exist, then false is returned.  Complex values
    * may be promoted from integers, floats, or doubles.
    */
   virtual bool
   isComplex(
      const std::string& key);

   /**
    * Create a complex scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putComplex(
      const std::string& key,
      const dcomplex& data);

   /**
    * Create a complex array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putComplexArray(
      const std::string& key,
      const Array<dcomplex>& data);

   /**
    * Create a complex array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putComplexArray(
      const std::string& key,
      const dcomplex * const data,
      const int nelements);

   /**
    * Get a complex entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * complex scalar, then an error message is printed and the program
    * exits.  Complex values may be promoted from integers, floats, or
    * doubles.
    */
   virtual dcomplex
   getComplex(
      const std::string& key);

   /**
    * Get a complex entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a complex scalar,
    * then an error message is printed and the program exits.  Complex
    * values may be promoted from integers, floats, or doubles.
    */
   virtual dcomplex
   getComplexWithDefault(
      const std::string& key,
      const dcomplex& defaultvalue);

   /**
    * Get a complex entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a complex array, then an error message is printed and
    * the program exits.  Complex values may be promoted from integers,
    * floats, or doubles.
    */
   virtual Array<dcomplex>
   getComplexArray(
      const std::string& key);

   /**
    * Get a complex entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a complex array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    * Complex values may be promoted from integers, floats, or doubles.
    */
   virtual void
   getComplexArray(
      const std::string& key,
      dcomplex * data,
      const int nelements);

   /**
    * Return whether the specified key represents a double entry.  If
    * the key does not exist, then false is returned.  Double values
    * may be promoted from integers or floats.
    */
   virtual bool
   isDouble(
      const std::string& key);

   /**
    * Create a double scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putDouble(
      const std::string& key,
      const double& data);

   /**
    * Create a double array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putDoubleArray(
      const std::string& key,
      const Array<double>& data);

   /**
    * Create a double array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putDoubleArray(
      const std::string& key,
      const double * const data,
      const int nelements);

   /**
    * Get a double entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * double scalar, then an error message is printed and the program
    * exits.  Double values may be promoted from integers or floats.
    */
   virtual double
   getDouble(
      const std::string& key);

   /**
    * Get a double entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a double scalar, then
    * an error message is printed and the program exits.  Double values may
    * be promoted from integers or floats.
    */
   virtual double
   getDoubleWithDefault(
      const std::string& key,
      const double& defaultvalue);

   /**
    * Get a double entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double array, then an error message is printed and
    * the program exits.  Double values may be promoted from integers
    * or floats.
    */
   virtual Array<double>
   getDoubleArray(
      const std::string& key);

   /**
    * Get a double entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    * Double values may be promoted from integers or floats.
    */
   virtual void
   getDoubleArray(
      const std::string& key,
      double * data,
      const int nelements);

   /**
    * Return whether the specified key represents a float entry.  If
    * the key does not exist, then false is returned.  Float values
    * may be promoted from integers or silently truncated from doubles.
    */
   virtual bool
   isFloat(
      const std::string& key);

   /**
    * Create a float scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putFloat(
      const std::string& key,
      const float& data);

   /**
    * Create a float array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putFloatArray(
      const std::string& key,
      const Array<float>& data);

   /**
    * Create a float array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putFloatArray(
      const std::string& key,
      const float * const data,
      const int nelements);

   /**
    * Get a float entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * float scalar, then an error message is printed and the program
    * exits.  Float values may be promoted from integers or silently
    * truncated from doubles.
    */
   virtual float
   getFloat(
      const std::string& key);

   /**
    * Get a float entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a float scalar, then
    * an error message is printed and the program exits.  Float values may
    * be promoted from integers or silently truncated from doubles.
    */
   virtual float
   getFloatWithDefault(
      const std::string& key,
      const float& defaultvalue);

   /**
    * Get a float entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float array, then an error message is printed and
    * the program exits.  Float values may be promoted from integers
    * or silently truncated from doubles.
    */
   virtual Array<float>
   getFloatArray(
      const std::string& key);

   /**
    * Get a float entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    * Float values may be promoted from integers or silently truncated
    * from doubles.
    */
   virtual void
   getFloatArray(
      const std::string& key,
      float * data,
      const int nelements);

   /**
    * Return whether the specified key represents an integer entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool
   isInteger(
      const std::string& key);

   /**
    * Create an integer scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putInteger(
      const std::string& key,
      const int& data);

   /**
    * Create an integer array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putIntegerArray(
      const std::string& key,
      const Array<int>& data);

   /**
    * Create an integer array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putIntegerArray(
      const std::string& key,
      const int * const data,
      const int nelements);

   /**
    * Get an integer entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * integer scalar, then an error message is printed and the program
    * exits.
    */
   virtual int
   getInteger(
      const std::string& key);

   /**
    * Get an integer entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not an integer scalar,
    * then an error message is printed and the program exits.
    */
   virtual int
   getIntegerWithDefault(
      const std::string& key,
      const int& defaultvalue);

   /**
    * Get an integer entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not an integer array, then an error message is printed and
    * the program exits.
    */
   virtual Array<int>
   getIntegerArray(
      const std::string& key);

   /**
    * Get an integer entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not an integer array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    */
   virtual void
   getIntegerArray(
      const std::string& key,
      int * data,
      const int nelements);

   /**
    * Return whether the specified key represents a std::string entry.  If
    * the key does not exist, then false is returned.
    */
   virtual bool
   isString(
      const std::string& key);

   /**
    * Create a string scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putString(
      const std::string& key,
      const std::string& data);

   /**
    * Create a string array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putStringArray(
      const std::string& key,
      const Array<std::string>& data);

   /**
    * Create a string array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    */
   virtual void
   putStringArray(
      const std::string& key,
      const std::string * const data,
      const int nelements);

   /**
    * Get a string entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * string scalar, then an error message is printed and the program
    * exits.
    */
   virtual std::string
   getString(
      const std::string& key);

   /**
    * Get a string entry in the database with the specified key name.
    * If the specified key does not exist in the database, then the default
    * value is returned.  If the key exists but is not a string scalar,
    * then an error message is printed and the program exits.
    */
   virtual std::string
   getStringWithDefault(
      const std::string& key,
      const std::string& defaultvalue);

   /**
    * Get a string entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a string array, then an error message is printed and
    * the program exits.
    */
   virtual Array<std::string>
   getStringArray(
      const std::string& key);

   /**
    * Get a string entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a string array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    */
   virtual void
   getStringArray(
      const std::string& key,
      std::string * data,
      const int nelements);

   /**
    * @brief Returns the name of this database.
    *
    * The name for the root of the database is the name supplied when creating
    * it.  Names for nested databases are the keyname of the database.
    */
   virtual std::string
   getName();

   /**
    * Return whether the specified key has been accessed by one of the
    * lookup member functions.  If the key does not exist in the database,
    * then false is returned.
    */
   bool
   keyAccessed(
      const std::string& key)
   {
      KeyData* keydata = findKeyData(key);
      return keydata ? keydata->d_accessed : false;
   }

   /**
    * Print the current database to the specified output stream.  After
    * each key, print whether that key came from the a file and was
    * used, came from the file but was not used (unused),
    * or came from a default key value (default).  If no output stream
    * is specified, then data is written to stream pout.
    *
    * NOTE:  under the g++ compiler libraries, printClassData has a
    * maximum output of 4096 characters per line.
    */
   virtual void
   printClassData(
      std::ostream& os = pout);

   /**
    * Print the database keys that were not used to the specified output
    * stream.
    */
   void
   printUnusedKeys(
      std::ostream& os = pout) const
   {
      printDatabase(os, 0, PRINT_UNUSED);
   }

   /**
    * Print the database keys that were set via default calls to the specified
    * output stream.
    */
   void
   printDefaultKeys(
      std::ostream& os = pout) const
   {
      printDatabase(os, 0, PRINT_DEFAULT);
   }

private:
   MemoryDatabase(
      const MemoryDatabase&);                   // not implemented
   void
   operator = (
      const MemoryDatabase&);                   // not implemented

   /*
    * The following structure holds the list of (key,value) pairs stored
    * in the database.  Note that only one of the arrays contains valid
    * data for any particular key.
    */
   struct KeyData {
      std::string d_key;                                // key name
      enum Database::DataType d_type;                   // type of entry
      int d_array_size;                                 // size of array data
      bool d_accessed;                                  // whether accessed
      bool d_from_default;                              // from default key
      boost::shared_ptr<Database> d_database;           // sub-database
      Array<bool> d_boolean;                    // boolean array value
      Array<DatabaseBox> d_box;                 // box array value
      Array<char> d_char;                       // char array value
      Array<dcomplex> d_complex;                // complex array value
      Array<double> d_double;                   // double array value
      Array<float> d_float;                     // float array value
      Array<int> d_integer;                     // integer array value
      Array<std::string> d_string;              // string array value
   };

   /*
    * Private utility routines for managing the database
    */
   bool
   deleteKeyIfFound(
      const std::string& key);
   KeyData *
   findKeyData(
      const std::string& key);
   KeyData *
   findKeyDataOrExit(
      const std::string& key);
   static void
   indentStream(
      std::ostream& os,
      const long indent)
   {
      for (int i = 0; i < indent; i++) {
         os << " ";
      }
   }
   void
   printDatabase(
      std::ostream& os,
      const int indent,
      const int toprint) const;

   /*
    * Private data members - name and a list of (key,value) pairs
    */
   std::string d_database_name;
   std::list<KeyData> d_keyvalues;

   static const int PRINT_DEFAULT;
   static const int PRINT_INPUT;
   static const int PRINT_UNUSED;
   static const int SSTREAM_BUFFER;
};

}
}

#endif
