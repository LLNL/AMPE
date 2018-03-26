/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Some simple generic database test functions
 *
 ************************************************************************/

#include "database_tests.h"
#include "database_values.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <vector>

/*
 * Write database and test contents.
 */
void setupTestData(
   void) {
   arraydb_boxArray0.setDim(tbox::Dimension(3));
   arraydb_boxArray1.setDim(tbox::Dimension(2));
   arraydb_boxArray2.setDim(tbox::Dimension(1));
}

/*
 * Write database and test contents.
 */
void writeTestData(
   boost::shared_ptr<tbox::Database> db)
{
   if (!db) {
      tbox::perr << "FAILED: - Test #0-write: database to write to is null"
                 << endl;
      tbox::SAMRAI_MPI::abort();
   }

   /*
    * Build database hierarchy and test.
    */

   boost::shared_ptr<tbox::Database> arraydb(db->putDatabase("Array Entries"));

   boost::shared_ptr<tbox::Database> scalardb(
      db->putDatabase("Scalar Entries"));
   boost::shared_ptr<tbox::Database> scalardb_empty(
      scalardb->putDatabase("Empty"));
   boost::shared_ptr<tbox::Database> scalardb_full(
      scalardb->putDatabase("Full"));
   boost::shared_ptr<tbox::Database> defaultdb(db->putDatabase("Default"));
   boost::shared_ptr<tbox::Database> namesdb(db->putDatabase("Name Entries"));
   boost::shared_ptr<tbox::Database> vectordb(db->putDatabase("stl_vector"));

   NULL_USE(defaultdb);

   if (!arraydb) {
      tbox::perr << "FAILED: - Test #1a-write: `arraydb' is null" << endl;
      tbox::SAMRAI_MPI::abort();
   }
   if (!scalardb) {
      tbox::perr << "FAILED: - Test #1b-write: `scalardb' is null" << endl;
      tbox::SAMRAI_MPI::abort();
   }
   if (!scalardb_empty) {
      tbox::perr << "FAILED: - Test #1c-write: `scalardb_empty' is null"
                 << endl;
      tbox::SAMRAI_MPI::abort();
   }
   if (!scalardb_full) {
      tbox::perr << "FAILED: - Test #1d-write: `scalardb_full' is null"
                 << endl;
      tbox::SAMRAI_MPI::abort();
   }
   if (!vectordb) {
      tbox::perr << "FAILED: - Test #1e-write: `vectordb' is null"
                 << endl;
      tbox::SAMRAI_MPI::abort();
   }

   /*
    * Set array values and write to database hierarchy.
    */

   tbox::Array<dcomplex> arraydb_dcomplexArray(3);
   arraydb_dcomplexArray[0] = arraydb_dcomplexArray0;
   arraydb_dcomplexArray[1] = arraydb_dcomplexArray1;
   arraydb_dcomplexArray[2] = arraydb_dcomplexArray2;

   tbox::Array<bool> arraydb_boolArray(3);
   arraydb_boolArray[0] = arraydb_boolArray0;
   arraydb_boolArray[1] = arraydb_boolArray1;
   arraydb_boolArray[2] = arraydb_boolArray2;

   tbox::Array<int> arraydb_intArray(5);
   arraydb_intArray[0] = arraydb_intArray0;
   arraydb_intArray[1] = arraydb_intArray1;
   arraydb_intArray[2] = arraydb_intArray2;
   arraydb_intArray[3] = arraydb_intArray3;
   arraydb_intArray[4] = arraydb_intArray4;

   tbox::Array<string> arraydb_stringArray(3);
   arraydb_stringArray[0] = arraydb_stringArray0;
   arraydb_stringArray[1] = arraydb_stringArray1;
   arraydb_stringArray[2] = arraydb_stringArray2;

   tbox::Array<float> arraydb_floatArray(5);
   arraydb_floatArray[0] = arraydb_floatArray0;
   arraydb_floatArray[1] = arraydb_floatArray1;
   arraydb_floatArray[2] = arraydb_floatArray2;
   arraydb_floatArray[3] = arraydb_floatArray3;
   arraydb_floatArray[4] = arraydb_floatArray4;

   tbox::Array<double> arraydb_doubleArray(6);
   arraydb_doubleArray[0] = arraydb_doubleArray0;
   arraydb_doubleArray[1] = arraydb_doubleArray1;
   arraydb_doubleArray[2] = arraydb_doubleArray2;
   arraydb_doubleArray[3] = arraydb_doubleArray3;
   arraydb_doubleArray[4] = arraydb_doubleArray4;
   arraydb_doubleArray[5] = arraydb_doubleArray5;

   tbox::Array<char> arraydb_charArray(2);
   arraydb_charArray[0] = arraydb_charArray0;
   arraydb_charArray[1] = arraydb_charArray1;

   tbox::Array<tbox::DatabaseBox> arraydb_boxArray(3);
   arraydb_boxArray[0] = arraydb_boxArray0;
   arraydb_boxArray[1] = arraydb_boxArray1;
   arraydb_boxArray[2] = arraydb_boxArray2;

   db->putFloat("float_val", db_float_val);
   db->putInteger("int_val", db_int_val);

   arraydb->putComplexArray("ComplexArray", arraydb_dcomplexArray);
   arraydb->putDatabaseBoxArray("BoxArray", arraydb_boxArray);
   arraydb->putBoolArray("BoolArray", arraydb_boolArray);
   arraydb->putIntegerArray("IntArray", arraydb_intArray);
   arraydb->putStringArray("StringArray", arraydb_stringArray);
   arraydb->putFloatArray("FloatArray", arraydb_floatArray);
   arraydb->putDoubleArray("DoubleArray", arraydb_doubleArray);
   arraydb->putCharArray("CharArray", arraydb_charArray);

   scalardb->putFloat("float1", scalardb_float1);
   scalardb->putFloat("float2", scalardb_float2);
   scalardb->putFloat("float3", scalardb_float3);

   scalardb_full->putDouble("thisDouble", scalardb_full_thisDouble);
   scalardb_full->putComplex("thisComplex", scalardb_full_thisComplex);
   scalardb_full->putInteger("thisInt", scalardb_full_thisInt);
   scalardb_full->putFloat("thisFloat", scalardb_full_thisFloat);
   scalardb_full->putBool("thisBool", scalardb_full_thisBool);
   scalardb_full->putString("thisString", scalardb_full_thisString);
   scalardb_full->putChar("thisChar", scalardb_full_thisChar);
   scalardb_full->putDatabaseBox("thisBox", scalardb_full_thisBox);

   namesdb->putDouble("Name with spaces", scalardb_full_thisDouble);
   namesdb->putDouble("Name-with-dashes", scalardb_full_thisDouble);
   namesdb->putDouble("Name-with-!@#$%^&*()_+-=", scalardb_full_thisDouble);

   std::vector<hier::IntVector> vector_IntVector(2);
   vector_IntVector[0] = intVector1;
   vector_IntVector[1] = intVector2;

   vectordb->putVector("vector_IntVector", vector_IntVector);

   testDatabaseContents(db, "write");
}

/*
 * Read database and test contents.
 */
void readTestData(
   boost::shared_ptr<tbox::Database> db)
{
   testDatabaseContents(db, "read");
}

/*
 * Test contents of database.
 */
void testDatabaseContents(
   boost::shared_ptr<tbox::Database> db,
   const string& tag)
{

   if (!db) {
      tbox::perr << "FAILED: - Test #0-" << tag
                 << ": database to read from is null" << endl;
      ++number_of_failures;
   }

   boost::shared_ptr<tbox::Database> arraydb(db->getDatabase("Array Entries"));

   boost::shared_ptr<tbox::Database> scalardb(
      db->getDatabase("Scalar Entries"));
   boost::shared_ptr<tbox::Database> scalardb_empty(
      scalardb->getDatabase("Empty"));
   boost::shared_ptr<tbox::Database> scalardb_full(
      scalardb->getDatabase("Full"));
   boost::shared_ptr<tbox::Database> defaultdb(db->getDatabase("Default"));

   boost::shared_ptr<tbox::Database> namesdb(db->getDatabase("Name Entries"));

   boost::shared_ptr<tbox::Database> vectordb(db->getDatabase("stl_vector"));

   if (!arraydb) {
      tbox::perr << "FAILED: - Test #1a-" << tag
                 << ": `arraydb' is null" << endl;
      ++number_of_failures;
   }
   if (!scalardb) {
      tbox::perr << "FAILED: - Test #1b-" << tag
                 << ": `scalardb' is null" << endl;
      ++number_of_failures;
   }
   if (!scalardb_empty) {
      tbox::perr << "FAILED: - Test #1c-" << tag
                 << ": `scalardb_empty' is null" << endl;
      ++number_of_failures;
   }
   if (!scalardb_full) {
      tbox::perr << "FAILED: - Test #1d-" << tag
                 << ": `scalardb_full' is null" << endl;
      ++number_of_failures;
   }

   if (!vectordb) {
      tbox::perr << "FAILED: - Test #1e-" << tag
                 << ": `vectordb' is null" << endl;
      ++number_of_failures;
   }

   tbox::Array<string> dbkeys = db->getAllKeys();
   tbox::Array<string> arraydbkeys = arraydb->getAllKeys();
   tbox::Array<string> scalardbkeys = scalardb->getAllKeys();
   tbox::Array<string> scalardb_emptykeys = scalardb_empty->getAllKeys();
   tbox::Array<string> scalardb_fullkeys = scalardb_full->getAllKeys();

   int i, nkeys;

   if (dbkeys.getSize() != 7) {
      tbox::perr << "FAILED: - Test #2a-" << tag
                 << ": # `db' keys wrong" << endl;
      ++number_of_failures;
   }
   nkeys = arraydbkeys.getSize();
   if (arraydbkeys.getSize() != 8) {
      tbox::perr << "FAILED: - Test #2b-" << tag
                 << ": # `arraydb' keys wrong"
                 << "\n\tFound " << nkeys << " keys:"
      ;
      ++number_of_failures;
      for (i = 0; i < nkeys; ++i) {
         tbox::pout << "\n\t\t" << i << ": '" << arraydbkeys[i] << "'";
      }
      tbox::pout << endl;
   }
   if (scalardbkeys.getSize() != 5) {
      tbox::perr << "FAILED: - Test #2c-" << tag
                 << ": # `scalardb' keys wrong" << endl;
      ++number_of_failures;
   }
   if (scalardb_emptykeys.getSize() != 0) {
      tbox::perr << "FAILED: - Test #2d-" << tag
                 << ": # `scalardb_empty' keys wrong" << endl;
      ++number_of_failures;
   }
   if (scalardb_fullkeys.getSize() != 8) {
      tbox::perr << "FAILED: - Test #2e-" << tag
                 << ": `scalardb_full' size is wrong" << endl
                 << " returned : " << scalardb_fullkeys.getSize() << endl
                 << " expected : " << 8 << endl;
      ++number_of_failures;
   }

   if (!db->isDatabase("Array Entries")) {
      tbox::perr << "FAILED: - #3a-" << tag
                 << ": `Array Entries' not a database" << endl;
      ++number_of_failures;
   }
   if (!db->isDatabase("Scalar Entries")) {
      tbox::perr << "FAILED: - #3b-" << tag
                 << ": `Scalar Entries' not a database" << endl;
      ++number_of_failures;
   }
   float tdb_float_val = db->getFloat("float_val");
   if (!tbox::MathUtilities<float>::equalEps(db_float_val, db_float_val)) {
      tbox::perr << "FAILED: - Test #3c-" << tag
                 << ": `RestartTester' database"
                 << "\n   Returned `float_val' = " << tdb_float_val
                 << "  , Expected = " << db_float_val << endl;
      ++number_of_failures;
   }
   int tdb_int_val = db->getInteger("int_val");
   if (tdb_int_val != db_int_val) {
      tbox::perr << "FAILED: - Test #3d-" << tag
                 << ": `RestartTester' database"
                 << "\n   Returned `int_val' = " << tdb_int_val
                 << "  , Expected = " << db_int_val << endl;
      ++number_of_failures;
   }

   /*
    * Set array values to test database.
    */

   tbox::Array<dcomplex> arraydb_dcomplexArray(3);
   arraydb_dcomplexArray[0] = arraydb_dcomplexArray0;
   arraydb_dcomplexArray[1] = arraydb_dcomplexArray1;
   arraydb_dcomplexArray[2] = arraydb_dcomplexArray2;

   tbox::Array<bool> arraydb_boolArray(3);
   arraydb_boolArray[0] = arraydb_boolArray0;
   arraydb_boolArray[1] = arraydb_boolArray1;
   arraydb_boolArray[2] = arraydb_boolArray2;

   tbox::Array<int> arraydb_intArray(5);
   arraydb_intArray[0] = arraydb_intArray0;
   arraydb_intArray[1] = arraydb_intArray1;
   arraydb_intArray[2] = arraydb_intArray2;
   arraydb_intArray[3] = arraydb_intArray3;
   arraydb_intArray[4] = arraydb_intArray4;

   tbox::Array<string> arraydb_stringArray(3);
   arraydb_stringArray[0] = arraydb_stringArray0;
   arraydb_stringArray[1] = arraydb_stringArray1;
   arraydb_stringArray[2] = arraydb_stringArray2;

   tbox::Array<float> arraydb_floatArray(5);
   arraydb_floatArray[0] = arraydb_floatArray0;
   arraydb_floatArray[1] = arraydb_floatArray1;
   arraydb_floatArray[2] = arraydb_floatArray2;
   arraydb_floatArray[3] = arraydb_floatArray3;
   arraydb_floatArray[4] = arraydb_floatArray4;

   tbox::Array<double> arraydb_doubleArray(6);
   arraydb_doubleArray[0] = arraydb_doubleArray0;
   arraydb_doubleArray[1] = arraydb_doubleArray1;
   arraydb_doubleArray[2] = arraydb_doubleArray2;
   arraydb_doubleArray[3] = arraydb_doubleArray3;
   arraydb_doubleArray[4] = arraydb_doubleArray4;
   arraydb_doubleArray[5] = arraydb_doubleArray5;

   tbox::Array<char> arraydb_charArray(2);
   arraydb_charArray[0] = arraydb_charArray0;
   arraydb_charArray[1] = arraydb_charArray1;

   tbox::Array<tbox::DatabaseBox> arraydb_boxArray(3);
   arraydb_boxArray[0] = arraydb_boxArray0;
   arraydb_boxArray[1] = arraydb_boxArray1;
   arraydb_boxArray[2] = arraydb_boxArray2;

   int tsize = 0;

   tbox::Array<dcomplex> tarraydb_dcomplexArray =
      arraydb->getComplexArray("ComplexArray");
   tsize = tarraydb_dcomplexArray.getSize();
   if (tsize != arraydb_dcomplexArray.getSize()) {
      tbox::perr << "FAILED: - Test #4a-" << tag
                 << ": `Array Entries' database"
                 << "\n   Returned `ComplexArray' size = " << tsize
                 << "  , Expected = " << arraydb_dcomplexArray.getSize()
                 << endl;
      ++number_of_failures;
   }
   for (i = 0; i < tsize; i++) {
      if (tarraydb_dcomplexArray[i] != arraydb_dcomplexArray[i]) {
         tbox::perr << "FAILED: - Test #4b-" << tag
                    << ": `Array Entries' database"
                    << "\n   `ComplexArray' entry " << i << " incorrect"
                    << endl;
         ++number_of_failures;
      }
   }
   tbox::Array<bool> tarraydb_boolArray =
      arraydb->getBoolArray("BoolArray");
   tsize = tarraydb_boolArray.getSize();
   if (tsize != arraydb_boolArray.getSize()) {
      tbox::perr << "FAILED: - Test #4c-" << tag
                 << ": `Array Entries' database"
                 << "\n   Returned `BoolArray' size = " << tsize
                 << "  , Expected = " << arraydb_boolArray.getSize() << endl;
      ++number_of_failures;
   }
   for (i = 0; i < tsize; i++) {
      if (tarraydb_boolArray[i] != arraydb_boolArray[i]) {
         tbox::perr << "FAILED: - Test #4d-" << tag
                    << ": `Array Entries' database"
                    << "\n   `BoolArray' entry " << i << " incorrect"
                    << "\n   " << tarraydb_boolArray[i] << " should be "
                    << arraydb_boolArray[i] << endl;
         ++number_of_failures;
      }
   }
   tbox::Array<int> tarraydb_intArray =
      arraydb->getIntegerArray("IntArray");
   tsize = tarraydb_intArray.getSize();
   if (tsize != arraydb_intArray.getSize()) {
      tbox::perr << "FAILED: - Test #4e-" << tag
                 << ": `Array Entries' database"
                 << "\n   Returned `IntArray' size = " << tsize
                 << "  , Expected = " << arraydb_intArray.getSize() << endl;
      ++number_of_failures;
   }
   for (i = 0; i < tsize; i++) {
      if (tarraydb_intArray[i] != arraydb_intArray[i]) {
         tbox::perr << "FAILED: - Test #4f-" << tag
                    << ": `Array Entries' database"
                    << "\n   `IntArray' entry " << i << " incorrect" << endl;
         ++number_of_failures;
      }
   }
   tbox::Array<string> tarraydb_stringArray =
      arraydb->getStringArray("StringArray");
   tsize = tarraydb_stringArray.getSize();
   if (tsize != arraydb_stringArray.getSize()) {
      tbox::perr << "FAILED: - Test #4g-" << tag
                 << ": `Array Entries' database"
                 << "\n   Returned `StringArray' size = " << tsize
                 << "  , Expected = " << arraydb_stringArray[i] << endl;
      ++number_of_failures;
   }
   for (i = 0; i < tsize; i++) {
      if (tarraydb_stringArray[i] != arraydb_stringArray[i]) {
         tbox::perr << "FAILED: - Test #4h-" << tag
                    << ": `Array Entries' database"
                    << "\n   `StringArray' entry " << i << " incorrect" << endl;
         ++number_of_failures;
      }
   }
   tbox::Array<float> tarraydb_floatArray =
      arraydb->getFloatArray("FloatArray");
   tsize = tarraydb_floatArray.getSize();
   if (tsize != arraydb_floatArray.getSize()) {
      tbox::perr << "FAILED: - Test #4i-" << tag
                 << ": `Array Entries' database"
                 << "\n   Returned `FloatArray' size = " << tsize
                 << "  , Expected = " << arraydb_floatArray.getSize() << endl;
      ++number_of_failures;
   }
   for (i = 0; i < tsize; i++) {
      if (!tbox::MathUtilities<float>::equalEps(tarraydb_floatArray[i],
             arraydb_floatArray[i])) {
         tbox::perr << "FAILED: - Test #4j-" << tag
                    << ": `Array Entries' database"
                    << "\n   `FloatArray' entry " << i << " incorrect" << endl;
         ++number_of_failures;
      }
   }

   tbox::Array<double> tarraydb_doubleArray =
      arraydb->getDoubleArray("DoubleArray");
   tsize = tarraydb_doubleArray.getSize();
   if (tsize != arraydb_doubleArray.getSize()) {
      tbox::perr << "FAILED: - Test #4k.b-" << tag
                 << ": `Array Entries' database"
                 << "\n   Returned `DoubleArray' size = " << tsize
                 << "  , Expected = " << arraydb_doubleArray.getSize() << endl;
      ++number_of_failures;
   }
   for (i = 0; i < tsize; i++) {
      if (!tbox::MathUtilities<double>::equalEps(tarraydb_doubleArray[i],
             arraydb_doubleArray[i])) {
         tbox::perr << "FAILED: - Test #4l-" << tag
                    << ": `Array Entries' database"
                    << "\n   `DoubleArray' entry " << i << " incorrect" << endl;
         ++number_of_failures;
      }
   }
   tbox::Array<char> tarraydb_charArray =
      arraydb->getCharArray("CharArray");
   tsize = tarraydb_charArray.getSize();
   if (tsize != arraydb_charArray.getSize()) {
      tbox::perr << "FAILED: - Test #4m-" << tag
                 << ": `Array Entries' database"
                 << "\n   Returned `CharArray' size = " << tsize
                 << "  , Expected = " << arraydb_charArray.getSize() << endl;
      ++number_of_failures;
   }
   for (i = 0; i < tsize; i++) {
      if (tarraydb_charArray[i] != arraydb_charArray[i]) {
         tbox::perr << "FAILED: - Test #4l-" << tag
                    << ": `Array Entries' database"
                    << "\n   `CharArray' entry " << i << " incorrect" << endl;
         ++number_of_failures;
      }
   }
   tbox::Array<tbox::DatabaseBox> tarraydb_boxArray =
      arraydb->getDatabaseBoxArray("BoxArray");
   tsize = tarraydb_boxArray.getSize();
   if (tsize != arraydb_boxArray.getSize()) {
      tbox::perr << "FAILED: - Test #4o-" << tag
                 << ": `Array Entries' database"
                 << "\n   Returned `BoxArray' size = " << tsize
                 << "  , Expected = " << arraydb_boxArray.getSize() << endl;
      ++number_of_failures;
   }
   for (i = 0; i < tsize; i++) {
      if (!(tarraydb_boxArray[i] == arraydb_boxArray[i])) {
         tbox::perr << "FAILED: - Test #4p-" << tag
                    << ": `Array Entries' database"
                    << "\n   `BoxArray' entry " << i << " incorrect" << endl;
         ++number_of_failures;
      }
   }

   if (!scalardb->isDatabase("Empty")) {
      tbox::perr << "FAILED: - #5a-" << tag
                 << ": `Empty' not a database" << endl;
      ++number_of_failures;
   }
   if (!scalardb->isDatabase("Full")) {
      tbox::perr << "FAILED: - #5b-" << tag
                 << ": `Full' not a database" << endl;
      ++number_of_failures;
   }
   float tscalardb_float1 = scalardb->getFloat("float1");
   if (!tbox::MathUtilities<float>::equalEps(tscalardb_float1,
          scalardb_float1)) {
      tbox::perr << "FAILED: - Test #5c-" << tag
                 << ": `Scalar Entries' database"
                 << "\n   Returned `float1' = " << tscalardb_float1
                 << "  , Expected = " << scalardb_float1 << endl;
      ++number_of_failures;
   }
   float tscalardb_float2 = scalardb->getFloat("float2");
   if (!tbox::MathUtilities<float>::equalEps(tscalardb_float2,
          scalardb_float2)) {
      tbox::perr << "FAILED: - Test #5d-" << tag
                 << ": `Scalar Entries' database"
                 << "\n   Returned `float2' = " << tscalardb_float2
                 << "  , Expected = " << scalardb_float2 << endl;
      ++number_of_failures;
   }
   float tscalardb_float3 = scalardb->getFloat("float3");
   if (!tbox::MathUtilities<float>::equalEps(tscalardb_float3,
          scalardb_float3)) {
      tbox::perr << "FAILED: - Test #5e-" << tag
                 << ": `Scalar Entries' database"
                 << "\n   Returned `float3' = " << tscalardb_float3
                 << "  , Expected = " << scalardb_float3 << endl;
      ++number_of_failures;
   }

   /*
    * Tests reading scalar
    */
   double tscalardb_full_thisDouble = scalardb_full->getDouble("thisDouble");
   if (tscalardb_full_thisDouble != scalardb_full_thisDouble) {
      tbox::perr << "FAILED: - Test #6a-" << tag
                 << ": `Full' database"
                 << "\n   Returned `thisDouble' = "
                 << tscalardb_full_thisDouble
                 << "  , Expected = " << scalardb_full_thisDouble << endl;
      ++number_of_failures;
   }
   dcomplex tscalardb_full_thisComplex =
      scalardb_full->getComplex("thisComplex");
   if (tscalardb_full_thisComplex != scalardb_full_thisComplex) {
      tbox::perr << "FAILED: - Test #6b-" << tag
                 << ": `Full' database"
                 << "\n   Returned `thisComplex' = "
                 << tscalardb_full_thisComplex
                 << "  , Expected = " << scalardb_full_thisComplex << endl;
      ++number_of_failures;
   }
   int tscalardb_full_thisInt = scalardb_full->getInteger("thisInt");
   if (tscalardb_full_thisInt != scalardb_full_thisInt) {
      tbox::perr << "FAILED: - Test #6c-" << tag
                 << ": `Full' database"
                 << "\n   Returned `thisInt' = " << tscalardb_full_thisInt
                 << "  , Expected = " << scalardb_full_thisInt << endl;
      ++number_of_failures;
   }
   float tscalardb_full_thisFloat = scalardb_full->getFloat("thisFloat");
   if (!tbox::MathUtilities<float>::equalEps(tscalardb_full_thisFloat,
          scalardb_full_thisFloat)) {
      tbox::perr << "FAILED: - Test #6d-" << tag
                 << ": `Full' database"
                 << "\n   Returned `thisFloat' = " << tscalardb_full_thisFloat
                 << "  , Expected = " << scalardb_full_thisFloat << endl;
      ++number_of_failures;
   }
   bool tscalardb_full_thisBool = scalardb_full->getBool("thisBool");
   if (tscalardb_full_thisBool != scalardb_full_thisBool) {
      tbox::perr << "FAILED: - Test #6e-" << tag
                 << ": `Full' database"
                 << "\n   Returned `thisBool' = " << tscalardb_full_thisBool
                 << "  , Expected = " << scalardb_full_thisBool << endl;
      ++number_of_failures;
   }
   string tscalardb_full_thisString = scalardb_full->getString("thisString");
   if (tscalardb_full_thisString != scalardb_full_thisString) {
      tbox::perr << "FAILED: - Test #6f-" << tag
                 << ": `Full' database"
                 << "\n   Returned `thisString' = "
                 << tscalardb_full_thisString
                 << "  , Expected = " << scalardb_full_thisString << endl;
      ++number_of_failures;
   }
   char tscalardb_full_thisChar = scalardb_full->getChar("thisChar");
   if (tscalardb_full_thisChar != scalardb_full_thisChar) {
      tbox::perr << "FAILED: - Test #6g-" << tag
                 << ": `Full' database"
                 << "\n   Returned `thisChar' = " << tscalardb_full_thisChar
                 << "  , Expected = " << scalardb_full_thisChar << endl;
      ++number_of_failures;
   }
   tbox::DatabaseBox tscalardb_full_thisBox =
      scalardb_full->getDatabaseBox("thisBox");
   if (!(tscalardb_full_thisBox == scalardb_full_thisBox)) {
      tbox::perr << "FAILED: - Test #6h-" << tag
                 << ": `Full' database"
                 << "\n   Returned `thisBox' does not match Expected value"
                 << endl;
      ++number_of_failures;
   }

   /*
    * Tests for special characters in names
    */
   tscalardb_full_thisDouble = namesdb->getDouble("Name with spaces");
   if (!tbox::MathUtilities<double>::equalEps(tscalardb_full_thisDouble,
          scalardb_full_thisDouble)) {
      tbox::perr << "FAILED: - Test #6i-" << tag
                 << ": `Full' database"
                 << "\n   Returned `Name with spaces' = "
                 << tscalardb_full_thisDouble
                 << "  , Expected = " << scalardb_full_thisDouble << endl;
      ++number_of_failures;
   }

   tscalardb_full_thisDouble = namesdb->getDouble("Name-with-dashes");
   if (!tbox::MathUtilities<double>::equalEps(tscalardb_full_thisDouble,
          scalardb_full_thisDouble)) {
      tbox::perr << "FAILED: - Test #6j-" << tag
                 << ": `Full' database"
                 << "\n   Returned `Name-with-dashes' = "
                 << tscalardb_full_thisDouble
                 << "  , Expected = " << scalardb_full_thisDouble << endl;
      ++number_of_failures;
   }

   tscalardb_full_thisDouble = namesdb->getDouble("Name-with-!@#$%^&*()_+-=");
   if (tscalardb_full_thisDouble != scalardb_full_thisDouble) {
      tbox::perr << "FAILED: - Test #6k-" << tag
                 << ": `Full' database"
                 << "\n   Returned `Name-with-!@#$%^&*()_+-=' = "
                 << tscalardb_full_thisDouble
                 << "  , Expected = " << scalardb_full_thisDouble << endl;
      ++number_of_failures;
   }

   /*
    * Tests for array size
    */
   int actual_size;

   tsize = arraydb->getArraySize("ComplexArray");
   actual_size = arraydb_dcomplexArray.getSize();
   if (tsize != actual_size) {
      tbox::perr << "FAILED: - Test #7a-" << tag
                 << ": `getArraySize'"
                 << "\n   Returned size = " << tsize
                 << "  , Expected = " << actual_size << endl;
      ++number_of_failures;
   }

   tsize = arraydb->getArraySize("BoolArray");
   actual_size = arraydb_boolArray.getSize();
   if (tsize != actual_size) {
      tbox::perr << "FAILED: - Test #7b-" << tag
                 << ": `getArraySize'"
                 << "\n   Returned size = " << tsize
                 << "  , Expected = " << actual_size << endl;
      ++number_of_failures;
   }

   tsize = arraydb->getArraySize("IntArray");
   actual_size = arraydb_intArray.getSize();
   if (tsize != actual_size) {
      tbox::perr << "FAILED: - Test #7c-" << tag
                 << ": `getArraySize'"
                 << "\n   Returned size = " << tsize
                 << "  , Expected = " << actual_size << endl;
      ++number_of_failures;
   }

   tsize = arraydb->getArraySize("StringArray");
   actual_size = arraydb_stringArray.getSize();
   if (tsize != actual_size) {
      tbox::perr << "FAILED: - Test #7d-" << tag
                 << ": `getArraySize'"
                 << "\n   Returned size = " << tsize
                 << "  , Expected = " << actual_size << endl;
      ++number_of_failures;
   }

   tsize = arraydb->getArraySize("FloatArray");
   actual_size = arraydb_floatArray.getSize();
   if (tsize != actual_size) {
      tbox::perr << "FAILED: - Test #7e-" << tag
                 << ": `getArraySize'"
                 << "\n   Returned size = " << tsize
                 << "  , Expected = " << actual_size << endl;
      ++number_of_failures;
   }

   tsize = arraydb->getArraySize("DoubleArray");
   actual_size = arraydb_doubleArray.getSize();
   if (tsize != actual_size) {
      tbox::perr << "FAILED: - Test #7f-" << tag
                 << ": `getArraySize'"
                 << "\n   Returned size = " << tsize
                 << "  , Expected = " << actual_size << endl;
      ++number_of_failures;
   }

   tsize = arraydb->getArraySize("CharArray");
   actual_size = arraydb_charArray.getSize();
   if (tsize != actual_size) {
      tbox::perr << "FAILED: - Test #7g-" << tag
                 << ": `getArraySize'"
                 << "\n   Returned size = " << tsize
                 << "  , Expected = " << actual_size << endl;
      ++number_of_failures;
   }

   tsize = arraydb->getArraySize("BoxArray");
   actual_size = arraydb_boxArray.getSize();
   if (tsize != actual_size) {
      tbox::perr << "FAILED: - Test #7h-" << tag
                 << ": `getArraySize'"
                 << "\n   Returned size = " << tsize
                 << "  , Expected = " << actual_size << endl;
      ++number_of_failures;
   }

   tsize = db->getArraySize("Array Entries");
   actual_size = 0;
   if (tsize != actual_size) {
      tbox::perr << "FAILED: - Test #7i-" << tag
                 << ": `getArraySize'"
                 << "\n   Returned size = " << tsize
                 << "  , Expected = " << actual_size << endl;
      ++number_of_failures;
   }

   /*
    * Tests for existance of each type
    */
   if (!arraydb->keyExists("ComplexArray")) {
      tbox::perr << "FAILED: - Test #8a-" << tag
                 << ": `keyExists'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }

   if (!arraydb->keyExists("BoolArray")) {
      tbox::perr << "FAILED: - Test #8b-" << tag
                 << ": `keyExists'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }
   if (!arraydb->keyExists("IntArray")) {
      tbox::perr << "FAILED: - Test #8c-" << tag
                 << ": `keyExists'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }
   if (!arraydb->keyExists("StringArray")) {
      tbox::perr << "FAILED: - Test #8d-" << tag
                 << ": `keyExists'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }
   if (!arraydb->keyExists("FloatArray")) {
      tbox::perr << "FAILED: - Test #8d-" << tag
                 << ": `keyExists'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }
   if (!arraydb->keyExists("DoubleArray")) {
      tbox::perr << "FAILED: - Test #8f-" << tag
                 << ": `keyExists'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }
   if (!arraydb->keyExists("CharArray")) {
      tbox::perr << "FAILED: - Test #8g-" << tag
                 << ": `keyExists'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }
   if (!arraydb->keyExists("BoxArray")) {
      tbox::perr << "FAILED: - Test #8h-" << tag
                 << ": `keyExists'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }
   if (!db->keyExists("Array Entries")) {
      tbox::perr << "FAILED: - Test #8i-" << tag
                 << ": `keyExists'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }

   /*
    * Tests for non-existant key
    */
   if (arraydb->keyExists("NonExistantKey")) {
      tbox::perr << "FAILED: - Test #10-" << tag
                 << ": `keyExists'"
                 << "\n   Returned true "
                 << "  , Expected false " << endl;
      ++number_of_failures;
   }

   /*
    * Tests for is methods
    */
   if (!arraydb->isComplex("ComplexArray")) {
      tbox::perr << "FAILED: - Test #11a-" << tag
                 << ": `isComplex'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }

   if (!arraydb->isBool("BoolArray")) {
      tbox::perr << "FAILED: - Test #11b-" << tag
                 << ": `isBool'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }
   if (!arraydb->isInteger("IntArray")) {
      tbox::perr << "FAILED: - Test #11c-" << tag
                 << ": `isInteger'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }
   if (!arraydb->isString("StringArray")) {
      tbox::perr << "FAILED: - Test #11d-" << tag
                 << ": `isString'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }
   if (!arraydb->isFloat("FloatArray")) {
      tbox::perr << "FAILED: - Test #11d-" << tag
                 << ": `isFloat'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }
   if (!arraydb->isDouble("DoubleArray")) {
      tbox::perr << "FAILED: - Test #11f-" << tag
                 << ": `isDouble'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }
   if (!arraydb->isChar("CharArray")) {
      tbox::perr << "FAILED: - Test #11g-" << tag
                 << ": `isChar'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }
   if (!arraydb->isDatabaseBox("BoxArray")) {
      tbox::perr << "FAILED: - Test #11h-" << tag
                 << ": `isBox'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }
   if (!db->isDatabase("Array Entries")) {
      tbox::perr << "FAILED: - Test #11i-" << tag
                 << ": `isDatabase'"
                 << "\n   Returned false "
                 << "  , Expected = true " << endl;
      ++number_of_failures;
   }

   /*
    * Tests for is methods on things that are not
    */
   if (arraydb->isComplex("nonComplexArray")) {
      tbox::perr << "FAILED: - Test #12a-" << tag
                 << ": `isComplex'"
                 << "\n   Returned true "
                 << "  , Expected false " << endl;
      ++number_of_failures;
   }

   if (arraydb->isBool("nonBoolArray")) {
      tbox::perr << "FAILED: - Test #12b-" << tag
                 << ": `isBool'"
                 << "\n   Returned true "
                 << "  , Expected false " << endl;
      ++number_of_failures;
   }
   if (arraydb->isInteger("nonIntArray")) {
      tbox::perr << "FAILED: - Test #12c-" << tag
                 << ": `isInteger'"
                 << "\n   Returned true "
                 << "  , Expected false " << endl;
      ++number_of_failures;
   }
   if (arraydb->isString("nonStringArray")) {
      tbox::perr << "FAILED: - Test #12d-" << tag
                 << ": `isString'"
                 << "\n   Returned true "
                 << "  , Expected false " << endl;
      ++number_of_failures;
   }
   if (arraydb->isFloat("nonFloatArray")) {
      tbox::perr << "FAILED: - Test #12e-" << tag
                 << ": `isFloat'"
                 << "\n   Returned true "
                 << "  , Expected false " << endl;
      ++number_of_failures;
   }
   if (arraydb->isDouble("nonDoubleArray")) {
      tbox::perr << "FAILED: - Test #12f-" << tag
                 << ": `isDouble'"
                 << "\n   Returned true "
                 << "  , Expected false " << endl;
      ++number_of_failures;
   }
   if (arraydb->isChar("nonCharArray")) {
      tbox::perr << "FAILED: - Test #12g-" << tag
                 << ": `isChar'"
                 << "\n   Returned true "
                 << "  , Expected false " << endl;
      ++number_of_failures;
   }
   if (arraydb->isDatabaseBox("nonBoxArray")) {
      tbox::perr << "FAILED: - Test #12h-" << tag
                 << ": `isBox'"
                 << "\n   Returned true "
                 << "  , Expected false " << endl;
      ++number_of_failures;
   }
   if (db->isDatabase("nonArray Entries")) {
      tbox::perr << "FAILED: - Test #12i-" << tag
                 << ": `isDatabase'"
                 << "\n   Returned true "
                 << "  , Expected false " << endl;
      ++number_of_failures;
   }

#if 0

   /*
    * Tests for getArrayType
    */
   if (arraydb->getArrayType("ComplexArray") != tbox::Database::COMPLEX) {
      tbox::perr << "FAILED: - Test #13a-" << tag
                 << ": `getArrayType returned the wrong type" << endl;
      ++number_of_failures;
   }
   if (arraydb->getArrayType("BoolArray") != tbox::Database::BOOL) {
      tbox::perr << "FAILED: - Test #13b-" << tag
                 << ": `getArrayType returned the wrong type" << endl;
      ++number_of_failures;
   }
   if (arraydb->getArrayType("IntArray") != tbox::Database::INT) {
      tbox::perr << "FAILED: - Test #13c-" << tag
                 << ": `getArrayType returned the wrong type" << endl;
      ++number_of_failures;
   }
   if (arraydb->getArrayType("StringArray") != tbox::Database::STRING) {
      tbox::perr << "FAILED: - Test #13d-" << tag
                 << ": `getArrayType returned the wrong type" << endl;
      ++number_of_failures;
   }
   if (arraydb->getArrayType("FloatArray") != tbox::Database::FLOAT) {
      tbox::perr << "FAILED: - Test #13e-" << tag
                 << ": `getArrayType returned the wrong type" << endl;
      ++number_of_failures;
   }
   if (arraydb->getArrayType("DoubleArray") != tbox::Database::DOUBLE) {
      tbox::perr << "FAILED: - Test #13f-" << tag
                 << ": `getArrayType returned the wrong type" << endl;
      ++number_of_failures;
   }
   if (arraydb->getArrayType("CharArray") != tbox::Database::CHAR) {
      tbox::perr << "FAILED: - Test #13g-" << tag
                 << ": `getArrayType returned the wrong type" << endl;
      ++number_of_failures;
   }
   if (arraydb->getArrayType("BoxArray") != tbox::Database::BOX) {
      tbox::perr << "FAILED: - Test #13h-" << tag
                 << ": `getArrayType returned the wrong type" << endl;
      ++number_of_failures;
   }
   if (db->getArrayType("Array Entries") != tbox::Database::DATABASE) {
      tbox::perr << "FAILED: - Test #13i-" << tag
                 << ": `getArrayType returned the wrong type" << endl;
      ++number_of_failures;
   }
   if (db->getArrayType("NonEntries") != tbox::Database::INVALID) {
      tbox::perr << "FAILED: - Test #13j-" << tag
                 << ": `getArrayType returned the wrong type" << endl;
      ++number_of_failures;
   }

#endif

   /*
    * Tests reading scalar with default
    */
   dcomplex test_scalar_complex(-9.0, -9.0);
   tscalardb_full_thisComplex = scalardb_full->getComplexWithDefault(
         "thisComplex",
         test_scalar_complex);
   if (tscalardb_full_thisComplex != scalardb_full_thisComplex) {
      tbox::perr << "FAILED: - Test #14a-" << tag
                 << ": `Full' database" << endl
                 << "   Returned `thisComplex' = "
                 << tscalardb_full_thisComplex
                 << "  , Expected = " << scalardb_full_thisComplex << endl;
      ++number_of_failures;
   }
   bool test_scalar_bool(false);
   tscalardb_full_thisBool = scalardb_full->getBoolWithDefault("thisBool",
         test_scalar_bool);
   if (tscalardb_full_thisBool != scalardb_full_thisBool) {
      tbox::perr << "FAILED: - Test #14b-" << tag
                 << ": `Full' database" << endl
                 << "   Returned `thisBool' = " << tscalardb_full_thisBool
                 << "  , Expected = " << scalardb_full_thisBool << endl;
      ++number_of_failures;
   }
   int test_scalar_int(-9);
   tscalardb_full_thisInt = scalardb_full->getIntegerWithDefault("thisInt",
         test_scalar_int);
   if (tscalardb_full_thisInt != scalardb_full_thisInt) {
      tbox::perr << "FAILED: - Test #14c-" << tag
                 << ": `Full' database" << endl
                 << "   Returned `thisInt' = " << tscalardb_full_thisInt
                 << "  , Expected = " << scalardb_full_thisInt << endl;
      ++number_of_failures;
   }
   string test_scalar_string("A fake string");
   tscalardb_full_thisString = scalardb_full->getStringWithDefault("thisString",
         test_scalar_string);
   if (tscalardb_full_thisString != scalardb_full_thisString) {
      tbox::perr << "FAILED: - Test #14d-" << tag
                 << ": `Full' database" << endl
                 << "   Returned `thisString' = " << tscalardb_full_thisString
                 << "  , Expected = " << scalardb_full_thisString << endl;
      ++number_of_failures;
   }
   float test_scalar_float(-9.0);
   tscalardb_full_thisFloat = scalardb_full->getFloatWithDefault("thisFloat",
         test_scalar_float);
   if (!tbox::MathUtilities<float>::equalEps(tscalardb_full_thisFloat,
          scalardb_full_thisFloat)) {
      tbox::perr << "FAILED: - Test #14e-" << tag
                 << ": `Full' database" << endl
                 << "   Returned `thisFloat' = " << tscalardb_full_thisFloat
                 << "  , Expected = " << scalardb_full_thisFloat << endl;
      ++number_of_failures;
   }
   double test_scalar_double(-9.0);
   tscalardb_full_thisDouble = scalardb_full->getDoubleWithDefault("thisDouble",
         test_scalar_double);
   if (!tbox::MathUtilities<double>::equalEps(tscalardb_full_thisDouble,
          scalardb_full_thisDouble)) {
      tbox::perr << "FAILED: - Test #14f-" << tag
                 << ": `Full' database" << endl
                 << "   Returned `thisDouble' = " << tscalardb_full_thisDouble
                 << "  , Expected = " << scalardb_full_thisDouble << endl;
      ++number_of_failures;
   }
   char test_scalar_char('h');
   tscalardb_full_thisChar = scalardb_full->getCharWithDefault("thisChar",
         test_scalar_char);
   if (tscalardb_full_thisChar != scalardb_full_thisChar) {
      tbox::perr << "FAILED: - Test #14h-" << tag
                 << ": `Full' database" << endl
                 << "   Returned `thisChar' = " << tscalardb_full_thisChar
                 << "  , Expected = " << scalardb_full_thisChar << endl;
      ++number_of_failures;
   }

   tscalardb_full_thisComplex = defaultdb->getComplexWithDefault("bogusComplex",
         test_scalar_complex);
   if (tscalardb_full_thisComplex != test_scalar_complex) {
      tbox::perr << "FAILED: - Test #15a-" << tag
                 << ": `Full' database" << endl
                 << "   Returned `thisComplex' = "
                 << tscalardb_full_thisComplex
                 << "  , Expected = " << test_scalar_complex << endl;
      ++number_of_failures;
   }
   tscalardb_full_thisBool = defaultdb->getBoolWithDefault("bogusBool",
         test_scalar_bool);
   if (tscalardb_full_thisBool != test_scalar_bool) {
      tbox::perr << "FAILED: - Test #15b-" << tag
                 << ": `Full' database" << endl
                 << "   Returned `thisBool' = " << tscalardb_full_thisBool
                 << "  , Expected = " << test_scalar_bool << endl;
      ++number_of_failures;
   }

   tscalardb_full_thisInt = defaultdb->getIntegerWithDefault("bogusInt",
         test_scalar_int);
   if (tscalardb_full_thisInt != test_scalar_int) {
      tbox::perr << "FAILED: - Test #15c-" << tag
                 << ": `Full' database" << endl
                 << "   Returned `thisInt' = " << tscalardb_full_thisInt
                 << "  , Expected = " << test_scalar_int << endl;
      ++number_of_failures;
   }
   tscalardb_full_thisString = defaultdb->getStringWithDefault("bogusString",
         test_scalar_string);
   if (tscalardb_full_thisString != test_scalar_string) {
      tbox::perr << "FAILED: - Test #15d-" << tag
                 << ": `Full' database" << endl
                 << "   Returned `thisString' = " << tscalardb_full_thisString
                 << "  , Expected = " << test_scalar_string << endl;
      ++number_of_failures;
   }
   tscalardb_full_thisFloat = defaultdb->getFloatWithDefault("bogusFloat",
         test_scalar_float);
   if (!tbox::MathUtilities<float>::equalEps(tscalardb_full_thisFloat,
          test_scalar_float)) {
      tbox::perr << "FAILED: - Test #15e-" << tag
                 << ": `Full' database" << endl
                 << "   Returned `thisFloat' = " << tscalardb_full_thisFloat
                 << "  , Expected = " << test_scalar_float << endl;
      ++number_of_failures;
   }
   tscalardb_full_thisDouble = defaultdb->getDoubleWithDefault("bogusDouble",
         test_scalar_double);
   if (!tbox::MathUtilities<double>::equalEps(tscalardb_full_thisDouble,
          test_scalar_double)) {
      tbox::perr << "FAILED: - Test #15f-" << tag
                 << ": `Full' database" << endl
                 << "   Returned `thisDouble' = " << tscalardb_full_thisDouble
                 << "  , Expected = " << test_scalar_double << endl;
      ++number_of_failures;
   }
   tscalardb_full_thisChar = defaultdb->getCharWithDefault("bogusChar",
         test_scalar_char);
   if (tscalardb_full_thisChar != test_scalar_char) {
      tbox::perr << "FAILED: - Test #15h-" << tag
                 << ": `Full' database" << endl
                 << "   Returned `thisChar' = " << tscalardb_full_thisChar
                 << "  , Expected = " << test_scalar_char << endl;
      ++number_of_failures;
   }

   /*
    * Tests for reading stl::vector
    */
   std::vector<hier::IntVector> vector_IntVector(2);

   vectordb->getVector("vector_IntVector", vector_IntVector);

   if (vector_IntVector[0] != intVector1) {
      tbox::perr << "FAILED: - Test #16a-" << tag
                 << ": stl::vector<IntVector> did not restore correctly" << endl
                 << "   Returned `IntVector' = " << vector_IntVector[0]
                 << "  , Expected = " << intVector1 << endl;
      ++number_of_failures;

   }

   if (vector_IntVector[1] != intVector2) {
      tbox::perr << "FAILED: - Test #16b-" << tag
                 << ": stl::vector<IntVector> did not restore correctly" << endl
                 << "   Returned `IntVector' = " << vector_IntVector[1]
                 << "  , Expected = " << intVector2 << endl;
      ++number_of_failures;
   }

}
