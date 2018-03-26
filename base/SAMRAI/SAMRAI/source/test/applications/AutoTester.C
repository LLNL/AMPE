/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   (c) 1997-2012 Lawrence Livermore National Security, LLC
 *                Description:   Class used for auto testing applications
 *
 ************************************************************************/

#include "AutoTester.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/MathUtilities.h"

AutoTester::AutoTester(
   const std::string& object_name,
   const tbox::Dimension& dim,
   boost::shared_ptr<tbox::Database> input_db):
   d_dim(dim)
#ifdef HAVE_HDF5
   ,
   d_hdf_db("AutoTesterDatabase")
#endif
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   d_object_name = object_name;
   d_test_fluxes = false;
   d_test_iter_num = 10;
   d_output_correct = false;

   d_write_patch_boxes = false;
   d_read_patch_boxes = false;
   d_test_patch_boxes_at_steps.resizeArray(0);
   d_test_patch_boxes_step_count = 0;

   getFromInput(input_db);

   const std::string hdf_filename =
      d_test_patch_boxes_filename
      + "." + tbox::Utilities::nodeToString(mpi.getSize())
      + "." + tbox::Utilities::processorToString(mpi.getRank());

#ifdef HAVE_HDF5
   if (d_read_patch_boxes) {
      d_hdf_db.open(hdf_filename);
      if (d_output_correct) {
         d_hdf_db.printClassData(tbox::pout);
      }

   }

   if (d_write_patch_boxes) {
      d_hdf_db.create(hdf_filename);
   }
#endif

}

AutoTester::~AutoTester()
{
}

/*
 ******************************************************************
 *
 *  Method "evalTestData" compares the result of the run with
 *  the correct result for runs with the TimeRefinementIntegrator
 *  and HyperbolicLevelIntegrator.
 *
 ******************************************************************
 */
int AutoTester::evalTestData(
   int iter,
   const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   const boost::shared_ptr<algs::TimeRefinementIntegrator> tri,
   const boost::shared_ptr<algs::HyperbolicLevelIntegrator> hli,
   const boost::shared_ptr<mesh::GriddingAlgorithm> ga)
{
   NULL_USE(ga);

   int num_failures = 0;

   /*
    * Compare "correct_result" array to the computed result on specified
    * iteration.
    */
   if (iter == d_test_iter_num && !d_test_fluxes) {

      /*
       * set precision of output stream.
       */
      tbox::plog.precision(12);

      /*
       * determine level.
       */
      int nlevels = hierarchy->getNumberOfLevels() - 1;
      boost::shared_ptr<hier::PatchLevel> level(
         hierarchy->getPatchLevel(nlevels));

      /*
       * Test 0: Time Refinement Integrator
       */
      double time = tri->getIntegratorTime();
      if (d_correct_result.getSize() > 0) {
         if (d_output_correct) {
            tbox::plog << "Test 0: Time Refinement Integrator "
                       << "\n   computed result: " << time;

            tbox::plog << "\n   specified result = "
                       << d_correct_result[0];
         }
         tbox::plog << std::endl;

         if (tbox::MathUtilities<double>::equalEps(time,
                d_correct_result[0])) {
            tbox::plog << "Test 0: Time Refinement check successful"
                       << std::endl;
         } else {
            tbox::perr << "Test 0 FAILED: Check Time Refinement Integrator"
                       << std::endl;
            num_failures++;
         }
      }

      /*
       * Test 1: Time Refinement Integrator
       */
      double dt = tri->getLevelDtMax(nlevels);
      if (d_correct_result.getSize() > 1) {
         if (d_output_correct) {
            tbox::plog << "Test 1: Time Refinement Integrator "
                       << "\n   computed result: " << dt;
            tbox::plog << "\n   specified result = "
                       << d_correct_result[1];
         }
         tbox::plog << std::endl;

         if (tbox::MathUtilities<double>::equalEps(dt, d_correct_result[1])) {
            tbox::plog << "Test 1: Time Refinement check successful"
                       << std::endl;
         } else {
            tbox::perr << "Test 1 FAILED: Check Time Refinement Integrator"
                       << std::endl;
            num_failures++;
         }
      }

      /*
       * Test 2: Hyperbolic Level Integrator
       */
      dt = hli->getLevelDt(level, time, false);
      if (d_correct_result.getSize() > 2) {
         if (d_output_correct) {
            tbox::plog << "Test 2: Hyperbolic Level Integrator "
                       << "\n   computed result: " << dt;

            tbox::plog << "\n   specified result = "
                       << d_correct_result[2];
         }
         tbox::plog << std::endl;

         if (tbox::MathUtilities<double>::equalEps(dt, d_correct_result[2])) {
            tbox::plog << "Test 2: Hyperbolic Level Int check successful"
                       << std::endl;
         } else {
            tbox::perr << "Test 2 FAILED: Check Hyperbolic Level Integrator"
                       << std::endl;
            num_failures++;
         }
      }

      /*
       * Test 3: Gridding Algorithm
       */
      int n = hierarchy->getMaxNumberOfLevels();
      if (d_output_correct) {
         tbox::plog << "Test 3: Gridding Algorithm "
                    << "\n   computed result: " << n;
         tbox::plog << "\n   correct result = " << nlevels + 1;
         tbox::plog << std::endl;
      }
      if (n == (nlevels + 1)) {
         tbox::plog << "Test 3: Gridding Algorithm check successful"
                    << std::endl;
      } else {
         tbox::perr << "Test 3 FAILED: Check Gridding Algorithm" << std::endl;
         num_failures++;
      }

   }

   if ((d_test_patch_boxes_at_steps.getSize() >
        d_test_patch_boxes_step_count) &&
       (d_test_patch_boxes_at_steps[d_test_patch_boxes_step_count] == iter)) {

      int num_levels = hierarchy->getNumberOfLevels();

#ifdef HAVE_HDF5
      if (d_read_patch_boxes) {

         if (d_output_correct) {
            d_hdf_db.printClassData(tbox::pout);
         }

         const std::string step_name =
            std::string("step_number_") + tbox::Utilities::intToString(
               d_test_patch_boxes_step_count,
               2);
         std::cout << std::endl;
         boost::shared_ptr<tbox::Database> step_db(
            d_hdf_db.getDatabase(step_name));

         /*
          * FIXME: This check give false positives!!!!!
          * It writes the same file regardless of the number of processors.
          * We should be checking against base runs with the same number of processors,
          * compare different data.
          */
         for (int ln = 0; ln < num_levels; ln++) {

            const std::string level_name =
               std::string("level_number_") + tbox::Utilities::levelToString(ln);
            boost::shared_ptr<tbox::Database> level_db(
               step_db->getDatabase(level_name));
            hier::BoxLevel correct_mapped_box_level(d_dim);
            boost::shared_ptr<const hier::BaseGridGeometry> grid_geometry(
               hierarchy->getGridGeometry());
            correct_mapped_box_level.getFromDatabase(*level_db,
               grid_geometry);

            num_failures += checkHierarchyBoxes(hierarchy,
                  ln,
                  correct_mapped_box_level,
                  iter);
         }

      }

      if (d_write_patch_boxes) {

         const std::string step_name =
            std::string("step_number_") + tbox::Utilities::intToString(
               d_test_patch_boxes_step_count,
               2);
         std::cout << std::endl;
         boost::shared_ptr<tbox::Database> step_db(
            d_hdf_db.putDatabase(step_name));

         for (int ln = 0; ln < num_levels; ln++) {
            boost::shared_ptr<hier::PatchLevel> level(
               hierarchy->getPatchLevel(ln));

            const std::string level_name =
               std::string("level_number_") + tbox::Utilities::levelToString(ln);
            boost::shared_ptr<tbox::Database> level_db(
               step_db->putDatabase(level_name));
            level->getBoxLevel()->putUnregisteredToDatabase(level_db);
         }

         if (d_output_correct) {
            d_hdf_db.printClassData(tbox::pout);
         }
      }
#endif

      d_test_patch_boxes_step_count++;

   }

   return num_failures;
}

/*
 ******************************************************************
 *
 *  Method "evalTestData" compares the result of the run with
 *  the correct result for runs with the MethodOfLinesIntegrator.
 *
 ******************************************************************
 */
int AutoTester::evalTestData(
   int iter,
   const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   double time,
   const boost::shared_ptr<algs::MethodOfLinesIntegrator> mol,
   const boost::shared_ptr<mesh::GriddingAlgorithm> ga)
{
   NULL_USE(ga);

   int num_failures = 0;

   /*
    * Compare "correct_result" array to the computed result on specified
    * iteration.
    */
   if (iter == d_test_iter_num && !d_test_fluxes) {

      /*
       * set precision of output stream.
       */
      tbox::plog.precision(12);

      /*
       * determine level.
       */
      int nlevels = hierarchy->getNumberOfLevels() - 1;
      boost::shared_ptr<hier::PatchLevel> level(
         hierarchy->getPatchLevel(nlevels));

      /*
       * Test 0: Time test
       */
      if (d_output_correct) {
         tbox::plog << "Test 0: Simulation Time: "
                    << "\n   computed result: " << time;
         if (d_correct_result.getSize() > 0) {
            tbox::plog << "\n   specified result = "
                       << d_correct_result[0];
         }
         tbox::plog << std::endl;
      }
      if (tbox::MathUtilities<double>::equalEps(time, d_correct_result[0])) {
         tbox::plog << "Test 0: Simulation Time check successful" << std::endl;
      } else {
         tbox::perr << "Test 0 FAILED: Simulation time incorrect" << std::endl;
         num_failures++;
      }

      /*
       * Test 1: MethodOfLinesIntegrator
       */
      double dt = mol->getTimestep(hierarchy, time);
      if (d_output_correct) {
         tbox::plog << "Test 1: Method of Lines Integrator "
                    << "\n   computed result: " << dt;
         if (d_correct_result.getSize() > 1) {
            tbox::plog << "\n   specified result = "
                       << d_correct_result[1];
         }
         tbox::plog << std::endl;
      }
      if (tbox::MathUtilities<double>::equalEps(dt, d_correct_result[1])) {
         tbox::plog << "Test 1: MOL Int check successful" << std::endl;
      } else {
         tbox::perr << "Test 1 FAILED: Check Method of Lines Integrator"
                    << std::endl;
         num_failures++;
      }

      /*
       * Test 2: Gridding Algorithm
       */
      int n = hierarchy->getMaxNumberOfLevels();
      if (d_output_correct) {
         tbox::plog << "Test 2: Gridding Algorithm "
                    << "\n   computed result: " << n;
         tbox::plog << "\n   correct result = " << nlevels + 1;
         tbox::plog << std::endl;
      }
      if (n == (nlevels + 1)) {
         tbox::plog << "Test 2: Gridding Alg check successful" << std::endl;
      } else {
         tbox::perr << "Test 2 FAILED: Check Gridding Algorithm" << std::endl;
         num_failures++;
      }

   }

   if ((d_test_patch_boxes_at_steps.getSize() > 0) &&
       (d_test_patch_boxes_at_steps[d_test_patch_boxes_step_count] == iter)) {

      int num_levels = hierarchy->getNumberOfLevels();

#ifdef HAVE_HDF5
      if (d_read_patch_boxes) {

         if (d_output_correct) {
            d_hdf_db.printClassData(tbox::pout);
         }

         const std::string step_name =
            std::string("step_number_") + tbox::Utilities::intToString(
               d_test_patch_boxes_step_count,
               2);
         std::cout << std::endl;
         boost::shared_ptr<tbox::Database> step_db(
            d_hdf_db.getDatabase(step_name));

         for (int ln = 0; ln < num_levels; ln++) {

            const std::string level_name =
               std::string("level_number_") + tbox::Utilities::levelToString(ln);
            boost::shared_ptr<tbox::Database> level_db(
               step_db->getDatabase(level_name));
            hier::BoxLevel correct_mapped_box_level(d_dim);
            boost::shared_ptr<const hier::BaseGridGeometry> grid_geometry(
               hierarchy->getGridGeometry());
            correct_mapped_box_level.getFromDatabase(*level_db,
               grid_geometry);

            num_failures += checkHierarchyBoxes(hierarchy,
                  ln,
                  correct_mapped_box_level,
                  iter);
         }

      }

      if (d_write_patch_boxes) {

         if (d_output_correct) {
            d_hdf_db.printClassData(tbox::pout);
         }

         const std::string step_name =
            std::string("step_number_") + tbox::Utilities::intToString(
               d_test_patch_boxes_step_count,
               2);
         std::cout << std::endl;
         boost::shared_ptr<tbox::Database> step_db(
            d_hdf_db.putDatabase(step_name));

         for (int ln = 0; ln < num_levels; ln++) {
            boost::shared_ptr<hier::PatchLevel> level(
               hierarchy->getPatchLevel(ln));

            const std::string level_name =
               std::string("level_number_") + tbox::Utilities::levelToString(ln);
            boost::shared_ptr<tbox::Database> level_db(
               step_db->putDatabase(level_name));
            level->getBoxLevel()->putUnregisteredToDatabase(level_db);
         }

      }
#endif

      d_test_patch_boxes_step_count++;

   }

   return num_failures;
}

/*
 ******************************************************************
 *
 *  Get test parameters from input.
 *
 ******************************************************************
 */

void AutoTester::getFromInput(
   boost::shared_ptr<tbox::Database> input_db)
{
   boost::shared_ptr<tbox::Database> tester_db(
      input_db->getDatabase(d_object_name));

   /*
    * Read testing parameters from testing_db
    */
   if (tester_db->keyExists("test_fluxes")) {
      d_test_fluxes = tester_db->getBool("test_fluxes");
   }

   if (tester_db->keyExists("test_iter_num")) {
      d_test_iter_num = tester_db->getInteger("test_iter_num");
   }

   if (tester_db->keyExists("write_patch_boxes")) {
      d_write_patch_boxes = tester_db->getBool("write_patch_boxes");
   }
   if (tester_db->keyExists("read_patch_boxes")) {
      d_read_patch_boxes = tester_db->getBool("read_patch_boxes");
   }
   if (d_read_patch_boxes && d_write_patch_boxes) {
      tbox::perr << "FAILED: - AutoTester " << d_object_name << "\n"
                 << "Cannot 'read_patch_boxes' and 'write_patch_boxes' \n"
                 << "at the same time." << std::endl;
   }
   if (d_read_patch_boxes || d_write_patch_boxes) {
      if (!tester_db->keyExists("test_patch_boxes_at_steps")) {
         tbox::perr << "FAILED: - AutoTester " << d_object_name << "\n"
                    << "Must provide 'test_patch_boxes_at_steps' data."
                    << std::endl;
      } else {
         d_test_patch_boxes_at_steps =
            tester_db->getIntegerArray("test_patch_boxes_at_steps");
      }
      if (!tester_db->keyExists("test_patch_boxes_filename")) {
         tbox::perr << "FAILED: - AutoTester " << d_object_name << "\n"
                    << "Must provide 'test_patch_boxes_filename' data."
                    << std::endl;
      } else {
         d_test_patch_boxes_filename =
            tester_db->getString("test_patch_boxes_filename");
      }
   }

   if (d_test_fluxes) {

      /*
       * Read expected result for flux test...
       * Fluxes not verified in this routine.  Rather, we let it
       * write the result and do a "diff" within the script
       */

      tbox::pout << "Do a diff on the resulting *.dat file to verify result."
                 << std::endl;

   } else {

      /*
       * Read correct_result array for timestep test...
       */
      if (tester_db->keyExists("correct_result")) {
         d_correct_result = tester_db->getDoubleArray("correct_result");
      } else {
         TBOX_WARNING("main.C: TESTING is on but no `correct_result' array"
            << "is given in input file." << std::endl);
      }

      /* Specify whether to output "correct_result" result */

      if (tester_db->keyExists("output_correct")) {
         d_output_correct = tester_db->getBool("output_correct");
      }

   }

}

int AutoTester::checkHierarchyBoxes(
   const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number,
   const hier::BoxLevel& correct_mapped_box_level,
   int iter)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   const boost::shared_ptr<hier::PatchLevel> patch_level(
      hierarchy->getPatchLevel(level_number));
   const hier::BoxLevel& mapped_box_level =
      *patch_level->getBoxLevel();

   const int local_exact_match =
      mapped_box_level == correct_mapped_box_level;

   int global_exact_match = local_exact_match;
   if (mpi.getSize() > 1) {
      mpi.AllReduce(&global_exact_match, 1, MPI_MIN);
   }

   /*
    * Check to make sure hierarchy's BoxLevel and
    * correct_mapped_box_level are identical.  If not, write an error
    * message.
    */

   int num_failures = 0;

   if (local_exact_match && global_exact_match) {
      tbox::plog << "Test 4: Level " << level_number
                 << " BoxLevel check successful for step " << iter
                 << std::endl << std::endl;
   } else {
      tbox::perr << "Test 4: FAILED: Level " << level_number
                 << " hier::BoxLevel configuration doesn't match at step " << iter
                 << std::endl << std::endl;
      num_failures++;
   }

   if (d_output_correct) {

      tbox::pout << "-------------------------------------------------------"
                 << std::endl;

      if (!local_exact_match) {
         tbox::pout << "LOCAL MAPPED BOX LEVEL DOES NOT MATCH "
                    << "ON LEVEL: " << level_number << std::endl;
      }

      if (!global_exact_match) {
         tbox::pout << "GLOBAL MAPPED BOX LEVEL DOES NOT MATCH "
                    << "ON LEVEL: " << level_number << std::endl;
      }
      tbox::pout << "BoxLevel: " << std::endl;
      mapped_box_level.recursivePrint(tbox::pout, "", 3);
      tbox::pout << "correct BoxLevel: " << std::endl;
      correct_mapped_box_level.recursivePrint(tbox::pout, "", 3);

      tbox::pout << "-------------------------------------------------------"
                 << std::endl << std::endl;

   }

   return num_failures;
}
