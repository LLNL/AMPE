/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   AMR communication tests for outerside-centered patch data
 *
 ************************************************************************/

#include "OutersideDataTest.h"

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "CommTester.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include "SAMRAI/pdat/SideIndex.h"
#include "SAMRAI/pdat/SideIterator.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/pdat/OutersideGeometry.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

namespace SAMRAI {

using namespace std;

OutersideDataTest::OutersideDataTest(
   const string& object_name,
   const tbox::Dimension& dim,
   boost::shared_ptr<tbox::Database> main_input_db,
   bool do_refine,
   bool do_coarsen,
   const string& refine_option):
   PatchDataTestStrategy(dim),
   d_dim(dim)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(main_input_db);
   TBOX_ASSERT(!refine_option.empty());
#endif

   d_object_name = object_name;

   d_do_refine = do_refine;
   d_do_coarsen = false;
   if (!do_refine) {
      d_do_coarsen = do_coarsen;
   }
   if (d_do_refine) {
      TBOX_ERROR("There is no refine test for Outerside data type, because\n"
         << "Outerside refinement does not exist at this time.");
      /*
       * The refine codes are still kept in this class in case we
       * implement Outerside refinement in the future.
       */
   }

   d_refine_option = refine_option;

   d_use_fine_value_at_interface.resizeArray(0);

   d_Acoef = 0.0;
   d_Bcoef = 0.0;
   d_Ccoef = 0.0;
   d_Dcoef = 0.0;

   d_finest_level_number = main_input_db->
      getDatabase("PatchHierarchy")->
      getInteger("max_levels") - 1;

   d_cart_grid_geometry.reset(
      new geom::CartesianGridGeometry(
         dim,
         "CartesianGridGeometry",
         main_input_db->getDatabase("CartesianGridGeometry")));

   setGridGeometry(d_cart_grid_geometry);
   readTestInput(main_input_db->getDatabase("OutersidePatchDataTest"));

}

OutersideDataTest::~OutersideDataTest()
{
}

void OutersideDataTest::readTestInput(
   boost::shared_ptr<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(db);
#endif

   /*
    * Base class reads variable parameters and boxes to refine.
    */

   readVariableInput(db->getDatabase("VariableData"));
   readRefinementInput(db->getDatabase("RefinementData"));

   boost::shared_ptr<tbox::Database> var_data(db->getDatabase("VariableData"));
   tbox::Array<string> var_keys = var_data->getAllKeys();
   int nkeys = var_keys.getSize();

   d_use_fine_value_at_interface.resizeArray(nkeys);

   for (int i = 0; i < nkeys; i++) {
      boost::shared_ptr<tbox::Database> var_db(
         var_data->getDatabase(var_keys[i]));

      if (var_db->keyExists("use_fine_value_at_interface")) {
         d_use_fine_value_at_interface[i] =
            var_db->getBool("use_fine_value_at_interface");
      } else {
         d_use_fine_value_at_interface[i] = true;
      }

   }

   if (db->keyExists("Acoef")) {
      d_Acoef = db->getDouble("Acoef");
   } else {
      TBOX_ERROR(d_object_name << " input error: No `Acoeff' found." << endl);
   }
   if (db->keyExists("Dcoef")) {
      d_Dcoef = db->getDouble("Dcoef");
   } else {
      TBOX_ERROR(d_object_name << " input error: No `Dcoef' found." << endl);
   }
   if (d_dim > tbox::Dimension(1)) {
      if (db->keyExists("Bcoef")) {
         d_Bcoef = db->getDouble("Bcoef");
      } else {
         TBOX_ERROR(d_object_name << " input error: No `Bcoef' found." << endl);
      }
   }
   if (d_dim > tbox::Dimension(2)) {
      if (db->keyExists("Ccoef")) {
         d_Ccoef = db->getDouble("Ccoef");
      } else {
         TBOX_ERROR(d_object_name << " input error: No `Ccoef' found." << endl);
      }
   }

}

void OutersideDataTest::registerVariables(
   CommTester* commtest)
{
   TBOX_ASSERT(commtest != (CommTester *)NULL);

   int nvars = d_variable_src_name.getSize();

   d_variables_src.resizeArray(nvars);
   d_variables_dst.resizeArray(nvars);

   for (int i = 0; i < nvars; i++) {
      d_variables_src[i].reset(
         new pdat::OutersideVariable<double>(
            d_dim,
            d_variable_src_name[i],
            d_variable_depth[i]));
      d_variables_dst[i].reset(
         new pdat::SideVariable<double>(
            d_dim,
            d_variable_dst_name[i],
            d_variable_depth[i],
            d_use_fine_value_at_interface[i]));

      if (d_do_refine) {
         commtest->registerVariable(d_variables_src[i],
            d_variables_dst[i],
            d_variable_src_ghosts[i],
            d_variable_dst_ghosts[i],
            d_cart_grid_geometry,
            d_variable_refine_op[i]);
      } else if (d_do_coarsen) {
         commtest->registerVariable(d_variables_src[i],
            d_variables_dst[i],
            d_variable_src_ghosts[i],
            d_variable_dst_ghosts[i],
            d_cart_grid_geometry,
            d_variable_coarsen_op[i]);
      }

   }

}

void OutersideDataTest::initializeDataOnPatch(
   const hier::Patch& patch,
   const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number,
   char src_or_dst)
{
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   hier::VariableDatabase* variable_db =
      hier::VariableDatabase::getDatabase();
   variable_db->printClassData();
   tbox::Array<boost::shared_ptr<hier::Variable> >& variables(
      src_or_dst == 's' ? d_variables_src : d_variables_dst);

   if (d_do_refine) {

      for (int i = 0; i < variables.getSize(); i++) {

         boost::shared_ptr<hier::PatchData> data(
            patch.getPatchData(variables[i], getDataContext()));

         TBOX_ASSERT(data);

         boost::shared_ptr<pdat::OutersideData<double> > oside_data(
            data,
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<double> > side_data(
            data,
            boost::detail::dynamic_cast_tag());

         hier::Box dbox = data->getBox();

         if (side_data) {
            setLinearData(side_data, dbox, patch);
         }
         if (oside_data) {
            setLinearData(oside_data, dbox, patch);
         }

      }

   } else if (d_do_coarsen) {

      for (int i = 0; i < variables.getSize(); i++) {

         boost::shared_ptr<hier::PatchData> data(
            patch.getPatchData(variables[i], getDataContext()));

         TBOX_ASSERT(data);

         boost::shared_ptr<pdat::OutersideData<double> > oside_data(
            data,
            boost::detail::dynamic_cast_tag());
         boost::shared_ptr<pdat::SideData<double> > side_data(
            data,
            boost::detail::dynamic_cast_tag());

         hier::Box dbox = data->getGhostBox();

         if (side_data) {
            setLinearData(side_data, dbox, patch);
         }
         if (oside_data) {
            setLinearData(oside_data, dbox, patch);
         }
      }

   }

}

void OutersideDataTest::checkPatchInteriorData(
   const boost::shared_ptr<pdat::OutersideData<double> >& data,
   const hier::Box& interior,
   const boost::shared_ptr<geom::CartesianPatchGeometry>& pgeom) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data);
#endif

   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x = 0., y = 0., z = 0.;

   const int depth = data->getDepth();

   for (int axis = 0; axis < d_dim.getValue(); axis++) {
      const pdat::SideIndex loweri(interior.lower(), axis, 0);
      pdat::SideIterator siend(interior, axis, false);
      for (pdat::SideIterator si(interior, axis, true); si != siend; ++si) {

         /*
          * Compute spatial location of face and
          * set data to linear profile.
          */

         if (axis == 0) {
            x = lowerx[0] + dx[0] * ((*si)(0) - loweri(0));
         } else {
            x = lowerx[0] + dx[0] * ((*si)(0) - loweri(0) + 0.5);
         }
         y = z = 0.;
         if (d_dim > tbox::Dimension(1)) {
            if (axis == 1) {
               y = lowerx[1] + dx[1] * ((*si)(1) - loweri(1));
            } else {
               y = lowerx[1] + dx[1] * ((*si)(1) - loweri(1) + 0.5);
            }
         }
         if (d_dim > tbox::Dimension(2)) {
            if (axis == 2) {
               z = lowerx[2] + dx[2] * ((*si)(2) - loweri(2));
            } else {
               z = lowerx[2] + dx[2] * ((*si)(2) - loweri(2) + 0.5);
            }
         }

         double value;
         for (int d = 0; d < depth; d++) {
            value = d_Dcoef + d_Acoef * x + d_Bcoef * y + d_Ccoef * z;

            if (!(tbox::MathUtilities<double>::equalEps((*data)(*si,
                                                                d), value))) {
               tbox::perr << "FAILED: -- patch interior not properly filled"
                          << endl;
            }
         }
      }
   }
}

void OutersideDataTest::setPhysicalBoundaryConditions(
   const hier::Patch& patch,
   const double time,
   const hier::IntVector& gcw) const
{
   NULL_USE(patch);
   NULL_USE(gcw);
   NULL_USE(time);
   TBOX_ERROR("Only coarsen operations can be done with this test.\n"
      << "Coarsen operations should not need physical bc.\n");
}

void OutersideDataTest::setLinearData(
   boost::shared_ptr<pdat::SideData<double> > data,
   const hier::Box& box,
   const hier::Patch& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data);
#endif

   boost::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      patch.getPatchGeometry(),
      boost::detail::dynamic_cast_tag());
   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x, y, z;

   const int depth = data->getDepth();

   const hier::Box sbox = data->getGhostBox() * box;

   hier::IntVector directions(data->getDirectionVector());

   for (int axis = 0; axis < d_dim.getValue(); axis++) {
      if (directions(axis)) {
         const pdat::SideIndex loweri(patch.getBox().lower(), axis, 0);
         pdat::SideIterator eiend(sbox, axis, false);
         for (pdat::SideIterator ei(sbox, axis, true); ei != eiend; ++ei) {

            /*
             * Compute spatial location of cell center and
             * set data to linear profile.
             */

            if (axis == 0) {
               x = lowerx[0] + dx[0] * ((*ei)(0) - loweri(0));
            } else {
               x = lowerx[0] + dx[0] * ((*ei)(0) - loweri(0) + 0.5);
            }
            y = z = 0.;
            if (d_dim > tbox::Dimension(1)) {
               if (axis == 1) {
                  y = lowerx[1] + dx[1] * ((*ei)(1) - loweri(1));
               } else {
                  y = lowerx[1] + dx[1] * ((*ei)(1) - loweri(1) + 0.5);
               }
            }
            if (d_dim > tbox::Dimension(2)) {
               if (axis == 2) {
                  z = lowerx[2] + dx[2] * ((*ei)(2) - loweri(2));
               } else {
                  z = lowerx[2] + dx[2] * ((*ei)(2) - loweri(2) + 0.5);
               }
            }

            for (int d = 0; d < depth; d++) {
               (*data)(*ei,
                       d) = d_Dcoef + d_Acoef * x + d_Bcoef * y + d_Ccoef * z;
            }

         }
      }
   }
}

void OutersideDataTest::setLinearData(
   boost::shared_ptr<pdat::OutersideData<double> > data,
   const hier::Box& box,
   const hier::Patch& patch) const
{
   NULL_USE(box);

#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data);
#endif

   boost::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      patch.getPatchGeometry(),
      boost::detail::dynamic_cast_tag());
   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x = 0., y = 0., z = 0.;

   const int depth = data->getDepth();

   for (int axis = 0; axis < d_dim.getValue(); axis++) {
      for (int s = 0; s < 2; s++) {
         const hier::Box databox = data->getArrayData(axis, s).getBox();
         const pdat::SideIndex loweri(patch.getBox().lower(), axis, 0);
         hier::Box::iterator biend(databox, false);
         for (hier::Box::iterator bi(databox, true); bi != biend; ++bi) {

            /*
             * Compute spatial location of cell center and
             * set data to linear profile.
             */

            if (axis == 0) {
               x = lowerx[0] + dx[0] * ((*bi)(0) - loweri(0));
            } else {
               x = lowerx[0] + dx[0] * ((*bi)(0) - loweri(0) + 0.5);
            }
            y = z = 0.;
            if (d_dim > tbox::Dimension(1)) {
               if (axis == 1) {
                  y = lowerx[1] + dx[1] * ((*bi)(1) - loweri(1));
               } else {
                  y = lowerx[1] + dx[1] * ((*bi)(1) - loweri(1) + 0.5);
               }
            }
            if (d_dim > tbox::Dimension(2)) {
               if (axis == 2) {
                  z = lowerx[2] + dx[2] * ((*bi)(2) - loweri(2));
               } else {
                  z = lowerx[2] + dx[2] * ((*bi)(2) - loweri(2) + 0.5);
               }
            }

            double value = d_Dcoef + d_Acoef * x + d_Bcoef * y + d_Ccoef * z;
            pdat::SideIndex si(*bi, axis, 0);
            for (int d = 0; d < depth; d++) {
               (*data)(si, s, d) = value;
            }
         }

      }
   }
}

/*
 *************************************************************************
 *
 * Verify results of communication operations.  This test must be
 * consistent with data initialization and boundary operations above.
 *
 *************************************************************************
 */

bool OutersideDataTest::verifyResults(
   const hier::Patch& patch,
   const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number)
{
   NULL_USE(hierarchy);
   bool test_failed = false;
   if (d_do_refine || d_do_coarsen) {

      tbox::plog << "\nEntering OutersideDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl;

      hier::IntVector tgcw(d_dim, 0);
      for (int i = 0; i < d_variables_dst.getSize(); i++) {
         tgcw.max(patch.getPatchData(d_variables_dst[i], getDataContext())->
            getGhostCellWidth());
      }
      hier::Box pbox = patch.getBox();

      boost::shared_ptr<pdat::SideData<double> > solution(
         new pdat::SideData<double>(pbox, 1, tgcw));

      hier::Box tbox(pbox);
      tbox.grow(tgcw);

      if (d_do_refine) {
         setLinearData(solution, tbox, patch);
      } else {
         setLinearData(solution, tbox, patch); //, hierarchy, level_number);
      }

      for (int i = 0; i < d_variables_dst.getSize(); i++) {

         boost::shared_ptr<pdat::SideData<double> > side_data(
            patch.getPatchData(d_variables_dst[i], getDataContext()),
            boost::detail::dynamic_cast_tag());
         int depth = side_data->getDepth();
         hier::Box dbox = side_data->getGhostBox();

         hier::IntVector directions(side_data->getDirectionVector());

         for (int id = 0; id < d_dim.getValue(); id++) {
            if (directions(id)) {
               pdat::SideIterator siend(dbox, id, false);
               for (pdat::SideIterator si(dbox, id, true);
                    si != siend; ++si) {
                  double correct = (*solution)(*si);
                  for (int d = 0; d < depth; d++) {
                     double result = (*side_data)(*si, d);
                     if (!tbox::MathUtilities<double>::equalEps(correct,
                            result)) {
                        tbox::perr << "Test FAILED: ...."
                                   << " : side_data index = " << *si << endl;
                        tbox::perr << "    hier::Variable = "
                                   << d_variable_src_name[i]
                                   << " : depth index = " << d << endl;
                        tbox::perr << "    result = " << result
                                   << " : correct = " << correct << endl;
                        test_failed = true;
                     }
                  }
               }
            }
         }

      }
      if (!test_failed) {
         tbox::plog << "Outerside test Successful!" << endl;
      }

      solution.reset();   // just to be anal...

      tbox::plog << "\nExiting OutersidedataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl << endl;

   }

   return !test_failed;

}

}
