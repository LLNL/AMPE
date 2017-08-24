/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2016 Lawrence Livermore National Security, LLC
 * Description:   AMR communication tests for node-centered patch data
 *
 ************************************************************************/

#include "OuternodeDataTest.h"

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/NodeData.h"
#include "SAMRAI/pdat/NodeIndex.h"
#include "SAMRAI/pdat/NodeIterator.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/OuternodeVariable.h"
#include "SAMRAI/hier/PatchData.h"
#include "CommTester.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"

namespace SAMRAI {

using namespace std;

OuternodeDataTest::OuternodeDataTest(
   const string& object_name,
   const tbox::Dimension& dim,
   boost::shared_ptr<tbox::Database> main_input_db,
   bool do_refine,
   bool do_coarsen,
   const string& refine_option):
   PatchDataTestStrategy(dim),
   d_dim(dim)
{
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(main_input_db);
   TBOX_ASSERT(!refine_option.empty());

   d_object_name = object_name;

   d_do_refine = do_refine;
   d_do_coarsen = false;
   if (!do_refine) {
      d_do_coarsen = do_coarsen;
   }
   if (d_do_refine) {
      TBOX_ERROR("There is no refine test for Outernode data type, because\n"
         << "Outernode refinement does not exist at this time.");
      /*
       * The refine codes are still kept in this class in case we
       * somehow define Outernode refinement in the future.
       */
   }

   d_refine_option = refine_option;

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

   readTestInput(main_input_db->getDatabase("OuternodePatchDataTest"));

}

OuternodeDataTest::~OuternodeDataTest()
{
}

void OuternodeDataTest::readTestInput(
   boost::shared_ptr<tbox::Database> db)
{
   TBOX_ASSERT(db);

   /*
    * Read coeeficients of linear profile to test interpolation.
    */
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

   /*
    * Base class reads variable parameters and boxes to refine.
    */

   readVariableInput(db->getDatabase("VariableData"));
}

void OuternodeDataTest::registerVariables(
   CommTester* commtest)
{
   TBOX_ASSERT(commtest != 0);

   int nvars = static_cast<int>(d_variable_src_name.size());

   d_variables_src.resize(nvars);
   d_variables_dst.resize(nvars);

   for (int i = 0; i < nvars; ++i) {
      d_variables_src[i].reset(
         new pdat::OuternodeVariable<double>(
            d_dim,
            d_variable_src_name[i],
            d_variable_depth[i]));
      d_variables_dst[i].reset(
         new pdat::NodeVariable<double>(
            d_dim,
            d_variable_dst_name[i],
            d_variable_depth[i]));

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

void OuternodeDataTest::setLinearData(
   boost::shared_ptr<pdat::NodeData<double> > data,
   const hier::Box& box,
   const hier::Patch& patch) const
{
   TBOX_ASSERT(data);

   boost::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(pgeom);
   const pdat::NodeIndex loweri(
      patch.getBox().lower(), (pdat::NodeIndex::Corner)0);
   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x, y, z;

   const int depth = data->getDepth();

   const hier::Box sbox = data->getGhostBox() * box;

   pdat::NodeIterator ciend(pdat::NodeGeometry::end(sbox));
   for (pdat::NodeIterator ci(pdat::NodeGeometry::begin(sbox));
        ci != ciend; ++ci) {

      /*
       * Compute spatial location of node center and
       * set data to linear profile.
       */

      x = lowerx[0] + dx[0] * ((*ci)(0) - loweri(0));
      y = z = 0.;
      if (d_dim > tbox::Dimension(1)) {
         y = lowerx[1] + dx[1] * ((*ci)(1) - loweri(1));
      }
      if (d_dim > tbox::Dimension(2)) {
         z = lowerx[2] + dx[2] * ((*ci)(2) - loweri(2));
      }

      for (int d = 0; d < depth; ++d) {
         (*data)(*ci, d) = d_Dcoef + d_Acoef * x + d_Bcoef * y + d_Ccoef * z;
      }

   }

}

void OuternodeDataTest::setLinearData(
   boost::shared_ptr<pdat::OuternodeData<double> > data,
   const hier::Box& box,
   const hier::Patch& patch) const
{
   NULL_USE(box);

   TBOX_ASSERT(data);
   TBOX_ASSERT(box.isSpatiallyEqual(patch.getBox()));
#ifdef DEBUG_CHECK_ASSERTIONS
   if (!box.isSpatiallyEqual(data->getBox())) {
      TBOX_ERROR("Box is not identical to data box, which is\n"
         << "required for testing Outernode communication.");
   }
#endif

   boost::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      BOOST_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
         patch.getPatchGeometry()));
   TBOX_ASSERT(pgeom);
   const pdat::NodeIndex loweri(
      patch.getBox().lower(), (pdat::NodeIndex::Corner)0);
   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x, y, z;

   const int depth = data->getDepth();

   for (tbox::Dimension::dir_t n = 0; n < d_dim.getValue(); ++n) {
      for (int s = 0; s < 2; ++s) {
         const hier::Box databox = data->getDataBox(n, s);
         hier::Box::iterator biend(databox.end());
         for (hier::Box::iterator bi(databox.begin()); bi != biend; ++bi) {

            /*
             * Compute spatial location of node center and
             * set data to linear profile.
             */

            x = lowerx[0] + dx[0] * ((*bi)(0) - loweri(0));
            y = z = 0.;
            if (d_dim > tbox::Dimension(1)) {
               y = lowerx[1] + dx[1] * ((*bi)(1) - loweri(1));
            }
            if (d_dim > tbox::Dimension(2)) {
               z = lowerx[2] + dx[2] * ((*bi)(2) - loweri(2));
            }

            pdat::NodeIndex ni(*bi, (pdat::NodeIndex::Corner)0);
            for (int d = 0; d < depth; ++d) {
               (*data)(ni,
                       d) = d_Dcoef + d_Acoef * x + d_Bcoef * y + d_Ccoef * z;
            }
         }
      }

   }

}

void OuternodeDataTest::initializeDataOnPatch(
   const hier::Patch& patch,
   const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number,
   char src_or_dst)
{
   NULL_USE(hierarchy);
   NULL_USE(level_number);
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();
   variable_db->printClassData();
   std::vector<boost::shared_ptr<hier::Variable> >& variables(
      src_or_dst == 's' ? d_variables_src : d_variables_dst);

   if (d_do_refine) {

      for (int i = 0; i < static_cast<int>(variables.size()); ++i) {

         boost::shared_ptr<hier::PatchData> data(
            patch.getPatchData(variables[i], getDataContext()));
         TBOX_ASSERT(data);

         boost::shared_ptr<pdat::OuternodeData<double> > onode_data(
            boost::dynamic_pointer_cast<pdat::OuternodeData<double>,
                                        hier::PatchData>(data));
         boost::shared_ptr<pdat::NodeData<double> > node_data(
            boost::dynamic_pointer_cast<pdat::NodeData<double>,
                                        hier::PatchData>(data));

         hier::Box dbox = data->getBox();

         if (node_data) {
            setLinearData(node_data, dbox, patch);
         }
         if (onode_data) {
            setLinearData(onode_data, dbox, patch);
         }

      }

   } else if (d_do_coarsen) {

      for (int i = 0; i < static_cast<int>(variables.size()); ++i) {

         boost::shared_ptr<hier::PatchData> data(
            patch.getPatchData(variables[i], getDataContext()));
         TBOX_ASSERT(data);
         boost::shared_ptr<pdat::OuternodeData<double> > onode_data(
            boost::dynamic_pointer_cast<pdat::OuternodeData<double>,
                                        hier::PatchData>(data));
         boost::shared_ptr<pdat::NodeData<double> > node_data(
            boost::dynamic_pointer_cast<pdat::NodeData<double>,
                                        hier::PatchData>(data));

         hier::Box dbox = data->getGhostBox();

         if (node_data) {
            setLinearData(node_data, dbox, patch);
         }
         if (onode_data) {
            setLinearData(onode_data, dbox, patch);
         }

      }

   }

}

void OuternodeDataTest::checkPatchInteriorData(
   const boost::shared_ptr<pdat::OuternodeData<double> >& data,
   const hier::Box& interior,
   const boost::shared_ptr<geom::CartesianPatchGeometry>& pgeom) const
{
   TBOX_ASSERT(data);

   const pdat::NodeIndex loweri(interior.lower(), (pdat::NodeIndex::Corner)0);
   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x, y, z;

   const int depth = data->getDepth();

   pdat::NodeIterator ciend(pdat::NodeGeometry::end(interior));
   for (pdat::NodeIterator ci(pdat::NodeGeometry::begin(interior));
        ci != ciend; ++ci) {

      /*
       * Compute spatial location of edge and
       * compare data to linear profile.
       */

      x = lowerx[0] + dx[0] * ((*ci)(0) - loweri(0));
      y = z = 0.;
      if (d_dim > tbox::Dimension(1)) {
         y = lowerx[1] + dx[1] * ((*ci)(1) - loweri(1));
      }
      if (d_dim > tbox::Dimension(2)) {
         z = lowerx[2] + dx[2] * ((*ci)(2) - loweri(2));
      }

      double value;
      for (int d = 0; d < depth; ++d) {
         value = d_Dcoef + d_Acoef * x + d_Bcoef * y + d_Ccoef * z;
         if (!(tbox::MathUtilities<double>::equalEps((*data)(*ci,
                                                             d), value))) {
            tbox::perr << "FAILED: -- patch interior not properly filled"
                       << endl;
         }
      }

   }

}

void OuternodeDataTest::setPhysicalBoundaryConditions(
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

/*
 *************************************************************************
 *
 * Verify results of communication operations.  This test must be
 * consistent with data initialization and boundary operations above.
 *
 *************************************************************************
 */
bool OuternodeDataTest::verifyResults(
   const hier::Patch& patch,
   const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number)
{
   NULL_USE(hierarchy);
   bool test_failed = false;
   if (d_do_refine || d_do_coarsen) {

      tbox::plog << "\nEntering OuternodeDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl;

      hier::IntVector tgcw(d_dim, 0);
      for (int i = 0; i < static_cast<int>(d_variables_dst.size()); ++i) {
         tgcw.max(patch.getPatchData(d_variables_dst[i], getDataContext())->
            getGhostCellWidth());
      }
      hier::Box pbox = patch.getBox();

      boost::shared_ptr<pdat::NodeData<double> > solution(
         new pdat::NodeData<double>(pbox, 1, tgcw));

      hier::Box tbox(pbox);
      tbox.grow(tgcw);

      if (d_do_refine) {
         setLinearData(solution, tbox, patch);
      } else {
         setLinearData(solution, tbox,
            patch);                 //, hierarchy, level_number);
      }

      for (int i = 0; i < static_cast<int>(d_variables_dst.size()); ++i) {

         boost::shared_ptr<pdat::NodeData<double> > node_data(
            BOOST_CAST<pdat::NodeData<double>, hier::PatchData>(
               patch.getPatchData(d_variables_dst[i], getDataContext())));
         TBOX_ASSERT(node_data);
         int depth = node_data->getDepth();
         hier::Box dbox = node_data->getGhostBox();

         pdat::NodeIterator ciend(pdat::NodeGeometry::end(dbox));
         for (pdat::NodeIterator ci(pdat::NodeGeometry::begin(dbox));
              ci != ciend; ++ci) {
            double correct = (*solution)(*ci);
            for (int d = 0; d < depth; ++d) {
               double result = (*node_data)(*ci, d);
               if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                  tbox::perr << "Test FAILED: ...."
                             << " : node index = " << *ci << endl;
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
      if (!test_failed) {
         tbox::plog << "Outernode test Successful!" << endl;
      }

      solution.reset();   // just to be anal...

      tbox::plog << "\nExiting OuternodeDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl << endl;

   }

   return !test_failed;
}

}
