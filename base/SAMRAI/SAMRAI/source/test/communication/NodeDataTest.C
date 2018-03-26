/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   AMR communication tests for node-centered patch data
 *
 ************************************************************************/

#include "NodeDataTest.h"

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include "SAMRAI/pdat/NodeIndex.h"
#include "SAMRAI/pdat/NodeIterator.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "CommTester.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"

namespace SAMRAI {

using namespace std;

NodeDataTest::NodeDataTest(
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

   readTestInput(main_input_db->getDatabase("NodePatchDataTest"));

}

NodeDataTest::~NodeDataTest()
{
}

void NodeDataTest::readTestInput(
   boost::shared_ptr<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(db);
#endif

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
   readRefinementInput(db->getDatabase("RefinementData"));
}

void NodeDataTest::registerVariables(
   CommTester* commtest)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(commtest != (CommTester *)NULL);
#endif

   int nvars = d_variable_src_name.getSize();

   d_variables.resizeArray(nvars);

   for (int i = 0; i < nvars; i++) {
      d_variables[i].reset(
         new pdat::NodeVariable<double>(d_dim,
                                        d_variable_src_name[i],
                                        d_variable_depth[i]));

      if (d_do_refine) {
         commtest->registerVariable(d_variables[i],
            d_variables[i],
            d_variable_src_ghosts[i],
            d_variable_dst_ghosts[i],
            d_cart_grid_geometry,
            d_variable_refine_op[i]);
      } else if (d_do_coarsen) {
         commtest->registerVariable(d_variables[i],
            d_variables[i],
            d_variable_src_ghosts[i],
            d_variable_dst_ghosts[i],
            d_cart_grid_geometry,
            d_variable_coarsen_op[i]);
      }

   }

}

void NodeDataTest::setLinearData(
   boost::shared_ptr<pdat::NodeData<double> > data,
   const hier::Box& box,
   const hier::Patch& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data);
#endif

   boost::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      patch.getPatchGeometry(),
      boost::detail::dynamic_cast_tag());
   const pdat::NodeIndex loweri(
      patch.getBox().lower(), (pdat::NodeIndex::Corner)0);
   const double* dx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x, y, z;

   const int depth = data->getDepth();

   const hier::Box sbox = data->getGhostBox() * box;

   pdat::NodeIterator ciend(sbox, false);
   for (pdat::NodeIterator ci(sbox, true); ci != ciend; ++ci) {

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

      for (int d = 0; d < depth; d++) {
         (*data)(*ci, d) = d_Dcoef + d_Acoef * x + d_Bcoef * y + d_Ccoef * z;
      }

   }

}

void NodeDataTest::setPeriodicData(
   boost::shared_ptr<pdat::NodeData<double> > data,
   const hier::Box& box,
   const hier::Patch& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data);
#endif
   NULL_USE(patch);

   const double* xlo = d_cart_grid_geometry->getXLower();
   const double* xup = d_cart_grid_geometry->getXUpper();
   std::vector<double> domain_len(d_dim.getValue());
   for (int d = 0; d < d_dim.getValue(); ++d) {
      domain_len[d] = xup[d] - xlo[d];
   }

   const boost::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
      patch.getPatchGeometry(),
      boost::detail::dynamic_cast_tag());
   const double* dx = patch_geom->getDx();

   const int depth = data->getDepth();

   const hier::Box sbox = data->getGhostBox() * box;

   pdat::NodeIterator niend(sbox, false);
   for (pdat::NodeIterator ni(sbox, true); ni != niend; ++ni) {

      double val = 1.0;
      for (int d = 0; d < d_dim.getValue(); ++d) {
         double tmpf = dx[d] * (*ni)(d) / domain_len[d];
         tmpf = sin(2 * M_PI * tmpf);
         val *= tmpf;
      }
      val = val + 20.0; // Shift function range to [1,3] to avoid bad floating point compares.
      for (int d = 0; d < depth; d++) {
         (*data)(*ni, d) = val;
      }

   }

}

void NodeDataTest::initializeDataOnPatch(
   const hier::Patch& patch,
   const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number,
   char src_or_dst)
{
   NULL_USE(src_or_dst);
   NULL_USE(level_number);
   NULL_USE(hierarchy);

   const hier::IntVector periodic_shift(
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(d_dim, 1)));
   bool is_periodic = periodic_shift.max() > 0;

   if (d_do_refine) {

      for (int i = 0; i < d_variables.getSize(); i++) {

         boost::shared_ptr<pdat::NodeData<double> > node_data(
            patch.getPatchData(d_variables[i], getDataContext()),
            boost::detail::dynamic_cast_tag());

         hier::Box dbox = node_data->getBox();

         if (is_periodic) {
            setPeriodicData(node_data, dbox, patch);
         } else {
            setLinearData(node_data, dbox, patch);
         }

      }

   } else if (d_do_coarsen) {

      for (int i = 0; i < d_variables.getSize(); i++) {

         boost::shared_ptr<pdat::NodeData<double> > node_data(
            patch.getPatchData(d_variables[i], getDataContext()),
            boost::detail::dynamic_cast_tag());

         hier::Box dbox = node_data->getGhostBox();

         if (is_periodic) {
            setPeriodicData(node_data, dbox, patch);
         } else {
            setLinearData(node_data, dbox, patch);
         }

      }

   }

}

void NodeDataTest::checkPatchInteriorData(
   const boost::shared_ptr<pdat::NodeData<double> >& data,
   const hier::Box& interior,
   const hier::Patch& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data);
#endif

   const bool is_periodic =
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(d_dim,
            1)).max() > 0;

   const int depth = data->getDepth();

   boost::shared_ptr<pdat::NodeData<double> > correct_data(
      new pdat::NodeData<double>(
         data->getBox(),
         depth,
         data->getGhostCellWidth()));
   if (is_periodic) {
      setPeriodicData(correct_data, correct_data->getGhostBox(), patch);
   } else {
      setLinearData(correct_data, correct_data->getGhostBox(), patch);
   }

   pdat::NodeIterator niend(interior, false);
   for (pdat::NodeIterator ni(interior, true); ni != niend; ++ni) {
      for (int d = 0; d < depth; d++) {
         if (!(tbox::MathUtilities<double>::equalEps((*data)(*ni, d),
                  (*correct_data)(*ni, d)))) {
            tbox::perr << "FAILED: -- patch interior not properly filled"
                       << endl;
         }
      }
   }

}

void NodeDataTest::setPhysicalBoundaryConditions(
   const hier::Patch& patch,
   const double time,
   const hier::IntVector& gcw) const
{
   NULL_USE(time);

   const hier::IntVector periodic_shift =
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(d_dim, 1));
   bool is_periodic = periodic_shift.max() > 0;

   boost::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      patch.getPatchGeometry(),
      boost::detail::dynamic_cast_tag());

   const tbox::Array<hier::BoundaryBox> node_bdry =
      pgeom->getCodimensionBoundaries(d_dim.getValue());
   const int num_node_bdry_boxes = node_bdry.getSize();

   tbox::Array<hier::BoundaryBox> edge_bdry;
   if (d_dim > tbox::Dimension(1)) {
      edge_bdry = pgeom->getCodimensionBoundaries(d_dim.getValue() - 1);
   }
   const int num_edge_bdry_boxes = d_dim > tbox::Dimension(1) ? edge_bdry.getSize() : -1;

   tbox::Array<hier::BoundaryBox> face_bdry;
   if (d_dim == tbox::Dimension(3)) {
      face_bdry = pgeom->getCodimensionBoundaries(d_dim.getValue() - 2);
   }
   const int num_face_bdry_boxes = d_dim == tbox::Dimension(3) ? face_bdry.getSize() : -1;

   for (int i = 0; i < d_variables.getSize(); i++) {

      boost::shared_ptr<pdat::NodeData<double> > node_data(
         patch.getPatchData(d_variables[i], getDataContext()),
         boost::detail::dynamic_cast_tag());

      hier::Box patch_interior = node_data->getBox();
      checkPatchInteriorData(node_data, patch_interior, patch);

      /*
       * Set node boundary data.
       */
      for (int ni = 0; ni < num_node_bdry_boxes; ni++) {

         hier::Box fill_box = pgeom->getBoundaryFillBox(node_bdry[ni],
               patch.getBox(),
               gcw);

         if (is_periodic) {
            setPeriodicData(node_data, fill_box, patch);
         } else {
            setLinearData(node_data, fill_box, patch);
         }
      }

      if (d_dim > tbox::Dimension(1)) {
         /*
          * Set edge boundary data.
          */
         for (int ei = 0; ei < num_edge_bdry_boxes; ei++) {

            hier::Box fill_box = pgeom->getBoundaryFillBox(edge_bdry[ei],
                  patch.getBox(),
                  gcw);

            if (is_periodic) {
               setPeriodicData(node_data, fill_box, patch);
            } else {
               setLinearData(node_data, fill_box, patch);
            }
         }
      }

      if (d_dim == tbox::Dimension(3)) {
         /*
          * Set face boundary data.
          */
         for (int fi = 0; fi < num_face_bdry_boxes; fi++) {

            hier::Box fill_box = pgeom->getBoundaryFillBox(face_bdry[fi],
                  patch.getBox(),
                  gcw);

            if (is_periodic) {
               setPeriodicData(node_data, fill_box, patch);
            } else {
               setLinearData(node_data, fill_box, patch);
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
bool NodeDataTest::verifyResults(
   const hier::Patch& patch,
   const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number)
{
   NULL_USE(hierarchy);

   bool test_failed = false;

   const hier::IntVector periodic_shift(
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(d_dim, 1)));
   bool is_periodic = periodic_shift.max() > 0;

   if (d_do_refine || d_do_coarsen) {

      tbox::plog << "\nEntering NodeDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl;

      hier::IntVector tgcw(periodic_shift.getDim(), 0);
      for (int i = 0; i < d_variables.getSize(); i++) {
         tgcw.max(patch.getPatchData(d_variables[i], getDataContext())->
            getGhostCellWidth());
      }
      hier::Box pbox = patch.getBox();

      boost::shared_ptr<pdat::NodeData<double> > solution(
         new pdat::NodeData<double>(pbox, 1, tgcw));

      hier::Box gbox(pbox);
      gbox.grow(tgcw);

      if (d_do_refine) {
         if (is_periodic) {
            setPeriodicData(solution, gbox, patch);
         } else {
            setLinearData(solution, gbox, patch);
         }
      } else {
         if (is_periodic) {
            setPeriodicData(solution, gbox, patch);
         } else {
            setLinearData(solution, gbox, patch);
         }
      }

      for (int i = 0; i < d_variables.getSize(); i++) {

         boost::shared_ptr<pdat::NodeData<double> > node_data(
            patch.getPatchData(d_variables[i], getDataContext()),
            boost::detail::dynamic_cast_tag());
         int depth = node_data->getDepth();
         hier::Box dbox = node_data->getGhostBox();

         pdat::NodeIterator ciend(dbox, false);
         for (pdat::NodeIterator ci(dbox, true); ci != ciend; ++ci) {
            double correct = (*solution)(*ci);
            for (int d = 0; d < depth; d++) {
               double result = (*node_data)(*ci, d);
               if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                  tbox::perr << "Test FAILED: ...."
                             << " : node index = " << *ci
                             << " of L" << level_number
                             << " P" << patch.getLocalId()
                             << " " << patch.getBox() << endl;
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
         tbox::plog << "Node test Successful!" << endl;
      }

      solution.reset();   // just to be anal...

      tbox::plog << "\nExiting NodeDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl << endl;

   }

   return !test_failed;

}

}
