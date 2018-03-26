/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   AMR communication tests for face-centered patch data
 *
 ************************************************************************/

#include "FaceDataTest.h"

#include "SAMRAI/pdat/ArrayData.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "CommTester.h"
#include "SAMRAI/pdat/FaceGeometry.h"
#include "SAMRAI/pdat/FaceIndex.h"
#include "SAMRAI/pdat/FaceIterator.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

namespace SAMRAI {

using namespace std;

FaceDataTest::FaceDataTest(
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

   readTestInput(main_input_db->getDatabase("FacePatchDataTest"));

}

FaceDataTest::~FaceDataTest()
{
}

void FaceDataTest::readTestInput(
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

void FaceDataTest::registerVariables(
   CommTester* commtest)
{
   TBOX_ASSERT(commtest != (CommTester *)NULL);

   int nvars = d_variable_src_name.getSize();

   d_variables.resizeArray(nvars);

   for (int i = 0; i < nvars; i++) {
      d_variables[i].reset(
         new pdat::FaceVariable<double>(d_dim,
                                        d_variable_src_name[i],
                                        d_variable_depth[i],
                                        d_use_fine_value_at_interface[i]));

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

void FaceDataTest::setConservativeData(
   boost::shared_ptr<pdat::FaceData<double> > data,
   const hier::Box& box,
   const hier::Patch& patch,
   const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data);
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= hierarchy->getFinestLevelNumber()));
#endif

   int i, j;
   boost::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(level_number));

   const hier::BoxContainer& domain =
      level->getPhysicalDomain(hier::BlockId::zero());
   int ncells = 0;
   for (hier::BoxContainer::const_iterator i(domain); i != domain.end(); ++i) {
      ncells += i->size();
   }

   const int depth = data->getDepth();

   const hier::Box sbox = data->getGhostBox() * box;

   if (level_number == 0) {

      /*
       * Set face value on level zero as follows:
       *
       *    u0(i,j,k) = (j + k)/ncells
       *    u1(i,j,k) = (i + k)/ncells
       *    u2(i,j,k) = (i + j)/ncells
       */

      for (int axis = 0; axis < d_dim.getValue(); axis++) {
         pdat::CellIterator ciend(sbox, false);
         for (pdat::CellIterator ci(sbox, true); ci != ciend; ++ci) {
            double value = 0.0;
            for (i = 0; i < d_dim.getValue(); i++) {
               if (i != axis) {
                  value += (double)((*ci)(i));
               }
            }
            value /= ncells;
            for (int face = pdat::FaceIndex::Lower;
                 face <= pdat::FaceIndex::Upper; face++) {
               pdat::FaceIndex si(*ci, axis, face);
               for (int d = 0; d < depth; d++) {
                  (*data)(si, d) = value;
               }
            }
         }
      }

   } else {

      /*
       * Set face value on level > 0 to
       *    u(i,j,k) = u_c + ci*del_i + cj*del_j + ck*del_k
       * where u_c is value on the underlying coarse face, (ci,cj,ck) is
       * the underlying coarse face index, and (del_i,del_j,del_k)
       * is the vector between the coarse and fine cell face centers.
       */

      hier::IntVector ratio(level->getRatioToLevelZero());
      const int max_ratio = ratio.max();

      boost::shared_ptr<geom::CartesianPatchGeometry> pgeom(
         patch.getPatchGeometry(),
         boost::detail::dynamic_cast_tag());
      const double* dx = pgeom->getDx();

      int coarse_ncells = ncells;
      double* delta = new double[max_ratio * d_dim.getValue()];
      for (j = 0; j < d_dim.getValue(); j++) {
         coarse_ncells /= ratio(j);
         double coarse_dx = dx[j] * ratio(j);
         for (i = 0; i < ratio(j); i++) {
            delta[j * max_ratio + i] = (i + 0.5) * dx[j] - coarse_dx * 0.5;
         }
      }

      for (int axis = 0; axis < d_dim.getValue(); axis++) {
         hier::IntVector ci(ratio.getDim());
         hier::IntVector del(ratio.getDim());
         pdat::CellIterator fiend(sbox, false);
         for (pdat::CellIterator fi(sbox, true); fi != fiend; ++fi) {
            double value = 0.0;
            for (i = 0; i < d_dim.getValue(); i++) {
               if (i != axis) {
                  int findx = (*fi)(i);
                  ci(i) = ((findx < 0) ? (findx + 1) / ratio(i) - 1
                           : findx / ratio(i));
                  del(i) = (int)delta[i * max_ratio + findx - ci(i) * ratio(i)];
                  value += (double)(ci(i));
               }
            }
            value /= coarse_ncells;

            for (j = 0; j < d_dim.getValue(); j++) {
               if (j != axis) {
                  value += ci(j) * del(j);
               }
            }

            for (int face = pdat::FaceIndex::Lower;
                 face <= pdat::FaceIndex::Upper; face++) {
               pdat::FaceIndex si(*fi, axis, face);
               for (int d = 0; d < depth; d++) {
                  (*data)(si, d) = value;
               }
            }
         }
      }
      delete[] delta;

   }

}

void FaceDataTest::initializeDataOnPatch(
   const hier::Patch& patch,
   const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number,
   char src_or_dst)
{
   NULL_USE(src_or_dst);

   if (d_do_refine) {

      for (int i = 0; i < d_variables.getSize(); i++) {

         boost::shared_ptr<pdat::FaceData<double> > face_data(
            patch.getPatchData(d_variables[i], getDataContext()),
            boost::detail::dynamic_cast_tag());

         hier::Box dbox = face_data->getBox();

         setLinearData(face_data, dbox, patch);
      }

   } else if (d_do_coarsen) {

      for (int i = 0; i < d_variables.getSize(); i++) {

         boost::shared_ptr<pdat::FaceData<double> > face_data(
            patch.getPatchData(d_variables[i], getDataContext()),
            boost::detail::dynamic_cast_tag());

         hier::Box dbox = face_data->getGhostBox();

         setConservativeData(face_data, dbox,
            patch, hierarchy, level_number);

      }

   }

}

void FaceDataTest::checkPatchInteriorData(
   const boost::shared_ptr<pdat::FaceData<double> >& data,
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
      const pdat::FaceIndex loweri(interior.lower(), axis, 0);
      pdat::FaceIterator fiend(interior, axis, false);
      for (pdat::FaceIterator fi(interior, axis, true); fi != fiend; ++fi) {

         /*
          * Compute spatial location of face and
          * set data to linear profile.
          */

         if (axis == 0) {
            x = lowerx[0] + dx[0] * ((*fi)(0) - loweri(0));
            if (d_dim > tbox::Dimension(1)) {
               y = lowerx[1] + dx[1] * ((*fi)(1) - loweri(1) + 0.5);
            }
            if (d_dim > tbox::Dimension(2)) {
               z = lowerx[2] + dx[2] * ((*fi)(2) - loweri(2) + 0.5);
            }
         } else if (axis == 1) {
            x = lowerx[0] + dx[0]
               * ((*fi)(d_dim.getValue() - 1) - loweri(d_dim.getValue() - 1) + 0.5);
            if (d_dim > tbox::Dimension(1)) {
               y = lowerx[1] + dx[1] * ((*fi)(0) - loweri(0));
            }
            if (d_dim > tbox::Dimension(2)) {
               z = lowerx[2] + dx[2] * ((*fi)(1) - loweri(1) + 0.5);
            }
         } else if (axis == 2) {
            x = lowerx[0] + dx[0] * ((*fi)(1) - loweri(1) + 0.5);
            if (d_dim > tbox::Dimension(1)) {
               y = lowerx[1] + dx[1] * ((*fi)(2) - loweri(2) + 0.5);
            }
            if (d_dim > tbox::Dimension(2)) {
               z = lowerx[2] + dx[2] * ((*fi)(0) - loweri(0));
            }
         }

         double value;
         for (int d = 0; d < depth; d++) {
            value = d_Dcoef + d_Acoef * x + d_Bcoef * y + d_Ccoef * z;
            if (!(tbox::MathUtilities<double>::equalEps((*data)(*fi,
                                                                d), value))) {
               tbox::perr << "FAILED: -- patch interior not properly filled"
                          << endl;
            }
         }
      }
   }
}

void FaceDataTest::setPhysicalBoundaryConditions(
   const hier::Patch& patch,
   const double time,
   const hier::IntVector& gcw) const
{
   NULL_USE(time);

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
   const int num_edge_bdry_boxes = d_dim.getValue() > 1 ? edge_bdry.getSize() : -1;

   tbox::Array<hier::BoundaryBox> face_bdry;
   if (d_dim == tbox::Dimension(3)) {
      face_bdry = pgeom->getCodimensionBoundaries(d_dim.getValue() - 2);
   }
   const int num_face_bdry_boxes = d_dim == tbox::Dimension(3) ? face_bdry.getSize() : -1;

   for (int i = 0; i < d_variables.getSize(); i++) {

      boost::shared_ptr<pdat::FaceData<double> > face_data(
         patch.getPatchData(d_variables[i], getDataContext()),
         boost::detail::dynamic_cast_tag());

      hier::Box patch_interior = face_data->getBox();
      checkPatchInteriorData(face_data, patch_interior, pgeom);

      /*
       * Set node boundary data.
       */
      for (int ni = 0; ni < num_node_bdry_boxes; ni++) {

         hier::Box fill_box = pgeom->getBoundaryFillBox(node_bdry[ni],
               patch.getBox(),
               gcw);

         setLinearData(face_data, fill_box, patch);
      }

      if (d_dim > tbox::Dimension(1)) {
         /*
          * Set edge boundary data.
          */
         for (int ei = 0; ei < num_edge_bdry_boxes; ei++) {

            hier::Box fill_box = pgeom->getBoundaryFillBox(edge_bdry[ei],
                  patch.getBox(),
                  gcw);

            setLinearData(face_data, fill_box, patch);
         }
      }

      if (d_dim > tbox::Dimension(2)) {
         /*
          * Set face boundary data.
          */
         for (int fi = 0; fi < num_face_bdry_boxes; fi++) {

            hier::Box fill_box = pgeom->getBoundaryFillBox(face_bdry[fi],
                  patch.getBox(),
                  gcw);

            setLinearData(face_data, fill_box, patch);
         }
      }

   }

}

void FaceDataTest::setLinearData(
   boost::shared_ptr<pdat::FaceData<double> > data,
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
   double x = 0., y = 0., z = 0.;

   const int depth = data->getDepth();

   const hier::Box sbox = data->getGhostBox() * box;

   for (int axis = 0; axis < d_dim.getValue(); axis++) {
      const pdat::FaceIndex loweri(patch.getBox().lower(), axis, 0);
      pdat::FaceIterator fiend(sbox, axis, false);
      for (pdat::FaceIterator fi(sbox, axis, true); fi != fiend; ++fi) {

         /*
          * Compute spatial location of cell center and
          * set data to linear profile.
          */

         if (axis == 0) {
            x = lowerx[0] + dx[0] * ((*fi)(0) - loweri(0));
            if (d_dim > tbox::Dimension(1)) {
               y = lowerx[1] + dx[1] * ((*fi)(1) - loweri(1) + 0.5);
            }
            if (d_dim > tbox::Dimension(2)) {
               z = lowerx[2] + dx[2] * ((*fi)(2) - loweri(2) + 0.5);
            }
         } else if (axis == 1) {
            x = lowerx[0] + dx[0]
               * ((*fi)(d_dim.getValue() - 1) - loweri(d_dim.getValue() - 1) + 0.5);
            if (d_dim > tbox::Dimension(1)) {
               y = lowerx[1] + dx[1] * ((*fi)(0) - loweri(0));
            }
            if (d_dim > tbox::Dimension(2)) {
               z = lowerx[2] + dx[2] * ((*fi)(1) - loweri(1) + 0.5);
            }
         } else if (axis == 2) {
            x = lowerx[0] + dx[0] * ((*fi)(1) - loweri(1) + 0.5);
            if (d_dim > tbox::Dimension(1)) {
               y = lowerx[1] + dx[1] * ((*fi)(2) - loweri(2) + 0.5);
            }
            if (d_dim > tbox::Dimension(2)) {
               z = lowerx[2] + dx[2] * ((*fi)(0) - loweri(0));
            }
         }

         for (int d = 0; d < depth; d++) {
            (*data)(*fi,
                    d) = d_Dcoef + d_Acoef * x + d_Bcoef * y + d_Ccoef * z;
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

bool FaceDataTest::verifyResults(
   const hier::Patch& patch,
   const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number)
{
   bool test_failed = false;
   if (d_do_refine || d_do_coarsen) {

      tbox::plog << "\nEntering FaceDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl;

      hier::IntVector tgcw(d_dim, 0);
      for (int i = 0; i < d_variables.getSize(); i++) {
         tgcw.max(patch.getPatchData(d_variables[i], getDataContext())->
            getGhostCellWidth());
      }
      hier::Box pbox = patch.getBox();

      boost::shared_ptr<pdat::FaceData<double> > solution(
         new pdat::FaceData<double>(pbox, 1, tgcw));

      hier::Box tbox(pbox);
      tbox.grow(tgcw);

      if (d_do_coarsen) {
         setConservativeData(solution, tbox,
            patch, hierarchy, level_number);
      }

      for (int i = 0; i < d_variables.getSize(); i++) {

         boost::shared_ptr<pdat::FaceData<double> > face_data(
            patch.getPatchData(d_variables[i], getDataContext()),
            boost::detail::dynamic_cast_tag());
         int depth = face_data->getDepth();
         hier::Box dbox = face_data->getGhostBox();

         if (d_do_refine) {
            setLinearData(solution, tbox, patch);
         }

         for (int id = 0; id < d_dim.getValue(); id++) {
            pdat::FaceIterator siend(dbox, id, false);
            for (pdat::FaceIterator si(dbox, id, true); si != siend; ++si) {
               double correct = (*solution)(*si);
               for (int d = 0; d < depth; d++) {
                  double result = (*face_data)(*si, d);
                  if (!tbox::MathUtilities<double>::equalEps(correct,
                         result)) {
                     test_failed = true;
                     tbox::perr << "Test FAILED: ...."
                                << " : face_data index = " << *si << endl;
                     tbox::perr << "    hier::Variable = "
                                << d_variable_src_name[i]
                                << " : depth index = " << d << endl;
                     tbox::perr << "    result = " << result
                                << " : correct = " << correct << endl;
                  }
               }
            }
         }

      }

      solution.reset();   // just to be anal...

      tbox::plog << "\nExiting FaceDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl << endl;

   }

   return !test_failed;
}

}
