/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   AMR communication tests for cell-centered patch data
 *
 ************************************************************************/

#include "CellDataTest.h"

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellIndex.h"
#include "SAMRAI/pdat/CellIterator.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "CommTester.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Database.h"

namespace SAMRAI {

using namespace std;

CellDataTest::CellDataTest(
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
      new geom::CartesianGridGeometry(dim,
         "CartesianGridGeometry",
         main_input_db->getDatabase("CartesianGridGeometry")));

   setGridGeometry(d_cart_grid_geometry);

   readTestInput(main_input_db->getDatabase("CellPatchDataTest"));

}

CellDataTest::~CellDataTest()
{
}

void CellDataTest::readTestInput(
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
   readRefinementInput(db->getDatabase("RefinementData"));
}

void CellDataTest::registerVariables(
   CommTester* commtest)
{
   TBOX_ASSERT(commtest != (CommTester *)NULL);

   int nvars = d_variable_src_name.getSize();

   d_variables.resizeArray(nvars);

   for (int i = 0; i < nvars; i++) {
      d_variables[i].reset(
         new pdat::CellVariable<double>(d_dim,
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

void CellDataTest::setLinearData(
   boost::shared_ptr<pdat::CellData<double> > data,
   const hier::Box& box,
   const hier::Patch& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(data);
#endif

   boost::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      patch.getPatchGeometry(),
      boost::detail::dynamic_cast_tag());
   const pdat::CellIndex loweri(patch.getBox().lower());
   const pdat::CellIndex upperi(patch.getBox().upper());
   const double* pdx = pgeom->getDx();
   const double* lowerx = pgeom->getXLower();
   double x, y, z;

   const int depth = data->getDepth();

   const hier::Box sbox = data->getGhostBox() * box;

   pdat::CellIterator ciend(sbox, false);
   for (pdat::CellIterator ci(sbox, true); ci != ciend; ++ci) {

      /*
       * Compute spatial location of cell center and
       * set data to linear profile.
       */

      x = lowerx[0] + pdx[0] * ((*ci)(0) - loweri(0) + 0.5);
      y = z = 0.;
      if (d_dim > tbox::Dimension(1)) {
         y = lowerx[1] + pdx[1] * ((*ci)(1) - loweri(1) + 0.5);
      }
      if (d_dim > tbox::Dimension(2)) {
         z = lowerx[2] + pdx[2] * ((*ci)(2) - loweri(2) + 0.5);
      }

      for (int d = 0; d < depth; d++) {
         (*data)(*ci, d) = d_Dcoef + d_Acoef * x + d_Bcoef * y + d_Ccoef * z;
      }

   }

}

void CellDataTest::setConservativeData(
   boost::shared_ptr<pdat::CellData<double> > data,
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
       * Set cell value on level zero to u(i,j,k) = (i + j + k)/ncells.
       */

      pdat::CellIterator fiend(sbox, false);
      for (pdat::CellIterator fi(sbox, true); fi != fiend; ++fi) {
         double value = 0.0;
         for (int d = 0; d < d_dim.getValue(); d++) {
            value += (double)((*fi)(d));
         }
         value /= ncells;
         for (int dep = 0; dep < depth; dep++) {
            (*data)(*fi, dep) = value;
         }
      }

   } else {

      /*
       * Set cell value on level > 0 to
       *    u(i,j,k) = u_c + ci*del_i + cj*del_j + ck*del_k
       * where u_c is underlying coarse value, (ci,cj,ck) is
       * the underlying coarse cell index, and (del_i,del_j,del_k)
       * is the vector between the coarse and fine cell centers.
       */

      hier::IntVector ratio(level->getRatioToLevelZero());

      boost::shared_ptr<geom::CartesianPatchGeometry> pgeom(
         patch.getPatchGeometry(),
         boost::detail::dynamic_cast_tag());
      const double* dx = pgeom->getDx();

      int coarse_ncells = ncells;
      std::vector<std::vector<double> > delta(d_dim.getValue());
      for (int d = 0; d < d_dim.getValue(); d++) {
         delta[d].resize(ratio(d), 0.0);
         coarse_ncells /= ratio(d);
         double coarse_dx = dx[d] * ratio(d);
         for (int i = 0; i < ratio(d); i++) {
            /*
             * delta[d][i] is the physical distance from i-th fine
             * cell centroid in d-direction to coarse cell centroid.
             * The distance is the d-th component of the displacement
             * vector.
             */
            delta[d][i] = (i + 0.5) * dx[d] - coarse_dx * 0.5;
         }
      }

      pdat::CellIterator fiend(sbox, false);
      for (pdat::CellIterator fi(sbox, true); fi != fiend; ++fi) {

         const hier::IntVector ci(hier::Index::coarsen(*fi, ratio));
         hier::IntVector del(ci.getDim());  // Index vector from ci to fi.
         double value = 0.0;
         for (int d = 0; d < d_dim.getValue(); d++) {
            del(d) = (int)delta[d][(*fi)(d) - ci(d) * ratio(d)];
            value += (double)(ci(d));
         }
         value /= coarse_ncells;

         for (int d = 0; d < d_dim.getValue(); d++) {
            value += ci(d) * del(d);
         }

         for (int dep = 0; dep < depth; dep++) {
            (*data)(*fi, dep) = value;
         }

      }

   }

}

void CellDataTest::setPeriodicData(
   boost::shared_ptr<pdat::CellData<double> > data,
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

   pdat::CellIterator ciend(sbox, false);
   for (pdat::CellIterator ci(sbox, true); ci != ciend; ++ci) {

      double val = 1.0;
      for (int d = 0; d < d_dim.getValue(); ++d) {
         double tmpf = dx[d] * ((*ci)(d) + 0.5) / domain_len[d];
         tmpf = sin(2 * M_PI * tmpf);
         val *= tmpf;
      }
      for (int d = 0; d < depth; d++) {
         (*data)(*ci, d) = val;
      }

   }

}

void CellDataTest::initializeDataOnPatch(
   const hier::Patch& patch,
   const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number,
   char src_or_dst)
{
   NULL_USE(src_or_dst);

   const hier::IntVector periodic_shift(d_cart_grid_geometry->getPeriodicShift(
                                           hier::IntVector(d_dim, 1)));
   bool is_periodic = periodic_shift.max() > 0;

   if (d_do_refine) {

      for (int i = 0; i < d_variables.getSize(); i++) {

         boost::shared_ptr<pdat::CellData<double> > cell_data(
            patch.getPatchData(d_variables[i], getDataContext()),
            boost::detail::dynamic_cast_tag());

         hier::Box dbox = cell_data->getBox();

         if (is_periodic) {
            setPeriodicData(cell_data, dbox, patch);
         } else {
            setLinearData(cell_data, dbox, patch);
         }

      }

   } else if (d_do_coarsen) {

      for (int i = 0; i < d_variables.getSize(); i++) {

         boost::shared_ptr<pdat::CellData<double> > cell_data(
            patch.getPatchData(d_variables[i], getDataContext()),
            boost::detail::dynamic_cast_tag());

         hier::Box dbox = cell_data->getGhostBox();

         if (is_periodic) {
            setPeriodicData(cell_data, dbox, patch);
         } else {
            setConservativeData(cell_data, dbox, patch, hierarchy, level_number);
         }

      }

   }

}

void CellDataTest::checkPatchInteriorData(
   const boost::shared_ptr<pdat::CellData<double> >& data,
   const hier::Box& interior,
   const hier::Patch& patch) const
{
   TBOX_ASSERT(data);

   const bool is_periodic =
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(d_dim,
            1)).max() > 0;

   const int depth = data->getDepth();

   boost::shared_ptr<pdat::CellData<double> > correct_data(
      new pdat::CellData<double>(
         data->getBox(),
         depth,
         data->getGhostCellWidth()));
   if (is_periodic) {
      setPeriodicData(correct_data, correct_data->getGhostBox(), patch);
   } else {
      setLinearData(correct_data, correct_data->getGhostBox(), patch);
   }

   pdat::CellIterator ciend(interior, false);
   for (pdat::CellIterator ci(interior, true); ci != ciend; ++ci) {
      for (int d = 0; d < depth; d++) {
         if (!(tbox::MathUtilities<double>::equalEps((*data)(*ci, d),
                  (*correct_data)(*ci, d)))) {
            tbox::perr << "FAILED: -- patch interior not properly filled"
                       << endl;
         }
      }
   }

}

void CellDataTest::setPhysicalBoundaryConditions(
   const hier::Patch& patch,
   const double time,
   const hier::IntVector& gcw_to_fill) const
{
   NULL_USE(time);

   const hier::IntVector periodic_shift(
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(d_dim, 1)));
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

      boost::shared_ptr<pdat::CellData<double> > cell_data(
         patch.getPatchData(d_variables[i], getDataContext()),
         boost::detail::dynamic_cast_tag());

      hier::Box patch_interior = cell_data->getBox();
      checkPatchInteriorData(cell_data, patch_interior, patch);

      /*
       * Set node boundary data.
       */
      for (int ni = 0; ni < num_node_bdry_boxes; ni++) {

         hier::Box fill_box = pgeom->getBoundaryFillBox(node_bdry[ni],
               patch.getBox(),
               gcw_to_fill);

         if (is_periodic) {
            setPeriodicData(cell_data, fill_box, patch);
         } else {
            setLinearData(cell_data, fill_box, patch);
         }
      }

      if (d_dim > tbox::Dimension(1)) {
         /*
          * Set edge boundary data.
          */
         for (int ei = 0; ei < num_edge_bdry_boxes; ei++) {

            hier::Box fill_box = pgeom->getBoundaryFillBox(edge_bdry[ei],
                  patch.getBox(),
                  gcw_to_fill);

            if (is_periodic) {
               setPeriodicData(cell_data, fill_box, patch);
            } else {
               setLinearData(cell_data, fill_box, patch);
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
                  gcw_to_fill);

            if (is_periodic) {
               setPeriodicData(cell_data, fill_box, patch);
            } else {
               setLinearData(cell_data, fill_box, patch);
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
bool CellDataTest::verifyResults(
   const hier::Patch& patch,
   const boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number)
{

   bool test_failed = false;

   const hier::IntVector periodic_shift(
      d_cart_grid_geometry->getPeriodicShift(hier::IntVector(d_dim, 1)));
   bool is_periodic = periodic_shift.max() > 0;

   if (d_do_refine || d_do_coarsen) {

      tbox::plog << "\nEntering CellDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl;

      hier::IntVector tgcw(periodic_shift.getDim(), 0);
      for (int i = 0; i < d_variables.getSize(); i++) {
         tgcw.max(patch.getPatchData(d_variables[i], getDataContext())->
            getGhostCellWidth());
      }
      hier::Box pbox = patch.getBox();

      boost::shared_ptr<pdat::CellData<double> > solution(
         new pdat::CellData<double>(pbox, 1, tgcw));

      hier::Box tbox(pbox);
      tbox.grow(tgcw);

      if (d_do_refine) {
         if (is_periodic) {
            setPeriodicData(solution, tbox, patch);
         } else {
            setLinearData(solution, tbox, patch);
         }
      } else {
         if (is_periodic) {
            setPeriodicData(solution, tbox, patch);
         } else {
            setConservativeData(solution, tbox, patch, hierarchy, level_number);
         }
      }

      for (int i = 0; i < d_variables.getSize(); i++) {

         boost::shared_ptr<pdat::CellData<double> > cell_data(
            patch.getPatchData(d_variables[i], getDataContext()),
            boost::detail::dynamic_cast_tag());
         int depth = cell_data->getDepth();
         hier::Box dbox = cell_data->getGhostBox();

         pdat::CellIterator ciend(dbox, false);
         for (pdat::CellIterator ci(dbox, true); ci != ciend; ++ci) {
            double correct = (*solution)(*ci);
            for (int d = 0; d < depth; d++) {
               double result = (*cell_data)(*ci, d);
               if (!tbox::MathUtilities<double>::equalEps(correct, result)) {
                  tbox::perr << "Test FAILED: ...."
                             << " : cell index = " << *ci
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
         tbox::plog << "CellDataTest Successful!" << endl;
      }

      solution.reset();   // just to be anal...

      tbox::plog << "\nExiting CellDataTest::verifyResults..." << endl;
      tbox::plog << "level_number = " << level_number << endl;
      tbox::plog << "Patch box = " << patch.getBox() << endl << endl;

   }

   return !test_failed;

}

}
