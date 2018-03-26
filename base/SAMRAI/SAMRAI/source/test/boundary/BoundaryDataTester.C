/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Class to test usage of boundary utilities
 *
 ************************************************************************/

#include "BoundaryDataTester.h"

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

//integer constants for boundary conditions
#define CHECK_BDRY_DATA (1)
#include "SAMRAI/appu/CartesianBoundaryDefines.h"

//integer constant for debugging improperly set boundary data
#define BOGUS_BDRY_DATA (-9999)

// routines for managing boundary data
#include "SAMRAI/appu/CartesianBoundaryUtilities2.h"
#include "SAMRAI/appu/CartesianBoundaryUtilities3.h"

/*
 *************************************************************************
 *
 * The constructor and destructor.
 *
 *************************************************************************
 */

BoundaryDataTester::BoundaryDataTester(
   const string& object_name,
   const tbox::Dimension& dim,
   boost::shared_ptr<tbox::Database> input_db,
   boost::shared_ptr<geom::CartesianGridGeometry> grid_geom):
   xfer::RefinePatchStrategy(dim),
   d_object_name(object_name),
   d_dim(dim),
   d_grid_geometry(grid_geom),
   d_variable_context(
      hier::VariableDatabase::getDatabase()->getContext("BOUNDARY_TEST"))
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!object_name.empty());
   TBOX_ASSERT(input_db);
#endif

   readVariableInputAndMakeVariables(input_db);

   setBoundaryDataDefaults();

   readBoundaryDataInput(input_db);

   postprocessBoundaryInput();

}

BoundaryDataTester::~BoundaryDataTester()
{
}

/*
 *************************************************************************
 *
 * Set physical boundary values for each variable acording to input data.
 *
 *************************************************************************
 */

void BoundaryDataTester::setPhysicalBoundaryConditions(
   hier::Patch& patch,
   const double fill_time,
   const hier::IntVector& ghost_width_to_fill)
{
   NULL_USE(fill_time);
   tbox::plog << "\n\nFilling boundary data on patch = " << patch.getBox()
              << endl;
   tbox::plog << "ghost_width_to_fill = " << ghost_width_to_fill << endl;

   for (int iv = 0; iv < d_variables.getSize(); iv++) {

      boost::shared_ptr<pdat::CellData<double> > cvdata(
         patch.getPatchData(d_variables[iv], d_variable_context),
         boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(cvdata);
#endif

      tbox::plog << "\n   iv = " << iv << " : " << d_variable_name[iv] << endl;
      tbox::plog << "   depth = " << cvdata->getDepth() << endl;

      hier::IntVector fill_gcw(hier::IntVector::min(cvdata->getGhostCellWidth(),
                                  ghost_width_to_fill));

      if (d_dim == tbox::Dimension(3)) {
         appu::CartesianBoundaryUtilities3::
         fillFaceBoundaryData(d_variable_name[iv], cvdata,
            patch,
            fill_gcw,
            ((cvdata->getDepth() > 1) ?
             d_vector_bdry_face_conds :
             d_scalar_bdry_face_conds),
            d_variable_bc_values[iv]);
         appu::CartesianBoundaryUtilities3::
         fillEdgeBoundaryData(d_variable_name[iv], cvdata,
            patch,
            fill_gcw,
            ((cvdata->getDepth() > 1) ?
             d_vector_bdry_edge_conds :
             d_scalar_bdry_edge_conds),
            d_variable_bc_values[iv]);

         appu::CartesianBoundaryUtilities3::
         fillNodeBoundaryData(d_variable_name[iv], cvdata,
            patch,
            fill_gcw,
            ((cvdata->getDepth() > 1) ?
             d_vector_bdry_node_conds :
             d_scalar_bdry_node_conds),
            d_variable_bc_values[iv]);
      }

      if (d_dim == tbox::Dimension(2)) {
         appu::CartesianBoundaryUtilities2::
         fillEdgeBoundaryData(d_variable_name[iv], cvdata,
            patch,
            fill_gcw,
            ((cvdata->getDepth() > 1) ?
             d_vector_bdry_edge_conds :
             d_scalar_bdry_edge_conds),
            d_variable_bc_values[iv]);

         appu::CartesianBoundaryUtilities2::
         fillNodeBoundaryData(d_variable_name[iv], cvdata,
            patch,
            fill_gcw,
            ((cvdata->getDepth() > 1) ?
             d_vector_bdry_node_conds :
             d_scalar_bdry_node_conds),
            d_variable_bc_values[iv]);
      }

   }

   if (d_dim == tbox::Dimension(2)) {
      checkBoundaryData(Bdry::EDGE2D, patch, ghost_width_to_fill);
      checkBoundaryData(Bdry::NODE2D, patch, ghost_width_to_fill);
   }
   if (d_dim == tbox::Dimension(3)) {
      checkBoundaryData(Bdry::FACE3D, patch, ghost_width_to_fill);
      checkBoundaryData(Bdry::EDGE3D, patch, ghost_width_to_fill);
      checkBoundaryData(Bdry::NODE3D, patch, ghost_width_to_fill);
   }

}

/*
 *************************************************************************
 *
 * Set data for each variable on patch interior acording to input data.
 *
 *************************************************************************
 */

void BoundaryDataTester::initializeDataOnPatchInteriors(
   boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(level_number == 0);

   boost::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(level_number));
   TBOX_ASSERT(level);

   level->allocatePatchData(d_patch_data_components);

   /*
    * Undefine the data so it is initialized to some value.
    */
   for (hier::PatchLevel::iterator ip(level->begin());
        ip != level->end(); ++ip) {
      const boost::shared_ptr<hier::Patch>& patch = *ip;

      for (int iv = 0; iv < d_variables.getSize(); iv++) {
         boost::shared_ptr<pdat::CellData<double> > cvdata(
            patch->getPatchData(d_variables[iv], d_variable_context),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(cvdata);
         cvdata->getArrayData().undefineData();
      }

   }

   for (hier::PatchLevel::iterator ip(level->begin());
        ip != level->end(); ++ip) {
      const boost::shared_ptr<hier::Patch>& patch = *ip;

      for (int iv = 0; iv < d_variables.getSize(); iv++) {
         boost::shared_ptr<pdat::CellData<double> > cvdata(
            patch->getPatchData(d_variables[iv], d_variable_context),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(cvdata);
         for (int id = 0; id < cvdata->getDepth(); id++) {
            cvdata->fill(d_variable_interior_values[iv][id],
               patch->getBox(),
               id);
         }
      }

   }

}

/*
 *************************************************************************
 *
 * Run boundary test:
 *
 *  1) register boundary filling operation for each variable
 *     with refine algorithm.
 *
 *  2) create communication schedule and fill data.
 *
 *  3) check all patch boundary values for correctness.
 *
 *************************************************************************
 */

int BoundaryDataTester::runBoundaryTest(
   boost::shared_ptr<hier::PatchHierarchy> hierarchy,
   int level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT(level_number == 0);
#endif

   int d_fail_count = 0;

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   xfer::RefineAlgorithm boundary_fill(d_dim);

   for (int iv = 0; iv < d_variables.getSize(); iv++) {
      int datid =
         variable_db->mapVariableAndContextToIndex(d_variables[iv],
            d_variable_context);

      boundary_fill.registerRefine(datid, datid, datid,
         boost::shared_ptr<hier::RefineOperator>());
   }

   boost::shared_ptr<hier::PatchLevel> level(
      hierarchy->getPatchLevel(level_number));
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(level);
#endif

   boundary_fill.createSchedule(level, this)->fillData(0.0);

   return d_fail_count;
}

/*
 *************************************************************************
 *
 * Read variable data from input, create variables,
 * and map variables into variable database.
 *
 *************************************************************************
 */

void BoundaryDataTester::readVariableInputAndMakeVariables(
   boost::shared_ptr<tbox::Database> db)
{
   TBOX_ASSERT(db);

   tbox::Array<string> var_keys = db->getAllKeys();
   int nkeys = var_keys.getSize();

   int var_cnt = 0;
   for (int i = 0; i < nkeys; i++) {
      if (db->getDatabase(var_keys[i])->keyExists("name")) {
         var_cnt++;
      }
   }

   d_variable_name.resizeArray(var_cnt);
   d_variable_depth.resizeArray(var_cnt);
   d_variable_num_ghosts.resizeArray(var_cnt, hier::IntVector(d_dim, 1));
   d_variable_interior_values.resizeArray(var_cnt);

   for (int i = 0; i < nkeys; i++) {

      boost::shared_ptr<tbox::Database> var_db(db->getDatabase(var_keys[i]));

      if (var_keys[i] != "Boundary_data" && var_db->keyExists("name")) {

         if (var_db->keyExists("name")) {
            d_variable_name[i] = var_db->getString("name");
         } else {
            TBOX_ERROR(d_object_name << ": "
                                     << "Variable input error: No 'name' string found for "
                                     << "key = " << var_keys[i] << endl);
         }

         if (var_db->keyExists("depth")) {
            d_variable_depth[i] = var_db->getInteger("depth");
         } else {
            d_variable_depth[i] = 1;
         }

         if (var_db->keyExists("num_ghosts")) {
            int* tmpg = &d_variable_num_ghosts[i][0];
            var_db->getIntegerArray("num_ghosts", tmpg, d_dim.getValue());
         }

         if (var_db->keyExists("interior_values")) {
            d_variable_interior_values[i].resizeArray(d_variable_depth[i]);
            var_db->getDoubleArray("interior_values",
               d_variable_interior_values[i].getPointer(),
               d_variable_depth[i]);
         } else {
            TBOX_ERROR(
               d_object_name << ": "
                             << "Variable input error: No 'interior_values' entry found for "
                             << "key = " << var_keys[i] << endl);
         }

      }

   }

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   d_variables.resizeArray(d_variable_name.getSize());

   for (int iv = 0; iv < d_variable_name.getSize(); iv++) {
      d_variables[iv].reset(
         new pdat::CellVariable<double>(d_dim, d_variable_name[iv],
                                        d_variable_depth[iv]));

      int datid =
         variable_db->registerVariableAndContext(d_variables[iv],
            d_variable_context,
            d_variable_num_ghosts[iv]);

      d_patch_data_components.setFlag(datid);
   }

}

/*
 *************************************************************************
 *
 * Set all boundary data to bogus default values for error checking.
 *
 *************************************************************************
 */

void BoundaryDataTester::setBoundaryDataDefaults()
{
   /*
    * Defaults for boundary conditions. Set to bogus values
    * for error checking.
    */

   if (d_dim == tbox::Dimension(2)) {
      d_master_bdry_edge_conds.resizeArray(NUM_2D_EDGES);
      d_scalar_bdry_edge_conds.resizeArray(NUM_2D_EDGES);
      d_vector_bdry_edge_conds.resizeArray(NUM_2D_EDGES);
      for (int ei = 0; ei < NUM_2D_EDGES; ei++) {
         d_master_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
         d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
         d_vector_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
      }

      d_master_bdry_node_conds.resizeArray(NUM_2D_NODES);
      d_scalar_bdry_node_conds.resizeArray(NUM_2D_NODES);
      d_vector_bdry_node_conds.resizeArray(NUM_2D_NODES);
      d_node_bdry_edge.resizeArray(NUM_2D_NODES);

      for (int ni = 0; ni < NUM_2D_NODES; ni++) {
         d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
         d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
         d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
         d_node_bdry_edge[ni] = BOGUS_BDRY_DATA;
      }
   }

   if (d_dim == tbox::Dimension(3)) {
      d_master_bdry_face_conds.resizeArray(NUM_3D_FACES);
      d_scalar_bdry_face_conds.resizeArray(NUM_3D_FACES);
      d_vector_bdry_face_conds.resizeArray(NUM_3D_FACES);
      for (int fi = 0; fi < NUM_3D_FACES; fi++) {
         d_master_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
         d_scalar_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
         d_vector_bdry_face_conds[fi] = BOGUS_BDRY_DATA;
      }

      d_master_bdry_edge_conds.resizeArray(NUM_3D_EDGES);
      d_scalar_bdry_edge_conds.resizeArray(NUM_3D_EDGES);
      d_vector_bdry_edge_conds.resizeArray(NUM_3D_EDGES);
      d_edge_bdry_face.resizeArray(NUM_3D_EDGES);
      for (int ei = 0; ei < NUM_3D_EDGES; ei++) {
         d_master_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
         d_scalar_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
         d_vector_bdry_edge_conds[ei] = BOGUS_BDRY_DATA;
         d_edge_bdry_face[ei] = BOGUS_BDRY_DATA;
      }

      d_master_bdry_node_conds.resizeArray(NUM_3D_NODES);
      d_scalar_bdry_node_conds.resizeArray(NUM_3D_NODES);
      d_vector_bdry_node_conds.resizeArray(NUM_3D_NODES);
      d_node_bdry_face.resizeArray(NUM_3D_NODES);

      for (int ni = 0; ni < NUM_3D_NODES; ni++) {
         d_master_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
         d_scalar_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
         d_vector_bdry_node_conds[ni] = BOGUS_BDRY_DATA;
         d_node_bdry_face[ni] = BOGUS_BDRY_DATA;
      }
   }

   d_variable_bc_values.resizeArray(d_variable_name.getSize());
   for (int iv = 0; iv < d_variable_name.getSize(); iv++) {
      if (d_dim == tbox::Dimension(2)) {
         d_variable_bc_values[iv].resizeArray(NUM_2D_EDGES
            * d_variable_depth[iv]);
      }
      if (d_dim == tbox::Dimension(3)) {
         d_variable_bc_values[iv].resizeArray(NUM_3D_FACES
            * d_variable_depth[iv]);
      }
      tbox::MathUtilities<double>::setArrayToSignalingNaN(d_variable_bc_values[
            iv]);
   }

}

/*
 *************************************************************************
 *
 * Functions to read boundary information from input database.
 *
 *************************************************************************
 */

void BoundaryDataTester::readDirichletBoundaryDataEntry(
   const boost::shared_ptr<tbox::Database>& db,
   string& db_name,
   int bdry_location_index)
{
   readBoundaryDataStateEntry(db, db_name, bdry_location_index);
}

void BoundaryDataTester::readNeumannBoundaryDataEntry(
   const boost::shared_ptr<tbox::Database>& db,
   string& db_name,
   int bdry_location_index)
{
   readBoundaryDataStateEntry(db, db_name, bdry_location_index);
}

void BoundaryDataTester::readBoundaryDataStateEntry(
   boost::shared_ptr<tbox::Database> db,
   string& db_name,
   int bdry_location_index)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(db);
   TBOX_ASSERT(!db_name.empty());
   TBOX_ASSERT(d_variable_bc_values.getSize() == d_variable_name.getSize());
#endif

   for (int iv = 0; iv < d_variable_name.getSize(); iv++) {

#ifdef DEBUG_CHECK_ASSERTIONS
      if (d_dim == tbox::Dimension(2)) {
         TBOX_ASSERT(d_variable_bc_values[iv].getSize() ==
            NUM_2D_EDGES * d_variable_depth[iv]);
      }
      if (d_dim == tbox::Dimension(3)) {
         TBOX_ASSERT(d_variable_bc_values[iv].getSize() ==
            NUM_3D_FACES * d_variable_depth[iv]);
      }
#endif

      if (db->keyExists(d_variable_name[iv])) {
         int depth = d_variable_depth[iv];
         tbox::Array<double> tmp_val(0);
         tmp_val = db->getDoubleArray(d_variable_name[iv]);
         if (tmp_val.getSize() < depth) {
            TBOX_ERROR(d_object_name << ": "
                                     << "Insufficient number of "
                                     << d_variable_name[iv] << " values given in "
                                     << db_name << " input database." << endl);
         }
         for (int id = 0; id < depth; id++) {
            d_variable_bc_values[iv][bdry_location_index * depth + id] =
               tmp_val[id];
         }
      } else {
         TBOX_ERROR(d_object_name << ": "
                                  << d_variable_name[iv]
                                  << " entry missing from " << db_name
                                  << " input database. " << endl);
      }

   }

}

void BoundaryDataTester::readBoundaryDataInput(
   boost::shared_ptr<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(db);
#endif

   hier::IntVector periodic(d_grid_geometry->getPeriodicShift(hier::IntVector(
                                  d_dim,
                                  1)));
   int num_per_dirs = 0;
   for (int id = 0; id < d_dim.getValue(); id++) {
      if (periodic(id)) num_per_dirs++;
   }

   if (num_per_dirs < d_dim.getValue()) {

      if (db->keyExists("Boundary_data")) {

         boost::shared_ptr<tbox::Database> bdry_db(
            db->getDatabase("Boundary_data"));

         if (d_dim == tbox::Dimension(2)) {
            appu::CartesianBoundaryUtilities2::
            readBoundaryInput(this,
               bdry_db,
               d_master_bdry_edge_conds,
               d_master_bdry_node_conds,
               periodic);
         }
         if (d_dim == tbox::Dimension(3)) {
            appu::CartesianBoundaryUtilities3::
            readBoundaryInput(this,
               bdry_db,
               d_master_bdry_face_conds,
               d_master_bdry_edge_conds,
               d_master_bdry_node_conds,
               periodic);
         }

      } else {
         TBOX_ERROR(
            d_object_name << ": "
                          << "Key data 'Boundary_data' not found in input. " << endl);
      }

   }

}

/*
 *************************************************************************
 *
 * Postprocess boundary data from input values
 * to make setting and checking easier.
 *
 *************************************************************************
 */

void BoundaryDataTester::postprocessBoundaryInput()
{
   if (d_dim == tbox::Dimension(2)) {
      for (int i = 0; i < NUM_2D_EDGES; i++) {
         d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
         d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];

         if (d_master_bdry_edge_conds[i] == BdryCond::REFLECT) {
            d_scalar_bdry_edge_conds[i] = BdryCond::FLOW;
         }
      }
      for (int i = 0; i < NUM_2D_NODES; i++) {
         d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
         d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];

         if (d_master_bdry_node_conds[i] == BdryCond::XREFLECT) {
            d_scalar_bdry_node_conds[i] = BdryCond::XFLOW;
         }
         if (d_master_bdry_node_conds[i] == BdryCond::YREFLECT) {
            d_scalar_bdry_node_conds[i] = BdryCond::YFLOW;
         }

         if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA) {
            d_node_bdry_edge[i] =
               appu::CartesianBoundaryUtilities2::getEdgeLocationForNodeBdry(
                  i, d_master_bdry_node_conds[i]);
         }
      }
   }
   if (d_dim == tbox::Dimension(3)) {
      for (int i = 0; i < NUM_3D_FACES; i++) {
         d_scalar_bdry_face_conds[i] = d_master_bdry_face_conds[i];
         d_vector_bdry_face_conds[i] = d_master_bdry_face_conds[i];

         if (d_master_bdry_face_conds[i] == BdryCond::REFLECT) {
            d_scalar_bdry_face_conds[i] = BdryCond::FLOW;
         }
      }

      for (int i = 0; i < NUM_3D_EDGES; i++) {
         d_scalar_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];
         d_vector_bdry_edge_conds[i] = d_master_bdry_edge_conds[i];

         if (d_master_bdry_edge_conds[i] == BdryCond::XREFLECT) {
            d_scalar_bdry_edge_conds[i] = BdryCond::XFLOW;
         }
         if (d_master_bdry_edge_conds[i] == BdryCond::YREFLECT) {
            d_scalar_bdry_edge_conds[i] = BdryCond::YFLOW;
         }
         if (d_master_bdry_edge_conds[i] == BdryCond::ZREFLECT) {
            d_scalar_bdry_edge_conds[i] = BdryCond::ZFLOW;
         }

         if (d_master_bdry_edge_conds[i] != BOGUS_BDRY_DATA) {
            d_edge_bdry_face[i] =
               appu::CartesianBoundaryUtilities3::getFaceLocationForEdgeBdry(
                  i, d_master_bdry_edge_conds[i]);
         }
      }

      for (int i = 0; i < NUM_3D_NODES; i++) {
         d_scalar_bdry_node_conds[i] = d_master_bdry_node_conds[i];
         d_vector_bdry_node_conds[i] = d_master_bdry_node_conds[i];

         if (d_master_bdry_node_conds[i] == BdryCond::XREFLECT) {
            d_scalar_bdry_node_conds[i] = BdryCond::XFLOW;
         }
         if (d_master_bdry_node_conds[i] == BdryCond::YREFLECT) {
            d_scalar_bdry_node_conds[i] = BdryCond::YFLOW;
         }
         if (d_master_bdry_node_conds[i] == BdryCond::ZREFLECT) {
            d_scalar_bdry_node_conds[i] = BdryCond::ZFLOW;
         }

         if (d_master_bdry_node_conds[i] != BOGUS_BDRY_DATA) {
            d_node_bdry_face[i] =
               appu::CartesianBoundaryUtilities3::getFaceLocationForNodeBdry(
                  i, d_master_bdry_node_conds[i]);
         }
      }
   }

}

/*
 *************************************************************************
 *
 * Check boundary values on patch for correctness.
 *
 *************************************************************************
 */

void BoundaryDataTester::checkBoundaryData(
   int btype,
   const hier::Patch& patch,
   const hier::IntVector& ghost_width_to_check)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_dim == tbox::Dimension(2)) {
      TBOX_ASSERT(btype == Bdry::EDGE2D ||
         btype == Bdry::NODE2D);
   }
   if (d_dim == tbox::Dimension(3)) {
      TBOX_ASSERT(btype == Bdry::FACE3D ||
         btype == Bdry::EDGE3D ||
         btype == Bdry::NODE3D);
   }
#endif

   const boost::shared_ptr<geom::CartesianPatchGeometry> pgeom(
      patch.getPatchGeometry(),
      boost::detail::dynamic_cast_tag());
   const tbox::Array<hier::BoundaryBox> bdry_boxes =
      pgeom->getCodimensionBoundaries(btype);

   for (int i = 0; i < bdry_boxes.getSize(); i++) {
      hier::BoundaryBox bbox = bdry_boxes[i];
#ifdef DEBUG_CHECK_ASSERTIONS
      TBOX_ASSERT(bbox.getBoundaryType() == btype);
#endif
      int bloc = bbox.getLocationIndex();

      for (int iv = 0; iv < d_variables.getSize(); iv++) {
         boost::shared_ptr<pdat::CellData<double> > cvdata(
            patch.getPatchData(d_variables[iv], d_variable_context),
            boost::detail::dynamic_cast_tag());
#ifdef DEBUG_CHECK_ASSERTIONS
         TBOX_ASSERT(cvdata);
#endif

         int depth = d_variable_depth[iv];

         int bscalarcase = 0;
         int bvectorcase = 0;
         int refbdryloc = 0;
         if (d_dim == tbox::Dimension(2)) {
            if (btype == Bdry::EDGE2D) {
               bscalarcase = d_scalar_bdry_edge_conds[bloc];
               bvectorcase = d_vector_bdry_edge_conds[bloc];
               refbdryloc = bloc;
            } else { // btype == Bdry::NODE2D
               bscalarcase = d_scalar_bdry_node_conds[bloc];
               bvectorcase = d_vector_bdry_node_conds[bloc];
               refbdryloc = d_node_bdry_edge[bloc];
            }
         }
         if (d_dim == tbox::Dimension(3)) {
            if (btype == Bdry::FACE3D) {
               bscalarcase = d_scalar_bdry_face_conds[bloc];
               bvectorcase = d_vector_bdry_face_conds[bloc];
               refbdryloc = bloc;
            } else if (btype == Bdry::EDGE3D) {
               bscalarcase = d_scalar_bdry_edge_conds[bloc];
               bvectorcase = d_vector_bdry_edge_conds[bloc];
               refbdryloc = d_edge_bdry_face[bloc];
            } else { // btype == Bdry::NODE3D
               bscalarcase = d_scalar_bdry_node_conds[bloc];
               bvectorcase = d_vector_bdry_node_conds[bloc];
               refbdryloc = d_node_bdry_face[bloc];
            }
         }

         int data_id = hier::VariableDatabase::getDatabase()->
            mapVariableAndContextToIndex(d_variables[iv], d_variable_context);

         int num_bad_values = 0;

         if (depth == 1) {

            if (d_dim == tbox::Dimension(2)) {
               num_bad_values =
                  appu::CartesianBoundaryUtilities2::
                  checkBdryData(d_variable_name[iv],
                     patch,
                     data_id,
                     0,
                     ghost_width_to_check,
                     bbox,
                     bscalarcase,
                     d_variable_bc_values[iv][refbdryloc]);
            }
            if (d_dim == tbox::Dimension(3)) {
               num_bad_values =
                  appu::CartesianBoundaryUtilities3::
                  checkBdryData(d_variable_name[iv],
                     patch,
                     data_id,
                     0,
                     ghost_width_to_check,
                     bbox,
                     bscalarcase,
                     d_variable_bc_values[iv][refbdryloc]);
            }
#if (TESTING == 1)
            if (num_bad_values > 0) {
               d_fail_count++;
               tbox::perr << "\nBoundary Test FAILED: \n"
                          << "     " << num_bad_values << " bad "
                          << d_variable_name[iv] << " values found for"
                          << "     boundary type " << btype
                          << " at location "
                          << bloc << endl;
            }
#endif

         } else {
            for (int id = 0; id < depth; id++) {
               int vbcase = bscalarcase;
               if (d_dim == tbox::Dimension(2)) {
                  if (btype == Bdry::EDGE2D) {
                     if ((id == 0 && (bloc == BdryLoc::XLO ||
                                      bloc == BdryLoc::XHI)) ||
                         (id == 1 && (bloc == BdryLoc::YLO ||
                                      bloc == BdryLoc::YHI))) {
                        vbcase = bvectorcase;
                     }
                  } else {
                     if ((id == 0 && bvectorcase == BdryCond::XREFLECT) ||
                         (id == 1 && bvectorcase == BdryCond::YREFLECT)) {
                        vbcase = bvectorcase;
                     }
                  }
               }
               if (d_dim == tbox::Dimension(3)) {
                  if (btype == Bdry::FACE3D) {
                     if ((id == 0 && (bloc == BdryLoc::XLO ||
                                      bloc == BdryLoc::XHI)) ||
                         (id == 1 && (bloc == BdryLoc::YLO ||
                                      bloc == BdryLoc::YHI)) ||
                         (id == 2 && (bloc == BdryLoc::ZLO ||
                                      bloc == BdryLoc::ZHI))) {
                        vbcase = bvectorcase;
                     }
                  } else {
                     if ((id == 0 && bvectorcase == BdryCond::XREFLECT) ||
                         (id == 1 && bvectorcase == BdryCond::YREFLECT) ||
                         (id == 2 && bvectorcase == BdryCond::ZREFLECT)) {
                        vbcase = bvectorcase;
                     }
                  }
               }

               if (d_dim == tbox::Dimension(2)) {
                  num_bad_values =
                     appu::CartesianBoundaryUtilities2::
                     checkBdryData(d_variable_name[iv],
                        patch,
                        data_id,
                        id,
                        ghost_width_to_check,
                        bbox,
                        vbcase,
                        d_variable_bc_values[iv][refbdryloc * depth + id]);
               }
               if (d_dim == tbox::Dimension(3)) {
                  num_bad_values =
                     appu::CartesianBoundaryUtilities3::
                     checkBdryData(d_variable_name[iv],
                        patch,
                        data_id,
                        id,
                        ghost_width_to_check,
                        bbox,
                        vbcase,
                        d_variable_bc_values[iv][refbdryloc * depth + id]);
               }
#if (TESTING == 1)
               if (num_bad_values > 0) {
                  d_fail_count++;
                  tbox::perr << "\nBoundary Test FAILED: \n"
                             << "     " << num_bad_values << " bad "
                             << d_variable_name[iv] << " values found for"
                             << "     boundary type " << btype
                             << " at location "
                             << bloc << endl;
               }
#endif

            }  // for (int id = 0; id < depth; id++)

         }  // else

      }   // for (int iv = 0; iv < d_variables.getSize(); iv++)

   }  // for (int i = 0; i < bdry_boxes.getSize(); i++ )

}

/*
 *************************************************************************
 *
 * Write all class data members to specified output stream.
 *
 *************************************************************************
 */

void BoundaryDataTester::printClassData(
   ostream& os) const
{
   int i, j;
   os << "\nBoundaryDataTester::printClassData..." << endl;
   os << "BoundaryDataTester: this = " << (BoundaryDataTester *)this << endl;
   os << "d_object_name = " << d_object_name << endl;
   os << "d_grid_geometry = "
      << d_grid_geometry.get() << endl;

   if (d_variable_context) {
      os << "d_variable_context = "
         << d_variable_context->getName() << endl;
   } else {
      os << "d_variable_context = NULL" << endl;
   }

   os << "\nVariables ...\n" << endl;
   for (i = 0; i < d_variable_name.getSize(); i++) {
      os << "Variable " << i << endl;
      os << "   name       = " << d_variable_name[i] << endl;
      os << "   depth      = " << d_variable_depth[i] << endl;
      os << "   num_ghosts = " << d_variable_num_ghosts[i] << endl;
      os << "   interior_values = " << d_variable_interior_values[i][0];
      for (j = 1; j < d_variable_depth[i]; j++) {
         os << " ,  " << d_variable_interior_values[i][j];
      }
      os << endl;
   }

   os << "\n   Boundary condition data... " << endl;

   if (d_dim == tbox::Dimension(2)) {
      for (j = 0; j < d_master_bdry_edge_conds.getSize(); j++) {
         os << "\n       d_master_bdry_edge_conds[" << j << "] = "
            << d_master_bdry_edge_conds[j] << endl;
         os << "       d_scalar_bdry_edge_conds[" << j << "] = "
            << d_scalar_bdry_edge_conds[j] << endl;
         os << "       d_vector_bdry_edge_conds[" << j << "] = "
            << d_vector_bdry_edge_conds[j] << endl;
         if (d_master_bdry_edge_conds[j] == BdryCond::DIRICHLET ||
             d_master_bdry_edge_conds[j] == BdryCond::NEUMANN) {
            for (i = 0; i < d_variable_name.getSize(); i++) {
               os << d_variable_name[i] << " bdry edge value[" << j << "] = "
                  << d_variable_bc_values[i][j * d_variable_depth[i]];
               for (int id = 1; id < d_variable_depth[i]; id++) {
                  os << " , "
                     << d_variable_bc_values[i][j * d_variable_depth[i] + id];
               }
               os << endl;
            }
         }
      }
      os << endl;
      for (j = 0; j < d_master_bdry_node_conds.getSize(); j++) {
         os << "\n       d_master_bdry_node_conds[" << j << "] = "
            << d_master_bdry_node_conds[j] << endl;
         os << "       d_scalar_bdry_node_conds[" << j << "] = "
            << d_scalar_bdry_node_conds[j] << endl;
         os << "       d_vector_bdry_node_conds[" << j << "] = "
            << d_vector_bdry_node_conds[j] << endl;
         os << "       d_node_bdry_edge[" << j << "] = "
            << d_node_bdry_edge[j] << endl;
      }
   }
   if (d_dim == tbox::Dimension(3)) {
      for (j = 0; j < d_master_bdry_face_conds.getSize(); j++) {
         os << "\n       d_master_bdry_face_conds[" << j << "] = "
            << d_master_bdry_face_conds[j] << endl;
         os << "       d_scalar_bdry_face_conds[" << j << "] = "
            << d_scalar_bdry_face_conds[j] << endl;
         os << "       d_vector_bdry_face_conds[" << j << "] = "
            << d_vector_bdry_face_conds[j] << endl;
         if (d_master_bdry_face_conds[j] == BdryCond::DIRICHLET) {
            for (i = 0; i < d_variable_name.getSize(); i++) {
               os << d_variable_name[i] << " bdry edge value[" << j << "] = "
                  << d_variable_bc_values[i][j * d_variable_depth[i]];
               for (int id = 1; id < d_variable_depth[i]; id++) {
                  os << " , "
                     << d_variable_bc_values[i][j * d_variable_depth[i] + id];
               }
               os << endl;
            }
         }
      }
      os << endl;
      for (j = 0; j < d_master_bdry_edge_conds.getSize(); j++) {
         os << "\n       d_master_bdry_edge_conds[" << j << "] = "
            << d_master_bdry_edge_conds[j] << endl;
         os << "       d_scalar_bdry_edge_conds[" << j << "] = "
            << d_scalar_bdry_edge_conds[j] << endl;
         os << "       d_vector_bdry_edge_conds[" << j << "] = "
            << d_vector_bdry_edge_conds[j] << endl;
         os << "       d_edge_bdry_face[" << j << "] = "
            << d_edge_bdry_face[j] << endl;
      }
      os << endl;
      for (j = 0; j < d_master_bdry_node_conds.getSize(); j++) {
         os << "\n       d_master_bdry_node_conds[" << j << "] = "
            << d_master_bdry_node_conds[j] << endl;
         os << "       d_scalar_bdry_node_conds[" << j << "] = "
            << d_scalar_bdry_node_conds[j] << endl;
         os << "       d_vector_bdry_node_conds[" << j << "] = "
            << d_vector_bdry_node_conds[j] << endl;
         os << "       d_node_bdry_face[" << j << "] = "
            << d_node_bdry_face[j] << endl;
      }
   }

}
