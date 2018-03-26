/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Class to test usage of boundary utilities
 *
 ************************************************************************/

#ifndef included_BoundaryDataTester
#define included_BoundaryDataTester

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/appu/BoundaryUtilityStrategy.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Array.h"

#include <boost/shared_ptr.hpp>
#include <string>

using namespace std;
using namespace SAMRAI;

class BoundaryDataTester:
   public xfer::RefinePatchStrategy,
   public appu::BoundaryUtilityStrategy
{
public:
   /**
    * The constructor reads variable data from input database.
    */
   BoundaryDataTester(
      const string& object_name,
      const tbox::Dimension& dim,
      boost::shared_ptr<tbox::Database> input_db,
      boost::shared_ptr<geom::CartesianGridGeometry> grid_geom);

   /**
    * Virtual destructor for BoundaryDataTester.
    */
   virtual ~BoundaryDataTester();

   /**
    * This routine is a concrete implementation of the virtual function
    * in the base class RefinePatchStrategy.  It sets the boundary
    * conditions for the variables.
    */
   void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double fill_time,
      const hier::IntVector& ghost_width_to_fill);

   /**
    * The next three functions are dummy implementations of the pure
    * virtual functions declared in the RefinePatchStrategy base class.
    * They are not needed for this example since we only have one level
    * in the hierarchy.
    */
   hier::IntVector getRefineOpStencilWidth() const {
      return hier::IntVector(d_dim, 0);
   }

   void preprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio)
   {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(fine_box);
      NULL_USE(ratio);
   }

   void postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio)
   {
      NULL_USE(fine);
      NULL_USE(coarse);
      NULL_USE(fine_box);
      NULL_USE(ratio);
   }

   /**
    * This routine is a concrete implementation of a virtual function
    * in the base class BoundaryUtilityStrategy.  It reads DIRICHLET
    * face or edge boundary state values from the given database with the
    * given name string idenifier.  The integer location index
    * indicates the face or edge to which the boundary condition applies.
    */
   void
   readDirichletBoundaryDataEntry(
      const boost::shared_ptr<tbox::Database>& db,
      string& db_name,
      int bdry_location_index);

   /**
    * This routine is a concrete implementation of a virtual function
    * in the base class BoundaryUtilityStrategy.  It reads NEUMANN
    * face or edge boundary state values from the given database with the
    * given name string idenifier.  The integer location index
    * indicates the face or edge to which the boundary condition applies.
    */
   void
   readNeumannBoundaryDataEntry(
      const boost::shared_ptr<tbox::Database>& db,
      string& db_name,
      int bdry_location_index);

   /**
    * Set data on patch interiors on given level in hierarchy.
    */
   void
   initializeDataOnPatchInteriors(
      boost::shared_ptr<hier::PatchHierarchy> hierarchy,
      int level_number);

   /**
    * Run boundary tests for given level in hierarchy and return integer
    * number of test failures.
    */
   int
   runBoundaryTest(
      boost::shared_ptr<hier::PatchHierarchy> hierarchy,
      int level_number);

   /**
    * Print all class data members to given output stream.
    */
   void
   printClassData(
      ostream& os) const;

private:
   /*
    * The object name is used for error/warning reporting.
    */
   string d_object_name;

   const tbox::Dimension d_dim;

   boost::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;

   /*
    * Arrays of information read from input file describing test variables
    */
   tbox::Array<string> d_variable_name;
   tbox::Array<int> d_variable_depth;
   tbox::Array<hier::IntVector> d_variable_num_ghosts;
   tbox::Array<tbox::Array<double> > d_variable_interior_values;

   /*
    * Items used to manage variables and data in test program.
    */
   tbox::Array<boost::shared_ptr<hier::Variable> > d_variables;
   boost::shared_ptr<hier::VariableContext> d_variable_context;
   hier::ComponentSelector d_patch_data_components;

   /*
    * Arrays of information read from input file for boundary conditions
    */
   tbox::Array<int> d_master_bdry_edge_conds;
   tbox::Array<int> d_scalar_bdry_edge_conds;
   tbox::Array<int> d_vector_bdry_edge_conds;

   tbox::Array<int> d_master_bdry_node_conds;
   tbox::Array<int> d_scalar_bdry_node_conds;
   tbox::Array<int> d_vector_bdry_node_conds;

   tbox::Array<int> d_master_bdry_face_conds; // Used only in 3D
   tbox::Array<int> d_scalar_bdry_face_conds; // Used only in 3D
   tbox::Array<int> d_vector_bdry_face_conds; // Used only in 3D

   tbox::Array<int> d_node_bdry_edge; // Used only in 2D
   tbox::Array<int> d_edge_bdry_face; // Used only in 3D
   tbox::Array<int> d_node_bdry_face; // Used only in 3D

   tbox::Array<tbox::Array<double> > d_variable_bc_values;

   int d_fail_count;

   /*
    * Private functions to perform tasks for boundary testing.
    */
   void
   readVariableInputAndMakeVariables(
      boost::shared_ptr<tbox::Database> db);
   void
   readBoundaryDataInput(
      boost::shared_ptr<tbox::Database> db);
   void
   readBoundaryDataStateEntry(
      boost::shared_ptr<tbox::Database> db,
      string& db_name,
      int bdry_location_index);
   void
   setBoundaryDataDefaults();
   void
   postprocessBoundaryInput();
   void
   checkBoundaryData(
      int btype,
      const hier::Patch& patch,
      const hier::IntVector& ghost_width_to_check);

};

#endif
