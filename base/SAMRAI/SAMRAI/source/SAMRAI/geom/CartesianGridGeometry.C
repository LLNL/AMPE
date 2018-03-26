/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Simple Cartesian grid geometry for an AMR hierarchy.
 *
 ************************************************************************/

#ifndef included_geom_CartesianGridGeometry_C
#define included_geom_CartesianGridGeometry_C

#include "SAMRAI/geom/CartesianGridGeometry.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"

// Cell data coarsen operators
#include "SAMRAI/geom/CartesianCellComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianCellDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianCellFloatWeightedAverage.h"

// Cell data refine operators
#include "SAMRAI/geom/CartesianCellComplexConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianCellComplexLinearRefine.h"
#include "SAMRAI/geom/CartesianCellDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianCellDoubleLinearRefine.h"
#include "SAMRAI/geom/CartesianCellFloatConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianCellFloatLinearRefine.h"

// Edge data coarsen operators
#include "SAMRAI/geom/CartesianEdgeComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianEdgeDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianEdgeFloatWeightedAverage.h"

// Edge data refine operators
#include "SAMRAI/geom/CartesianEdgeDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianEdgeFloatConservativeLinearRefine.h"

// Face data coarsen operators
#include "SAMRAI/geom/CartesianFaceComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianFaceDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianFaceFloatWeightedAverage.h"

// Face data refine operators
#include "SAMRAI/geom/CartesianFaceDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianFaceFloatConservativeLinearRefine.h"

// Node data refine operators
#include "SAMRAI/geom/CartesianNodeComplexLinearRefine.h"
#include "SAMRAI/geom/CartesianNodeDoubleLinearRefine.h"
#include "SAMRAI/geom/CartesianNodeFloatLinearRefine.h"

// Outerface data coarsen operators
#include "SAMRAI/geom/CartesianOuterfaceComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianOuterfaceDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianOuterfaceFloatWeightedAverage.h"

// Outerside data coarsen operators
#include "SAMRAI/geom/CartesianOutersideDoubleWeightedAverage.h"

// Side data coarsen operators
#include "SAMRAI/geom/CartesianSideComplexWeightedAverage.h"
#include "SAMRAI/geom/CartesianSideDoubleWeightedAverage.h"
#include "SAMRAI/geom/CartesianSideFloatWeightedAverage.h"

// Side data refine operators
#include "SAMRAI/geom/CartesianSideDoubleConservativeLinearRefine.h"
#include "SAMRAI/geom/CartesianSideFloatConservativeLinearRefine.h"

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/EdgeVariable.h"
#include "SAMRAI/pdat/FaceVariable.h"
#include "SAMRAI/pdat/NodeVariable.h"
#include "SAMRAI/pdat/OuterfaceVariable.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/pdat/SideVariable.h"

#include "SAMRAI/hier/BoundaryLookupTable.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <boost/make_shared.hpp>
#include <cstdlib>
#include <fstream>
#include <typeinfo>

namespace SAMRAI {
namespace geom {

const int CartesianGridGeometry::GEOM_CARTESIAN_GRID_GEOMETRY_VERSION = 2;

// using namespace std;

/*
 *************************************************************************
 *
 * Constructors for CartesianGridGeometry.  Both set up operator
 * handlers and register the geometry object with the RestartManager.
 * However, one initializes data members based on arguments.
 * The other initializes the object based on input file information.
 *
 *************************************************************************
 */
CartesianGridGeometry::CartesianGridGeometry(
   const tbox::Dimension& dim,
   const std::string& object_name,
   const boost::shared_ptr<tbox::Database>& input_db,
   bool register_for_restart):
   GridGeometry(dim, object_name),
   d_domain_box(dim)
{
   TBOX_ASSERT(input_db);

   d_registered_for_restart = register_for_restart;

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
      registerRestartItem(getObjectName(), this);
   }

   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();
   if (is_from_restart && d_registered_for_restart) {
      getFromRestart();
   }

   getFromInput(input_db, is_from_restart);
}

CartesianGridGeometry::CartesianGridGeometry(
   const std::string& object_name,
   const double* x_lo,
   const double* x_up,
   const hier::BoxContainer& domain,
   bool register_for_restart):
   GridGeometry(domain.front().getDim(), object_name),
   d_domain_box(domain.front().getDim())
{
   TBOX_ASSERT(domain.size() > 0);
   TBOX_ASSERT(!(x_lo == (double *)NULL));
   TBOX_ASSERT(!(x_up == (double *)NULL));

   d_registered_for_restart = register_for_restart;

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
      registerRestartItem(getObjectName(), this);
   }

   setGeometryData(x_lo, x_up, domain);
}

CartesianGridGeometry::CartesianGridGeometry(
   const std::string& object_name,
   const double* x_lo,
   const double* x_up,
   const hier::BoxContainer& domain,
   const boost::shared_ptr<hier::TransferOperatorRegistry>& op_reg,
   bool register_for_restart) :
   GridGeometry(domain.front().getDim(), object_name, op_reg),
   d_domain_box(domain.front().getDim())
{
   TBOX_ASSERT(domain.size() > 0);
   TBOX_ASSERT(!(x_lo == (double *)NULL));
   TBOX_ASSERT(!(x_up == (double *)NULL));

   d_registered_for_restart = register_for_restart;

   if (d_registered_for_restart) {
      tbox::RestartManager::getManager()->
      registerRestartItem(getObjectName(), this);
   }

   setGeometryData(x_lo, x_up, domain);
}

/*
 *************************************************************************
 *
 * Destructor for CartesianGridGeometry deallocates grid storage.
 *
 *************************************************************************
 */

CartesianGridGeometry::~CartesianGridGeometry()
{
}

/*
 *************************************************************************
 *
 * Create and return pointer to refined version of this Cartesian
 * grid geometry object refined by the given ratio.
 *
 *************************************************************************
 */

boost::shared_ptr<hier::BaseGridGeometry>
CartesianGridGeometry::makeRefinedGridGeometry(
   const std::string& fine_geom_name,
   const hier::IntVector& refine_ratio,
   bool register_for_restart) const
{
   const tbox::Dimension dim(getDim());

   TBOX_ASSERT(!fine_geom_name.empty());
   TBOX_ASSERT(fine_geom_name != getObjectName());
   TBOX_ASSERT(refine_ratio > hier::IntVector::getZero(dim));

   hier::BoxContainer fine_domain(getPhysicalDomain());
   fine_domain.refine(refine_ratio);

   boost::shared_ptr<hier::BaseGridGeometry> fine_geometry(
      boost::make_shared<CartesianGridGeometry>(fine_geom_name,
         d_x_lo,
         d_x_up,
         fine_domain,
         d_transfer_operator_registry,
         register_for_restart));

   fine_geometry->initializePeriodicShift(getPeriodicShift(hier::
         IntVector::getOne(dim)));

   return fine_geometry;
}

/*
 *************************************************************************
 *
 * Create and return pointer to coarsened version of this Cartesian
 * grid geometry object coarsened by the given ratio.
 *
 *************************************************************************
 */

boost::shared_ptr<hier::BaseGridGeometry>
CartesianGridGeometry::makeCoarsenedGridGeometry(
   const std::string& coarse_geom_name,
   const hier::IntVector& coarsen_ratio,
   bool register_for_restart) const
{
   const tbox::Dimension& dim(getDim());

   TBOX_ASSERT(!coarse_geom_name.empty());
   TBOX_ASSERT(coarse_geom_name != getObjectName());
   TBOX_ASSERT(coarsen_ratio > hier::IntVector::getZero(dim));

   hier::BoxContainer coarse_domain(getPhysicalDomain());
   coarse_domain.coarsen(coarsen_ratio);

   /*
    * Need to check that domain can be coarsened by given ratio.
    */
   const hier::BoxContainer& fine_domain = getPhysicalDomain();
   const int nboxes = fine_domain.size();
   hier::BoxContainer::const_iterator fine_domain_itr(fine_domain);
   hier::BoxContainer::iterator coarse_domain_itr(coarse_domain);
   for (int ib = 0; ib < nboxes; ib++, ++fine_domain_itr, ++coarse_domain_itr) {
      hier::Box testbox = hier::Box::refine(*coarse_domain_itr, coarsen_ratio);
      if (!testbox.isSpatiallyEqual(*fine_domain_itr)) {
#ifdef DEBUG_CHECK_ASSERTIONS
         tbox::plog
         << "CartesianGridGeometry::makeCoarsenedGridGeometry : Box # "
         << ib << std::endl;
         tbox::plog << "      fine box = " << *fine_domain_itr << std::endl;
         tbox::plog << "      coarse box = " << *coarse_domain_itr << std::endl;
         tbox::plog << "      refined coarse box = " << testbox << std::endl;
#endif
         TBOX_ERROR(
            "geom::CartesianGridGeometry::makeCoarsenedGridGeometry() error...\n"
            << "    geometry object with name = " << getObjectName()
            << "\n    Cannot be coarsened by ratio " << coarsen_ratio
            << std::endl);
      }
   }

   boost::shared_ptr<hier::BaseGridGeometry> coarse_geometry(
      boost::make_shared<CartesianGridGeometry>(coarse_geom_name,
         d_x_lo,
         d_x_up,
         coarse_domain,
         d_transfer_operator_registry,
         register_for_restart));

   coarse_geometry->initializePeriodicShift(getPeriodicShift(hier::
         IntVector::getOne(dim)));

   return coarse_geometry;
}

/*
 *************************************************************************
 *
 * Set data members for this geometry object based on arguments.
 *
 *************************************************************************
 */

void
CartesianGridGeometry::setGeometryData(
   const double* x_lo,
   const double* x_up,
   const hier::BoxContainer& domain)
{
   const tbox::Dimension& dim(getDim());

   TBOX_ASSERT(!(x_lo == (double *)NULL));
   TBOX_ASSERT(!(x_up == (double *)NULL));

   for (int id = 0; id < dim.getValue(); id++) {
      d_x_lo[id] = x_lo[id];
      d_x_up[id] = x_up[id];
   }

   setPhysicalDomain(domain, 1);

   hier::Box bigbox(dim);
   const hier::BoxContainer& block_domain = getPhysicalDomain();
   for (hier::BoxContainer::const_iterator k(block_domain);
        k != block_domain.end(); ++k) {
      bigbox += *k;
   }

   d_domain_box = bigbox;

   hier::IntVector ncells = d_domain_box.numberCells();
   for (int id2 = 0; id2 < dim.getValue(); id2++) {
      double length = d_x_up[id2] - d_x_lo[id2];
      d_dx[id2] = length / ((double)ncells(id2));
   }
}

/*
 *************************************************************************
 *
 * Create CartesianPatchGeometry geometry object, initializing its
 * boundary and grid information and assign it to the given patch.
 *
 *************************************************************************
 */

void
CartesianGridGeometry::setGeometryDataOnPatch(
   hier::Patch& patch,
   const hier::IntVector& ratio_to_level_zero,
   const TwoDimBool& touches_regular_bdry,
   const TwoDimBool& touches_periodic_bdry) const
{
   const tbox::Dimension& dim(getDim());

   TBOX_DIM_ASSERT_CHECK_DIM_ARGS2(dim, patch, ratio_to_level_zero);

#ifdef DEBUG_CHECK_ASSERTIONS
   /*
    * All components of ratio must be nonzero.  Additionally,
    * all components not equal to 1 must have the same sign.
    */
   TBOX_ASSERT(ratio_to_level_zero != hier::IntVector::getZero(dim));

   if (dim > tbox::Dimension(1)) {
      for (int i = 0; i < dim.getValue(); i++) {
         TBOX_ASSERT((ratio_to_level_zero(i)
                      * ratio_to_level_zero((i + 1) % dim.getValue()) > 0)
            || (ratio_to_level_zero(i) == 1)
            || (ratio_to_level_zero((i + 1) % dim.getValue()) == 1));
      }
   }
#endif

   double dx[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   double x_lo[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   double x_up[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   bool coarsen = false;
   if (ratio_to_level_zero(0) < 0) coarsen = true;
   hier::IntVector tmp_rat = (ratio_to_level_zero);
   for (int id2 = 0; id2 < dim.getValue(); id2++) {
      tmp_rat(id2) = abs(ratio_to_level_zero(id2));
   }

   hier::Box index_box = d_domain_box;
   hier::Box box = patch.getBox();

   if (coarsen) {
      index_box.coarsen(tmp_rat);
      for (int id3 = 0; id3 < dim.getValue(); id3++) {
         dx[id3] = d_dx[id3] * ((double)tmp_rat(id3));
      }
   } else {
      index_box.refine(tmp_rat);
      for (int id4 = 0; id4 < dim.getValue(); id4++) {
         dx[id4] = d_dx[id4] / ((double)tmp_rat(id4));
      }
   }

   for (int id5 = 0; id5 < dim.getValue(); id5++) {
      x_lo[id5] = d_x_lo[id5]
         + ((double)(box.lower(id5) - index_box.lower(id5))) * dx[id5];
      x_up[id5] = x_lo[id5] + ((double)box.numberCells(id5)) * dx[id5];
   }

   boost::shared_ptr<CartesianPatchGeometry> geom(
      boost::make_shared<CartesianPatchGeometry>(ratio_to_level_zero,
         touches_regular_bdry,
         touches_periodic_bdry,
         dx, x_lo, x_up));

   patch.setPatchGeometry(geom);

}

void
CartesianGridGeometry::buildOperators()
{
   GridGeometry::buildOperators();

   // CartesianGridGeometry specific Coarsening Operators
   addCoarsenOperator(
      typeid(pdat::CellVariable<dcomplex>).name(),
      boost::make_shared<CartesianCellComplexWeightedAverage>(d_dim));
   addCoarsenOperator(
      typeid(pdat::CellVariable<double>).name(),
      boost::make_shared<CartesianCellDoubleWeightedAverage>(d_dim));
   addCoarsenOperator(
      typeid(pdat::CellVariable<float>).name(),
      boost::make_shared<CartesianCellFloatWeightedAverage>(d_dim));
   addCoarsenOperator(
      typeid(pdat::EdgeVariable<dcomplex>).name(),
      boost::make_shared<CartesianEdgeComplexWeightedAverage>(d_dim));
   addCoarsenOperator(
      typeid(pdat::EdgeVariable<double>).name(),
      boost::make_shared<CartesianEdgeDoubleWeightedAverage>(d_dim));
   addCoarsenOperator(
      typeid(pdat::EdgeVariable<float>).name(),
      boost::make_shared<CartesianEdgeFloatWeightedAverage>(d_dim));
   addCoarsenOperator(
      typeid(pdat::FaceVariable<dcomplex>).name(),
      boost::make_shared<CartesianFaceComplexWeightedAverage>(d_dim));
   addCoarsenOperator(
      typeid(pdat::FaceVariable<double>).name(),
      boost::make_shared<CartesianFaceDoubleWeightedAverage>(d_dim));
   addCoarsenOperator(
      typeid(pdat::FaceVariable<float>).name(),
      boost::make_shared<CartesianFaceFloatWeightedAverage>(d_dim));
   addCoarsenOperator(
      typeid(pdat::OuterfaceVariable<dcomplex>).name(),
      boost::make_shared<CartesianOuterfaceComplexWeightedAverage>(d_dim));
   addCoarsenOperator(
      typeid(pdat::OuterfaceVariable<double>).name(),
      boost::make_shared<CartesianOuterfaceDoubleWeightedAverage>(d_dim));
   addCoarsenOperator(
      typeid(pdat::OuterfaceVariable<float>).name(),
      boost::make_shared<CartesianOuterfaceFloatWeightedAverage>(d_dim));
   addCoarsenOperator(
      typeid(pdat::OutersideVariable<double>).name(),
      boost::make_shared<CartesianOutersideDoubleWeightedAverage>(d_dim));
   addCoarsenOperator(
      typeid(pdat::SideVariable<dcomplex>).name(),
      boost::make_shared<CartesianSideComplexWeightedAverage>(d_dim));
   addCoarsenOperator(
      typeid(pdat::SideVariable<double>).name(),
      boost::make_shared<CartesianSideDoubleWeightedAverage>(d_dim));
   addCoarsenOperator(
      typeid(pdat::SideVariable<float>).name(),
      boost::make_shared<CartesianSideFloatWeightedAverage>(d_dim));

   // CartesianGridGeometry specific Refinement Operators
   addRefineOperator(
      typeid(pdat::CellVariable<dcomplex>).name(),
      boost::make_shared<CartesianCellComplexConservativeLinearRefine>(d_dim));
   addRefineOperator(
      typeid(pdat::CellVariable<double>).name(),
      boost::make_shared<CartesianCellDoubleConservativeLinearRefine>(d_dim));
   addRefineOperator(
      typeid(pdat::CellVariable<float>).name(),
      boost::make_shared<CartesianCellFloatConservativeLinearRefine>(d_dim));
   addRefineOperator(
      typeid(pdat::EdgeVariable<double>).name(),
      boost::make_shared<CartesianEdgeDoubleConservativeLinearRefine>(d_dim));
   addRefineOperator(
      typeid(pdat::EdgeVariable<float>).name(),
      boost::make_shared<CartesianEdgeFloatConservativeLinearRefine>(d_dim));
   addRefineOperator(
      typeid(pdat::FaceVariable<double>).name(),
      boost::make_shared<CartesianFaceDoubleConservativeLinearRefine>(d_dim));
   addRefineOperator(
      typeid(pdat::FaceVariable<float>).name(),
      boost::make_shared<CartesianFaceFloatConservativeLinearRefine>(d_dim));
   addRefineOperator(
      typeid(pdat::SideVariable<double>).name(),
      boost::make_shared<CartesianSideDoubleConservativeLinearRefine>(d_dim));
   addRefineOperator(
      typeid(pdat::SideVariable<float>).name(),
      boost::make_shared<CartesianSideFloatConservativeLinearRefine>(d_dim));
   addRefineOperator(
      typeid(pdat::CellVariable<dcomplex>).name(),
      boost::make_shared<CartesianCellComplexLinearRefine>(d_dim));
   addRefineOperator(
      typeid(pdat::CellVariable<double>).name(),
      boost::make_shared<CartesianCellDoubleLinearRefine>(d_dim));
   addRefineOperator(
      typeid(pdat::CellVariable<float>).name(),
      boost::make_shared<CartesianCellFloatLinearRefine>(d_dim));
   addRefineOperator(
      typeid(pdat::NodeVariable<dcomplex>).name(),
      boost::make_shared<CartesianNodeComplexLinearRefine>(d_dim));
   addRefineOperator(
      typeid(pdat::NodeVariable<double>).name(),
      boost::make_shared<CartesianNodeDoubleLinearRefine>(d_dim));
   addRefineOperator(
      typeid(pdat::NodeVariable<float>).name(),
      boost::make_shared<CartesianNodeFloatLinearRefine>(d_dim));
}

/*
 *************************************************************************
 *
 * Print CartesianGridGeometry class data.
 *
 *************************************************************************
 */

void
CartesianGridGeometry::printClassData(
   std::ostream& os) const
{
   const tbox::Dimension& dim(getDim());

   os << "Printing CartesianGridGeometry data: this = "
      << (CartesianGridGeometry *)this << std::endl;
   os << "d_object_name = " << getObjectName() << std::endl;

   int id;
   os << "d_x_lo = ";
   for (id = 0; id < dim.getValue(); id++) {
      os << d_x_lo[id] << "   ";
   }
   os << std::endl;
   os << "d_x_up = ";
   for (id = 0; id < dim.getValue(); id++) {
      os << d_x_up[id] << "   ";
   }
   os << std::endl;
   os << "d_dx = ";
   for (id = 0; id < dim.getValue(); id++) {
      os << d_dx[id] << "   ";
   }
   os << std::endl;

   os << "d_domain_box = " << d_domain_box << std::endl;

   GridGeometry::printClassData(os);
}

/*
 *************************************************************************
 *
 * Write class version number and object state to database.
 *
 *************************************************************************
 */

void
CartesianGridGeometry::putToDatabase(
   const boost::shared_ptr<tbox::Database>& db) const
{
   TBOX_ASSERT(db);

   const tbox::Dimension& dim(getDim());

   db->putInteger("GEOM_CARTESIAN_GRID_GEOMETRY_VERSION",
      GEOM_CARTESIAN_GRID_GEOMETRY_VERSION);
   tbox::Array<tbox::DatabaseBox> temp_box_array = getPhysicalDomain();
   db->putDatabaseBoxArray("d_physical_domain", temp_box_array);

   db->putDoubleArray("d_dx", d_dx, dim.getValue());
   db->putDoubleArray("d_x_lo", d_x_lo, dim.getValue());
   db->putDoubleArray("d_x_up", d_x_up, dim.getValue());

   hier::IntVector level0_shift(
      getPeriodicShift(hier::IntVector::getOne(dim)));
   int* temp_shift = &level0_shift[0];
   db->putIntegerArray("d_periodic_shift", temp_shift, dim.getValue());

}

/*
 *************************************************************************
 *
 * Data is read from input only if the simulation is not from restart.
 * Otherwise, all values specifed in the input database are ignored.
 * In this method data from the database are read to local
 * variables and the setGeometryData() method is called to
 * initialize the data members.
 *
 *************************************************************************
 */

void
CartesianGridGeometry::getFromInput(
   const boost::shared_ptr<tbox::Database>& db,
   bool is_from_restart)
{
   TBOX_ASSERT(db);

   const tbox::Dimension& dim(getDim());

   if (!is_from_restart) {

      hier::BoxContainer domain;
      if (db->keyExists("domain_boxes")) {
         hier::BoxContainer input_domain(
            db->getDatabaseBoxArray("domain_boxes"));
         if (input_domain.size() == 0) {
            TBOX_ERROR(
               "CartesianGridGeometry::getFromInput() error...\n"
               << "    geometry object with name = " << getObjectName()
               << "\n    Empty `domain_boxes' array found in input."
               << std::endl);
         }
         hier::LocalId local_id(0);
         for (hier::BoxContainer::iterator itr = input_domain.begin();
              itr != input_domain.end(); ++itr) {
            itr->setBlockId(hier::BlockId(0));
            domain.pushBack(hier::Box(*itr, local_id++, 0));
         }
      } else {
         TBOX_ERROR("CartesianGridGeometry::getFromInput() error...\n"
            << "    geometry object with name = " << getObjectName()
            << "\n    Key data `domain_boxes' not found in input." << std::endl);
      }

      double x_lo[tbox::Dimension::MAXIMUM_DIMENSION_VALUE],
             x_up[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

      if (db->keyExists("x_lo")) {
         db->getDoubleArray("x_lo", x_lo, dim.getValue());
      } else {
         TBOX_ERROR("CartesianGridGeometry::getFromInput() error...\n"
            << "    geometry object with name = " << getObjectName()
            << "\n    Key data `x_lo' not found in input." << std::endl);
      }
      if (db->keyExists("x_up")) {
         db->getDoubleArray("x_up", x_up, dim.getValue());
      } else {
         TBOX_ERROR("CartesianGridGeometry::getFromInput() error...\n"
            << "    geometry object with name = " << getObjectName()
            << "\n   Key data `x_up' not found in input." << std::endl);
      }

      int pbc[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
      hier::IntVector per_bc(dim, 0);
      if (db->keyExists("periodic_dimension")) {
         db->getIntegerArray("periodic_dimension", pbc, dim.getValue());
         for (int i = 0; i < dim.getValue(); i++) {
            per_bc(i) = ((pbc[i] == 0) ? 0 : 1);
         }
      }


      initializePeriodicShift(per_bc);

      setGeometryData(x_lo, x_up, domain);

   }
}

/*
 *************************************************************************
 *
 * Checks to see if the version number for the class is the same as
 * as the version number of the restart file.
 * If they are equal, then the data from the database are read to local
 * variables and the setGeometryData() method is called to
 * initialize the data members.
 *
 *************************************************************************
 */
void
CartesianGridGeometry::getFromRestart()
{
   const tbox::Dimension& dim(getDim());

   boost::shared_ptr<tbox::Database> restart_db(
      tbox::RestartManager::getManager()->getRootDatabase());

   if (!restart_db->isDatabase(getObjectName())) {
      TBOX_ERROR("CartesianGridGeometry::getFromRestart() error...\n"
         << "    database with name " << getObjectName()
         << " not found in the restart file" << std::endl);
   }
   boost::shared_ptr<tbox::Database> db(
      restart_db->getDatabase(getObjectName()));

   int ver = db->getInteger("GEOM_CARTESIAN_GRID_GEOMETRY_VERSION");
   if (ver != GEOM_CARTESIAN_GRID_GEOMETRY_VERSION) {
      TBOX_ERROR("CartesianGridGeometry::getFromRestart() error...\n"
         << "    geometry object with name = " << getObjectName()
         << "Restart file version is different than class version" << std::endl);
   }
   hier::BoxContainer restart_domain(
      db->getDatabaseBoxArray("d_physical_domain"));
   double x_lo[tbox::Dimension::MAXIMUM_DIMENSION_VALUE],
          x_up[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];
   db->getDoubleArray("d_x_lo", x_lo, dim.getValue());
   db->getDoubleArray("d_x_up", x_up, dim.getValue());

   hier::BoxContainer domain;
   hier::LocalId local_id(0);
   for (hier::BoxContainer::iterator itr = restart_domain.begin();
        itr != restart_domain.end(); ++itr) {
      itr->setBlockId(hier::BlockId(0));
      domain.pushBack(hier::Box(*itr, local_id++, 0));
   }

   setGeometryData(x_lo, x_up, domain);

   hier::IntVector periodic_shift(dim);
   int* temp_shift = &periodic_shift[0];
   db->getIntegerArray("d_periodic_shift", temp_shift, dim.getValue());
   initializePeriodicShift(periodic_shift);

}

}
}
#endif
