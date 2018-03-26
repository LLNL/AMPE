/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Simple Cartesian grid geometry for an AMR hierarchy.
 *
 ************************************************************************/

#ifndef included_geom_CartesianGridGeometry
#define included_geom_CartesianGridGeometry

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/geom/GridGeometry.h"
#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Serializable.h"

#include <boost/shared_ptr.hpp>
#include <string>

namespace SAMRAI {
namespace geom {

/**
 * Class CartesianGridGeometry provides simple Cartesian mesh geometry
 * management on an AMR hierarchy.  The mesh coordinates on each hierarchy
 * level are limited to mesh increments specified as DIM-tuple
 * (dx[0],...,dx[DIM-1]) and spatial coordinates of the lower and upper
 * corners of the smallest parallelepiped bounding the entire computational
 * domain.  The mesh increments on each level are defined with respect to
 * the coarsest hierarchy level and multiplying those values by the proper
 * refinement ratio.  This class sets geometry information on each patch in
 * an AMR hierarchy.  This class is derived from the geom::GridGeometry
 * base class.
 *
 * An object of this class requires numerous parameters to be read from
 * input.  Also, data must be written to and read from files for restart.
 * The input and restart data are summarized as follows:
 *
 * Required input keys and data types:
 *
 *
 *
 *
 *    - \b    domain_boxes
 *       tbox::Array of boxes representing the index space for the entire
 *       domain (on the coarsest refinement level).
 *
 *    - \b    x_lo
 *       tbox::Array of double values representing the spatial coordinates of
 *       the lower corner of the physical domain.
 *
 *    - \b    x_up
 *       tbox::Array of double values representing the spatial coordinates of
 *       the upper corner of the physical domain.
 *
 *
 *
 *
 *
 * Optional input keys, data types, and defaults:
 *
 *
 *
 *
 *    - \b    periodic_dimension
 *       tbox::Array of integer values representing the directions in which
 *       the physical domain is periodic.  A non-zero value indicates
 *       that the direction is periodic.  A zero value indicates that
 *       the direction is not periodic.  If no values are specified, then
 *       the array is initialized to all zeros (no periodic directions).
 *
 *
 *
 * No input values can overwrite restart values.
 *
 * A sample input file for a two-dimensional problem might look like:
 *
 * @verbatim
 *
 *    domain_boxes = [(0,0) , (49,39)]
 *    x_lo = 0.0 , 0.0
 *    x_up = 50.0 , 40.0
 *    periodic_dimension = 0, 1  // periodic in y only
 *
 * @endverbatim
 *
 * This generates a two-dimensional rectangular domain periodic in the
 * y-direction, and having 50 cells in the x-direction and 40 cells in
 * the y-direction, with the cell size 1 unit in each direction.
 *
 * @see geom::GridGeometry
 */

class CartesianGridGeometry:
   public geom::GridGeometry
{
   friend class TransferOperatorRegistry;

   typedef hier::PatchGeometry::TwoDimBool TwoDimBool;

public:
   /**
    * Constructor for CartesianGridGeometry initializes data
    * members based on parameters read from the specified input database
    * or from the restart database corresponding to the specified
    * object name.  The constructor also registers this object
    * for restart using the specified object name when the boolean
    * argument is true.  Whether object will write its state to restart
    * files during program execution is determined by this argument.
    * Note that it has a default state of true.
    *
    * Errors: passing in a null database pointer or an empty std::string
    * will result in an unrecoverable assertion.
    */
   CartesianGridGeometry(
      const tbox::Dimension& dim,
      const std::string& object_name,
      const boost::shared_ptr<tbox::Database>& input_db,
      bool register_for_restart = true);

   /**
    * Constructor for CartesianGridGeometry sets data members
    * based on arguments.  The constructor also registers this object
    * for restart using the specified object name when the boolean
    * argument is true.  Whether object will write its state to restart
    * files during program execution is determined by this argument.
    * Note that it has a default state of true.
    *
    * Errors: passing in an empty std::string, or null data pointers will
    * result in an unrecoverable assertion.
    */
   CartesianGridGeometry(
      const std::string& object_name,
      const double* x_lo,
      const double* x_up,
      const hier::BoxContainer& domain,
      bool register_for_restart = true);

   /*!
    * @brief Construct a new coarsened/refined CartesianGridGeometry object
    * with the supplied domain.
    *
    * This method is intended to be called only by boost::make_shared from the
    * make[Coarsened, Refined]GridGeometry methods to make a coarsened or
    * refined version of a given CartesianGridGeometry.
    *
    * @param[in] object_name The same name as the uncoarsened/unrefined grid
    *            geometry.
    * @param[in] x_lo The same lower corner as the uncoarsened/unrefined grid
    *            geometry.
    * @param[in] x_up The same upper corner as the uncoarsened/unrefined grid
    *            geometry.
    * @param[in] domain The coarsened/refined domain.
    * @param[in] op_reg The same operator registry as the uncoarsened/unrefined
    *            grid geometry.
    * @param[in] register_for_restart Flag indicating whether this instance
    *            should be registered for restart.
    */
   CartesianGridGeometry(
      const std::string& object_name,
      const double* x_lo,
      const double* x_up,
      const hier::BoxContainer& domain,
      const boost::shared_ptr<hier::TransferOperatorRegistry>& op_reg,
      bool register_for_restart);

   /**
    * Destructor for CartesianGridGeometry deallocates
    * data describing grid geometry and unregisters the object with
    * the restart manager if previously registered.
    */
   virtual ~CartesianGridGeometry();

   /**
    * Create and return a pointer to a refined version of this Cartesian grid
    * geometry object.
    */
   boost::shared_ptr<hier::BaseGridGeometry>
   makeRefinedGridGeometry(
      const std::string& fine_geom_name,
      const hier::IntVector& refine_ratio,
      bool register_for_restart) const;

   /**
    * Create and return a pointer to a coarsened version of this Cartesian grid
    * geometry object.
    */
   boost::shared_ptr<hier::BaseGridGeometry>
   makeCoarsenedGridGeometry(
      const std::string& coarse_geom_name,
      const hier::IntVector& coarsen_ratio,
      bool register_for_restart) const;

   /*
    * Compute grid data for patch and assign new geom_CartesianPatchGeometry
    * object to patch.
    */
   void
   setGeometryDataOnPatch(
      hier::Patch& patch,
      const hier::IntVector& ratio_to_level_zero,
      const TwoDimBool& touches_regular_bdry,
      const TwoDimBool& touches_periodic_bdry) const;

   /**
    * Set data members for this CartesianGridGeometry object.
    */
   void
   setGeometryData(
      const double* x_lo,
      const double* x_up,
      const hier::BoxContainer& domain);

   /**
    * Return const pointer to dx array for reference level in hierarchy.
    */
   const double *
   getDx() const
   {
      return d_dx;
   }

   /**
    * Return const pointer to lower spatial coordinate for reference
    * level in hierarchy.
    */
   const double *
   getXLower() const
   {
      return d_x_lo;
   }

   /**
    * Return const pointer to upper spatial coordinate for reference
    * level in hierarchy.
    */
   const double *
   getXUpper() const
   {
      return d_x_up;
   }

   /**
    * Print class data representation.
    */
   virtual void
   printClassData(
      std::ostream& os) const;

   /**
    * Writes the state of the CartesianGridGeometry object to the
    * database.
    *
    * When assertion checking is active, db cannot be a null database pointer.
    */
   virtual void
   putToDatabase(
      const boost::shared_ptr<tbox::Database>& db) const;

protected:
   /*!
    * @brief Build operators appropriate for a CartesianGridGeometry.
    */
   virtual void
   buildOperators();

private:
   /*
    * Static integer constant describing class's version number.
    */
   static const int GEOM_CARTESIAN_GRID_GEOMETRY_VERSION;

   /*
    * Reads in domain_boxes, x_lo, and x_up from the input database.
    * Data is read from input only if the simulation is not from restart.
    * Otherwise, all values specified in the input database are ignored.
    *
    * Arguments: is_from_restart is true when simulation is from restart
    * Assertions: db must not be a null pointer.
    */
   void
   getFromInput(
      const boost::shared_ptr<tbox::Database>& db,
      bool is_from_restart);

   /*
    * Read object state from the restart file and initialize class data
    * members.  The database from which the restart data is read is
    * determined by the object_name specified in the constructor.
    *
    * Unrecoverable Errors:
    *
    *    -The database corresponding to object_name is not found
    *     in the restart file.
    *
    *    -The class version number and restart version number do not
    *     match.
    *
    */
   void
   getFromRestart();

   /*
    * Flag to determine whether this instance is registered for restart.
    */
   bool d_registered_for_restart;

   double d_dx[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];                // mesh increments for level 0.
   double d_x_lo[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];              // spatial coordinates of lower corner
   // (i.e., box corner) of problem domain.
   double d_x_up[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];              // spatial coordinates of upper corner
   // (i.e., box corner) of problem domain.

   hier::Box d_domain_box;           // smallest box covering coarsest level
                                     // (i.e., reference level) index space.

};

}
}

#endif
