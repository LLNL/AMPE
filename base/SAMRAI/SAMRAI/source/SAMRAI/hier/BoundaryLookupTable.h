/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Lookup table to aid in BoundaryBox construction
 *
 ************************************************************************/

#ifndef included_hier_BoundaryLookupTable
#define included_hier_BoundaryLookupTable

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/Array.h"

namespace SAMRAI {
namespace hier {

/*!
 * @brief Singleton class to organize patch boundary information.
 *
 * Class BoundaryLookupTable is a singleton class that maintains
 * a table that organizes all of the possible boundary region cases
 * for a patch.  It is used primarily by the grid geometry during the
 * construction of physical boundary boxes for patches and by the
 * PatchGeometry class to determine box regions to be filled during
 * a physical boundary fill.
 *
 * This class is useful for any situation where enumerating the
 * cases for boundary regions around a box is needed. The main advantage
 * of using this class is that such calculations can be programmed in
 * a dimension-independent way.
 *
 * @see hier::BoundaryBox
 * @see hier::BaseGridGeometry
 * @see hier::PatchGeometry
 */

class BoundaryLookupTable
{
public:
   /*!
    * @brief Return pointer to singleton instance of the boundary
    * lookup table.
    *
    * Following the Singleton design pattern, users of this class
    * do not explicitly allocate or deallocate the Singleton instance.
    *
    * @param   dim   Dimension of the object.
    *
    * @return  Pointer to lookup table instance.
    */
   static BoundaryLookupTable *
   getLookupTable(
      const tbox::Dimension& dim)
   {
     int idx = dim.getValue() - 1;
      if (!s_lookup_table_instance[idx]) {
         s_lookup_table_instance[idx] = new BoundaryLookupTable(dim);
      }
      return s_lookup_table_instance[idx];
   }

   /*!
    * @brief Get array of active directions for specific boundary
    * location and codimension case.
    *
    * The active directions are those coordinate direction in which
    * the boundary region would have to be shifted in order to be
    * contained in the corresponding patch box (whose boundary
    * we are interested in).
    *
    * @param loc   integer location index of boundary region
    * @param codim integer codimension of boundary region
    *
    * @return  const reference to integer array of length codim
    *          containing the active directions for this boundary case.
    */
   const tbox::Array<int>&
   getDirections(
      int loc,
      int codim) const
   {
      TBOX_ASSERT((codim > 0) && (codim <= d_dim.getValue()));
      TBOX_ASSERT((loc >= 0) && (loc < d_max_li[codim - 1]));
      int iloc = loc / (1 << codim);
      return d_table[codim - 1][iloc];
   }

   /*!
    * @brief Get array of maximum number of locations for each
    * codimension boundary case.
    *
    * For example, a 2D patch has 4 possible locations for codimension 1
    * (edges) and 4 possible locations for codimension 2 (nodes).  A 3D
    * patch has 6 possible locations for codimension 1 (faces), 12 for
    * codimension 2 (edges), and 8 for codimension 3 (nodes).
    *
    * @return integer array of length dim, each entry of which indicates
    *         the maximum number of boundary locations for each
    *         codimension
    */
   const tbox::Array<int>&
   getMaxLocationIndices() const
   {
      return d_max_li;
   }

   /*!
    * @brief Determines if given boundary information indicates a
    * a lower boundary region.
    *
    * A boundary region is a "lower" boundary region if the associated patch
    * box contains higher values along the axis in the coordinate
    * direction than the boundary region.
    *
    * @param loc       integer location index of boundary region
    * @param codim     integer codimension of boundary region
    * @param dir_index integer spatial direction identifier, an index into the
    *                  array returned by getDirections().
    *
    * @return bool true if the boundary type of codimension codim indexed
    * by loc is a lower boundary in the specified direction;
    * return false if the boundary is an upper boundary.
    */
   bool
   isLower(
      int loc,
      int codim,
      int dir_index) const
   {
      TBOX_ASSERT((codim > 0) && (codim <= d_dim.getValue()));
      TBOX_ASSERT((loc >= 0) && (loc < d_max_li[codim - 1]));
      TBOX_ASSERT((dir_index >= 0) && (dir_index < codim));
      return !isUpper(loc, codim, dir_index);
   }

   /*!
    * @brief Determines if given boundary information indicates a
    * an upper boundary region.
    *
    * A boundary region is an "upper" boundary region if the associated patch
    * box contains lower values along the axis in the coordinate
    * direction than the boundary region.
    *
    * @param loc       integer location index of boundary region
    * @param codim     integer codimension of boundary region
    * @param dir_index integer spatial direction identifier, an index into the
    *                  array returned by getDirections().
    *
    * @return bool true if the boundary type of codimension codim indexed
    * by loc is an upper boundary in the specified direction;
    * return false if the boundary is a lower boundary.
    */
   bool
   isUpper(
      int loc,
      int codim,
      int dir_index) const
   {
      TBOX_ASSERT((codim > 0) && (codim <= d_dim.getValue()));
      TBOX_ASSERT((loc >= 0) && (loc < d_max_li[codim - 1]));
      TBOX_ASSERT((dir_index >= 0) && (dir_index < codim));
      return (loc % (1 << codim)) & (1 << (dir_index));
   }

   /*!
    * @brief Get array of boundary direction IntVectors.
    *
    * For any codimension, there is a particular number of valid boundary
    * locations.  This function returns an array of IntVectors that provide
    * information about where each boundary location lies in relation to
    * a patch.  The array's length is the number of valid locations for
    * the given codimension, and the array is indexed by the location id's
    * that are set up by this BoundaryLookupTable class.
    *
    * For a particular location, each element of the IntVector tells whether
    * the location is on the lower or upper side, or neither, of the patch in
    * a specific coordinate direction.  -1 indicates lower, 0 indicates
    * neither, and 1 indicates upper.
    *
    * @param codim  codimension
    * @return       Array of IntVectors, one element for each valid location
    */
   const tbox::Array<IntVector>&
   getBoundaryDirections(
      int codim) const
   {
      TBOX_ASSERT((codim > 0) && (codim <= d_dim.getValue()));
      return d_bdry_dirs[codim - 1];
   }

protected:
   /*!
    * @brief protected constructor
    *
    * The constructor for BoundaryLookupTable is protected.
    * Consistent with the definition of a Singleton class, only the
    * lookup table has access to the constructor for the class.
    *
    * The constructor initializes the state of lookup table contents.
    *
    * @param dim  Dimension
    */
   explicit BoundaryLookupTable(
      const tbox::Dimension& dim);

   /*!
    * @brief protected destructor
    *
    * The destructor for BoundaryLookupTable is protected. See the
    * comments for the constructor.
    *
    * The destructor deallocates lookup table contents.
    */
   ~BoundaryLookupTable();

private:
   /*!
    * @brief Build table by recursively computing the entries in the
    * lookup table for a given codimension.
    *
    * TODO:  Document the parameters.
    */
   void
   buildTable(int* table,
              int codim,
              int ibeg,
              int(&work)[tbox::Dimension::MAXIMUM_DIMENSION_VALUE],
              int& lvl,
              int * & ptr);

   /*!
    * @brief Build table of direction IntVectors
    *
    * Private member function that builds table of intvectors that stores
    * whether a boundary location is upper, lower, or neither in each
    * coordinate direction
    */
   void
   buildBoundaryDirectionVectors();

   /*!
    * @brief Deallocate the BoundaryLookupTable instance.
    *
    * It is not necessary to call this function at program termination,
    * since it is automatically called by the StartupShutdownManager class.
    */
   static void
   finalizeCallback();

   /*!
    * @brief Static data members used to control access to and destruction of
    * singleton variable database instance.
    */
   static BoundaryLookupTable* s_lookup_table_instance[tbox::Dimension::
                                                       MAXIMUM_DIMENSION_VALUE];

   /*!
    * @brief Dimension of the object
    */
   const tbox::Dimension d_dim;

   /*!
    * @brief Array used to store the number of combinations for
    * each codimension.
    */
   tbox::Array<int> d_ncomb;

   /*!
    * @brief Array used to store the number of possible location indices
    * for each codimension.
    */
   tbox::Array<int> d_max_li;

   /*!
    * @brief Data member used to store the lookup table.
    */
   tbox::Array<tbox::Array<int> >
   d_table[tbox::Dimension::MAXIMUM_DIMENSION_VALUE];

   /*!
    * @brief Array to hold information about possible directions for each
    * codimension.
    */
   tbox::Array<tbox::Array<IntVector> > d_bdry_dirs;

   static tbox::StartupShutdownManager::Handler
      s_finalize_handler;

};

}
}

#endif
