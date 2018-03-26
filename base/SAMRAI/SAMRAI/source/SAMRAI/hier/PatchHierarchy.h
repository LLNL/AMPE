/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   An AMR hierarchy of patch levels
 *
 ************************************************************************/

#ifndef included_hier_PatchHierarchy
#define included_hier_PatchHierarchy

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/BaseGridGeometry.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/PatchDescriptor.h"
#include "SAMRAI/hier/PatchFactory.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/PatchLevelFactory.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Serializable.h"
#include "SAMRAI/tbox/Utilities.h"

#include <boost/shared_ptr.hpp>
#include <string>

namespace SAMRAI {
namespace hier {

/*!
 * @brief Class PatchHierarchy maintains the patch levels that
 * define the AMR hierarchy.
 *
 * The following describes the input file keys and data types. For each
 * input key that involves input for levels in the hierarchy, assume that
 * we have N levels, numbered coarsest to finest as 0,..., N-1.  For such
 * input, the value for a level must be given as ``level_n = value'' where
 * n is the level number.  When more values are given than needed for the
 * maximum number of levels, extra values are ignored.  When fewer values
 * are given, the last value provided will be used on each level without a
 * specified input value.  See example input below.
 *
 *   - \b    max_levels
 *      Integer value specifying maximum number of levels
 *      allowed in the AMR patch hierarchy.
 *
 *   - \b    ratio_to_coarser
 *      A set of max_levels - 1 integer vectors (each has length = DIM),
 *      each of which indicates the ratio of the index space of a patch
 *      level to that of the next coarser level in the hierarchy.  The
 *      input is given for each level n, where n (= 1, 2,..., N-1) is the
 *      level number.
 *
 *   - \b    largest_patch_size
 *      A set of max_levels integer vectors (each has length = DIM) each of
 *      which indicates the dimensions of largest patch allowed on the level
 *      in the hierarchy.  Non-positive values for patch size corresponds to
 *      no upper limit on patch size.  The input is given for each level n,
 *      where n (= 0, 1,..., N-1) is the level number.
 *
 *   - \b    smallest_patch_size
 *      A set of max_levels integer vectors (each has length = DIM) each of
 *      which indicates the dimensions of smallest patch allowed on the level
 *      in the hierarchy.  The smallest patch allowed must be at least as
 *      large as the maximum ghost cell width for all variables in the problem.
 *      If some smaller patch size is given in input, then it will be
 *      overridden by the maximum ghost width.  If no input is given, a
 *      default of the maximum ghost cell width over all variables is used.
 *      The input is given for each level n, where n (= 0, 1,..., N-1) is
 *      the level number.
 *
 *   - \b    proper_nesting_buffer
 *      A set of max_levels - 1 integer values specifying the number of
 *      coarse cells by which the next finer level is nested within the
 *      interior of the union of the patches on the next coarser level.
 *      The input is given for each level n, where n (= 0, 1,..., N-2) is
 *      the level number.
 *
 *    - \b    allow_patches_smaller_than_ghostwidth
 *      A boolean value ("TRUE" or "FALSE") indicating whether patches
 *      are allowed that are smaller than the maximum variable ghost width
 *      along some coordinate direction.  Recall that when a smallest
 *      patch size provided in the input file is smaller than the maximum
 *      ghost width of all the registered variables, then by default the
 *      smallest patch size will be set to the maximum ghost width.  Set this
 *      flag to TRUE to override this default behavior and to allow the
 *      smallest patch size given in the input to remain in effect.
 *      The default value is FALSE.
 *
 *    - \b    allow_patches_smaller_than_minimum_size_to_prevent_overlaps
 *      A boolean value ("TRUE" or "FALSE") indicating whether patches
 *      are allowed to be smaller than the minimum patch size to prevent
 *      overlapping patches.  In order to enforce minimum patch size
 *      restrictions, boxes may be grown during adaptive gridding operations.
 *      This may lead to patches whose boxes overlap.  This may be a problem
 *      for some applications. If overlaps are undesirable and you are
 *      willing to relax the minimum size constraints, set this parameter
 *      TRUE.  By default, it is FALSE.
 *
 * Note that when continuing from restart, the input values in the
 * input file override all values read in from the restart database.
 *
 * The following represents sample input data for a three-dimensional problem:
 *
 * @code
 *
 *   // Required input: maximum number of levels in patch hierarchy
 *   max_levels = 4
 *
 *   // Required input: vector ratio between each finer level and next coarser
 *   ratio_to_coarser {
 *      level_1 = 2, 2, 2
 *      level_2 = 2, 2, 2
 *      level_3 = 4, 4, 4
 *   }
 *
 *   // Optional input: int vector for largest patch size on each level.
 *   largest_patch_size {
 *      level_0 = 40, 40, 40
 *      level_1 = 30, 30, 30
 *      // all finer levels will use same values as level_1...
 *   }
 *
 *   // Optional input: int vector for smallest patch size on each level.
 *   smallest_patch_size {
 *      level_0 = 16, 16, 16
 *      // all finer levels will use same values as level_0...
 *   }
 *
 *   // Optional input:  buffer of one cell used on each level
 *   proper_nesting_buffer = 1
 *
 * @endcode
 *
 * @see hier::PatchLevel
 * @see hier::PatchDescriptor
 */

class PatchHierarchy:public tbox::Serializable
{
public:
   /*
    *  TODO: There must be a better way to do what this base class
    *  provides.  If not, then we need to make this documentation clearer.
    *  Specifically, it is not clear how one actually determines the
    *  correct width information to provide.  Also, the documentation
    *  alludes to usage that is internal to SAMRAI (e.g., the stuff about
    *  GriddingAlgorithm) so it is not clear what a user needs to do
    *  and what she does not.
    *
    *  The notion of "connector width" does not appear to be defined
    *  anywhere.  If is is, where is the definition?
    */
   /*!
    * @brief Abstract base class defining a Strategy pattern interface
    * for providing Connector width information to the PatchHierarchy.
    *
    * The hierarchy needs to know that Connector width constraints
    * must be applied so it can provide correct connector width information
    * to mesh generation operations so that they can generate Connectors
    * with sufficient widths.  It is most efficient (scales better) to
    * generate overlap Connectors when the mesh is built rather than searching
    * for overlaps later.  If the required widths are not set up this way,
    * or if the code building the hierarchy does not provide Connectors
    * of sufficient widths, SAMRAI will still work but will not scale
    * well.  The mechanism for generating missing Connectors as needed
    * is in the class PersistentOverlapConnector.
    *
    * The user of a PatchHierarchy object registers implementations of
    * this strategy class with the PatchHierarchy.  See
    * PatchHierarchy::registerConnectorWidthRequestor().
    * By default, requesters implemented by the SAMRAI library will be
    * automatically registered. Unless you want to provide your own
    * requesters, no special action is required.
    *
    * TODO: The following paragraph is more confusing than helpful.
    * It exposes internal SAMRAI implementation details and confuses
    * what a user needs to do.
    *
    * When the hierarchy is first asked for any required width, it
    * uses the registered implementations to compute all required
    * widths.  After that, no further registration is allowed.
    * Implementations of this class usually requires the max ghost
    * data width and the max stencil widths, so it is best to make
    * sure that those have been registered before doing something that
    * calls getRequiredConnectorWidth(), such as using
    * mesh::GriddingAlgorithm to populate the hierarchy with a level.
    *
    * Example of a typical scenario: GriddingAlgorithm is used to
    * build the hierarchy.  Two classes that will want the hierarchy's
    * Connectors to have some specified widths are GriddingAlgorithm
    * and RefineSchedule.  This strategy should therefore be
    * implemented by GriddingAlgorithm and RefineSchedule.
    */
   class ConnectorWidthRequestorStrategy
   {
public:
      /*
       * TODO: As a general rule in SAMRAI, we provide a default ctor
       * explicitly (at least).
       */

      /*!
       * @brief Destructor
       */
      virtual ~ConnectorWidthRequestorStrategy();

      /*
       * TODO: How is a developer supposed to know what IntVector values
       * are needed when she implements this routine?
       */
      /*!
       * @brief Provide Connector widths the child class requires in
       * order to work properly on a given hierarchy.
       *
       * The Connector widths are computed to the maximum that could be
       * requested between levels in the given PatchHierarchy. The two vector
       * parameters will contain the computed widths.
       * @par Assumptions
       * <ul>
       * <li> On completion of the function call, self_connector_width must be
       *      of length @c patch_hierarchy.getMaxNumberOfLevels().
       * <li> On completion of the function call, fine_connector_width must be
       *      of length @c patch_hierarchy.getMaxNumberOfLevels() - 1.
       * </ul>
       *
       * @param[out] self_connector_widths Array of widths for Connectors
       * from a level to itself.
       *
       * @param[out] fine_connector_widths Array of widths for Connectors
       * from a level to the next finer level.
       *
       * @param[in]  patch_hierarchy
       */
      virtual void
      computeRequiredConnectorWidths(
         std::vector<IntVector>& self_connector_widths,
         std::vector<IntVector>& fine_connector_widths,
         const PatchHierarchy& patch_hierarchy) const = 0;
   };

/*
 * TODO: How does the hierarchy get properly initialized if the input
 * database argument is a null pointer?
 */
/*!
 * @brief Constructor initializing the state of PatchHierarchy.
 *
 * The constructor for the PatchHierarchy initializes the number of
 * levels to zero, sets the geometry for the PatchHierarchy, and
 * registers the PatchHierarchy for restart with the specified name
 * by default.
 *
 * @par Errors/Assertions
 * Passing in an empty std::string or a null grid geometry pointer
 * will result in an unrecoverable assertion when assertion checking is
 * active.
 *
 * @param[in]  object_name
 * @param[in]  geometry
 * @param[in]  database Database specifying hierarchy parameters.
 * @param[in]  register_for_restart @b Default: true
 */
   PatchHierarchy(
      const std::string& object_name,
      const boost::shared_ptr<BaseGridGeometry>& geometry,
      const boost::shared_ptr<tbox::Database>& database =
         boost::shared_ptr<tbox::Database>(),
      bool register_for_restart = true);

   /*!
    * @brief Destructor for PatchHierarchy.
    */
   ~PatchHierarchy();

   /*!
    * @brief Create a refined version of this patch hierarchy.
    *
    * Create and return pointer to a patch hierarchy that is a refined
    * version of this patch hierarchy object.  That is, the data
    * members of the returned patch hierarchy are set by refining
    * information on this hierarchy using the given ratio.  The refined
    * hierarchy will cover the same physical space as this hierarchy
    * and will have the same number of levels and same mapping of
    * patches to processors on each level.  However, the index space
    * of each level will be refined by the specified ratio.
    * @note
    * This function does not allocate patch data so this must be done
    * before any data operations can be performed on the new
    * hierarchy.
    *
    * @return The refined patch hierarchy.
    *
    * @param[in]  fine_hierarchy_name
    * @param[in]  refine_ratio
    * @param[in]  register_for_restart
    */
   boost::shared_ptr<PatchHierarchy>
   makeRefinedPatchHierarchy(
      const std::string& fine_hierarchy_name,
      const IntVector& refine_ratio,
      bool register_for_restart) const;

   /*!
    * @brief Create a coarsened version of this patch hierarchy.
    *
    * Create and return pointer to a patch hierarchy that is a
    * coarsened version of this patch hierarchy object.  That is, the
    * data members of the returned patch hierarchy are set by
    * coarsening information on this hierarchy by the given ratio.
    * The coarsened hierarchy will cover the same physical space as
    * this hierarchy and will have the same number of levels and same
    * mapping of patches to processors on each level.  However, the
    * index space of each level will be coarsened by the specified
    * ratio.
    * @note
    * This function does not allocate patch data so
    * this must be done before any data operations can be performed on
    * the new hierarchy.
    *
    * @return boost::shared_ptr to the coarsened patch hierarchy.
    *
    * @param[in]  coarse_hierarchy_name
    * @param[in]  coarsen_ratio
    * @param[in]  register_for_restart
    */
   boost::shared_ptr<PatchHierarchy>
   makeCoarsenedPatchHierarchy(
      const std::string& coarse_hierarchy_name,
      const IntVector& coarsen_ratio,
      bool register_for_restart) const;

/*
 * TODO: Is it an error to call this method when a level with the given
 * level number already exists?  Are some preconditions assumed?
 */
/*!
 * @brief Construct new PatchLevel in hierarchy at given level number.
 *
 * Boxes, their mappings and the refinement ratio are obtained from
 * @c new_mapped_box_level.
 *
 * @param[in]  level_number
 * @param[in]  new_mapped_box_level
 */
   void
   makeNewPatchLevel(
      const int level_number,
      const BoxLevel& new_mapped_box_level);

   /*!
    * @brief Remove a patch level
    *
    * Remove PatchLevel from the hierarchy and adjust number of levels
    * accordingly.
    *
    * @param[in]  level
    */
   void
   removePatchLevel(
      const int level);

   /*!
    * @brief Get a PatchLevel from the hierarchy.
    *
    * @return a pointer to the specified patch level.
    *
    * @param[in]  level
    */
   boost::shared_ptr<PatchLevel>
   getPatchLevel(
      const int level) const
   {
      TBOX_ASSERT((level >= 0) && (level < d_number_levels));
      return d_patch_levels[level];
   }

   /*!
    * @brief Get the patch descriptor.
    *
    * @return a pointer to the patch descriptor used for the patches in
    * the patch hierarchy.
    */
   boost::shared_ptr<PatchDescriptor>
   getPatchDescriptor() const
   {
      return d_patch_descriptor;
   }

   /*!
    * @brief Check if the level exists
    *
    * @return True if the hierarchy contains a patch level with the given
    * level number; otherwise false.
    *
    * @param[in]  level
    */
   bool
   levelExists(
      const int level) const
   {
      return (level < d_number_levels) && d_patch_levels[level];
   }

   /*!
    * @brief Check if a finer level exists.
    *
    * @return True if the hierarchy contains a patch level with level
    * number level + 1 (i.e., finer than level with given number);
    * otherwise, false.
    *
    * @param[in]  level
    */
   bool
   finerLevelExists(
      const int level) const
   {
      return (level + 1 < d_number_levels) && d_patch_levels[level + 1];
   }

   /*!
    * @brief Get the number of levels in the hierarchy.
    *
    * @return The number of levels that currently exist in the hierarchy.
    */
   int
   getNumberOfLevels() const
   {
      return d_number_levels;
   }

   /*!
    * @brief Get the finest level in the hierarchy.
    *
    * @return The level number of the finest resolution patch level currently
    * existing in the hierarchy.
    */
   int
   getFinestLevelNumber() const
   {
      return d_number_levels - 1;
   }

   /*!
    * @brief Check whether specified level can be refined.
    *
    * @return true if level associated with the specified level number can
    * be refined; i.e., the level number is less than that of the finest
    * level allowed in the hierarchy.  Otherwise, false is returned.
    */
   bool
   levelCanBeRefined(
      const int level_number) const
   {
      TBOX_ASSERT(level_number >= 0);
      return level_number < getMaxNumberOfLevels() - 1;
   }

   /*!
    * @brief Return a pointer to the specified BoxLevel.
    *
    * @return The BoxLevel owned by PatchLevel number @c level.
    *
    * @param[in]  level
    */
   const boost::shared_ptr<BoxLevel>&
   getBoxLevel(
      const int level) const
   {
      return d_patch_levels[level]->getBoxLevel();
   }

   /*!
    * @brief Get the connector between two levels
    *
    * Get const access to the Connector between two given BoxLevels
    * between two given levels in the hierarchy.
    *
    * @return Connector between the two given level numbers.
    *
    * @param[in]  base_ln The base level indicating one end of the
    *             connector.
    * @param[in]  head_ln The head level indicating the other end of
    *             the connector.
    */
   const Connector&
   getConnector(
      const int base_ln,
      const int head_ln) const;

   //@{

   //! @name Connector width coordination between hierarchy builders and users.

   /*!
    * @brief Register the ConnectorWidthRequestorStrategy with the hierarchy.
    *
    * Used to adjust the results of getRequiredConnectorWidth().
    * This method can be called multiple times to provide as much
    * connector width information for any object that requires a connector
    * width constraint.  The largest Connector width over all provided
    * will be used in the hierarchy.
    *
    * @par Important
    * All registrations must occur before the first call to
    * getRequiredConnectorWidth().
    *
    * @par Errors/Assertions
    * Calling this method after the first getRequiredConnectorWidth()
    * results in an unrecoverable assertion.
    *
    * @param[in]  cwrs The connector width requester strategy instance.
    */
   void
   registerConnectorWidthRequestor(
      const ConnectorWidthRequestorStrategy& cwrs);

   /*!
    * @brief Get the width required of the Connector from the base to
    * the head levels, as specified by the given level numbers.
    *
    * The two level numbers must differ by no more than one.
    *
    * @return The widths required from the base to head levels.
    *
    * @param[in]  base_ln
    * @param[in]  head_ln
    */
   IntVector
   getRequiredConnectorWidth(
      int base_ln,
      int head_ln) const;

   /*!
    * @brief Add a ConnectorWidthRequestorStrategy implementation to
    * the set of implementations that, by default, gets registered to
    * all PatchHierarchy objects.
    *
    * @param[in] cwrs   ConnectorWidthRequestorStrategy to be registered.
    */
   static void
   registerAutoConnectorWidthRequestorStrategy(
      const ConnectorWidthRequestorStrategy& cwrs);

   //@}

/*
 * TODO: This DomainBoxLevel, etc. stuff (if it is really needed) should
 * be moved to the BaseGridGeometry class.  It makes the role of this class
 * too divergent by having it here.
 */
/*!
 * @brief Access the domain description as a BoxLevel.
 *
 * The domain BoxLevel is maintained in Globalized mode with
 * processor 0 owning all mapped boxes.
 *
 * @return The domain description as a BoxLevel
 */
   const BoxLevel&
   getDomainBoxLevel() const
   {
      return d_domain_mapped_box_level;
   }

   /*!
    * @brief Returns the SAMRAI_MPI communicator over which the domain
    * BoxLevel is distributed.
    */
   const tbox::SAMRAI_MPI&
   getMPI() const
   {
      return d_domain_mapped_box_level.getMPI();
   }

   //@{

   //! @name Methods related to the hierarchy parameters

   /*!
    * @brief Set the maximum number of levels allowed on the hierarchy.
    *
    * This method can only be called when the hierarchy has no PatchLevels.
    *
    * If this increases the max number of levels, the missing smallest
    * patch size, largest patch size and ratio to coarser level are
    * taken from the values for the current finest level.  These values can
    * be changed afterward using set methods below.
    *
    * @param[in]  max_levels
    */
   void
   setMaxNumberOfLevels(
      int max_levels)
   {
      d_max_levels = max_levels;
      if (d_max_levels != int(d_ratio_to_coarser.size())) {
         d_ratio_to_coarser.resize(d_max_levels, d_ratio_to_coarser.back());
         d_smallest_patch_size.resize(
            d_max_levels,
            d_smallest_patch_size.back());
         d_largest_patch_size.resize(
            d_max_levels,
            d_largest_patch_size.back());
         d_proper_nesting_buffer.resize(
            d_max_levels - 1,
            d_proper_nesting_buffer.empty() ?
               1 : d_proper_nesting_buffer.back());
      }
   }

   /*!
    * @brief Get the maximum number of levels allowed on the hierarchy.
    *
    * This number is set by setMaxNumberOfLevels() or in the input database.
    *
    * @return The maximum number of levels allowed on the hierarchy.
    */
   int
   getMaxNumberOfLevels() const
   {
      return d_max_levels;
   }

   /*!
    * @brief Set the ratio to coarser level.
    *
    * @param[in]  ratio   Refinement ratio in each direction
    * @param[in]  level
    */
   void
   setRatioToCoarserLevel(
      const IntVector& ratio,
      int level)
   {
      TBOX_ASSERT(level > 0 && level < d_max_levels);
      d_ratio_to_coarser[level] = ratio;
   }

   /*!
    * @brief Get the ratio between specified level and next coarser level
    *
    * @return The ratio between specified @c level and @ level-1
    *
    * @param[in]  level
    */
   const IntVector&
   getRatioToCoarserLevel(
      int level) const
   {
      TBOX_ASSERT(level < d_max_levels);
      return d_ratio_to_coarser[level];
   }

   /*!
    * @brief Set the smallest patch size on the given level.
    *
    * @param[in]  size   Smallest size in each direction
    * @param[in]  level
    */
   void
   setSmallestPatchSize(
      const IntVector& size,
      int level)
   {
      TBOX_ASSERT(level >= 0 && level < d_max_levels);
      d_smallest_patch_size[level] = size;
   }

   /*!
    * @brief Get the smallest patch size on the given level.
    *
    * @return The smallest patch size allowed on the given level.
    *
    * @param[in]  level
    */
   const IntVector&
   getSmallestPatchSize(
      int level) const
   {
      TBOX_ASSERT(level >= 0 && level < d_max_levels);
      return d_smallest_patch_size[level];
   }

   /*!
    * @brief Set the largest patch size on the given level.
    *
    * @param[in]  size   Largest size in each direction.
    * @param[in]  level
    */
   void
   setLargestPatchSize(
      const IntVector& size,
      int level)
   {
      TBOX_ASSERT(level >= 0 && level < d_max_levels);
      d_largest_patch_size[level] = size;
   }

   /*!
    * @brief Get the largest patch size.
    *
    * @return The largest patch size allowed on the given level.
    *
    * @param[in]  level
    */
   const IntVector&
   getLargestPatchSize(
      int level) const
   {
      TBOX_ASSERT(level >= 0 && level < d_max_levels);
      return d_largest_patch_size[level];
   }

   /*!
    * @brief Get the proper nesting buffer for a specific level.
    *
    * The proper nesting requirement is that level ln+1 must nest in
    * level ln by getProperNestingBuffer(ln) cells.  The nesting
    * requirement does not apply where level ln+1 touches a physical
    * boundary.
    *
    * If @c ln >= getMaxNumberOfLevels()-1, the return value will be -1.
    *
    * @return The proper nesting buffer
    *
    * @param[in]  ln
    */
   int
   getProperNestingBuffer(
      int ln) const
   {
      TBOX_ASSERT(ln >= 0 && ln < d_max_levels);
      return (ln < d_max_levels - 1) ? d_proper_nesting_buffer[ln] : -1;
   }

   /*!
    * @brief Get flag for allowing patches smaller than ghost width.
    */
   bool
   allowPatchesSmallerThanGhostWidth() const
   {
      return d_allow_patches_smaller_than_ghostwidth;
   }

   /*!
    * @brief Get flag for allowing patches smaller than user-provided minimum
    * size.
    */
   bool
   allowPatchesSmallerThanMinimumSize() const
   {
      return d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps;
   }

   //@}

   /*
    * TODO: Since we have really never used the patch factory and patch
    * level factory concepts beyond their defaults, should we remove them?
    */
   /*!
    * @brief Set the factory used to create patch objects.
    *
    * If a factory is not specified, then the default factory will create
    * patch objects of type Patch.
    *
    * @param[in]  factory
    */
   void
   setPatchFactory(
      const boost::shared_ptr<PatchFactory>& factory)
   {
      d_patch_factory = factory;
   }

   /*!
    * @brief Set the factory used to create patch level objects.
    *
    * If a factory is not specified, then the default factory will
    * create patch level objects of type PatchLevel.
    *
    * @param[in]  factory
    */
   void
   setPatchLevelFactory(
      const boost::shared_ptr<PatchLevelFactory>& factory)
   {
      d_patch_level_factory = factory;
   }

   /*!
    * @brief Get the grid geometry.
    *
    * @return a pointer to the grid geometry object.
    */
   boost::shared_ptr<BaseGridGeometry>
   getGridGeometry() const
   {
      return d_grid_geometry;
   }

   /*!
    * @brief Writes the state of the PatchHierarchy object and the PatchLevels
    * it contains to the database.
    *
    * @note
    * Only those patch data which have been registered for restart with
    * the VariableDatabase will be written to the database.
    * This method implements the pure virtual method in tbox::Serializable
    * class which is used by the tbox::RestartManager for writing the
    * PatchHierarchy to a restart file.
    * @par Assertions
    * When assertion checking is active, the database pointer must be non-null.
    *
    * @param[out]  database
    */
   void
   putToDatabase(
      const boost::shared_ptr<tbox::Database>& database) const;

   /*!
    * @brief Read in the entire hierarchy from the restart file.
    *
    * The database from which the restart data is read is determined by the
    * object_name specified in the constructor.
    *
    * @note
    * <ul>
    * <li> This method handles the memory allocation for each PatchLevel
    *      it reads in.
    * <li> The number of levels read in is the minimum of the d_max_levels
    *      and the number of levels stored in the database.
    * </ul>
    *
    * @par Assertions
    * When assertion checking is active, @c d_max_levels must be
    * greater than zero.  An unrecoverable assertion will result if the
    * database cannot be found in the restart file or the data in the
    * restart file is erroneous.
    */
   void
   getFromRestart();

   /*!
    * @brief Read in the entire hierarchy from the specified database.
    *
    * The component_selector specifies which patch data components
    * to read in from the database.
    *
    * @note
    * <ul>
    * <li>   Warning messages will be printed to the log file if
    *        any patch data component specified in the
    *        component_selector cannot be found in the database.
    * <li>   This method handles the memory allocation for each PatchLevel
    *        it reads in.
    * </ul>
    *
    * @par Assertions
    *
    * <ul>
    * <li>   The database argument must not be null.
    * <li>   The number of levels (if given) must be greater than zero.
    * </ul>
    *
    * @param[in]  database
    * @param[in]  component_selector
    */
   void
   getFromDatabase(
      const boost::shared_ptr<tbox::Database>& database,
      const ComponentSelector& component_selector);

   /*!
    * @brief Get the dimension of this object.
    *
    * @return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const
   {
      return d_dim;
   }

   /*!
    * @brief Print a patch hierarchy to a specified degree of detail.
    *
    * The amount of detail is controlled by the value of depth.
    *
    * depth <= 0 means only the number of levels will be printed.
    * depth == 1 means information about each level of the hierarchy will be
    * printed, as well as the total number of patches and cells in the
    * hierarchy.  depth >= 2 will print everything that is printed for
    * depth == 1, plus details about each patch on each level.
    *
    * If depth >= 1, global communication steps will be invoked, which can
    * have negative effects on scaling performance.
    *
    * @return 0.  Always.
    *
    * @param[in,out]    os The output stream
    * @param[in]        border string around output text
    * @param[in]     depth
    */
   int
   recursivePrint(
      std::ostream& os,
      const std::string& border = std::string(),
      int depth = 0);

   /*!
    * @brief Get the name of this object.
    *
    * @return The name of this object.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

private:
   /*
    * Static integer constant describing class's version number.
    */
   static const int HIER_PATCH_HIERARCHY_VERSION;

   /*!
    * @brief Writes the state of the PatchHierarchy object and the PatchLevels
    * it contains to the database.
    *
    * Only those patch data indicated in the ComponentSelector are written to
    * the specified database.
    *
    * @par Assertions
    * When assertion checking is active, the database pointer must be non-null.
    *
    * @param[out]  database
    * @param[in]  patchdata_write_table
    */
   void
   putToDatabase(
      const boost::shared_ptr<tbox::Database>& database,
      const ComponentSelector& patchdata_write_table) const;

   /*!
    * @brief Read input data from specified database and initialize
    * class members.
    *
    * When assertion checking is active, the database pointer must be
    * non-null.
    *
    * @param[in]  database   Input database
    */
   void
   getFromInput(
      const boost::shared_ptr<tbox::Database>& database);

   /*!
    * @brief Set up things for the entire class.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback()
   {
      /*
       * No-op.  This class doesn't
       */
   }

   /*!
    * @brief Free static timers.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   /*!
    * @brief String identifier of the object
    */
   std::string d_object_name;

   /*!
    * @brief Dimension of the object
    */
   const tbox::Dimension d_dim;

   /*!
    * @brief Boolean telling if object is registered for restart
    */
   bool d_registered_for_restart;

   /*!
    * @brief Number of levels currently in the hierarchy
    */
   int d_number_levels;

   int d_number_blocks;

   /*!
    * @brief Array of pointers to PatchLevels that make up the hierarchy
    */
   tbox::Array<boost::shared_ptr<PatchLevel> > d_patch_levels;

   /*!
    * @brief BaseGridGeometry that was used to construct the hierarchy
    */
   boost::shared_ptr<BaseGridGeometry> d_grid_geometry;

   /*!
    * @brief PatchDescriptor that is shared by every patch on the hierarchy
    */
   boost::shared_ptr<PatchDescriptor> d_patch_descriptor;

/*
 * TODO: Since we have really never used the patch factory and patch
 * level factory concepts beyond their defaults, should we remove them?
 */
/*!
 * @brief Factory used to create patches on the hierarchy
 */
   boost::shared_ptr<PatchFactory> d_patch_factory;

   /*!
    * @brief Factory used to create levels on the hierarchy
    */
   boost::shared_ptr<PatchLevelFactory> d_patch_level_factory;

   //@{
   //! @name Parameters for setting up the hierarchy.

   /*!
    * @brief Max number of levels that can be in the hierarchy.
    */
   int d_max_levels;

   /*!
    * @brief Ratio to coarser level.  d_ratio_to_coarser[0] is unused.
    *
    * The vector will be sized to d_max_levels.  d_ratio_to_coarser[n] is
    * the refinement ratio betwee level @c n and level @c n-1.
    */
   std::vector<IntVector> d_ratio_to_coarser;

   /*
    * @brief Proper nesting buffer for each level.
    *
    * The proper nesting buffer specifies the number of coarse cells
    * by which the next finer level is nested within the interior of
    * the domain of the next coarser level.  These buffer values
    * are used to compute the proper nesting boxes on each level.  They
    * can be used to guarantee that:
    *
    * (1) Adjacent cells on the composite grid differ by no more than
    *    one refinement level.
    *
    * (2) Data on the interior of fine patches can be filled by
    *    interpolation from the next coarser level if the buffer
    *    width is at least as large as the maximum stencil width
    *    of all interpolation operators.
    *
    * The proper nesting buffer size for each level may be specified in the
    * input file.  If not, a nesting buffer of 1 is used.  This value should
    * be suitable for most problems.  It is not recommended that the buffer
    * be set to 0 since this may cause features of the solution to move off
    * of refined mesh regions before subsequent regridding occurs.  If
    * additional buffering is required (e.g., if the stencil width for some
    * interpolation operator is greater than one), then it may be necessary
    * to increase this value.
    */
   std::vector<int> d_proper_nesting_buffer;

   /*!
    * @brief Smallest patch sizes for each level.
    */
   std::vector<IntVector> d_smallest_patch_size;

   /*!
    * @brief Largest patch sizes for each level.  A negative value means
    * unlimited.
    */
   std::vector<IntVector> d_largest_patch_size;

   /*!
    * @brief Whether to normally allow patches smaller than the max
    * ghost width.
    */
   bool d_allow_patches_smaller_than_ghostwidth;

   /*!
    * @brief Whether to allow patches smaller than d_smallest_patch_size
    * in order to prevent overlaping patches.
    */
   bool d_allow_patches_smaller_than_minimum_size_to_prevent_overlaps;

   /*!
    * @brief Required Connector width for self connectors.
    *
    * This is mutable because it may have to be updated by
    * getRequiredConnectorWidth(), which is a const method.
    *
    * See registerConnectorWidthRequestor() and getRequiredConnectorWidth().
    */
   mutable std::vector<IntVector> d_self_connector_widths;

   /*!
    * @brief Required Connector width for fine connectors.
    *
    * The width for coarse Connectors are the equivalent width for the
    * transpose Connector, times refinement ratio between the two
    * levels.
    *
    * This is mutable because it may have to be updated by
    * getRequiredConnectorWidth(), which is a const method.
    *
    * See registerConnectorWidthRequestor() and getRequiredConnectorWidth().
    */
   mutable std::vector<IntVector> d_fine_connector_widths;

   /*!
    * @brief Whether Connector widths have been computed.
    *
    * This is mutable because it is computed as required by the const
    * method getConnectorWidth().
    */
   mutable bool d_connector_widths_are_computed;

   /*!
    * @brief Vector of all ConnectorWidthRequestorStrategy objects registered
    * with this particular object.
    *
    * Note that this contains strategies registered with this object
    * using registerConnectorWidthRequestor().  It is
    * different from the strategies in s_class_cwrs.
    */
   std::vector<const ConnectorWidthRequestorStrategy *> d_individual_cwrs;

   /*!
    * @brief All ConnectorWidthRequestorStrategy registered with
    * PatchHierarchy's auto-registry.
    *
    * These are the implementations registered for the entire class,
    * using registerAutoConnectorWidthRequestorStrategy().
    */
   static std::vector<const ConnectorWidthRequestorStrategy *> s_class_cwrs;

   /*!
    * @brief Shutdown handler for clearing out static registry.
    */
   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;

   //@}

   //@{

   //! @name Domain-related objects.

/*
 * TODO: These things (if really needed) should be moved to the
 * BaseGridGeometry class.  However, the BaseGridGeometry object cannot own a
 * MappedBoxLevel because the BaseGridGeometry object is incapable of creating
 * a boost::shared_ptr to itself.  Might need to change BoxLevel to take a raw
 * pointer to BaseGridGeometry.
 */

   /*!
    * @brief Physical domain BoxLevel.
    *
    * All mapped_boxes in the domain
    * BoxLevel are owned by process 0.
    *
    * The physical domain BoxLevel is maintained in GLOBALIZED
    * mode with processor 0 owning all mapped boxes.
    */
   BoxLevel d_domain_mapped_box_level;

   //@}

};

}
}

#endif
