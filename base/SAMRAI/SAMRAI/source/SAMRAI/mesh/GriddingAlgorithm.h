/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   AMR hierarchy generation and regridding routines.
 *
 ************************************************************************/

#ifndef included_mesh_GriddingAlgorithm
#define included_mesh_GriddingAlgorithm

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/mesh/GriddingAlgorithmStrategy.h"
#include "SAMRAI/mesh/BoxGeneratorStrategy.h"
#include "SAMRAI/mesh/LoadBalanceStrategy.h"
#include "SAMRAI/mesh/GriddingAlgorithmConnectorWidthRequestor.h"
#include "SAMRAI/mesh/MultiblockGriddingTagger.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/hier/Connector.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/Utilities.h"

#include <boost/shared_ptr.hpp>
#include <iostream>
#include <string>
#include <vector>

#define GA_RECORD_STATS
// #undef GA_RECORD_STATS

#ifdef GA_RECORD_STATS
#include "SAMRAI/tbox/Statistic.h"
#include "SAMRAI/tbox/Statistician.h"
#endif

namespace SAMRAI {
namespace mesh {

/*!
 * @brief Class GriddingAlgorithm manages AMR patch hierarchy construction
 * operations in SAMRAI.  Specifically, it provides AMR patch hierarchy
 * generation and regridding routines that may be used with a variety
 * of AMR solution algorithms and application codes.
 *
 * The three main functions provided by this class are:
 *   - @b    makeCoarsestLevel()
 *      This routine constructs or repartitions
 *      the coarsest hierarchy level (level 0).
 *
 *   - @b    makeFinerLevel()
 *      This routine will attempt to add a new
 *      finest level to the hierarchy if the
 *      maximum number of levels allows it and
 *      cells on the current finest level are
 *      tagged for refinement.
 *
 *   - @b    regridAllFinerLevels()
 *      This routine will regrid all levels finer
 *      than some specified level based on cells
 *      that are tagged for refinement on each
 *      level finer than and including the given
 *      level.  This routine may add a new finest
 *      hierarchy level if the maximum number of
 *      levels allows it and cells on the current
 *      finest level are tagged for refinement.
 *      Levels may also be removed from the
 *      hierarchy if no cells are tagged.
 *
 *
 * These basic AMR operations are used to generate levels in
 * the AMR patch hierarchy at the beginning of a simulation, and regridding
 * collections of levels during an adaptive calculation.  More details are
 * found in the comments accompanying each member function below.
 *
 * Other objects passed to the
 * constructor provide the gridding algorithm with particular operations
 * needed during meshing operations.  Operations that tag cells for
 * refinement on a patch level and initialize data on new levels are
 * provided by the TagAndInitializeStrategy argument.  Operations that
 * cluster tagged cells into a boxes are provided by the
 * BoxGeneratorStrategy argument.  Routines that load balance patches
 * on a level are provided by the LoadBalanceStrategy constructor argument.
 *
 * Initialization of a GriddingAlgorithm object is performed via a
 * combination of default parameters and values read from an input
 * database.  Data read from input is summarized as follows:
 *
 * Optional input keys, data types, and defaults:
 *
 *   - \b    efficiency_tolerance
 *      An array of double values, each of which specifies the minimum
 *      fraction of tagged cells to total cells in boxes used to
 *      construct patches on a new level.  If the ratio is below the
 *      tolerance value, the box may be split into smaller boxes and
 *      pieces removed until the ratio becomes greater than or equal to
 *      the tolerance.  This tolerance helps users control the amount
 *      of extra refined cells created (beyond those tagged explicitly)
 *      that is typical in patch-based AMR computations.
 *      If no input values are given, a default of 0.8 is
 *      used.  See sample input below for input file format.
 *      The index of the value in the array corresponds to the number
 *      of the level to which the tolerance value applies.  If more
 *      values are given than the maximum number of levels allowed in
 *      the hierarchy, extra values will be ignored.  If fewer values
 *      are given, then the last value given will be used for each
 *      level without a specified input value.  For example, if only a
 *      single value is specified, then that value will be used on all
 *      levels.
 *
 *   - \b    combine_efficiency
 *      An array of double values, each of which serves as a threshold
 *      for the ratio of the total number of cells in two boxes into which
 *      a box may be split and the number of cells in the original box.
 *      If that ratio is greater than the combine efficiency, the box will not
 *      be split.  This tolerance helps users avoids splitting up portions
 *      of the domain into into very small patches which can increase
 *      the overhead of AMR operations.
 *      If no input values are given, a default of 0.8 is
 *      used.  See sample input below for input file format.
 *      Multiple values in the array are handled similar to
 *      efficiency_tolerance.
 *
 *   - @b check_nonnesting_user_boxes
 *      A flag to control how user-specified refinement boxes that violate
 *      proper nesting are handled.
 *      Set to one of these strings:
 *      @b "IGNORE" - nesting violations will be quietly disregarded.
 *      @b "WARN" - nesting violations will cause a warning but the
 *      code will continue anyway.
 *      @b "ERROR" - nesting violations will cause an unrecoverable
 *      assertion.
 *      The default is "ERROR".  We highly recommend making nesting violation
 *      an error.  The code may work anyway, but there are no guarantees.
 *
 *   - @b check_boundary_proximity_violation
 *      A flag to control how to resolve refinement boxes that violate
 *      boundary proximity (are less than the max ghost cell width of
 *      physical boundaries without touching the boundary).
 *      Set to one of these strings:
 *      @b "IGNORE" - violations will be quietly disregarded.
 *      @b "WARN" - violations will cause a warning but the
 *      code will continue anyway.
 *      @b "ERROR" - violations will cause an unrecoverable
 *      assertion.
 *      The default is "ERROR".  We highly recommend making boundary
 *      proximity violation an error.  The code may work anyway, but there
 *      are no guarantees.
 *
 *   - \b    check_nonrefined_tags
 *      A flag to control how to resolve user-specified tags that violate
 *      proper nesting.
 *
 *      If a tag violates the nesting requirements, its location in index space
 *      will not be refined when creating the finer level.  This flag allows the
 *      user to determine what to do when this occurs
 *
 *      Set to one of these strings:
 *      @b "IGNORE" - violating tags will be quietly disregarded.
 *      @b "WARN" - violating tags will cause a warning and be
 *      disregarded.
 *      @b "ERROR" - violating tags will cause an unrecoverable
 *      assertion.
 *      The default is "WARN".  It is fastest to ignore non-nesting tags
 *      because no checking has to be done.
 *
 *   - \b    check_overlapping_patches
 *      A flag to control checking for overlapping patches on a new level.
 *      Set to one of these strings:
 *      @b "IGNORE" - there is no check for overlapping patches,
 *      and they will be quietly disregarded.
 *      @b "WARN" - overlapping patches will cause a warning and be
 *      disregarded.
 *      @b "ERROR" - violating tags will cause an unrecoverable
 *      assertion.
 *      The default is "WARN".  The check for overlapping patches may be
 *      and should be bypassed by application that can tolerate overlaps.
 *      To prevent the creation of levels with overlapping patches, see
 *      the input flag
 *      "allow_patches_smaller_than_minimum_size_to_prevent_overlaps"
 *
 *   - \b   sequentialize_patch_indices
 *      A flag to specify whether patch indices will be globally sequentialized.
 *      This is not scalable, but is required for writing correct VisIt files.
 *      Due to the current VisIt requirement, this is currently true by default.
 *      It will evetually be set back to false after we remove the VisIt
 *      requirement.
 *
 *   - @b log_metadata_statistics = FALSE
 *      Whether to log metadata statistics after generating a new level.
 *      This flag writes out data that would be of interest to analyzing
 *      how metadata statistics affects performance.
 *
 *
 * Note that when continuing from restart, the input values in the
 * input file override all values read in from the restart database.
 *
 * The following represents sample input data for a three-level problem:
 *
 * \verbatim
 *
 *   // Optional input: different efficiency tolerance for each coarser level
 *   efficiency_tolerance = 0.80e0, 0.85e0, 0.90e0
 *
 *   // Optional input: combine efficiency is same for all levels.
 *   combine_efficiency = 0.95e0
 *
 * \endverbatim
 *
 * @see mesh::TagAndInitializeStrategy
 * @see mesh::LoadBalanceStrategy
 * @see mesh::BoxGeneratorStrategy
 */

class GriddingAlgorithm:
   public GriddingAlgorithmStrategy,
   public tbox::Serializable
{
public:
   /*!
    * @brief The constructor for GriddingAlgorithm configures the
    * gridding algorithm with the patch hierarchy and concrete algorithm
    * strategy objects in the argument list.
    *
    * Gridding parameters are initialized from values provided in the
    * specified input and in the restart database corresponding to the
    * specified object_name argument.  The constructor also registers
    * this object for restart using the specified object name when the
    * boolean argument is true (default).
    *
    * If assertion checking is turned on, an unrecoverable assertion
    * will result if any of the required pointer arguments is null.
    * Assertions may also be thrown if any checks for consistency
    * among input parameters fail.
    *
    * @param[in] hierarchy The hierarchy that this GriddingAlgorithm will
    * work on.  The pointer is cached.  All hierarchy operations will
    * be on this hierarchy.
    *
    * @param[in] object_name For registering the object in the restart
    * database.
    *
    * @param[in] input_db
    *
    * @param[in] level_strategy
    *
    * @param[in] generator
    *
    * @param[in] balancer Load balancer
    *
    * @param[in] balancer_zero Special load balancer to use for level
    * zero.  If omitted, will use @c balancer instead.
    *
    * @param[in] register_for_restart
    */
   GriddingAlgorithm(
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const std::string& object_name,
      const boost::shared_ptr<tbox::Database>& input_db,
      const boost::shared_ptr<TagAndInitializeStrategy>& level_strategy,
      const boost::shared_ptr<BoxGeneratorStrategy>& generator,
      const boost::shared_ptr<LoadBalanceStrategy>& balancer,
      const boost::shared_ptr<LoadBalanceStrategy>& balancer_zero =
         boost::shared_ptr<LoadBalanceStrategy>(),
      bool register_for_restart = true);

   /*!
    * @brief Destructor
    */
   virtual ~GriddingAlgorithm();

   /*!
    * @brief Create or rebalance the coarsest level.
    *
    * This is an implementation of interface defined in GriddingAlgorithmStrategy.
    *
    * This routine will attempt to construct the coarsest level in the AMR
    * patch hierarchy (i.e., level 0).  If level 0 does not already exist,
    * then the domain specification is checked against the constraints of
    * the grid generation procedures.  Recall that the domain specification
    * is maintained by the grid geometry object associated with the hierarchy.
    * Generally, an unrecoverable assertion will result if the constraints
    * are not satisfied.
    *
    * If level 0 already exists in the hierarchy, then the routine
    * will generate a new level zero by re-applying the load balancing
    * procedure to the existing level.  Data will be moved from the
    * old level to the new level and the pre-existing level 0 will be
    * discarded.  Note that this routine is different than the routine
    * makeFinerLevel() below, which is used to construct levels finer
    * than level zero.  In particular, this routine does not select
    * cells for refinement, whereas the other routine does.
    *
    * @param[in] level_time Simulation time.
    */
   void
   makeCoarsestLevel(
      const double level_time);

   /*!
    * @brief Attempt to create a new level in the hierarchy finer than
    * the finest level currently residing in the hierarchy.
    *
    * This is an implementation of interface method
    * GriddingAlgorithmStrategy::makeFinerLevel().
    *
    * The tag buffer indicates the number of cells by which cells
    * selected for refinement should be buffered before new finer
    * level boxes are constructed.  All tagged cells should be refined
    * except where refinement would violate proper nesting.  The
    * buffer is meant to keep phenomena of interest on refined regions
    * of the mesh until adaptive regridding occurs next.  Callers of
    * this method should take into account how the simulation may
    * evolve before regridding occurs (e.g., number of timesteps
    * taken) when calculating the tag_buffer.
    *
    * @param[in] level_time See above text.
    *
    * @param[in] initial_time See above text.
    *
    * @param[in] tag_buffer See above text.
    *
    * @param[in] regrid_start_time The simulation time when the
    * regridding operation began (this parameter is ignored except
    * when using Richardson extrapolation)
    */
   void
   makeFinerLevel(
      const double level_time,
      const bool initial_time,
      const int tag_buffer,
      const double regrid_start_time = 0.0);

   /*!
    * @brief Attempt to regrid each level in the PatchHierarchy
    * that is finer than the specified level.
    *
    * This method implements the virtual interface
    * GriddingAlgorithmStrategy::regridAllFinerLevels().
    *
    * Note that the current algorithm permits at most one new finest level
    * to be added to the hierarchy with each invocation of the regridding
    * process.  This constraint, though seemingly restrictive makes the
    * process of maintaining properly nested levels much easier.
    *
    * Note that the current algorithm permits at most one new finest level
    * to be added to the hierarchy with each call to this method.
    * This constraint, though seemingly restrictive makes the
    * process of maintaining properly nested levels much easier.
    *
    * Important note: If assertion checking is activated, several
    * checks are applied to the functions arguments.  If any check is
    * violated, an unrecoverable assertion will result.  In
    * particular, the given level number must match that of of some
    * level in the hierarchy.  Also, the tag buffer array must contain
    * a positive value for each level in the hierarchy.
    *
    * @param[in] level_number Coarsest level on which cells will be
    * tagged for refinement
    *
    * @param[in] regrid_time Simulaition time when regridding occurs
    *
    * @param[in] tag_buffer Size of buffer on each level around tagged
    * cells that will be covered by the next finer level
    *
    * @param[in] regrid_start_time The simulation time when the
    * regridding operation began on each level (this parameter is
    * ignored except when using Richardson extrapolation)
    *
    * @param[in] level_is_coarsest_to_sync Level is the coarsest to sync
    */
   void
   regridAllFinerLevels(
      const int level_number,
      const double regrid_time,
      const tbox::Array<int>& tag_buffer,
      tbox::Array<double> regrid_start_time = tbox::Array<double>(),
      const bool level_is_coarsest_to_sync = true);

   /*!
    * @brief Return pointer to level gridding strategy data member.
    *
    * Access to this member is useful when an integrator implementation
    * needs to know if the error estimator uses time integration.
    *
    * @return pointer to TagAndInitializeStrategy data member.
    */
   virtual
   boost::shared_ptr<TagAndInitializeStrategy>
   getTagAndInitializeStrategy() const;

   /*!
    * @brief Return pointer to load balance strategy data member.
    */
   virtual
   boost::shared_ptr<LoadBalanceStrategy>
   getLoadBalanceStrategy() const;

   /*!
    * @brief Return pointer to load balance strategy specialized for
    * balancing level zero.
    *
    * @return pointer to load balance strategy specialized for the case
    * where one processor owns all the initial loads.
    */
   virtual
   boost::shared_ptr<LoadBalanceStrategy>
   getLoadBalanceStrategyZero() const;

   /*!
    * @brief Set efficiency tolerance for clustering tags on level.
    *
    * @param[in] efficiency_tolerance
    * @param[in] level_number
    */
   void
   setEfficiencyTolerance(
      const double efficiency_tolerance,
      const int level_number)
   {
      TBOX_ASSERT((level_number >= 0) &&
                  (level_number < d_hierarchy->getMaxNumberOfLevels()));
      TBOX_ASSERT((efficiency_tolerance >= 0) &&
                  (efficiency_tolerance <= 1.0));
      d_efficiency_tolerance[level_number] = efficiency_tolerance;
   }

   /*!
    * @brief Return efficiency tolerance for clustering tags on level.
    *
    * @return efficiency tolerance for clustering tags on level.
    */
   double
   getEfficiencyTolerance(
      const int level_number) const
   {
      TBOX_ASSERT((level_number >= 0) &&
                  (level_number < d_hierarchy->getMaxNumberOfLevels()));
      int size = d_efficiency_tolerance.getSize();
      return (level_number < size) ?
         d_efficiency_tolerance[level_number] :
         d_efficiency_tolerance[size - 1];
   }

   /*!
    * @brief Set combine efficiency for clustering tags on level.
    *
    * @param[in] combine_efficiency
    * @param[in] level_number
    */
   void
   setCombineEfficiency(
      const double combine_efficiency,
      const int level_number)
   {
      TBOX_ASSERT((level_number >= 0) &&
                  (level_number < d_hierarchy->getMaxNumberOfLevels()));
      TBOX_ASSERT((combine_efficiency >= 0) && (combine_efficiency <= 1.0));
      d_combine_efficiency[level_number] = combine_efficiency;
   }

   /*!
    * @brief Return combine efficiency for clustering tags on level.
    *
    * @return combine efficiency for clustering tags on level.
    */
   double
   getCombineEfficiency(
      const int level_number) const
   {
      TBOX_ASSERT((level_number >= 0) &&
                  (level_number < d_hierarchy->getMaxNumberOfLevels()));
      int size = d_combine_efficiency.getSize();
      return (level_number < size) ?
         d_combine_efficiency[level_number] :
         d_combine_efficiency[size - 1];
   }

   /*!
    * @brief Print all data members of the class instance to given output stream.
    */
   void
   printClassData(
      std::ostream& os) const;

   /*!
    * @brief Write object state out to the gien database.
    *
    * When assertion checking is active, the database pointer must be non-null.
    */
   void
   putToDatabase(
      const boost::shared_ptr<tbox::Database>& db) const;

   /*
    * @brief Write out statistics recorded on numbers of cells and patches generated.
    */
   void
   printStatistics(
      std::ostream& s = tbox::plog) const;

private:
   /*
    * Static integer constant describing this class's version number.
    */
   static const int ALGS_GRIDDING_ALGORITHM_VERSION;

   //! @brief Shorthand typedef.
   typedef hier::Connector::NeighborSet NeighborSet;

   /*
    * @brief Read input data from specified database and initialize class members.
    *
    * When assertion checking is active, the database pointer must be non-null.
    *
    * @param[in] db
    *
    * @param[in] is_from_restart Should be set to true if the
    * simulation is from restart.  Otherwise, it should be set to
    * false.
    */
   void
   getFromInput(
      const boost::shared_ptr<tbox::Database>& db,
      bool is_from_restart);

   /*
    * @brief Read object state from the restart file and initialize
    * class data members.
    *
    * The database from which the restart data is read is determined
    * by the object_name specified in the constructor.
    *
    * Unrecoverable Errors:
    *
    *   -The database corresponding to object_name is not found
    *    in the restart file.
    *
    *   -The class version number and restart version number do not
    *    match.
    */
   void
   getFromRestart();

   /*
    * @brief Recursively regrid the hierarchy level and all finer
    * levels in the hierarchy.
    *
    * This private member function is invoked by the
    * regridAllFinerLevels() routine.  It may invoke recursively
    * invoke itself.
    */
   void
   regridFinerLevel(
      const int level_number,
      const double regrid_time,
      const int finest_level_not_regridded,
      const bool level_is_coarsest_to_sync,
      const tbox::Array<int>& tag_buffer,
      const tbox::Array<double>& regrid_start_time = tbox::Array<double>(0));

   /*
    * @brief Tagging stuff before recursive regrid, called from regridFinerLevel.
    */
   void
   regridFinerLevel_doTaggingBeforeRecursiveRegrid(
      const int tag_ln,
      const bool level_is_coarsest_sync_level,
      const tbox::Array<double>& regrid_start_time,
      const double regrid_time);

   /*
    * @brief Tagging stuff after recursive regrid, called from regridFinerLevel.
    *
    * A side-effect of this method is setting the overlap Connectors
    * between the tag level (number tag_ln) and the finer level
    * (number tag_ln+2) if the finer level exists in the hierarchy.
    */
   void
   regridFinerLevel_doTaggingAfterRecursiveRegrid(
      hier::Connector& tag_to_finer,
      hier::Connector& finer_to_tag,
      const int tag_ln,
      const tbox::Array<int>& tag_buffer);

   /*
    * @brief Given the metadata describing the new level, this method
    * creates and installs new PatchLevel in the hierarchy.
    */
   void
   regridFinerLevel_createAndInstallNewLevel(
      const int tag_ln,
      const double regrid_time,
      hier::Connector* tag_to_new,
      hier::Connector* new_to_tag,
      const hier::Connector& tag_to_finer,
      const hier::Connector& finer_to_tag);

   /*
    * @brief Set all tags on a level to a given value.
    *
    * @param[in] tag_value
    *
    * @param[in] level
    *
    * @param[in] tag_index
    */
   void
   fillTags(
      const int tag_value,
      const boost::shared_ptr<hier::PatchLevel>& tag_level,
      const int tag_index) const;

   /*
    * @brief Set integer tags to specified value on intersection
    * between patch level and a BoxLevel
    *
    * The index value corresponds to the patch descriptor entry of the
    * cell-centered integer tag array.  The boolean flag indicates
    * whether tags are to be set on the regions corresponding to the
    * interior of the level only, if the tag arrays contain ghosts.
    *
    * @param[in] tag_value
    *
    * @param[in] level
    *
    * @param[in] index
    *
    * @param[in] level_to_fill_mapped_box_level Connector from the
    * level with the tags to the BoxLevel describing where to
    * fill.
    *
    * @param[in] interior_only
    *
    * @param fill_box_growth How much to grow fill boxes before using
    * them to tag.  Must be in index space of level holding fill boxes.
    */
   void
   fillTagsFromBoxLevel(
      const int tag_value,
      const boost::shared_ptr<hier::PatchLevel>& tag_level,
      const int index,
      const hier::Connector& tag_level_to_fill_mapped_box_level,
      const bool interior_only,
      const hier::IntVector& fill_box_growth) const;

   /*
    * @brief Make a map that, when applied to an improperly nested
    * BoxLevel, removes the nonnesting parts.
    *
    * @param[out] nested_mapped_box_level  The nesting parts of the
    * unnested BoxLevel.
    *
    * @param[out] unnested_to_nested  The mapping from the unnested
    * BoxLevel to the nested BoxLevel.
    *
    * @param[in] unnested_mapped_box_level
    *
    * @param[in] tag_to_unnested  Overlap Connector from the tag level
    * to unnested_mapped_box_level.
    *
    * @param[in] unnested_to_tag  Transpose of tag_to_unnested.
    *
    * @param[in] unnested_ln Level number of PatchLevel being
    * generated (one more than the tag level number).
    */
   void
   makeProperNestingMap(
      hier::BoxLevel& nested_mapped_box_level,
      hier::Connector& unnested_to_nested,
      const hier::BoxLevel& unnested_mapped_box_level,
      const hier::Connector& tag_to_unnested,
      const hier::Connector& unnested_to_tag,
      const int unnested_ln) const;

   /*
    * @brief Make a map that, when applied to a BoxLevel that
    * extends outside of a reference BoxLevel, removes those
    * outside parts.
    *
    * @param[out] nested_mapped_box_level  The nesting parts of the
    * unnested BoxLevel.
    *
    * @param[out] unnested_to_nested  The mapping from the unnested
    * BoxLevel to the nested BoxLevel.
    *
    * @param[in] unnested_mapped_box_level
    *
    * @param[in] unnested_to_reference
    */
   void
   makeOverflowNestingMap(
      hier::BoxLevel& nested_mapped_box_level,
      hier::Connector& unnested_to_nested,
      const hier::BoxLevel& unnested_mapped_box_level,
      const hier::Connector& unnested_to_reference) const;

   /*!
    * @brief Make a map from a BoxLevel to parts of that BoxLevel
    * that violate proper nesting.
    *
    * @param[in] candidate BoxLevel being examined for nesting violation.
    *
    * @param[out] violator BoxLevel containing violating parts of candidate.
    *
    * @param[in] tag_ln Level number of the level that candidate should nest in.
    *
    * @param[in] candidate_to_hierarchy Connector to mapped_box_level number
    *       tag_ln in the hierarchy.
    *
    * @param[in] hierarchy_to_candidate Connector from mapped_box_level number
    *       tag_ln in the hierarchy.
    */
   void
   computeNestingViolator(
      hier::BoxLevel& violator,
      hier::Connector& candidate_to_violator,
      const hier::BoxLevel& candidate,
      const hier::Connector& candidate_to_hierarchy,
      const hier::Connector& hierarchy_to_candidate,
      const int tag_ln) const;

   /*
    * @brief Extend Boxes to domain boundary if they they are too close.
    *
    * See hier::BoxUtilities::extendBoxToDomainBoundary().
    *
    * @param[in,out] new_mapped_box_level
    *
    * @param[in,out] tag_to_new
    *
    * @param[in,out] new_to_tag
    *
    * @param[in] physical_domain_array  Array holding domain for each block
    *
    * @param[in] extend_ghosts Extend the boxes to the boundary if
    * they less than extend_ghosts cells from the boundary.
    */
   void
   extendBoxesToDomainBoundary(
      hier::BoxLevel& new_mapped_box_level,
      hier::Connector& tag_to_new,
      hier::Connector& new_to_tag,
      const tbox::Array<hier::BoxContainer>& physical_domain_array,
      const hier::IntVector& extend_ghosts) const;

   /*!
    * @brief Precompute data used to define proper nesting.
    *
    * Data is associated with level number ln, to be used for
    * constructing level number ln+1.
    *
    * @param[in] ln
    */
   void
   computeProperNestingData(
      const int ln);

   /*!
    * @brief Attempt to fix boxes that are too small by growing them
    * within the nesting-restricted domain.
    *
    * @param[in,out] new_mapped_box_level BoxLevel being formed
    * into the new new level.
    *
    * @param[in,out] tag_to_new  Connector to be updated.
    *
    * @param[in,out] new_to_tag  Connector to be updated.
    *
    * @param[in] min_size
    *
    * @param[in] tag_ln Level number of the tag level.
    */
   void
   growBoxesWithinNestingDomain(
      hier::BoxLevel& new_mapped_box_level,
      hier::Connector& tag_to_new,
      hier::Connector& new_to_tag,
      const hier::IntVector& min_size,
      const int tag_ln) const;

   /*
    * @brief Refine a BoxLevel from the resolution of the tag
    * level to the resolution of the level being created.
    *
    * @param[in,out] new_mapped_box_level BoxLevel being formed
    * into the new new level.
    *
    * @param[in,out] tag_to_new  Connector to be updated.
    *
    * @param[in,out] new_to_tag  Connector to be updated.
    *
    * @param[in] ratio
    */
   void
   refineNewBoxLevel(
      hier::BoxLevel& new_mapped_box_level,
      hier::Connector& tag_to_new,
      hier::Connector& new_to_tag,
      const hier::IntVector& ratio) const;

   /*
    * @brief Renumber Boxes in a BoxLevel.
    *
    * This method renumbers Boxes to give them globally
    * sequential LocalIds and/or sorts the Boxes by their
    * coordinates.  Sequentializing the global indices numbers them
    * sequentially, like Patch numbering in the traditional SAMR
    * approach.  Sorting by box coordinates removes the randomness in
    * the ordering of boxes.
    *
    * @param[in,out] new_mapped_box_level
    *
    * @param[in,out] tag_to_new
    *
    * @param[in,out] new_to_tag
    *
    * @param[in] sort_by_corners
    *
    * @param[in] sequentialize_global_indices
    */
   void
   renumberBoxes(
      hier::BoxLevel& new_mapped_box_level,
      hier::Connector& tag_to_new,
      hier::Connector& new_to_tag,
      bool sort_by_corners,
      bool sequentialize_global_indices) const;

   /*
    * @brief Buffer each integer tag on patch level matching given tag
    * value with a border of matching tags.
    */
   void
   bufferTagsOnLevel(
      const int tag_value,
      const boost::shared_ptr<hier::PatchLevel>& level,
      const int buffer_size) const;

   /*
    * @brief Set the new level boxes using information stored in a file.
    *
    * If cell tagging is not performed, the new level boxes may
    * be specified either from user input (by specifying "REFINE_BOXES"
    * input) or from output from a previous run stored in an
    * HDF file.  This method sets the "new_level_boxes" based on
    * the information in the file.  It also sets the boolean
    * parameter "remove_old_fine_level" which indicates whether
    * the level box configuration has changed and, consequently,
    * whether we need to remove the old level.
    */
   void
   readLevelBoxes(
      hier::BoxLevel& new_mapped_box_level,
      hier::Connector& coarser_to_new,
      hier::Connector& new_to_coarser,
      const int level_number,
      const double regrid_time,
      bool& remove_old_fine_level);

   /*
    * @brief Given the number for the level where cells are tagged for
    * refinement, compute a BoxLevel from which a refinement of
    * the level may be constructed.
    *
    * It is assumed that the integer tags that identify cells for
    * refinement have been set on the level before this routine is
    * called.  At the end of this function, new_mapped_box_level will
    * represent the box extents of a new fine level on the given
    * block.
    */
   void
   findRefinementBoxes(
      hier::BoxLevel& new_mapped_box_level,
      hier::Connector& tag_to_new,
      hier::Connector& new_to_tag,
      const int tag_ln) const;

   /*
    * @brief Set patch size and ghost cell information needed to
    * create new refinement levels.
    *
    * This method applies to levels that are being used to build new
    * finer levels (i.e. level_number is a coarser level in the
    * hierarchy) and to levels which are simply being reconstructed
    * (i.e. the same level in the hierarchy).  The boolean value
    * "for_building_finer" controls the logic for the two cases - in
    * the former case, it is true while in the latter case it is
    * false.
    *
    * When a finer level is being constructed, the maximum number of ghost
    * cells needed in any variable is used to compute the smallest patch
    * size allowed and the extent to which patches may be extended to touch
    * the physical boundary.  This avoids problems in setting ghost cell
    * data that may occur when ghost cell regions intersect the physical
    * boundary in strange ways.
    *
    * This routine sets the smallest and largest patch sizes for the specified
    * level, the smallest box to refine on the next coarser level, and the
    * number of cells that a patch may be extended to the boundary if it
    * sufficiently close to the boundary (extend_ghosts).
    */
   void
   getGriddingParameters(
      hier::IntVector& smallest_patch,
      hier::IntVector& smallest_box_to_refine,
      hier::IntVector& largest_patch,
      hier::IntVector& extend_ghosts,
      const int level_number,
      const bool for_building_finer) const;

   /*!
    * @brief Check domain boxes for violations of certain constraints.
    */
   void
   checkDomainBoxes(
      const hier::BoxContainer& domain_boxes) const;

   /*!
    * @brief Check for non-nesting user-specified boxes.
    */
   void
   checkNonnestingUserBoxes(
      const hier::Connector& new_to_tag,
      const hier::IntVector& nesting_buffer) const;

   /*!
    * @brief Check for boxes that are too close to the physical
    * boundary without touching it.
    *
    * Boxes close to the physical boundaries causes ghost boxes to
    * intersect the boundary in weird ways, so we disallow it.  This
    * method writes a warning describing each violation found.
    *
    * @param[in] mapped_box_level
    *
    * @param[in] extend_ghosts
    *
    * Return the number of violations found.
    */
   size_t
   checkBoundaryProximityViolation(
      const hier::BoxLevel& mapped_box_level,
      const hier::IntVector& extend_ghosts) const;

   /*!
    * @brief Check for boxes that are too close to the physical
    * boundary without touching it.
    *
    * Compute the allowable distance from boxes in new_mapped_box_level
    * to domain boundary and delegate the checking.
    */
   void
   checkBoundaryProximityViolation(
      const int tag_ln,
      const hier::BoxLevel& new_mapped_box_level) const;

   /*!
    * @brief Check for user tags that violate proper nesting.
    *
    * @param tag_level  Tag level
    *
    * @param tag_ln  Tag level number
    */
   void
   checkNonrefinedTags(
      const hier::PatchLevel& tag_level,
      int tag_ln) const;

   /*
    * @brief Reset data that handles tag buffering.
    *
    * Resets the tag buffering data so that it will be able to handle a
    * a tag buffer of the given size.
    *
    * @param tag_buffer  The size of the buffer
    */
   void
   resetTagBufferingData(const int tag_buffer);

   /*!
    * @brief Check for overlapping patches within a level.
    *
    * @param[in] mapped_box_level_to_self
    */
   void
   checkOverlappingPatches(
      const hier::Connector& mapped_box_level_to_self) const;

   /*!
    * @brief Warn if the domain is too small any periodic direction.
    */
   void
   warnIfDomainTooSmallInPeriodicDir() const;

   /*!
    * @brief Allocate timers.
    */
   void
   allocateTimers();

   /*!
    * @brief Log metadata statistics after generating a new level.
    */
   void logMetadataStatistics(
      const char *caller_name,
      int ln,
      bool log_fine_connector,
      bool log_coarse_connector) const;

   /*!
    * @brief Initialize static objects and register shutdown routine.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   startupCallback()
   {
      s_tag_indx = new tbox::Array<int>(
         tbox::Dimension::MAXIMUM_DIMENSION_VALUE,
         -1);
      s_buf_tag_indx = new tbox::Array<int>(
         tbox::Dimension::MAXIMUM_DIMENSION_VALUE,
         -1);
   }

   /*!
    * @brief Method registered with ShutdownRegister to cleanup statics.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   shutdownCallback()
   {
      delete s_tag_indx;
      delete s_buf_tag_indx;
   }

   /*!
    * @brief Initialize static objects and register shutdown routine.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback()
   {
   }

   /*!
    * @brief Method registered with ShutdownRegister to cleanup statics.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback()
   {
   }

   /*
    * @brief Record statistics on how many patches and cells were generated.
    */
   void
   recordStatistics(
      double current_time);

   /*!
    * @brief The hierarchy that this GriddingAlgorithm works on.
    */
   boost::shared_ptr<hier::PatchHierarchy> d_hierarchy;

   /*!
    * @brief Implementation registered with the hierarchy, telling the
    * hierarchy what width the GriddingAlgorithm will be requesting.
    */
   GriddingAlgorithmConnectorWidthRequestor d_connector_width_requestor;

   /*
    * Static members for managing shared tag data among multiple
    * GriddingAlgorithm objects.
    */
   static tbox::Array<int> * s_tag_indx;
   static tbox::Array<int> * s_buf_tag_indx;

   const tbox::Dimension d_dim;

   hier::IntVector d_buf_tag_ghosts;

   /*
    * The object name is used for error reporting and accessing
    * restart file information.  Whether the object writes restart
    * file data depends on the value of this boolean which is
    * set in the constructor.
    */
   std::string d_object_name;
   bool d_registered_for_restart;

   /*
    * Data members that manage application-specific level initialization
    * and cell tagging, clustering of tagged cells into boxes, and load
    * balancing of patches to processors, respectively.
    */
   boost::shared_ptr<TagAndInitializeStrategy> d_tag_init_strategy;
   boost::shared_ptr<BoxGeneratorStrategy> d_box_generator;
   boost::shared_ptr<LoadBalanceStrategy> d_load_balancer;
   boost::shared_ptr<LoadBalanceStrategy> d_load_balancer0;

   /*
    * MultiblockGriddingTagger is the RefinePatchStrategy
    * implementation provided to the RefineSchedule d_bdry_sched_tags
    * inside GriddingAlgorithm.
    *
    * @see setMultiblockGriddingTagger().
    */
   MultiblockGriddingTagger * d_mb_tagger_strategy;

   /*
    * Cell-centered integer variables use to tag cells for refinement.
    * The descriptor index d_tag_indx is used to obtain tag information
    * from user-defined routines on patch interiors.  The descriptor index
    * d_buf_tag_indx is used to buffer tags on patches that may be
    * distributed across processors.  The refine algorithm and schedule are
    * used for interprocessor communication.
    */
   boost::shared_ptr<pdat::CellVariable<int> > d_tag;
   boost::shared_ptr<pdat::CellVariable<int> > d_buf_tag;
   int d_tag_indx;
   int d_buf_tag_indx;

   boost::shared_ptr<xfer::RefineAlgorithm> d_bdry_fill_tags;
   tbox::Array<boost::shared_ptr<xfer::RefineSchedule> > d_bdry_sched_tags;

   /*
    * True and false integer tag values set in constructor and used to
    * set and compare integer tag values on the patch hierarchy.  Generally,
    * these variables are easier to read in the code than integer constants,
    * such as 0 and 1.
    */
   const int d_true_tag;
   const int d_false_tag;

   /*!
    * @brief Finest level not changed during regrid.
    *
    * This member is temporary, used to coordinate with private methods
    * during mesh generation and regridding.  It has a value of -1 when
    * not being used.
    */
   int d_base_ln;

   /*
    * @brief Efficiency tolerance during clustering.
    *
    * See input parameter efficiency_tolerance.
    */
   tbox::Array<double> d_efficiency_tolerance;

   /*
    * @brief Combine efficiency during clustering.
    *
    * See input parameter combine_efficiency.
    */
   tbox::Array<double> d_combine_efficiency;

   /*
    * @brief When regridding level ln+1, the new level ln must not flow into
    * d_proper_nesting_complement[ln].
    *
    * Has length d_hierarchy->getMaxNumberOfLevels().  The objects are
    * initialized only during gridding/regridding.
    */
   std::vector<hier::BoxLevel> d_proper_nesting_complement;

   /*
    * @brief Connectors from the hierarchy to d_proper_nesting_complement.
    *
    * d_to_nesting_complement[ln] goes from level ln to
    * d_proper_nesting_complelemt[ln].
    */
   std::vector<hier::Connector> d_to_nesting_complement;

   /*
    * @brief Connectors from d_proper_nesting_complement to the hierarchy.
    *
    * d_from_nesting_complement[ln] goes from
    * d_proper_nesting_complement[ln] to level ln.
    */
   std::vector<hier::Connector> d_from_nesting_complement;

   /*!
    * @brief How to resolve user tags that violate nesting requirements.
    *
    * See input parameter "check_nonrefined_tags".
    */
   char d_check_nonrefined_tags;

   /*!
    * @brief Whether to check for overlapping patches on a level.
    *
    * See input parameter "check_overlapping_patches".
    */
   char d_check_overlapping_patches;

   /*!
    * @brief Whether to check for non-nesting user-specified refine
    * boxes.
    *
    * See input parameter "check_nonnesting_user_boxes".
    */
   char d_check_nonnesting_user_boxes;

   /*!
    * @brief Whether to check for user-specified refine boxes that
    * violate boundary proximity.
    *
    * See input parameter "check_boundary_proximity_violation".
    */
   char d_check_boundary_proximity_violation;

   /*!
    * @brief Whether to globally sequentialize patch indices on every level.
    */
   bool d_sequentialize_patch_indices;

   /*!
    * @brief Whether to log metadata statistics after generating a new
    * level.
    *
    * See input parameter description.
    */
   bool d_log_metadata_statistics;

   /*
    * Switches for massaging boxes after clustering.  Should be on for
    * most AMR applications.  Turning off is mainly for debugging
    * purposes.
    */
   bool d_enforce_proper_nesting;
   bool d_extend_to_domain_boundary;
   bool d_load_balance;

   //@{
   //! @name Used for evaluating peformance.
   bool d_barrier_and_time;
   //@}

   /*
    * Timers interspersed throughout the class.
    */
   boost::shared_ptr<tbox::Timer> t_find_domain_complement;
   boost::shared_ptr<tbox::Timer> t_load_balance;
   boost::shared_ptr<tbox::Timer> t_load_balance0;
   boost::shared_ptr<tbox::Timer> t_load_balance_setup;
   boost::shared_ptr<tbox::Timer> t_bdry_fill_tags_create;
   boost::shared_ptr<tbox::Timer> t_make_coarsest;
   boost::shared_ptr<tbox::Timer> t_make_finer;
   boost::shared_ptr<tbox::Timer> t_make_finer_setup;
   boost::shared_ptr<tbox::Timer> t_make_finer_tagging;
   boost::shared_ptr<tbox::Timer> t_make_finer_create;
   boost::shared_ptr<tbox::Timer> t_regrid_all_finer;
   boost::shared_ptr<tbox::Timer> t_regrid_finer_create;
   boost::shared_ptr<tbox::Timer> t_bridge_links;
   boost::shared_ptr<tbox::Timer> t_fill_tags;
   boost::shared_ptr<tbox::Timer> t_tag_cells_for_refinement;
   boost::shared_ptr<tbox::Timer> t_buffer_tags;
   boost::shared_ptr<tbox::Timer> t_bdry_fill_tags_comm;
   boost::shared_ptr<tbox::Timer> t_second_finer_tagging;
   boost::shared_ptr<tbox::Timer> t_find_refinement;
   boost::shared_ptr<tbox::Timer> t_bridge_new_to_new;
   boost::shared_ptr<tbox::Timer> t_find_new_to_new;
   boost::shared_ptr<tbox::Timer> t_bridge_new_to_coarser;
   boost::shared_ptr<tbox::Timer> t_bridge_new_to_finer;
   boost::shared_ptr<tbox::Timer> t_bridge_new_to_old;
   boost::shared_ptr<tbox::Timer> t_find_boxes_containing_tags;
   boost::shared_ptr<tbox::Timer> t_enforce_nesting;
   boost::shared_ptr<tbox::Timer> t_make_nesting_map;
   boost::shared_ptr<tbox::Timer> t_make_nesting_map_compute;
   boost::shared_ptr<tbox::Timer> t_make_nesting_map_convert;
   boost::shared_ptr<tbox::Timer> t_use_nesting_map;
   boost::shared_ptr<tbox::Timer> t_make_overflow_map;
   boost::shared_ptr<tbox::Timer> t_make_overflow_map_compute;
   boost::shared_ptr<tbox::Timer> t_make_overflow_map_convert;
   boost::shared_ptr<tbox::Timer> t_use_overflow_map;
   boost::shared_ptr<tbox::Timer> t_compute_external_parts;
   boost::shared_ptr<tbox::Timer> t_compute_nesting_violator;
   boost::shared_ptr<tbox::Timer> t_extend_to_domain_boundary;
   boost::shared_ptr<tbox::Timer> t_extend_within_domain;
   boost::shared_ptr<tbox::Timer> t_grow_boxes_within_domain;
   boost::shared_ptr<tbox::Timer> t_sort_nodes;
   boost::shared_ptr<tbox::Timer> t_modify_connector;
   boost::shared_ptr<tbox::Timer> t_misc1;
   boost::shared_ptr<tbox::Timer> t_misc2;
   boost::shared_ptr<tbox::Timer> t_misc3;
   boost::shared_ptr<tbox::Timer> t_misc4;
   boost::shared_ptr<tbox::Timer> t_misc5;
   boost::shared_ptr<tbox::Timer> t_make_domain;
   boost::shared_ptr<tbox::Timer> t_get_balance;
   boost::shared_ptr<tbox::Timer> t_use_balance;
   boost::shared_ptr<tbox::Timer> t_make_new;
   boost::shared_ptr<tbox::Timer> t_process_error;
   boost::shared_ptr<tbox::Timer> t_limit_overflow;
   boost::shared_ptr<tbox::Timer> t_reset_hier;
   boost::shared_ptr<tbox::Timer> t_box_massage;

#ifdef GA_RECORD_STATS
   /*
    * Statistics on number of cells and patches generated.
    */
   tbox::Array<boost::shared_ptr<tbox::Statistic> > d_boxes_stat;
   tbox::Array<boost::shared_ptr<tbox::Statistic> > d_cells_stat;
   tbox::Array<boost::shared_ptr<tbox::Statistic> > d_timestamp_stat;
#endif

   // The following are not yet implemented:
   GriddingAlgorithm(
      const GriddingAlgorithm&);
   void
   operator = (
      const GriddingAlgorithm&);

   // Verbose flags.
   bool d_check_overflow_nesting;
   bool d_check_proper_nesting;
   bool d_check_connectors;
   bool d_print_steps;

   /*
    * Static initialization and cleanup handler.
    */
   static tbox::StartupShutdownManager::Handler
      s_initialize_handler;

   /*
    *
    */
   static tbox::StartupShutdownManager::Handler
      s_startup_shutdown_handler;

};

}
}

#endif
