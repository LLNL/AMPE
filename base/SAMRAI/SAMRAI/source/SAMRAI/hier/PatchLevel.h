/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   A collection of patches at one level of the AMR hierarchy
 *
 ************************************************************************/

#ifndef included_hier_PatchLevel
#define included_hier_PatchLevel

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/BoxContainerSingleBlockIterator.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/hier/PatchFactory.h"
#include "SAMRAI/hier/ProcessorMapping.h"
#include "SAMRAI/tbox/Utilities.h"

#include <boost/shared_ptr.hpp>
#include <map>

namespace SAMRAI {
namespace hier {

class BaseGridGeometry;

/*!
 * @brief Container class for patches defined at a single level of the
 * AMR hierarchy.
 *
 * The patches in a patch level are distributed across the processors
 * of a parallel machine, so not all patches reside on the local processor.
 * (However, each patch is assigned to one and only one processor.)
 *
 * To iterate over the local patches in a patch level, use the patch
 * level iterator class (PatchLevel::Iterator).
 *
 * @see hier::BasePatchLevel
 * @see hier::Patch
 * @see hier::PatchDescriptor
 * @see hier::PatchFactory
 * @see hier::PatchLevelFactory
 * @see hier::PatchLevel::Iterator
 */

class PatchLevel
{
public:
   /*!
    * @brief Default constructor.  PatchLevel must be initialized before it can
    * be used.
    */
   explicit PatchLevel(
      const tbox::Dimension& dim);

   /*!
    * @brief Construct a new patch level given a BoxLevel.
    *
    * The BoxLevel provides refinement ratio information, establishing
    * the ratio between the index space of the new level and some reference
    * level (typically level zero) in some patch hierarchy.
    *
    * The ratio information provided by the BoxLevel is also used
    * by the grid geometry instance to initialize geometry information
    * of both the level and the patches on that level.
    *
    * @par Error conditions
    * When assertion checking is active, an unrecoverable assertion results
    * if either the grid geometry pointer or patch descriptor pointer is
    * null, or if the number of boxes in the array does not match the
    * mapping array.
    *
    * @param[in]  mapped_box_level
    * @param[in]  grid_geometry
    * @param[in]  descriptor The PatchDescriptor used to allocate patch data
    *             on the local processor
    * @param[in]  factory Optional PatchFactory.  If none specified, a default
    *             (standard) patch factory will be used.
    * @param[in]  defer_boundary_box_creation Flag to indicate suppressing
    *             construction of the boundary boxes.
    *
    */
   PatchLevel(
      const BoxLevel& mapped_box_level,
      const boost::shared_ptr<BaseGridGeometry>& grid_geometry,
      const boost::shared_ptr<PatchDescriptor>& descriptor,
      const boost::shared_ptr<PatchFactory>& factory =
         boost::shared_ptr<PatchFactory>(),
      bool defer_boundary_box_creation = false);

   /*!
    * @brief Construct a new patch level from the specified PatchLevel database.
    *
    * The box, mapping, and ratio to level zero data which are normally
    * passed in during the construction of a new patch level are
    * retrieved from the specified database.  The component_selector
    * argument specifies which patch data components should be allocated
    * and read in from the level_database.  By default, all bits in the
    * component selector are set to false so that no patch data are
    * allocated.
    *
    * @par Error conditions
    * When assertion checking is turned on, the level_database,
    * grid_geometry, and descriptor are checked to make sure that
    * they are not null.  If null, an unrecoverable assertion will result.
    *
    * @param[in]  level_database
    * @param[in]  grid_geometry
    * @param[in]  descriptor The PatchDescriptor used to allocate patch
    *             data.
    * @param[in]  factory
    * @param[in]  component_selector Optional ComponentSelector. @b Default:
    *             a ComponentSelector with all elements set to false
    * @param[in]  defer_boundary_box_creation Flag to indicate suppressing
    *             construction of the boundary boxes.  @b Default: false
    */
   PatchLevel(
      const boost::shared_ptr<tbox::Database>& level_database,
      const boost::shared_ptr<BaseGridGeometry>& grid_geometry,
      const boost::shared_ptr<PatchDescriptor>& descriptor,
      const boost::shared_ptr<PatchFactory>& factory,
      const ComponentSelector& component_selector =
         *(new ComponentSelector(false)),
      bool defer_boundary_box_creation = false);

   /*!
    * @brief The virtual destructor for patch level deallocates all patches.
    */
   virtual ~PatchLevel();

   /*!
    * @brief Get the level number
    *
    * @return the number of this level in a hierarchy, or the number of
    * a hierarchy level matching the index space of this level. If this
    * level does not align with the index space of a level in the hierarchy,
    * then this value is -1.  When the level is in a hierarchy, the return
    * value of the number of the level in the hierarchy.
    *
    * @see inHierarchy()
    */
   int
   getLevelNumber() const
   {
      return d_level_number;
   }

   /*!
    * @brief Set the number of this level to the level in the hierarchy
    * aligning with the index space of this level.
    *
    * The default value is -1 meaning the level index space does not align
    * with that of any hierarchy level.
    *
    * @param[in]  level
    */
   void
   setLevelNumber(
      const int level)
   {
      d_level_number = level;
      for (Iterator p(begin()); p != end(); p++) {
         p->setPatchLevelNumber(d_level_number);
      }
   }

   /*!
    * @brief Convenience method to get the next coarser level
    * number in a hierarchy.
    *
    * Used for data interpolation from coarser levels.  If the
    * level is in a hierarchy, then this value is getLevelNumber() - 1.
    *
    * @see inHierarchy()
    *
    * @return The next coarser level in the hierarchy or -1 if the level
    * does not exist in the hierarchy.
    */
   int
   getNextCoarserHierarchyLevelNumber() const
   {
      return d_next_coarser_level_number;
   }

   /*!
    * @brief Convenience method to set the number of of the next coarser
    * level in a hierarchy.
    *
    * For the purposes of data interpolation from coarser levels, set the
    * next coarser level in a hierarchy.  The default of -1 means the
    * level does not relate to any hierarchy.
    *
    * @param[in]  level
    */
   void
   setNextCoarserHierarchyLevelNumber(
      const int level)
   {
      d_next_coarser_level_number = level;
   }

   /*!
    * @brief Determine if this level resides in a hierarchy.
    *
    * @return true if this level resides in a hierarchy, otherwise false.
    */
   bool
   inHierarchy() const
   {
      return d_in_hierarchy;
   }

   /*!
    * @brief Setting to indicate whether this level resides in a hierarchy.
    *
    * @param[in]  in_hierarchy Flag to indicate whether this level resides
    *             in a hierarchy.  @b Default: false
    */
   void
   setLevelInHierarchy(
      bool in_hierarchy)
   {
      d_in_hierarchy = in_hierarchy;
      for (Iterator p(begin()); p != end(); p++) {
         p->setPatchInHierarchy(d_in_hierarchy);
      }
   }

   /*!
    * @brief Get the number of patches.
    *
    * This is equivalent to calling PatchLevel::getGlobalNumberOfPatches().
    */
   int
   getNumberOfPatches() const
   {
      return getGlobalNumberOfPatches();
   }

   /*!
    * @brief Get the local number of patches
    */
   int
   getLocalNumberOfPatches() const
   {
      return static_cast<int>(d_mapped_box_level->getLocalNumberOfBoxes());
   }

   /*!
    * @brief Get the global number of patches.
    */
   int
   getGlobalNumberOfPatches() const
   {
      return d_mapped_box_level->getGlobalNumberOfBoxes();
   }

   /*!
    * @brief Get the local number of Cells.
    */
   int
   getLocalNumberOfCells() const
   {
      return static_cast<int>(d_mapped_box_level->getLocalNumberOfCells());
   }

   /*!
    * @brief Get the global number of cells
    */
   int
   getGlobalNumberOfCells() const
   {
      return d_mapped_box_level->getGlobalNumberOfCells();
   }

   /*!
    * @brief Get a Patch based on its GlobalId.
    *
    * @param[in]  gid
    *
    * @return A boost::shared_ptr to the Patch indicated by the GlobalId.
    */
   const boost::shared_ptr<Patch>&
   getPatch(
      const GlobalId& gid) const
   {
      BoxId mbid(gid);
      PatchContainer::const_iterator it = d_patches.find(mbid);
      TBOX_ASSERT(it != d_patches.end());
      return it->second;
   }

   /*!
    * @brief Get a Patch based on its BoxId.
    *
    * @param[in]  mbid
    *
    * @return A boost::shared_ptr to the Patch indicated by the BoxId.
    */
   boost::shared_ptr<Patch>
   getPatch(
      const BoxId& mbid) const
   {
      const PatchContainer::const_iterator mi = d_patches.find(mbid);
#ifdef DEBUG_CHECK_ASSERTIONS
      if (mi == d_patches.end()) {
         TBOX_ERROR("PatchLevel::getPatch(" << mbid
            << "): patch does not exist locally.");
      }
#endif
      return (*mi).second;
   }

   /*!
    * @brief Get the PatchDescriptor
    *
    * @return pointer to the patch descriptor for the hierarchy.
    */
   boost::shared_ptr<PatchDescriptor>
   getPatchDescriptor() const
   {
      return d_descriptor;
   }

   /*!
    * @brief Get the PatchFactory
    *
    * @return the factory object used to created patches in the level.
    */
   boost::shared_ptr<PatchFactory>
   getPatchFactory() const
   {
      return d_factory;
   }

   /*!
    * @brief Get the grid geometry
    *
    * @return A boost::shared_ptr to the grid geometry description.
    */
   boost::shared_ptr<BaseGridGeometry>
   getGridGeometry() const
   {
      return d_geometry;
   }

   /*!
    * @brief Update this patch level through refining.
    *
    * The data members of this patch level are updated by refining the
    * information on a given coarse level using the given ratio between
    * the two levels.  The fine level will cover the same physical space as
    * the coarse level and will have the same number of patches with the
    * same mapping of those patches to processors.  However, the index
    * space of the level will be refined by the specified ratio.
    *
    * @par Assumptions
    * If the fine grid geometry is null (default case), then it is assumed
    * that this level is to use the same grid geometry as the given coarse
    * level and the ratio to level zero is set relative to the given coarse
    * level.  Otherwise, we use the given grid geometry (assumed to be a proper
    * refinement of the grid geometry used on the given coarse level) and copy
    * ratio to level zero from given coarse level.  In other words, the function
    * can be used to produce two different results.
    *
    * <ol>
    *   <li> When passed a null grid geometry pointer, the refined patch level
    *        can be used for data exchange operations with the AMR hierarchy
    *        in which the coarse level resides -- both levels are defined with
    *        respect to the index space of the grid geometry object which they
    *        share.  Thus, the refined patch level can be used in data
    *        exchanges with the AMR hierarchy of the coarse level
    *        automatically.
    *   <li> Second, when passed a non-null fine grid geometry pointer, the
    *        level is defined relative to that geometry and the refined patch
    *        level cannot be used in data exchanges with the AMR hierarchy
    *        of the coarse level automatically in general.  This mode is
    *        used to construct a refined copy of an entire patch hierarchy,
    *        typically.
    * </ol>
    *
    * @param[in]  coarse_level
    * @param[in]  refine_ratio
    * @param[in]  fine_grid_geometry @b Default: boost::shared_ptr to a null
    *             grid geometry
    * @param[in]  defer_boundary_box_creation @b Default: false
    */
   void
   setRefinedPatchLevel(
      const boost::shared_ptr<PatchLevel>& coarse_level,
      const IntVector& refine_ratio,
      const boost::shared_ptr<BaseGridGeometry>& fine_grid_geometry =
         boost::shared_ptr<BaseGridGeometry>(),
      bool defer_boundary_box_creation = false);

   /*!
    * @brief Update this patch through coarsening.
    *
    * The data members of this patch level are updated by coarsening the
    * information on a given fine level using the given ratio between
    * the two levels.  The coarse level will cover the same physical space as
    * the fine level and will have the patches with the same
    * GlobalIndices.  However, the index space of the level will be coarsened
    * by the specified ratio.
    * @par Assumptions
    * If the coarse grid geometry is null (default case), then it is assumed
    * that this level is to use the same grid geometry as the given fine
    * level and the ratio to level zero is set relative to the given fine
    * level.  Otherwise, we use the given grid geometry (assumed to be a proper
    * coarsening of the grid geometry used on the given fine level) and copy
    * ratio to level zero from given fine level.  In other words, the function
    * can be used to produce two different results.
    *
    * <ol>
    *   <li> When passed a null grid geometry pointer, the coarsened
    *        patch level can be used for data exchange operations with the
    *        AMR hierarchy in which the fine level resides -- both levels
    *        are defined with respect to the index space of the grid geometry
    *        object which they share.  Thus, the coarsened patch level can be
    *        used in data exchanges with the AMR hierarchy of the fine level
    *        automatically.
    *   <li> When passed a non-null coarse grid geometry pointer, the level is
    *        defined relative to that geometry and the coarsened patch level
    *        cannot be used in data exchanges with the AMR hierarchy of the
    *        fine level automatically in general.  This mode is used to
    *        construct a coarsened copy of an entire patch hierarchy,
    *        typically.
    * </ol>
    *
    * @param[in]  fine_level
    * @param[in]  coarsen_ratio
    * @param[in]  coarse_grid_geom @b Default: boost::shared_ptr to a null
    *             grid geometry
    * @param[in]  defer_boundary_box_creation @b Default: false
    */
   void
   setCoarsenedPatchLevel(
      const boost::shared_ptr<PatchLevel>& fine_level,
      const IntVector& coarsen_ratio,
      const boost::shared_ptr<BaseGridGeometry>& coarse_grid_geom =
         boost::shared_ptr<BaseGridGeometry>(),
      bool defer_boundary_box_creation = false);

   /*!
    * @brief Create and store the boundary boxes for this level.
    *
    * If boundary boxes have already been constructed, this function
    * does nothing.
    * @note
    * If the level is constructed with boundary box creation deferred,
    * this method must be called before any attempt at filling data at
    * physical boundaries.  This function is called from
    * xfer::RefineSchedule prior to any physical boundary operations.
    */
   void
   setBoundaryBoxes()
   {
      if (!d_boundary_boxes_created) {
         d_geometry->setBoundaryBoxes(*this);
         d_boundary_boxes_created = true;
      }
   }

   /*!
    * @brief Get the physical domain.
    *
    * @return A const reference to the box array that defines
    * the extent of the index space on the level.
    */
   const tbox::Array<BoxContainer>&
   getPhysicalDomainArray() const
   {
      return d_physical_domain;
   }

   const BoxContainer&
   getPhysicalDomain(
      const BlockId& block_id) const
   {
      return d_physical_domain[block_id.getBlockValue()];
   }

   /*!
    * @brief Get the box defining the patches on the level.
    *
    * The internal state of PatchLevel (where boxes are concern) is
    * dependent on the BoxLevel associated with it, and computed
    * only if getBoxes() is called.  The first call to getBoxes() must be
    * done by all processors as it requires communication.
    *
    * @return a const reference to the box array that defines
    * the patches on the level.
    */
   const BoxContainer&
   getBoxes() const
   {
      if (!d_has_globalized_data) {
         initializeGlobalizedBoxLevel();
      }
      return d_boxes;
   }

   /*!
    * @brief Get boxes for a particular block.
    *
    * The boxes from only the specified block will be stored in the
    * output BoxContainer.
    *
    * @param[out] boxes
    * @param[in] block_id
    */
   void
   getBoxes(
      BoxContainer& boxes,
      const BlockId& block_id) const;

   /*!
    * @brief Get the BoxLevel associated with the PatchLevel.
    *
    * @return a reference to a boost::shared_ptr to the BoxLevel
    * associated with the PatchLevel.
    */
   const boost::shared_ptr<BoxLevel>&
   getBoxLevel() const
   {
      return d_mapped_box_level;
   }

   /*!
    * @brief Get the globalized version of the BoxLevel associated
    * with the PatchLevel.
    * @note
    * The first time this method is used, a global communication is
    * done.  Thus all processors must use this method the first time
    * any processor uses it.
    *
    * @return The globalized version of the BoxLevel associated
    * with the PatchLevel.
    */
   const BoxLevel&
   getGlobalizedBoxLevel() const
   {
      if (!d_has_globalized_data) {
         initializeGlobalizedBoxLevel();
      }
      return d_mapped_box_level->getGlobalizedVersion();
   }

   /*!
    * @brief Get the mapping of patches to processors.
    *
    * @return A const reference to the mapping of patches to processors.
    */
   const ProcessorMapping&
   getProcessorMapping() const
   {
      if (!d_has_globalized_data) {
         initializeGlobalizedBoxLevel();
      }
      return d_mapping;
   }

   /*!
    * @brief Get the ratio between the index space of this PatchLevel and
    * the reference level in the AMR hierarchy.
    *
    * @return A const reference to the vector ratio between the index
    * space of this patch level and that of a reference level in AMR
    * hierarchy (that is, level zero).
    */
   const IntVector&
   getRatioToLevelZero() const
   {
      return d_ratio_to_level_zero;
   }

   /*!
    * @brief Get the ratio between this level and the next coarser
    * level in the patch hierarchy.
    *
    * This vector is set with the setRatioToCoarserLevel() function.
    * If the level is not in a hierarchy, a default ratio of zero is returned.
    *
    * @return the vector ratio between this level and the next coarser
    * level in the patch hierarchy.
    */
   const IntVector&
   getRatioToCoarserLevel() const
   {
      return d_ratio_to_coarser_level;
   }

   /*!
    * @brief Set the ratio between this level and the next coarser
    * level in the patch hierarchy.
    *
    * This is required only when level resides in a hierarchy.
    *
    * @param[in] ratio
    */
   void
   setRatioToCoarserLevel(
      const IntVector& ratio)
   {
      d_ratio_to_coarser_level = ratio;
   }

   /*!
    * @brief Get the processor mapping for the patch.
    *
    * @return the processor that owns the specified patch.  The patches
    * are numbered starting at zero.
    *
    * @param[in] mapped_box_id Patch's BoxId
    */
   int
   getMappingForPatch(
      const BoxId& mapped_box_id) const
   {
      // Note: p is required to be a local index.
      /*
       * This must be for backward compatability, because if p is a local
       * index, the mapping is always to d_mapped_box_level->getRank().
       * Here is the old code:
       *
       * return d_mapped_box_level->getBoxStrict(p)->getOwnerRank();
       */
      NULL_USE(mapped_box_id);
      return d_mapped_box_level->getMPI().getRank();
   }

   /*!
    * @brief Get the box for the specified patch
    *
    * @return The box for the specified patch.
    *
    * @param[in] mapped_box_id Patch's BoxId
    */
   const Box&
   getBoxForPatch(
      const BoxId& mapped_box_id) const
   {
      TBOX_ASSERT(mapped_box_id.getOwnerRank() ==
                  d_mapped_box_level->getMPI().getRank());
      return getPatch(mapped_box_id)->getBox();
   }

   /*!
    * @brief Determine if the patch is adjacent to a non-periodic
    * physical domain boundary.
    *
    * @param[in] mapped_box_id Patch's BoxId
    *
    * @return True if patch with given number is adjacent to a non-periodic
    * physical domain boundary.  Otherwise, false.
    */
   bool
   patchTouchesRegularBoundary(
      const BoxId& mapped_box_id) const
   {
      TBOX_ASSERT(mapped_box_id.getOwnerRank() ==
                  d_mapped_box_level->getMPI().getRank());
      return getPatch(mapped_box_id)->getPatchGeometry()->getTouchesRegularBoundary();
   }

   /*!
    * @brief Allocate the specified component on all patches.
    *
    * @param[in]  id
    * @param[in]  timestamp @b Default: zero (0.0)
    */
   void
   allocatePatchData(
      const int id,
      const double timestamp = 0.0)
   {
      for (Iterator ip(begin()); ip != end(); ip++) {
         ip->allocatePatchData(id, timestamp);
      }
   }

   /*!
    * @brief Allocate the specified components on all patches.
    *
    * @param[in]  components The componentSelector indicating
    *             which elements to allocate
    * @param[in]  timestamp @b Default: zero (0.0)
    */
   void
   allocatePatchData(
      const ComponentSelector& components,
      const double timestamp = 0.0)
   {
      for (Iterator ip(begin()); ip != end(); ip++) {
         ip->allocatePatchData(components, timestamp);
      }
   }

   /*!
    * @brief Determine if the patch data has been allocated.
    *
    * @return True if (1) there are no patches in this patch level or  if
    *         (2) all of the patches have allocated the patch data component,
    *         otherwise false.
    *
    * @param[in] id The patch identifier.
    */
   bool
   checkAllocated(
      const int id) const
   {
      bool allocated = true;
      for (PatchContainer::const_iterator mi = d_patches.begin();
           mi != d_patches.end(); ++mi) {
         allocated &= (*mi).second->checkAllocated(id);
      }
      return allocated;
   }

   /*!
    * @brief  Deallocate the specified component on all patches.
    *
    * This component will need to be reallocated before its next use.
    *
    * @param[in]  id The patch identifier
    */
   void
   deallocatePatchData(
      const int id)
   {
      for (Iterator ip(begin()); ip != end(); ip++) {
         ip->deallocatePatchData(id);
      }
   }

   /*!
    * @brief Deallocate the specified components on all patches.
    *
    * Components will need to be reallocated before their next use.
    *
    * @param[in]  components The ComponentSelector indicating which
    *             components to deallocate.
    */
   void
   deallocatePatchData(
      const ComponentSelector& components)
   {
      for (Iterator ip(begin()); ip != end(); ip++) {
         ip->deallocatePatchData(components);
      }
   }

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
    * @brief Set the simulation time for the specified patch component.
    *
    * @param[in]  timestamp
    * @param[in]  id The patch identifier
    */
   void
   setTime(
      const double timestamp,
      const int id)
   {
      for (Iterator ip(begin()); ip != end(); ip++) {
         ip->setTime(timestamp, id);
      }
   }

   /*!
    * @brief Set the simulation time for the specified patch components.
    *
    * @param[in] timestamp
    * @param[in] components The ComponentSelector indicating on which
    *            components to set the simulation time.
    */
   void
   setTime(
      const double timestamp,
      const ComponentSelector& components)
   {
      for (Iterator ip(begin()); ip != end(); ip++) {
         ip->setTime(timestamp, components);
      }
   }

   /*!
    * @brief Set the simulation time for all allocated patch components.
    *
    * @param[in]  timestamp
    */
   void
   setTime(
      const double timestamp)
   {
      for (Iterator ip(begin()); ip != end(); ip++) {
         ip->setTime(timestamp);
      }
   }

   /*!
    * @brief Use the PatchLevel database to set the state of the PatchLevel
    * and to create all patches on the local processor.
    *
    * @par Assertions
    * Assertions will check that database is a non-null boost::shared_ptr,
    * that the data being retrieved from the database are of
    * the type expected.  Also checked is the number of patches is positive,
    * and the number of patches and size of processor mapping array are the
    * same, and that the number of patches and the number of boxes on the
    * level are equal.
    *
    * @param[in,out] database
    * @param[in]     component_selector
    */
   void
   getFromDatabase(
      const boost::shared_ptr<tbox::Database>& database,
      const ComponentSelector& component_selector);

   /*!
    * @brief Write data to the database.
    *
    * Writes the data from the PatchLevel to the database.
    * Also tells all local patches to write out their state to
    * the database.
    *
    * @par Assertions
    * Check that database is a non-null boost::shared_ptr.
    *
    * @param[in,out]  database
    * @param[in]      patchdata_write_table The ComponentSelector specifying
    *                 which patch data to write to the database
    */
   void
   putUnregisteredToDatabase(
      const boost::shared_ptr<tbox::Database>& database,
      const ComponentSelector& patchdata_write_table) const;

   /*!
    * @brief Print a patch level to varying details.
    *
    * If depth>0, print function will be called for each patch in the level.
    *
    * @param[in]  os The std::ostream in which to print to
    * @param[in]  border @b Default: empty string
    * @param[in]  depth @b Default: zero (0).
    *
    * @return 0.  Always.
    */
   int
   recursivePrint(
      std::ostream& os,
      const std::string& border = std::string(),
      int depth = 0);

private:
   /*
    * Static integer constant describing class's version number.
    */
   static const int HIER_PATCH_LEVEL_VERSION;

   /*
    * @brief Container of distributed patches on level.
    *
    */
   typedef std::map<BoxId, boost::shared_ptr<Patch> > PatchContainer;

public:
   /*!
    * @brief Iterator for looping through local patches.
    */
   class Iterator
   {
friend class PatchLevel;

public:
      /*!
       * @brief Copy constructor.
       *
       * @param[in]  other
       */
      Iterator(
         const Iterator& other);

      /*!
       * @brief Assignment operator
       */
      Iterator&
      operator = (
         const Iterator& rhs)
      {
         d_iterator = rhs.d_iterator;
         d_patches = rhs.d_patches;
         return *this;
      }

      /*!
       * @brief Dereference operator.
       */
      const boost::shared_ptr<Patch>&
      operator * () const
      {
         return d_iterator->second;
      }

      /*!
       * @brief Delegation operations to the Patch pointer.
       */
      const boost::shared_ptr<Patch>&
      operator -> () const
      {
         return d_iterator->second;
      }

      /*!
       * @brief Equality comparison.
       */
      bool
      operator == (
         const Iterator& rhs) const
      {
         return d_iterator == rhs.d_iterator;
      }

      /*!
       * @brief Inequality operator.
       */
      bool
      operator != (
         const Iterator& rhs) const
      {
         return d_iterator != rhs.d_iterator;
      }

      /*!
       * @brief Pre-increment.
       */
      Iterator&
      operator ++ ()
      {
         ++d_iterator;
         return *this;
      }

      /*!
       * @brief Post-increment.
       */
      Iterator
      operator ++ (
         int)
      {
         Iterator tmp_iterator = *this;
         ++d_iterator;
         return tmp_iterator;
      }

private:
      /*
       * Unimplemented default constructor.
       */
      Iterator();

      /*!
       * @brief Construct from a PatchLevel.
       *
       * @param[in]  patch_level
       * @param[in]  begin
       */
      explicit Iterator(
         const PatchLevel* patch_level,
         bool begin);

      /*!
       * @brief The real iterator (this class is basically a wrapper).
       */
      PatchContainer::const_iterator d_iterator;

      /*!
       * @brief For supporting backward-compatible interface.
       */
      const PatchContainer* d_patches;

   };

   /*!
    * @brief Construct an iterator pointing to the first Patch in the
    * PatchLevel.
    */
   Iterator
   begin() const
   {
      return Iterator(this, true);
   }

   /*!
    * @brief Construct an iterator pointing to the last Patch in the
    * PatchLevel.
    */
   Iterator
   end() const
   {
      return Iterator(this, false);
   }

   /*!
    * @brief Typdef PatchLevel::Iterator to standard iterator nomenclature.
    */
   typedef Iterator iterator;

private:
   /**
    * @brief Static initialization to be done at startup.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /**
    * @brief Static cleanup to be done at shutdown.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   void
   initializeGlobalizedBoxLevel() const;

   /*!
    * @brief Dimension of the object
    */
   const tbox::Dimension d_dim;

   /*!
    * @brief Number of blocks that can be represented by this level.
    */
   int d_number_blocks;

   /*!
    * Primary metadata describing the PatchLevel.
    */
   boost::shared_ptr<BoxLevel> d_mapped_box_level;

   /*
    * Whether we have a globalized version of d_mapped_box_level.
    */
   mutable bool d_has_globalized_data;
   /*
    * Boxes for all level patches.
    *
    * d_boxes is slave to d_mapped_box_level and computed only if getBoxes() is called.
    * This means that the first getBoxes() has to be called by all processors,
    * because it requires communication.
    */
   mutable BoxContainer d_boxes;

   /*
    * Patch mapping to processors.
    */
   mutable ProcessorMapping d_mapping;

   /*
    * ratio to reference level
    */
   IntVector d_ratio_to_level_zero;

   /*
    * Grid geometry description.
    */
   boost::shared_ptr<BaseGridGeometry> d_geometry;
   /*
    * PatchDescriptor - patch data info shared by all patches in the hierarchy
    */
   boost::shared_ptr<PatchDescriptor> d_descriptor;
   /*
    * Factory for creating patches.
    */
   boost::shared_ptr<PatchFactory> d_factory;

   /*
    * Local number of patches on the level.
    */
   int d_local_number_patches;

   /*
    * Extent of the index space.
    */
   tbox::Array<BoxContainer> d_physical_domain;

   /*
    * The ratio to coarser level applies only when the level resides
    * in a hierarchy.  The level number is that of the hierarchy level
    * that aligns with the index space of the level; if level aligns with
    * no such level then the value is -1 (default value).  The next coarser
    * level number is the next coarser level in the hierarchy for the
    * purposes of filling data from coarser levels.   It is -1 by default
    * but is usually a valid level number more often than level number.
    * The boolean is true when the level is in a hierarchy, false otherwise.
    */
   IntVector d_ratio_to_coarser_level;

   /*
    * Level number in the hierarchy.
    */
   int d_level_number;

   /*
    * Aligning with the index space of the next coarser level number.
    */
   int d_next_coarser_level_number;

   /*
    * Flag indicating the level is in a hierarchy.
    */
   bool d_in_hierarchy;

   /*
    * Container for patches.
    */
   PatchContainer d_patches;

   /*
    * Flag to indicate boundary boxes are created.
    */
   bool d_boundary_boxes_created;

   /*!
    * @brief Has shutdown handler been initialized.
    *
    * This should be checked and set in every ctor.
    */
   static bool s_initialized;

   /*!
    * @brief Initialize static state
    */
   static bool
   initialize();

   static tbox::StartupShutdownManager::Handler
      s_initialize_finalize_handler;
};

}
}

#endif
