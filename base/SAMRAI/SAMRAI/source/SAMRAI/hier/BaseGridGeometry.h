/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Base class for geometry management in AMR hierarchy
 *
 ************************************************************************/

#ifndef included_hier_BaseGridGeometry
#define included_hier_BaseGridGeometry

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/Serializable.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchGeometry.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/hier/TimeInterpolateOperator.h"
#include "SAMRAI/hier/TransferOperatorRegistry.h"
#include "SAMRAI/hier/Transformation.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <list>

namespace SAMRAI {
namespace hier {

class BoxLevel;
class PatchLevel;
class BoxTree;

/*!
 * @brief Class BaseGridGeometry manages the index space that determines the
 * extent of the coarse-level domain of a SAMRAI hierarchy.
 *
 * A BaseGridGeometry object may be used as a base class to derive child
 * classes that manage particular grid types (%e.g., Cartesian, cylindrical,
 * etc.).
 *
 * The grid geometry class is responsible for maintaining information
 * about the index space describing the physical domain and computing this
 * information for patches in an AMR hierarchy.
 *
 * Required input file keys and data types (only required when using first
 * constructor which takes an input database):
 *
 *    - @b    num_blocks
 *       Integer value specifying the number of blocks in the multiblock
 *       mesh configuration.
 *
 *    - @b    domain_boxes_
 *       For each block, an array of boxes representing the index space for the
 *       entire domain within a block (on the coarsest mesh level; i.e., level
 *       zero).  The key must have an integer value as a suffix
 *       (domain_boxes_0, domain_boxes_1, etc.), and there must be an entry
 *       for every block from 0 to num_blocks-1.
 *
 * Optional input file keys and data types (for first constructor which
 * takes an input database):
 *
 *    - @b    periodic_dimension
 *       An array of integer values (expected number of values is equal to
 *       the spatial dimension of the mesh) representing the directions in
 *       which the physical domain is periodic.  A non-zero value indicates
 *       that the direction is periodic.  A zero value indicates that
 *       the direction is not periodic.  If no values are specified, then
 *       the array is initialized to all zeros (no periodic directions).
 *       This key should only be used when the number of blocks is 1 (a
 *       single block mesh), as periodic boundaries are not supported for
 *       multiblock meshes.
 *
 *    - @b   BlockNeighbors
 *
 *       For multiblock grids, a BlockNeighbors entry must be given for
 *       every pair of blocks that touch each other in any way.  The key
 *       for this entry must include a unique trailing integer, and the
 *       integers for the full set of BlockNeighbors keys must be a
 *       continuous sequence beginning with 0.
 *
 *    - @b  Singularity
 *
 *       When there is a reduced or enhanced connectivity singularity,
 *       this key must be used to identify which blocks touch the
 *       singularity and the position of the singularity in relation to
 *       each block's index space.  Like BlockNeighbors, each entry must
 *       have a trailing integer beginning with 0.
 *
 * A description of the input format for BlockNeighbors* and Singularity
 * is included in the Multiblock.pdf document in the docs/userdocs
 * directory of the SAMRAI distribution.
 *
 * @par Additional Functionality
 * Operations performed by this class include determining which patches are
 * adjacent to the physical domain boundary and computing boundary boxes
 * for patches which describe how the patch touches the domain boundary
 * (useful for filling ghost cell data for physical boundary conditions).
 *
 * @see hier::BoundaryBox
 */

class BaseGridGeometry:
   public tbox::Serializable
{
   friend class TransferOperatorRegistry;

public:
   typedef  PatchGeometry::TwoDimBool TwoDimBool;

   /*!
    * @brief Construct a new BaseBaseGridGeometry object and initialize from
    * input.
    *
    * This constructor for BaseGridGeometry initializes data members
    * based on parameters read from the specified input database.
    * The constructor also registers this object for restart using
    * the specified object name, when the boolean argument is true.
    * Whether the object will write its state to restart files during
    * program execution is determined by this argument.
    *
    * This constructor is intended for use when directly constructing a
    * BaseGridGeometry without using a derived child class.  The object will
    * contain all index space grid information for a mesh, but nothing about
    * the physical coordinates of the mesh.
    *
    * @note
    * @b Errors: passing in a null database pointer or an empty string
    * will result in an unrecoverable assertion.
    *
    * @param[in]  dim
    * @param[in]  object_name
    * @param[in]  input_db
    * @param[in]  register_for_restart Flag indicating whether this instance
    *             should be registered for restart.  @b Default: true
    */
   BaseGridGeometry(
      const tbox::Dimension& dim,
      const std::string& object_name,
      const boost::shared_ptr<tbox::Database>& input_db,
      bool register_for_restart = true);

   /*!
    * @brief Construct a new BaseGridGeometry object based on arguments.
    *
    * This constructor creates a new BaseGridGeometry object based on the
    * arguments, rather than relying on input or restart data.  The
    * constructor also registers this object for restart using
    * the specified object name, when the boolean argument is true.
    * Whether the object will write its state to restart files during
    * program execution is determined by this argument.
    *
    * @param[in]  object_name
    * @param[in]  domain      Each element of the array describes the index
    *                         space for a block.
    * @param[in]  register_for_restart Flag indicating whether this instance
    *             should be registered for restart.  @b Default: true
    */
   BaseGridGeometry(
      const std::string& object_name,
      const BoxContainer& domain,
      bool register_for_restart = true);

   /*!
    * @brief Virtual destructor
    */
   virtual ~BaseGridGeometry();

   //! @{
   /*!
    * @name Functions for computing boundary boxes
    */
   /*!
    * @brief For every patch on a level, find all ways a patch touches
    * a physical or periodic boundary.
    *
    * For every patch on the given PatchLevel, this routine determines which
    * kinds of boundaries are touched (regular, periodic, both, or neither).
    *
    * @param[out]   touches_regular_boundary map to store which patches touch
    *               non-periodic boundaries.
    * @param[out]   touches_periodic_boundary map to store which patches touch
    *               periodic boundaries.
    * @param[in]    level containing the patches to be checked
    */
   void
   findPatchesTouchingBoundaries(
      std::map<BoxId, TwoDimBool>& touches_regular_boundary,
      std::map<BoxId, TwoDimBool>& touches_periodic_boundary,
      const PatchLevel& level) const;

   /*!
    * @brief Version of findPatchTouchingBoundaries for a single box.
    *
    * @param[out]   touches_regular_bdry TwoDimBool to store which patches touch
    *               non-periodic boundaries.
    * @param[out]   touches_periodic_bdry TwoDimBool to store which patches touch
    *               periodic boundaries.
    * @param[in]    box to be checked
    *
    * @param[in] refinement_ratio Refinement ratio of given box
    *
    * @param[in] refined_periodic_domain_tree
    */
   void
   computeBoxTouchingBoundaries(
      TwoDimBool& touches_regular_bdry,
      TwoDimBool& touches_periodic_bdry,
      const Box& box,
      const IntVector &refinement_ratio,
      const BoxContainer& refined_periodic_domain_tree) const;

   /*!
    * @brief Sets geometry data for patches on a level.
    *
    * Using the boundary information previously computed, this method
    * will pass the information to the concrete implementation of the
    * geometry class, and construct boundary boxes if required.
    *
    * @param[in]   level containing the patches to be checked.
    * @param[in]   ratio_to_level_zero ratio to the coarsest level.
    * @param[in]   touches_regular_bdry Array storing which patches touch
    *              non-periodic boundaries.
    * @param[in]   touches_periodic_bdry Array storing which patches touch
    *              periodic boundaries.
    * @param[in]   defer_boundary_box_creation Flag to indicate if boundary
    *              boxes should be created
    */
   /*
    * TODO:  Are all these parameters really only input params?  If so, these
    * should be const references rather than non-const.
    *
    * TODO:  The confusing coordination and sequence of calls between
    *        grid geometry classes and patch geometry classes should
    *        be reworked, including simplification of the data we're
    *        passing around.  Part of the issue is tied to the fact that
    *        we are inheriting implementation here as well as providing
    *        a concrete implementation of an interface.
    */
   virtual void
   setGeometryOnPatches(
      PatchLevel& level,
      const IntVector& ratio_to_level_zero,
      std::map<BoxId, TwoDimBool>& touches_regular_bdry,
      std::map<BoxId, TwoDimBool>& touches_periodic_bdry,
      bool defer_boundary_box_creation);

   /*!
    * @brief Construct and set the boundary boxes for each patch.
    *
    * Construct the boundary boxes for every patch in the level.
    * Once constructed, the boundary boxes are set on each patch's
    * PatchGeometry object.
    *
    * @param[in] level The level for which boundary boxes are constructed.
    */
   /*
    * TODO:  really input param?  If so, should be const qualified.
    *
    * TODO:  See the second TODO item for the previous method.
    */
   void
   setBoundaryBoxes(
      PatchLevel& level);

   //@}

   /*!
    * @brief Compute the box array describing the index space for a
    * given block of the physical domain.
    *
    * The domain description does not include periodic images.
    *
    * If any entry of the ratio vector is negative, the index space is
    * coarsened with respect to the physical domain description.
    * Otherwise, the index space is refined.
    *
    * @param[out]    domain The BoxContainer to be computed
    * @param[in]     ratio_to_level_zero ratio to the coarsest level
    * @param[in]     block_id
    */
   void
   computePhysicalDomain(
      BoxContainer& domain,
      const IntVector& ratio_to_level_zero,
      const BlockId& block_id) const;

   void
   computePhysicalDomainWithPeriodicImages(
      BoxContainer& domain,
      const IntVector& ratio_to_level_zero,
      const BlockId& block_id) const;

   /*!
    * @brief Compute the BoxLevel describing the index space for a
    * given block of the physical domain.
    *
    * The domain description includes periodic images, if any exist.
    *
    * Unlike the BoxContainer version of this function, the domain computed
    * contains periodic image boxes.  If any entry of ratio vector is
    * negative, the index space is coarsened with respect to the physical
    * domain description.  Otherwise, the index space is refined.
    *
    * @param[out]    box_level The BoxLevel containing all
    *                Boxes describing the index space
    * @param[in]     ratio_to_level_zero ratio to the coarsest level
    * @param[in]     block_id
    */
   void
   computePhysicalDomain(
      BoxLevel& box_level,
      const IntVector& ratio_to_level_zero,
      const BlockId& block_id) const;

   /*!
    * @brief Compute the BoxContainer describing the complete physical
    * domain for all blocks.
    *
    * The domain description includes periodic images, if any exist.
    *
    * If any entry of the ratio vector is negative, the index space is
    * coarsened with respect to the physical domain description.
    * Otherwise, the index space is refined.
    *
    * @param[out]    domain_mapped_boxes The BoxContainer containing all
    *                Boxes describing the physical domain
    * @param[in]     ratio_to_level_zero ratio to the coarsest level
    */
   void
   computePhysicalDomain(
      BoxContainer& domain_mapped_boxes,
      const IntVector& ratio_to_level_zero) const;

   /*!
    * @brief Compute the BoxLevel describing the complete physical
    * domain for all blocks.
    *
    * The domain description includes periodic images, if any exist.
    *
    * If any entry of the ratio vector is negative, the index space is
    * coarsened with respect to the physical domain description.
    * Otherwise, the index space is refined.
    *
    * @param[out]    box_level The BoxLevel containing all
    *                Boxes describing the physical domain
    * @param[in]     ratio_to_level_zero ratio to the coarsest level
    */
   void
   computePhysicalDomain(
      BoxLevel& box_level,
      const IntVector& ratio_to_level_zero) const;

   /*!
    * @brief Set the physical domain (for level zero)
    *
    * Each entry in the array of box arrays represents the physical domain
    * for a single block
    *
    * The extents of the input domain boxes are used, but their
    * LocalId's are disregarded.  BaseGridGeometry will assign new and
    * unique LocalId's to the domain box description.  Subsequent
    * calls to getPhysicalDomain() will return boxes with the new
    * LocalId's.
    *
    * @param[in]     domain The input array of BoxContainer
    * @param[in]     number_blocks
    */
   void
   setPhysicalDomain(
      const BoxContainer& domain,
      const int number_blocks);

   /*!
    * @brief Get the physical domain description on level zero.
    *
    * @return const reference to physical domain description for level 0.
    */
   const BoxContainer&
   getPhysicalDomain() const
   {
      return d_physical_domain;
   }

   /*!
    * @brief Access the multiblock domain description as a tree
    * without periodic images.
    *
    * This tree does not contain periodic images, even if there is
    * only one block and the domain is periodic.
    *
    * @return The multiblock domain description as a search tree.
    */
   const BoxContainer&
   getDomainSearchTree() const
   {
      return d_physical_domain;
   }

   /*!
    * @brief Access the multiblock domain description with periodic
    * images (if any)
    *
    * @return The domain description as a search tree with periodic
    * images (if any).
    */
   const BoxContainer& 
   getPeriodicDomainSearchTree() const
   {
      return d_domain_with_images;
   }

   /*!
    * @brief returns whether the physical domain for a block managed by this
    * geometry object can be represented as a single box.
    *
    * @return true if the physical domain can be represented as a single box,
    *         otherwise false.
    *
    * @param[in]   block_id
    */
   bool
   getDomainIsSingleBox(
      const BlockId& block_id) const
   {
      return d_domain_is_single_box[block_id.getBlockValue()];
   }

   /*!
    * @brief Initialize the periodic shift on the coarsest level.
    *
    * @param[in]  directions an array indicating periodic directions(1) or
    *             all others (0).
    *
    * @note
    * The IntVector argument should be set to 1 for periodic directions
    * and 0 for all other directions.  The shift will be calculated to
    * the number of cells in the periodic direction and zero in all
    * other directions.
    */
   void
   initializePeriodicShift(
      const IntVector& directions);

   /*!
    * @brief Get the periodic shift in each direction for the physical domain
    * managed by this geometry object.
    *
    * The returned IntVector contains the periodic shift in each direction
    * for a domain represented by a refinement of the reference physical
    * domain (i.e. level zero) by the given ratio vector.  Entries
    * will be zero for non-periodic directions.
    *
    * @param[in]     ratio_to_level_zero ratio to the coarsest level.
    *
    * @return        The periodic shift in each direction for a domain
    *                represented by a refinement of the reference physical
    *                domain.
    */
   IntVector
   getPeriodicShift(
      const IntVector& ratio_to_level_zero) const;

   /*!
    * @brief Get the number of blocks in the geometry.
    */
   int
   getNumberBlocks() const
   {
      return d_number_blocks;
   }

   /*!
    * @brief Get the max stencil width of all transfer operators.
    *
    * The max stencil width is required by the DLBG to determine when
    * two nearby boxes are defined as neighbors.  The DLBG in turn
    * provides the neighbor information for various operations such as
    * schedule construction.
    *
    * @return The max stencil width of all transfer operators.
    */
   IntVector
   getMaxTransferOpStencilWidth()
   {
      return d_transfer_operator_registry->getMaxTransferOpStencilWidth();
   }

   /*!
    * @brief Create a pointer to a refined version of this grid geometry
    *        object.
    *
    * Virtual method -- should be overridden in specialized grid geometry
    * classes
    *
    * @param[in]     fine_geom_name std::string name of the geometry object
    * @param[in]     refine_ratio the refinement ratio.
    * @param[in]     register_for_restart Flag to indicate whether to register
    *                for restart.
    *
    * @return The pointer to the grid geometry object.
    */
   virtual boost::shared_ptr<BaseGridGeometry>
   makeRefinedGridGeometry(
      const std::string& fine_geom_name,
      const IntVector& refine_ratio,
      bool register_for_restart) const = 0;

   /*!
    * @brief Create a pointer to a coarsened version of this grid geometry
    *        object.
    *
    * Virtual method -- should be overridden in specialized grid geometry
    * classes
    *
    * @param[in]     coarse_geom_name std::string name of the geometry object
    * @param[in]     coarsen_ratio the coasening ratio
    * @param[in]     register_for_restart Flag to indicate whether to register
    *                for restart.
    *
    * @return The pointer to a coarsened version of this grid geometry object.
    */
   virtual boost::shared_ptr<BaseGridGeometry>
   makeCoarsenedGridGeometry(
      const std::string& coarse_geom_name,
      const IntVector& coarsen_ratio,
      bool register_for_restart) const = 0;

   /*!
    * @brief Compute and set grid data for patch.
    *
    * Virtual method -- should be overridden in specialized grid geometry
    * classes
    *
    * @param[in,out]    patch The patch on which to set grid data and the new
    *                   concrete patch geometry object.
    * @param[in]        ratio_to_level_zero ratio to coarsest level
    * @param[in]        touches_regular_bdry Array storing which patches touch
    *                   non-periodic boundaries.
    * @param[in]        touches_periodic_bdry Array storing which patches touch
    *                   periodic boundaries.
    */
   /*
    * TODO:  See the second TODO item for the setGeometryOnPatches() method.
    */
   virtual void
   setGeometryDataOnPatch(
      Patch& patch,
      const IntVector& ratio_to_level_zero,
      const TwoDimBool& touches_regular_bdry,
      const TwoDimBool& touches_periodic_bdry) const;

   /*!
    * @brief Compute boundary boxes for each patch in patch level.
    *
    * Boundary boxes for each patch in a patch level will be computed and
    * assign to boundary box arrays.  These arrays are assumed to be
    * of length
    * @code
    *    DIM * num_patches
    * @endcode.
    *
    * The DIM arrays of boundary boxes for each patch will be stored in groups
    * of DIM. For example in 3D, with @c n patches on the level, the array
    * of boundary box arrays will be ordered as follows:
    *
    * @code
    * (patch 0 face array, patch 0 edge array, patch 0 mapped_box array,
    *  patch 1 face array, patch 1 edge array, patch 1 mapped_box array, . . . ,
    *  patch n-1 face array, patch n-1 edge array, patch n-1 mapped_box array)
    * @endcode
    *
    * @note
    * The optional argument @c do_all_patches defaults to false:
    * the boundary box computation is executed only on patches that touch a
    * non-periodic boundary.  When this routine is called during patch
    * level construction to describe a physical boundary, it is known that
    * only patches that touch a non-periodic boundary will have non-empty
    * sets of boundary boxes, so for efficiency's sake the boundary box
    * box computation is supressed for all other patches.
    *
    * When this routine is called to create boundary boxes that describe a
    * coarse-fine boundary, the computation must occur for every patch, so
    * @c do_all_patches must be set to true.
    *
    * @param[out]    boundaries output boundary description
    * @param[in]     level level on which to generate boundaries
    * @param[in]     periodic_shift periodic shift for the level
    *                (see getPeriodicShift())
    * @param[in]     ghost_width ghost width to compute geometry for the domain
    * @param[in]     domain Physical domain (in index space of level) for
    *                computing boundary boxes.
    * @param[in]     do_all_patches Flag to indicate boundary box computation
    *                on all patches, even those known to not touch a boundary
    */
   void
   computeBoundaryBoxesOnLevel(
      std::map<BoxId, PatchBoundaries>& boundaries,
      const PatchLevel& level,
      const IntVector& periodic_shift,
      const IntVector& ghost_width,
      const tbox::Array<BoxContainer>& domain,
      bool do_all_patches = false) const;

   /*!
    * @brief Compute boundary boxes for patch
    *
    * Decompose patch boundary region into pieces depending on spatial
    * dimensions. Boxes are extended along the boundary to the edge
    * of the ghost mapped_box_level if necessary.
    *
    * @param[out]    patch_boundaries output boundaries
    * @param[in]     box
    * @param[in]     domain_boxes
    * @param[in]     ghosts
    * @param[in]     periodic_shift
    */
   void
   getBoundaryBoxes(
      PatchBoundaries& patch_boundaries,
      const Box& box,
      const BoxContainer& domain_boxes,
      const IntVector& ghosts,
      const IntVector& periodic_shift) const;

   /*!
    * @brief Adjust boundary data of a level to be consistent with the
    * multiblock nature of the domain.
    *
    * In a multiblock problem, the PatchLevels contain patches that were
    * constructed independent of any knowledge of the multiblock nature of the
    * complete domain.  Thus the patches will contain boundary data that
    * recognizes no difference between a physical domain boundary and a block
    * boundary.  Calling this method will adjust the boundary data on all
    * patches in the given level such that the true boundaries of the domain
    * are represented.
    *
    * @param[in,out] patch_level Level where boundaries need to be adjusted.
    * TODO:  Incorporate into regular boundary box computation once
    * PatchLevel is multiblock-aware.
    */
   void
   adjustMultiblockPatchLevelBoundaries(
      PatchLevel& patch_level);

   /*!
    * @brief Add a concrete spatial coarsening operator.
    *
    * @param[in]  var_type_name The type name of the variable with which
    *             coarsen_op is associated.
    * @param[in]  coarsen_op The concrete coarsening operator to add to the
    *             lookup list.
    */
   void
   addCoarsenOperator(
      const char* var_type_name,
      const boost::shared_ptr<CoarsenOperator>& coarsen_op)
   {
      d_transfer_operator_registry->addCoarsenOperator(
         var_type_name,
         coarsen_op);
   }

   /*!
    * @brief Add a concrete spatial refinement operator.
    *
    * @param[in]  var_type_name The type name of the variable with which
    *             refine_op is associated.
    * @param[in]  refine_op The concrete refinement operator to add to the
    *             lookup list.
    */
   void
   addRefineOperator(
      const char* var_type_name,
      const boost::shared_ptr<RefineOperator>& refine_op)
   {
      d_transfer_operator_registry->addRefineOperator(
         var_type_name,
         refine_op);
   }

   /*!
    * @brief Add a concrete time interpolation operator.
    *
    * @param[in]  var_type_name The type name of the variable with which
    *             time_op is associated.
    * @param[in]  time_op The concrete time interpolation operator to add
    *             to the lookup list.
    */
   void
   addTimeInterpolateOperator(
      const char* var_type_name,
      const boost::shared_ptr<TimeInterpolateOperator>& time_op)
   {
      d_transfer_operator_registry->addTimeInterpolateOperator(
         var_type_name,
         time_op);
   }

   /*!
    * @brief Lookup function for coarsening operator.
    *
    * Search list for the spatial coarsening operator matching the
    * request for the given variable.  If the operator is found, a
    * pointer to it will be returned.  Otherwise, an unrecoverable
    * error will result and the program will abort.
    *
    * @param[in]     var The Variable for which the corresponding coarsening
    *                operator should match.
    * @param[in]     op_name The string identifier of the coarsening operator.
    */
   boost::shared_ptr<CoarsenOperator>
   lookupCoarsenOperator(
      const boost::shared_ptr<Variable>& var,
      const std::string& op_name)
   {
      return d_transfer_operator_registry->lookupCoarsenOperator(
         *this, var, op_name);
   }

   /*!
    * @brief Lookup function for refinement operator.
    *
    * Search list for the spatial refinement operator matching the
    * request for the given variable.  If the operator is found, a
    * pointer to it will be returned.  Otherwise, an unrecoverable
    * error will result and the program will abort.
    *
    * @param[in]     var The Variable for which the corresponding refinement
    *                operator should match.
    * @param[in]     op_name The string identifier of the refinement operator.
    */
   boost::shared_ptr<RefineOperator>
   lookupRefineOperator(
      const boost::shared_ptr<Variable>& var,
      const std::string& op_name)
   {
      return d_transfer_operator_registry->lookupRefineOperator(
         *this, var, op_name);
   }

   /*!
    * @brief Lookup function for time interpolation operator.
    *
    * Search list for the time interpolation operator matching the
    * request for the given variable.  If the operator is found, a
    * pointer to it will be returned.  Otherwise, an unrecoverable
    * error will result and the program will abort.
    *
    * @param[in]     var The Variable for which the corresponding time
    *                interpolation operator should match.
    * @param[in]     op_name The string identifier of the time interpolation
    *                operator.  \b Default: STD_LINEAR_TIME_INTERPOLATE
    */
   boost::shared_ptr<TimeInterpolateOperator>
   lookupTimeInterpolateOperator(
      const boost::shared_ptr<Variable>& var,
      const std::string& op_name =
         "STD_LINEAR_TIME_INTERPOLATE")
   {
      return d_transfer_operator_registry->lookupTimeInterpolateOperator(
         *this, var, op_name);
   }

   /*!
    * @brief Set a minimum value on the value returned by
    * getMaxTransferOpStencilWidth().
    *
    * This method allows users to specify a minimum value returned by
    * getMaxTransferOpStencilWidth().  The default minimum is zero.
    * This value can be used as a substitute for data that is not yet
    * registered with the Geometry and therefore cannot be reflected
    * in getMaxTransferOpStencilWidth().
    *
    * @param[in]  min_value The minimum value to set.
    */
   void
   setMinTransferOpStencilWidth(
      const IntVector& min_value)
   {
      d_transfer_operator_registry->setMinTransferOpStencilWidth(min_value);
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
    * @brief Get the name of this object.
    *
    * @return The name of this object.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

   /*!
    * @brief Class to represent the neighbor of a given block.
    */

   class Neighbor
   {
public:
      /*!
       * @brief Constructor
       *
       * @param[in] block_id   The block id of the neighboring block
       * @param[in] domain     The neighboring block's domain in the current
       *                       block's index space
       * @param[in] transformation The transformation needed to align the
       *                           neighboring index spaces
       * @param[in] is_singularity True if the current block and the
       *                           neighboring block abut at a reduced
       *                           or enhanced connectivity singularity
       */
      Neighbor(
         const BlockId& block_id,
         const BoxContainer& domain,
         const Transformation& transformation,
         const bool is_singularity);

      /*!
       * @brief Get the block number of the neighboring block.
       */
      const BlockId&
      getBlockId() const
      {
         return d_block_id;
      }

      /*!
       * @brief Get the neighboring block's domain in the current block's
       * index space.
       */
      const BoxContainer&
      getTransformedDomain() const
      {
         return d_transformed_domain;
      }

      /*!
       * @brief Get the Transformation for the neighbor relationship.
       */
      const Transformation&
      getTransformation() const
      {
         return d_transformation;
      }

      /*!
       * @brief Get the rotation identifier for the neighbor relationship.
       */
      Transformation::RotationIdentifier
      getRotationIdentifier() const
      {
         return d_transformation.getRotation();
      }

      /*!
       * @brief Get the shift for the neighbor relationship.
       */
      const IntVector&
      getShift() const
      {
         return d_transformation.getOffset();
      }

      /*!
       * @brief Tell if the neighboring block touch each other at an
       * enhanced connectivity singularity.
       */
      bool
      isSingularity() const
      {
         return d_is_singularity;
      }

private:
      /*!
       * @brief The block number of the neighboring block
       */
      BlockId d_block_id;

      /*!
       * @brief The neighboring block's domain in the current block's
       * index space.
       */
      BoxContainer d_transformed_domain;

      /*!
       * @brief The transformation to transform the neighboring block's
       * coordinate system into the current coordinate system.
       */
      Transformation d_transformation;

      /*!
       * True if the current block and the neighboring block abut at a
       * reduced or enhanced connectivity singularity
       */
      bool d_is_singularity;

   };

   /*!
    * @brief Register a relationship between two neighboring blocks of a
    * multiblock domain.
    *
    * @param[in] block_a         One block in the relationship
    * @param[in] block_b         The other block
    * @param[in] rotation_b_to_a The rotation that aligns block b's index space
    *                            with block a's
    * @param[in] shift_b_to_a    The post-rotation shift to move b into its
    *                            correct location within a's index space
    * @param[in] neighbor_type   The type (codimension) of the neighbor
    *                            relationship
    */
   void
   registerNeighbors(
      const BlockId& block_a,
      const BlockId& block_b,
      const Transformation::RotationIdentifier rotation_b_to_a,
      const IntVector& shift_b_to_a,
      const int neighbor_type);

   /*!
    * @brief Get a BoxContainer that contains all of the index space of all other
    * blocks in the multiblock domain.
    *
    * A BoxContainer will be constructed that contains the full set of the
    * coarse level domains of all blocks except the one identified by
    * block_id.  The domains will all be transformed into the index space
    * represented by block_id.
    *
    * @param[out] domain_outside_block
    * @param[in] block_id
    *
    */
   void
   getDomainOutsideBlock(
      BoxContainer& domain_outside_block,
      const BlockId& block_id) const;

   /*!
    * @brief Return the number of block singularities in the block
    * configuration.
    */
   int
   getNumberOfBlockSingularities() const
   {
      return d_number_of_block_singularities;
   }

   /*!
    * @brief Return a BoxContainer that describes all of the singularities
    * touched by the block indicated by block_id.
    *
    * @return For every singularity point the block touches, the BoxContainer will
    * contain a single-cell box that lies just outside the block domain,
    * touching the block only at the singularity point.  For line singularities,
    * the BoxContainer will contain boxes of width 1 in all dimensions except one,
    * lying outside the block's coarse-level domain and touching the domain
    * only along the line of singularity.
    *
    * @param[in] block_id
    */
   const BoxContainer&
   getSingularityBoxContainer(
      const BlockId& block_id) const
   {
      return d_singularity[block_id.getBlockValue()];
   }

   /*!
    * @brief Return a list of integers indicating all of the
    * singularities touched by the block indicated by block_id.
    *
    * @return For every singularity point the block touches, the
    * vector<int> will contain the index of that singularity.
    *
    * @param[in] block_id
    */
   const std::vector<int>&
   getSingularityIndices(
      const BlockId& block_id) const
   {
      return d_singularity_indices[block_id.getBlockValue()];
   }

   /*!
    * @brief Tell if block represented by block_id touches
    * a reduced-connectivity singularity
    *
    * @return True if the block touches reduced connectivity singularity,
    * false if not.
    *
    * @param[in] block_id
    */
   bool
   reducedConnectivityExists(
      const BlockId& block_id) const
   {
      return d_reduced_connect[block_id.getBlockValue()];
   }

   /*!
    * @brief Modify a box by rotating and shifting from the index space of
    * the transformed_block to the index space of the base_block at the
    * resolution level defined by ratio_to_level_zero.
    *
    * @param[in,out] box The boxes will be transformed from the
    *                      transformed_block index space to the base_block
    *                      index space.
    * @param[in] ratio
    * @param[in] output_block Integer identifier of the block whose index space
    *                       will be represented in the boxes at output
    * @param[in] input_block Integer identifier of the block whose index
    *                             space is represented in the boxes at input
    *
    * @return Whether the box has been transformed.  True if there is a
    * relationship between input_block and output_block.  False if
    * there is no such relationship.
    */
   bool
   transformBox(
      Box& box,
      const IntVector& ratio,
      const BlockId& output_block,
      const BlockId& input_block) const;

   /*!
    * @brief Modify boxes by rotating and shifting from the index space of
    * the input_block to the index space of the output_block at the
    * resolution level defined by ratio_to_level_zero.
    *
    * @param[in,out] boxes The boxes will be transformed from the
    *                      input_block index space to the output_block
    *                      index space.
    * @param[in] ratio
    * @param[in] output_block Integer identifier of the block whose index space
    *                       will be represented in the boxes at output
    * @param[in] input_block Integer identifier of the block whose index
    *                             space is represented in the boxes at input
    *
    * @return Whether the boxes have been transformed.  True if there
    * is a relationship between input_block and output_block.  False
    * if there is no such relationship.
    */
   bool
   transformBoxContainer(
      BoxContainer& boxes,
      const IntVector& ratio,
      const BlockId& output_block,
      const BlockId& input_block) const;

   /*!
    * @brief Get a box array that describes the coarse-level domain of the
    * transformed_block in terms of the index space of base_block.
    *
    * @param[out] block_boxes The coarse-level domain of the block
    *                         identified by transformed_block, represented
    *                         in the index space of the block identified by
    *                         base_block
    * @param[in] base_block  The block whose index space will be used for
    *                        the output boxes
    * @param[in] transformed_block ID of another block whose
    *                             domain will be represented in the index space
    *                             of the base block
    */
   void
   getTransformedBlock(
      BoxContainer& block_boxes,
      const BlockId& base_block,
      const BlockId& transformed_block);

   /*!
    * @brief Return a list of Neighbor objects describing all of the neighbors
    * of the block indicated by the block_id.
    *
    * @return The list of neighbors
    *
    * @param[in] block_id
    */
   const std::list<Neighbor>&
   getNeighbors(
      const BlockId& block_id) const
   {
      return d_block_neighbors[block_id.getBlockValue()];
   }

   /*!
    * @brief Return the number of neighbors a specific block of the Multiblock
    * domain has.
    *
    * A block is the neighbor of another block if the two blocks abut in
    * any way, whether at a point, a 1D line, or a 2D plane.
    *
    * @return The number of neighbors
    *
    * @param[in] block_id
    */
   int
   getNumberOfNeighbors(
      const BlockId& block_id) const
   {
      return static_cast<int>(
         d_block_neighbors[block_id.getBlockValue()].size());
   }

   /*!
    * @brief Tell if the given BlockIds represent neighboring blocks.
    */
   bool
   areNeighbors(
      const BlockId& block_a,
      const BlockId& block_b) const;

   /*!
    * @brief Tell if the given BlockIds represent neighboring blocks.
    */
   bool
   areSingularityNeighbors(
      const BlockId& block_a,
      const BlockId& block_b) const;

   /*!
    * @brief Get the rotation identifier to rotate from src to dst.
    */
   Transformation::RotationIdentifier
   getRotationIdentifier(
      const BlockId& dst,
      const BlockId& src) const;

   /*!
    * @brief Get the offset to shift from src to dst after rotation.
    */
   const IntVector&
   getOffset(
      const BlockId& dst,
      const BlockId& src) const;

   /*!
    * @brief Query if the geometry has enhanced connectivity.
    */
   bool
   hasEnhancedConnectivity() const
   {
      return d_has_enhanced_connectivity;
   }

   /*!
    * @brief Print object data to the specified output stream.
    *
    * @param[out] stream The output stream (as a std::ostream&) to print to.
    */
   void
   printClassData(
      std::ostream& stream) const;

   /*!
    * @brief Writes the state of the BaseGridGeometry object to the
    * database.
    *
    * When assertion checking is active, db cannot be a null database pointer.
    *
    * @param[in,out]    db The database to write to write/serialize.
    */
   virtual void
   putToDatabase(
      const boost::shared_ptr<tbox::Database>& db) const;

protected:
   /*!
    * @brief Construct a new BaseGridGeometry object in its default state.
    *
    * This constructor is intended to be called from a child class derived
    * from BaseGridGeometry.  It will not register for restart nor read any
    * input data, as it is expected that the child class will handle those
    * operations.
    *
    * @param[in]  dim
    * @param[in]  object_name
    * @param[in]  op_reg
    */
   BaseGridGeometry(
      const tbox::Dimension& dim,
      const std::string& object_name,
      const boost::shared_ptr<TransferOperatorRegistry>& op_reg);

   /*!
    * @brief Construct a new BaseGridGeometry object in its default state.
    *
    * This constructor is intended to be called from a child class derived
    * from BaseGridGeometry.  It will not register for restart nor read any
    * input data, as it is expected that the child class will handle those
    * operations.
    *
    * @param[in]  dim
    * @param[in]  object_name
    */
   BaseGridGeometry(
      const tbox::Dimension& dim,
      const std::string& object_name);

   /*!
    * @brief Construct a new BaseGridGeometry object based on arguments.
    *
    * @param[in]  object_name
    * @param[in]  domain      Each element of the array describes the index
    *                         space for a block.
    * @param[in]  op_reg
    * @param[in]  register_for_restart Flag indicating whether this instance
    *             should be registered for restart.  @b Default: true
    */
   BaseGridGeometry(
      const std::string& object_name,
      const BoxContainer& domain,
      const boost::shared_ptr<TransferOperatorRegistry>& op_reg,
      bool register_for_restart);

   /*!
    * @brief Read multiblock metadata input from the input database
    *
    * @param[in] input_db
    */
   void
   readBlockDataFromInput(
      const boost::shared_ptr<tbox::Database>& input_db);

   /*!
    * The holder of all the transfer operators.
    */
   boost::shared_ptr<TransferOperatorRegistry> d_transfer_operator_registry;

   tbox::Dimension d_dim;

private:
   /*!
    * @brief Virtual method to build operators appropriate to a specific
    * grid geometry.  This must be defined by all sub-classes of
    * BaseGridGeometry.
    */
   virtual void
   buildOperators() = 0;

   /*!
    * @brief Reset domain BoxContainer after data it depends on has changed.
    */
   void
   resetDomainBoxContainer();

   /*!
    * @brief Check that the domain is valid for periodic boundary conditions
    */
   bool
   checkPeriodicValidity(
      const BoxContainer& domain);

   /*!
    * @brief Check on each BoundaryBox when it is created.
    *
    * This is a check performed on each BoundaryBox when it is created.
    * It returns true when a BoundaryBox has a width of 1 in at least
    * one direction, is adjacent to the patch boundary (possible extended
    * into the patch's ghost region) and is outside the physical domain.
    */
   bool
   checkBoundaryBox(
      const BoundaryBox& boundary_box,
      const Patch& patch,
      const BoxContainer& domain,
      const int num_per_dirs,
      const IntVector& max_data_ghost_width) const;

   /*!
    * @brief Compute shifts for a box on periodic boundary.
    *
    * If box is located on a periodic boundary, all of its possible shifts
    * will be computed and stored in shifts.  If box is not on a periodic
    * boundary, shifts will be an empty list.
    *
    * @param[out]    shifts vector storing the periodic shifts for a box
    * @param[in]     box Box on which to compute periodic shifts
    * @param[in]     domain_search_tree search tree for the domain
    * @param[in]     periodic_shift
    */
   void
   computeShiftsForBox(
      std::vector<IntVector>& shifts,
      const Box& box,
      const BoxContainer& domain_search_tree,
      const IntVector& periodic_shift) const;

   /*!
    * @brief Adjust the BoundaryBoxes held by the Patch so that they take
    * into account the full multiblock configuration.
    *
    * @param[in] patch
    * @param[in] geometry
    * @param[in] pseudo_domain The full multiblock domain represented as if
    *            every block were in the index space where the patch exists
    * @param[in] gcw         The maximum patch data ghost width
    * @param[in] singularity BoxContainer obtained by getSingularityBoxContainer()
    */

   void
   adjustBoundaryBoxesOnPatch(
      const Patch& patch,
      const BoxContainer& pseudo_domain,
      const IntVector& gcw,
      const BoxContainer& singularity);

   /*!
    * @brief Reads in data from the specified input database.
    *
    * If the simulation is from restart, these values are taken from restart
    * and newly specified values in the input database are ignored.
    *
    * @param[in] db  input database, must not be NULL pointer
    * @param[in] is_from_restart  set to true if simulation is from restart
    */
   void
   getFromInput(
      const boost::shared_ptr<tbox::Database>& db,
      bool is_from_restart);

   /*!
    * @brief Read object state from the restart file and initialize class data
    * members.
    *
    * The database from which the restart data is read is
    * determined by the object_name specified in the constructor.
    *
    * Unrecoverable Errors:
    *
    *    -The database corresponding to object_name is not found
    *     in the restart file.
    *
    *    -The class version number and restart version number do not
    *     match.
    */
   void
   getFromRestart();

   /*!
    * @brief Initialize static objects and register shutdown routine.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   initializeCallback();

   /*!
    * @brief Method registered with ShutdownRegister to cleanup statics.
    *
    * Only called by StartupShutdownManager.
    */
   static void
   finalizeCallback();

   /*!
    * Object name used for error reporting purposes.
    */
   std::string d_object_name;

   /*!
    * BoxContainer defining computational domain on coarsest level.
    */
   BoxContainer d_physical_domain;

   /*!
    * Boolean array telling for each block whether the domain of that block
    * can be represented by a single box.
    */
   tbox::Array<bool> d_domain_is_single_box;

   /*!
    * @brief BoxContainer representation of the physical domain, including
    * its periodic image Boxes.
    */
   BoxContainer d_domain_with_images;

   /*!
    * Integer array vector describing periodic shift coarsest level.
    * An entry of zero means direction is not periodic.
    */
   IntVector d_periodic_shift;

   /*!
    * Current maximum ghost cell width over all patch data objects
    * known to the patch descriptor.  This is used to compute
    * boundary boxes.
    */
   IntVector d_max_data_ghost_width;

   /*!
    * The number of blocks represented in the BaseGridGeometry.
    */
   int d_number_blocks;

   /*!
    * The number of blocks singularities in the block configuration.
    */
   int d_number_of_block_singularities;

   /*!
    * @brief Associated with each block is a list of Neighbors that
    * it shares a block boundary with.
    */
   tbox::Array<std::list<Neighbor> > d_block_neighbors;

   /*!
    * @brief An array of BoxContainers defining the singularities of a multiblock
    * domain.  Each BoxContainer element defines the singularities that a single
    * block touches.
    */
   tbox::Array<BoxContainer> d_singularity;

   /*!
    * @brief An array of singularity indices of a multiblock
    * domain.  d_singularity_indices[bn] is a list of singularity indices
    * touched by block bn.
    */
   tbox::Array<std::vector<int> > d_singularity_indices;

   /*!
    * @brief Tell whether each block touches a reduced-connectivity
    * singularity.
    */
   tbox::Array<bool> d_reduced_connect;

   /*!
    * @brief Tell whether there is enhanced connectivity anywhere in the
    * geometry.
    */
   bool d_has_enhanced_connectivity;

   /*!
    * Flag to determine whether this instance is registered for restart.
    */
   bool d_registered_for_restart;

   static boost::shared_ptr<tbox::Timer> t_find_patches_touching_boundaries;
   static boost::shared_ptr<tbox::Timer> t_touching_boundaries_init;
   static boost::shared_ptr<tbox::Timer> t_touching_boundaries_loop;
   static boost::shared_ptr<tbox::Timer> t_set_geometry_on_patches;
   static boost::shared_ptr<tbox::Timer> t_set_boundary_boxes;
   static boost::shared_ptr<tbox::Timer> t_set_geometry_data_on_patches;
   static boost::shared_ptr<tbox::Timer> t_compute_boundary_boxes_on_level;
   static boost::shared_ptr<tbox::Timer> t_get_boundary_boxes;

   /*
    * Static initialization and cleanup handler.
    */
   static tbox::StartupShutdownManager::Handler
      s_initialize_handler;

};

}
}

#endif
