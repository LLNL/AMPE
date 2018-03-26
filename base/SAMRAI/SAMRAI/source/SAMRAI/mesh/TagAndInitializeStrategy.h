/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Strategy interface for params, tagging, init for gridding.
 *
 ************************************************************************/

#ifndef included_mesh_TagAndInitializeStrategy
#define included_mesh_TagAndInitializeStrategy

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/BoxLevel.h"
#include "SAMRAI/tbox/Array.h"

#include <boost/shared_ptr.hpp>

namespace SAMRAI {
namespace mesh {

/*!
 * Class TagAndInitializeStrategy is a base class that defines a
 * Strategy pattern interface for level initialization and cell tagging
 * routines that are needed by the adaptive meshing algorithms provided
 * by the class GriddingAlgorithm.  The class
 * maintains functionality to construct refined regions based on a
 * user-supplied set of boxes, but its main role is to provide interfaces
 * for level initialization and cell tagging operations.
 *
 * The operations that identify mesh cells for refinement or initialize
 * data on a new hierarchy level are problem-specific and must be
 * supplied by a concrete sub-class of this base class.
 *
 * If user supplied refine boxes are used, they may be supplied through
 * input.  Alternatively, they may be supplied through the "resetRefineBoxes()"
 * method.  If they are supplied through input, the format is as follows:
 *
 *    - \b    RefineBoxes
 *      input section describing the refine boxes for each level.
 *          - \b  level_0
 *             input section providing the hier::Box arrays
 *             describing where user-specified refinement is to occur on
 *             Level 0.
 *          - \b  level_1
 *             input section providing the hier::Box arrays
 *             describing where user-specified refinement is to occur on
 *             Level 1.
 *            \b . . .
 *          - \b  level_n
 *             input section providing the hier::Box arrays
 *             describing where user-specified refinement is to occur on
 *             Level N.
 *
 *       For each level, the input section can have the following entries:
 *
 *               - \b times = optional entry, a double array specifying times
 *                            at which a particular set of boxes is to be used
 *                            as a region of refinement.
 *               - \b cycles = optional entry, an integer array specifying
 *                             regrid cycles at which a particular set of boxes
 *                             is to be used as a region of refinement.
 *               - \b boxes_0 = box array specifying refine boxes for element
 *                              0 of the times or cycles array.
 *               - \b boxes_1 = box array specifying refine boxes for element
 *                              1 of the times or cycles array.
 *               - \b boxes_n = box array specifying refine boxes for element
 *                              n of the times or cycles array.
 *
 *       The @b times and @b cycles entries are optional.  If neither
 *       is provided, a uniform set of refine boxes specified in the
 *       boxes_0 entry will be used over the entire calculation.  If no
 *       boxes_0 entry is provided, no refinement will occur on that level.
 *
 *       If both @b times or @b cycles entries are supplied, the times entry
 *       takes precedence so the cycles entry is ignored. The particular
 *       box array chosen during regridding is determined by a ``greater-than''
 *       convention.  That is, if boxes are accessed at regridding time t,
 *       where t is greater-than the specified times[n] entry, then the array
 *       given for boxes_n is used.  Otherwise, the corresponding previous
 *       box array that satisfies the criteria is used.  The same convention is
 *       followed for regridding cycles.  To avoid errant behavior, the times
 *       and cycles entries should always be supplied in increasing order.
 *
 *       The hier::BoxContainer entries withing each level's input section
 *       must be of the form ``boxes_n'' (where n is the corresponds to the
 *       elements of the times or cycles array), or the input parser will
 *       ignore the entry.  If there is no ``boxes_n'' entry corresponding
 *       to element n of the times or cycles array, then no refinement will
 *       occur on that level at the given time or cycle.
 *
 * A sample input file entry might look like:
 *
 * \verbatim
 *
 *    RefineBoxes {
 *       level_0 {
 *          cycles = 0, 10
 *          boxes_0 = [(5,5),(9,9)],[(12,15),(18,19)]
 *          boxes_1 = [(7,7),(11,11)],[(14,17),(20,21)]
 *       }
 *       level_1 {
 *          times  = 0., 0.05, 0.10
 *          boxes_0 = [(25,30),(29,35)]
 *          boxes_1 = [(30,35),(34,40)]
 *          boxes_2 = [(35,40),(39,45)]
 *       }
 *       level_2 {
 *          boxes_0 = [(60,70),(70,80)]
 *       }
 *    }
 * \endverbatim
 *
 * The virtual methods in this class may place constraints on the patch
 * hierarchy by the particluar error estimation procedure in use.  Those
 * constraints and operations must be honored in the concrete subclass
 * implementations of these methods.  The constraints are discussed in
 * the method descriptions below.
 *
 * @see mesh::GriddingAlgorithm
 */

class TagAndInitializeStrategy
{
public:
   /*!
    * Empty constructor for TagAndInitializeStrategy.
    */
   TagAndInitializeStrategy(
      const tbox::Dimension& dim,
      const std::string& object_name);

   /*!
    * Empty destructor for TagAndInitializeStrategy.
    */
   virtual ~TagAndInitializeStrategy();

   /*!
    * Return user supplied set of refine boxes for specified level number
    * and time.  The boolean return value specifies whether the boxes
    * have been reset from the last time this method was called.  If they
    * have been reset, it returns true.  If they are unchanged, it returns
    * false.
    */
   bool
   getUserSuppliedRefineBoxes(
      hier::BoxContainer& refine_boxes,
      const int level_number,
      const double time);

   /*!
    * Reset the static refine boxes for the specified level number in the
    * hierarchy.  The level number must be greater than or equal to zero.
    */
   void
   resetRefineBoxes(
      const hier::BoxContainer& refine_boxes,
      const int level_number);

   /*!
    * Initialize data on a new level after it is inserted into an AMR patch
    * hierarchy by the gridding algorithm.  The level number indicates
    * that of the new level.  The old_level pointer corresponds to
    * the level that resided in the hierarchy before the level with the
    * specified number was introduced.  If this pointer is null, there was
    * no level in the hierarchy prior to the call and the data on the new
    * level is set by interpolating data from coarser levels in the hierarchy.
    * Otherwise, the the new level is initialized by interpolating data from
    * coarser levels and copying data from the old level before it is
    * destroyed.
    *
    * The boolean argument initial_time indicates whether the integration
    * time corresponds to the initial simulation time.  If true, the level
    * should be initialized with initial simulation values.  Otherwise, it
    * should be assumed that the simulation time is at some point after the
    * start of the simulation.  This information is provided since the
    * initialization of the data may be different in each of those
    * circumstances.  In any case, the double "time" value is the current
    * simulation time for the level.  The can_be_refined boolean argument
    * indicates whether the level is the finest allowable level in the
    * hierarchy.  This flag is included since data management on the finest
    * level may be different than other levels in the hierarchy in some cases.
    *
    * The last two (optional) arguments specify an old level from which the
    * data may be used to initialize data on this level, and a flag that
    * indicates whether data on the initialized level must first be allocated.
    * The allocate_data argument is used in cases where one wishes to
    * simply reset data to an initialized state on a level that has already
    * been allocated.
    */
   virtual void
   initializeLevelData(
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double init_data_time,
      const bool can_be_refined,
      const bool initial_time,
      const boost::shared_ptr<hier::PatchLevel>& old_level =
         boost::shared_ptr<hier::PatchLevel>(),
      const bool allocate_data = true) = 0;

   /*!
    * After hierarchy levels have changed and data has been initialized on
    * the new levels, this routine can be used to reset any information
    * needed by the solution method that is particular to the hierarchy
    * configuration.  For example, the solution procedure may cache
    * communication schedules to amortize the cost of data movement on the
    * AMR patch hierarchy.  This function will be called by the gridding
    * algorithm after the initialization occurs so that the algorithm-specific
    * subclass can reset such things.  Also, if the solution method must
    * make the solution consistent across multiple levels after the hierarchy
    * is changed, this process may be invoked by this routine.  Of course the
    * details of these processes are determined by the particular solution
    * methods in use.
    *
    * The level number arguments indicate the coarsest and finest levels
    * in the current hierarchy configuration that have changed.  It should
    * be assumed that all intermediate levels have changed as well.
    */
   virtual void
   resetHierarchyConfiguration(
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level) = 0;

   /*!
    * Set integer tags to "one" on the given level to identify
    * where refinement of that level should occur.  The index is that of the
    * cell-centered integer tag array on each patch.  The boolean argument
    * initial_time indicates whether cells are being tagged at
    * initialization time, or at some later time during the calculation.
    * If it is false, it should be assumed that the error estimation process
    * is being invoked at some later time after the AMR hierarchy was
    * initially constructed.  This information is provided since application
    * of the error estimator may be different in each of those circumstances.
    *
    * The cell-tagging operation may use time advancement to determine
    * tag regions. The argument coarsest_sync_level provides information
    * for the tagging method to coordinate time advance with an integrator.
    * When time integration is used during regridding, this value is true
    * if the level is the coarsest level involved in level synchronization
    * immediately preceeding the regrid process; otherwise it is false.
    * If time advancement is not used, this argument are ignored.
    *
    * The boolean can_be_refined is used to coordinate data reset operations
    * with the time integrator when time-dependent regridding is used.  This
    * is provided since data may be managed differently on the finest hierarchy
    * level than on coarser levels.
    */
   virtual void
   tagCellsForRefinement(
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double error_data_time,
      const int tag_index,
      const bool initial_time,
      const bool coarsest_sync_level,
      const bool can_be_refined = true,
      const double regrid_start_time = 0.) = 0;

   /*!
    * Certain cases may require pre-processing of error estimation data
    * before tagging cells, which is handled by this method. For example,
    * Richardson extrapolation may require advances of data in time before
    * the error estimation procedure is implemented.
    *
    * The level number indicates the level in which pre-process steps
    * are applied, time is the time at which the operation is performed
    * (generally the regrid time), and the boolean argument indicates
    * whether the operation is performed at the initial time.
    */
   virtual void
   preprocessErrorEstimation(
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double regrid_time,
      const double regrid_start_time,
      const bool initial_time) = 0;

   /*!
    * Return true if regridding process advances the data using some time
    * integration procedure; otherwise, return false.
    */
   virtual bool
   usesTimeIntegration() const = 0;

   /*!
    * Return true if boxes for coarsest hierarchy level are not appropriate
    * for gridding strategy.  Otherwise, return false.  If false is returned,
    * it is useful to provide a detailed explanatory message describing the
    * problems with the boxes.
    */
   virtual bool
   coarsestLevelBoxesOK(
      const hier::BoxContainer& boxes) const = 0;

   /*!
    * Return ratio by which level may be coarsened during the error
    * estimation process.  Generally, this is needed by the gridding
    * algorithm class so that the new patch levels that it constructs can
    * be coarsened properly (if needed) during the error estimation process.
    */
   virtual int
   getErrorCoarsenRatio() const = 0;

   /*!
    * Check ratios between hierarchy levels against any constraints that
    * may be required for the error estimation scheme.
    */
   virtual void
   checkCoarsenRatios(
      const tbox::Array<hier::IntVector>& ratio_to_coarser) = 0;

   /*!
    * Return whether refinement is being performed using ONLY
    * user-supplied refine boxes.  If any method is used that invokes
    * tagging, this will return false.
    */
   virtual bool
   refineUserBoxInputOnly() const = 0;

   /*!
    * Read user supplied refine boxes from the provided database.  The
    * database must be non-null, or an unrecoverable assertion will be thrown.
    */
   void
   getFromInput(
      const boost::shared_ptr<tbox::Database>& db);

   /*!
    * Return the dimension of this object.
    */
   const tbox::Dimension&
   getDim() const
   {
      return d_dim;
   }

   /*!
    * Returns the object name.
    */
   const std::string&
   getObjectName() const
   {
      return d_object_name;
   }

private:
   const tbox::Dimension d_dim;

   std::string d_object_name;

   /*
    * Arrays of data for user-specified refinement.  The user controls
    * the particular boxes to be used for refinement by specifying
    * "cycles" -or- "times" and "refine_boxes".  The arrays below hold
    * hold entries for each level and seq number.  The boolean array
    * specifies whether to use time or cycles as the guiding criteria
    * (by default, time is used).  The integer cycle counter holds
    * internally the number of times the getRefineBoxes() method has
    * been accessed for each level.
    */
   tbox::Array<tbox::Array<hier::BoxContainer> > d_refine_boxes;
   tbox::Array<tbox::Array<int> > d_refine_boxes_cycles;
   tbox::Array<tbox::Array<double> > d_refine_boxes_times;
   tbox::Array<bool> d_refine_boxes_use_times;
   tbox::Array<int> d_refine_boxes_cycle_counter;

   /*
    * Arrays to hold boxes that are specifically reset by the user (via the
    * resetRefineBoxes() method).  The boolean array specifies which levels
    * have been reset while the box array specifies the new set of refine
    * boxes for the level. The int array holds the sequence number from
    * the last time getRefineBoxes() was called, allowing us to
    * determine when refine boxes change between steps.
    */
   tbox::Array<bool> d_refine_boxes_reset;
   tbox::Array<hier::BoxContainer> d_reset_refine_boxes;
   tbox::Array<int> d_refine_boxes_old_seq_num;

   // The following are not implemented:
   TagAndInitializeStrategy(
      const TagAndInitializeStrategy&);
   void
   operator = (
      const TagAndInitializeStrategy&);

};

}
}

#endif
