/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Gridding routines and params for Richardson Extrapolation.
 *
 ************************************************************************/

#ifndef included_mesh_StandardTagAndInitialize
#define included_mesh_StandardTagAndInitialize

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/mesh/StandardTagAndInitializeConnectorWidthRequestor.h"
#include "SAMRAI/mesh/TagAndInitializeStrategy.h"

#include <boost/shared_ptr.hpp>

namespace SAMRAI {
namespace mesh {

/*!
 * Class StandardTagAndInitialize defines an implementation
 * for level initialization and cell tagging routines needed by
 * the GriddingAlgorithm class.  This class is derived from
 * the abstract base class TagAndInitializeStrategy. It invokes
 * problem-specific level initialization routines after AMR patch
 * hierarchy levels change and routines for tagging cells for refinement
 * using one (or more) of the following methods:
 *
 *   - Gradient Detection
 *   - Richardson Extrapolation
 *   - Explicitly defined refine boxes
 *
 * It is possible to use combinations of these three methods (e.g.,
 * use gradient detection, Richardson extrapolation, and static refine boxes
 * at the same time). The order in which they are executed is fixed (
 * Richardson extrapolation first, gradient detection second, and refine
 * boxes third).  An input entry for this class is optional.
 * If none is provided, the class will by default not use any criteria
 * to tag cells for refinement.
 *
 * Required input keys and data types: NONE
 *
 * Optional input keys, data types, and defaults:
 *
 *
 *    - \b    tagging_method
 *       std::string array specification of the type of cell-tagging used.  Valid
 *       choices include:
 *          - ``GRADIENT_DETECTOR''
 *          - ``RICHARDSON_EXTRAPOLATION''
 *          - ``REFINE_BOXES''
 *          - ``RICHARDSON_EXTRAPOLATION'', ``GRADIENT_DETECTOR'',
 *                ``REFINE_BOXES''
 *       (i.e. a combination of any or all of the above - the choices may
 *       be placed in any order). If no input is given, no tagging will be
 *       performed.
 *
 *    - \b    input section describing the refine boxes for each level.
 *      (@see mesh::TagAndInitializeStrategy for details on format)
 *
 * A sample input file entry might look like:
 *
 * \verbatim
 *
 *    tagging_method = "GRADIENT_DETECTOR", "REFINE_BOXES"
 *    <refine boxes input> (@see mesh::TagAndInitializeStrategy)
 *
 * \endverbatim
 *
 * This class supplies the routines for tagging cells
 * and invoking problem-specific level initialization routines after AMR
 * patch hierarchy levels change.  A number of problem-specific operations
 * are performed in the StandardTagAndInitStrategy
 * data member, for which methods are specified in a derived subclass.
 *
 * @see mesh::TagAndInitializeStrategy
 * @see mesh::GriddingAlgorithm
 * @see mesh::StandardTagAndInitStrategy
 */

class StandardTagAndInitialize:
   public TagAndInitializeStrategy
{
public:
   /*!
    * Constructor for StandardTagAndInitialize which
    * may read inputs from the provided input_db.  If no input
    * database is provided, the class interprets that no tagging
    * is desired so no cell-tagging will be performed.
    */
   StandardTagAndInitialize(
      const tbox::Dimension& dim,
      const std::string& object_name,
      StandardTagAndInitStrategy* tag_strategy,
      const boost::shared_ptr<tbox::Database>& input_db =
         boost::shared_ptr<tbox::Database>());

   /*!
    * Virtual destructor for StandardTagAndInitialize.
    */
   virtual ~StandardTagAndInitialize();

   /*!
    * Specifies whether the chosen method advances the solution data
    * in the regridding process (Richardson extrapolation does, the
    * others will not).
    */
   bool
   usesTimeIntegration() const;

   /*!
    * Return coarsen ratio used for applying cell tagging. An error
    * coarsen ratio other than 2 or 3 will throw an error.
    */
   int
   getErrorCoarsenRatio() const;

   /*!
    * Some restrictions may be placed on the coarsen ratio used for
    * cell tagging.  Check these here.
    */
   void
   checkCoarsenRatios(
      const tbox::Array<hier::IntVector>& ratio_to_coarser);

   /*!
    * Pass the request to initialize the data on a new level in the
    * hierarchy to the StandardTagAndInitStrategy data member. Required
    * arguments specify the grid hierarchy, level number being initialized,
    * simulation time at which the data is initialized, whether the level
    * can be refined, and whether it is the initial time.  Optional arguments
    * include an old level, from which data may be used to initialize this
    * level, and a flag that indicates whether data on the initialized level
    * must first be allocated.  For more information on the operations that
    * must be performed, see the
    * TagAndInitializeStrategy::initializeLevelData() method.
    */
   void
   initializeLevelData(
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double init_data_time,
      const bool can_be_refined,
      const bool initial_time,
      const boost::shared_ptr<hier::PatchLevel>& old_level =
         boost::shared_ptr<hier::PatchLevel>(),
      const bool allocate_data = true);

   /*!
    * Pass the request to reset information that depends on the hierarchy
    * configuration to the StandardTagAndInitStrategy data member.
    * For more information on the operations that must be performed, see
    * the TagAndInitializeStrategy::resetHierarchyConfiguration()
    * method.
    */
   void
   resetHierarchyConfiguration(
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int coarsest_level,
      const int finest_level);

   /*!
    * Certain cases may require pre-processing of error estimation data
    * before tagging cells, which is handled by this method.  For more
    * information on the operations that must be performed, see the
    * TagAndInitializeStrategy::preprocessErrorEstimation()
    * method
    */
   void
   preprocessErrorEstimation(
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double regrid_time,
      const double regrid_start_time,
      const bool initial_time);

   /*!
    * Pass the request to set tags on the given level where refinement of
    * that level should occur.  Gradient detection, Richardson extrapolation,
    * and tagging on static refine boxes is performed here.
    *
    * For more information on the operations that must be performed, see the
    * TagAndInitializeStrategy::tagCellsForRefinement() routine.
    */
   void
   tagCellsForRefinement(
      const boost::shared_ptr<hier::PatchHierarchy>& level,
      const int level_number,
      const double regrid_time,
      const int tag_index,
      const bool initial_time,
      const bool coarsest_sync_level,
      const bool can_be_refined,
      const double regrid_start_time = 0);

   /*!
    * Return true if boxes for coarsest hierarchy level are not appropriate
    * for gridding strategy.  Otherwise, return false.  If false is returned,
    * it is useful to provide a detailed explanatory message describing the
    * problems with the boxes.
    */
   bool
   coarsestLevelBoxesOK(
      const hier::BoxContainer& boxes) const;

   /*!
    * Return whether refinement is being performed using ONLY
    * user-supplied refine boxes.  If any method is used that invokes
    * tagging, this will return false.
    */
   bool
   refineUserBoxInputOnly() const;

   /*!
    * Turn on gradient detector to tag cells for refinement.
    */
   void
   turnOnGradientDetector()
   {
      d_use_gradient_detector = true;
   }

   /*!
    * Turn off gradient detector.
    */
   void
   turnOffGradientDetector()
   {
      d_use_gradient_detector = false;
   }

   /*!
    * Turn on Richardson extrapolation to tag cells for refinement.
    */
   void
   turnOnRichardsonExtrapolation()
   {
      d_use_richardson_extrapolation = true;
   }

   /*!
    * Turn off Richardson extrapolation.
    */
   void
   turnOffRichardsonExtrapolation()
   {
      d_use_richardson_extrapolation = false;
   }

   /*!
    * Turn on static refine box regions where refinement should occur.
    */
   void
   turnOnRefineBoxes()
   {
      d_use_refine_boxes = true;
   }

   /*!
    * Turn off static refine box regions.
    */
   void
   turnOffRefineBoxes()
   {
      d_use_refine_boxes = false;
   }

   /*!
    * Read input values, indicated above, from given database.
    *
    * When assertion checking is active, the database pointer must be non-null.
    */
   void
   getFromInput(
      const boost::shared_ptr<tbox::Database>& db);

   const StandardTagAndInitializeConnectorWidthRequestor&
   getConnectorWidthRequestor() const
   {
      return d_staicwri;
   }

private:
   /*
    * Apply preprocessing for Richardson extrapolation.
    */
   void
   preprocessRichardsonExtrapolation(
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double regrid_time,
      const double regrid_start_time,
      const bool initial_time);

   /*
    * Apply Richardson extrapolation algorithm.
    */
   void
   tagCellsUsingRichardsonExtrapolation(
      const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
      const int level_number,
      const double regrid_time,
      const double regrid_start_time,
      const int tag_index,
      const bool initial_time,
      const bool coarsest_sync_level,
      const bool can_be_refined);

   /*
    * Booleans specifying the tagging method.  Any combination of the
    * three methods may be used.
    */
   bool d_use_refine_boxes;
   bool d_use_gradient_detector;
   bool d_use_richardson_extrapolation;

   /*
    * Concrete object that supplies problem-specific initialization
    * and regridding operations.
    */
   StandardTagAndInitStrategy* d_tag_strategy;

   /*
    * The error_coarsen_ratio used for all levels in the hierarchy.
    * If Richardson extrapolation is not used, the error coarsen ratio
    * is 1.  If Richardson extrapolation is used, the error coarsen ratio
    * is set in the method coarsestLevelBoxesOK().
    */
   int d_error_coarsen_ratio;

   /*
    * tbox::Array of patch levels containing coarsened versions of the patch
    * levels, for use with Richardson extrapolation.
    */
   tbox::Array<boost::shared_ptr<hier::PatchLevel> > d_rich_extrap_coarsened_levels;

   StandardTagAndInitializeConnectorWidthRequestor d_staicwri;

};

}
}

#endif
