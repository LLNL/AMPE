/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Routines for performing cell-tagging and initializing
 *                a new level.
 *
 ************************************************************************/

#ifndef included_mesh_StandardTagAndInitialize_C
#define included_mesh_StandardTagAndInitialize_C

#include "SAMRAI/mesh/StandardTagAndInitialize.h"

#include "SAMRAI/pdat/CellIntegerConstantRefine.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/OverlapConnectorAlgorithm.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <stdio.h>

/*
 *************************************************************************
 *
 * External declarations for FORTRAN 77 routines used in Richardson
 * extrapolation algorithm to coarsen tagged cells from fine to coarse
 * level.
 *
 *************************************************************************
 */

extern "C" {
#ifdef __INTEL_COMPILER
#pragma warning (disable:1419)
#endif
// in coarsentags1d.f:
void F77_FUNC(coarsentags1d, COARSENTAGS1D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&,
   const int *,
   const int *, int *);
// in coarsentags2d.f:
void F77_FUNC(coarsentags2d, COARSENTAGS2D) (const int&, const int&,
   const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int&, const int&, const int&, const int&,
   const int *,
   const int *, int *);
// in coarsentags3d.f:
void F77_FUNC(coarsentags3d, COARSENTAGS3D) (const int&, const int&,
   const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int&, const int&, const int&,
   const int *,
   const int *, int *);
}

#define DEBUG_TIMES
//#undef DEBUG_TIMES

namespace SAMRAI {
namespace mesh {

// using namespace std;

/*
 *************************************************************************
 *
 * Static function that computes greatest common divisor.
 *
 *************************************************************************
 */

static int
GCD(
   const int a,
   const int b);

/*
 *************************************************************************
 *
 * Constructors and destructor for StandardTagAndInitialize.
 *
 *************************************************************************
 */

StandardTagAndInitialize::StandardTagAndInitialize(
   const tbox::Dimension& dim,
   const std::string& object_name,
   StandardTagAndInitStrategy* tag_strategy,
   const boost::shared_ptr<tbox::Database>& input_db):
   TagAndInitializeStrategy(dim, object_name)
{
   TBOX_ASSERT(!object_name.empty());

   d_tag_strategy = tag_strategy;

   d_error_coarsen_ratio = 1;

   d_use_gradient_detector = false;
   d_use_richardson_extrapolation = false;
   d_use_refine_boxes = false;

   /*
    * If no input database is provided, no criteria is set to tag cells
    * so cell-tagging will not occur.  Print a warning to indicate if
    * this is the case.
    */
   if (!input_db) {
      TBOX_WARNING(
         getObjectName() << ":constructor \n"
                         << "no input database specified - NO METHOD IS SPECIFIED TO TAG \n"
                         << "CELLS FOR REFINEMENT so no tagging is performed.");
   } else {
      getFromInput(input_db);
   }

   /*
    * If the user wishes to only use the REFINE_BOXES tagging option,
    * the registered strategy class may be null.  In order to use
    * the GRADIENT_DETECTOR or RICHARDSON_EXTRAPOLATION options, the
    * registered StandardTagAndInitStrategy must be non-NULL.
    */
   if (d_use_gradient_detector || d_use_richardson_extrapolation) {
      if (tag_strategy == ((StandardTagAndInitStrategy *)NULL)) {
         TBOX_ERROR(
            getObjectName() << ":constructor "
                            << "\nThe supplied implementation of the "
                            << "\nStandardTagAndInitStrategy is NULL.  It must be"
                            << "\nnon-NULL to use the GRADIENT_DETECTOR or"
                            << "\nRICHARDSON_EXTRAPOLATION tagging options." << std::endl);
      }
   }
}

StandardTagAndInitialize::~StandardTagAndInitialize()
{
}

/*
 *************************************************************************
 *
 * Pass requests to initialize level data, reset hierarchy information,
 * and apply an application-specific gradient detector to
 * the subclass of the StandardTagAndInitStrategX data member.
 *
 *************************************************************************
 */

void
StandardTagAndInitialize::initializeLevelData(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const double init_data_time,
   const bool can_be_refined,
   const bool initial_time,
   const boost::shared_ptr<hier::PatchLevel>& old_level,
   const bool allocate_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= hierarchy->getFinestLevelNumber()));
   if (old_level) {
      TBOX_ASSERT(level_number == old_level->getLevelNumber());
   }
   TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
#endif
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(getDim(), *hierarchy);

   if (d_tag_strategy != ((StandardTagAndInitStrategy *)NULL)) {
      d_tag_strategy->initializeLevelData(hierarchy,
         level_number,
         init_data_time,
         can_be_refined,
         initial_time,
         old_level,
         allocate_data);
   }

}

/*
 *************************************************************************
 *
 * Reset hierarchy configuration information where the range of new
 * hierarchy levels is specified.
 *
 *************************************************************************
 */

void
StandardTagAndInitialize::resetHierarchyConfiguration(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int coarsest_level,
   const int finest_level)
{

   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((coarsest_level >= 0)
      && (coarsest_level <= finest_level)
      && (finest_level <= hierarchy->getFinestLevelNumber()));
#ifdef DEBUG_CHECK_ASSERTIONS
   for (int ln0 = 0; ln0 <= finest_level; ln0++) {
      TBOX_ASSERT(hierarchy->getPatchLevel(ln0));
   }
#endif
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(getDim(), *hierarchy);

   if (d_tag_strategy != ((StandardTagAndInitStrategy *)NULL)) {
      d_tag_strategy->resetHierarchyConfiguration(hierarchy,
         coarsest_level,
         finest_level);
   }

}

/*
 *************************************************************************
 *
 * Tag cells on level where refinement should occur.   The method can
 * tag cells using either of three options:
 *
 *    1) Richardson extrapolation
 *    2) gradient detection
 *    3) user supplied refine boxes.
 *
 * These options may be used individually or in combination.  If used in
 * combination,  it is IMPORTANT TO PRESERVE THE ORDER of the calls
 * (Richardson extrapolation 1st, gradient detection 2nd, user-supplied
 * refine boxes 3rd) in this method because users may have logic in
 * their code to compare how cells are tagged and changing the order
 * could destroy this logic.
 *
 *************************************************************************
 */

void
StandardTagAndInitialize::tagCellsForRefinement(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const double regrid_time,
   const int tag_index,
   const bool initial_time,
   const bool coarsest_sync_level,
   const bool can_be_refined,
   const double regrid_start_time)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
   TBOX_ASSERT(tag_index >= 0);
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(getDim(), *hierarchy);

   if (d_use_richardson_extrapolation) {
      tagCellsUsingRichardsonExtrapolation(hierarchy,
         level_number,
         regrid_time,
         regrid_start_time,
         tag_index,
         initial_time,
         coarsest_sync_level,
         can_be_refined);
   }

   if (d_use_gradient_detector) {

      NULL_USE(regrid_start_time);
      NULL_USE(can_be_refined);
      NULL_USE(coarsest_sync_level);

      TBOX_ASSERT(d_tag_strategy != ((StandardTagAndInitStrategy *)NULL));

      d_tag_strategy->applyGradientDetector(hierarchy,
         level_number,
         regrid_time,
         tag_index,
         initial_time,
         d_use_richardson_extrapolation);
   }

   /*
    * If user-supplied refine boxes are to be used, get refine box information
    * from the TagAndInitializeStrategy class, from which this class is
    * derived.
    */
   if (d_use_refine_boxes) {

      hier::BoxContainer refine_boxes;
      getUserSuppliedRefineBoxes(refine_boxes, level_number, regrid_time);

      boost::shared_ptr<hier::PatchLevel> level(
         hierarchy->getPatchLevel(level_number));

      for (hier::PatchLevel::iterator ip(level->begin());
           ip != level->end(); ++ip) {
         const boost::shared_ptr<hier::Patch>& patch = *ip;

         boost::shared_ptr<pdat::CellData<int> > tag_data(
            patch->getPatchData(tag_index),
            boost::detail::dynamic_cast_tag());

         TBOX_ASSERT(tag_data);

         for (hier::BoxContainer::iterator ib(refine_boxes);
              ib != refine_boxes.end(); ++ib) {
            hier::Box intersection = *ib * tag_data->getBox();
            if (!(intersection.empty())) {
               tag_data->fill(1, intersection);
            }
         }
      }
   }

}

/*
 *************************************************************************
 *
 * The Richardson extrapolation error estimation tags cells according
 * to differences in the solution computed on two different levels of
 * refinement.  The preprocessRichardsonExtrapolation method advanced
 * data on a COARSENED VERSION of the hierarchy level where regridding
 * is applied.  This method advances data on the level itself and
 * compares the solutions on the level and the coarsened version of the
 * level at the advanced time.
 *
 * The steps are summarized as follows:
 *
 * 0) Advance data on the level
 *    0a) if (initial_time) {
 *           advance level for ErrorCoarsenRatio steps with time
 *           increment dt
 *        } else {
 *           advance level by 1 step (see discussion
 *           under 2b for reasons for the difference).
 *        }
 *        NOTE: The "first_step" argument in the
 *        tag_strategy->advanceLevel() method sets the conditions of
 *        the advance.  The following conditions define the state of
 *        "first_step":
 *           - first_step is always true at the initial time.
 *           - at subsequent times, it is true only if:
 *                    level <=  finest level that has not been regridded
 *              -AND- the level is the coarsest used in flux syncs
 *    0b) reset the time dependent data for all but the last step.  The
 *        reason we don't do it on the last step is that we are going to
 *        use the allocated space for step 4.
 *
 * 1) Coarsen data computed on the level to the coarsened version of the
 *    level.  Apply the
 *    tag_strategy->coarsenDataForRichardsonExtrapolation() method again.
 *    This time, before_advance set to false since data has already been
 *    advanced on both levels.
 *
 * 2) Allocate tags and apply the
 *    tag_strategy->applyRichardsonExtrapolation() method.  This method
 *    sets tags on the COARSENED VERSION of the level according to the
 *    criteria set in the tag_strategy.
 *
 * 3) Refine tags from the coarsened version of the level to the level
 *    where tagging is performed.
 *
 * 4) Put data back into the state it was before calling the RE routine.
 *    4a) Apply the d_tag_strategy->resetDataToPreAdvanceState() to take
 *        care of resetting data after the last advance step (see 3b
 *        discussion).
 *    4b) Reset the timestamp of the data if we are at the initial time
 *        by calling the initializeLevelData() method with
 *        "allocate_data"
 *        set to false. We have already allocated and operated on the
 *        initialized data, we just want to reset the timestamp of it.
 *
 *************************************************************************
 */

void
StandardTagAndInitialize::tagCellsUsingRichardsonExtrapolation(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const double regrid_time,
   const double regrid_start_time,
   const int tag_index,
   const bool initial_time,
   const bool coarsest_sync_level,
   const bool can_be_refined)
{
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(getDim(), *hierarchy);
   TBOX_ASSERT(regrid_start_time <= regrid_time);
   TBOX_ASSERT(d_tag_strategy != ((StandardTagAndInitStrategy *)NULL));

   const tbox::Dimension& dim(getDim());

   boost::shared_ptr<hier::PatchLevel> patch_level(
      hierarchy->getPatchLevel(level_number));

   /*
    * Determine the level timestep.  If the error coarsen ratio is 2
    * (i.e. even refinement ratio) only a single step has been taken
    * between the regrid_start_time and regrid_time, so the timestep is
    * the difference between the regrid_time and the regrid_start_time.
    * If the error coarsen ratio is 3 (i.e. factor 3 refine ratio)
    * then two timesteps have been taken between the regrid_start_time
    * and the regrid_time, so the timestep will be one-half the difference
    * between the regrid_start_time and the regrid_time.
    */
   double dt = (regrid_time - regrid_start_time)
      / (double)(d_error_coarsen_ratio - 1);

   /*
    * Determine number of advance steps for time integration on the level.
    */
   int n_steps;
   if (initial_time) {
      n_steps = d_error_coarsen_ratio;
   } else {
      n_steps = 1;
      d_tag_strategy->resetTimeDependentData(patch_level,
         regrid_time,
         can_be_refined);
   }

   /*
    * Arguments to the advanceLevel() method
    * are set as follows:
    *   first_step - if the level is NOT the coarsest sync level (that is
    *                the coarsest level to synchronize at the regridding
    *                time), first step is true.  Otherwise, it is false.
    *   last_step - true or false: depends on step count
    *   regrid_advance - true: this is a time-dependent regrid advance
    */

   bool first_step = !coarsest_sync_level;
   bool regrid_advance = true;

   double start_time = regrid_time;
   double end_time = 0.0;
   for (int step_cnt = 0; step_cnt < n_steps; step_cnt++) {

      end_time = start_time + dt;
      bool last_step = (step_cnt == (n_steps - 1));

#ifdef DEBUG_TIMES
      tbox::plog << "\nAdvancing Data on level in Rich. Extrap" << std::endl;
      tbox::plog << "level number = " << patch_level->getLevelNumber()
                 << std::endl;
      tbox::plog << "level in hierarchy? " << patch_level->inHierarchy()
                 << std::endl;
      tbox::plog << "start time = " << start_time << std::endl;
      tbox::plog << "end time = " << end_time << std::endl;
      tbox::plog << "first step? = " << first_step << std::endl;
      tbox::plog << "last step? = " << last_step << std::endl;
      tbox::plog << "regrid advance? = " << regrid_advance << std::endl;
#endif

      (void)d_tag_strategy->advanceLevel(patch_level,
         hierarchy,
         start_time,
         end_time,
         first_step,
         last_step,
         regrid_advance);

      if (step_cnt < (n_steps - 1)) {
         d_tag_strategy->resetTimeDependentData(patch_level,
            end_time,
            can_be_refined);
      }

      start_time = end_time;
   }

   boost::shared_ptr<hier::PatchLevel> coarser_level(
      d_rich_extrap_coarsened_levels[level_number]);

   /*
    * Coarsen data from hierarchy level to coarser level.
    */
   bool before_advance = false;
   d_tag_strategy->coarsenDataForRichardsonExtrapolation(hierarchy,
      level_number,
      coarser_level,
      end_time,
      before_advance);

   coarser_level->allocatePatchData(tag_index, end_time);

   /*
    * Coarsen tags from level to coarser level.
    */
   hier::IntVector coarsen_ratio(dim, d_error_coarsen_ratio);
   for (hier::PatchLevel::iterator ip(coarser_level->begin());
        ip != coarser_level->end(); ++ip) {
      const boost::shared_ptr<hier::Patch>& coarse_patch = *ip;
      boost::shared_ptr<hier::Patch> fine_patch(
         patch_level->getPatch(coarse_patch->getGlobalId()));
      boost::shared_ptr<pdat::CellData<int> > ftags(
         fine_patch->getPatchData(tag_index),
         boost::detail::dynamic_cast_tag());
      boost::shared_ptr<pdat::CellData<int> > ctags(
         coarse_patch->getPatchData(tag_index),
         boost::detail::dynamic_cast_tag());

      TBOX_ASSERT(ftags);
      TBOX_ASSERT(ctags);
      TBOX_ASSERT(ctags->getDepth() == ftags->getDepth());

      const hier::Index filo = ftags->getGhostBox().lower();
      const hier::Index fihi = ftags->getGhostBox().upper();
      const hier::Index cilo = ctags->getGhostBox().lower();
      const hier::Index cihi = ctags->getGhostBox().upper();

      const hier::Index ifirstc = coarse_patch->getBox().lower();
      const hier::Index ilastc = coarse_patch->getBox().upper();

      for (int d = 0; d < ctags->getDepth(); d++) {
         if (dim == tbox::Dimension(1)) {
            F77_FUNC(coarsentags1d, COARSENTAGS1D) (ifirstc(0), ilastc(0),
               filo(0), fihi(0),
               cilo(0), cihi(0),
               &coarsen_ratio[0],
               ftags->getPointer(d),
               ctags->getPointer(d));
         } else if ((dim == tbox::Dimension(2))) {
            F77_FUNC(coarsentags2d, COARSENTAGS2D) (ifirstc(0), ifirstc(1),
               ilastc(0), ilastc(1),
               filo(0), filo(1), fihi(0), fihi(1),
               cilo(0), cilo(1), cihi(0), cihi(1),
               &coarsen_ratio[0],
               ftags->getPointer(d),
               ctags->getPointer(d));
         } else if ((dim == tbox::Dimension(3))) {
            F77_FUNC(coarsentags3d, COARSENTAGS3D) (ifirstc(0), ifirstc(1),
               ifirstc(2),
               ilastc(0), ilastc(1), ilastc(2),
               filo(0), filo(1), filo(2),
               fihi(0), fihi(1), fihi(2),
               cilo(0), cilo(1), cilo(2),
               cihi(0), cihi(1), cihi(2),
               &coarsen_ratio[0],
               ftags->getPointer(d),
               ctags->getPointer(d));

         } else {
            TBOX_ERROR("StandardTagAndInitialize error...\n"
               << "DIM > 3 not supported." << std::endl);

         }

      }
   }

   /*
    * Tag cells on coarser level.
    */
   d_tag_strategy->applyRichardsonExtrapolation(coarser_level,
      end_time,
      tag_index,
      dt,
      d_error_coarsen_ratio,
      initial_time,
      d_use_gradient_detector);

   /*
    * Refine tags from coarser level to level.
    */
   pdat::CellIntegerConstantRefine copytags(dim);
   for (hier::PatchLevel::iterator ip(coarser_level->begin());
        ip != coarser_level->end(); ++ip) {

      const boost::shared_ptr<hier::Patch>& coarse_patch = *ip;
      boost::shared_ptr<hier::Patch> fine_patch(
         patch_level->getPatch(coarse_patch->getGlobalId()));
      copytags.refine(*fine_patch, *coarse_patch,
         tag_index, tag_index,
         fine_patch->getBox(), coarsen_ratio);
   }

   /*
    * Final cleanup.  Reset data to initial state before entering this routine.
    */
   d_tag_strategy->resetDataToPreadvanceState(patch_level);

   if (initial_time) {
      bool allocate_data = false;
      initializeLevelData(hierarchy, level_number, regrid_time,
         can_be_refined, initial_time,
         boost::shared_ptr<hier::PatchLevel>(),
         allocate_data);
   }

}

/*
 *************************************************************************
 *
 * Preprocess data before cell tagging, if appropriate.  For the options
 * provided in this class, only Richardson extrapolation requires
 * any pre-processing.
 *
 *************************************************************************
 */

void
StandardTagAndInitialize::preprocessErrorEstimation(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const double regrid_time,
   const double regrid_start_time,
   const bool initial_time)
{
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((level_number >= 0)
      && (level_number <= hierarchy->getFinestLevelNumber()));
   TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(getDim(), *hierarchy);

   if (d_use_richardson_extrapolation) {
      preprocessRichardsonExtrapolation(hierarchy,
         level_number,
         regrid_time,
         regrid_start_time,
         initial_time);
   }
}

/*
 *************************************************************************
 *
 * The preprocess method for Richardson extrapolation error estimation
 * creates a coarsened version of a level in the hierarchy and advances
 * data on the coarsened level to a prescribed time.
 *
 * The steps are summarized as follows:
 *
 * 0) Create a coarser version of the patch level where tagging is
 *    being performed. The coarser level is coarsened by the
 *    "error coarsen ratio", which is the greatest common divisor of the
 *    refinement ratio (e.g. GCD of ratio 4 refinement would be 2).
 *
 * 1) Initialize data on the coarser level by applying the
 *    tag_strategy->coarsenDataForRichardsonExtrapolation() method.
 *    Note that "before_advance" is set true in this call since we have
 *    not yet advanced data on the coarser level.
 *
 * 2) Advance data on the coarsened version of the level:
 *    2a) get timestep (dt) for level where tagging is performed
 *    2b) if (initial_time) {
 *          Advance coarse level by ErrorCoarsenRatio*dt
 *        } else {
 *          Advance coarse level by dt
 *        }
 *    2c) reset the time dependent data on the coarser level by calling
 *        the tag_strategy's resetTimeDependentData() function.
 *
 * The constructed coarsened levels are stored and used in the
 * tagCellsUsingRichardsonExtrapolation() method.
 *
 *************************************************************************
 */

void
StandardTagAndInitialize::preprocessRichardsonExtrapolation(
   const boost::shared_ptr<hier::PatchHierarchy>& hierarchy,
   const int level_number,
   const double regrid_time,
   const double regrid_start_time,
   const bool initial_time)
{
   TBOX_ASSERT(hierarchy);
   TBOX_DIM_ASSERT_CHECK_DIM_ARGS1(getDim(), *hierarchy);
   TBOX_ASSERT(regrid_start_time <= regrid_time);
   TBOX_ASSERT(d_tag_strategy != ((StandardTagAndInitStrategy *)NULL));

   const tbox::Dimension& dim(getDim());

   boost::shared_ptr<hier::PatchLevel> patch_level(
      hierarchy->getPatchLevel(level_number));

   /*
    * Determine the level timestep.  If the error coarsen ratio is 2
    * (i.e. even refinement ratio) only a single step has been taken
    * between the regrid_start_time and regrid_time, so the timestep is
    * the difference between the regrid_time and the regrid_start_time.
    * If the error coarsen ratio is 3 (i.e. factor 3 refine ratio)
    * then two timesteps have been taken between the regrid_start_time
    * and the regrid_time, so the timestep will be one-half the difference
    * between the regrid_start_time and the regrid_time.
    */
   double dt = (regrid_time - regrid_start_time)
      / (double)(d_error_coarsen_ratio - 1);

   /*
    * Determine start and end times for integration on the coarsened level.
    *
    * At the initial time, the start time is the regrid_time and the end
    * time is the regrid_time + error_coarsen_ratio*dt.
    *
    * In subsequent to the initial time, the start time is the
    * regrid_start_time and the end time is the regrid_time + dt.
    */
   double coarse_start_time;
   double coarse_end_time;

   if (initial_time) {
      coarse_start_time = regrid_time;
      coarse_end_time =
         regrid_time + dt * d_error_coarsen_ratio;
   } else {
      coarse_start_time = regrid_start_time;
      coarse_end_time = regrid_time + dt;
   }

#ifdef DEBUG_TIMES
   tbox::plog << "\nRegridding on level using Rich. Extrap" << std::endl;
   tbox::plog << "level number = " << patch_level->getLevelNumber()
              << std::endl;
   tbox::plog << "level in hierarchy? " << patch_level->inHierarchy()
              << std::endl;
   tbox::plog << "coarsened hier level = " << level_number - 1 << std::endl;
   tbox::plog << "coarsened start time = " << coarse_start_time << std::endl;
   tbox::plog << "coarsened end time = " << coarse_end_time << std::endl;
#endif

   /*
    * Generate coarsened version and initialize data on it.  If coarsened
    * level aligns with next coarsened level in hierarchy, set level number
    * so user routines can use this information.
    */

   boost::shared_ptr<hier::PatchLevel> coarsened_level(
      boost::make_shared<hier::PatchLevel>(dim));
   hier::IntVector coarsen_ratio(dim, d_error_coarsen_ratio);
   coarsened_level->setCoarsenedPatchLevel(patch_level, coarsen_ratio);

   if ((level_number > 0)
       && (hierarchy->getPatchLevel(level_number - 1)->getRatioToLevelZero() ==
           coarsened_level->getRatioToLevelZero())) {
      coarsened_level->setLevelNumber(level_number - 1);
   }
   coarsened_level->setNextCoarserHierarchyLevelNumber(level_number - 1);

   /*
    * Generate Connector patch_level<==>coarsened_level and
    * coarsened_level--->coarsened_level.  To support recursive data
    * transfer, these Connectors should have gcw equivalent to
    * patch_level<==>patch_level.
    */
   const hier::IntVector level_to_level_gcw =
      hierarchy->getRequiredConnectorWidth(level_number,
         level_number);

   const hier::Connector level_to_level =
      patch_level->getBoxLevel()->getPersistentOverlapConnectors().
      findConnector(
         *patch_level->getBoxLevel(),
         level_to_level_gcw);

   hier::Connector coarsened_to_level = level_to_level;
   coarsened_to_level.setBase(*coarsened_level->getBoxLevel());
   coarsened_to_level.setHead(*patch_level->getBoxLevel());
   coarsened_to_level.setWidth(
      hier::IntVector::ceilingDivide(level_to_level_gcw, coarsen_ratio),
      true);

   hier::Connector tmp_coarsened(level_to_level);
   tmp_coarsened.setBase( *patch_level->getBoxLevel());
   tmp_coarsened.setHead(*coarsened_level->getBoxLevel(), true);
   tmp_coarsened.coarsenLocalNeighbors(coarsen_ratio);
   tmp_coarsened.setConnectorType(hier::Connector::COMPLETE_OVERLAP);

   const hier::Connector& level_to_coarsened =
      patch_level->getBoxLevel()->getPersistentOverlapConnectors().
      createConnector(
         *coarsened_level->getBoxLevel(),
         level_to_level_gcw,
         tmp_coarsened);

   coarsened_level->getBoxLevel()->getPersistentOverlapConnectors().
   createConnector(
      *coarsened_level->getBoxLevel(),
      hier::IntVector::ceilingDivide(level_to_level_gcw, coarsen_ratio),
      tmp_coarsened);

   if (level_number > 0) {
      /*
       * Get Connectors coarsened<==>coarser, which are used for recursive
       * refinement filling of the coarsened level's ghosts.
       */
      hier::Connector* coarsened_to_coarser = new hier::Connector;
      hier::Connector* coarser_to_coarsened = new hier::Connector;
      boost::shared_ptr<hier::PatchLevel> coarser_level(
         hierarchy->getPatchLevel(level_number - 1));
      const hier::Connector& level_to_coarser =
         patch_level->getBoxLevel()->getPersistentOverlapConnectors().
         findConnector(
            *coarser_level->getBoxLevel(),
            hierarchy->getRequiredConnectorWidth(
               level_number, level_number - 1));
      const hier::Connector& coarser_to_level =
         coarser_level->getBoxLevel()->getPersistentOverlapConnectors()
         .findConnector(
            *patch_level->getBoxLevel(),
            hierarchy->getRequiredConnectorWidth(
               level_number - 1, level_number));
      hier::OverlapConnectorAlgorithm oca;
      oca.bridge(*coarsened_to_coarser,
         *coarser_to_coarsened,
         coarsened_to_level,
         level_to_coarser,
         coarser_to_level,
         level_to_coarsened);
      coarsened_to_coarser->setConnectorType(hier::Connector::COMPLETE_OVERLAP);
      coarser_to_coarsened->setConnectorType(hier::Connector::COMPLETE_OVERLAP);
      coarsened_level->getBoxLevel()->getPersistentOverlapConnectors().
      cacheConnector(
         *coarser_level->getBoxLevel(),
         coarsened_to_coarser);
      coarser_level->getBoxLevel()->getPersistentOverlapConnectors().
      cacheConnector(
         *coarsened_level->getBoxLevel(),
         coarser_to_coarsened);
   }

   bool before_advance = true;
   d_tag_strategy->coarsenDataForRichardsonExtrapolation(hierarchy,
      level_number,
      coarsened_level,
      coarse_start_time,
      before_advance);

   /*
    * Advance data on coarsened level.  Arguments to the advanceLevel() method
    * are set as follows:
    *   first_step - true: this is the first step on the coarsened level
    *                so it is necessary to do any required
    *                setup in the advance method.
    *   last_step - true: only one step will occur on the the coarsened level
    *   regrid_advance - true: this is a time-dependent regrid advance
    */
   bool first_step = true;
   bool last_step = true;
   bool regrid_advance = true;

#ifdef DEBUG_TIMES
   tbox::plog << "\nAdvancing Data on coarsened in Rich. Extrap" << std::endl;
   tbox::plog << "level number = " << coarsened_level->getLevelNumber()
              << std::endl;
   tbox::plog << "level in hierarchy? " << patch_level->inHierarchy()
              << std::endl;
   tbox::plog << "start time = " << coarse_start_time << std::endl;
   tbox::plog << "end time = " << coarse_end_time << std::endl;
   tbox::plog << "first step? = " << first_step << std::endl;
   tbox::plog << "last step? = " << last_step << std::endl;
   tbox::plog << "regrid advance? = " << regrid_advance << std::endl;
#endif

   (void)d_tag_strategy->advanceLevel(coarsened_level,
      hierarchy,
      coarse_start_time,
      coarse_end_time,
      first_step,
      last_step,
      regrid_advance);

   /*
    * Reset data on the coarse level.  Since this level does not reside
    * in the hierarchy, it cannot be refined (and no special storage
    * manipulation is needed; so we set this to false.
    */
   bool level_can_be_refined = false;
   d_tag_strategy->resetTimeDependentData(coarsened_level,
      coarse_end_time,
      level_can_be_refined);

   /*
    * Add the constructed coarsened level to the array of maintained
    * coarsened levels for Richardson extrapolation.
    */
   if (d_rich_extrap_coarsened_levels.getSize() < level_number + 1) {
      d_rich_extrap_coarsened_levels.resizeArray(level_number + 1);
   }

   d_rich_extrap_coarsened_levels[level_number] = coarsened_level;

}

/*
 *************************************************************************
 *
 * Boxes on the coarsest level must be able to be coarsened by the
 * error coarsen ratio to apply Richardson Extrapolation.  This method
 * simply checks that this is the case.
 *
 *************************************************************************
 */

bool
StandardTagAndInitialize::coarsestLevelBoxesOK(
   const hier::BoxContainer& boxes) const
{
   TBOX_ASSERT(boxes.size() > 0);

   bool boxes_ok = true;
   if (d_use_richardson_extrapolation) {

      for (hier::BoxContainer::const_iterator ib(boxes);
           ib != boxes.end(); ++ib) {
         hier::IntVector n_cells = ib->numberCells();
         for (int i = 0; i < getDim().getValue(); i++) {
            int error_coarsen_ratio = getErrorCoarsenRatio();
            if (!((n_cells(i) % error_coarsen_ratio) == 0)) {
               tbox::perr << "Bad domain box: " << *ib << std::endl;
               TBOX_WARNING(
                  getObjectName() << "At least one box on the \n"
                                  << "coarsest level could not be coarsened by the ratio: "
                                  << error_coarsen_ratio);
               boxes_ok = false;
            }
         }
      }
   }
   return boxes_ok;

}

/*
 *************************************************************************
 *
 * Compute Error coarsen ratio for Richardson extrapolation. For a given
 * level, the error coarsen ratio should be the greatest common divisor
 * (GCD) of the refinement ratio applied to the level.  This value
 * should generally be 2 or 3 (e.g. refinement ratio=2 gives GCD=2;
 * rr=3 gives GCD=3; rr=4 gives GCD=2; etc.).
 *
 *************************************************************************
 */

void
StandardTagAndInitialize::checkCoarsenRatios(
   const tbox::Array<hier::IntVector>& ratio_to_coarser)
{
   if (d_use_richardson_extrapolation) {

      /*
       * Compute GCD on first dimension of level 1
       */
      int error_coarsen_ratio = 0;
      int gcd_level1 = ratio_to_coarser[1](0);
      if ((gcd_level1 % 2) == 0) {
         error_coarsen_ratio = 2;
      } else if ((gcd_level1 % 3) == 0) {
         error_coarsen_ratio = 3;
      } else {
         TBOX_ERROR("Unable to perform Richardson extrapolation algorithm "
            << "with ratio_to_coarser[1](0) = " << gcd_level1);
      }

      /*
       * Iterate through levels and check the coarsen ratios to make sure the
       * error coarsen ratios computed in every dimension on every
       * level are between the supported 2 or 3, and that the error coarsen
       * ratios are constant over the hierarchy.
       */
      for (int ln = 1; ln < ratio_to_coarser.getSize(); ln++) {

         for (int d = 0; d < getDim().getValue(); d++) {
            int gcd = GCD(error_coarsen_ratio, ratio_to_coarser[ln](d));
            if ((gcd % error_coarsen_ratio) != 0) {
               gcd = ratio_to_coarser[ln](d);
               TBOX_ERROR(
                  getObjectName() << "\n"
                                  << "Unable to perform Richardson extrapolation because\n"
                                  << "the error coarsen ratio computed from the\n"
                                  << "ratio_to_coarser entries is not constant across all\n"
                                  << "levels, in all dimensions, of the hierarchy. In\n"
                                  << "order to use Richardson extrapolation, the minimum\n"
                                  << "divisor (> 1) of all the ratio_to_coarser entries must\n"
                                  << "be 2 -or- 3:\n"
                                  << "   level 1(0): minimum divisor: "
                                  << error_coarsen_ratio
                                  << "\n   level " << ln << "(" << d
                                  << "):"
                                  << ": ratio_to_coarser = " << gcd);
            }
         }
      }

      d_error_coarsen_ratio = error_coarsen_ratio;

   }

}

bool
StandardTagAndInitialize::usesTimeIntegration() const
{
   return d_use_richardson_extrapolation;
}

int
StandardTagAndInitialize::getErrorCoarsenRatio() const
{
   return d_error_coarsen_ratio;
}

bool
StandardTagAndInitialize::refineUserBoxInputOnly() const
{
   bool use_only_refine_boxes = false;
   if (d_use_refine_boxes) {
      use_only_refine_boxes = true;
      if (d_use_gradient_detector || d_use_richardson_extrapolation) {
         use_only_refine_boxes = false;
      }
   }
   return use_only_refine_boxes;
}

/*
 *************************************************************************
 *
 * Read cell tagging option and, if required, specified refinement boxes.
 *
 *************************************************************************
 */

void
StandardTagAndInitialize::getFromInput(
   const boost::shared_ptr<tbox::Database>& db)
{
   TBOX_ASSERT(db);

   tbox::Array<std::string> tagging_method;
   if (db->keyExists("tagging_method")) {
      tagging_method = db->getStringArray("tagging_method");
   }

   if (tagging_method.getSize() > 3) {
      TBOX_ERROR(getObjectName() << ":getFromInput\n"
                                 << tagging_method.getSize()
                                 << "entries specified"
                                 << "in `tagging_method' input.  Maximum allowable is 3.");
   }

   d_use_gradient_detector = false;
   d_use_richardson_extrapolation = false;
   d_use_refine_boxes = false;

   /*
    * Check tagging method input.
    */

   bool found_method = false;
   for (int i = 0; i < tagging_method.getSize(); i++) {

      if (tagging_method[i] == "GRADIENT_DETECTOR") {

         d_use_gradient_detector = true;
         found_method = true;

      }

      if (tagging_method[i] == "RICHARDSON_EXTRAPOLATION") {

         d_use_richardson_extrapolation = true;
         found_method = true;

      }

      if (tagging_method[i] == "REFINE_BOXES") {

         d_use_refine_boxes = true;
         found_method = true;

      }
   }

   /*
    * Check for valid entries
    */
   if (!found_method) {
      TBOX_WARNING(
         getObjectName() << ":getFromInput \n"
                         << "No `tagging_method' entry specified, so cell tagging \n"
                         << "will NOT be performed.  If you wish to invoke cell \n"
                         << "tagging, you must enter one or more valid tagging \n"
                         << "methods, of type GRADIENT_DETECTOR, "
                         << "RICHARDSON_EXTRAPOLATION, or REFINE_BOXES\n"
                         << "See class header for details.\n");
   }

   /*
    * If user-supplied refine boxes are to be used, get refine box information
    * from input using the TagAndInitializeStrategy class, from which
    * this class is derived.
    */
   if (d_use_refine_boxes) {
      TagAndInitializeStrategy::getFromInput(db);
   }

}

static int GCD(
   const int a,
   const int b)
{
   int at = tbox::MathUtilities<int>::Min(a, b);
   int bt = tbox::MathUtilities<int>::Max(a, b);

   if (at == 0 || bt == 0) return bt;

   at = (at > 0 ? at : -at);
   bt = (bt > 0 ? bt : -bt);

   int r0 = at;
   int r1 = at;
   int r2 = bt;

   while (!(r2 == 0)) {
      r0 = r1;
      r1 = r2;
      int q = r0 / r1;
      r2 = r0 - r1 * q;
   }

   return r1;
}

}
}
#endif
