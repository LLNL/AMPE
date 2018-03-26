/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Interface to application-specific patch functions in support
 *                Method of Lines integration algorithm
 *
 ************************************************************************/

#ifndef included_algs_MethodOfLinesPatchStrategy
#define included_algs_MethodOfLinesPatchStrategy

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/xfer/CoarsenPatchStrategy.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"

#include <boost/shared_ptr.hpp>

namespace SAMRAI {
namespace algs {

class MethodOfLinesIntegrator;

/**
 * Class MethodOfLinesPatchStrategy is an abstract type defining the
 * interface for operations invoked during the integration routines defined
 * in the MethodOfLinesIntegrator class.  This class is derived from
 * the xfer::RefinePatchStrategy and xfer::CoarsenPatchStrategy abstract
 * base classes.  These base classes define the interfaces for user-defined
 * interlevel data refining and coarsening operations and the specification
 * of physical boundary conditions.
 *
 * @see algs::MethodOfLinesIntegrator
 * @see xfer::RefinePatchStrategy
 * @see xfer::CoarsenPatchStrategy
 */

class MethodOfLinesPatchStrategy:
   public xfer::RefinePatchStrategy,
   public xfer::CoarsenPatchStrategy
{
public:
   /*!
    * Blank constructor for MethodOfLinesPatchStrategy.
    */
   explicit MethodOfLinesPatchStrategy(
      const tbox::Dimension& dim);

   /*!
    * Virtual destructor for MethodOfLinesPatchStrategy.
    */
   virtual ~MethodOfLinesPatchStrategy() = 0;

   /*!
    * Register variables specific to the problem to be solved with the
    * integrator using the registerVariable function.  This
    * defines the way data for each quantity will be manipulated on the
    * patches.  For more information, refer to
    * MethodOfLinesIntegrator::registerVariable.
    */
   virtual void
   registerModelVariables(
      MethodOfLinesIntegrator* integrator) = 0;

   /*!
    * Set the initial data on a patch interior (i.e., NO GHOST CELLS).
    * Setting "initial_time" true will initialize data at time. Setting it
    * false will interpolate data from the appropriate coarser patch.
    */
   virtual void
   initializeDataOnPatch(
      hier::Patch& patch,
      const double time,
      const bool initial_time) const = 0;

   /*!
    * Compute the stable time increment for a patch.
    */
   virtual double
   computeStableDtOnPatch(
      hier::Patch& patch,
      const double time) const = 0;

   /*!
    * Advance a single Runge Kutta step.
    *
    * @param patch patch that RK step is being applied
    * @param dt    timestep
    * @param alpha_1 first coefficient applied in the RK step
    * @param alpha_2 second coefficient
    * @param beta    third coefficient
    */
   virtual void
   singleStep(
      hier::Patch& patch,
      const double dt,
      const double alpha_1,
      const double alpha_2,
      const double beta) const = 0;

   /*!
    * Using a user-specified gradient detection scheme, determine cells which
    * have high gradients and, consequently, should be refined.
    */
   virtual void tagGradientDetectorCells(
      hier::Patch& patch,
      const double regrid_time,
      const bool initial_error,
      const int tag_index,
      const bool uses_richardson_extrapolation_too) = 0;

   /*!
    * Set user-defined boundary conditions at the physical domain boundary.
    */
   virtual void
   setPhysicalBoundaryConditions(
      hier::Patch& patch,
      const double fill_time,
      const hier::IntVector& ghost_width_to_fill) = 0;

   /*!
    * Return maximum stencil width needed for user-defined
    * data interpolation operations.  Default is to return
    * zero, assuming no user-defined operations provided.
    */
   virtual hier::IntVector
   getRefineOpStencilWidth() const;

   /*!
    * Pre- and post-processing routines for implementing user-defined
    * spatial interpolation routines applied to variables.  The
    * interpolation routines are used in the MOL AMR algorithm
    * for filling patch ghost cells before advancing data on a level
    * and after regridding a level to fill portions of the new level
    * from some coarser level.  These routines are called automatically
    * from within patch boundary filling schedules; thus, some concrete
    * function matching these signatures must be provided in the user's
    * patch model.  However, the routines only need to perform some
    * operations when "USER_DEFINED_REFINE" is given as the interpolation
    * method for some variable when the patch model registers variables
    * with the MOL integration algorithm, typically.  If the
    * user does not provide operations that refine such variables in either
    * of these routines, then they will not be refined.
    *
    * The order in which these operations are used in each patch
    * boundary filling schedule is:
    *
    * - \b (1) {Call user's preprocessRefine() routine.}
    * - \b (2) {Refine all variables with standard interpolation operators.}
    * - \b (3) {Call user's postprocessRefine() routine.}
    *
    *
    * Also, user routines that implement these functions must use
    * data corresponding to the d_scratch context on both coarse and
    * fine patches.
    */
   virtual void
   preprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio) = 0;

   ///
   virtual void
   postprocessRefine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const hier::Box& fine_box,
      const hier::IntVector& ratio) = 0;

   /*!
    * Return maximum stencil width needed for user-defined
    * data coarsen operations.  Default is to return
    * zero, assuming no user-defined operations provided.
    */
   virtual hier::IntVector
   getCoarsenOpStencilWidth() const;

   /*!
    * Pre- and post-processing routines for implementing user-defined
    * spatial coarsening routines applied to variables.  The coarsening
    * routines are used in the MOL AMR algorithm synchronizing
    * coarse and fine levels when they have been integrated to the same
    * point.  These routines are called automatically from within the
    * data synchronization coarsen schedules; thus, some concrete
    * function matching these signatures must be provided in the user's
    * patch model.  However, the routines only need to perform some
    * operations when "USER_DEFINED_COARSEN" is given as the coarsening
    * method for some variable when the patch model registers variables
    * with the MOL level integration algorithm, typically.  If the
    * user does not provide operations that coarsen such variables in either
    * of these routines, then they will not be coarsened.
    *
    * The order in which these operations are used in each coarsening
    * schedule is:
    *
    * - \b (1) {Call user's preprocessCoarsen() routine.}
    * - \b (2) {Coarsen all variables with standard coarsening operators.}
    * - \b (3) {Call user's postprocessCoarsen() routine.}
    *
    *
    * Also, user routines that implement these functions must use
    * corresponding to the d_new context on both coarse and fine patches
    * for time-dependent quantities.
    */
   virtual void
   preprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio) = 0;

   ///
   virtual void
   postprocessCoarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio) = 0;

   /*!
    * The method of lines integrator controls the context for the data to
    * be used in the numerical routines implemented in the concrete patch
    * strategy. These accessor methods allow the patch strategy to access
    * the particular data contexts used in the integrator.
    *
    * Return pointer to data context with ghost cells.
    */
   boost::shared_ptr<hier::VariableContext>
   getInteriorWithGhostsContext() const
   {
      return d_interior_with_ghosts;
   }

   /*!
    * Return pointer to data context with NO ghosts.
    */
   boost::shared_ptr<hier::VariableContext>
   getInteriorContext() const
   {
      return d_interior;
   }

   /*!
    * Set pointer to data context with ghosts.
    */
   void
   setInteriorWithGhostsContext(
      const boost::shared_ptr<hier::VariableContext>& context)
   {
      d_interior_with_ghosts = context;
   }

   /*!
    * Set pointer to data context with NO ghosts.
    */
   void
   setInteriorContext(
      const boost::shared_ptr<hier::VariableContext>& context)
   {
      d_interior = context;
   }

   const tbox::Dimension&
   getDim() const
   {
      return d_dim;
   }

private:
   const tbox::Dimension d_dim;

   boost::shared_ptr<hier::VariableContext> d_interior_with_ghosts;
   boost::shared_ptr<hier::VariableContext> d_interior;
};

}
}

#endif
