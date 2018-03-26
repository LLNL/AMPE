/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Interface to patch routines for hyperbolic integration scheme.
 *
 ************************************************************************/

#ifndef included_algs_HyperbolicPatchStrategy_C
#define included_algs_HyperbolicPatchStrategy_C

#include "SAMRAI/algs/HyperbolicPatchStrategy.h"

#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI {
namespace algs {

HyperbolicPatchStrategy::HyperbolicPatchStrategy(
   const tbox::Dimension& dim):
   xfer::RefinePatchStrategy(dim),
   xfer::CoarsenPatchStrategy(dim),
   d_dim(dim),
   d_data_context()
{
}

HyperbolicPatchStrategy::~HyperbolicPatchStrategy()
{
}

/*
 *************************************************************************
 *
 * Default virtual function implementations.
 *
 *************************************************************************
 */

void
HyperbolicPatchStrategy::tagGradientDetectorCells(
   hier::Patch& patch,
   const double regrid_time,
   const bool initial_error,
   const int tag_index,
   const bool uses_richardson_extrapolation_too)
{
   NULL_USE(patch);
   NULL_USE(regrid_time);
   NULL_USE(initial_error);
   NULL_USE(tag_index);
   NULL_USE(uses_richardson_extrapolation_too);
   TBOX_ERROR("HyperbolicPatchStrategy::tagGradientDetectorCells()"
      << "\nNo derived class supplies a concrete implementation for "
      << "\nthis method." << std::endl);
}

void
HyperbolicPatchStrategy::tagRichardsonExtrapolationCells(
   hier::Patch& patch,
   const int error_level_number,
   const boost::shared_ptr<hier::VariableContext>& coarsened_fine,
   const boost::shared_ptr<hier::VariableContext>& advanced_coarse,
   const double regrid_time,
   const double deltat,
   const int error_coarsen_ratio,
   const bool initial_error,
   const int tag_index,
   const bool uses_gradient_detector_too)
{
   NULL_USE(patch);
   NULL_USE(error_level_number);
   NULL_USE(coarsened_fine);
   NULL_USE(advanced_coarse);
   NULL_USE(regrid_time);
   NULL_USE(deltat);
   NULL_USE(error_coarsen_ratio);
   NULL_USE(initial_error);
   NULL_USE(tag_index);
   NULL_USE(uses_gradient_detector_too);
   TBOX_ERROR("HyperbolicPatchStrategy::tagRichardsonExtrapolationCells()"
      << "\nNo derived class supplies a concrete implementation for "
      << "\nthis method." << std::endl);
}

void
HyperbolicPatchStrategy::setupLoadBalancer(
   HyperbolicLevelIntegrator* integrator,
   mesh::GriddingAlgorithm* gridding_algorithm)
{
   NULL_USE(integrator);
   NULL_USE(gridding_algorithm);
}

void
HyperbolicPatchStrategy::preprocessAdvanceLevelState(
   const boost::shared_ptr<hier::PatchLevel>& level,
   double current_time,
   double dt,
   bool first_step,
   bool last_step,
   bool regrid_advance)
{
   NULL_USE(level);
   NULL_USE(current_time);
   NULL_USE(dt);
   NULL_USE(first_step);
   NULL_USE(last_step);
   NULL_USE(regrid_advance);
}

void
HyperbolicPatchStrategy::postprocessAdvanceLevelState(
   const boost::shared_ptr<hier::PatchLevel>& level,
   double current_time,
   double dt,
   bool first_step,
   bool last_step,
   bool regrid_advance)
{
   NULL_USE(level);
   NULL_USE(current_time);
   NULL_USE(dt);
   NULL_USE(first_step);
   NULL_USE(last_step);
   NULL_USE(regrid_advance);
}

hier::IntVector
HyperbolicPatchStrategy::getRefineOpStencilWidth() const
{
   return hier::IntVector::getZero(d_dim);
}

void
HyperbolicPatchStrategy::preprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio)
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}

void
HyperbolicPatchStrategy::postprocessRefine(
   hier::Patch& fine,
   const hier::Patch& coarse,
   const hier::Box& fine_box,
   const hier::IntVector& ratio)
{
   NULL_USE(fine);
   NULL_USE(coarse);
   NULL_USE(fine_box);
   NULL_USE(ratio);
}

hier::IntVector
HyperbolicPatchStrategy::getCoarsenOpStencilWidth() const
{
   return hier::IntVector::getZero(d_dim);
}

void
HyperbolicPatchStrategy::preprocessCoarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio)
{
   NULL_USE(coarse);
   NULL_USE(fine);
   NULL_USE(coarse_box);
   NULL_USE(ratio);
}

void
HyperbolicPatchStrategy::postprocessCoarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio)
{
   NULL_USE(coarse);
   NULL_USE(fine);
   NULL_USE(coarse_box);
   NULL_USE(ratio);
}

}
}
#endif
