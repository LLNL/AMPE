/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2012 Lawrence Livermore National Security, LLC
 * Description:   Interface to patch routines for hyperbolic integration scheme.
 *
 ************************************************************************/

#include "MblkHyperbolicPatchStrategy.h"

#include "SAMRAI/tbox/Utilities.h"

using namespace std;
using namespace SAMRAI;

MblkHyperbolicPatchStrategy::MblkHyperbolicPatchStrategy(
   const tbox::Dimension& dim):
   xfer::RefinePatchStrategy(dim),
   xfer::CoarsenPatchStrategy(dim),
   d_dim(dim)
{
   d_data_context.reset();
}

MblkHyperbolicPatchStrategy::~MblkHyperbolicPatchStrategy()
{
}

/*
 *************************************************************************
 *
 * Default virtual function implementations.
 *
 *************************************************************************
 */

void MblkHyperbolicPatchStrategy::tagGradientDetectorCells(
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
   TBOX_WARNING("MblkHyperbolicPatchStrategy::tagGradientDetectorCells()"
      << "\nNo class supplies a concrete implementation for "
      << "\nthis method.  The default abstract method (which "
      << "\ndoes no cell tagging) is executed" << endl);
}

void MblkHyperbolicPatchStrategy::tagRichardsonExtrapolationCells(
   hier::Patch& patch,
   const int error_level_number,
   const boost::shared_ptr<hier::VariableContext> coarsened_fine,
   const boost::shared_ptr<hier::VariableContext> advanced_coarse,
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
   TBOX_WARNING(
      "MblkHyperbolicPatchStrategy::tagRichardsonExtrapolationCells()"
      << "\nNo class supplies a concrete implementation for "
      << "\nthis method.  The default abstract method (which "
      << "\ndoes no cell tagging) is executed" << endl);
}
