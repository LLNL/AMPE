// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
#include "EllipticFACOps.h"
#include "CellPoissonHypreSolver.h"

#include "SAMRAI/hier/BoundaryBoxUtils.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/pdat/CellDoubleConstantRefine.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/pdat/OutersideData.h"
#include "SAMRAI/pdat/OutersideVariable.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/tbox/MathUtilities.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/PatchLevelFullFillPattern.h"
#include "SAMRAI/tbox/Database.h"

#include "fc_samrai_mangle.h"

#include <cassert>


extern "C" {

void EFO_COMPFLUXVARDC2D(double *xflux, double *yflux, const int *fluxgi,
                         const int *fluxgj, const double *xdiff_coef,
                         const double *ydiff_coef, const int *dcgi,
                         const int *dcgj, const double *soln, const int *solngi,
                         const int *solngj, const int *ifirst, const int *ilast,
                         const int *jfirst, const int *jlast, const double *dx);
void SAMRAI_F77_FUNC(compfluxcondc2d, COMPFLUXCONDC2D)(
    double *xflux, double *yflux, const int *fluxgi, const int *fluxgj,
    const double &diff_coef, const double *soln, const int *solngi,
    const int *solngj, const int *ifirst, const int *ilast, const int *jfirst,
    const int *jlast, const double *dx);
void SAMRAI_F77_FUNC(compfluxcondc3d, COMPFLUXCONDC3D)(
    double *xflux, double *yflux, double *zflux, const int *fluxgi,
    const int *fluxgj, const int *fluxgk, const double &diff_coef,
    const double *soln, const int *solngi, const int *solngj, const int *solngk,
    const int *ifirst, const int *ilast, const int *jfirst, const int *jlast,
    const int *kfirst, const int *klast, const double *dx);
void EFO_RBGSWITHFLUXMAXVARDCVARSF2D(
    const double *xflux, const double *yflux, const int *fluxgi,
    const int *fluxgj, const double *xdiff_coef, const double *ydiff_coef,
    const int *dcgi, const int *dcgj, const double *rhs, const int *rhsgi,
    const int *rhsgj, const double *scalar_field, const int *scalar_field_gi,
    const int *scalar_field_gj, const double *m, const int *m_gi,
    const int *m_gj, double *soln, const int *solngi, const int *solngj,
    const int *ifirst, const int *ilast, const int *jfirst, const int *jlast,
    const double *dx, const int *offset, const double *maxres);
void EFO_RBGSWITHFLUXMAXCONDCVARSF2D(
    const double *xflux, const double *yflux, const int *fluxgi,
    const int *fluxgj, const double &dc, const double *rhs, const int *rhsgi,
    const int *rhsgj, const double *scalar_field, const int *scalar_field_gi,
    const int *scalar_field_gj, const double *m, const int *m_gi,
    const int *m_gj, double *soln, const int *solngi, const int *solngj,
    const int *ifirst, const int *ilast, const int *jfirst, const int *jlast,
    const double *dx, const int *offset, const double *maxres);
void EFO_RBGSWITHFLUXMAXVARDCCONSF2D(
    const double *xflux, const double *yflux, const int *fluxgi,
    const int *fluxgj, const double *xdiff_coef, const double *ydiff_coef,
    const int *dcgi, const int *dcgj, const double *rhs, const int *rhsgi,
    const int *rhsgj, const double &scalar_field, const double *m,
    const int *m_gi, const int *m_gj, double *soln, const int *solngi,
    const int *solngj, const int *ifirst, const int *ilast, const int *jfirst,
    const int *jlast, const double *dx, const int *offset,
    const double *maxres);
void EFO_RBGSWITHFLUXMAXCONDCCONSF2D(
    const double *xflux, const double *yflux, const int *fluxgi,
    const int *fluxgj, const double &dc, const double *rhs, const int *rhsgi,
    const int *rhsgj, const double &scalar_field, const double *m,
    const int *m_gi, const int *m_gj, double *soln, const int *solngi,
    const int *solngj, const int *ifirst, const int *ilast, const int *jfirst,
    const int *jlast, const double *dx, const int *offset,
    const double *maxres);
void EFO_COMPRESVARSCA2D(const double *xflux, const double *yflux,
                         const int *fluxgi, const int *fluxgj,
                         const double *rhs, const int *rhsgi, const int *rhsgj,
                         double *residual, const int *residualgi,
                         const int *residualgj, const double *scalar_field,
                         const int *scalar_field_gi, const int *scalar_field_gj,
                         const double *m, const int *mgi, const int *mgj,
                         const double *soln, const int *solngi,
                         const int *solngj, const int *ifirst, const int *ilast,
                         const int *jfirst, const int *jlast, const double *dx);
void EFO_COMPRESCONSCA2D(const double *xflux, const double *yflux,
                         const int *fluxgi, const int *fluxgj,
                         const double *rhs, const int *rhsgi, const int *rhsgj,
                         double *residual, const int *residualgi,
                         const int *residualgj, const double &scalar_field,
                         const double *m, const int *mgi, const int *mgj,
                         const double *soln, const int *solngi,
                         const int *solngj, const int *ifirst, const int *ilast,
                         const int *jfirst, const int *jlast, const double *dx);
void ACCUMOPVARSCA2D(const double *xflux, const double *yflux,
                     const int *fluxgi, const int *fluxgj, const double *accum,
                     const int *accumgi, const int *accumgj,
                     const double *scalar_field, const int *scalar_field_gi,
                     const int *scalar_field_gj, const double *m,
                     const int *mgi, const int *mgj, const double *soln,
                     const int *solngi, const int *solngj, const int *ifirst,
                     const int *ilast, const int *jfirst, const int *jlast,
                     const double *dx);
void ACCUMOPCONSCA2D(const double *xflux, const double *yflux,
                     const int *fluxgi, const int *fluxgj, const double *accum,
                     const int *accumgi, const int *accumgj,
                     const double &scalar_field, const double *m,
                     const int *mgi, const int *mgj, const double *soln,
                     const int *solngi, const int *solngj, const int *ifirst,
                     const int *ilast, const int *jfirst, const int *jlast,
                     const double *dx);
void EFO_EWINGFIXFLUXVARDC2D(
    const double *xflux, const double *yflux, const int *fluxgi,
    const int *fluxgj, const double *xdiff_coef, const double *ydiff_coef,
    const int *dcgi, const int *dcgj, const double *soln, const int *solngi,
    const int *solngj, const int *ifirst, const int *ilast, const int *jfirst,
    const int *jlast, const int *location_index, const int *ratio_to_coarser,
    const int *blower, const int *bupper, const double *dx);
void EFO_EWINGFIXFLUXCONDC2D(const double *xflux, const double *yflux,
                             const int *fluxgi, const int *fluxgj,
                             const double &diff_coef, const double *soln,
                             const int *solngi, const int *solngj,
                             const int *ifirst, const int *ilast,
                             const int *jfirst, const int *jlast,
                             const int *location_index,
                             const int *ratio_to_coarser, const int *blower,
                             const int *bupper, const double *dx);

void EFO_COMPFLUXVARDC3D(double *xflux, double *yflux, double *zflux,
                         const int *fluxgi, const int *fluxgj,
                         const int *fluxgk, const double *xdiff_coef,
                         const double *ydiff_coef, const double *zdiff_coef,
                         const int *dcgi, const int *dcgj, const int *dcgk,
                         const double *soln, const int *solngi,
                         const int *solngj, const int *solngk,
                         const int *ifirst, const int *ilast, const int *jfirst,
                         const int *jlast, const int *kfirst, const int *klast,
                         const double *dx);
void EFO_RBGSWITHFLUXMAXVARDCVARSF3D(
    const double *xflux, const double *yflux, const double *zflux,
    const int *fluxgi, const int *fluxgj, const int *fluxgk,
    const double *xdiff_coef, const double *ydiff_coef,
    const double *zdiff_coef, const int *dcgi, const int *dcgj, const int *dcgk,
    const double *rhs, const int *rhsgi, const int *rhsgj, const int *rhsgk,
    const double *scalar_field, const int *scalar_field_gi,
    const int *scalar_field_gj, const int *scalar_field_gk, const double *m,
    const int *m_gi, const int *m_gj, const int *m_gk, double *soln,
    const int *solngi, const int *solngj, const int *solngk, const int *ifirst,
    const int *ilast, const int *jfirst, const int *jlast, const int *kfirst,
    const int *klast, const double *dx, const int *offset,
    const double *maxres);
void EFO_RBGSWITHFLUXMAXCONDCVARSF3D(
    const double *xflux, const double *yflux, const double *zflux,
    const int *fluxgi, const int *fluxgj, const int *fluxgk, const double &dc,
    const double *rhs, const int *rhsgi, const int *rhsgj, const int *rhsgk,
    const double *scalar_field, const int *scalar_field_gi,
    const int *scalar_field_gj, const int *scalar_field_gk, const double *m,
    const int *m_gi, const int *m_gj, const int *m_gk, double *soln,
    const int *solngi, const int *solngj, const int *solngk, const int *ifirst,
    const int *ilast, const int *jfirst, const int *jlast, const int *kfirst,
    const int *klast, const double *dx, const int *offset,
    const double *maxres);
void EFO_RBGSWITHFLUXMAXVARDCCONSF3D(
    const double *xflux, const double *yflux, const double *zflux,
    const int *fluxgi, const int *fluxgj, const int *fluxgk,
    const double *xdiff_coef, const double *ydiff_coef,
    const double *zdiff_coef, const int *dcgi, const int *dcgj, const int *dcgk,
    const double *rhs, const int *rhsgi, const int *rhsgj, const int *rhsgk,
    const double &scalar_field, const double *m, const int *m_gi,
    const int *m_gj, const int *m_gk, double *soln, const int *solngi,
    const int *solngj, const int *solngk, const int *ifirst, const int *ilast,
    const int *jfirst, const int *jlast, const int *kfirst, const int *klast,
    const double *dx, const int *offset, const double *maxres);
void EFO_RBGSWITHFLUXMAXCONDCCONSF3D(
    const double *xflux, const double *yflux, const double *zflux,
    const int *fluxgi, const int *fluxgj, const int *fluxgk, const double &dc,
    const double *rhs, const int *rhsgi, const int *rhsgj, const int *rhsgk,
    const double &scalar_field, const double *m, const int *m_gi,
    const int *m_gj, const int *m_gk, double *soln, const int *solngi,
    const int *solngj, const int *solngk, const int *ifirst, const int *ilast,
    const int *jfirst, const int *jlast, const int *kfirst, const int *klast,
    const double *dx, const int *offset, const double *maxres);
void EFO_COMPRESVARSCA3D(
    const double *xflux, const double *yflux, const double *zflux,
    const int *fluxgi, const int *fluxgj, const int *fluxgk, const double *rhs,
    const int *rhsgi, const int *rhsgj, const int *rhsgk, double *residual,
    const int *residualgi, const int *residualgj, const int *residualgk,
    const double *scalar_field, const int *scalar_field_gi,
    const int *scalar_field_gj, const int *scalar_field_gk, const double *m,
    const int *m_gi, const int *m_gj, const int *m_gk, const double *soln,
    const int *solngi, const int *solngj, const int *solngk, const int *ifirst,
    const int *ilast, const int *jfirst, const int *jlast, const int *kfirst,
    const int *klast, const double *dx);
void EFO_COMPRESCONSCA3D(
    const double *xflux, const double *yflux, const double *zflux,
    const int *fluxgi, const int *fluxgj, const int *fluxgk, const double *rhs,
    const int *rhsgi, const int *rhsgj, const int *rhsgk, double *residual,
    const int *residualgi, const int *residualgj, const int *residualgk,
    const double &scalar_field, const double *m, const int *m_gi,
    const int *m_gj, const int *m_gk, const double *soln, const int *solngi,
    const int *solngj, const int *solngk, const int *ifirst, const int *ilast,
    const int *jfirst, const int *jlast, const int *kfirst, const int *klast,
    const double *dx);
void ACCUMOPVARSCA3D(const double *xflux, const double *yflux,
                     const double *zflux, const int *fluxgi, const int *fluxgj,
                     const int *fluxgk, const double *accum, const int *accumgi,
                     const int *accumgj, const int *accumgk,
                     const double *scalar_field, const int *scalar_field_gi,
                     const int *scalar_field_gj, const int *scalar_field_gk,
                     const double *m, const int *m_gi, const int *m_gj,
                     const int *m_gk, const double *soln, const int *solngi,
                     const int *solngj, const int *solngk, const int *ifirst,
                     const int *ilast, const int *jfirst, const int *jlast,
                     const int *kfirst, const int *klast, const double *dx);
void ACCUMOPCONSCA3D(const double *xflux, const double *yflux,
                     const double *zflux, const int *fluxgi, const int *fluxgj,
                     const int *fluxgk, const double *accum, const int *accumgi,
                     const int *accumgj, const int *accumgk,
                     const double &scalar_field, const double *m,
                     const int *m_gi, const int *m_gj, const int *m_gk,
                     const double *soln, const int *solngi, const int *solngj,
                     const int *solngk, const int *ifirst, const int *ilast,
                     const int *jfirst, const int *jlast, const int *kfirst,
                     const int *klast, const double *dx);
void EFO_EWINGFIXFLUXVARDC3D(
    const double *xflux, const double *yflux, const double *zflux,
    const int *fluxgi, const int *fluxgj, const int *fluxgk,
    const double *xdiff_coef, const double *ydiff_coef,
    const double *zdiff_coef, const int *dcgi, const int *dcgj, const int *dcgk,
    const double *soln, const int *solngi, const int *solngj, const int *solngk,
    const int *ifirst, const int *ilast, const int *jfirst, const int *jlast,
    const int *kfirst, const int *klast, const int *location_index,
    const int *ratio_to_coarser, const int *blower, const int *bupper,
    const double *dx);
void EFO_EWINGFIXFLUXCONDC3D(
    const double *xflux, const double *yflux, const double *zflux,
    const int *fluxgi, const int *fluxgj, const int *fluxgk,
    const double &diff_coef, const double *soln, const int *solngi,
    const int *solngj, const int *solngk, const int *ifirst, const int *ilast,
    const int *jfirst, const int *jlast, const int *kfirst, const int *klast,
    const int *location_index, const int *ratio_to_coarser, const int *blower,
    const int *bupper, const double *dx);
}

/*
********************************************************************
* Constructor.                                                     *
********************************************************************
*/
EllipticFACOps::EllipticFACOps(const tbox::Dimension &dim,
                               const std::string &object_name,
                               const std::shared_ptr<tbox::Database> &database,
                               const int depth)
    : d_dim(dim),
      d_object_name(object_name),
      d_ln_min(-1),
      d_ln_max(-1),
      d_cf_boundary(),
      d_smoothing_choice("redblack"),
      d_coarse_solver_choice("hypre"),
      d_cf_discretization("Ewing"),
      d_prolongation_method("CONSTANT_REFINE"),
      d_coarse_solver_tolerance(1.e-2),
      d_coarse_solver_max_iterations(20),
      d_residual_tolerance_during_smoothing(-1.0),
      d_flux_id(-1),
      d_context(hier::VariableDatabase::getDatabase()->getContext(
          object_name + "::PRIVATE_CONTEXT")),
      d_cell_scratch_id(-1),
      d_flux_scratch_id(-1),
      d_oflux_scratch_id(-1),
      d_bc_helper(tbox::Dimension(NDIM), d_object_name + "::bc helper"),
      d_enable_logging(false),
      d_preconditioner(NULL),
      d_hopscell(),
      d_hopsside(),
      d_physical_bc_coef(NULL),
      d_depth(depth)
{
   if (NDIM == 1) {
      TBOX_ERROR(d_object_name << ": 1D not implemented yet.\n");
   }

   for (int i = 0; i < depth; i++) {
      CellPoissonHypreSolver *hypre_solver =
          new CellPoissonHypreSolver(object_name + "::hypre_solver" +
                                         std::to_string(i),
                                     database && database->isDatabase("hypre_"
                                                                      "solver")
                                         ? database->getDatabase("hypre_solver")
                                         : std::shared_ptr<tbox::Database>());
      d_hypre_solver.push_back(hypre_solver);

      PoissonSpecifications *poisson_spec =
          new PoissonSpecifications(object_name + "::Poisson specs");
      d_poisson_spec.push_back(*poisson_spec);

      d_C_is_set.push_back(false);
      d_D_is_set.push_back(false);
   }
   d_M_is_set = false;

   t_restrict_solution = tbox::TimerManager::getManager()->getTimer(
       "AMPE::EllipticFACOps::restrictSolution()");
   t_restrict_residual = tbox::TimerManager::getManager()->getTimer(
       "AMPE::EllipticFACOps::restrictResidual()");
   t_prolong = tbox::TimerManager::getManager()->getTimer(
       "AMPE::EllipticFACOps::prolongErrorAndCorrect()");
   t_smooth_error = tbox::TimerManager::getManager()->getTimer(
       "AMPE::EllipticFACOps::smoothError()");
   t_solve_coarsest = tbox::TimerManager::getManager()->getTimer(
       "AMPE::EllipticFACOps::solveCoarsestLevel()");
   t_compute_composite_residual = tbox::TimerManager::getManager()->getTimer(
       "AMPE::EllipticFACOps::computeCompositeResidualOnLevel()");
   t_accumulate_operator = tbox::TimerManager::getManager()->getTimer(
       "AMPE::EllipticFACOps::accumulateOperatorOnLevel()");
   t_compute_residual_norm = tbox::TimerManager::getManager()->getTimer(
       "AMPE::EllipticFACOps::computeResidualNorm()");
   t_compute_rhs = tbox::TimerManager::getManager()->getTimer(
       "AMPE::EllipticFACOps::evaluateRHS()");
   t_finalizecoeffs = tbox::TimerManager::getManager()->getTimer(
       "AMPE::EllipticFACOps::finalizeCoefficients()");

   TBOX_ASSERT(!d_cell_scratch_var);

   d_cell_scratch_var.reset(new pdat::CellVariable<double>(
       tbox::Dimension(NDIM), object_name + "::private_cell_scratch", d_depth));
   d_flux_scratch_var.reset(new pdat::SideVariable<double>(
       tbox::Dimension(NDIM), object_name + "::private_flux_scratch", d_depth));
   d_oflux_scratch_var.reset(
       new pdat::OutersideVariable<double>(tbox::Dimension(NDIM),
                                           object_name + "::private_oflux_"
                                                         "scratch",
                                           d_depth));

   d_m_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                object_name + "EllipticFACOps::"
                                                              "private_m",
                                                1));
   for (int i = 0; i < depth; i++) {
      std::shared_ptr<pdat::SideVariable<double> > d_var;
      d_var.reset(new pdat::SideVariable<double>(
          tbox::Dimension(NDIM),
          object_name + "EllipticFACOps::privateD" + std::to_string(i), 1));
      d_d_var.push_back(d_var);

      std::shared_ptr<pdat::CellVariable<double> > c_var;
      c_var.reset(new pdat::CellVariable<double>(
          tbox::Dimension(NDIM),
          object_name + "EllipticFACOps::privateC" + std::to_string(i), 1));
      d_c_var.push_back(c_var);
   }

   /*
    * Some variables initialized by default are overriden by input.
    */
   getFromInput(database);

   hier::VariableDatabase *vdb = hier::VariableDatabase::getDatabase();
   d_cell_scratch_id =
       vdb->registerVariableAndContext(d_cell_scratch_var, d_context,
                                       hier::IntVector(tbox::Dimension(NDIM),
                                                       1));
   d_flux_scratch_id =
       vdb->registerVariableAndContext(d_flux_scratch_var, d_context,
                                       hier::IntVector(tbox::Dimension(NDIM),
                                                       0));
   d_oflux_scratch_id =
       vdb->registerVariableAndContext(d_oflux_scratch_var, d_context,
                                       hier::IntVector(tbox::Dimension(NDIM),
                                                       0));
   d_m_id = vdb->registerVariableAndContext(
       d_m_var, d_context, hier::IntVector(tbox::Dimension(NDIM), 1));

   for (int i = 0; i < depth; i++) {
      assert(i < static_cast<int>(d_d_var.size()));
      d_d_id.push_back(vdb->registerVariableAndContext(
          d_d_var[i], d_context, hier::IntVector(tbox::Dimension(NDIM), 0)));

      d_c_id.push_back(vdb->registerVariableAndContext(
          d_c_var[i], d_context, hier::IntVector(tbox::Dimension(NDIM), 0)));
   }

   /*
    * Check input validity and correctness.
    */
   checkInputPatchDataIndices();
}

void EllipticFACOps::getFromInput(
    const std::shared_ptr<tbox::Database> &input_db)
{
   if (input_db) {
      d_coarse_solver_choice =
          input_db->getStringWithDefault("coarse_solver_choice",
                                         d_coarse_solver_choice);
      if (!(d_coarse_solver_choice == "hypre" ||
            d_coarse_solver_choice == "redblack")) {
         INPUT_VALUE_ERROR("coarse_solver_choice");
      }

      d_coarse_solver_tolerance =
          input_db->getDoubleWithDefault("coarse_solver_tolerance",
                                         d_coarse_solver_tolerance);

      d_coarse_solver_max_iterations =
          input_db->getIntegerWithDefault("coarse_solver_max_iterations",
                                          d_coarse_solver_max_iterations);
      if (!(d_coarse_solver_max_iterations >= 1)) {
         INPUT_RANGE_ERROR("coarse_solver_max_iterations");
      }

      d_cf_discretization =
          input_db->getStringWithDefault("cf_discretization", "Ewing");
      if (!(d_cf_discretization == "Ewing" ||
            d_cf_discretization == "CONSTANT_REFINE" ||
            d_cf_discretization == "LINEAR_REFINE" ||
            d_cf_discretization == "CONSERVATIVE_LINEAR_REFINE")) {
         INPUT_VALUE_ERROR("cf_discretization");
      }

      d_prolongation_method =
          input_db->getStringWithDefault("prolongation_method",
                                         "CONSTANT_REFINE");
      if (!(d_prolongation_method == "CONSTANT_REFINE" ||
            d_prolongation_method == "LINEAR_REFINE" ||
            d_prolongation_method == "CONSERVATIVE_LINEAR_REFINE")) {
         INPUT_VALUE_ERROR("prolongation_method");
      }

      d_enable_logging = input_db->getBoolWithDefault("enable_logging", false);
   }
}


/*
************************************************************************
* FACOperatorStrategy virtual initializeOperatorState function.  *
*                                                                      *
* Set internal variables to correspond to the solution passed in.      *
* Look up transfer operators.                                          *
************************************************************************
*/
void EllipticFACOps::initializeOperatorState(
    const solv::SAMRAIVectorReal<double> &solution,
    const solv::SAMRAIVectorReal<double> &rhs)
{
   deallocateOperatorState();
   hier::VariableDatabase *vdb = hier::VariableDatabase::getDatabase();

   d_hierarchy = solution.getPatchHierarchy();
   d_ln_min = solution.getCoarsestLevelNumber();
   d_ln_max = solution.getFinestLevelNumber();
   d_hopscell.reset(new math::HierarchyCellDataOpsReal<double>(d_hierarchy,
                                                               d_ln_min,
                                                               d_ln_max));
   d_hopsside.reset(new math::HierarchySideDataOpsReal<double>(d_hierarchy,
                                                               d_ln_min,
                                                               d_ln_max));

#ifdef DEBUG_CHECK_ASSERTIONS
   if (d_physical_bc_coef == NULL) {
      /*
       * It's an error not to have bc object set.
       * Note that the bc object cannot be passed in through
       * the argument because the interface is inherited.
       */
      TBOX_ERROR(d_object_name << ": No physical bc object in\n"
                               << "EllipticFACOps::initializeOperatorState\n"
                               << "You must use "
                               << "EllipticFACOps::setPhysicalBcCoefObject\n"
                               << "to set one before calling "
                                  "initializeOperatorState\n");
   }

   /*
    * Make sure that solution and rhs data
    *   are of correct type
    *   are allocated
    *   has sufficient ghost width
    */
   std::shared_ptr<hier::Variable> var;
   {
      vdb->mapIndexToVariable(rhs.getComponentDescriptorIndex(0), var);
      TBOX_ASSERT(var);
      std::shared_ptr<pdat::CellVariable<double> > cell_var(
          SAMRAI_SHARED_PTR_CAST<pdat::CellVariable<double>, hier::Variable>(
              var));
      TBOX_ASSERT(cell_var);
   }
   {
      vdb->mapIndexToVariable(solution.getComponentDescriptorIndex(0), var);
      TBOX_ASSERT(var);
      std::shared_ptr<pdat::CellVariable<double> > cell_var(
          SAMRAI_SHARED_PTR_CAST<pdat::CellVariable<double>, hier::Variable>(
              var));
      TBOX_ASSERT(cell_var);
   }
   for (int ln = d_ln_min; ln <= d_ln_max; ++ln) {
      std::shared_ptr<hier::PatchLevel> level_ptr(
          d_hierarchy->getPatchLevel(ln));
      hier::PatchLevel &level = *level_ptr;
      for (hier::PatchLevel::iterator pi(level.begin()); pi != level.end();
           ++pi) {
         hier::Patch &patch = **pi;
         std::shared_ptr<hier::PatchData> fd(
             patch.getPatchData(rhs.getComponentDescriptorIndex(0)));
         if (fd) {
            /*
             * Some data checks can only be done if the data already exists.
             */
            std::shared_ptr<pdat::CellData<double> > cd(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    fd));
            TBOX_ASSERT(cd);
         }
         std::shared_ptr<hier::PatchData> ud(
             patch.getPatchData(solution.getComponentDescriptorIndex(0)));
         if (ud) {
            /*
             * Some data checks can only be done if the data already exists.
             */
            std::shared_ptr<pdat::CellData<double> > cd(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    ud));
            TBOX_ASSERT(cd);
            if (cd->getGhostCellWidth() <
                hier::IntVector::getOne(tbox::Dimension(NDIM))) {
               TBOX_ERROR(d_object_name << ": Solution data has insufficient "
                                           "ghost width\n");
            }
         }
      }
   }

   /*
    * Solution and rhs must have some similar properties.
    */
   if (rhs.getPatchHierarchy() != d_hierarchy ||
       rhs.getCoarsestLevelNumber() != d_ln_min ||
       rhs.getFinestLevelNumber() != d_ln_max) {
      TBOX_ERROR(d_object_name << ": solution and rhs do not have\n"
                               << "the same set of patch levels.\n");
   }

#endif

   /*
    * Initialize the coarse-fine boundary description for the
    * hierarchy.
    */
   d_cf_boundary.resize(d_hierarchy->getNumberOfLevels());

   hier::IntVector max_gcw(tbox::Dimension(NDIM), 1);
   for (int ln = d_ln_min; ln <= d_ln_max; ++ln) {
      d_cf_boundary[ln].reset(
          new hier::CoarseFineBoundary(*d_hierarchy, ln, max_gcw));
   }

   /*
    * Get the transfer operators.
    * Flux coarsening is conservative.
    * Cell (solution, error, etc) coarsening is conservative.
    * Cell refinement from same level is constant refinement.
    * Cell refinement from coarser level is chosen by the
    *   choice of coarse-fine discretization, d_cf_discretization,
    *   which should be set to either "Ewing" or one of the
    *   acceptable std::strings for looking up the refine operator.
    */
   std::shared_ptr<geom::CartesianGridGeometry> geometry(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry,
                              hier::BaseGridGeometry>(
           d_hierarchy->getGridGeometry()));
   TBOX_ASSERT(geometry);

   std::shared_ptr<hier::Variable> variable;

   vdb->mapIndexToVariable(d_cell_scratch_id, variable);
   d_prolongation_refine_operator =
       geometry->lookupRefineOperator(variable, d_prolongation_method);

   vdb->mapIndexToVariable(d_cell_scratch_id, variable);
   d_urestriction_coarsen_operator = d_rrestriction_coarsen_operator =
       geometry->lookupCoarsenOperator(variable, "CONSERVATIVE_COARSEN");

   vdb->mapIndexToVariable(d_oflux_scratch_id, variable);
   d_flux_coarsen_operator =
       geometry->lookupCoarsenOperator(variable, "CONSERVATIVE_COARSEN");

   vdb->mapIndexToVariable(d_cell_scratch_id, variable);
   d_ghostfill_refine_operator =
       geometry->lookupRefineOperator(variable, d_cf_discretization == "Ewing"
                                                    ? "CONSTANT_REFINE"
                                                    : d_cf_discretization);

   vdb->mapIndexToVariable(d_cell_scratch_id, variable);
   d_ghostfill_nocoarse_refine_operator =
       geometry->lookupRefineOperator(variable, "CONSTANT_REFINE");

#ifdef DEBUG_CHECK_ASSERTIONS
   if (!d_prolongation_refine_operator) {
      TBOX_ERROR(d_object_name << ": Cannot find prolongation refine operator");
   }
   if (!d_urestriction_coarsen_operator) {
      TBOX_ERROR(d_object_name << ": Cannot find restriction coarsening "
                                  "operator");
   }
   if (!d_rrestriction_coarsen_operator) {
      TBOX_ERROR(d_object_name << ": Cannot find restriction coarsening "
                                  "operator");
   }
   if (!d_flux_coarsen_operator) {
      TBOX_ERROR(d_object_name << ": Cannot find flux coarsening operator");
   }
   if (!d_ghostfill_refine_operator) {
      TBOX_ERROR(d_object_name << ": Cannot find ghost filling refinement "
                                  "operator");
   }
   if (!d_ghostfill_nocoarse_refine_operator) {
      TBOX_ERROR(d_object_name << ": Cannot find ghost filling refinement "
                                  "operator");
   }
#endif

   for (int ln = d_ln_min + 1; ln <= d_ln_max; ++ln) {
      d_hierarchy->getPatchLevel(ln)->allocatePatchData(d_oflux_scratch_id);
   }

   for (int ln = d_ln_min; ln <= d_ln_max; ++ln) {
      d_hierarchy->getPatchLevel(ln)->allocatePatchData(d_m_id);
   }

   for (int ln = d_ln_min; ln <= d_ln_max; ++ln) {
      for (int i = 0; i < d_depth; i++)
         d_hierarchy->getPatchLevel(ln)->allocatePatchData(d_c_id[i]);
   }

   for (int ln = d_ln_min; ln <= d_ln_max; ++ln) {
      for (int i = 0; i < d_depth; i++)
         d_hierarchy->getPatchLevel(ln)->allocatePatchData(d_d_id[i]);
   }

   /*
     Make space for saving communication schedules.
     There is no need to delete the old schedules first
     because we have deallocated the solver state above.
   */
   d_prolongation_refine_schedules.resize(d_ln_max + 1);
   d_ghostfill_refine_schedules.resize(d_ln_max + 1);
   d_ghostfill_nocoarse_refine_schedules.resize(d_ln_max + 1);
   d_urestriction_coarsen_schedules.resize(d_ln_max + 1);
   d_rrestriction_coarsen_schedules.resize(d_ln_max + 1);
   d_flux_coarsen_schedules.resize(d_ln_max + 1);

   d_prolongation_refine_algorithm.reset(new xfer::RefineAlgorithm());
   d_urestriction_coarsen_algorithm.reset(
       new xfer::CoarsenAlgorithm(tbox::Dimension(NDIM)));
   d_rrestriction_coarsen_algorithm.reset(
       new xfer::CoarsenAlgorithm(tbox::Dimension(NDIM)));
   d_flux_coarsen_algorithm.reset(
       new xfer::CoarsenAlgorithm(tbox::Dimension(NDIM)));
   d_ghostfill_refine_algorithm.reset(new xfer::RefineAlgorithm());
   d_ghostfill_nocoarse_refine_algorithm.reset(new xfer::RefineAlgorithm());

   d_prolongation_refine_algorithm->registerRefine(
       d_cell_scratch_id, solution.getComponentDescriptorIndex(0),
       d_cell_scratch_id, d_prolongation_refine_operator);
   d_urestriction_coarsen_algorithm->registerCoarsen(
       solution.getComponentDescriptorIndex(0),
       solution.getComponentDescriptorIndex(0),
       d_urestriction_coarsen_operator);
   d_rrestriction_coarsen_algorithm->registerCoarsen(
       rhs.getComponentDescriptorIndex(0), rhs.getComponentDescriptorIndex(0),
       d_rrestriction_coarsen_operator);
   d_ghostfill_refine_algorithm->registerRefine(
       solution.getComponentDescriptorIndex(0),
       solution.getComponentDescriptorIndex(0),
       solution.getComponentDescriptorIndex(0), d_ghostfill_refine_operator);
   d_flux_coarsen_algorithm->registerCoarsen(
       ((d_flux_id != -1) ? d_flux_id : d_flux_scratch_id), d_oflux_scratch_id,
       d_flux_coarsen_operator);
   d_ghostfill_nocoarse_refine_algorithm->registerRefine(
       solution.getComponentDescriptorIndex(0),
       solution.getComponentDescriptorIndex(0),
       solution.getComponentDescriptorIndex(0),
       d_ghostfill_nocoarse_refine_operator);

   for (int dest_ln = d_ln_min + 1; dest_ln <= d_ln_max; ++dest_ln) {

      std::shared_ptr<xfer::PatchLevelFullFillPattern> fill_pattern(
          std::make_shared<xfer::PatchLevelFullFillPattern>());
      d_prolongation_refine_schedules[dest_ln] =
          d_prolongation_refine_algorithm->createSchedule(
              fill_pattern, d_hierarchy->getPatchLevel(dest_ln),
              std::shared_ptr<hier::PatchLevel>(), dest_ln - 1, d_hierarchy,
              &d_bc_helper);
      if (!d_prolongation_refine_schedules[dest_ln]) {
         TBOX_ERROR(d_object_name << ": Cannot create a refine schedule for "
                                     "prolongation!\n");
      }
      d_ghostfill_refine_schedules[dest_ln] =
          d_ghostfill_refine_algorithm->createSchedule(
              d_hierarchy->getPatchLevel(dest_ln), dest_ln - 1, d_hierarchy,
              &d_bc_helper);
      if (!d_ghostfill_refine_schedules[dest_ln]) {
         TBOX_ERROR(d_object_name << ": Cannot create a refine schedule for "
                                     "ghost filling!\n");
      }
      d_ghostfill_nocoarse_refine_schedules[dest_ln] =
          d_ghostfill_nocoarse_refine_algorithm->createSchedule(
              d_hierarchy->getPatchLevel(dest_ln), &d_bc_helper);
      if (!d_ghostfill_nocoarse_refine_schedules[dest_ln]) {
         TBOX_ERROR(d_object_name << ": Cannot create a refine schedule for "
                                     "ghost filling on bottom level!\n");
      }
   }
   for (int dest_ln = d_ln_min; dest_ln < d_ln_max; ++dest_ln) {
      d_urestriction_coarsen_schedules[dest_ln] =
          d_urestriction_coarsen_algorithm->createSchedule(
              d_hierarchy->getPatchLevel(dest_ln),
              d_hierarchy->getPatchLevel(dest_ln + 1));
      if (!d_urestriction_coarsen_schedules[dest_ln]) {
         TBOX_ERROR(d_object_name << ": Cannot create a coarsen schedule for U "
                                     "restriction!\n");
      }
      d_rrestriction_coarsen_schedules[dest_ln] =
          d_rrestriction_coarsen_algorithm->createSchedule(
              d_hierarchy->getPatchLevel(dest_ln),
              d_hierarchy->getPatchLevel(dest_ln + 1));
      if (!d_rrestriction_coarsen_schedules[dest_ln]) {
         TBOX_ERROR(d_object_name << ": Cannot create a coarsen schedule for R "
                                     "restriction!\n");
      }
      d_flux_coarsen_schedules[dest_ln] =
          d_flux_coarsen_algorithm->createSchedule(
              d_hierarchy->getPatchLevel(dest_ln),
              d_hierarchy->getPatchLevel(dest_ln + 1));
      if (!d_flux_coarsen_schedules[dest_ln]) {
         TBOX_ERROR(d_object_name << ": Cannot create a coarsen schedule for "
                                     "flux transfer!\n");
      }
   }
   d_ghostfill_nocoarse_refine_schedules[d_ln_min] =
       d_ghostfill_nocoarse_refine_algorithm->createSchedule(
           d_hierarchy->getPatchLevel(d_ln_min), &d_bc_helper);
   if (!d_ghostfill_nocoarse_refine_schedules[d_ln_min]) {
      TBOX_ERROR(d_object_name << ": Cannot create a refine schedule for ghost "
                                  "filling on bottom level!\n");
   }
}

/*
********************************************************************
* FACOperatorStrategy virtual deallocateOperatorState        *
* function.  Deallocate internal hierarchy-dependent data.         *
* State is allocated iff hierarchy is set.                         *
********************************************************************
*/
void EllipticFACOps::deallocateOperatorState()
{
   if (d_hierarchy) {

      int ln;
      for (ln = d_ln_min; ln <= d_ln_max; ++ln) {
         d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_m_id);
      }
      for (ln = d_ln_min; ln <= d_ln_max; ++ln) {
         for (int i = 0; i < d_depth; i++)
            d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_c_id[i]);
      }
      for (ln = d_ln_min; ln <= d_ln_max; ++ln) {
         for (int i = 0; i < d_depth; i++)
            d_hierarchy->getPatchLevel(ln)->deallocatePatchData(d_d_id[i]);
      }
      for (ln = d_ln_min + 1; ln <= d_ln_max; ++ln) {
         d_hierarchy->getPatchLevel(ln)->deallocatePatchData(
             d_oflux_scratch_id);
      }
      d_cf_boundary.resize(0);
      std::vector<CellPoissonHypreSolver *>::iterator it(
          d_hypre_solver.begin());
      for (; it != d_hypre_solver.end(); ++it)
         (*it)->deallocateSolverState();
      d_hierarchy.reset();
      d_ln_min = -1;
      d_ln_max = -1;

      d_prolongation_refine_algorithm.reset();
      d_urestriction_coarsen_algorithm.reset();
      d_rrestriction_coarsen_algorithm.reset();
      d_flux_coarsen_algorithm.reset();
      d_ghostfill_refine_algorithm.reset();
      d_ghostfill_nocoarse_refine_algorithm.reset();
   }
   d_C_is_set.clear();
   d_D_is_set.clear();
}


/*
********************************************************************
* Set the object specifying the parameters of the equation *
********************************************************************
*/
void EllipticFACOps::setM(const int m_id)
{
   assert(d_m_id >= 0);
   assert(m_id >= 0);
   assert(d_ln_min >= 0);

   // tbox::pout<<"EllipticFACOps::setM()"<<endl;
   // copy one patch at a time to include ghost values
   SAMRAI::math::PatchCellDataNormOpsReal<double> ops;
   for (int ln = d_ln_min; ln <= d_ln_max; ++ln) {
      std::shared_ptr<hier::PatchLevel> level_ptr(
          d_hierarchy->getPatchLevel(ln));
      hier::PatchLevel &level = *level_ptr;
      for (hier::PatchLevel::iterator pi(level.begin()); pi != level.end();
           ++pi) {
         hier::Patch &patch = **pi;
         std::shared_ptr<pdat::CellData<double> > src(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch.getPatchData(m_id)));
         std::shared_ptr<pdat::CellData<double> > dst(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch.getPatchData(d_m_id)));
         assert(src->getDepth() == dst->getDepth());
#ifdef DEBUG_CHECK_ASSERTIONS
         double l2 = ops.L2Norm(src, src->getBox());
         assert(l2 == l2);
         double l2g = ops.L2Norm(src, src->getGhostBox());
         assert(l2g == l2g);
#endif
         dst->copy(*src);
      }
   }

   for (std::vector<PoissonSpecifications>::iterator it =
            d_poisson_spec.begin();
        it != d_poisson_spec.end(); ++it) {
      it->setMPatchDataId(d_m_id);
   }
   d_M_is_set = true;

   return;
}

void EllipticFACOps::finalizeCoefficients()
{
   tbox::SAMRAI_MPI::getSAMRAIWorld().Barrier();

   t_finalizecoeffs->start();

   size_t n = d_C_is_set.size();
   for (size_t i = 0; i < n; i++)
      if (!d_C_is_set[i] || !d_D_is_set[i] || !d_M_is_set) {
         TBOX_ERROR(d_object_name << ": Operator coefficients have not be "
                                     "set;\n cannot finalize\n");
      }

   if (d_coarse_solver_choice == "hypre") {

      int counter = 0;
      std::vector<CellPoissonHypreSolver *>::iterator it(
          d_hypre_solver.begin());
      for (; it != d_hypre_solver.end(); ++it) {
         (*it)->deallocateSolverState();
         (*it)->initializeSolverState(d_hierarchy, d_ln_min);

         /*
          * Share the boundary condition object with the hypre solver
          * to make sure that boundary condition settings are consistent
          * between the two objects.
          */
         (*it)->setPhysicalBcCoefObject(d_physical_bc_coef);
         (*it)->setMatrixCoefficients(d_poisson_spec[counter]);

         counter++;
      }
   }
   t_finalizecoeffs->stop();
}


/*
********************************************************************
*eFACOperatorStrategy virtual postprocessOneCycle function.  *
l
********************************************************************
*/
void EllipticFACOps::postprocessOneCycle(
    int fac_cycle_num, const solv::SAMRAIVectorReal<double> &current_soln,
    const solv::SAMRAIVectorReal<double> &residual)
{
   NULL_USE(current_soln);
   NULL_USE(residual);

   if (d_enable_logging) {
      if (d_preconditioner) {
         /*
          * Output convergence progress.  This is probably only appropriate
          * if the solver is NOT being used as a preconditioner.
          */
         double avg_factor, final_factor;
         d_preconditioner->getConvergenceFactors(avg_factor, final_factor);
         tbox::plog << "iter=" << std::setw(4) << fac_cycle_num
                    << " resid=" << d_preconditioner->getResidualNorm()
                    << " net conv="
                    << d_preconditioner->getNetConvergenceFactor()
                    << " final conv="
                    << d_preconditioner->getNetConvergenceFactor()
                    << " avg conv="
                    << d_preconditioner->getAvgConvergenceFactor() << std::endl;
      }
   }
   return;
}

/*
********************************************************************
* FACOperatorStrategy virtual restrictSolution function.     *
* After restricting solution, update ghost cells of the affected   *
* level.                                                           *
********************************************************************
*/
void EllipticFACOps::restrictSolution(const solv::SAMRAIVectorReal<double> &s,
                                      solv::SAMRAIVectorReal<double> &d,
                                      int dest_ln)
{
   t_restrict_solution->start();

   xeqScheduleURestriction(d.getComponentDescriptorIndex(0),
                           s.getComponentDescriptorIndex(0), dest_ln);

   d_bc_helper.setHomogeneousBc(false);
   d_bc_helper.setTargetDataId(d.getComponentDescriptorIndex(0));

   if (dest_ln == d_ln_min) {
      xeqScheduleGhostFillNoCoarse(d.getComponentDescriptorIndex(0), dest_ln);
   } else {
      xeqScheduleGhostFill(d.getComponentDescriptorIndex(0), dest_ln);
   }

   t_restrict_solution->stop();
}

/*
********************************************************************
* FACOperatorStrategy virtual restrictresidual function.     *
********************************************************************
*/
void EllipticFACOps::restrictResidual(const solv::SAMRAIVectorReal<double> &s,
                                      solv::SAMRAIVectorReal<double> &d,
                                      int dest_ln)
{
   t_restrict_residual->start();

   xeqScheduleRRestriction(d.getComponentDescriptorIndex(0),
                           s.getComponentDescriptorIndex(0), dest_ln);

   t_restrict_residual->stop();
}

/*
***********************************************************************
* FACOperatorStrategy virtual prolongErrorAndCorrect function.  *
* After the prolongation, we set the physical boundary condition      *
* for the correction, which is zero.  Other ghost cell values,        *
* which are preset to zero, need not be set.                          *
***********************************************************************
*/
void EllipticFACOps::prolongErrorAndCorrect(
    const solv::SAMRAIVectorReal<double> &s, solv::SAMRAIVectorReal<double> &d,
    int dest_ln)
{
   t_prolong->start();

#ifdef DEBUG_CHECK_ASSERTIONS
   if (s.getPatchHierarchy() != d_hierarchy ||
       d.getPatchHierarchy() != d_hierarchy) {
      TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
                                  "internal state hierarchy.");
   }
#endif

   std::shared_ptr<hier::PatchLevel> fine_level =
       d_hierarchy->getPatchLevel(dest_ln);

   /*
    * Data is prolonged into the scratch space corresponding
    * to index d_cell_scratch_id and allocated here.
    */
   fine_level->allocatePatchData(d_cell_scratch_id);

   /*
    * Refine solution into scratch space to fill the fine level
    * interior in the scratch space, then use that refined data
    * to correct the fine level error.
    */
   d_bc_helper.setTargetDataId(d_cell_scratch_id);
   d_bc_helper.setHomogeneousBc(true);
   const int src_index = s.getComponentDescriptorIndex(0);
   xeqScheduleProlongation(d_cell_scratch_id, src_index, d_cell_scratch_id,
                           dest_ln);

   /*
    * Add the refined error in the scratch space
    * to the error currently residing in the destination level.
    */
   math::HierarchyCellDataOpsReal<double> hierarchy_math_ops(d_hierarchy,
                                                             dest_ln, dest_ln);
   const int dst_index = d.getComponentDescriptorIndex(0);
   hierarchy_math_ops.add(dst_index, dst_index, d_cell_scratch_id);

   fine_level->deallocatePatchData(d_cell_scratch_id);

   t_prolong->stop();
}


/*
 ********************************************************************
 ********************************************************************
 */
void EllipticFACOps::smoothError(solv::SAMRAIVectorReal<double> &data,
                                 const solv::SAMRAIVectorReal<double> &residual,
                                 int ln, int num_sweeps)
{
   t_smooth_error->start();

   checkInputPatchDataIndices();
   smoothErrorByRedBlack(data, residual, ln, num_sweeps,
                         d_residual_tolerance_during_smoothing);

   t_smooth_error->stop();
}


/*
********************************************************************
* Workhorse function to smooth error using red-black               *
* Gauss-Seidel iterations.                                         *
********************************************************************
*/
void EllipticFACOps::smoothErrorByRedBlack(
    solv::SAMRAIVectorReal<double> &data,
    const solv::SAMRAIVectorReal<double> &residual, int ln, int num_sweeps,
    double residual_tolerance)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   if (data.getPatchHierarchy() != d_hierarchy ||
       residual.getPatchHierarchy() != d_hierarchy) {
      TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
                                  "internal hierarchy.");
   }
#endif
   std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);

   const int data_id = data.getComponentDescriptorIndex(0);

   const int flux_id = (d_flux_id != -1) ? d_flux_id : d_flux_scratch_id;

   d_bc_helper.setTargetDataId(data_id);
   d_bc_helper.setHomogeneousBc(true);
   xeqScheduleGhostFillNoCoarse(data_id, ln);

   if (ln > d_ln_min) {
      /*
       * Perform a one-time transfer of data from coarser level,
       * to fill ghost boundaries that will not change through
       * the smoothing loop.
       */
      xeqScheduleGhostFill(data_id, ln);
   }

   /*
    * Smooth the number of sweeps specified or until
    * the convergence is satisfactory.
    */
   int isweep;
   double red_maxres, blk_maxres, maxres = 0;
   red_maxres = blk_maxres = residual_tolerance + 1;
   /*
    * Instead of checking residual convergence globally,
    * we check the not_converged flag.  This avoids possible
    * round-off errors affecting different processes differently,
    * leading to disagreement on whether to continue smoothing.
    */
   int not_converged = 1;
   for (isweep = 0; isweep < num_sweeps && not_converged; ++isweep) {
      // tbox::plog<<"isweep="<<isweep
      //          <<", num_sweep="<<num_sweeps
      //          <<", residual_tolerance="<<residual_tolerance<<std::endl;
      red_maxres = blk_maxres = 0;

      // Red sweep.
      xeqScheduleGhostFillNoCoarse(data_id, ln);
      for (hier::PatchLevel::Iterator pi(level->begin()); pi != level->end();
           pi++) {
         std::shared_ptr<hier::Patch> patch = *pi;

         bool deallocate_flux_data_when_done = false;
         if (flux_id == d_flux_scratch_id) {
            /*
             * Using internal temporary storage for flux.
             * For each patch, make sure the internal
             * side-centered data is allocated and note
             * whether that data should be deallocated when done.
             */
            if (!patch->checkAllocated(flux_id)) {
               patch->allocatePatchData(flux_id);
               deallocate_flux_data_when_done = true;
            }
         }

         std::shared_ptr<pdat::CellData<double> > err_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 data.getComponentPatchData(0, *patch)));
         std::shared_ptr<pdat::CellData<double> > residual_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 residual.getComponentPatchData(0, *patch)));
         std::shared_ptr<pdat::SideData<double> > flux_data(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(flux_id)));

         TBOX_ASSERT(err_data);
         TBOX_ASSERT(residual_data);
         TBOX_ASSERT(flux_data);
         TBOX_ASSERT(err_data->getDepth() == flux_data->getDepth());

         for (int depth = 0; depth < residual_data->getDepth(); depth++) {
            computeFluxOnPatch(*patch, level->getRatioToCoarserLevel(),
                               *err_data, *flux_data, depth);

            redOrBlackSmoothingOnPatch(*patch, *flux_data, *residual_data,
                                       *err_data, 'r', depth, &red_maxres);
         }
         if (deallocate_flux_data_when_done) {
            patch->deallocatePatchData(flux_id);
         }
      }  // End patch number pn
      xeqScheduleGhostFillNoCoarse(data_id, ln);

      // Black sweep.
      for (hier::PatchLevel::Iterator pi(level->begin()); pi != level->end();
           pi++) {
         std::shared_ptr<hier::Patch> patch = *pi;

         bool deallocate_flux_data_when_done = false;
         if (flux_id == d_flux_scratch_id) {
            /*
             * Using internal temporary storage for flux.
             * For each patch, make sure the internal
             * side-centered data is allocated and note
             * whether that data should be deallocated when done.
             */
            if (!patch->checkAllocated(flux_id)) {
               patch->allocatePatchData(flux_id);
               deallocate_flux_data_when_done = true;
            }
         }

         std::shared_ptr<pdat::CellData<double> > err_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 data.getComponentPatchData(0, *patch)));
         std::shared_ptr<pdat::CellData<double> > residual_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 residual.getComponentPatchData(0, *patch)));
         std::shared_ptr<pdat::SideData<double> > flux_data(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(flux_id)));

         TBOX_ASSERT(err_data);
         TBOX_ASSERT(residual_data);
         TBOX_ASSERT(flux_data);
         TBOX_ASSERT(residual_data->getDepth() == d_depth);

         for (int depth = 0; depth < residual_data->getDepth(); depth++) {
            computeFluxOnPatch(*patch, level->getRatioToCoarserLevel(),
                               *err_data, *flux_data, depth);

            redOrBlackSmoothingOnPatch(*patch, *flux_data, *residual_data,
                                       *err_data, 'b', depth, &blk_maxres);
         }
         if (deallocate_flux_data_when_done) {
            patch->deallocatePatchData(flux_id);
         }
      }  // End patch number pn
      xeqScheduleGhostFillNoCoarse(data_id, ln);
      if (residual_tolerance >= 0.0) {
         /*
           Check for early end of sweeps due to convergence
           only if it is numerically possible (user gave a
           non negative value for residual tolerance).
          */
         maxres = tbox::MathUtilities<double>::Max(red_maxres, blk_maxres);
         not_converged = maxres > residual_tolerance;
         const tbox::SAMRAI_MPI &mpi(d_hierarchy->getMPI());
         if (mpi.getSize() > 1) {
            mpi.AllReduce(&not_converged, 1, MPI_MAX);
         }
      }
   }  // End sweep number isweep
   if (d_enable_logging)
      tbox::plog << d_object_name << " RBGS smoothing maxres = " << maxres
                 << "\n"
                 << "  after " << isweep << " sweeps.\n";
}


/*
********************************************************************
* Fix flux on coarse-fine boundaries computed from a               *
* constant-refine interpolation of coarse level data.              *
********************************************************************
*/
void EllipticFACOps::ewingFixFlux(const hier::Patch &patch,
                                  const pdat::CellData<double> &soln_data,
                                  pdat::SideData<double> &flux_data,
                                  const hier::IntVector &ratio_to_coarser,
                                  const int depth) const
{
   TBOX_ASSERT_DIM_OBJDIM_EQUALITY4(tbox::Dimension(NDIM), patch, soln_data,
                                    flux_data, ratio_to_coarser);

   const int patch_ln = patch.getPatchLevelNumber();
   const hier::GlobalId id = patch.getGlobalId();
   std::shared_ptr<geom::CartesianGridGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry,
                              hier::BaseGridGeometry>(
           d_hierarchy->getGridGeometry()));
   const double *dx = patch_geom->getDx();
   const hier::Box &patch_box(patch.getBox());
   const hier::Index &plower = patch_box.lower();
   const hier::Index &pupper = patch_box.upper();

   const std::vector<hier::BoundaryBox> &bboxes =
       d_cf_boundary[patch_ln]->getBoundaries(id, 1);
   int bn;
   int nboxes = static_cast<int>(bboxes.size());

   const double *sol = soln_data.getPointer(depth);
   const double *flux0 = flux_data.getPointer(0, depth);
   const double *flux1 = flux_data.getPointer(1, depth);
#if NDIM == 3
   const double *flux2 = flux_data.getPointer(2, depth);
#endif

   if (d_poisson_spec[0].dIsVariable()) {

      std::shared_ptr<pdat::SideData<double> > diffcoef_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch.getPatchData(d_poisson_spec[0].getDPatchDataId())));

      for (bn = 0; bn < nboxes; ++bn) {
         const hier::BoundaryBox &boundary_box = bboxes[bn];
         TBOX_ASSERT(boundary_box.getBoundaryType() == 1);
         const hier::Box &bdry_box = boundary_box.getBox();
         const hier::Index &blower = bdry_box.lower();
         const hier::Index &bupper = bdry_box.upper();
         const int location_index = boundary_box.getLocationIndex();
#if NDIM == 2
         EFO_EWINGFIXFLUXVARDC2D(flux0, flux1,
                                 &flux_data.getGhostCellWidth()[0],
                                 &flux_data.getGhostCellWidth()[1],
                                 diffcoef_data->getPointer(0),
                                 diffcoef_data->getPointer(1),
                                 &diffcoef_data->getGhostCellWidth()[0],
                                 &diffcoef_data->getGhostCellWidth()[1], sol,
                                 &soln_data.getGhostCellWidth()[0],
                                 &soln_data.getGhostCellWidth()[1], &plower[0],
                                 &pupper[0], &plower[1], &pupper[1],
                                 &location_index, &ratio_to_coarser[0],
                                 &blower[0], &bupper[0], dx);
#else
         EFO_EWINGFIXFLUXVARDC3D(
             flux0, flux1, flux2, &flux_data.getGhostCellWidth()[0],
             &flux_data.getGhostCellWidth()[1],
             &flux_data.getGhostCellWidth()[2], diffcoef_data->getPointer(0),
             diffcoef_data->getPointer(1), diffcoef_data->getPointer(2),
             &diffcoef_data->getGhostCellWidth()[0],
             &diffcoef_data->getGhostCellWidth()[1],
             &diffcoef_data->getGhostCellWidth()[2], sol,
             &soln_data.getGhostCellWidth()[0],
             &soln_data.getGhostCellWidth()[1],
             &soln_data.getGhostCellWidth()[2], &plower[0], &pupper[0],
             &plower[1], &pupper[1], &plower[2], &pupper[2], &location_index,
             &ratio_to_coarser[0], &blower[0], &bupper[0], dx);
#endif
      }
   } else {

      const double diffcoef_constant = d_poisson_spec[0].getDConstant();

      for (bn = 0; bn < nboxes; ++bn) {
         const hier::BoundaryBox &boundary_box = bboxes[bn];
         TBOX_ASSERT(boundary_box.getBoundaryType() == 1);
         const hier::Box &bdry_box = boundary_box.getBox();
         const hier::Index &blower = bdry_box.lower();
         const hier::Index &bupper = bdry_box.upper();
         const int location_index = boundary_box.getLocationIndex();
#if NDIM == 2
         EFO_EWINGFIXFLUXCONDC2D(flux0, flux1,
                                 &flux_data.getGhostCellWidth()[0],
                                 &flux_data.getGhostCellWidth()[1],
                                 diffcoef_constant, sol,
                                 &soln_data.getGhostCellWidth()[0],
                                 &soln_data.getGhostCellWidth()[1], &plower[0],
                                 &pupper[0], &plower[1], &pupper[1],
                                 &location_index, &ratio_to_coarser[0],
                                 &blower[0], &bupper[0], dx);
#else
         EFO_EWINGFIXFLUXCONDC3D(
             flux0, flux1, flux2, &flux_data.getGhostCellWidth()[0],
             &flux_data.getGhostCellWidth()[1],
             &flux_data.getGhostCellWidth()[2], diffcoef_constant, sol,
             &soln_data.getGhostCellWidth()[0],
             &soln_data.getGhostCellWidth()[1],
             &soln_data.getGhostCellWidth()[2], &plower[0], &pupper[0],
             &plower[1], &pupper[1], &plower[2], &pupper[2], &location_index,
             &ratio_to_coarser[0], &blower[0], &bupper[0], dx);
#endif
      }
   }
}


/*
********************************************************************
* FACOperatorStrategy virtual solveCoarsestLevel             *
* function                                                         *
********************************************************************
*/
int EllipticFACOps::solveCoarsestLevel(
    solv::SAMRAIVectorReal<double> &data,
    const solv::SAMRAIVectorReal<double> &residual, int coarsest_ln)
{
   t_solve_coarsest->start();

   checkInputPatchDataIndices();

   int return_value = 0;

   if (d_coarse_solver_choice == "redblack") {
      d_residual_tolerance_during_smoothing = d_coarse_solver_tolerance;
      smoothError(data, residual, coarsest_ln, d_coarse_solver_max_iterations);
      d_residual_tolerance_during_smoothing = -1.0;
   } else if (d_coarse_solver_choice == "hypre") {
      return_value = solveCoarsestLevel_HYPRE(data, residual, coarsest_ln);
   } else {
      TBOX_ERROR(d_object_name << ": Bad coarse level solver choice '"
                               << d_coarse_solver_choice
                               << "' in "
                                  "scapCellPoissonOps::solveCoarsestLevel.");
   }

   xeqScheduleGhostFillNoCoarse(data.getComponentDescriptorIndex(0),
                                coarsest_ln);

   t_solve_coarsest->stop();

   return return_value;
}


/*******************************************************************
 * Solve coarsest level using Hypre                                 *
 * We only solve for the error, so we always use homogeneous bc.    *
 ********************************************************************/
int EllipticFACOps::solveCoarsestLevel_HYPRE(
    solv::SAMRAIVectorReal<double> &data,
    const solv::SAMRAIVectorReal<double> &residual, int coarsest_ln)
{
   NULL_USE(coarsest_ln);

   checkInputPatchDataIndices();

   int solver_ret = 1;
   /*
    * We use a different matrix for each depth.
    * We need to the the solver about which depth to use
    * out of "data" and "residual"
    */
   for (int i = 0; i < d_depth; i++) {
      d_hypre_solver[i]->setSolnIdDepth(i);
      d_hypre_solver[i]->setRhsIdDepth(i);
      d_hypre_solver[i]->setStoppingCriteria(d_coarse_solver_tolerance);
      int iret = d_hypre_solver[i]->solveSystem(
          data.getComponentDescriptorIndex(0),
          residual.getComponentDescriptorIndex(0), true, true);
      solver_ret = iret > 0 ? solver_ret : 0;
   }

   /*
    * Present data on the solve.
    * The Hypre solver returns 0 if converged.
    */
   if (d_enable_logging) {
      std::vector<CellPoissonHypreSolver *>::iterator it(
          d_hypre_solver.begin());
      for (; it != d_hypre_solver.end(); ++it)
         tbox::plog << d_object_name << " Hypre solve "
                    << (solver_ret ? "" : "NOT ") << "converged\n"
                    << "\titerations: " << (*it)->getNumberOfIterations()
                    << "\n"
                    << "\tresidual: " << (*it)->getRelativeResidualNorm()
                    << "\n";
   }
   return !solver_ret;
}


void EllipticFACOps::accumulateOperatorOnLevel(const int soln_id,
                                               const int accum_id, int ln,
                                               bool error_equation_indicator)
{
   t_accumulate_operator->start();

   std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);

   /*
    * Set up the bc helper so that when we use a refine schedule
    * to fill ghosts, the correct data is operated on.
    */

   d_bc_helper.setTargetDataId(soln_id);
   d_bc_helper.setHomogeneousBc(error_equation_indicator);

   const int flux_id = (d_flux_id != -1) ? d_flux_id : d_flux_scratch_id;

   /*
    * Assumptions:
    * 1. Data does not yet exist in ghost boundaries.
    * 2. Operator image data on next finer grid (if any)
    *    has been computed already.
    * 3. Flux data from next finer grid (if any) has
    *    been computed but has not been coarsened to
    *    this level.
    *
    * Steps:
    * S1. Fill solution ghost data by refinement
    *     or setting physical boundary conditions.
    *     This also brings in information from coarser
    *     to form the composite grid flux.
    * S2. Compute flux on ln.
    * S3. If next finer is available,
    *     Coarsen flux data on next finer level,
    *     overwriting flux computed from coarse data.
    * S4. Compute operaor image data from flux.
    */

   /* S1. Fill solution ghost data. */
   {
      if (ln > d_ln_min) {
         /* Fill from current, next coarser level and physical boundary */
         xeqScheduleGhostFill(soln_id, ln);
      } else {
         /* Fill from current and physical boundary */
         xeqScheduleGhostFillNoCoarse(soln_id, ln);
      }
   }

   /*
    * For the whole level, make sure the internal
    * side-centered data is allocated and note
    * whether that data should be deallocated when done.
    * We do this for the whole level because the data
    * undergoes transfer operations which require the
    * whole level data.
    */
   bool deallocate_flux_data_when_done = false;
   if (flux_id == d_flux_scratch_id) {
      if (!level->checkAllocated(flux_id)) {
         level->allocatePatchData(flux_id);
         deallocate_flux_data_when_done = true;
      }
   }

   /*
    * S2. Compute flux on patches in level.
    */
   for (hier::PatchLevel::Iterator pi(level->begin()); pi != level->end();
        pi++) {
      const std::shared_ptr<hier::Patch> patch = *pi;

      std::shared_ptr<pdat::CellData<double> > soln_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(soln_id)));
      std::shared_ptr<pdat::SideData<double> > flux_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(flux_id)));

      TBOX_ASSERT(soln_data);
      TBOX_ASSERT(flux_data);

      for (int depth = 0; depth < soln_data->getDepth(); depth++)
         computeFluxOnPatch(*patch, level->getRatioToCoarserLevel(), *soln_data,
                            *flux_data, depth);
   }

   /*
    * S3. Coarsen oflux data from next finer level so that
    * the computed flux becomes the composite grid flux.
    */
   if (ln < d_ln_max) {
      xeqScheduleFluxCoarsen(flux_id, d_oflux_scratch_id, ln);
   }

   /*
    * S4. Accumulate operator image on level.
    */
   for (hier::PatchLevel::Iterator pi(level->begin()); pi != level->end();
        pi++) {
      std::shared_ptr<hier::Patch> patch = *pi;
      std::shared_ptr<pdat::SideData<double> > flux_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(flux_id)));
      std::shared_ptr<pdat::CellData<double> > m_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_m_id)));
      std::shared_ptr<pdat::CellData<double> > soln_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(soln_id)));
      std::shared_ptr<pdat::CellData<double> > accum_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(accum_id)));

      TBOX_ASSERT(soln_data);
      TBOX_ASSERT(m_data);
      TBOX_ASSERT(accum_data);
      TBOX_ASSERT(flux_data);

      TBOX_ASSERT(flux_data->getDepth() == accum_data->getDepth());
      TBOX_ASSERT(flux_data->getDepth() == soln_data->getDepth());

      for (int depth = 0; depth < soln_data->getDepth(); depth++)
         accumulateOperatorOnPatch(*patch, *flux_data, *m_data, *soln_data,
                                   *accum_data, depth);

      if (ln > d_ln_min) {
         /*
          * Save outerflux data so that next coarser level
          *  can compute its coarse-fine composite flux.
          *  This is not strictly needed in this "compute residual"
          *  loop through the patches, but we put it here to
          *  avoid writing another loop for it.
          */
         std::shared_ptr<pdat::OutersideData<double> > oflux_data(
             SAMRAI_SHARED_PTR_CAST<pdat::OutersideData<double>,
                                    hier::PatchData>(
                 patch->getPatchData(d_oflux_scratch_id)));
         TBOX_ASSERT(oflux_data);
         oflux_data->copy(*flux_data);
      }
   }


   if (deallocate_flux_data_when_done) {
      level->deallocatePatchData(flux_id);
   }

   t_accumulate_operator->stop();
}

/*
********************************************************************
* FACOperatorStrategy virtual                                *
* computeCompositeResidualOnLevel function                         *
********************************************************************
*/
void EllipticFACOps::computeCompositeResidualOnLevel(
    solv::SAMRAIVectorReal<double> &residual,
    const solv::SAMRAIVectorReal<double> &solution,
    const solv::SAMRAIVectorReal<double> &rhs, int ln,
    bool error_equation_indicator)
{
   t_compute_composite_residual->start();

   checkInputPatchDataIndices();
#ifdef DEBUG_CHECK_ASSERTIONS
   if (residual.getPatchHierarchy() != d_hierarchy ||
       solution.getPatchHierarchy() != d_hierarchy ||
       rhs.getPatchHierarchy() != d_hierarchy) {
      TBOX_ERROR(d_object_name << ": Vector hierarchy does not match\n"
                                  "internal hierarchy.");
   }
#endif
   std::shared_ptr<hier::PatchLevel> level = d_hierarchy->getPatchLevel(ln);

   /*
    * Set up the bc helper so that when we use a refine schedule
    * to fill ghosts, the correct data is operated on.
    */
   const int soln_id = solution.getComponentDescriptorIndex(0);
   d_bc_helper.setTargetDataId(soln_id);
   d_bc_helper.setHomogeneousBc(error_equation_indicator);

   const int flux_id = (d_flux_id != -1) ? d_flux_id : d_flux_scratch_id;

   /*
    * Assumptions:
    * 1. Data does not yet exist in ghost boundaries.
    * 2. Residual data on next finer grid (if any)
    *    has been computed already.
    * 3. Flux data from next finer grid (if any) has
    *    been computed but has not been coarsened to
    *    this level.
    *
    * Steps:
    * S1. Fill solution ghost data by refinement
    *     or setting physical boundary conditions.
    *     This also brings in information from coarser
    *     to form the composite grid flux.
    * S2. Compute flux on ln.
    * S3. If next finer is available,
    *     Coarsen flux data on next finer level,
    *     overwriting flux computed from coarse data.
    * S4. Compute residual data from flux.
    */

   /* S1. Fill solution ghost data. */
   {
      if (ln > d_ln_min) {
         /* Fill from current, next coarser level and physical boundary */
         xeqScheduleGhostFill(soln_id, ln);
      } else {
         /* Fill from current and physical boundary */
         xeqScheduleGhostFillNoCoarse(soln_id, ln);
      }
   }

   /*
    * For the whole level, make sure the internal
    * side-centered data is allocated and note
    * whether that data should be deallocated when done.
    * We do this for the whole level because the data
    * undergoes transfer operations which require the
    * whole level data.
    */
   bool deallocate_flux_data_when_done = false;
   if (flux_id == d_flux_scratch_id) {
      if (!level->checkAllocated(flux_id)) {
         level->allocatePatchData(flux_id);
         deallocate_flux_data_when_done = true;
      }
   }

   /*
    * S2. Compute flux on patches in level.
    */
   for (hier::PatchLevel::Iterator pi(level->begin()); pi != level->end();
        pi++) {
      std::shared_ptr<hier::Patch> patch = *pi;

      std::shared_ptr<pdat::CellData<double> > soln_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              solution.getComponentPatchData(0, *patch)));
      std::shared_ptr<pdat::SideData<double> > flux_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(flux_id)));
      for (int depth = 0; depth < soln_data->getDepth(); depth++)
         computeFluxOnPatch(*patch, level->getRatioToCoarserLevel(), *soln_data,
                            *flux_data, depth);
   }

   /*
    * S3. Coarsen oflux data from next finer level so that
    * the computed flux becomes the composite grid flux.
    */
   if (ln < d_ln_max) {
      xeqScheduleFluxCoarsen(flux_id, d_oflux_scratch_id, ln);
   }

   /*
    * S4. Compute residual on patches in level.
    */
   for (hier::PatchLevel::Iterator pi(level->begin()); pi != level->end();
        pi++) {
      std::shared_ptr<hier::Patch> patch = *pi;
      std::shared_ptr<pdat::CellData<double> > soln_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              solution.getComponentPatchData(0, *patch)));
      std::shared_ptr<pdat::CellData<double> > m_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_m_id)));
      std::shared_ptr<pdat::CellData<double> > rhs_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              rhs.getComponentPatchData(0, *patch)));
      std::shared_ptr<pdat::CellData<double> > residual_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              residual.getComponentPatchData(0, *patch)));
      std::shared_ptr<pdat::SideData<double> > flux_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(flux_id)));

      TBOX_ASSERT(soln_data);
      TBOX_ASSERT(m_data);
      TBOX_ASSERT(rhs_data);
      TBOX_ASSERT(residual_data);
      TBOX_ASSERT(flux_data);

      for (int depth = 0; depth < soln_data->getDepth(); depth++)
         computeResidualOnPatch(*patch, *flux_data, *m_data, *soln_data,
                                *rhs_data, *residual_data, depth);

      if (ln > d_ln_min) {
         /*
          * Save outerflux data so that next coarser level
          *  can compute its coarse-fine composite flux.
          *  This is not strictly needed in this "compute residual"
          *  loop through the patches, but we put it here to
          *  avoid writing another loop for it.
          */
         std::shared_ptr<pdat::OutersideData<double> > oflux_data(
             SAMRAI_SHARED_PTR_CAST<pdat::OutersideData<double>,
                                    hier::PatchData>(
                 patch->getPatchData(d_oflux_scratch_id)));
         TBOX_ASSERT(oflux_data);
         oflux_data->copy(*flux_data);
      }
   }


   if (deallocate_flux_data_when_done) {
      level->deallocatePatchData(flux_id);
   }

   t_compute_composite_residual->stop();
}


/*
********************************************************************
* FACOperatorStrategy virtual computeResidualNorm             *
* function                                                         *
********************************************************************
*/
double EllipticFACOps::computeResidualNorm(
    const solv::SAMRAIVectorReal<double> &residual, int fine_ln, int coarse_ln)
{

   if (coarse_ln != residual.getCoarsestLevelNumber() ||
       fine_ln != residual.getFinestLevelNumber()) {
      TBOX_ERROR("EllipticFACOps::computeResidualNorm() is not\n"
                 << "set up to compute residual except on the range of\n"
                 << "levels defining the std::vector.\n");
   }
   t_compute_residual_norm->start();
   /*
    * The residual std::vector was cloned from std::vectors that has
    * the proper weights associated with them, so we do not
    * have to explicitly weight the residuals.
    *
    * maxNorm: not good to use because Hypre's norm does not
    *   correspond to it.  Also maybe too sensitive to spikes.
    * L2Norm: maybe good.  But does not correspond to the
    *   scale of the quantity.
    * L1Norm: maybe good.  Correspond to scale of quantity,
    *   but may be too insensitive to spikes.
    * RMSNorm: maybe good.
    */
   double norm = residual.RMSNorm();
   t_compute_residual_norm->stop();
   return norm;
}


/*
********************************************************************
* Compute the std::vector weight and put it at a specified patch data   *
* index.                                                           *
********************************************************************
*/
void EllipticFACOps::computeVectorWeights(
    std::shared_ptr<hier::PatchHierarchy> hierarchy, int weight_id,
    int coarsest_ln, int finest_ln) const
{
   TBOX_ASSERT(hierarchy);

   if (coarsest_ln == -1) coarsest_ln = 0;
   if (finest_ln == -1) finest_ln = hierarchy->getFinestLevelNumber();
   if (finest_ln < coarsest_ln) {
      TBOX_ERROR(d_object_name << ": Illegal level number range.  finest_ln < "
                                  "coarsest_ln.");
   }

   int ln;
   for (ln = finest_ln; ln >= coarsest_ln; --ln) {

      /*
       * On every level, first assign cell volume to std::vector weight.
       */

      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));
      for (hier::PatchLevel::iterator p(level->begin()); p != level->end();
           ++p) {
         const std::shared_ptr<hier::Patch> &patch = *p;
         std::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
             SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                                    hier::PatchGeometry>(
                 patch->getPatchGeometry()));
         const double *dx = patch_geometry->getDx();
         double cell_vol = dx[0];
         if (NDIM > 1) {
            cell_vol *= dx[1];
         }

         if (NDIM > 2) {
            cell_vol *= dx[2];
         }

         std::shared_ptr<pdat::CellData<double> > w(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(weight_id)));
         if (!w) {
            TBOX_ERROR(d_object_name << ": weight id must refer to a "
                                        "pdat::CellVariable");
         }
         w->fillAll(cell_vol);
      }

      /*
       * On all but the finest level, assign 0 to std::vector
       * weight to cells covered by finer cells.
       */

      if (ln < finest_ln) {

         /*
          * First get the boxes that describe index space of the next finer
          * level and coarsen them to describe corresponding index space
          * at this level.
          */

         std::shared_ptr<hier::PatchLevel> next_finer_level(
             hierarchy->getPatchLevel(ln + 1));
         hier::BoxContainer coarsened_boxes = next_finer_level->getBoxes();
         hier::IntVector coarsen_ratio(next_finer_level->getRatioToLevelZero());
         coarsen_ratio /= level->getRatioToLevelZero();
         coarsened_boxes.coarsen(coarsen_ratio);

         /*
          * Then set std::vector weight to 0 wherever there is
          * a nonempty intersection with the next finer level.
          * Note that all assignments are local.
          */

         for (hier::PatchLevel::iterator p(level->begin()); p != level->end();
              ++p) {

            const std::shared_ptr<hier::Patch> &patch = *p;
            for (hier::BoxContainer::iterator i = coarsened_boxes.begin();
                 i != coarsened_boxes.end(); ++i) {

               hier::Box intersection = *i * (patch->getBox());
               if (!intersection.empty()) {
                  std::shared_ptr<pdat::CellData<double> > w(
                      SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>,
                                             hier::PatchData>(
                          patch->getPatchData(weight_id)));
                  w->fillAll(0.0, intersection);

               }  // assignment only in non-empty intersection
            }     // loop over coarsened boxes from finer level
         }        // loop over patches in level
      }           // all levels except finest
   }              // loop over levels
}


/*
********************************************************************
* Check the validity and correctness of input data for this class. *
********************************************************************
*/
void EllipticFACOps::checkInputPatchDataIndices(const int depth) const
{
   /*
    * Check input validity and correctness.
    */
   hier::VariableDatabase &vdb(*hier::VariableDatabase::getDatabase());

   if (!d_poisson_spec[depth].dIsConstant() &&
       d_poisson_spec[depth].getDPatchDataId() != -1) {
      std::shared_ptr<hier::Variable> var;
      vdb.mapIndexToVariable(d_poisson_spec[depth].getDPatchDataId(), var);
      std::shared_ptr<pdat::SideVariable<double> > diffcoef_var(
          SAMRAI_SHARED_PTR_CAST<pdat::SideVariable<double>, hier::Variable>(
              var));
      if (!diffcoef_var) {
         TBOX_ERROR(d_object_name << ": Bad diffusion coefficient patch data "
                                     "index.");
      }
   }

   if (d_poisson_spec[depth].cIsVariable()) {
      std::shared_ptr<hier::Variable> var;
      vdb.mapIndexToVariable(d_poisson_spec[depth].getCPatchDataId(), var);
      std::shared_ptr<pdat::CellVariable<double> > scalar_field_var(
          SAMRAI_SHARED_PTR_CAST<pdat::CellVariable<double>, hier::Variable>(
              var));
      if (!scalar_field_var) {
         TBOX_ERROR(d_object_name << ": Bad linear term patch data index.");
      }
   }

   if (d_flux_id != -1) {
      std::shared_ptr<hier::Variable> var;
      vdb.mapIndexToVariable(d_flux_id, var);
      std::shared_ptr<pdat::SideVariable<double> > flux_var(
          SAMRAI_SHARED_PTR_CAST<pdat::SideVariable<double>, hier::Variable>(
              var));

      TBOX_ASSERT(flux_var);
   }
}

/*
*******************************************************************
*                                                                 *
* AMR-unaware patch-centered computational kernels.               *
*                                                                 *
*******************************************************************
*/
void EllipticFACOps::computeFluxOnPatch(
    const hier::Patch &patch, const hier::IntVector &ratio_to_coarser_level,
    const pdat::CellData<double> &w_data, pdat::SideData<double> &Dgradw_data,
    const int depth) const
{
   TBOX_ASSERT(patch.inHierarchy());
   TBOX_ASSERT(w_data.getGhostCellWidth() >=
               hier::IntVector::getOne(ratio_to_coarser_level.getDim()));
   TBOX_ASSERT(w_data.getDepth() > depth);
   TBOX_ASSERT(Dgradw_data.getDepth() > depth);

   std::shared_ptr<geom::CartesianGridGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianGridGeometry,
                              hier::BaseGridGeometry>(
           d_hierarchy->getGridGeometry()));
   TBOX_ASSERT(patch_geom);
   const hier::Box &box = patch.getBox();
   const int *lower = &box.lower()[0];
   const int *upper = &box.upper()[0];
   const double *dx = patch_geom->getDx();

   double *flux0 = Dgradw_data.getPointer(0, depth);
   double *flux1 = Dgradw_data.getPointer(1, depth);
#if NDIM == 3
   double *flux2 = Dgradw_data.getPointer(2, depth);
#endif

   if (d_poisson_spec[depth].dIsConstant()) {
      double D_value = d_poisson_spec[depth].getDConstant();
#if NDIM == 2
      SAMRAI_F77_FUNC(compfluxcondc2d, COMPFLUXCONDC2D)
      (flux0, flux1, &Dgradw_data.getGhostCellWidth()[0],
       &Dgradw_data.getGhostCellWidth()[1], D_value, w_data.getPointer(depth),
       &w_data.getGhostCellWidth()[0], &w_data.getGhostCellWidth()[1],
       &lower[0], &upper[0], &lower[1], &upper[1], dx);
#else
      SAMRAI_F77_FUNC(compfluxcondc3d, COMPFLUXCONDC3D)
      (flux0, flux1, flux2, &Dgradw_data.getGhostCellWidth()[0],
       &Dgradw_data.getGhostCellWidth()[1], &Dgradw_data.getGhostCellWidth()[2],
       D_value, w_data.getPointer(depth), &w_data.getGhostCellWidth()[0],
       &w_data.getGhostCellWidth()[1], &w_data.getGhostCellWidth()[2],
       &lower[0], &upper[0], &lower[1], &upper[1], &lower[2], &upper[2], dx);
#endif
   } else {
      std::shared_ptr<pdat::SideData<double> > D_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch.getPatchData(d_poisson_spec[depth].getDPatchDataId())));
      TBOX_ASSERT(D_data);
      TBOX_ASSERT(D_data->getDepth() == 1);
#if NDIM == 2
      EFO_COMPFLUXVARDC2D(flux0, flux1, &Dgradw_data.getGhostCellWidth()[0],
                          &Dgradw_data.getGhostCellWidth()[1],
                          D_data->getPointer(0), D_data->getPointer(1),
                          &D_data->getGhostCellWidth()[0],
                          &D_data->getGhostCellWidth()[1],
                          w_data.getPointer(depth),
                          &w_data.getGhostCellWidth()[0],
                          &w_data.getGhostCellWidth()[1], &lower[0], &upper[0],
                          &lower[1], &upper[1], dx);
#else
      EFO_COMPFLUXVARDC3D(
          flux0, flux1, flux2, &Dgradw_data.getGhostCellWidth()[0],
          &Dgradw_data.getGhostCellWidth()[1],
          &Dgradw_data.getGhostCellWidth()[2], D_data->getPointer(0),
          D_data->getPointer(1), D_data->getPointer(2),
          &D_data->getGhostCellWidth()[0], &D_data->getGhostCellWidth()[1],
          &D_data->getGhostCellWidth()[2], w_data.getPointer(depth),
          &w_data.getGhostCellWidth()[0], &w_data.getGhostCellWidth()[1],
          &w_data.getGhostCellWidth()[2], &lower[0], &upper[0], &lower[1],
          &upper[1], &lower[2], &upper[2], dx);
#endif
   }

   const int patch_ln = patch.getPatchLevelNumber();

   if (d_cf_discretization == "Ewing" && patch_ln > d_ln_min) {
      ewingFixFlux(patch, w_data, Dgradw_data, ratio_to_coarser_level, depth);
   }
}

/*
 * accum_data -= m*div(flux)-D*soln_data
 */
void EllipticFACOps::accumulateOperatorOnPatch(
    const hier::Patch &patch, const pdat::SideData<double> &flux_data,
    const pdat::CellData<double> &m_data,
    const pdat::CellData<double> &soln_data, pdat::CellData<double> &accum_data,
    const int depth) const
{
   assert(flux_data.getDepth() > depth);
   assert(soln_data.getDepth() > depth);
   assert(accum_data.getDepth() > depth);

   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));
   const hier::Box &box = patch.getBox();
   const hier::Index &lower = box.lower();
   const hier::Index &upper = box.upper();
   const double *dx = patch_geom->getDx();

   const double *sol = soln_data.getPointer(depth);

   const double *flux0 = flux_data.getPointer(0, depth);
   const double *flux1 = flux_data.getPointer(1, depth);
#if NDIM == 3
   const double *flux2 = flux_data.getPointer(2, depth);
#endif

   std::shared_ptr<pdat::CellData<double> > scalar_field_data;
   double scalar_field_constant;
   if (d_poisson_spec[0].cIsVariable()) {
      scalar_field_data =
          std::dynamic_pointer_cast<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(d_poisson_spec[0].getCPatchDataId()));
#if NDIM == 2
      ACCUMOPVARSCA2D(
          flux0, flux1, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], accum_data.getPointer(depth),
          &accum_data.getGhostCellWidth()[0],
          &accum_data.getGhostCellWidth()[1], scalar_field_data->getPointer(),
          &scalar_field_data->getGhostCellWidth()[0],
          &scalar_field_data->getGhostCellWidth()[1], m_data.getPointer(),
          &m_data.getGhostCellWidth()[0], &m_data.getGhostCellWidth()[1], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &lower[0], &upper[0], &lower[1], &upper[1], dx);
#else
      ACCUMOPVARSCA3D(
          flux0, flux1, flux2, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], &flux_data.getGhostCellWidth()[2],
          accum_data.getPointer(depth), &accum_data.getGhostCellWidth()[0],
          &accum_data.getGhostCellWidth()[1],
          &accum_data.getGhostCellWidth()[2], scalar_field_data->getPointer(),
          &scalar_field_data->getGhostCellWidth()[0],
          &scalar_field_data->getGhostCellWidth()[1],
          &scalar_field_data->getGhostCellWidth()[2], m_data.getPointer(),
          &m_data.getGhostCellWidth()[0], &m_data.getGhostCellWidth()[1],
          &m_data.getGhostCellWidth()[2], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &soln_data.getGhostCellWidth()[2], &lower[0], &upper[0], &lower[1],
          &upper[1], &lower[2], &upper[2], dx);
#endif
   } else if (d_poisson_spec[0].cIsConstant()) {
      scalar_field_constant = d_poisson_spec[0].getCConstant();
#if NDIM == 2
      ACCUMOPCONSCA2D(flux0, flux1, &flux_data.getGhostCellWidth()[0],
                      &flux_data.getGhostCellWidth()[1],
                      accum_data.getPointer(depth),
                      &accum_data.getGhostCellWidth()[0],
                      &accum_data.getGhostCellWidth()[1], scalar_field_constant,
                      m_data.getPointer(), &m_data.getGhostCellWidth()[0],
                      &m_data.getGhostCellWidth()[1], sol,
                      &soln_data.getGhostCellWidth()[0],
                      &soln_data.getGhostCellWidth()[1], &lower[0], &upper[0],
                      &lower[1], &upper[1], dx);
#else
      ACCUMOPCONSCA3D(
          flux0, flux1, flux2, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], &flux_data.getGhostCellWidth()[2],
          accum_data.getPointer(depth), &accum_data.getGhostCellWidth()[0],
          &accum_data.getGhostCellWidth()[1],
          &accum_data.getGhostCellWidth()[2], scalar_field_constant,
          m_data.getPointer(), &m_data.getGhostCellWidth()[0],
          &m_data.getGhostCellWidth()[1], &m_data.getGhostCellWidth()[2], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &soln_data.getGhostCellWidth()[2], &lower[0], &upper[0], &lower[1],
          &upper[1], &lower[2], &upper[2], dx);
#endif
   } else {
      scalar_field_constant = 0.0;
#if NDIM == 2
      ACCUMOPCONSCA2D(flux0, flux1, &flux_data.getGhostCellWidth()[0],
                      &flux_data.getGhostCellWidth()[1],
                      accum_data.getPointer(depth),
                      &accum_data.getGhostCellWidth()[0],
                      &accum_data.getGhostCellWidth()[1], 0.0,
                      m_data.getPointer(), &m_data.getGhostCellWidth()[0],
                      &m_data.getGhostCellWidth()[1], sol,
                      &soln_data.getGhostCellWidth()[0],
                      &soln_data.getGhostCellWidth()[1], &lower[0], &upper[0],
                      &lower[1], &upper[1], dx);
#else
      ACCUMOPCONSCA3D(
          flux0, flux1, flux2, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], &flux_data.getGhostCellWidth()[2],
          accum_data.getPointer(depth), &accum_data.getGhostCellWidth()[0],
          &accum_data.getGhostCellWidth()[1],
          &accum_data.getGhostCellWidth()[2], 0.0, m_data.getPointer(),
          &m_data.getGhostCellWidth()[0], &m_data.getGhostCellWidth()[1],
          &m_data.getGhostCellWidth()[2], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &soln_data.getGhostCellWidth()[2], &lower[0], &upper[0], &lower[1],
          &upper[1], &lower[2], &upper[2], dx);
#endif
   }
}


void EllipticFACOps::computeResidualOnPatch(
    const hier::Patch &patch, const pdat::SideData<double> &flux_data,
    const pdat::CellData<double> &m_data,
    const pdat::CellData<double> &soln_data,
    const pdat::CellData<double> &rhs_data,
    pdat::CellData<double> &residual_data, const int depth) const
{
   TBOX_ASSERT(flux_data.getDepth() > depth);
   TBOX_ASSERT(soln_data.getDepth() > depth);
   TBOX_ASSERT(rhs_data.getDepth() > depth);

   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));
   const hier::Box &box = patch.getBox();
   const hier::Index &lower = box.lower();
   const hier::Index &upper = box.upper();
   const double *dx = patch_geom->getDx();

   const double *flux0 = flux_data.getPointer(0, depth);
   const double *flux1 = flux_data.getPointer(1, depth);
#if NDIM == 3
   const double *flux2 = flux_data.getPointer(2, depth);
#endif

   double *res = residual_data.getPointer(depth);
   const double *rhs = rhs_data.getPointer(depth);
   const double *sol = soln_data.getPointer(depth);

   std::shared_ptr<pdat::CellData<double> > scalar_field_data;
   double scalar_field_constant;
   if (d_poisson_spec[0].cIsVariable()) {
      scalar_field_data =
          std::dynamic_pointer_cast<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(d_poisson_spec[0].getCPatchDataId()));
      TBOX_ASSERT(scalar_field_data);

#if NDIM == 2
      EFO_COMPRESVARSCA2D(flux0, flux1, &flux_data.getGhostCellWidth()[0],
                          &flux_data.getGhostCellWidth()[1], rhs,
                          &rhs_data.getGhostCellWidth()[0],
                          &rhs_data.getGhostCellWidth()[1], res,
                          &residual_data.getGhostCellWidth()[0],
                          &residual_data.getGhostCellWidth()[1],
                          scalar_field_data->getPointer(),
                          &scalar_field_data->getGhostCellWidth()[0],
                          &scalar_field_data->getGhostCellWidth()[1],
                          m_data.getPointer(), &m_data.getGhostCellWidth()[0],
                          &m_data.getGhostCellWidth()[1], sol,
                          &soln_data.getGhostCellWidth()[0],
                          &soln_data.getGhostCellWidth()[1], &lower[0],
                          &upper[0], &lower[1], &upper[1], dx);
#else
      EFO_COMPRESVARSCA3D(
          flux0, flux1, flux2, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], &flux_data.getGhostCellWidth()[2],
          rhs, &rhs_data.getGhostCellWidth()[0],
          &rhs_data.getGhostCellWidth()[1], &rhs_data.getGhostCellWidth()[2],
          res, &residual_data.getGhostCellWidth()[0],
          &residual_data.getGhostCellWidth()[1],
          &residual_data.getGhostCellWidth()[2],
          scalar_field_data->getPointer(),
          &scalar_field_data->getGhostCellWidth()[0],
          &scalar_field_data->getGhostCellWidth()[1],
          &scalar_field_data->getGhostCellWidth()[2], m_data.getPointer(),
          &m_data.getGhostCellWidth()[0], &m_data.getGhostCellWidth()[1],
          &m_data.getGhostCellWidth()[2], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &soln_data.getGhostCellWidth()[2], &lower[0], &upper[0], &lower[1],
          &upper[1], &lower[2], &upper[2], dx);
#endif
   } else if (d_poisson_spec[0].cIsConstant()) {
      scalar_field_constant = d_poisson_spec[0].getCConstant();
#if NDIM == 2
      EFO_COMPRESCONSCA2D(flux0, flux1, &flux_data.getGhostCellWidth()[0],
                          &flux_data.getGhostCellWidth()[1], rhs,
                          &rhs_data.getGhostCellWidth()[0],
                          &rhs_data.getGhostCellWidth()[1], res,
                          &residual_data.getGhostCellWidth()[0],
                          &residual_data.getGhostCellWidth()[1],
                          scalar_field_constant, m_data.getPointer(),
                          &m_data.getGhostCellWidth()[0],
                          &m_data.getGhostCellWidth()[1], sol,
                          &soln_data.getGhostCellWidth()[0],
                          &soln_data.getGhostCellWidth()[1], &lower[0],
                          &upper[0], &lower[1], &upper[1], dx);
#else
      EFO_COMPRESCONSCA3D(
          flux0, flux1, flux2, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], &flux_data.getGhostCellWidth()[2],
          rhs, &rhs_data.getGhostCellWidth()[0],
          &rhs_data.getGhostCellWidth()[1], &rhs_data.getGhostCellWidth()[2],
          res, &residual_data.getGhostCellWidth()[0],
          &residual_data.getGhostCellWidth()[1],
          &residual_data.getGhostCellWidth()[2], scalar_field_constant,
          m_data.getPointer(), &m_data.getGhostCellWidth()[0],
          &m_data.getGhostCellWidth()[1], &m_data.getGhostCellWidth()[2], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &soln_data.getGhostCellWidth()[2], &lower[0], &upper[0], &lower[1],
          &upper[1], &lower[2], &upper[2], dx);
#endif
   } else {
      scalar_field_constant = 0.0;
#if NDIM == 2
      EFO_COMPRESCONSCA2D(flux0, flux1, &flux_data.getGhostCellWidth()[0],
                          &flux_data.getGhostCellWidth()[1], rhs,
                          &rhs_data.getGhostCellWidth()[0],
                          &rhs_data.getGhostCellWidth()[1], res,
                          &residual_data.getGhostCellWidth()[0],
                          &residual_data.getGhostCellWidth()[1], 0.0,
                          m_data.getPointer(), &m_data.getGhostCellWidth()[0],
                          &m_data.getGhostCellWidth()[1], sol,
                          &soln_data.getGhostCellWidth()[0],
                          &soln_data.getGhostCellWidth()[1], &lower[0],
                          &upper[0], &lower[1], &upper[1], dx);
#else
      EFO_COMPRESCONSCA3D(
          flux0, flux1, flux2, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], &flux_data.getGhostCellWidth()[2],
          rhs, &rhs_data.getGhostCellWidth()[0],
          &rhs_data.getGhostCellWidth()[1], &rhs_data.getGhostCellWidth()[2],
          res, &residual_data.getGhostCellWidth()[0],
          &residual_data.getGhostCellWidth()[1],
          &residual_data.getGhostCellWidth()[2], 0.0, m_data.getPointer(),
          &m_data.getGhostCellWidth()[0], &m_data.getGhostCellWidth()[1],
          &m_data.getGhostCellWidth()[2], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &soln_data.getGhostCellWidth()[2], &lower[0], &upper[0], &lower[1],
          &upper[1], &lower[2], &upper[2], dx);
#endif
   }
}


void EllipticFACOps::evaluateRHS(const int soln_id, const int rhs_id)
{
   t_compute_rhs->start();

   // Initialize the output array
   d_hopscell->setToScalar(rhs_id, 0., false);

   for (int ln = d_ln_max; ln >= d_ln_min; ln--) {
      accumulateOperatorOnLevel(soln_id, rhs_id, ln, false);
   }

   t_compute_rhs->stop();
}


void EllipticFACOps::redOrBlackSmoothingOnPatch(
    const hier::Patch &patch, const pdat::SideData<double> &flux_data,
    const pdat::CellData<double> &rhs_data, pdat::CellData<double> &soln_data,
    char red_or_black, const int depth, double *p_maxres) const
{
   TBOX_ASSERT(red_or_black == 'r' || red_or_black == 'b');
   assert(d_m_id >= 0);

   const int offset = red_or_black == 'r' ? 0 : 1;
   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));
   const hier::Box &box = patch.getBox();
   const hier::Index &lower = box.lower();
   const hier::Index &upper = box.upper();
   const double *dx = patch_geom->getDx();

   std::shared_ptr<pdat::CellData<double> > scalar_field_data;
   double scalar_field_constant;
   std::shared_ptr<pdat::SideData<double> > diffcoef_data;
   double diffcoef_constant;

   if (d_poisson_spec[depth].cIsVariable()) {
      scalar_field_data =
          std::dynamic_pointer_cast<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(d_poisson_spec[depth].getCPatchDataId()));
   } else if (d_poisson_spec[depth].cIsConstant()) {
      scalar_field_constant = d_poisson_spec[depth].getCConstant();
   } else {
      scalar_field_constant = 0.0;
   }
   if (d_poisson_spec[depth].dIsVariable()) {
      diffcoef_data =
          std::dynamic_pointer_cast<pdat::SideData<double>, hier::PatchData>(
              patch.getPatchData(d_poisson_spec[depth].getDPatchDataId()));
   } else {
      diffcoef_constant = d_poisson_spec[depth].getDConstant();
   }
   std::shared_ptr<pdat::CellData<double> > m_data(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(d_m_id)));

   double maxres = 0.0;
   const double *rhs = rhs_data.getPointer(depth);
   const double *pm = m_data->getPointer();
   double *sol = soln_data.getPointer(depth);

   const double *flux0 = flux_data.getPointer(0, depth);
   const double *flux1 = flux_data.getPointer(1, depth);
#if NDIM == 3
   const double *flux2 = flux_data.getPointer(2, depth);
#endif

   if (d_poisson_spec[0].dIsVariable() && d_poisson_spec[depth].cIsVariable()) {
#if NDIM == 2
      EFO_RBGSWITHFLUXMAXVARDCVARSF2D(
          flux0, flux1, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], diffcoef_data->getPointer(0),
          diffcoef_data->getPointer(1), &diffcoef_data->getGhostCellWidth()[0],
          &diffcoef_data->getGhostCellWidth()[1], rhs,
          &rhs_data.getGhostCellWidth()[0], &rhs_data.getGhostCellWidth()[1],
          scalar_field_data->getPointer(),
          &scalar_field_data->getGhostCellWidth()[0],
          &scalar_field_data->getGhostCellWidth()[1], pm,
          &m_data->getGhostCellWidth()[0], &m_data->getGhostCellWidth()[1], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &lower[0], &upper[0], &lower[1], &upper[1], dx, &offset, &maxres);
#else
      EFO_RBGSWITHFLUXMAXVARDCVARSF3D(
          flux0, flux1, flux2, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], &flux_data.getGhostCellWidth()[2],
          diffcoef_data->getPointer(0), diffcoef_data->getPointer(1),
          diffcoef_data->getPointer(2), &diffcoef_data->getGhostCellWidth()[0],
          &diffcoef_data->getGhostCellWidth()[1],
          &diffcoef_data->getGhostCellWidth()[2], rhs,
          &rhs_data.getGhostCellWidth()[0], &rhs_data.getGhostCellWidth()[1],
          &rhs_data.getGhostCellWidth()[2], scalar_field_data->getPointer(),
          &scalar_field_data->getGhostCellWidth()[0],
          &scalar_field_data->getGhostCellWidth()[1],
          &scalar_field_data->getGhostCellWidth()[2], pm,
          &m_data->getGhostCellWidth()[0], &m_data->getGhostCellWidth()[1],
          &m_data->getGhostCellWidth()[2], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &soln_data.getGhostCellWidth()[2], &lower[0], &upper[0], &lower[1],
          &upper[1], &lower[2], &upper[2], dx, &offset, &maxres);
#endif
   } else if (d_poisson_spec[0].dIsVariable() &&
              d_poisson_spec[0].cIsConstant()) {
#if NDIM == 2
      EFO_RBGSWITHFLUXMAXVARDCCONSF2D(
          flux0, flux1, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], diffcoef_data->getPointer(0),
          diffcoef_data->getPointer(1), &diffcoef_data->getGhostCellWidth()[0],
          &diffcoef_data->getGhostCellWidth()[1], rhs,
          &rhs_data.getGhostCellWidth()[0], &rhs_data.getGhostCellWidth()[1],
          scalar_field_constant, pm, &m_data->getGhostCellWidth()[0],
          &m_data->getGhostCellWidth()[1], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &lower[0], &upper[0], &lower[1], &upper[1], dx, &offset, &maxres);
#else
      EFO_RBGSWITHFLUXMAXVARDCCONSF3D(
          flux0, flux1, flux2, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], &flux_data.getGhostCellWidth()[2],
          diffcoef_data->getPointer(0), diffcoef_data->getPointer(1),
          diffcoef_data->getPointer(2), &diffcoef_data->getGhostCellWidth()[0],
          &diffcoef_data->getGhostCellWidth()[1],
          &diffcoef_data->getGhostCellWidth()[2], rhs,
          &rhs_data.getGhostCellWidth()[0], &rhs_data.getGhostCellWidth()[1],
          &rhs_data.getGhostCellWidth()[2], scalar_field_constant, pm,
          &m_data->getGhostCellWidth()[0], &m_data->getGhostCellWidth()[1],
          &m_data->getGhostCellWidth()[2], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &soln_data.getGhostCellWidth()[2], &lower[0], &upper[0], &lower[1],
          &upper[1], &lower[2], &upper[2], dx, &offset, &maxres);
#endif
   } else if (d_poisson_spec[0].dIsVariable() && d_poisson_spec[0].cIsZero()) {
#if NDIM == 2
      EFO_RBGSWITHFLUXMAXVARDCCONSF2D(
          flux0, flux1, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], diffcoef_data->getPointer(0),
          diffcoef_data->getPointer(1), &diffcoef_data->getGhostCellWidth()[0],
          &diffcoef_data->getGhostCellWidth()[1], rhs,
          &rhs_data.getGhostCellWidth()[0], &rhs_data.getGhostCellWidth()[1],
          0.0, pm, &m_data->getGhostCellWidth()[0],
          &m_data->getGhostCellWidth()[1], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &lower[0], &upper[0], &lower[1], &upper[1], dx, &offset, &maxres);
#else
      EFO_RBGSWITHFLUXMAXVARDCCONSF3D(
          flux0, flux1, flux2, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], &flux_data.getGhostCellWidth()[2],
          diffcoef_data->getPointer(0), diffcoef_data->getPointer(1),
          diffcoef_data->getPointer(2), &diffcoef_data->getGhostCellWidth()[0],
          &diffcoef_data->getGhostCellWidth()[1],
          &diffcoef_data->getGhostCellWidth()[2], rhs,
          &rhs_data.getGhostCellWidth()[0], &rhs_data.getGhostCellWidth()[1],
          &rhs_data.getGhostCellWidth()[2], 0.0, pm,
          &m_data->getGhostCellWidth()[0], &m_data->getGhostCellWidth()[1],
          &m_data->getGhostCellWidth()[2], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &soln_data.getGhostCellWidth()[2], &lower[0], &upper[0], &lower[1],
          &upper[1], &lower[2], &upper[2], dx, &offset, &maxres);
#endif
   } else if (!d_poisson_spec[0].dIsVariable() &&
              d_poisson_spec[0].cIsVariable()) {
#if NDIM == 2
      EFO_RBGSWITHFLUXMAXCONDCVARSF2D(
          flux0, flux1, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], diffcoef_constant, rhs,
          &rhs_data.getGhostCellWidth()[0], &rhs_data.getGhostCellWidth()[1],
          scalar_field_data->getPointer(),
          &scalar_field_data->getGhostCellWidth()[0],
          &scalar_field_data->getGhostCellWidth()[1], pm,
          &m_data->getGhostCellWidth()[0], &m_data->getGhostCellWidth()[1], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &lower[0], &upper[0], &lower[1], &upper[1], dx, &offset, &maxres);
#else
      EFO_RBGSWITHFLUXMAXCONDCVARSF3D(
          flux0, flux1, flux2, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], &flux_data.getGhostCellWidth()[2],
          diffcoef_constant, rhs, &rhs_data.getGhostCellWidth()[0],
          &rhs_data.getGhostCellWidth()[1], &rhs_data.getGhostCellWidth()[2],
          scalar_field_data->getPointer(),
          &scalar_field_data->getGhostCellWidth()[0],
          &scalar_field_data->getGhostCellWidth()[1],
          &scalar_field_data->getGhostCellWidth()[2], pm,
          &m_data->getGhostCellWidth()[0], &m_data->getGhostCellWidth()[1],
          &m_data->getGhostCellWidth()[2], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &soln_data.getGhostCellWidth()[2], &lower[0], &upper[0], &lower[1],
          &upper[1], &lower[2], &upper[2], dx, &offset, &maxres);
#endif
   } else if (!d_poisson_spec[0].dIsVariable() &&
              d_poisson_spec[0].cIsConstant()) {
#if NDIM == 2
      EFO_RBGSWITHFLUXMAXCONDCCONSF2D(
          flux0, flux1, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], diffcoef_constant, rhs,
          &rhs_data.getGhostCellWidth()[0], &rhs_data.getGhostCellWidth()[1],
          scalar_field_constant, pm, &m_data->getGhostCellWidth()[0],
          &m_data->getGhostCellWidth()[1], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &lower[0], &upper[0], &lower[1], &upper[1], dx, &offset, &maxres);
#else
      EFO_RBGSWITHFLUXMAXCONDCCONSF3D(
          flux0, flux1, flux2, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], &flux_data.getGhostCellWidth()[2],
          diffcoef_constant, rhs, &rhs_data.getGhostCellWidth()[0],
          &rhs_data.getGhostCellWidth()[1], &rhs_data.getGhostCellWidth()[2],
          scalar_field_constant, pm, &m_data->getGhostCellWidth()[0],
          &m_data->getGhostCellWidth()[1], &m_data->getGhostCellWidth()[2], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &soln_data.getGhostCellWidth()[2], &lower[0], &upper[0], &lower[1],
          &upper[1], &lower[2], &upper[2], dx, &offset, &maxres);
#endif
   } else if (!d_poisson_spec[0].dIsVariable() && d_poisson_spec[0].cIsZero()) {
#if NDIM == 2
      EFO_RBGSWITHFLUXMAXCONDCCONSF2D(
          flux0, flux1, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], diffcoef_constant, rhs,
          &rhs_data.getGhostCellWidth()[0], &rhs_data.getGhostCellWidth()[1],
          0.0, pm, &m_data->getGhostCellWidth()[0],
          &m_data->getGhostCellWidth()[1], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &lower[0], &upper[0], &lower[1], &upper[1], dx, &offset, &maxres);
#else
      EFO_RBGSWITHFLUXMAXCONDCCONSF3D(
          flux0, flux1, flux2, &flux_data.getGhostCellWidth()[0],
          &flux_data.getGhostCellWidth()[1], &flux_data.getGhostCellWidth()[2],
          diffcoef_constant, rhs, &rhs_data.getGhostCellWidth()[0],
          &rhs_data.getGhostCellWidth()[1], &rhs_data.getGhostCellWidth()[2],
          0.0, pm, &m_data->getGhostCellWidth()[0],
          &m_data->getGhostCellWidth()[1], &m_data->getGhostCellWidth()[2], sol,
          &soln_data.getGhostCellWidth()[0], &soln_data.getGhostCellWidth()[1],
          &soln_data.getGhostCellWidth()[2], &lower[0], &upper[0], &lower[1],
          &upper[1], &lower[2], &upper[2], dx, &offset, &maxres);
#endif
   }

   *p_maxres = maxres;
}


void EllipticFACOps::xeqScheduleProlongation(int dst_id, int src_id, int scr_id,
                                             int dest_ln)
{
   if (!d_prolongation_refine_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
   }
   xfer::RefineAlgorithm refiner;

   refiner.registerRefine(dst_id, src_id, scr_id,
                          d_prolongation_refine_operator);
   refiner.resetSchedule(d_prolongation_refine_schedules[dest_ln]);
   d_prolongation_refine_schedules[dest_ln]->fillData(0.0);
   d_prolongation_refine_algorithm->resetSchedule(
       d_prolongation_refine_schedules[dest_ln]);
   return;
}


void EllipticFACOps::xeqScheduleURestriction(int dst_id, int src_id,
                                             int dest_ln)
{
   if (!d_urestriction_coarsen_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
   }
   xfer::CoarsenAlgorithm coarsener(tbox::Dimension(NDIM));
   coarsener.registerCoarsen(dst_id, src_id, d_urestriction_coarsen_operator);
   coarsener.resetSchedule(d_urestriction_coarsen_schedules[dest_ln]);
   d_urestriction_coarsen_schedules[dest_ln]->coarsenData();
   d_urestriction_coarsen_algorithm->resetSchedule(
       d_urestriction_coarsen_schedules[dest_ln]);
}


void EllipticFACOps::xeqScheduleRRestriction(int dst_id, int src_id,
                                             int dest_ln)
{
   if (!d_rrestriction_coarsen_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
   }
   xfer::CoarsenAlgorithm coarsener(tbox::Dimension(NDIM));
   coarsener.registerCoarsen(dst_id, src_id, d_rrestriction_coarsen_operator);
   coarsener.resetSchedule(d_rrestriction_coarsen_schedules[dest_ln]);
   d_rrestriction_coarsen_schedules[dest_ln]->coarsenData();
   d_rrestriction_coarsen_algorithm->resetSchedule(
       d_rrestriction_coarsen_schedules[dest_ln]);
}


void EllipticFACOps::xeqScheduleFluxCoarsen(int dst_id, int src_id, int dest_ln)
{
   if (!d_flux_coarsen_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
   }
   xfer::CoarsenAlgorithm coarsener(tbox::Dimension(NDIM));
   coarsener.registerCoarsen(dst_id, src_id, d_flux_coarsen_operator);
   coarsener.resetSchedule(d_flux_coarsen_schedules[dest_ln]);
   d_flux_coarsen_schedules[dest_ln]->coarsenData();
   d_flux_coarsen_algorithm->resetSchedule(d_flux_coarsen_schedules[dest_ln]);
}


void EllipticFACOps::xeqScheduleGhostFill(int dst_id, int dest_ln)
{
   if (!d_ghostfill_refine_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
   }
   xfer::RefineAlgorithm refiner;

   refiner.registerRefine(dst_id, dst_id, dst_id, d_ghostfill_refine_operator);
   refiner.resetSchedule(d_ghostfill_refine_schedules[dest_ln]);
   d_ghostfill_refine_schedules[dest_ln]->fillData(0.0);
   d_ghostfill_refine_algorithm->resetSchedule(
       d_ghostfill_refine_schedules[dest_ln]);
}

void EllipticFACOps::xeqScheduleGhostFillNoCoarse(int dst_id, int dest_ln)
{
   if (!d_ghostfill_nocoarse_refine_schedules[dest_ln]) {
      TBOX_ERROR("Expected schedule not found.");
   }
   xfer::RefineAlgorithm refiner;
   refiner.registerRefine(dst_id, dst_id, dst_id,
                          d_ghostfill_nocoarse_refine_operator);
   refiner.resetSchedule(d_ghostfill_nocoarse_refine_schedules[dest_ln]);
   d_ghostfill_nocoarse_refine_schedules[dest_ln]->fillData(0.0);
   d_ghostfill_nocoarse_refine_algorithm->resetSchedule(
       d_ghostfill_nocoarse_refine_schedules[dest_ln]);
}
