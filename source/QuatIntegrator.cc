// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#include "QuatIntegrator.h"

// external definitions for Fortran numerical routines
#include "QuatFort.h"
#include "ConcFort.h"

#include "PhaseFACSolver.h"
#include "EtaFACSolver.h"
#include "ConcFACSolver.h"
#include "TemperatureFACSolver.h"
#include "PhaseFACOps.h"
#include "EtaFACOps.h"
#include "ConcFACOps.h"
#include "TemperatureFACOps.h"
#include "QuatModel.h"
#include "FreeEnergyStrategy.h"
#include "QuatGradStrategy.h"
#include "QuatRefinePatchStrategy.h"
#include "QuatMobilityStrategy.h"
#include "CompositionRHSStrategy.h"
#include "PhaseFluxStrategy.h"
#include "PhaseConcentrationsStrategy.h"
#include "DeltaTemperatureFreeEnergyStrategy.h"
#include "UniformNoise.h"
#include "toolsSAMRAI.h"
#include "EBSCompositionRHSStrategy.h"
#include "PhaseRHSStrategyWithQ.h"
#include "FolchPlappRHSStrategy.h"
#include "SimpleTemperatureRHSStrategy.h"
#include "KKSCompositionRHSStrategy.h"
#include "GradientTemperatureStrategy.h"
#include "MultiOrderRHSStrategy.h"
#include "SinteringUWangRHSStrategy.h"
#include "WangSinteringCompositionRHSStrategy.h"

#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/math/PatchSideDataNormOpsReal.h"

#define BDF CV_BDF
#define ONE_STEP CV_ONE_STEP
#include "cvode/cvode_spils.h"

#define NUM_FACES 2 * NDIM

#ifndef TRUE
#define TRUE (1)
#endif

//-----------------------------------------------------------------------

void setBChomogeneous(solv::LocationIndexRobinBcCoefs* bc_coefs)
{
   tbox::plog << "BC for " << bc_coefs->getObjectName() << " set to homogeneous"
              << std::endl;
   for (int n = 0; n < 2 * NDIM; n++) {
      double a, b, g;
      bc_coefs->getCoefficients(n, a, b, g);
      // tbox::plog<<"BC for linear solver:"<<endl;
      // tbox::plog<<"old values: "<<a<<","<<b<<","<<g<<endl;
      g = 0.;
      // tbox::plog<<"new values: "<<a<<","<<b<<","<<g<<endl;
      bc_coefs->setRawCoefficients(n, a, b, g);
   }
}

//-----------------------------------------------------------------------

QuatIntegrator::QuatIntegrator(
    const std::string& name, const QuatModelParameters& model_parameters,
    QuatModel* model, const std::shared_ptr<hier::VariableContext> current,
    const std::shared_ptr<hier::VariableContext> scratch, const int qlen,
    const int ncompositions, std::shared_ptr<tbox::Database> db,
    std::shared_ptr<geom::CartesianGridGeometry> grid_geom,
    std::shared_ptr<tbox::Database> bc_db, const bool with_phase,
    const bool with_concentration, const bool with_heat_equation,
    const bool with_steady_temperature, const bool with_gradT,
    const bool with_antitrapping, const bool with_partition_coeff,
    const bool use_warm_start, const bool symmetry_aware,
    const bool use_gradq_for_flux)
    :  // protected data members
      d_phase_component_index(-1),
      d_eta_component_index(-1),
      d_quat_component_index(-1),
      d_conc_component_index(-1),
      d_temperature_component_index(-1),
      d_phase_id(-1),
      d_phase_scratch_id(-1),
      d_temperature_id(-1),
      d_temperature_scratch_id(-1),
      d_quat_id(-1),
      d_quat_scratch_id(-1),
      d_conc_id(-1),
      d_conc_scratch_id(-1),
      d_weight_id(-1),
      d_phase_mobility_id(-1),
      d_phase_temperature_mobility_id(-1),
      d_ncompositions(ncompositions),
      d_with_phase(with_phase),
      d_with_concentration(with_concentration),
      d_with_orientation(model_parameters.with_orientation()),
      d_evolve_quat(model_parameters.evolveQuat()),
      d_precond_has_dquatdphi(true),
      d_precond_has_dTdphi(false),
      d_precond_has_dPhidT(false),
      d_grid_geometry(grid_geom),
      d_lag_quat_sidegrad(true),
      d_uniform_diffusion_time_threshold(tbox::IEEE::getSignalingNaN()),
      d_conc_mobility(tbox::IEEE::getSignalingNaN()),
      d_show_conc_sys_stats(false),
      d_conc_bc_coefs(nullptr),
      d_dphidt_bc_coefs(nullptr),
      // private data members
      d_name(name),
      d_model_parameters(model_parameters),
      d_qlen(qlen),
      d_compute_velocity(false),
      d_current(current),
      d_scratch(scratch),
      d_quat_grad_strategy(nullptr),
      d_partition_coeff_strategy(nullptr),
      d_current_time(tbox::IEEE::getSignalingNaN()),
      d_previous_timestep(0.),
      d_eta_id(-1),
      d_eta_scratch_id(-1),
      d_cp_id(-1),
      d_quat_grad_cell_id(-1),
      d_quat_grad_side_id(-1),
      d_quat_grad_modulus_id(-1),
      d_eta_mobility_id(-1),
      d_quat_mobility_id(-1),
      d_quat_diffusion_deriv_id(-1),
      d_quat_symm_rotation_id(-1),
      d_conc_diffusion_id(-1),
      d_conc_phase_coupling_diffusion_id(-1),
      d_conc_eta_coupling_diffusion_id(-1),
      d_quat_diffs_id(-1),
      d_f_l_id(-1),
      d_f_a_id(-1),
      d_f_b_id(-1),
      d_phase_rhs_visit_id(-1),
      d_driving_force_visit_id(-1),
      d_q_rhs_visit_id(-1),
      d_modulus_q_rhs_visit_id(-1),
      d_temperature_rhs_visit_id(-1),
      d_conc_rhs_visit_id(-1),
      d_q_rhs1_visit_id(-1),
      d_phase_sol_id(-1),
      d_phase_rhs_id(-1),
      d_eta_sol_id(-1),
      d_eta_rhs_id(-1),
      d_temperature_sol_id(-1),
      d_temperature_rhs_id(-1),
      d_quat_sol_id(-1),
      d_quat_rhs_id(-1),
      d_u_sol_id(-1),
      d_u_rhs_id(-1),
      d_conc_sol_id(-1),
      d_conc_rhs_id(-1),
      d_dphidt_scratch_id(-1),
      d_quat_mobility_deriv_id(-1),
      d_flux_id(-1),
      d_flux_conc_id(-1),
      d_velocity_id(-1),
      d_quat_grad_side_copy_id(-1),
      d_tmp1_id(-1),
      d_tmp2_id(-1),
      d_conc_diffusion_coarsen(tbox::Dimension(NDIM)),
      d_flux_coarsen_algorithm(tbox::Dimension(NDIM)),
      d_flux_conc_coarsen_algorithm(tbox::Dimension(NDIM)),
      d_coarsen_alg(tbox::Dimension(NDIM)),
      d_H_parameter(tbox::IEEE::getSignalingNaN()),
      d_epsilon_phase(tbox::IEEE::getSignalingNaN()),
      d_epsilon_eta(tbox::IEEE::getSignalingNaN()),
      d_phase_well_scale(tbox::IEEE::getSignalingNaN()),
      d_eta_well_scale(tbox::IEEE::getSignalingNaN()),
      d_thermal_diffusivity(tbox::IEEE::getSignalingNaN()),
      d_latent_heat(tbox::IEEE::getSignalingNaN()),
      d_epsilon_q(tbox::IEEE::getSignalingNaN()),
      d_quat_grad_floor(tbox::IEEE::getSignalingNaN()),
      d_quat_smooth_floor_type("max"),
      d_sundials_solver(nullptr),
      d_end_time(tbox::IEEE::getSignalingNaN()),
      d_boundary_cond_db(bc_db),
      d_all_refine_patch_strategy(nullptr),
      d_show_integrator_stats(false),
      d_show_solver_stats(false),
      d_show_phase_sys_stats(false),
      d_show_eta_sys_stats(false),
      d_show_quat_sys_stats(false),
      d_show_temperature_sys_stats(false),
      d_use_preconditioner(true),
      d_max_precond_steps(1),
      d_cum_newton_iter(0),
      d_cum_lin_iter(0),
      d_cum_newton_fail(0),
      d_cum_lin_fail(0),
      d_cum_err_test_fail(0),
      d_cum_f_eval(0),
      d_cum_p_setup(0),
      d_cum_p_apply(0),
      d_phase_bc_coefs(nullptr),
      d_eta_bc_coefs(nullptr),
      d_temperature_bc_coefs(nullptr),
      d_quat_bc_coefs(nullptr)
{
   assert(db);
   assert(grid_geom);
   if (with_concentration) {
      assert(d_ncompositions <= 2);
      assert(d_ncompositions > 0);
   }

   d_quat_model = model;

   d_with_third_phase = false;
   d_with_heat_equation = with_heat_equation;
   d_with_steady_temperature = with_steady_temperature;
   d_with_gradT = with_gradT;
   d_with_antitrapping = with_antitrapping;
   d_with_partition_coeff = with_partition_coeff;

   d_with_unsteady_temperature =
       (d_with_heat_equation && !d_with_steady_temperature);

   if ((d_with_partition_coeff &&
        d_model_parameters.with_Aziz_partition_coeff()) ||
       d_model_parameters.with_velocity())
      d_compute_velocity = true;

   d_use_warm_start = use_warm_start;
   d_use_cvode = true;

   d_symmetry_aware = symmetry_aware;
   d_use_gradq_for_flux = use_gradq_for_flux;

   d_frame_velocity = d_model_parameters.movingVelocity();

   // Get timers
   tbox::TimerManager* tman = tbox::TimerManager::getManager();
   t_advance_timer = tman->getTimer("AMPE::QuatIntegrator::Advance()");
   t_rhs_timer = tman->getTimer("AMPE::QuatIntegrator::evaluateRHSFunction()");
   t_eta_rhs_timer = tman->getTimer("AMPE::QuatIntegrator::evaluateEtaRHS()");
   t_conc_rhs_timer =
       tman->getTimer("AMPE::QuatIntegrator::evaluateConcentrationRHS()");
   t_symm_rhs_timer =
       tman->getTimer("AMPE::QuatIntegrator::correctRhsForSymmetry()");
   t_set_coeff_timer =
       tman->getTimer("AMPE::QuatIntegrator::setCoefficients()");
   t_set_diffcoeff_conc_timer = tman->getTimer(
       "AMPE::QuatIntegrator::setDiffusionCoeffForConcentration()");
   t_psolve_setup_timer =
       tman->getTimer("AMPE::QuatIntegrator::CVSpgmrPrecondSet()");
   t_psolve_solve_timer =
       tman->getTimer("AMPE::QuatIntegrator::CVSpgmrPrecondSolve()");
   t_phase_precond_timer =
       tman->getTimer("AMPE::QuatIntegrator::PhasePrecondSolve()");
   t_conc_precond_timer =
       tman->getTimer("AMPE::QuatIntegrator::ConcPrecondSolve()");
   t_quat_grad_timer =
       tman->getTimer("AMPE::QuatIntegrator::computeQuatGradients()");
   t_phase_rhs_timer =
       tman->getTimer("AMPE::QuatIntegrator::evaluatePhaseRHS()");

   std::shared_ptr<tbox::Database> integrator_db =
       db->getDatabase("Integrator");

   d_show_integrator_stats =
       integrator_db->getBoolWithDefault("verbose", false);

   d_max_step_size = integrator_db->getDoubleWithDefault("max_step_size", 0.);

   d_use_cvode = integrator_db->getBoolWithDefault("use_cvode", true);
   d_atol = integrator_db->getDoubleWithDefault("atol", 3.e-4);
   if (integrator_db->keyExists("tolerance")) {
      d_atol = integrator_db->getDouble("tolerance");
   }
   d_rtol = integrator_db->getDoubleWithDefault("rtol", d_atol * 1.e-2);

   if (d_use_cvode)
      d_max_order = integrator_db->getIntegerWithDefault("max_order", 2);
   else
      d_max_order = integrator_db->getIntegerWithDefault("max_order", 4);

   d_uniform_diffusion_time_threshold =
       integrator_db->getDoubleWithDefault("uniform_diffusion_time_threshold",
                                           0.);

   d_lag_quat_sidegrad =
       integrator_db->getBoolWithDefault("lag_quat_sidegrad", true);

   d_max_krylov_dimension =
       integrator_db->getIntegerWithDefault("max_krylov_dimension", 5);

   d_integrator_db = integrator_db;


   tbox::Dimension dim(NDIM);
   const hier::IntVector& one_vec = hier::IntVector::getOne(dim);
   const hier::IntVector periodic = d_grid_geometry->getPeriodicShift(one_vec);

   d_all_periodic = true;
   for (int dd = 0; dd < NDIM; dd++) {
      d_all_periodic = d_all_periodic && periodic[dd];
   }
}

//-----------------------------------------------------------------------

QuatIntegrator::~QuatIntegrator()
{
   d_phase_sys_solver.reset();
   d_eta_sys_solver.reset();
   d_conc_sys_solver.reset();
   d_quat_sys_solver.reset();
   d_temperature_sys_solver.reset();
   d_grid_geometry.reset();
}

//-----------------------------------------------------------------------

void QuatIntegrator::setupPreconditionersPhase(
    std::shared_ptr<tbox::Database> integrator_db)
{
   std::shared_ptr<tbox::Database> phase_sys_solver_database;
   if (integrator_db->isDatabase("PhaseSysSolver")) {
      phase_sys_solver_database = integrator_db->getDatabase("PhaseSysSolver");
      d_show_phase_sys_stats =
          phase_sys_solver_database->getBoolWithDefault("verbose", false);
   }

   std::shared_ptr<PhaseFACOps> d_phase_fac_ops(
       new PhaseFACOps(d_name + "_QIPhaseFACOps", phase_sys_solver_database));

   d_phase_sys_solver.reset(new PhaseFACSolver(d_name + "QIPhaseSysSolver",
                                               d_phase_fac_ops,
                                               phase_sys_solver_database));
}

//-----------------------------------------------------------------------

void QuatIntegrator::setupPreconditionersConcentration(
    std::shared_ptr<tbox::Database> integrator_db)
{
   std::shared_ptr<tbox::Database> conc_sys_solver_database;
   if (integrator_db->isDatabase("ConcentrationSysSolver")) {
      conc_sys_solver_database =
          integrator_db->getDatabase("ConcentrationSysSolver");
      d_show_conc_sys_stats =
          conc_sys_solver_database->getBoolWithDefault("verbose", false);
   }

   std::shared_ptr<ConcFACOps> fac_ops(
       new ConcFACOps(d_name + "_QIConcFACOps", d_ncompositions,
                      conc_sys_solver_database));

   d_conc_sys_solver.reset(new ConcFACSolver(d_name + "_QIConcSysSolver",
                                             fac_ops,
                                             conc_sys_solver_database));
}

//-----------------------------------------------------------------------

void QuatIntegrator::setupPreconditionersEta(
    std::shared_ptr<tbox::Database> integrator_db)
{
   std::shared_ptr<tbox::Database> eta_sys_solver_database;
   if (integrator_db->isDatabase("EtaSysSolver")) {
      eta_sys_solver_database = integrator_db->getDatabase("EtaSysSolver");
      d_show_eta_sys_stats =
          eta_sys_solver_database->getBoolWithDefault("verbose", false);
   }

   std::shared_ptr<EtaFACOps> fac_ops(
       new EtaFACOps(d_name + "_QIEtaFACOps", eta_sys_solver_database));

   d_eta_sys_solver.reset(new EtaFACSolver(d_name + "_QIEtaSysSolver", fac_ops,
                                           eta_sys_solver_database));
}

//-----------------------------------------------------------------------

void QuatIntegrator::setupPreconditionersTemperature(
    std::shared_ptr<tbox::Database> integrator_db)
{
   std::shared_ptr<tbox::Database> temperature_sys_solver_database;
   if (integrator_db->isDatabase("TemperatureSysSolver")) {
      temperature_sys_solver_database =
          integrator_db->getDatabase("TemperatureSysSolver");
      d_show_temperature_sys_stats =
          temperature_sys_solver_database->getBoolWithDefault("verbose", false);
   }

   std::shared_ptr<TemperatureFACOps> d_temperature_fac_ops(
       new TemperatureFACOps(d_name + "_QITemperatureFACOps",
                             temperature_sys_solver_database));

   d_temperature_sys_solver.reset(
       new TemperatureFACSolver(d_name + "_QITemperatureSysSolver",
                                d_temperature_fac_ops,
                                temperature_sys_solver_database));

   if (d_precond_has_dTdphi) {
      d_phase_temperature_fac_ops.reset(
          new PhaseTemperatureFACOps(d_name + "_QIPhaseTemperatureFACOps",
                                     temperature_sys_solver_database));
      d_phase_temperature_sys_solver.reset(
          new EllipticFACSolver(d_name + "_QIPhaseTemperatureFACOSolver",
                                d_phase_temperature_fac_ops,
                                temperature_sys_solver_database));
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::setupPreconditioners()
{
   tbox::pout << "QuatIntegrator::setupPreconditioners()" << std::endl;

   d_use_preconditioner = true;
   bool precondition_phase = true;
   d_precondition_quat = true;
   bool precondition_eta = true;
   bool precondition_conc = true;
   bool precondition_temperature = d_with_steady_temperature ? false : true;

   if (d_integrator_db->isDatabase("Preconditioner")) {
      std::shared_ptr<tbox::Database> precond_db =
          d_integrator_db->getDatabase("Preconditioner");

      d_use_preconditioner = precond_db->getBoolWithDefault("enabled", true);

      if (d_use_preconditioner) {
         d_max_precond_steps =
             precond_db->getIntegerWithDefault("max_steps", 1);

         precondition_phase =
             precond_db->getBoolWithDefault("precondition_phase",
                                            precondition_phase);
         d_precondition_quat =
             precond_db->getBoolWithDefault("precondition_quat", true);
         precondition_eta = precond_db->getBoolWithDefault("precondition_eta",
                                                           precondition_eta);
         precondition_conc =
             precond_db->getBoolWithDefault("precondition_concentration",
                                            precondition_conc);
         precondition_temperature =
             precond_db->getBoolWithDefault("precondition_temperature",
                                            precondition_temperature);

         // Check to see if the preconditioner includes the Jacobian block
         // corresponding to the derivative of the quaternion components with
         // respect to the phase variable.
         bool default_value =
             d_model_parameters.concRHSstrategyIsKKS() ? true : false;
         d_precond_has_dquatdphi =
             precond_db->getBoolWithDefault("precond_has_dquatdphi",
                                            default_value);

         d_precond_has_dTdphi =
             precond_db->getBoolWithDefault("precond_has_dTdphi", false);
         d_precond_has_dPhidT =
             precond_db->getBoolWithDefault("precond_has_dPhidT", false);
      } else {
         precondition_phase = false;
         d_precondition_quat = false;
         precondition_eta = false;
         precondition_conc = false;
         precondition_temperature = false;
      }
   }

   if (!d_with_phase) d_precond_has_dquatdphi = false;

   // d_quat_sys_solver used not only for preconditioning, but also to evaluate
   // RHS
   if (d_evolve_quat) {

      std::shared_ptr<tbox::Database> quatsys_db;
      if (d_integrator_db->isDatabase("QuatSysSolver")) {
         quatsys_db = d_integrator_db->getDatabase("QuatSysSolver");
         d_show_quat_sys_stats =
             quatsys_db->getBoolWithDefault("verbose", false);
      }
      d_quat_face_coeff_strategy.reset(
          new QuatFaceCoeff(d_qlen, d_epsilon_q, d_quat_grad_floor,
                            d_quat_smooth_floor_type, d_H_parameter,
                            d_orient_interp_func_type1,
                            d_orient_interp_func_type2, d_avg_func_type));

      d_quat_sys_solver.reset(new QuatSysSolver(d_qlen,
                                                d_quat_face_coeff_strategy,
                                                d_name + "_QuatIntegratorQuatSy"
                                                         "sSolver",
                                                quatsys_db));

      assert(d_quat_sys_solver);
   } else {
      d_quat_sys_solver.reset();
   }

   if (d_use_preconditioner) {

      if (d_with_phase && precondition_phase) {

         setupPreconditionersPhase(d_integrator_db);
      } else {
         d_phase_sys_solver.reset();
      }

      if (d_with_third_phase && precondition_eta) {

         setupPreconditionersEta(d_integrator_db);
      } else {
         d_eta_sys_solver.reset();
      }

      if (d_with_concentration && precondition_conc) {

         setupPreconditionersConcentration(d_integrator_db);

      } else {
         d_conc_sys_solver.reset();
      }

      if (d_with_unsteady_temperature && precondition_temperature) {

         setupPreconditionersTemperature(d_integrator_db);
      } else {
         d_temperature_sys_solver.reset();
      }

   }  // d_use_preconditioner
   else {
      d_phase_sys_solver.reset();
      d_eta_sys_solver.reset();
      d_conc_sys_solver.reset();
      d_temperature_sys_solver.reset();
   }
}

//-----------------------------------------------------------------------
// Virtual function from QuatIntegrator (from StandardTagAndInitStrategy)

void QuatIntegrator::resetHierarchyConfiguration(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int coarsest_level, const int finest_level)
{
   tbox::plog << "QuatIntegrator::resetHierarchyConfiguration()" << std::endl;

   d_flux_coarsen_schedule.resize(hierarchy->getNumberOfLevels());
   d_flux_conc_coarsen_schedule.resize(hierarchy->getNumberOfLevels());

   d_conc_diffusion_coarsen_schedule.resize(hierarchy->getNumberOfLevels());

   int ln_beg = coarsest_level - (coarsest_level > 0);
   int ln_end = finest_level;
   for (int ln = ln_beg; ln < ln_end; ++ln) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
      std::shared_ptr<hier::PatchLevel> finer_level =
          hierarchy->getPatchLevel(ln + 1);

      d_conc_diffusion_coarsen_schedule[ln] =
          d_conc_diffusion_coarsen.createSchedule(level, finer_level);

      d_flux_coarsen_schedule[ln] =
          d_flux_coarsen_algorithm.createSchedule(level, finer_level);
      d_flux_conc_coarsen_schedule[ln] =
          d_flux_conc_coarsen_algorithm.createSchedule(level, finer_level);
   }

   if (!d_use_warm_start) {
      resetIntegrator(hierarchy, coarsest_level, finest_level);
#ifdef USE_CPODE
   } else {
      if (d_current_time > 0.0) {
         d_sundials_solver->reinitializeAfterRegrid();
         resetAfterRegrid(hierarchy, coarsest_level, finest_level);
      }
#endif
   }
}

//-----------------------------------------------------------------------
// Pure virtual function from QuatIntegrator
void QuatIntegrator::RegisterFreeEnergyVariables(
    const std::shared_ptr<pdat::CellVariable<double> > f_l_var,
    const std::shared_ptr<pdat::CellVariable<double> > f_a_var,
    const std::shared_ptr<pdat::CellVariable<double> > f_b_var)
{
   d_f_l_var = f_l_var;
   d_f_a_var = f_a_var;
   d_f_b_var = f_b_var;

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   d_f_l_id = variable_db->registerVariableAndContext(
       d_f_l_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 0));
   assert(d_f_l_id >= 0);

   d_f_a_id = variable_db->registerVariableAndContext(
       d_f_a_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 0));
   assert(d_f_a_id >= 0);

   if (d_f_b_var) {
      d_f_b_id = variable_db->registerVariableAndContext(
          d_f_b_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(d_f_b_id >= 0);
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::RegisterConcentrationVariables(
    const std::shared_ptr<pdat::CellVariable<double> > conc_var,
    const std::vector<std::shared_ptr<pdat::SideVariable<double> > >
        conc_pfm_diffusion_var,
    const std::shared_ptr<pdat::SideVariable<double> >
        conc_phase_coupling_diffusion_var,
    const std::shared_ptr<pdat::SideVariable<double> >
        conc_eta_coupling_diffusion_var,
    const std::shared_ptr<pdat::SideVariable<double> > conc_diffusion_var)
{
   assert(d_with_concentration);
   assert(d_ncompositions > 0);

   d_conc_var = conc_var;
   d_conc_pfm_diffusion_var = conc_pfm_diffusion_var;
   d_conc_diffusion_var = conc_diffusion_var;
   d_conc_phase_coupling_diffusion_var = conc_phase_coupling_diffusion_var;
   d_conc_eta_coupling_diffusion_var = conc_eta_coupling_diffusion_var;

   if (conc_var) {
      hier::VariableDatabase* variable_db =
          hier::VariableDatabase::getDatabase();

      d_conc_id = variable_db->registerVariableAndContext(
          d_conc_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(d_conc_id >= 0);

      d_conc_scratch_id =
          variable_db->mapVariableAndContextToIndex(d_conc_var, d_scratch);
      assert(d_conc_scratch_id >= 0);

      if (conc_diffusion_var) {
         d_conc_diffusion_id = variable_db->registerVariableAndContext(
             d_conc_diffusion_var, d_current,
             hier::IntVector(tbox::Dimension(NDIM), 0));
         assert(d_conc_diffusion_id >= 0);
      }
      for (std::vector<std::shared_ptr<pdat::SideVariable<double> > >::iterator
               it = d_conc_pfm_diffusion_var.begin();
           it != d_conc_pfm_diffusion_var.end(); ++it) {
         d_conc_pfm_diffusion_id.push_back(
             variable_db->mapVariableAndContextToIndex(*it, d_current));
         assert(d_conc_pfm_diffusion_id[0] >= 0);
      }

      if (d_conc_phase_coupling_diffusion_var) {
         d_conc_phase_coupling_diffusion_id =
             variable_db->registerVariableAndContext(
                 d_conc_phase_coupling_diffusion_var, d_current,
                 hier::IntVector(tbox::Dimension(NDIM), 0));
         assert(d_conc_phase_coupling_diffusion_id >= 0);
      }

      if (d_with_third_phase) {
         d_conc_eta_coupling_diffusion_id =
             variable_db->registerVariableAndContext(
                 d_conc_eta_coupling_diffusion_var, d_current,
                 hier::IntVector(tbox::Dimension(NDIM), 0));
         assert(d_conc_eta_coupling_diffusion_id >= 0);
      }

      if (!conc_pfm_diffusion_var.empty()) {
         // schedules
         std::shared_ptr<hier::CoarsenOperator> diff_coarsen_op =
             d_grid_geometry->lookupCoarsenOperator(d_conc_pfm_diffusion_var[0],
                                                    "CONSERVATIVE_COARSEN");

         for (std::vector<int>::iterator it = d_conc_pfm_diffusion_id.begin();
              it != d_conc_pfm_diffusion_id.end(); ++it) {
            d_conc_diffusion_coarsen.registerCoarsen(*it, *it, diff_coarsen_op);
         }

         if (d_conc_phase_coupling_diffusion_var)
            d_conc_diffusion_coarsen.registerCoarsen(
                d_conc_phase_coupling_diffusion_id,
                d_conc_phase_coupling_diffusion_id, diff_coarsen_op);

         if (d_with_third_phase) {
            d_conc_diffusion_coarsen.registerCoarsen(
                d_conc_eta_coupling_diffusion_id,
                d_conc_eta_coupling_diffusion_id, diff_coarsen_op);
         }
      }

      RegisterLocalConcentrationVariables();
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::RegisterQuatVariables(
    const std::shared_ptr<pdat::CellVariable<double> > quat_var,
    const std::shared_ptr<pdat::CellVariable<double> > quat_grad_cell_var,
    const std::shared_ptr<pdat::SideVariable<double> > quat_grad_side_var,
    const std::shared_ptr<pdat::CellVariable<double> > quat_grad_modulus_var,
    const std::shared_ptr<pdat::CellVariable<double> > quat_mobility_var,
    const std::shared_ptr<pdat::SideVariable<double> > quat_diffs_var,
    const std::shared_ptr<pdat::SideVariable<int> > quat_symm_rotation_var)
{
   assert(quat_var);

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   d_quat_id = variable_db->mapVariableAndContextToIndex(d_quat_var, d_current);
   d_quat_scratch_id =
       variable_db->mapVariableAndContextToIndex(d_quat_var, d_scratch);

   if (d_evolve_quat) {
      assert(quat_grad_cell_var);
      assert(quat_grad_side_var);
      assert(quat_grad_modulus_var);
      assert(quat_mobility_var);

      d_quat_grad_cell_var = quat_grad_cell_var;
      d_quat_grad_cell_id =
          variable_db->mapVariableAndContextToIndex(d_quat_grad_cell_var,
                                                    d_current);

      d_quat_grad_side_var = quat_grad_side_var;
      d_quat_grad_side_id =
          variable_db->mapVariableAndContextToIndex(d_quat_grad_side_var,
                                                    d_current);

      d_quat_grad_modulus_var = quat_grad_modulus_var;
      d_quat_grad_modulus_id =
          variable_db->mapVariableAndContextToIndex(d_quat_grad_modulus_var,
                                                    d_current);

      d_quat_mobility_var = quat_mobility_var;
      d_quat_mobility_id =
          variable_db->mapVariableAndContextToIndex(d_quat_mobility_var,
                                                    d_current);

      d_quat_diffs_var = quat_diffs_var;
      d_quat_diffs_id =
          variable_db->mapVariableAndContextToIndex(d_quat_diffs_var,
                                                    d_current);

      if (d_symmetry_aware) {
         assert(quat_symm_rotation_var);

         d_quat_symm_rotation_var = quat_symm_rotation_var;
         d_quat_symm_rotation_id =
             variable_db->mapVariableAndContextToIndex(d_quat_symm_rotation_var,
                                                       d_current);

         assert(d_quat_symm_rotation_id >= 0);
      }

      assert(d_quat_diffs_id >= 0);
      assert(d_quat_mobility_id >= 0);
      assert(d_quat_grad_cell_id >= 0);
      assert(d_quat_grad_side_id >= 0);
      assert(d_quat_grad_modulus_id >= 0);
   }
   assert(d_quat_scratch_id >= 0);
   assert(d_quat_id >= 0);
}

//-----------------------------------------------------------------------

void QuatIntegrator::RegisterVariables(
    const std::shared_ptr<pdat::CellVariable<double> > phase_var,
    const std::shared_ptr<pdat::CellVariable<double> > eta_var,
    const std::shared_ptr<pdat::CellVariable<double> > quat_var,
    const std::shared_ptr<pdat::CellVariable<double> > quat_grad_cell_var,
    const std::shared_ptr<pdat::SideVariable<double> > quat_grad_side_var,
    const std::shared_ptr<pdat::CellVariable<double> > quat_grad_modulus_var,
    const std::shared_ptr<pdat::CellVariable<double> > phase_mobility_var,
    const std::shared_ptr<pdat::CellVariable<double> > eta_mobility_var,
    const std::shared_ptr<pdat::CellVariable<double> > quat_mobility_var,
    const std::shared_ptr<pdat::SideVariable<double> > quat_diffs_var,
    const std::shared_ptr<pdat::SideVariable<int> > quat_symm_rotation_var,
    const std::shared_ptr<pdat::CellVariable<double> > weight_var,
    const std::shared_ptr<pdat::CellVariable<double> > temperature_var,
    const std::shared_ptr<pdat::CellVariable<double> > cp_var)
{
   tbox::pout << "QuatIntegrator::RegisterVariables()" << std::endl;

   if (d_with_phase) {
      assert(phase_var);
      assert(phase_mobility_var);
   }
   assert(weight_var);

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   // Variables owned by QuatModel
   if (phase_var) {
      d_phase_var = phase_var;
      d_phase_id =
          variable_db->mapVariableAndContextToIndex(d_phase_var, d_current);
      d_phase_scratch_id =
          variable_db->mapVariableAndContextToIndex(d_phase_var, d_scratch);
      assert(d_phase_id >= 0);
      assert(d_phase_scratch_id >= 0);
   }

   d_eta_var = eta_var;
   if (d_with_third_phase) {
      assert(eta_var);
      d_eta_id =
          variable_db->mapVariableAndContextToIndex(d_eta_var, d_current);
      d_eta_scratch_id =
          variable_db->mapVariableAndContextToIndex(d_eta_var, d_scratch);
      assert(d_eta_id >= 0);
      assert(d_eta_scratch_id >= 0);
   }

   d_temperature_var = temperature_var;
   d_temperature_id =
       variable_db->mapVariableAndContextToIndex(d_temperature_var, d_current);
   assert(d_temperature_id >= 0);
   d_temperature_scratch_id =
       variable_db->mapVariableAndContextToIndex(d_temperature_var, d_scratch);
   assert(d_temperature_scratch_id >= 0);

   d_cp_var = cp_var;
   if (d_with_unsteady_temperature) {
      assert(d_cp_var);
      d_cp_id = variable_db->mapVariableAndContextToIndex(d_cp_var, d_current);
      assert(d_cp_id >= 0);
   }

   d_quat_var = quat_var;
   if (d_with_orientation) {
      RegisterQuatVariables(quat_var, quat_grad_cell_var, quat_grad_side_var,
                            quat_grad_modulus_var, quat_mobility_var,
                            quat_diffs_var, quat_symm_rotation_var);
   }

   if (d_with_phase) {
      d_phase_mobility_var = phase_mobility_var;
      // 1 ghost layer needed for non-periodic BC
      d_phase_mobility_id =
          variable_db->mapVariableAndContextToIndex(d_phase_mobility_var,
                                                    d_current);
      assert(d_phase_mobility_id >= 0);

      if (d_precond_has_dTdphi) {
         d_phase_temperature_mobility_var.reset(
             new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                            d_name + "_QI_phase_Temperature_"
                                                     "mobility",
                                            1));
         d_phase_temperature_mobility_id =
             variable_db->registerVariableAndContext(
                 d_phase_temperature_mobility_var, d_current,
                 hier::IntVector(tbox::Dimension(NDIM), 1));
         d_local_data.setFlag(d_phase_temperature_mobility_id);
      }
   }

   if (d_with_third_phase) {
      d_eta_mobility_var = eta_mobility_var;
      d_eta_mobility_id = variable_db->registerVariableAndContext(
          d_eta_mobility_var, d_current,
          hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(d_eta_mobility_id >= 0);
   }

   d_weight_var = weight_var;
   d_weight_id = variable_db->registerVariableAndContext(
       d_weight_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 0));

   // Variables owned locally
   RegisterLocalVisitVariables();

   if (d_with_phase) {
      RegisterLocalPhaseVariables();
   }

   const int nphases = d_model_parameters.norderp();
   d_flux_var.reset(new pdat::SideVariable<double>(tbox::Dimension(NDIM),
                                                   d_name + "_QUI_flux_",
                                                   nphases));
   d_flux_id = variable_db->registerVariableAndContext(
       d_flux_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 0));
   assert(d_flux_id >= 0);
   d_local_data.setFlag(d_flux_id);

   if (d_with_concentration) {
      d_flux_conc_var.reset(new pdat::SideVariable<double>(
          tbox::Dimension(NDIM), d_name + "_QUI_conc_flux_", d_ncompositions));
      int nghosts = d_model_parameters.concRHSneeds2Ghosts() ? 1 : 0;
      d_flux_conc_id = variable_db->registerVariableAndContext(
          d_flux_conc_var, d_current,
          hier::IntVector(tbox::Dimension(NDIM), nghosts));
      assert(d_flux_conc_id >= 0);
      d_local_data.setFlag(d_flux_conc_id);
   }

   if (d_compute_velocity) {
      d_velocity_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                          d_name + "_velocity_",
                                                          1));
      d_velocity_id = variable_db->registerVariableAndContext(
          d_velocity_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(d_velocity_id >= 0);
      d_local_data.setFlag(d_velocity_id);
      d_tmp1_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                      d_name + "_tmp1_", 1));
      d_tmp1_id = variable_db->registerVariableAndContext(
          d_tmp1_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(d_tmp1_id >= 0);
      d_local_data.setFlag(d_tmp1_id);

      d_tmp2_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                      d_name + "_tmp2_", 1));
      d_tmp2_id = variable_db->registerVariableAndContext(
          d_tmp2_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(d_tmp2_id >= 0);
      d_local_data.setFlag(d_tmp2_id);
   }

   if (d_with_third_phase) {
      RegisterLocalEtaVariables();
   }

   d_flux_coarsen_op = d_grid_geometry->lookupCoarsenOperator(d_flux_var,
                                                              "CONSERVATIVE_"
                                                              "COARSEN");

   d_flux_coarsen_algorithm.registerCoarsen(d_flux_id, d_flux_id,
                                            d_flux_coarsen_op);

   if (d_with_concentration) {
      d_flux_conc_coarsen_algorithm.registerCoarsen(d_flux_conc_id,
                                                    d_flux_conc_id,
                                                    d_flux_coarsen_op);
   }
   if (d_with_orientation) {
      RegisterLocalQuatVariables();
   }

   if (d_with_unsteady_temperature) {
      RegisterLocalUnsteadyTemperatureVariables();
   }

   if (d_evolve_quat) {
      d_quat_grad_side_copy_var.reset(
          new pdat::SideVariable<double>(tbox::Dimension(NDIM),
                                         d_name + "_quat_grad_side_copy",
                                         NDIM * d_qlen));
      d_quat_grad_side_copy_id = variable_db->registerVariableAndContext(
          d_quat_grad_side_copy_var, d_current,
          hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(d_quat_grad_side_copy_id >= 0);
      d_local_data.setFlag(d_quat_grad_side_copy_id);
   }

   if (d_model_parameters.noise_amplitude() > 0.) {
      d_noise_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                       d_name + "_noise", 1));
      d_noise_id = variable_db->registerVariableAndContext(
          d_noise_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 0));
      d_local_data.setFlag(d_noise_id);
   }
}


void QuatIntegrator::RegisterLocalVisitVariables()
{
   if (d_model_parameters.with_rhs_visit_output()) {
      hier::VariableDatabase* variable_db =
          hier::VariableDatabase::getDatabase();

      if (d_with_phase) {
         const int nphases = d_model_parameters.norderp();
         d_phase_rhs_visit_var.reset(new pdat::CellVariable<double>(
             tbox::Dimension(NDIM), d_name + "_phase_rhs_visit_", nphases));
         d_phase_rhs_visit_id = variable_db->registerVariableAndContext(
             d_phase_rhs_visit_var, d_current,
             hier::IntVector(tbox::Dimension(NDIM), 0));
         assert(d_phase_rhs_visit_id >= 0);
         d_local_data.setFlag(d_phase_rhs_visit_id);

         d_driving_force_visit_var.reset(new pdat::CellVariable<double>(
             tbox::Dimension(NDIM), d_name + "_driving_force_visit_", nphases));
         d_driving_force_visit_id = variable_db->registerVariableAndContext(
             d_driving_force_visit_var, d_current,
             hier::IntVector(tbox::Dimension(NDIM), 0));
         assert(d_driving_force_visit_id >= 0);
         d_local_data.setFlag(d_driving_force_visit_id);
      }
      if (d_with_concentration) {
         d_conc_rhs_visit_var.reset(
             new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                            d_name + "_conc_rhs_visit_", 1));
         d_conc_rhs_visit_id = variable_db->registerVariableAndContext(
             d_conc_rhs_visit_var, d_current,
             hier::IntVector(tbox::Dimension(NDIM), 0));
         assert(d_conc_rhs_visit_id >= 0);
         d_local_data.setFlag(d_conc_rhs_visit_id);
      }
      if (d_with_heat_equation) {
         d_temperature_rhs_visit_var.reset(new pdat::CellVariable<double>(
             tbox::Dimension(NDIM), d_name + "_temperature_rhs_visit_", 1));
         d_temperature_rhs_visit_id = variable_db->registerVariableAndContext(
             d_temperature_rhs_visit_var, d_current,
             hier::IntVector(tbox::Dimension(NDIM), 0));
         assert(d_temperature_rhs_visit_id >= 0);
         d_local_data.setFlag(d_temperature_rhs_visit_id);
      }
      if (d_evolve_quat) {
         d_modulus_q_rhs_visit_var.reset(new pdat::CellVariable<double>(
             tbox::Dimension(NDIM), d_name + "_modulus_q_rhs_visit_", 1));
         d_modulus_q_rhs_visit_id = variable_db->registerVariableAndContext(
             d_modulus_q_rhs_visit_var, d_current,
             hier::IntVector(tbox::Dimension(NDIM), 0));
         assert(d_modulus_q_rhs_visit_id >= 0);
         d_local_data.setFlag(d_modulus_q_rhs_visit_id);

         d_q_rhs_visit_var.reset(
             new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                            d_name + "_q_rhs_visit_", d_qlen));
         d_q_rhs_visit_id = variable_db->registerVariableAndContext(
             d_q_rhs_visit_var, d_current,
             hier::IntVector(tbox::Dimension(NDIM), 0));
         assert(d_q_rhs_visit_id >= 0);
         d_local_data.setFlag(d_q_rhs_visit_id);

         d_q_rhs1_visit_var.reset(
             new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                            d_name + "_q_rhs1_visit_", 1));
         d_q_rhs1_visit_id = variable_db->registerVariableAndContext(
             d_q_rhs1_visit_var, d_current,
             hier::IntVector(tbox::Dimension(NDIM), 0));
         assert(d_q_rhs1_visit_id >= 0);
         d_local_data.setFlag(d_q_rhs1_visit_id);
      }
   }
}

void QuatIntegrator::RegisterLocalEtaVariables()
{
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   d_eta_sol_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                      d_name + "_QI_eta_sol_",
                                                      1));
   d_eta_sol_id = variable_db->registerVariableAndContext(
       d_eta_sol_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 1));
   d_local_data.setFlag(d_eta_sol_id);

   d_eta_rhs_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                      d_name + "_QI_eta_rhs_",
                                                      1));
   d_eta_rhs_id = variable_db->registerVariableAndContext(
       d_eta_rhs_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 0));
   d_local_data.setFlag(d_eta_rhs_id);

   assert(d_eta_rhs_id >= 0);
   assert(d_eta_sol_id >= 0);
}

void QuatIntegrator::RegisterLocalPhaseVariables()
{
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   const int nphases = d_model_parameters.norderp();

   // d_phase_sol and d_phase_rhs are used as temporary arrays for
   // preconditioner preconditioning is done one depth at a time, so a depth of
   // 1 is sufficient
   d_phase_sol_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                        d_name + "_QI_phase_"
                                                                 "sol_",
                                                        1));
   d_phase_sol_id = variable_db->registerVariableAndContext(
       d_phase_sol_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 1));
   assert(d_phase_sol_id >= 0);
   d_local_data.setFlag(d_phase_sol_id);

   d_phase_rhs_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                        d_name + "_QI_phase_"
                                                                 "rhs_",
                                                        1));
   d_phase_rhs_id = variable_db->registerVariableAndContext(
       d_phase_rhs_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 0));
   assert(d_phase_rhs_id >= 0);
   d_local_data.setFlag(d_phase_rhs_id);

   // dphidt is needed in ghost cells to compute
   // antitrapping fluxes
   if (needDphiDt()) {
      std::string name(d_name + "_QI_dphiddphidt");
      d_dphidt_scratch_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM), name, nphases));
      d_dphidt_scratch_id = variable_db->registerVariableAndContext(
          d_dphidt_scratch_var, d_scratch,
          hier::IntVector(tbox::Dimension(NDIM), 1));
      assert(d_dphidt_scratch_id >= 0);
      d_local_data.setFlag(d_dphidt_scratch_id);
   }
}

void QuatIntegrator::RegisterLocalConcentrationVariables()
{
   assert(d_with_concentration);
   assert(d_ncompositions > 0);
   assert(d_conc_var->getDepth() > 0);

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   d_conc_rhs_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                       d_name + "_QI_conc_rhs_",
                                                       d_ncompositions));
   d_conc_rhs_id = variable_db->registerVariableAndContext(
       d_conc_rhs_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 0));
   assert(d_conc_rhs_id >= 0);
   d_local_data.setFlag(d_conc_rhs_id);

   d_conc_sol_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                       d_name + "_QI_conc_sol_",
                                                       d_ncompositions));
   d_conc_sol_id = variable_db->registerVariableAndContext(
       d_conc_sol_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 1));
   assert(d_conc_sol_id >= 0);
   d_local_data.setFlag(d_conc_sol_id);

   assert(d_conc_rhs_var->getDepth() == d_conc_var->getDepth());
}

void QuatIntegrator::RegisterLocalQuatVariables()
{
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   d_quat_sol_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                       d_name + "_QI_quat_sol_",
                                                       d_qlen));
   d_quat_sol_id = variable_db->registerVariableAndContext(
       d_quat_sol_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 1));
   assert(d_quat_sol_id >= 0);
   d_local_data.setFlag(d_quat_sol_id);

   d_quat_rhs_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                       d_name + "_QI_quat_rhs_",
                                                       d_qlen));
   d_quat_rhs_id = variable_db->registerVariableAndContext(
       d_quat_rhs_var, d_current, hier::IntVector(tbox::Dimension(NDIM), 0));
   assert(d_quat_rhs_id >= 0);
   d_local_data.setFlag(d_quat_rhs_id);

   if (d_precond_has_dquatdphi) {

      d_quat_mobility_deriv_var.reset(new pdat::CellVariable<double>(
          tbox::Dimension(NDIM), d_name + "_QI_quat_mobility_deriv_", 1));
      d_quat_mobility_deriv_id = variable_db->registerVariableAndContext(
          d_quat_mobility_deriv_var, d_current,
          hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(d_quat_mobility_deriv_id >= 0);
      d_local_data.setFlag(d_quat_mobility_deriv_id);

      d_quat_diffusion_deriv_var.reset(
          new pdat::SideVariable<double>(tbox::Dimension(NDIM),
                                         d_name + "_QI_quat_diffusion_deriv_",
                                         2 * d_qlen));
      d_quat_diffusion_deriv_id = variable_db->registerVariableAndContext(
          d_quat_diffusion_deriv_var, d_current,
          hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(d_quat_diffusion_deriv_id >= 0);
      d_local_data.setFlag(d_quat_diffusion_deriv_id);
   }
}


void QuatIntegrator::RegisterLocalUnsteadyTemperatureVariables()
{
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   d_temperature_sol_var.reset(
       new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                      d_name + "_QI_temperature_sol_", 1));
   d_temperature_sol_id = variable_db->registerVariableAndContext(
       d_temperature_sol_var, d_current,
       hier::IntVector(tbox::Dimension(NDIM), 1));
   assert(d_temperature_sol_id >= 0);
   d_local_data.setFlag(d_temperature_sol_id);

   d_temperature_rhs_var.reset(
       new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                      d_name + "_QI_temperature_rhs_", 1));
   d_temperature_rhs_id = variable_db->registerVariableAndContext(
       d_temperature_rhs_var, d_current,
       hier::IntVector(tbox::Dimension(NDIM), 0));
   assert(d_temperature_rhs_id >= 0);
   d_local_data.setFlag(d_temperature_rhs_id);
}

// setup BC for linear solvers
// Values of g are set to 0 since these linear solvers compute corrections
void QuatIntegrator::setupBC()
{
   // boundary conditions
   if (!d_all_periodic) {
      if (d_with_phase) {
         std::shared_ptr<tbox::Database> phase_bc_db =
             d_boundary_cond_db->getDatabase("Phase");
         d_phase_bc_coefs =
             new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
                                                 "PhaseBcCoefs", phase_bc_db);
         setBChomogeneous(d_phase_bc_coefs);
      }
      if (needDphiDt()) {
         std::shared_ptr<tbox::Database> phase_bc_db =
             d_boundary_cond_db->getDatabase("Phase");

         assert(d_dphidt_scratch_id >= 0);
         d_dphidt_bc_coefs =
             new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
                                                 "DphiDtBcCoefs", phase_bc_db);
         // set BC values to 0 to have a 0 antitrapping flux at the boundary
         for (int i = 0; i < 2 * NDIM; i++) {
            d_dphidt_bc_coefs->setBoundaryValue(i, 0.);
         }
         d_dphidt_bc_helper =
             new solv::CartesianRobinBcHelper(tbox::Dimension(NDIM),
                                              "DphiDtBcHelper");
         d_dphidt_bc_helper->setTargetDataId(d_dphidt_scratch_id);
         d_dphidt_bc_helper->setCoefImplementation(d_dphidt_bc_coefs);
      }
      if (d_with_concentration) {
         std::shared_ptr<tbox::Database> conc_bc_db =
             d_boundary_cond_db->getDatabase("Conc");
         d_conc_bc_coefs =
             new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
                                                 "ConcBcCoefs", conc_bc_db);
         setBChomogeneous(d_conc_bc_coefs);
      }
      if (d_with_orientation) {
         std::shared_ptr<tbox::Database> quat_bc_db =
             d_boundary_cond_db->getDatabase("Quat");
         d_quat_bc_coefs =
             new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
                                                 "QuatBcCoefs", quat_bc_db);
         setBChomogeneous(d_quat_bc_coefs);
      }
      if (d_with_third_phase) {
         std::shared_ptr<tbox::Database> eta_bc_db =
             d_boundary_cond_db->getDatabase("Eta");
         d_eta_bc_coefs =
             new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
                                                 "EtaBcCoefs", eta_bc_db);
         setBChomogeneous(d_eta_bc_coefs);
      }
      if (d_with_unsteady_temperature) {
         std::shared_ptr<tbox::Database> temp_bc_db =
             d_boundary_cond_db->isDatabase("TemperatureCorrections")
                 ? d_boundary_cond_db->getDatabase("TemperatureCorrections")
                 : d_boundary_cond_db->getDatabase("Temperature");
         d_temperature_bc_coefs =
             new solv::LocationIndexRobinBcCoefs(tbox::Dimension(NDIM),
                                                 "LSTemperatureBcCoefs",
                                                 temp_bc_db);
         setBChomogeneous(d_temperature_bc_coefs);
      }
   }
}

//-----------------------------------------------------------------------
// Pure virtual function from QuatIntegrator

void QuatIntegrator::initializeCoarseRefineOperators(
    std::shared_ptr<mesh::GriddingAlgorithm> gridding_alg,
    std::shared_ptr<hier::RefineOperator> quat_refine_op,
    std::shared_ptr<hier::CoarsenOperator> quat_coarsen_op)
{
   // tbox::pout<<"QuatIntegrator::InitializeOperators()"<<endl;

   assert(gridding_alg);

   d_quat_refine_op = quat_refine_op;
   d_quat_coarsen_op = quat_coarsen_op;
   d_gridding_algorithm = gridding_alg;

   if (d_phase_id >= 0) {
      d_phase_refine_op =
          d_grid_geometry->lookupRefineOperator(d_phase_var, "LINEAR_REFINE");

      d_phase_coarsen_op = d_grid_geometry->lookupCoarsenOperator(d_phase_var,
                                                                  "CONSERVATIVE"
                                                                  "_COARSEN");
   }

   d_temperature_refine_op =
       d_grid_geometry->lookupRefineOperator(d_temperature_var,
                                             "LINEAR_"
                                             "REFINE");

   d_temperature_coarsen_op =
       d_grid_geometry->lookupCoarsenOperator(d_temperature_var,
                                              "CONSERVATIVE_COARSEN");

   if (d_with_third_phase) {
      d_eta_refine_op =
          d_grid_geometry->lookupRefineOperator(d_eta_var, "LINEAR_REFINE");
      d_eta_coarsen_op = d_grid_geometry->lookupCoarsenOperator(d_eta_var,
                                                                "CONSERVATIVE_"
                                                                "COARSEN");
   }

   if (d_conc_id >= 0) {
      d_conc_refine_op =
          d_grid_geometry->lookupRefineOperator(d_conc_var, "LINEAR_REFINE");

      d_conc_coarsen_op = d_grid_geometry->lookupCoarsenOperator(d_conc_var,
                                                                 "CONSERVATIVE_"
                                                                 "COARSEN");
   }

   setSolversBoundaries();

   if (!d_all_periodic) {
      double factor = d_model_parameters.with_rescaled_temperature()
                          ? 1. / d_model_parameters.rescale_factorT()
                          : -1.;
      const int phase_id = d_with_phase ? d_phase_scratch_id : -1;
      d_all_refine_patch_strategy =
          new QuatRefinePatchStrategy("QuatRefinePatchStrategy",
                                      d_boundary_cond_db, phase_id,
                                      d_eta_scratch_id, d_quat_scratch_id,
                                      d_conc_scratch_id,
                                      d_temperature_scratch_id, factor);
   }
}

void QuatIntegrator::setSolversBoundaries()
{
   // This should be okay with periodic boundaries
   // to satisfy solvers
   if (d_all_periodic) {
      if (d_phase_sys_solver) {
         d_phase_sys_solver->setBoundaries("Dirichlet");
      }
      if (d_eta_sys_solver && d_with_third_phase) {
         d_eta_sys_solver->setBoundaries("Dirichlet");
      }
      setConcentrationSolverBoundaries();
      if (d_evolve_quat) {
         assert(d_quat_sys_solver);
         d_quat_sys_solver->setBoundaries("Dirichlet");
      }
      if (d_temperature_sys_solver)
         if (d_with_unsteady_temperature) {
            d_temperature_sys_solver->setBoundaries("Dirichlet");
            if (d_precond_has_dTdphi)
               d_phase_temperature_sys_solver->setBoundaries("Dirichlet");
         }
   }
}

void QuatIntegrator::setConcentrationSolverBoundaries()
{
   if (d_conc_sys_solver && d_with_concentration) {
      d_conc_sys_solver->setBoundaries("Dirichlet");
   }
}

//-----------------------------------------------------------------------
// Pure virtual function from QuatIntegrator

void QuatIntegrator::setModelParameters(
    const double current_time, const double end_time, const double h_parameter,
    const double epsilon_phase, const double epsilon_eta,
    const double epsilon_q, const double quat_grad_floor,
    const std::string quat_smooth_floor_type, const double phase_well_scale,
    const double eta_well_scale, const std::string orient_interp_func_type1,
    const std::string orient_interp_func_type2, const std::string avg_func_type,
    const std::string phase_well_func_type,
    const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
    const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
    const std::string eta_well_func_type)
{
   tbox::plog << "QuatIntegrator::setModelParameters()" << std::endl;

   d_current_time = current_time;
   d_end_time = end_time;

   d_epsilon_phase = epsilon_phase;
   d_epsilon_eta = epsilon_eta;
   d_H_parameter = h_parameter;
   d_epsilon_q = epsilon_q;
   d_quat_grad_floor = quat_grad_floor;
   d_quat_smooth_floor_type = quat_smooth_floor_type;
   d_phase_well_scale = phase_well_scale;
   d_eta_well_scale = eta_well_scale;

   d_orient_interp_func_type1 = orient_interp_func_type1;
   d_orient_interp_func_type2 = orient_interp_func_type2;
   d_avg_func_type = avg_func_type;
   d_phase_well_func_type = phase_well_func_type;
   d_energy_interp_func_type = energy_interp_func_type;
   d_conc_interp_func_type = conc_interp_func_type;
   d_eta_well_func_type = eta_well_func_type;

   d_thermal_diffusivity = d_model_parameters.thermal_diffusivity();
   d_latent_heat = d_model_parameters.latent_heat();
   d_T_source = d_model_parameters.T_source();

   d_alpha_AT = d_epsilon_phase / sqrt(32. * d_phase_well_scale);

   setupPreconditioners();
}

//-----------------------------------------------------------------------

void QuatIntegrator::setConcentrationModelParameters(const double conc_mobility)
{
   d_conc_mobility = conc_mobility;
}

//-----------------------------------------------------------------------
// Pure virtual function from QuatIntegrator

void QuatIntegrator::RegisterWithVisit(
    std::shared_ptr<appu::VisItDataWriter> visit_data_writer)
{

   if (d_compute_velocity) {
#if 0
      assert( d_tmp1_id>=0 );
      std::string visit_name1("tmp1");
      visit_data_writer->registerPlotQuantity(
         visit_name1, "SCALAR", d_tmp1_id, 0 );
      assert( d_tmp2_id>=0 );
      std::string visit_name2("tmp2");
      visit_data_writer->registerPlotQuantity(
         visit_name2, "SCALAR", d_tmp2_id, 0 );
#endif
   }

   if (d_model_parameters.with_rhs_visit_output()) {
      if (d_with_phase) {
         assert(d_phase_rhs_visit_id >= 0);
         assert(d_driving_force_visit_id >= 0);
         const int n = d_model_parameters.norderp();
         for (int i = 0; i < n; i++) {
            visit_data_writer->registerPlotQuantity("phase_rhs" +
                                                        std::to_string(i),
                                                    "SCALAR",
                                                    d_phase_rhs_visit_id, i);

            visit_data_writer->registerPlotQuantity("driving_force" +
                                                        std::to_string(i),
                                                    "SCALAR",
                                                    d_driving_force_visit_id,
                                                    i);
         }
      }
      if (d_with_concentration) {
         assert(d_conc_rhs_visit_id >= 0);
         visit_data_writer->registerPlotQuantity("conc_rhs", "SCALAR",
                                                 d_conc_rhs_visit_id, 0);
      }
      if (d_with_heat_equation) {
         assert(d_temperature_rhs_visit_id > 0);
         visit_data_writer->registerPlotQuantity("temperature_rhs", "SCALAR",
                                                 d_temperature_rhs_visit_id, 0);
      }
      if (d_evolve_quat) {
         assert(d_modulus_q_rhs_visit_id > 0);
         std::string visit_nameq("modulus_q_rhs");
         visit_data_writer->registerPlotQuantity(visit_nameq, "SCALAR",
                                                 d_modulus_q_rhs_visit_id, 0);

         assert(d_q_rhs_visit_id > 0);
         for (int q = 0; q < d_qlen; q++) {
            std::string visit_namem("q_rhs" + std::to_string(q));
            visit_data_writer->registerPlotQuantity(visit_namem, "SCALAR",
                                                    d_q_rhs_visit_id, q);
         }

         assert(d_q_rhs1_visit_id > 0);
         std::string visit_name1("modulus_q_rhs1");
         visit_data_writer->registerPlotQuantity(visit_name1, "SCALAR",
                                                 d_q_rhs1_visit_id, 0);
      }
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::setQuatGradStrategy(QuatGradStrategy* quat_grad_strategy)
{
   assert(quat_grad_strategy != nullptr);
   d_quat_grad_strategy = quat_grad_strategy;
}

//-----------------------------------------------------------------------

void QuatIntegrator::setMobilityStrategy(
    std::shared_ptr<QuatMobilityStrategy>& mobility_strategy)
{
   d_mobility_strategy = mobility_strategy;
}

//-----------------------------------------------------------------------

void QuatIntegrator::setFreeEnergyStrategy(
    std::shared_ptr<FreeEnergyStrategy> free_energy_strategy)
{
   d_free_energy_strategy = free_energy_strategy;
}

//-----------------------------------------------------------------------

void QuatIntegrator::setCompositionRHSStrategy(
    std::shared_ptr<CompositionRHSStrategy> composition_rhs_strategy)
{
   assert(composition_rhs_strategy);
   d_composition_rhs_strategy = composition_rhs_strategy;
}

void QuatIntegrator::setCompositionDiffusionStrategy(
    std::shared_ptr<CompositionDiffusionStrategy>
        composition_diffusion_strategy)
{
   d_composition_diffusion_strategy = composition_diffusion_strategy;
}

void QuatIntegrator::setPhaseFluxStrategy(
    std::shared_ptr<PhaseFluxStrategy> phase_flux_strategy)
{
   assert(phase_flux_strategy);
   d_phase_flux_strategy = phase_flux_strategy;
}

//-----------------------------------------------------------------------

void QuatIntegrator::setTemperatureStrategy(
    std::shared_ptr<TemperatureStrategy> temperature_strategy)
{
   assert(temperature_strategy);
   d_temperature_strategy = temperature_strategy;
}

void QuatIntegrator::setHeatCapacityStrategy(
    HeatCapacityStrategy* heat_capacity_strategy)
{
   d_heat_capacity_strategy = heat_capacity_strategy;
}

//-----------------------------------------------------------------------
// Set solver options
void QuatIntegrator::setSundialsOptions()
{
   if (d_arkode_solver)
      setARKODEOptions();
   else {
      d_sundials_solver->setMaxPrecondSteps(d_max_precond_steps);
      d_sundials_solver->setLinearMultistepMethod(BDF);
      d_sundials_solver->setSteppingMethod(ONE_STEP);
      d_sundials_solver->setInitialStepSize(d_previous_timestep);
      d_sundials_solver->setRelativeTolerance(d_rtol);
      d_sundials_solver->setAbsoluteTolerance(d_atol);
      d_sundials_solver->setInitialValueOfIndependentVariable(d_current_time);
      d_sundials_solver->setMaximumLinearMultistepMethodOrder(d_max_order);

      if (d_use_preconditioner) {
         d_sundials_solver->setPreconditioningType(PREC_LEFT);
      }
      if (d_max_step_size > 0.) {
         d_sundials_solver->setMaximumAbsoluteStepSize(d_max_step_size);
      }

      bool needs_initialization = true;
      d_sundials_solver->setFinalValueOfIndependentVariable(
          d_end_time, needs_initialization);
      d_sundials_solver->setProjectionFunction(true);
      d_sundials_solver->setJTimesRhsFunction(true);
   }
}

void QuatIntegrator::setARKODEOptions()
{
   d_arkode_solver->setSteppingMethod(ARK_ONE_STEP);
   d_arkode_solver->setInitialStepSize(d_previous_timestep);
   d_arkode_solver->setRelativeTolerance(d_rtol);
   d_arkode_solver->setAbsoluteTolerance(d_atol);
   d_arkode_solver->setInitialValueOfIndependentVariable(d_current_time);
   d_arkode_solver->setMethodOrder(d_max_order);

   if (d_max_step_size > 0.) {
      d_arkode_solver->setMaximumAbsoluteStepSize(d_max_step_size);
   }

   bool needs_initialization = true;
   d_arkode_solver->setFinalValueOfIndependentVariable(d_end_time,
                                                       needs_initialization);
}

//-----------------------------------------------------------------------

void QuatIntegrator::createSundialsSolver()
{
   // Make a new solver
   if (d_use_cvode) {
      d_sundials_solver.reset(new CVODESolver(d_name + "_cvode_solver", this,
                                              d_use_preconditioner));
   } else {
      const int im_ex = 0;  // explicit
      const bool uses_preconditioner = false;
      d_arkode_solver.reset(new ARKODESolver(d_name + "_arkodesolver", this,
                                             uses_preconditioner, im_ex));
   };
}

void QuatIntegrator::resetSolutionVector(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   if (d_solution_vec) d_solution_vec.reset();

   d_solution_vec.reset(new solv::SAMRAIVectorReal<double>(
       d_name + "_solution_vector", hierarchy,
       0,  // even if all the levels have not changed, we still must do this for
           // all levels
       hierarchy->getFinestLevelNumber()));
}

//-----------------------------------------------------------------------
// Make a new solution std::vector
void QuatIntegrator::createSolutionvector(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   std::shared_ptr<hier::PatchLevel> level(
       hierarchy->getPatchLevel(hierarchy->getFinestLevelNumber()));
   assert(level);

   resetSolutionVector(hierarchy);

   int ncomponents = 0;
   if (d_with_phase) {
      assert(d_phase_id > -1);
      assert(level->checkAllocated(d_phase_id));
      d_solution_vec->addComponent(d_phase_var, d_phase_id, d_weight_id);
      d_phase_component_index = ncomponents;
      ncomponents++;
   }

   if (d_with_third_phase) {
      assert(d_eta_id > -1);
      d_solution_vec->addComponent(d_eta_var, d_eta_id, d_weight_id);
      d_eta_component_index = ncomponents;
      ncomponents++;
   }

   if (d_with_unsteady_temperature) {
      assert(d_temperature_id > -1);
      assert(level->checkAllocated(d_temperature_id));
      d_solution_vec->addComponent(d_temperature_var, d_temperature_id,
                                   d_weight_id);
      d_temperature_component_index = ncomponents;
      ncomponents++;
   }

   if (d_with_concentration) {
      assert(d_conc_id > -1);
      d_solution_vec->addComponent(d_conc_var, d_conc_id, d_weight_id);
      d_conc_component_index = ncomponents;
      ncomponents++;
   }

   if (d_evolve_quat) {
      assert(d_quat_id > -1);
      d_solution_vec->addComponent(d_quat_var, d_quat_id, d_weight_id);
      d_quat_component_index = ncomponents;
      ncomponents++;
   }
}

//-----------------------------------------------------------------------
/*
 * This function resets the integrator.
 */
void QuatIntegrator::resetIntegrator(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int coarsest_level, const int finest_level)
{
   tbox::plog << "QuatIntegrator::resetIntegrator()" << std::endl;

   assert(d_weight_id != -1);

   std::shared_ptr<hier::PatchLevel> level(
       hierarchy->getPatchLevel(finest_level));
   assert(level);

   // tbox::pout << "QuatIntegrator::resetIntegrator()" << std::endl;

   // Reset the linear solvers since the hierarchy has been changed
   resetSolversState(hierarchy, coarsest_level, finest_level);

   createSundialsSolver();

   createSolutionvector(hierarchy);

   setSundialsOptions();

   // Convert the SAMRAI solution std::vector to a Sundials std::vector
   solv::Sundials_SAMRAIVector* y = (solv::Sundials_SAMRAIVector*)
       solv::Sundials_SAMRAIVector::createSundialsVector(d_solution_vec);

   if (d_sundials_solver)
      d_sundials_solver->setInitialConditionVector(y);
   else
      d_arkode_solver->setInitialConditionVector(y);

   /*
     Complete the initialization of the integrator having now set the desired
     options.

     NB: The argument y is the std::vector in which we want
     the integrator to store its current solution, which may or may not
     be the same std::vector passed to the preceding initial condition
     setting function.  y must have the same structure as the
     initial condition std::vector, however, since the integrator uses the
     initial condition std::vector as a template to clone its
     internal std::vectors.
   */

   if (d_sundials_solver)
      d_sundials_solver->initialize(y);
   else
      d_arkode_solver->initialize(y);

   // Reset the cumulative counters
   d_cum_newton_iter = 0;
   d_cum_lin_iter = 0;
   d_cum_newton_fail = 0;
   d_cum_lin_fail = 0;
   d_cum_err_test_fail = 0;
   d_cum_f_eval = 0;
   d_cum_p_setup = 0;
   d_cum_p_apply = 0;
}

//-----------------------------------------------------------------------

void QuatIntegrator::resetAfterRegrid(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int coarsest_level, const int finest_level)
{
   // tbox::pout<<"QuatIntegrator::resetAfterRegrid()"<<endl;

   assert(d_weight_id != -1);

   // tbox::pout << "QuatIntegrator::resetAfterRegrid()" << std::endl;
   /*
    * This function resets the parts of the integrator needed after regrid
    */

   // Reset the linear solvers since the hierarchy has been changed
   resetSolversState(hierarchy, coarsest_level, finest_level);
}

//-----------------------------------------------------------------------

void QuatIntegrator::resetSolversStateConcentration(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   if (d_with_concentration && d_conc_sys_solver) {
      d_conc_sys_solver->resetSolverState(d_conc_sol_id, d_conc_rhs_id,
                                          hierarchy);
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::resetSolversState(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int coarsest_level, const int finest_level)
{
   tbox::plog << "QuatIntegrator::resetSolversState()" << std::endl;

   if (d_with_phase && d_phase_sys_solver) {
      assert(d_phase_sys_solver);
      d_phase_sys_solver->resetSolverState(d_phase_sol_id, d_phase_rhs_id,
                                           hierarchy);
   }

   if (d_with_third_phase && d_eta_sys_solver) {
      d_eta_sys_solver->resetSolverState(d_eta_sol_id, d_eta_rhs_id, hierarchy);
   }

   if (d_with_unsteady_temperature && d_temperature_sys_solver) {
      d_temperature_sys_solver->resetSolverState(d_temperature_sol_id,
                                                 d_temperature_rhs_id,
                                                 hierarchy);

      if (d_precond_has_dTdphi)
         d_phase_temperature_sys_solver->resetSolverState(d_temperature_sol_id,
                                                          d_temperature_rhs_id,
                                                          hierarchy);
   }

   resetSolversStateConcentration(hierarchy);

   if (d_evolve_quat) {
      assert(d_quat_sys_solver);
      d_quat_sys_solver->resetSolverState(d_quat_sol_id, d_quat_rhs_id,
                                          hierarchy);
   }


   if (d_with_steady_temperature && d_temperature_strategy) {
      d_temperature_strategy->resetSolversState(hierarchy, coarsest_level,
                                                finest_level);
   }
}

//-----------------------------------------------------------------------
// Inherited from StandardTagAndInitStrategy

void QuatIntegrator::initializeLevelData(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int level_number, const double time, const bool can_be_refined,
    const bool initial_time, const std::shared_ptr<hier::PatchLevel>& old_level,
    const bool allocate_data)
{
   // tbox::pout<<"QuatIntegrator::initializeLevelData()"<<endl;

   std::shared_ptr<hier::PatchLevel> level =
       hierarchy->getPatchLevel(level_number);

   if (level_number > 0) {
      if (old_level) old_level->deallocatePatchData(d_local_data);
   }

   level->allocatePatchData(d_local_data);

   if (!initial_time && d_use_warm_start) {

      // allocate CPODES internal std::vectors for warm start

      std::set<int> cpodes_id_set;
      std::set<int> phase_id_set;
      std::set<int> eta_id_set;
      std::set<int> orient_id_set;
      std::set<int> conc_id_set;
      std::set<int> temp_id_set;
#ifdef USE_CPODE
      getCPODESIdsRequiringRegrid(cpodes_id_set, phase_id_set, eta_id_set,
                                  orient_id_set, conc_id_set, temp_id_set);
#endif
      std::set<int>::iterator it;

      for (it = cpodes_id_set.begin(); it != cpodes_id_set.end(); it++) {
         level->allocatePatchData(*it, time);
      }

   }  // if ( !initial_time )
}

//-----------------------------------------------------------------------

void QuatIntegrator::initializeConcentrationNonPeriodicBC()
{
   if (d_conc_sys_solver) {
      assert(d_conc_bc_coefs != nullptr);
      d_conc_sys_solver->setBcObject(d_conc_bc_coefs);
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::initializeNonPeriodicBC()
{
   if (!d_all_periodic) {
      if (d_phase_sys_solver) {
         assert(d_phase_bc_coefs != nullptr);
         d_phase_sys_solver->setBcObject(d_phase_bc_coefs);
      }
      if (d_eta_sys_solver) {
         assert(d_eta_bc_coefs != nullptr);
         d_eta_sys_solver->setBcObject(d_eta_bc_coefs);
      }
      initializeConcentrationNonPeriodicBC();
      if (d_evolve_quat) {
         assert(d_quat_bc_coefs != nullptr);
         d_quat_sys_solver->setBcObject(d_quat_bc_coefs);
      }
      if (d_temperature_sys_solver) {
         assert(d_temperature_bc_coefs != nullptr);
         d_temperature_sys_solver->setBcObject(d_temperature_bc_coefs);
         if (d_precond_has_dTdphi)
            d_phase_temperature_sys_solver->setBcObject(d_temperature_bc_coefs);
      }
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::initializeConcentrationSolver(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
   if (d_with_concentration && d_conc_sys_solver) {
      int finest = hierarchy->getFinestLevelNumber();
      d_conc_sys_solver->initializeSolverState(d_conc_sol_id, d_conc_rhs_id,
                                               hierarchy, 0, finest);
   }
}
//-----------------------------------------------------------------------

void QuatIntegrator::initializeSolvers(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
   int finest = hierarchy->getFinestLevelNumber();

   if (d_with_phase && d_phase_sys_solver) {
      d_phase_sys_solver->initializeSolverState(d_phase_sol_id, d_phase_rhs_id,
                                                hierarchy, 0, finest);
   }

   if (d_with_third_phase && d_eta_sys_solver) {
      d_eta_sys_solver->initializeSolverState(d_eta_sol_id, d_eta_rhs_id,
                                              hierarchy, 0, finest);
   }

   if (d_with_heat_equation && d_temperature_sys_solver) {
      assert(d_temperature_sol_id >= 0);
      assert(d_temperature_rhs_id >= 0);
      d_temperature_sys_solver->initializeSolverState(d_temperature_sol_id,
                                                      d_temperature_rhs_id,
                                                      hierarchy, 0, finest);

      if (d_precond_has_dTdphi)
         d_phase_temperature_sys_solver->initializeSolverState(
             d_temperature_sol_id, d_temperature_rhs_id, hierarchy, 0, finest);
   }

   initializeConcentrationSolver(hierarchy);

   if (d_evolve_quat) {
      assert(d_quat_sys_solver);
      d_quat_sys_solver->initializeSolverState(d_quat_sol_id, d_quat_rhs_id,
                                               hierarchy);
   }
}

//-----------------------------------------------------------------------

// The way things are now (MEW, 08-08-04) the ::resetIntegrator method
// will have been called from resetHierarchyConfiguration before this
// is ever called.

void QuatIntegrator::initialize(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy)
{
   // tbox::pout<<"QuatIntegrator::initialize()"<<endl;

   initializeNonPeriodicBC();

   int finest = hierarchy->getFinestLevelNumber();

   resetIntegrator(hierarchy, 0, finest);

   if (d_evolve_quat) {
      if (d_precond_has_dquatdphi) {
         d_diffusion4quatderiv.reset(new DerivDiffusionCoeffForQuat(
             d_grid_geometry, d_H_parameter, d_quat_grad_floor, d_qlen,
             d_orient_interp_func_type1, d_avg_func_type,
             d_quat_smooth_floor_type, d_quat_diffusion_deriv_id));
         d_diffusion4quatderiv->setup(hierarchy);
      }
   }

   if (d_temperature_strategy) d_temperature_strategy->initialize(hierarchy);

   initializeSolvers(hierarchy);

   if (d_sundials_solver)
      d_sundials_solver->setInitialConditionVector(
          (solv::Sundials_SAMRAIVector*)
              solv::Sundials_SAMRAIVector::createSundialsVector(
                  d_solution_vec));
   else
      d_arkode_solver->setInitialConditionVector(
          (solv::Sundials_SAMRAIVector*)
              solv::Sundials_SAMRAIVector::createSundialsVector(
                  d_solution_vec));

   if (d_model_parameters.use_FolchPlapp())
      d_phase_rhs_strategy.reset(new FolchPlappRHSStrategy(
          d_model_parameters, d_phase_scratch_id, d_conc_scratch_id,
          d_temperature_scratch_id, d_f_l_id, d_f_a_id, d_f_b_id,
          d_phase_mobility_id, d_flux_id, this, d_free_energy_strategy,
          d_grid_geometry, d_phase_flux_strategy));
   else if (d_with_concentration) {
      if (d_model_parameters.isConcentrationModelWangSintering()) {
         tbox::plog << "Use SinteringUWangRHSStrategy..." << std::endl;
         d_phase_rhs_strategy.reset(new SinteringUWangRHSStrategy(
             d_model_parameters, d_phase_scratch_id, d_conc_scratch_id,
             d_temperature_scratch_id, d_phase_mobility_id, d_flux_id, this,
             d_grid_geometry, d_phase_flux_strategy));
      }
   }

   if (!d_phase_rhs_strategy) {
      if (d_model_parameters.norderp() > 1) {
         tbox::plog << "Use MultiOrderRHSStrategy..." << std::endl;
         d_phase_rhs_strategy.reset(new MultiOrderRHSStrategy(
             d_model_parameters, d_phase_scratch_id, d_conc_scratch_id,
             d_quat_scratch_id, d_temperature_scratch_id, d_f_l_id, d_f_a_id,
             d_f_b_id, d_phase_mobility_id, d_flux_id, d_phase_rhs_visit_id,
             this, d_free_energy_strategy, d_grid_geometry,
             d_phase_flux_strategy));
      } else {
         d_phase_rhs_strategy.reset(new PhaseRHSStrategyWithQ(
             d_model_parameters, d_phase_scratch_id, d_conc_scratch_id,
             d_quat_scratch_id, d_temperature_scratch_id, d_eta_scratch_id,
             d_f_l_id, d_f_a_id, d_f_b_id, d_phase_mobility_id, d_flux_id,
             d_quat_grad_modulus_id, d_noise_id, d_phase_rhs_visit_id, this,
             d_free_energy_strategy, d_grid_geometry, d_phase_flux_strategy));
      }
   }

   d_phase_rhs_strategy->setup(hierarchy);

   if (d_with_unsteady_temperature)
      d_temperature_rhs_strategy.reset(
          new SimpleTemperatureRHSStrategy(d_thermal_diffusivity, d_latent_heat,
                                           d_temperature_scratch_id, d_cp_id));

   if (d_model_parameters.inMovingFrame()) {
      d_movingframe_phi.reset(
          new MovingFrameRHS(d_phase_scratch_id,
                             d_model_parameters.useUpwindScheme()));
      d_adapt_moving_frame.reset(new AdaptMovingFrame(d_grid_geometry,
                                                      d_quat_model,
                                                      d_frame_velocity));
      if (d_with_concentration) {
         d_movingframe_conc.reset(
             new MovingFrameRHS(d_conc_scratch_id,
                                d_model_parameters.useUpwindScheme()));
      }
      if (d_with_unsteady_temperature) {
         d_movingframe_temp.reset(
             new MovingFrameRHS(d_temperature_scratch_id,
                                d_model_parameters.useUpwindScheme()));
      }
   }

   if (d_model_parameters.withRBmotion()) {
      d_rb_motion.reset(
          new RigidBodyMotionRHS(d_phase_scratch_id,
                                 d_model_parameters.norderpA() +
                                     d_model_parameters.norderpB(),
                                 d_weight_id, d_model_parameters.rbMobility()));
      d_rb_motion_conc.reset(new RigidBodyMotionConcRHS(
          d_phase_scratch_id,
          d_model_parameters.norderpA() + d_model_parameters.norderpB(),
          d_weight_id, d_model_parameters.rbMobility()));
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::updateSolution(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int coarsest_level, const int finest_level)
{
   if (!d_use_warm_start) {
      resetIntegrator(hierarchy, coarsest_level, finest_level);
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::updateDependentVariables(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const std::shared_ptr<hier::VariableContext> src_context,
    const std::shared_ptr<hier::VariableContext> dst_context)
{
   // tbox::pout<<"QuatIntegrator::updateDependentVariables()..."<<
   // std::endl;
   int finest_hiera_level = hierarchy->getFinestLevelNumber();

   for (int ln = 0; ln <= finest_hiera_level; ln++) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln));

      updateDependentVariables(level, src_context, dst_context);
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::updateDependentVariables(
    const std::shared_ptr<hier::PatchLevel> level,
    const std::shared_ptr<hier::VariableContext> src_context,
    const std::shared_ptr<hier::VariableContext> dst_context)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(level);
   TBOX_ASSERT(src_context);
   TBOX_ASSERT(dst_context);
#endif

   for (hier::PatchLevel::iterator ip(level->begin()); ip != level->end();
        ++ip) {
      const std::shared_ptr<hier::Patch>& patch = *ip;

      std::shared_ptr<hier::PatchData> src_data(
          patch->getPatchData(d_temperature_var, src_context));
      std::shared_ptr<hier::PatchData> dst_data(
          patch->getPatchData(d_temperature_var, dst_context));

      dst_data->copy(*src_data);
   }
}

//-----------------------------------------------------------------------

double QuatIntegrator::Advance(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   //   tbox::pout<<"QuatIntegrator::advance()"<<endl;

   t_advance_timer->start();

   // Take one step
   int return_code = d_sundials_solver ? d_sundials_solver->solve()
                                       : d_arkode_solver->solve();

   updateDependentVariables(hierarchy, d_scratch, d_current);

   // Check the return code
   if (return_code != 0) {
      tbox::pout << "   SUNDIALS solver return code          " << return_code
                 << ": ";

      switch (return_code) {
         case -1:
            tbox::pout << "The cvode_mem argument was nullptr" << std::endl;
            break;
         case -2:
            tbox::pout << "One of the inputs to the integrator is illegal"
                       << std::endl;
            break;
         case -3:
            tbox::pout << "The solver took maxstep internal steps but "
                       << "could not reach t_f" << std::endl;
            break;
         case -4:
            tbox::pout << "The solver could not satisfy the accuracy "
                          "demanded "
                       << "by the user for some internal step" << std::endl;
            break;
         case -5:
            tbox::pout << "Error test failures occurred too many times "
                       << "during one internal time step or occurred"
                       << "with |h| = hmin" << std::endl;
            break;
         case -6:
            tbox::pout << "Convergence test failures occurred too many "
                          "times "
                       << "during one internal time step or occurred "
                       << "with |h| = hmin" << std::endl;
            break;
         case -7:
            tbox::pout << "The linear solver's setup routine failed in an "
                       << "unrecoverable manner" << std::endl;
            break;
         case -8:
            tbox::pout << "The linear solver's solve routine failed in an "
                       << "unrecoverable manner" << std::endl;
            break;
         case -22:
            tbox::pout << "One of the function inputs is illegal." << std::endl;
            break;
         default: tbox::pout << "Unknown return code" << std::endl;
      }

      TBOX_ERROR("Integrator failure\n");
   }

   // Collect and (if requested) print the integrator statistics
   double time = currentTime();
   double dt = currentDT();
   int newton_iter =
       d_sundials_solver ? d_sundials_solver->getNumberOfNewtonIterations() : 0;
   int lin_iter =
       d_sundials_solver ? d_sundials_solver->getNumberOfLinearIterations() : 0;
   int newton_fail =
       d_sundials_solver
           ? d_sundials_solver->getNumberOfNonlinearConvergenceFailures()
           : 0;
   int lin_fail =
       d_sundials_solver
           ? d_sundials_solver->getNumberOfLinearConvergenceFailures()
           : 0;
   int err_test_fail =
       d_sundials_solver
           ? d_sundials_solver->getNumberOfLocalErrorTestFailures()
           : 0;
   int order = d_sundials_solver
                   ? d_sundials_solver->getOrderUsedDuringLastInternalStep()
                   : 4;
   int f_eval = numberRHSeval();
   int p_setup = d_sundials_solver
                     ? d_sundials_solver->getNumberOfPreconditionerEvaluations()
                     : 0;
   int p_apply = d_sundials_solver
                     ? d_sundials_solver->getNumberOfPrecondSolveCalls()
                     : 0;

   if (d_show_integrator_stats) {
      tbox::pout << "Integrator statistics:"
                 << "\n   Time step            " << dt
                 << "\n   Newton iterations    "
                 << newton_iter - d_cum_newton_iter
                 << "\n   Newton failures      "
                 << newton_fail - d_cum_newton_fail
                 << "\n   Linear iterations    " << lin_iter - d_cum_lin_iter
                 << "\n   Linear failures      " << lin_fail - d_cum_lin_fail
                 << "\n   Error test failures  "
                 << err_test_fail - d_cum_err_test_fail
                 << "\n   Integration order    " << order
                 << "\n   Function evaluations " << f_eval - d_cum_f_eval
                 << "\n   Precond setups       " << p_setup - d_cum_p_setup
                 << "\n   Precond applications " << p_apply - d_cum_p_apply;
      tbox::pout << std::endl;
   }

   // Update the timestamp of the current solution
   updateSolutionTime(hierarchy, time);

   // Coarsen updated data to have consistent data on all levels when
   // regridding
   coarsenData(d_phase_id, d_eta_id, d_quat_id, d_conc_id, d_temperature_id,
               hierarchy);

   if (d_model_parameters.adaptMovingFrame()) {
      static int astep = 0;
      if (astep % 25 == 0) {
         d_frame_velocity =
             d_adapt_moving_frame->adaptVelocity(time, hierarchy,
                                                 d_frame_velocity,
                                                 d_phase_scratch_id);

         std::shared_ptr<GradientTemperatureStrategy> temperature_strategy =
             std::dynamic_pointer_cast<GradientTemperatureStrategy>(
                 d_temperature_strategy);
         if (temperature_strategy) {
            temperature_strategy->resetVelocity(time, d_frame_velocity);
            tbox::pout << "resetVelocity to " << d_frame_velocity << std::endl;
         }
      }
      astep++;
   }

   // Save the cumulative counters
   d_previous_timestep = dt;
   d_cum_newton_iter = newton_iter;
   d_cum_lin_iter = lin_iter;
   d_cum_newton_fail = newton_fail;
   d_cum_lin_fail = lin_fail;
   d_cum_err_test_fail = err_test_fail;
   d_cum_f_eval = f_eval;
   d_cum_p_setup = p_setup;
   d_cum_p_apply = p_apply;

   t_advance_timer->stop();

   return dt;
}

//-----------------------------------------------------------------------

void QuatIntegrator::setVerbosity(const int v)
{
   if (v == 0) {
      d_show_integrator_stats = false;
      d_show_phase_sys_stats = false;
      d_show_eta_sys_stats = false;
      d_show_conc_sys_stats = false;
      d_show_quat_sys_stats = false;
      d_show_temperature_sys_stats = false;
   }
   if (v > 3) {
      d_show_integrator_stats = true;
   }
   if (v > 4) {
      d_show_phase_sys_stats = true;
      d_show_eta_sys_stats = true;
      d_show_conc_sys_stats = true;
      d_show_quat_sys_stats = true;
      d_show_temperature_sys_stats = true;
   }
}

//-----------------------------------------------------------------------

// If this routine were to be called at any time other than at
// initialization or reset of the integrator, this would potentially
// not be adequate.  Rather, a reset of the integrator may be needed.

void QuatIntegrator::setTimestep(const double dt)
{
   d_previous_timestep = dt;
   if (d_sundials_solver) {
      d_sundials_solver->setInitialStepSize(d_previous_timestep);
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::printSolverTotals(void)
{
   tbox::pout << "Integrator Totals:"
              << "\n   Newton iterations    " << d_cum_newton_iter
              << "\n   Newton failures      " << d_cum_newton_fail
              << "\n   Linear iterations    " << d_cum_lin_iter
              << "\n   Linear failures      " << d_cum_lin_fail
              << "\n   Error test failures  " << d_cum_err_test_fail
              << "\n   Function evaluations " << d_cum_f_eval
              << "\n   Precond setups       " << d_cum_p_setup
              << "\n   Precond applications " << d_cum_p_apply;
   tbox::pout << std::endl;
}

//-----------------------------------------------------------------------

void QuatIntegrator::coarsenData(
    const int phase_id, const int eta_id, const int quat_id, const int conc_id,
    const int temperature_id,
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   xfer::CoarsenAlgorithm coarsen_alg(tbox::Dimension(NDIM));
   if (d_with_phase) {
      coarsen_alg.registerCoarsen(phase_id, phase_id, d_phase_coarsen_op);
   }

   if (d_with_third_phase) {
      coarsen_alg.registerCoarsen(eta_id, eta_id, d_eta_coarsen_op);
   }

   if (d_evolve_quat) {
      coarsen_alg.registerCoarsen(quat_id, quat_id, d_quat_coarsen_op);
   }

   if (d_with_unsteady_temperature) {
      coarsen_alg.registerCoarsen(temperature_id, temperature_id,
                                  d_temperature_coarsen_op);
   }

   if (d_with_concentration) {
      std::shared_ptr<hier::CoarsenOperator> conc_coarsen_op =
          d_grid_geometry->lookupCoarsenOperator(d_conc_var,
                                                 "CONSERVATIVE_COARSEN");

      coarsen_alg.registerCoarsen(conc_id, conc_id, conc_coarsen_op);
   }

   for (int amr_level = hierarchy->getFinestLevelNumber(); amr_level > 0;
        amr_level--) {
      std::shared_ptr<hier::PatchLevel> fine_level =
          hierarchy->getPatchLevel(amr_level);
      std::shared_ptr<hier::PatchLevel> coarse_level =
          hierarchy->getPatchLevel(amr_level - 1);

      std::shared_ptr<xfer::CoarsenSchedule> schedule =
          coarsen_alg.createSchedule(coarse_level, fine_level, nullptr);

      schedule->coarsenData();
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::updateSolutionTime(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time)
{
   //   tbox::pout<<"QuatIntegrator::updateSolutionTime()"<<endl;

   /*
     Update time to the new solution time.
   */

   d_current_time = time;

   if (d_with_phase) {
      updateTimeForPatchID(hierarchy, d_phase_id);
   }

   if (d_with_third_phase) {
      updateTimeForPatchID(hierarchy, d_eta_id);
   }

   if (d_with_orientation) {
      updateTimeForPatchID(hierarchy, d_quat_id);
   }

   if (d_with_concentration) {
      updateTimeForPatchID(hierarchy, d_conc_id);
   }

   if (d_with_unsteady_temperature) {
      updateTimeForPatchID(hierarchy, d_temperature_id);
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::updateTimeForPatchID(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const int id)
{
   assert(id >= 0);

   int nlevels = hierarchy->getNumberOfLevels();

   for (int amr_level = 0; amr_level < nlevels; amr_level++) {

      std::shared_ptr<hier::PatchLevel> patch_level =
          hierarchy->getPatchLevel(amr_level);

      patch_level->setTime(d_current_time, id);
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::setDiffusionCoeffForConcentration(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time)
{
   assert(hierarchy);

   t_set_diffcoeff_conc_timer->start();
   // tbox::plog << "QuatIntegrator::setDiffusionCoeffForConcentration"
   //           << std::endl;

   std::shared_ptr<EBSCompositionRHSStrategy> ebs_rhs =
       std::dynamic_pointer_cast<EBSCompositionRHSStrategy>(
           d_composition_rhs_strategy);
   if (ebs_rhs) {
      ebs_rhs->setDiffusionCoeffForPreconditioner(hierarchy);
      // set diffusion coefficients which is a function of the free energy
      // form
      assert(d_composition_diffusion_strategy);
      d_composition_diffusion_strategy->setDiffusion(hierarchy,
                                                     d_temperature_scratch_id,
                                                     d_phase_scratch_id);
   }

   std::shared_ptr<KKSCompositionRHSStrategy> kks_rhs =
       std::dynamic_pointer_cast<KKSCompositionRHSStrategy>(
           d_composition_rhs_strategy);
   if (kks_rhs) {
      kks_rhs->setDiffusionCoeff(hierarchy, time);
   }

   std::shared_ptr<WangSinteringCompositionRHSStrategy> ws_rhs =
       std::dynamic_pointer_cast<WangSinteringCompositionRHSStrategy>(
           d_composition_rhs_strategy);
   if (ws_rhs) {
      ws_rhs->setDiffusionCoeff(hierarchy, time);
   }

   for (int amr_level = hierarchy->getFinestLevelNumber() - 1; amr_level >= 0;
        amr_level--) {
      d_conc_diffusion_coarsen_schedule[amr_level]->coarsenData();
   }

   d_quat_model->setDiffusionVisit(hierarchy);

   t_set_diffcoeff_conc_timer->stop();
}

//-----------------------------------------------------------------------
// set grad q at side to constant value
// ( s.t. norm2(grad q)=1 ) to avoid problems with
// initial data that may not be smooth enough
// grad q at side is later used by linear solver in denominator
// of diffusion coefficient
// Note: value for diffusion may be problem dependent and may need to be
// tuned jlf
void QuatIntegrator::setUniformDiffusionCoeffForQuat(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   // tbox::pout<<"QuatIntegrator::setUniformDiffusionCoeffForQuat"<<endl;

   assert(hierarchy);

   const int maxl = hierarchy->getNumberOfLevels();

   // set value of alpha s.t. 1./|grad q| is 10.
   const double alpha = 0.1 / sqrt((double)(d_qlen * NDIM));
   // tbox::pout<<"alpha="<<alpha<<endl;

   std::shared_ptr<hier::PatchLevel> patch_level =
       hierarchy->getPatchLevel(maxl - 1);
   std::shared_ptr<hier::Patch> patch = *(patch_level->begin());

   // evaluate vel at finest level
   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch->getPatchGeometry()));
   const double* dx = patch_geom->getDx();
   double vel = 1.;
   for (int i = 0; i < NDIM; i++)
      vel *= dx[i];
   // tbox::pout<<"vel="<<vel<<endl;

   // vel/max(Mq)
   math::HierarchyCellDataOpsReal<double> mathops(hierarchy);
   const double maxmobility = mathops.max(d_quat_mobility_id);
   const double dval = 0.2 * vel / (maxmobility);
   tbox::pout << "dval=" << dval << std::endl;

   // set diffusion coefficients
   math::HierarchySideDataOpsReal<double> smathops(hierarchy);
   smathops.setToScalar(d_quat_grad_side_copy_id, alpha);
}

//-----------------------------------------------------------------------

void QuatIntegrator::evaluatePhaseRHS(
    const double time, std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int phase_id, const int phase_rhs_id, const bool eval_flag)
{
   assert(phase_id >= 0);
   assert(phase_rhs_id >= 0);

   t_phase_rhs_timer->start();

   d_phase_rhs_strategy->evaluateRHS(time, hierarchy, phase_rhs_id, eval_flag);

   // save dphidt if needed for other purposes
   if (needDphiDt()) fillDphiDt(hierarchy, time, phase_rhs_id);

   if (d_model_parameters.withRBmotion()) {
      d_rb_forces.clear();
      const int nop =
          d_model_parameters.norderpA() + d_model_parameters.norderpB();
      for (int i = 0; i < nop; i++) {
         std::array<double, NDIM> f;
         // external force is the same on all grains
         d_model_parameters.rbExternalForce(f);
         d_rb_forces.push_back(f);
      }

      assert(d_rigid_body_forces);
      d_rigid_body_forces->evaluatePairForces(hierarchy);
      for (int i = 0; i < nop; i++) {
         // sum up forces contributions from all other particles
         for (int j = 0; j < nop; j++) {
            if (j != i) {
               std::array<double, NDIM> f;
               d_rigid_body_forces->getPairForce(i, j, f);
               for (short d = 0; d < NDIM; d++)
                  d_rb_forces[i][d] += f[d];
            }
         }
      }

      d_rb_motion->addRHS(hierarchy, phase_rhs_id, d_rb_forces);
   }

   if (d_model_parameters.inMovingFrame()) {
      assert(d_movingframe_phi);
      d_movingframe_phi->addRHS(hierarchy, phase_rhs_id, d_frame_velocity);
   }

   if (d_model_parameters.with_rhs_visit_output() && eval_flag &&
       d_free_energy_strategy) {
      assert(d_driving_force_visit_id >= 0);
      // recompute driving force just for visualization
      d_free_energy_strategy->computeDrivingForce(time, hierarchy,
                                                  d_temperature_scratch_id,
                                                  d_phase_scratch_id,
                                                  d_conc_scratch_id, d_f_l_id,
                                                  d_f_a_id, d_f_b_id,
                                                  d_driving_force_visit_id);
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   math::HierarchyCellDataOpsReal<double> cellops(hierarchy);
   double l2rhs = cellops.L2Norm(phase_rhs_id);
   assert(l2rhs == l2rhs);
#endif

   t_phase_rhs_timer->stop();
}

//-----------------------------------------------------------------------

void QuatIntegrator::evaluateTemperatureRHS(
    std::shared_ptr<hier::PatchHierarchy> hierarchy, const int temperature_id,
    const int phase_rhs_id, const int temperature_rhs_id, const bool visit_flag)
{
   // tbox::pout << "QuatIntegrator::evaluateTemperatureRHS()..." <<
   // std::endl;
   assert(temperature_id >= 0);
   assert(temperature_rhs_id >= 0);
   assert(d_temperature_rhs_strategy);

   d_heat_capacity_strategy->setCurrentValue(hierarchy);
   d_temperature_rhs_strategy->evaluateRHS(hierarchy, temperature_rhs_id,
                                           phase_rhs_id);

   // add component related to moving frame if moving velocity!=0
   if (d_model_parameters.inMovingFrame())
      d_movingframe_temp->addRHS(hierarchy, temperature_rhs_id,
                                 d_frame_velocity);

   if (d_model_parameters.with_rhs_visit_output() && visit_flag) {
      assert(d_temperature_rhs_visit_id >= 0);
      math::HierarchyCellDataOpsReal<double> cellops(hierarchy);
      cellops.copyData(d_temperature_rhs_visit_id, temperature_rhs_id, false);
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::evaluateConcentrationRHS(
    std::shared_ptr<hier::PatchHierarchy> hierarchy, const int phase_id,
    const int conc_id, const int conc_rhs_id, const int temperature_id,
    const bool visit_flag)
{
   assert(conc_rhs_id >= 0);
   assert(d_conc_mobility > 0.);
   assert(temperature_id >= 0);

   t_conc_rhs_timer->start();

   // tbox::pout<<"QuatIntegrator::evaluateConcentrationRHS()"<<endl;
   if (d_with_antitrapping)
      assert(checkForNans(hierarchy, d_dphidt_scratch_id) == 0);

   math::PatchCellDataBasicOps<double> mathops;

   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
           ++ip) {
         std::shared_ptr<hier::Patch> patch = *ip;

         assert(d_composition_rhs_strategy);
         d_composition_rhs_strategy->computeFluxOnPatch(*patch, d_flux_conc_id);
         if (d_with_gradT)
            d_composition_rhs_strategy->addFluxFromGradTonPatch(*patch,
                                                                temperature_id,
                                                                d_flux_conc_id);
         if (d_with_antitrapping) {
            assert(phase_id >= 0);
            d_composition_rhs_strategy->addFluxFromAntitrappingonPatch(
                *patch, phase_id, d_dphidt_scratch_id, d_alpha_AT,
                d_flux_conc_id);
         }
      }

      // Coarsen flux data from next finer level so that
      // the computed flux becomes the composite grid flux.
      if (ln < hierarchy->getFinestLevelNumber()) {
         d_flux_conc_coarsen_schedule[ln]->coarsenData();
      }

      for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
           ++ip) {
         std::shared_ptr<hier::Patch> patch = *ip;

         const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
             SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                                    hier::PatchGeometry>(
                 patch->getPatchGeometry()));
         const double* dx = patch_geom->getDx();

         std::shared_ptr<pdat::CellData<double> > conc(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(conc_id)));
         assert(conc);

         std::shared_ptr<pdat::SideData<double> > flux(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(d_flux_conc_id)));
         assert(flux);
         assert(flux->getDepth() == d_ncompositions);

         std::shared_ptr<pdat::CellData<double> > conc_rhs(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(conc_rhs_id)));
         assert(conc_rhs);
         assert(conc_rhs->getDepth() == d_ncompositions);
         assert(conc_rhs->getGhostCellWidth() ==
                hier::IntVector(tbox::Dimension(NDIM), 0));

         const hier::Box& pbox = patch->getBox();
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast = pbox.upper();

         // now compute r.h.s., one species at a time
         for (int ic = 0; ic < d_ncompositions; ic++) {
            COMPUTERHSCONCENTRATION(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                                    ifirst(2), ilast(2),
#endif
                                    dx, flux->getPointer(0, ic),
                                    flux->getPointer(1, ic),
#if (NDIM == 3)
                                    flux->getPointer(2, ic),
#endif
                                    flux->getGhostCellWidth()[0],
                                    d_conc_mobility, conc_rhs->getPointer(ic),
                                    0);
         }
      }
   }

   if (d_model_parameters.withRBmotion()) {
      // d_rb_forces set when computing r.h.s. for phase variable
      d_rb_motion_conc->addRHS(hierarchy, conc_rhs_id, d_rb_forces);
   }

   // add component related to moving frame if moving velocity!=0
   if (d_model_parameters.inMovingFrame())
      d_movingframe_conc->addRHS(hierarchy, conc_rhs_id, d_frame_velocity);

   if (d_model_parameters.with_rhs_visit_output() && visit_flag) {
      assert(d_conc_rhs_visit_id >= 0);
      math::HierarchyCellDataOpsReal<double> cellops(hierarchy);
      cellops.copyData(d_conc_rhs_visit_id, conc_rhs_id, false);
   }

   t_conc_rhs_timer->stop();
}

//-----------------------------------------------------------------------

void QuatIntegrator::evaluateQuatRHS(
    std::shared_ptr<hier::PatchHierarchy> hierarchy, const int quat_rhs_id,
    const bool visit_flag)
{
   assert(d_quat_sys_solver);
   assert(d_quat_grad_side_id >= 0);

   // tbox::plog<<"QuatIntegrator::evaluateQuatRHS()"<<std::endl;
   int quat_symm_rotation_id =
       d_use_gradq_for_flux ? d_quat_symm_rotation_id : -1;

   // compute RHS using the gradient of q at sides (d_quat_grad_side_id)
   // computed with physical BC
   d_quat_sys_solver->evaluateRHS(d_phase_scratch_id, d_temperature_scratch_id,
                                  d_quat_grad_side_id, d_quat_grad_side_copy_id,
                                  quat_symm_rotation_id, d_quat_mobility_id,
                                  d_quat_scratch_id, quat_rhs_id, true);


   if (visit_flag && d_model_parameters.with_rhs_visit_output())
      for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++) {
         std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

         for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
              ++ip) {
            std::shared_ptr<hier::Patch> patch = *ip;
            const hier::Box box(patch->getBox());

            std::shared_ptr<pdat::CellData<double> > rhs(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(quat_rhs_id)));
            std::shared_ptr<pdat::CellData<double> > nrhs(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_q_rhs1_visit_id)));

            pdat::CellIterator cend(pdat::CellGeometry::end(box));
            for (pdat::CellIterator c(pdat::CellGeometry::begin(box));
                 c != cend; ++c) {
               pdat::CellIndex cell = *c;
               (*nrhs)(cell) = 0.;
               for (int m = 0; m < d_qlen; m++)
                  (*nrhs)(cell) += (*rhs)(cell, m) * (*rhs)(cell, m);
               (*nrhs)(cell) = sqrt((*nrhs)(cell));
            }
         }
      }

   // correction needed if d_quat_sys_solver
   // doesn't know about symmetry!
   if (d_symmetry_aware && !d_use_gradq_for_flux) {
      // tbox::pout<<"Correct r.h.s. for symmetry..."<<endl;
      correctRhsForSymmetry(hierarchy, d_quat_scratch_id, quat_rhs_id);
   }

   if (visit_flag && d_model_parameters.with_rhs_visit_output()) {
      for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++) {
         std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

         for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
              ++ip) {
            std::shared_ptr<hier::Patch> patch = *ip;
            const hier::Box box(patch->getBox());

            std::shared_ptr<pdat::CellData<double> > rhs(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(quat_rhs_id)));
            std::shared_ptr<pdat::CellData<double> > nrhs(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_modulus_q_rhs_visit_id)));

            pdat::CellIterator cend(pdat::CellGeometry::end(box));
            for (pdat::CellIterator c(pdat::CellGeometry::begin(box));
                 c != cend; c++) {
               pdat::CellIndex cell = *c;
               (*nrhs)(cell) = 0.;
               for (int m = 0; m < d_qlen; m++) {
                  (*nrhs)(cell) += (*rhs)(cell, m) * (*rhs)(cell, m);
               }
               (*nrhs)(cell) = sqrt((*nrhs)(cell));
            }
         }
      }
      math::HierarchyCellDataOpsReal<double> cellops(hierarchy);
      cellops.copyData(d_q_rhs_visit_id, quat_rhs_id, false);
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::computeQuatGradients(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, double time,
    const bool recompute_quat_sidegrad)
{
   // tbox::pout<<"QuatIntegrator::computeQuatGradients()..."<<endl;
   t_quat_grad_timer->start();

   int diff_id = -1;
   d_quat_grad_strategy->computeDiffs(hierarchy, d_quat_scratch_id, diff_id,
                                      time, QuatGradStrategy::FORCE);
   assert(diff_id >= 0);

   // Compute cell-centered gradients
   d_quat_grad_strategy->computeGradCell(hierarchy, diff_id,
                                         d_quat_grad_cell_id, time,
                                         QuatGradStrategy::FORCE);

   // Compute gradients on cell faces
   d_quat_grad_strategy->computeGradSide(hierarchy, diff_id,
                                         d_quat_grad_side_id, time,
                                         QuatGradStrategy::FORCE);

   if (recompute_quat_sidegrad) {
      math::HierarchySideDataOpsReal<double> sideops(hierarchy);

      // update d_quat_grad_side_copy_id to be used for D_q(phi)/|nabla q|
      sideops.copyData(d_quat_grad_side_copy_id, d_quat_grad_side_id, false);
   }

   // Compute modulus of cell-centered gradients
   // to be used in phi r.h.s.
   if (d_model_parameters.quat_grad_modulus_from_cells()) {
      // tbox::pout<<"computeGradModulus()..."<<endl;
      d_quat_grad_strategy->computeGradModulus(hierarchy, d_quat_grad_cell_id,
                                               d_quat_grad_modulus_id, time,
                                               QuatGradStrategy::FORCE);

   } else {
      // tbox::pout<<"computeGradModulusFromSides()..."<<endl;
      d_quat_grad_strategy->computeGradModulusFromSides(
          hierarchy, d_quat_grad_side_id, d_quat_grad_modulus_id, time,
          QuatGradStrategy::FORCE);
   }

   t_quat_grad_timer->stop();
}

//-----------------------------------------------------------------------

void QuatIntegrator::fillScratchComposition(
    double time, std::shared_ptr<solv::SAMRAIVectorReal<double> > y,
    xfer::RefineAlgorithm& copy_to_scratch)
{
   (void)time;

   int y_conc_id = y->getComponentDescriptorIndex(d_conc_component_index);

   copy_to_scratch.registerRefine(d_conc_scratch_id,  // destination
                                  y_conc_id,          // source
                                  d_conc_scratch_id,  // temporary
                                  d_conc_refine_op);
}

void QuatIntegrator::fillDphiDt(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                                const double time, const int phase_rhs_id)
{
   assert(d_dphidt_scratch_id >= 0);
   assert(checkForNans(hierarchy, phase_rhs_id) == 0);
   if (!d_all_periodic) assert(d_dphidt_bc_helper != 0);

   xfer::RefineAlgorithm copy_to_scratch;

   copy_to_scratch.registerRefine(d_dphidt_scratch_id,  // destination
                                  phase_rhs_id,         // source
                                  d_dphidt_scratch_id,  // temporary
                                  d_phase_refine_op);

   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      // no physical BC filling performed
      copy_to_scratch
          .createSchedule(level, ln - 1, hierarchy, d_dphidt_bc_helper)
          ->fillData(time);
   }

   assert(checkForNans(hierarchy, d_dphidt_scratch_id) == 0);
}

//-----------------------------------------------------------------------

void QuatIntegrator::fillScratch(
    double time, std::shared_ptr<solv::SAMRAIVectorReal<double> > y)
{
   assert(time >= 0.);

   std::shared_ptr<hier::PatchHierarchy> hierarchy = y->getPatchHierarchy();
#ifdef DEBUG_CHECK_ASSERTIONS
   math::HierarchyCellDataOpsReal<double> mathops(hierarchy);
#endif
   // Copy the input phase/quat/conc std::vectors to the phase/quat/conc
   // scratch arrays and fill the ghost cells

   xfer::RefineAlgorithm copy_to_scratch;

   if (d_with_phase) {
      int y_phase_id = y->getComponentDescriptorIndex(d_phase_component_index);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(y_phase_id > -1);
      const double norm_y_phi = mathops.L2Norm(y_phase_id);
      assert(norm_y_phi == norm_y_phi);
#endif
      copy_to_scratch.registerRefine(
          d_phase_scratch_id,  // destination
          y_phase_id,          // source
          d_phase_scratch_id,  // temporary work space
          d_phase_refine_op);
   }

   if (d_with_third_phase) {
      const int y_eta_id =
          y->getComponentDescriptorIndex(d_eta_component_index);
      copy_to_scratch.registerRefine(d_eta_scratch_id,  // destination
                                     y_eta_id,          // source
                                     d_eta_scratch_id,  // temporary work space
                                     d_eta_refine_op);
   }

   if (d_evolve_quat) {
      const int y_quat_id =
          y->getComponentDescriptorIndex(d_quat_component_index);
      assert(y_quat_id > -1);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(y_quat_id > -1);
      const double norm_y_q = mathops.L2Norm(y_quat_id);
      assert(norm_y_q == norm_y_q);
#endif
      // tbox::pout<<"Fill scratch from "<<y_quat_id<<" into
      // "<<d_quat_scratch_id<<endl;
      copy_to_scratch.registerRefine(d_quat_scratch_id,  // destination
                                     y_quat_id,          // source
                                     d_quat_scratch_id,  // temporary work space
                                     d_quat_refine_op);
   }

   if (d_with_concentration) {
      fillScratchComposition(time, y, copy_to_scratch);
   }

   if (d_with_unsteady_temperature) {
      assert(d_temperature_scratch_id >= 0);
      int y_temp_id =
          y->getComponentDescriptorIndex(d_temperature_component_index);
      assert(d_temperature_scratch_id > -1);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert(y_temp_id > -1);
      const double norm_y_temp = mathops.L2Norm(y_temp_id);
      assert(norm_y_temp == norm_y_temp);
#endif
      copy_to_scratch.registerRefine(
          d_temperature_scratch_id,  // destination
          y_temp_id,                 // source
          d_temperature_scratch_id,  // temporary work space
          d_temperature_refine_op);
   }

   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      copy_to_scratch
          .createSchedule(level, ln - 1, hierarchy, d_all_refine_patch_strategy)
          ->fillData(time);
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::computeMobilities(
    double time, std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   if (d_with_phase) {
      assert(d_mobility_strategy);
      d_mobility_strategy->computePhaseMobility(hierarchy, d_phase_scratch_id,
                                                d_phase_mobility_id, time,
                                                QuatMobilityStrategy::FORCE);
   }
   if (d_with_third_phase) {
      d_mobility_strategy->computeEtaMobility(hierarchy, d_eta_scratch_id,
                                              d_eta_mobility_id, time,
                                              QuatMobilityStrategy::FORCE);
   }
   if (d_evolve_quat) {
      d_mobility_strategy->computeQuatMobility(hierarchy, d_phase_scratch_id,
                                               d_quat_mobility_id, time,
                                               QuatMobilityStrategy::FORCE);

      if (d_precond_has_dquatdphi) {
         d_mobility_strategy->computeQuatMobilityDeriv(
             hierarchy, d_phase_scratch_id, d_quat_mobility_deriv_id, time,
             QuatMobilityStrategy::FORCE);
      }
   }
   if (d_precond_has_dTdphi) {
      assert(d_phase_temperature_mobility_id >= 0);
      d_mobility_strategy->computePhaseTemperatureMobility(
          hierarchy, d_phase_mobility_id, d_cp_id,
          d_phase_temperature_mobility_id);
   }
}

//-----------------------------------------------------------------------

void QuatIntegrator::setCoefficients(
    double time, std::shared_ptr<solv::SAMRAIVectorReal<double> > y,
    const bool recompute_quat_sidegrad)
{
   // tbox::pout << "Entering QuatIntegrator::setCoefficients()" <<
   // std::endl;
   t_set_coeff_timer->start();

#ifdef DEBUG_CHECK_ASSERTIONS
   int n = 0;
   if (d_with_phase) n++;
   if (d_with_third_phase) n++;
   if (d_evolve_quat) n++;
   if (d_with_concentration) n++;
   if (d_with_unsteady_temperature) n++;
   assert(y->getNumberOfComponents() == n);
#endif

   std::shared_ptr<hier::PatchHierarchy> hierarchy = y->getPatchHierarchy();

   int phase_id = d_with_phase
                      ? y->getComponentDescriptorIndex(d_phase_component_index)
                      : -1;
   int eta_id = d_with_third_phase
                    ? y->getComponentDescriptorIndex(d_eta_component_index)
                    : -1;
   int quat_id = (d_evolve_quat)
                     ? y->getComponentDescriptorIndex(d_quat_component_index)
                     : -1;
   int conc_id = d_with_concentration
                     ? y->getComponentDescriptorIndex(d_conc_component_index)
                     : -1;
   int temperature_id =
       d_with_unsteady_temperature
           ? y->getComponentDescriptorIndex(d_temperature_component_index)
           : -1;
#ifdef DEBUG_CHECK_ASSERTIONS
   if (temperature_id >= 0) {
      math::HierarchyCellDataOpsReal<double> mathops(hierarchy);
      const double norm_y_temp = mathops.L2Norm(temperature_id);
      assert(norm_y_temp == norm_y_temp);
   }
   if (conc_id >= 0) {
      assert(checkForNans(hierarchy, conc_id) == 0);
   }
#endif

   coarsenData(phase_id, eta_id, quat_id, conc_id, temperature_id, hierarchy);
   fillScratch(time, y);

   if (d_evolve_quat) {
      computeQuatGradients(hierarchy, time, recompute_quat_sidegrad);
   }

   if (d_with_phase) {
      if (d_compute_velocity) {
         int phase_diffs_id = -1;
         int phase_grad_cell_id = -1;
         d_quat_model->computePhaseDiffs(hierarchy, d_phase_scratch_id,
                                         phase_diffs_id, time);
         d_quat_model->computePhaseGradCell(hierarchy, phase_diffs_id,
                                            phase_grad_cell_id, time);
      }
   }

   if (d_evolve_quat) {

      // Compute quaternion diffusion coefficient
      if (time < d_uniform_diffusion_time_threshold) {
         setUniformDiffusionCoeffForQuat(hierarchy);
      } else {
         if (d_precond_has_dquatdphi)
            d_diffusion4quatderiv->setDerivDiffusion(hierarchy,
                                                     d_phase_scratch_id,
                                                     d_temperature_scratch_id,
                                                     d_quat_grad_side_copy_id);
      }
   }

   if (d_with_concentration &&
       d_model_parameters.concentrationModelNeedsPhaseConcentrations()) {
      d_quat_model->computePhaseConcentrations(hierarchy);
   }

   // mobilities may depend on cl and cs, thus they should be computed after
   // cl and cs
   computeMobilities(time, hierarchy);

   t_set_coeff_timer->stop();
}

//-----------------------------------------------------------------------
// Virtual function from CVODEAbstractFunction
// jlf (7/11): fd_flag is set to 1 when the function is called to evaluate
// F(u+sigma*v) and calculate a finite difference with F(u) (call with
// fd_flag=0)

int QuatIntegrator::evaluateRHSFunction(double time,
                                        solv::SundialsAbstractVector* y,
                                        solv::SundialsAbstractVector* y_dot,
                                        int fd_flag)
{
   if (d_with_unsteady_temperature) assert(d_temperature_sys_solver);

   // tbox::pout << "Entering QuatIntegrator::evaluateRHSFunction" <<
   // std::endl; tbox::pout << "QuatIntegrator::evaluateRHSFunction with
   // fd_flag="<<fd_flag
   // << std::endl;

   t_rhs_timer->start();

   // Convert the Sundials std::vectors to SAMRAI std::vectors
   std::shared_ptr<solv::SAMRAIVectorReal<double> > y_samvect =
       solv::Sundials_SAMRAIVector::getSAMRAIVector(y);
   std::shared_ptr<solv::SAMRAIVectorReal<double> > y_dot_samvect =
       solv::Sundials_SAMRAIVector::getSAMRAIVector(y_dot);

   std::shared_ptr<hier::PatchHierarchy> hierarchy =
       y_samvect->getPatchHierarchy();

   // #ifdef DEBUG_CHECK_ASSERTIONS
   int n = 0;
   if (d_with_phase) n++;
   if (d_with_third_phase) n++;
   if (d_evolve_quat) n++;
   if (d_with_concentration) n++;
   if (d_with_unsteady_temperature) n++;
   assert(y_dot_samvect->getNumberOfComponents() == n);
   // #endif

#ifdef DEBUG_CHECK_ASSERTIONS
   int temperature_id = d_with_unsteady_temperature
                            ? y_samvect->getComponentDescriptorIndex(
                                  d_temperature_component_index)
                            : -1;
   if (temperature_id >= 0) {
      math::HierarchyCellDataOpsReal<double> mathops(hierarchy);
      const double norm_y_temp = mathops.L2Norm(temperature_id);
      assert(norm_y_temp == norm_y_temp);
   }
   if (d_with_concentration) {
      int conc_id =
          y_samvect->getComponentDescriptorIndex(d_conc_component_index);
      assert(checkForNans(hierarchy, conc_id) == 0);
   }
#endif

   /*
      If fd_flag != 0, the integrator is calling this function to compute a
      finite difference approximation of the system Jacobian. In this case,
      if d_lag_quat_sidegrad is true, we lag the computation of the
      quaternion side gradients.
   */
   const bool recompute_quat_sidegrad = (fd_flag == 0) || !d_lag_quat_sidegrad;

   setTemperatureField(hierarchy, time);
#ifdef DEBUG_CHECK_ASSERTIONS
   if (temperature_id >= 0) {
      assert(checkForNans(hierarchy, temperature_id) == 0);
   }
#endif

   setCoefficients(time, y_samvect, recompute_quat_sidegrad);

   // Set the phase component of the RHS
   int ydot_phase_id = -1;
   if (d_with_phase) {
      ydot_phase_id =
          y_dot_samvect->getComponentDescriptorIndex(d_phase_component_index);

      bool need_iterate = false;  // is a fixed point iteration needed?
      math::HierarchyCellDataOpsReal<double> cellops(hierarchy);
      do {
         // tbox::pout<<"Evaluate phase rhs..."<<endl;
         evaluatePhaseRHS(time, hierarchy, y_dot_samvect, fd_flag);

         if (d_with_partition_coeff &&
             d_model_parameters.with_Aziz_partition_coeff()) {

            assert(d_tmp1_id >= 0);
            assert(d_tmp2_id >= 0);
            assert(d_velocity_id >= 0);

            need_iterate = true;  // c_l, c_s depend on ydot_phase_id, which
                                  // depend on c_l, c_s

            cellops.copyData(d_tmp1_id, d_velocity_id, false);
         }

         if (d_compute_velocity && (fd_flag == 0)) {
            d_quat_model->computeVelocity(hierarchy, ydot_phase_id);
         }

         if (d_with_partition_coeff &&
             d_model_parameters.with_Aziz_partition_coeff()) {
            assert(d_velocity_id >= 0);
            assert(d_tmp1_id);
            assert(d_tmp2_id);

            cellops.subtract(d_tmp2_id, d_tmp1_id, d_velocity_id);

            double norm_diff = cellops.L1Norm(d_tmp2_id);
            if (norm_diff < 1.e-8)
               need_iterate = false;
            else
               tbox::pout << "Norm diff = " << norm_diff << std::endl;

            assert(d_partition_coeff_strategy != nullptr);

            // compute phase concentrations again if they depend on velocity
            // tbox::pout<<"Evaluate c_L, c_S..."<<endl;
            d_quat_model->computePhaseConcentrations(hierarchy);
         }

      } while (need_iterate);
   }


   // math::HierarchyCellDataOpsReal<double> mathops(hierarchy);
   // const double norm_ydot_phase_id = mathops.L2Norm( ydot_phase_id );
   // tbox::plog<<"L2 Norm ydot_phase_id="<<norm_ydot_phase_id<<endl;

   // Set the quaternion component of the RHS
   if (d_evolve_quat) {

      evaluateQuatRHS(hierarchy, y_dot_samvect, fd_flag);
   }

   // Set the concentration component of the RHS
   if (d_with_concentration) {
      assert(y_dot_samvect->getNumberOfComponents() > 0);

      if (recompute_quat_sidegrad)
         setDiffusionCoeffForConcentration(hierarchy, time);

      const int ydot_conc_id =
          y_dot_samvect->getComponentDescriptorIndex(d_conc_component_index);
      evaluateConcentrationRHS(hierarchy, d_phase_scratch_id, d_conc_scratch_id,
                               ydot_conc_id, d_temperature_scratch_id,
                               fd_flag == 0);

      assert(checkForNans(hierarchy, ydot_conc_id) == 0);
   }

   // Set the temperature component of the RHS
   if (d_with_unsteady_temperature) {

      evaluateTemperatureRHS(hierarchy, y_dot_samvect, fd_flag);

#ifdef DEBUG_CHECK_ASSERTIONS
      int Tdot_id = y_dot_samvect->getComponentDescriptorIndex(
          d_temperature_component_index);
      assert(checkForNans(hierarchy, Tdot_id) == 0);
#endif
   }

   t_rhs_timer->stop();

   return 0;  // Always successful
}

//-----------------------------------------------------------------------
// Virtual function from CPODESAbstractFunction or CVODEAbstractFunction

int QuatIntegrator::CVSpgmrPrecondSet(double t, solv::SundialsAbstractVector* y,
                                      solv::SundialsAbstractVector* fy, int jok,
                                      int* jcurPtr, double gamma)
{
   (void)fy;
   (void)jok;

   // tbox::pout << "QuatIntegrator::CVSpgmrPrecondSet, jok = " << jok <<
   // std::endl;
   t_psolve_setup_timer->start();

   // Ignore jok flag
   // If (jok == TRUE), Jacobian-related data saved
   // from previous call could be reused with new gamma

   // Convert passed-in std::vector into SAMRAI std::vector
   std::shared_ptr<solv::SAMRAIVectorReal<double> > y_samvect =
       solv::Sundials_SAMRAIVector::getSAMRAIVector(y);

   setCoefficients(t, y_samvect, true);

   if (d_with_phase && d_phase_sys_solver) {
      const double phase_well_scale =
          d_model_parameters.norderp() > 1 ? 0. : d_phase_well_scale;
      d_phase_sys_solver->setOperatorCoefficients(d_phase_scratch_id,
                                                  d_phase_mobility_id,
                                                  d_epsilon_phase, gamma,
                                                  phase_well_scale,
                                                  d_phase_well_func_type);
   }

   if (d_with_third_phase && d_eta_sys_solver) {
      d_eta_sys_solver->setOperatorCoefficients(
          d_phase_scratch_id, d_eta_scratch_id, d_eta_mobility_id,
          d_epsilon_eta, gamma, d_energy_interp_func_type, d_eta_well_scale,
          d_eta_well_func_type);
   }

   if (d_with_unsteady_temperature)
      if (d_temperature_sys_solver) {
         const double m = 1.;
         const double d = (-gamma * d_thermal_diffusivity);
         const double c = 1.;

         d_temperature_sys_solver->setOperatorCoefficients(m, c, d);

         if (d_precond_has_dTdphi) {
            TBOX_ASSERT(d_phase_temperature_fac_ops);
            d_phase_temperature_fac_ops->setOperatorCoefficients(
                d_phase_scratch_id, d_phase_temperature_mobility_id,
                d_epsilon_phase, d_latent_heat, d_phase_well_scale,
                d_phase_well_func_type);
         }
      }

   if (d_with_concentration) setCompositionOperatorCoefficients(gamma);

   if (d_evolve_quat) {
      assert(d_quat_sys_solver);

      d_quat_sys_solver->setOperatorCoefficients(gamma, d_quat_mobility_id,
                                                 d_quat_mobility_deriv_id,
                                                 d_phase_scratch_id,
                                                 d_temperature_scratch_id,
                                                 d_quat_diffusion_deriv_id,
                                                 d_quat_grad_side_copy_id,
                                                 d_quat_scratch_id);
   }

   // Tell the integrator that the Jacobian data was recomputed
   *jcurPtr = TRUE;

   t_psolve_setup_timer->stop();

   return 0;
}

//-----------------------------------------------------------------------
void QuatIntegrator::setCompositionOperatorCoefficients(const double gamma)
{
   TBOX_ASSERT(d_conc_mobility > 0.);

   if (d_with_concentration && d_conc_sys_solver) {
      // Set concentration block coefficients
      d_conc_sys_solver->setOperatorCoefficients(gamma, d_conc_pfm_diffusion_id,
                                                 d_conc_mobility);
   }
}

//-----------------------------------------------------------------------
// returns 0 if converged
int QuatIntegrator::PhasePrecondSolve(
    std::shared_ptr<hier::PatchHierarchy> hierarchy, int r_phase_id,
    int z_phase_id, const double delta, const double gamma)
{
   t_phase_precond_timer->start();

   if (d_show_phase_sys_stats) {
      tbox::plog << "Preconditioner for Phase block with tol " << delta
                 << std::endl;
   }

   math::HierarchyCellDataOpsReal<double> cellops(hierarchy);

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(cellops.max(r_phase_id) == cellops.max(r_phase_id));
#endif

   if (d_precond_has_dPhidT) {
      assert(d_model_parameters.norderp() == 1);
      std::shared_ptr<DeltaTemperatureFreeEnergyStrategy> free_energy_strategy =
          std::dynamic_pointer_cast<DeltaTemperatureFreeEnergyStrategy>(
              d_free_energy_strategy);
      TBOX_ASSERT(free_energy_strategy);
      free_energy_strategy->applydPhidTBlock(
          hierarchy, d_temperature_sol_id, z_phase_id, d_phase_rhs_id,
          d_model_parameters.phase_mobility());

      // Add -1.*gamma times the just computed product to the right-hand side
      cellops.axpy(d_phase_rhs_id, -1. * gamma, d_phase_rhs_id, r_phase_id,
                   false);
   }

   // Set the tolerance for the FAC solve as requested by integrator
   d_phase_sys_solver->setResidualTolerance(delta);

   int nconverged = 0;
   // solve for one depth at a time
   for (int depth = 0; depth < d_model_parameters.norderp(); depth++) {
      if (!d_precond_has_dPhidT)
         copyDepthCellData(hierarchy, d_phase_rhs_id, 0, r_phase_id, depth);

      // Zero out the initial guess in the temporary solution array
      cellops.setToScalar(d_phase_sol_id, 0., false);

      // Solve the phase block system
      bool converged =
          d_phase_sys_solver->solveSystem(d_phase_sol_id, d_phase_rhs_id);
      if (converged) nconverged++;

      // Copy solution from the local temporary to the output array
      // (z_phase_id)
      copyDepthCellData(hierarchy, z_phase_id, depth, d_phase_sol_id, 0);
   }
   int retcode = (nconverged == d_model_parameters.norderp()) ? 0 : 1;

   t_phase_precond_timer->stop();

   return retcode;
}

//-----------------------------------------------------------------------
int QuatIntegrator::EtaPrecondSolve(
    std::shared_ptr<hier::PatchHierarchy> hierarchy, int r_eta_id, int z_eta_id,
    const double delta)
{
   if (d_show_eta_sys_stats) {
      tbox::pout << "Preconditioner for Eta block" << std::endl;
   }

   math::HierarchyCellDataOpsReal<double> cellops(hierarchy);

   // Zero out the initial guess in the temporary solution array
   cellops.setToScalar(d_eta_sol_id, 0., false);

   // Copy the right-hand side to the temporary right-hand side array
   cellops.copyData(d_eta_rhs_id, r_eta_id, false);

   // Set the tolerance for the FAC solve as requested by integrator
   d_eta_sys_solver->setResidualTolerance(delta);

   // Solve the eta block system
   bool converged = d_eta_sys_solver->solveSystem(d_eta_sol_id, d_eta_rhs_id);

   int retcode = converged ? 0 : 1;

   // Copy solution from the local temporary to the output array
   cellops.copyData(z_eta_id, d_eta_sol_id, false);

   return retcode;
}

//-----------------------------------------------------------------------
int QuatIntegrator::TemperaturePrecondSolve(
    std::shared_ptr<hier::PatchHierarchy> hierarchy, int r_temperature_id,
    int z_temperature_id, const double delta, const double gamma)
{
   if (d_show_temperature_sys_stats) {
      tbox::plog << "Preconditioner for temperature block with tol " << delta
                 << std::endl;
   }

   math::HierarchyCellDataOpsReal<double> cellops(hierarchy);

   if (d_precond_has_dTdphi) {
      TBOX_ASSERT(d_phase_temperature_fac_ops);
      // Compute the product of DTDPhi block of the Jacobian with the
      // just computed phi correction
      // double norm_phase=cellops.L2Norm(d_phase_sol_id);
      // tbox::pout << "Off-diagonal Preconditioner for temperature block,
      // norm phase="<<norm_phase << std::endl;
      d_phase_temperature_fac_ops->multiplyDTDPhiBlock(d_phase_sol_id,
                                                       d_temperature_rhs_id);

      // Add gamma times the just computed product to the right-hand side
      cellops.axpy(d_temperature_rhs_id, gamma, d_temperature_rhs_id,
                   r_temperature_id, false);
   } else {
      // Copy the right-hand side to the temporary right-hand side array
      cellops.copyData(d_temperature_rhs_id, r_temperature_id, false);
   }

   // Zero out the initial guess in the temporary solution array
   cellops.setToScalar(d_temperature_sol_id, 0., false);

   // Set the tolerance for the FAC solve as requested by integrator
   d_temperature_sys_solver->setResidualTolerance(delta);

   // Solve the temperature block system
   bool converged = d_temperature_sys_solver->solveSystem(d_temperature_sol_id,
                                                          d_temperature_rhs_id);
#if 0
   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++ ) {

      std::shared_ptr<hier::PatchLevel > level = hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
         std::shared_ptr<hier::Patch > patch = *p;

         std::shared_ptr< pdat::CellData<double> > y (
            SAMRAI_SHARED_PTR_CAST< pdat::CellData<double>, hier::PatchData>(patch->getPatchData( d_temperature_sol_id ) ) );
         const hier::Box& pbox = patch->getBox();

         y->print(y->getGhostBox() );
      }
   }
#endif
   int retcode = converged ? 0 : 1;

   // Copy solution from the local temporary temperature_sol to the output
   // array, z_temperature, including ghost values
   cellops.copyData(z_temperature_id, d_temperature_sol_id, false);

   return retcode;
}

//-----------------------------------------------------------------------
int QuatIntegrator::ConcentrationPrecondSolve(
    std::shared_ptr<hier::PatchHierarchy> hierarchy,
    std::shared_ptr<solv::SAMRAIVectorReal<double> > r_samvect,
    std::shared_ptr<solv::SAMRAIVectorReal<double> > z_samvect,
    const double delta)
{
   t_conc_precond_timer->start();

   if (d_show_conc_sys_stats) {
      tbox::plog << "Preconditioner for Concentration block with tol " << delta
                 << std::endl;
   }

   math::HierarchyCellDataOpsReal<double> cellops(hierarchy);

   const int z_conc_id =
       z_samvect->getComponentDescriptorIndex(d_conc_component_index);
   const int r_conc_id =
       r_samvect->getComponentDescriptorIndex(d_conc_component_index);

   assert(z_conc_id >= 0);
   assert(r_conc_id >= 0);

   // Copy the right-hand side to the temporary right-hand side array
   cellops.copyData(d_conc_rhs_id, r_conc_id, false);

   // Set the tolerance for the FAC solve as requested by integrator
   d_conc_sys_solver->setResidualTolerance(delta);

   // Zero out the initial guess in the temporary solution array
   cellops.setToScalar(d_conc_sol_id, 0., false);

   // Solve the concentration block system
   bool converged =
       d_conc_sys_solver->solveSystem(d_conc_sol_id, d_conc_rhs_id);

   int retcode = converged ? 0 : 1;

   // Copy solution from the local temporary to the output array
   cellops.copyData(z_conc_id, d_conc_sol_id, false);

   t_conc_precond_timer->stop();

   return retcode;
}

//-----------------------------------------------------------------------
int QuatIntegrator::QuatPrecondSolve(
    std::shared_ptr<hier::PatchHierarchy> hierarchy, int r_quat_id,
    int z_quat_id, const double delta, const double gamma)
{
   if (d_show_quat_sys_stats) {
      tbox::plog << "Preconditioner for Quaternion block with tol " << delta
                 << std::endl;
   }

   math::HierarchyCellDataOpsReal<double> cellops(hierarchy);
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(cellops.max(r_quat_id) == cellops.max(r_quat_id));
#endif
   if (d_precond_has_dquatdphi) {

      assert(!d_use_gradq_for_flux);

      // Compute the product of DQuatDPhi block of the Jacobian with the
      // just computed phi correction
      d_quat_sys_solver->multiplyDQuatDPhiBlock(d_phase_sol_id, d_quat_rhs_id);

      // Add gamma times the just computed product to the right-hand side
      cellops.axpy(d_quat_rhs_id, gamma, d_quat_rhs_id, r_quat_id, false);
   } else {

#ifdef DEBUG_CHECK_ASSERTIONS
      assert(cellops.max(r_quat_id) == cellops.max(r_quat_id));
#endif
      // Copy the right-hand side to the temporary right-hand side array
      cellops.copyData(d_quat_rhs_id, r_quat_id, false);

#ifdef DEBUG_CHECK_ASSERTIONS
      assert(cellops.max(d_quat_rhs_id) == cellops.max(d_quat_rhs_id));
#endif
   }

   // Set the tolerance for the FAC solve as requested by integrator
   d_quat_sys_solver->setResidualTolerance(delta);

   // Zero out the initial guess in the temporary solution array
   cellops.setToScalar(d_quat_sol_id, 0., false);

   // Solve the quaternion block system
   bool converged =
       d_quat_sys_solver->solveSystem(d_quat_sol_id, d_quat_rhs_id);

   int retcode = converged ? 0 : 1;

   // Copy solution from the local temporary to the output array, z_quat_id
   cellops.copyData(z_quat_id, d_quat_sol_id, false);

   return retcode;
}

//-----------------------------------------------------------------------
// Virtual function from CVODEAbstractFunction

/*!
  Solve the system Mz = r for z given r, where M is an
  approximation of I - gamma J, where J is the Jacobian
  of the ODE right-hand side f (whose evaluation is
  performed by the evaluateRHS() member of this class).
  The input parameter gamma is the reciprocal of the
  BDF coefficient of the current time solution computed
  by CPODE and is proportional to the current time step.
  The parameter delta is the solver tolerance requested
  by CPODE.  This routine assumes that M has already been
  setup by a prior call to CPSpgmrPrecondSet().

  Returns 0 if all block solves converged.
*/

int QuatIntegrator::CVSpgmrPrecondSolve(double t,
                                        solv::SundialsAbstractVector* y,
                                        solv::SundialsAbstractVector* fy,
                                        solv::SundialsAbstractVector* r,
                                        solv::SundialsAbstractVector* z,
                                        double gamma, double delta, int lr)
{
   (void)y;
   (void)fy;

   assert(d_use_preconditioner);
   // tbox::pout << "QuatIntegrator::CVSpgmrPrecondSolve" << std::endl;

   t_psolve_solve_timer->start();

   // Convert passed-in std::vectors into SAMRAI std::vectors
   std::shared_ptr<solv::SAMRAIVectorReal<double> > r_samvect =
       solv::Sundials_SAMRAIVector::getSAMRAIVector(r);
   std::shared_ptr<solv::SAMRAIVectorReal<double> > z_samvect =
       solv::Sundials_SAMRAIVector::getSAMRAIVector(z);

   int retcode = 0;

   if (lr == 1 || lr == 2) {  // Applying left or right preconditioner

      std::shared_ptr<hier::PatchHierarchy> hierarchy =
          r_samvect->getPatchHierarchy();
      math::HierarchyCellDataOpsReal<double> cellops(hierarchy);

      if (d_with_unsteady_temperature && d_precond_has_dPhidT) {
         int converged =
             applyTemperaturePreconditioner(hierarchy, t, r_samvect, z_samvect,
                                            delta, gamma);
         retcode = (converged == 0 && retcode == 0) ? 0 : 1;
      }
      if (d_with_phase) {
         // Apply the preconditioner phase block
         int converged = applyPhasePreconditioner(hierarchy, t, r_samvect,
                                                  z_samvect, delta, gamma);
         retcode = (converged == 0 && retcode == 0) ? 0 : 1;
      }

      if (d_with_third_phase) {

         int z_eta_id =
             z_samvect->getComponentDescriptorIndex(d_eta_component_index);
         int r_eta_id =
             r_samvect->getComponentDescriptorIndex(d_eta_component_index);

         assert(z_eta_id >= 0);
         assert(r_eta_id >= 0);

         if (!d_eta_sys_solver) {
            cellops.copyData(z_eta_id, r_eta_id, false);
         } else {
            int converged =
                EtaPrecondSolve(hierarchy, r_eta_id, z_eta_id, delta);
            retcode = (converged == 0 && retcode == 0) ? 0 : 1;
         }
      }

      // Apply the preconditioner quaternion block
      if (d_evolve_quat) {

         int z_quat_id =
             z_samvect->getComponentDescriptorIndex(d_quat_component_index);
         int r_quat_id =
             r_samvect->getComponentDescriptorIndex(d_quat_component_index);

         assert(z_quat_id >= 0);
         assert(r_quat_id >= 0);

         if (d_precondition_quat) {
            int converged =
                QuatPrecondSolve(hierarchy, r_quat_id, z_quat_id, delta, gamma);
            retcode = (converged == 0 && retcode == 0) ? 0 : 1;
         } else {  // !d_precondition_quat
            cellops.copyData(z_quat_id, r_quat_id, false);
         }
      }

      // Apply the preconditioner temperature block
      if (d_with_unsteady_temperature && !d_precond_has_dPhidT) {

         int converged =
             applyTemperaturePreconditioner(hierarchy, t, r_samvect, z_samvect,
                                            delta, gamma);
         retcode = (converged == 0 && retcode == 0) ? 0 : 1;
      }

      // Apply the preconditioner concentration block
      if (d_with_concentration) {

         int converged = applyConcentrationPreconditioner(hierarchy, r_samvect,
                                                          z_samvect, delta);
         retcode = (converged == 0 && retcode == 0) ? 0 : 1;
      }
   } else {  // Identity (no) preconditioner
      z_samvect->copyVector(r_samvect);
   }

   t_psolve_solve_timer->stop();

   return retcode;
}

//-----------------------------------------------------------------------

int QuatIntegrator::applyPhasePreconditioner(
    std::shared_ptr<hier::PatchHierarchy> hierarchy, const double t,
    std::shared_ptr<solv::SAMRAIVectorReal<double> > r_samvect,
    std::shared_ptr<solv::SAMRAIVectorReal<double> > z_samvect,
    const double delta, const double gamma)
{
   int retcode = 0;

   const int z_phase_id =
       z_samvect->getComponentDescriptorIndex(d_phase_component_index);
   const int r_phase_id =
       r_samvect->getComponentDescriptorIndex(d_phase_component_index);

   assert(z_phase_id >= 0);
   assert(r_phase_id >= 0);

   if (!d_phase_sys_solver) {
      math::HierarchyCellDataOpsReal<double> cellops(hierarchy);
      cellops.copyData(z_phase_id, r_phase_id, false);

      if (d_precond_has_dquatdphi || d_precond_has_dTdphi) {

         // Copy the phase correction to an array with ghost cells and fill
         // them
         xfer::RefineAlgorithm copy_with_ghosts;

         copy_with_ghosts.registerRefine(
             d_phase_sol_id,  // destination
             z_phase_id,      // source
             d_phase_sol_id,  // temporary work space
             d_phase_refine_op);

         for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++) {
            std::shared_ptr<hier::PatchLevel> level =
                hierarchy->getPatchLevel(ln);

            std::shared_ptr<xfer::RefineSchedule> schedule(
                copy_with_ghosts.createSchedule(level, ln - 1, hierarchy,
                                                d_all_refine_patch_strategy));
            schedule->fillData(t);
         }
      }
   } else {

      int converged =
          PhasePrecondSolve(hierarchy, r_phase_id, z_phase_id, delta, gamma);
      retcode = (converged == 0 && retcode == 0) ? 0 : 1;
   }

   return retcode;
}

//-----------------------------------------------------------------------

int QuatIntegrator::applyTemperaturePreconditioner(
    std::shared_ptr<hier::PatchHierarchy> hierarchy, const double t,
    std::shared_ptr<solv::SAMRAIVectorReal<double> > r_samvect,
    std::shared_ptr<solv::SAMRAIVectorReal<double> > z_samvect,
    const double delta, const double gamma)
{
   int retcode = 0;

   const int z_temperature_id =
       z_samvect->getComponentDescriptorIndex(d_temperature_component_index);
   const int r_temperature_id =
       r_samvect->getComponentDescriptorIndex(d_temperature_component_index);

   assert(z_temperature_id >= 0);
   assert(r_temperature_id >= 0);

   if (!d_temperature_sys_solver) {
      math::HierarchyCellDataOpsReal<double> cellops(hierarchy);
      cellops.copyData(z_temperature_id, r_temperature_id, false);

      // save computed T correction in d_temperature_sol_id to use in
      // Phase off-diagonal block preconditioner
      //(done by TemperaturePrecondSolve if called)
      if (d_precond_has_dPhidT) {
         // Copy the T correction to an array with ghost cells and fill them
         xfer::RefineAlgorithm copy_with_ghosts;

         copy_with_ghosts.registerRefine(
             d_temperature_sol_id,  // destination
             z_temperature_id,      // source
             d_temperature_sol_id,  // temporary work space
             d_temperature_refine_op);

         for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++) {
            std::shared_ptr<hier::PatchLevel> level =
                hierarchy->getPatchLevel(ln);

            std::shared_ptr<xfer::RefineSchedule> schedule(
                copy_with_ghosts.createSchedule(level, ln - 1, hierarchy,
                                                d_all_refine_patch_strategy));
            schedule->fillData(t);
         }
      }
   } else {
      int converged = TemperaturePrecondSolve(hierarchy, r_temperature_id,
                                              z_temperature_id, delta, gamma);
      retcode = (converged == 0 && retcode == 0) ? 0 : 1;
   }

   return retcode;
}

//-----------------------------------------------------------------------

int QuatIntegrator::applyConcentrationPreconditioner(
    std::shared_ptr<hier::PatchHierarchy> hierarchy,
    std::shared_ptr<solv::SAMRAIVectorReal<double> > r_samvect,
    std::shared_ptr<solv::SAMRAIVectorReal<double> > z_samvect,
    const double delta)
{
   int retcode = 0;

   if (!d_conc_sys_solver) {
      const int z_conc_id =
          z_samvect->getComponentDescriptorIndex(d_conc_component_index);
      const int r_conc_id =
          r_samvect->getComponentDescriptorIndex(d_conc_component_index);

      assert(z_conc_id >= 0);
      assert(r_conc_id >= 0);

      math::HierarchyCellDataOpsReal<double> cellops(hierarchy);
      cellops.copyData(z_conc_id, r_conc_id, false);
   } else {
      retcode =
          ConcentrationPrecondSolve(hierarchy, r_samvect, z_samvect, delta);
   }

   return retcode;
}

//-----------------------------------------------------------------------
// Virtual function from CVODEAbstractFunction

int QuatIntegrator::applyProjection(double time,
                                    solv::SundialsAbstractVector* y,
                                    solv::SundialsAbstractVector* corr,
                                    double epsProj,
                                    solv::SundialsAbstractVector* err)
{
   (void)time;

   // Zero all components of the correction
   corr->setToScalar(0.);

   if (d_qlen > 1 && d_evolve_quat) {
      // tbox::pout<<"QuatIntegrator::applyProjection()"<<endl;
      // Convert the Sundials std::vectors to SAMRAI std::vectors
      std::shared_ptr<solv::SAMRAIVectorReal<double> > y_samvect =
          solv::Sundials_SAMRAIVector::getSAMRAIVector(y);
      std::shared_ptr<solv::SAMRAIVectorReal<double> > corr_samvect =
          solv::Sundials_SAMRAIVector::getSAMRAIVector(corr);
      std::shared_ptr<solv::SAMRAIVectorReal<double> > err_samvect =
          solv::Sundials_SAMRAIVector::getSAMRAIVector(err);

      const int q_id =
          y_samvect->getComponentDescriptorIndex(d_quat_component_index);
      const int corr_id =
          corr_samvect->getComponentDescriptorIndex(d_quat_component_index);
      const int err_id =
          err_samvect->getComponentDescriptorIndex(d_quat_component_index);
      // tbox::pout<<"q_id="<<q_id<<endl;
      d_quat_sys_solver->applyProjection(q_id, corr_id, err_id);
   }

   if (d_model_parameters.use_FolchPlapp()) {
      // tbox::pout << "3 phases projection..." << std::endl;
      std::shared_ptr<solv::SAMRAIVectorReal<double> > y_samvect =
          solv::Sundials_SAMRAIVector::getSAMRAIVector(y);
      std::shared_ptr<solv::SAMRAIVectorReal<double> > corr_samvect =
          solv::Sundials_SAMRAIVector::getSAMRAIVector(corr);
      std::shared_ptr<solv::SAMRAIVectorReal<double> > err_samvect =
          solv::Sundials_SAMRAIVector::getSAMRAIVector(err);

      const int phi_id =
          y_samvect->getComponentDescriptorIndex(d_phase_component_index);
      const int corr_id =
          corr_samvect->getComponentDescriptorIndex(d_phase_component_index);
      const int err_id =
          err_samvect->getComponentDescriptorIndex(d_phase_component_index);
      std::shared_ptr<FolchPlappRHSStrategy> phase_rhs_strategy =
          std::dynamic_pointer_cast<FolchPlappRHSStrategy>(
              d_phase_rhs_strategy);
      assert(phase_rhs_strategy);
      phase_rhs_strategy->projectPhases(phi_id, corr_id, err_id);
   }

   return 0;  // Always successful
}

//=======================================================================

void QuatIntegrator::correctRhsForSymmetry(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const int quat_id,
    const int quat_rhs_id)
{
   t_symm_rhs_timer->start();

   assert(quat_rhs_id != -1);
   assert(d_quat_diffs_id != -1);

   // tbox::pout<<"QuatIntegrator::correctRhsForSymmetry()"<<endl;
   const int face_coeff_id = d_quat_sys_solver->getFaceDiffCoeffScratchId();
   assert(face_coeff_id >= 0);

   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
             SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                                    hier::PatchGeometry>(
                 patch->getPatchGeometry()));
         const double* dx = patch_geom->getDx();

         const hier::Box& box = patch->getBox();
         const hier::Index& lower = box.lower();
         const hier::Index& upper = box.upper();

         std::shared_ptr<pdat::CellData<double> > q_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(quat_id)));
         assert(q_data);

         std::shared_ptr<pdat::SideData<double> > face_coeff(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(face_coeff_id)));
         assert(face_coeff);
         assert(face_coeff->getDepth() == d_qlen);
         assert(face_coeff->getGhostCellWidth() ==
                hier::IntVector(tbox::Dimension(NDIM), 0));

         std::shared_ptr<pdat::SideData<double> > quat_diffs(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(d_quat_diffs_id)));
         assert(quat_diffs);
         assert(quat_diffs->getDepth() == 2 * d_qlen);

         std::shared_ptr<pdat::CellData<double> > mobility(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_quat_mobility_id)));
         assert(mobility);
         assert(mobility->getDepth() == 1);
         assert(mobility->getGhostCellWidth()[0] > 0);

         std::shared_ptr<pdat::CellData<double> > rhs(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(quat_rhs_id)));
         assert(rhs);
         assert(rhs->getDepth() == d_qlen);

         std::shared_ptr<pdat::SideData<int> > rotation_index(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<int>, hier::PatchData>(
                 patch->getPatchData(d_quat_symm_rotation_id)));
         assert(rotation_index);

         const int symm_offset = 0 * d_qlen;
         const int nonsymm_offset = 1 * d_qlen;

         CORRECTRHSQUATFORSYMMETRY(
             lower[0], upper[0], lower[1], upper[1],
#if (NDIM == 3)
             lower[2], upper[2],
#endif
             d_qlen, dx, quat_diffs->getPointer(0, nonsymm_offset),
             quat_diffs->getPointer(1, nonsymm_offset),
#if (NDIM == 3)
             quat_diffs->getPointer(2, nonsymm_offset),
#endif
             quat_diffs->getPointer(0, symm_offset),
             quat_diffs->getPointer(1, symm_offset),
#if (NDIM == 3)
             quat_diffs->getPointer(2, symm_offset),
#endif
             quat_diffs->getGhostCellWidth()[0], rhs->getPointer(),
             rhs->getGhostCellWidth()[0], q_data->getPointer(),
             q_data->getGhostCellWidth()[0], face_coeff->getPointer(0),
             face_coeff->getPointer(1),
#if (NDIM == 3)
             face_coeff->getPointer(2),
#endif
             face_coeff->getGhostCellWidth()[0], mobility->getPointer(),
             mobility->getGhostCellWidth()[0], rotation_index->getPointer(0),
             rotation_index->getPointer(1),
#if (NDIM == 3)
             rotation_index->getPointer(2),
#endif
             rotation_index->getGhostCellWidth()[0]);
      }
   }

   t_symm_rhs_timer->stop();
}

//=======================================================================
#ifdef USE_CPODE
std::vector<std::shared_ptr<solv::SAMRAIVectorReal<double> > >* QuatIntegrator::
    getCPODESVectorsRequiringRegrid(void)
{
   assert(d_sundials_solver != nullptr);

   std::vector<std::shared_ptr<solv::SAMRAIVectorReal<double> > >* cpodes_vec =
       new std::vector<std::shared_ptr<solv::SAMRAIVectorReal<double> > >;

   std::vector<solv::SundialsAbstractVector*>* sundials_vec =
       d_sundials_solver->getVectorsRequiringRegrid();

   std::vector<solv::SundialsAbstractVector*>::iterator it;

   for (it = sundials_vec->begin(); it < sundials_vec->end(); it++) {
      std::shared_ptr<solv::SAMRAIVectorReal<double> > samvec =
          solv::Sundials_SAMRAIVector::getSAMRAIVector(*it);

      cpodes_vec->push_back(samvec);
   }

   delete sundials_vec;

   return cpodes_vec;
}
#endif

//-----------------------------------------------------------------------

#ifdef USE_CPODE
void QuatIntegrator::getCPODESIdsRequiringRegrid(std::set<int>& cpode_id_set,
                                                 std::set<int>& phase_id_set,
                                                 std::set<int>& eta_id_set,
                                                 std::set<int>& orient_id_set,
                                                 std::set<int>& conc_id_set,
                                                 std::set<int>& temp_id_set)
{
   std::vector<std::shared_ptr<solv::SAMRAIVectorReal<double> > >* cpodes_vec =
       getCPODESVectorsRequiringRegrid();

   for (std::vector<std::shared_ptr<solv::SAMRAIVectorReal<double> > >::iterator
            it = cpodes_vec->begin();
        it < cpodes_vec->end(); it++) {

      if (d_with_phase) {
         assert(d_phase_component_index != -1);
         int id = (*it)->getComponentDescriptorIndex(d_phase_component_index);

         if (id != d_phase_id && phase_id_set.find(id) == phase_id_set.end()) {
            phase_id_set.insert(id);
            cpode_id_set.insert(id);
         }
      }

      if (d_with_third_phase) {
         assert(d_eta_component_index != -1);
         int id = (*it)->getComponentDescriptorIndex(d_eta_component_index);

         if (id != d_eta_id && eta_id_set.find(id) == eta_id_set.end()) {
            eta_id_set.insert(id);
            cpode_id_set.insert(id);
         }
      }

      if (d_evolve_quat) {
         int id = (*it)->getComponentDescriptorIndex(d_quat_component_index);

         if (id != d_quat_id && orient_id_set.find(id) == orient_id_set.end()) {
            orient_id_set.insert(id);
            cpode_id_set.insert(id);
         }
      }

      if (d_with_concentration) {
         assert(d_conc_component_index != -1);
         int id = (*it)->getComponentDescriptorIndex(d_conc_component_index);

         if (id != d_conc_id && conc_id_set.find(id) == conc_id_set.end()) {
            conc_id_set.insert(id);
            cpode_id_set.insert(id);
         }
      }

      if (d_with_unsteady_temperature) {
         assert(d_temperature_component_index != -1);
         int id =
             (*it)->getComponentDescriptorIndex(d_temperature_component_index);

         if (id != d_temperature_id &&
             temp_id_set.find(id) == temp_id_set.end()) {
            temp_id_set.insert(id);
            cpode_id_set.insert(id);
         }
      }
   }

   delete cpodes_vec;
}
#endif
