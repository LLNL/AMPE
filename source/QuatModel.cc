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
#include "QuatModel.h"
#include "QuatIntegrator.h"
#include "SimpleGradStrategy.h"
#include "SimpleQuatGradStrategy.h"
#include "QuatGradModulusStrategy.h"
#include "TemperatureFreeEnergyStrategy.h"
#include "CALPHADFreeEnergyStrategyBinary.h"
#include "CALPHADFreeEnergyStrategyWithPenalty.h"
#include "ConstantTemperatureStrategy.h"
#include "SteadyStateTemperatureStrategy.h"
#include "QuatFort.h"
#include "QuatLinearRefine.h"
#include "QuatWeightedAverage.h"
#include "MinIntCoarsen.h"
#include "ConcFort.h"
#include "ConstantMolarVolumeStrategy.h"
#include "TemperatureStrategyFactory.h"
#include "ConstantHeatCapacityStrategy.h"
#include "NKRHeatCapacityStrategy.h"
#include "PhaseFluxStrategyFactory.h"
#include "AzizPartitionCoefficientStrategy.h"
#include "UniformPartitionCoefficientStrategy.h"
#include "ConstantMeltingTemperatureStrategy.h"
#include "LinearMeltingTemperatureStrategy.h"
#include "QuatIntegratorFactory.h"
#include "CompositionStrategyMobilities.h"
#include "TbasedCompositionDiffusionStrategy.h"
#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "toolsSAMRAI.h"
#include "tools.h"
#include "MobilityFactory.h"
#include "CompositionDiffusionStrategyFactory.h"
#include "CompositionRHSStrategyFactory.h"
#include "PhaseConcentrationsStrategyFactory.h"
#include "FreeEnergyStrategyFactory.h"
#include "FuncFort.h"
#include "diagnostics.h"
#include "FieldsWriter.h"
#include "TwoPhasesEnergyEvaluationStrategy.h"
#include "computeQDiffs.h"
#include "HierarchyStencilOps.h"
#include "WangRigidBodyForces.h"

#include "Database2JSON.h"
namespace pt = boost::property_tree;

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/math/PatchCellDataBasicOps.h"
#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/hier/PatchDataRestartManager.h"

#include "PhysicalConstants.h"

#include <set>
#include <map>

const double gas_constant_R_JpKpmol = GASCONSTANT_R_JPKPMOL;

const double um2tom2 = 1.e-12;

QuatModel::QuatModel(int ql) : d_qlen(ql), d_ncompositions(-1)
{
   d_symmetry_aware = false;

   d_integrator.reset();
   d_integrator_quat_only.reset();
   d_quat_grad_strategy = nullptr;
   d_mobility_strategy = nullptr;
   d_all_refine_patch_strategy = nullptr;
   d_partition_coeff_refine_patch_strategy = nullptr;
   d_phase_conc_strategy = nullptr;
   d_partition_coeff_strategy = nullptr;

   d_meltingT_strategy = nullptr;

   d_composition_strategy_mobilities = nullptr;

   d_heat_capacity_strategy = nullptr;


   d_test_interval.reset();
   d_fundamental_interval.reset();
   d_scalar_diag_interval.reset();
   d_grain_extend_interval.reset();

   d_phase_id = -1;
   d_eta_id = -1;
   d_quat_id = -1;
   d_quat_relax_id = -1;
   d_conc_id = -1;
   d_conc_l_ref_id = -1;
   d_conc_a_ref_id = -1;
   d_conc_b_ref_id = -1;
   d_phase_scratch_id = -1;
   d_eta_scratch_id = -1;
   d_quat_scratch_id = -1;
   d_conc_scratch_id = -1;
   d_conc_l_id = -1;
   d_conc_a_id = -1;
   d_conc_b_id = -1;
   d_phase_grad_cell_id = -1;
   d_phase_grad_side_id = -1;
   d_phase_diffs_id = -1;
   d_phase_diffs_cell_id = -1;
   d_eta_grad_cell_id = -1;
   d_eta_grad_side_id = -1;
   d_eta_diffs_id = -1;
   d_quat_diffs_id = -1;
   d_quat_diffs_cell_id = -1;
   d_quat_nonsymm_diffs_cell_id = -1;
   d_quat_grad_cell_id = -1;
   d_quat_grad_side_id = -1;
   d_quat_grad_modulus_id = -1;
   d_quat_symm_rotation_id = -1;
   d_quat_symm_rotation_cell_id = -1;
   d_phase_mobility_id = -1;
   d_eta_mobility_id = -1;
   d_quat_mobility_id = -1;
   d_quat_norm_error_id = -1;
   d_weight_id = -1;
   d_weight_diagnostics_id = -1;
   d_work_id = -1;
   d_conc_diffusion_id = -1;
   d_conc_phase_coupling_diffusion_id = -1;
   d_conc_eta_coupling_diffusion_id = -1;
   d_conc_pfm_diffusion_l_id = -1;
   d_conc_pfm_diffusion_a_id = -1;
   d_conc_pfm_diffusion_b_id = -1;
   d_conc_diffusion_coeff_l_id = -1;
   d_conc_diffusion_coeff_a_id = -1;
   d_conc_diffusion_coeff_b_id = -1;
   d_velocity_id = -1;
   d_partition_coeff_id = -1;
   d_partition_coeff_scratch_id = -1;
   d_temperature_id = -1;
   d_temperature_scratch_id = -1;
   d_temperature_rhs_steady_id = -1;
   d_f_l_id = -1;
   d_f_a_id = -1;
   d_f_b_id = -1;
   d_cp_id = -1;
   d_conc_Mq_id = -1;
   d_visit_conc_pfm_diffusion_id = -1;

   d_number_of_grains = -1;
   d_phase_threshold = 0.85;

   d_nghosts = -1;
   d_nghosts_aux_conc = -1;

   d_use_warm_start = false;

   d_conc_l_var.reset();
   d_conc_l_ref_var.reset();

   double def_val = tbox::IEEE::getSignalingNaN();

   d_time = def_val;

   tbox::RestartManager::getManager()->registerRestartItem("QuatModel", this);

   d_fenergy_diag_filename = "";
   d_grains.reset();

   d_verbosity = new QuatVerbosity();
   PFModel::setVerbosity(d_verbosity);

   tbox::TimerManager* tman = tbox::TimerManager::getManager();
   t_resetGrains_timer = tman->getTimer("AMPE::QuatModel::resetGrains()");

   t_phase_diffs_timer = tman->getTimer("AMPE::QuatModel::phaseDiffs");
   t_phase_conc_timer =
       tman->getTimer("AMPE::QuatIntegrator::computePhaseConcentrations()");
}

//=======================================================================

QuatModel::~QuatModel()
{
   delete d_quat_grad_strategy;
   if (d_heat_capacity_strategy) delete d_heat_capacity_strategy;
   if (d_partition_coeff_strategy) delete d_partition_coeff_strategy;
   d_integrator.reset();
   d_integrator_quat_only.reset();

   delete d_verbosity;
}


//=======================================================================

void QuatModel::initializeTemperature(
    std::shared_ptr<tbox::Database> model_db,
    std::shared_ptr<tbox::Database> integrator_db)
{
   if (d_model_parameters.with_heat_equation()) {
      if (d_model_parameters.with_concentration()) {
         tbox::pout << "QuatModel::initializeTemperature() "
                    << "Using NKR model for heat capacity..." << std::endl;
         d_heat_capacity_strategy =
             new NKRHeatCapacityStrategy(d_model_parameters.cp(), d_cp_id,
                                         d_conc_id, d_temperature_id);
      } else {
         const std::map<short, double>& cp(d_model_parameters.cp(0));
         std::map<short, double>::const_iterator p = cp.find(0);
         const double cpval = p->second;
         d_heat_capacity_strategy =
             new ConstantHeatCapacityStrategy(cpval, d_cp_id);
      }
   }

   TemperatureStrategyFactory factory(d_temperature_id,
                                      d_temperature_scratch_id, d_conc_id,
                                      d_weight_id, d_temperature_rhs_steady_id,
                                      d_cp_id,
                                      d_model_parameters.molar_volume_liquid(),
                                      d_model_parameters.with_concentration(),
                                      d_grid_geometry,
                                      d_heat_capacity_strategy);

   d_temperature_strategy =
       factory.create(model_db, integrator_db, d_model_parameters);

   d_temperature_strategy_quat_only.reset(
       new ConstantTemperatureStrategy(d_temperature_id,
                                       d_temperature_scratch_id));
}

//=======================================================================

void QuatModel::initializeAmr(std::shared_ptr<tbox::Database> amr_db)
{
   d_use_warm_start = amr_db->getBoolWithDefault("use_warm_start", false);

   if (amr_db->isDatabase("TaggingCriteria")) {
      std::shared_ptr<tbox::Database> tag_db =
          amr_db->getDatabase("TaggingCriteria");

      if (tag_db->isDatabase("Phi")) {
         d_tag_phase = true;

         std::shared_ptr<tbox::Database> p_db = tag_db->getDatabase("Phi");

         d_phase_threshold_tagged = p_db->getDouble("threshold_tagged");
         d_phase_threshold_untagged = p_db->getDouble("threshold_untagged");
      }

      if (d_model_parameters.with_third_phase() && tag_db->isDatabase("Eta")) {
         d_tag_eta = true;
         if (!d_model_parameters.with_third_phase()) {
            d_tag_eta = false;
         }

         std::shared_ptr<tbox::Database> p_db = tag_db->getDatabase("Eta");

         d_eta_threshold_tagged = p_db->getDouble("threshold_tagged");
         d_eta_threshold_untagged = p_db->getDouble("threshold_untagged");
      }

      if (d_model_parameters.with_orientation() && tag_db->isDatabase("Orien"
                                                                      "t")) {
         d_tag_quat = true;

         std::shared_ptr<tbox::Database> q_db = tag_db->getDatabase("Orient");

         d_quat_threshold_tagged = q_db->getDouble("threshold_tagged");
         d_quat_threshold_untagged = q_db->getDouble("threshold_untagged");
      }
   }
}

//=======================================================================

void QuatModel::initializeRHSandEnergyStrategies(
    std::shared_ptr<tbox::MemoryDatabase>& input_db)
{
   tbox::plog << "QuatModel::initializeRHSandEnergyStrategies()" << std::endl;

   assert(d_ncompositions >= 0);

   const double Tref = d_model_parameters.with_rescaled_temperature()
                           ? d_model_parameters.meltingT() /
                                 d_model_parameters.rescale_factorT()
                           : d_model_parameters.meltingT();

   std::shared_ptr<tbox::Database> model_db =
       input_db->getDatabase("ModelParameters");

   if (d_model_parameters.with_phase())
      d_phase_flux_strategy =
          PhaseFluxStrategyFactory::create(d_model_parameters);

   std::shared_ptr<tbox::MemoryDatabase> calphad_db;

   if (d_model_parameters.with_concentration()) {
      d_conc_db = model_db->getDatabase("ConcentrationModel");

      if (d_model_parameters.isConcentrationModelCALPHAD() ||
          d_model_parameters.isConcentrationModelKKSdilute() ||
          d_model_parameters.isConcentrationModelParabolic()) {
         d_mvstrategy = new ConstantMolarVolumeStrategy(
             d_model_parameters.molar_volume_liquid(),
             d_model_parameters.molar_volume_solid_A(),
             d_model_parameters.molar_volume_solid_B());
      }

      if (d_model_parameters.isConcentrationModelCALPHAD()) {
         tbox::pout << "QuatModel: "
                    << "Using CALPHAD model for concentration" << std::endl;
         d_calphad_db = d_conc_db->getDatabase("Calphad");
      }
   }

   if (d_model_parameters.with_concentration()) {
      if (d_model_parameters.with_bias_well()) {
         if (d_model_parameters.isConcentrationModelLinear()) {
            d_meltingT_strategy = new LinearMeltingTemperatureStrategy(
                Tref, d_model_parameters.average_concentration(),
                d_model_parameters.liquidus_slope(), d_conc_l_id,
                d_equilibrium_temperature_id);

         } else {
            d_meltingT_strategy = new ConstantMeltingTemperatureStrategy(
                Tref, d_equilibrium_temperature_id);
         }
      }
   } else if (d_model_parameters.with_heat_equation()) {
      if (d_model_parameters.with_bias_well()) {
         d_meltingT_strategy = new ConstantMeltingTemperatureStrategy(
             Tref, d_equilibrium_temperature_id);
      }
   }

   d_free_energy_strategy =
       FreeEnergyStrategyFactory::create(d_model_parameters, d_ncompositions,
                                         d_conc_l_id, d_conc_a_id, d_conc_b_id,
                                         d_mvstrategy, d_meltingT_strategy,
                                         Tref, d_conc_db);

   if (d_model_parameters.with_concentration()) {
      if (d_model_parameters.concentrationModelNeedsPhaseConcentrations()) {

         // Newton solver parameters to solve KKS equations
         std::shared_ptr<tbox::Database> newton_db;
         if (d_conc_db->isDatabase("NewtonSolver")) {
            newton_db = d_conc_db->getDatabase("NewtonSolver");
         }

         d_phase_conc_strategy = PhaseConcentrationsStrategyFactory::create(
             d_model_parameters, d_conc_l_id, d_conc_a_id, d_conc_b_id,
             d_conc_l_ref_id, d_conc_a_ref_id, d_conc_b_ref_id,
             d_partition_coeff_id, d_ncompositions, d_conc_db, newton_db);
      }

      d_free_energy_strategy_for_diffusion = d_free_energy_strategy;

      // tbox::pout << "Initialize cl, ca, cb with c..." << std::endl;
      // math::HierarchyCellDataOpsReal<double> mathops(d_patch_hierarchy);
      // mathops.copyData(d_conc_l_id, d_conc_id);
      // mathops.copyData(d_conc_a_id, d_conc_id);
      // if (d_conc_b_id > -1) mathops.copyData(d_conc_b_id, d_conc_id);

      std::shared_ptr<tbox::Database> integrator_db =
          input_db->getDatabase("Integrator");
      if (d_model_parameters.concRHSstrategyIsEBS()) {
         if (d_model_parameters.conDiffusionStrategyIsCTD()) {
            d_composition_strategy_mobilities =
                new CompositionStrategyMobilities(d_calphad_db,
                                                  (d_eta_scratch_id > -1),
                                                  static_cast<unsigned short>(
                                                      d_ncompositions),
                                                  d_free_energy_strategy);
         }
      }
      if (d_model_parameters.concRHSstrategyIsEBS() ||
          d_model_parameters.concRHSstrategyIsWangSintering()) {
         d_diffusion_for_conc_in_phase =
             CompositionDiffusionStrategyFactory::create(
                 this, d_model_parameters, d_temperature_scratch_id,
                 static_cast<unsigned short>(d_ncompositions),
                 d_conc_scratch_id, d_conc_l_id, d_conc_a_id, d_conc_b_id,
                 d_conc_pfm_diffusion_id, d_conc_pfm_diffusion_l_id,
                 d_conc_pfm_diffusion_a_id, d_conc_pfm_diffusion_b_id,
                 d_conc_diffusion_coeff_l_id, d_conc_diffusion_coeff_a_id,
                 d_conc_diffusion_coeff_b_id, d_composition_strategy_mobilities,
                 d_free_energy_strategy);
      }

      d_composition_rhs_strategy = CompositionRHSStrategyFactory::create(
          this, d_free_energy_strategy_for_diffusion, d_model_parameters,
          d_ncompositions, d_conc_scratch_id, d_phase_scratch_id,
          d_temperature_scratch_id, d_conc_pfm_diffusion_id,
          d_conc_pfm_diffusion_l_id, d_conc_pfm_diffusion_a_id,
          d_conc_pfm_diffusion_b_id, d_conc_phase_coupling_diffusion_id,
          d_eta_scratch_id, d_conc_eta_coupling_diffusion_id, d_conc_l_id,
          d_conc_a_id, d_conc_b_id, d_partition_coeff_scratch_id, d_conc_Mq_id,
          d_composition_strategy_mobilities, d_diffusion_for_conc_in_phase);

   }  // d_model_parameters.with_concentration()

   if (d_model_parameters.with_Aziz_partition_coeff()) {
      setupAziz();
   }

   if (d_model_parameters.with_uniform_partition_coeff()) {

      // use d_temperature_scratch_id since that field will be uptodate when
      // integrator needs partion function
      d_partition_coeff_strategy =
          new UniformPartitionCoefficientStrategy(d_velocity_id,
                                                  d_temperature_scratch_id,
                                                  d_partition_coeff_id,
                                                  d_model_parameters.keq());
   }
}

void QuatModel::setupAziz()
{
   // use d_temperature_scratch_id since that field will be uptodate when
   // integrator needs partion function
   if (d_model_parameters.isConcentrationModelCALPHAD()) {
      if (d_ncompositions > 1)
         d_partition_coeff_strategy = new AzizPartitionCoefficientStrategy<
             Thermo4PFM::CALPHADFreeEnergyFunctionsTernary>(
             d_velocity_id, d_temperature_scratch_id, d_partition_coeff_id,
             d_model_parameters.vd(), d_model_parameters.keq(),
             d_model_parameters.energy_interp_func_type(),
             d_model_parameters.conc_interp_func_type(), d_conc_db);
      else {
         d_partition_coeff_strategy = new AzizPartitionCoefficientStrategy<
             Thermo4PFM::CALPHADFreeEnergyFunctionsBinary>(
             d_velocity_id, d_temperature_scratch_id, d_partition_coeff_id,
             d_model_parameters.vd(), d_model_parameters.keq(),
             d_model_parameters.energy_interp_func_type(),
             d_model_parameters.conc_interp_func_type(), d_conc_db);
      }
   } else
      d_partition_coeff_strategy = new AzizPartitionCoefficientStrategy<
          Thermo4PFM::KKSFreeEnergyFunctionDiluteBinary>(
          d_velocity_id, d_temperature_scratch_id, d_partition_coeff_id,
          d_model_parameters.vd(), d_model_parameters.keq(),
          d_model_parameters.energy_interp_func_type(),
          d_model_parameters.conc_interp_func_type(), d_conc_db);
}

//=======================================================================

void QuatModel::Initialize(std::shared_ptr<tbox::MemoryDatabase>& input_db,
                           const std::string& run_name,
                           const bool is_from_restart,
                           const std::string& restart_read_dirname,
                           const int restore_num)
{
   tbox::plog << "QuatModel::Initialize()" << std::endl;
   d_is_from_restart = is_from_restart;

   if (input_db->isDatabase("FreeEnergyDiagnostics")) {
      std::shared_ptr<tbox::Database> fenergy_diag_db =
          input_db->getDatabase("FreeEnergyDiagnostics");

      d_fenergy_diag_filename = fenergy_diag_db->getString("filename");
   }

   std::shared_ptr<tbox::Database> model_db =
       input_db->getDatabase("ModelParameters");

   d_model_parameters.readModelParameters(model_db);

   d_nghosts = d_model_parameters.nghosts_required();
   tbox::pout << "QuatModel using " << d_nghosts
              << " ghost cells for main field" << std::endl;
   assert(d_nghosts > 0);

   // we need internal composition with ghost values for EBS r.h.s.
   // in particular
   d_nghosts_aux_conc = d_nghosts;

   d_ncompositions = d_model_parameters.ncompositions();

   d_tag_phase = false;
   d_tag_eta = false;
   d_tag_quat = false;

   if (input_db->isDatabase("Amr")) {
      std::shared_ptr<tbox::Database> amr_db = input_db->getDatabase("Amr");
      d_amr_enabled = amr_db->getBoolWithDefault("enabled", true);

      initializeAmr(amr_db);
   } else {
      d_amr_enabled = false;
   }

   d_test_interval.reset(new EventInterval(input_db, "Test", 0.0, "step"));

   d_symmetry_aware = false;
   if (input_db->isDatabase("Symmetry")) {
      std::shared_ptr<tbox::Database> symm_db =
          input_db->getDatabase("Symmetry");

      if (symm_db->keyExists("enabled")) {
         d_symmetry_aware = symm_db->getBoolWithDefault("enabled", true);
      }

      d_fundamental_interval.reset(
          new EventInterval(symm_db, "Fundamental", 0.0, "step"));
   }

   d_scalar_diag_interval.reset(
       new EventInterval(input_db, "ScalarDiagnostics", 0.0, "step"));

   for (int i = 0; i < NDIM; i++)
      d_subdomain_diagnostics[i] = 0.;
   if (d_scalar_diag_interval->isActive()) {

      std::shared_ptr<tbox::Database> tmp_db =
          input_db->getDatabase("ScalarDiagnostics");

      d_extra_energy_detail =
          tmp_db->getBoolWithDefault("extra_energy_detail", false);

      // default value of 0 will lead to no diagnostics
      if (tmp_db->keyExists("domain_fraction"))
         tmp_db->getDoubleArray("domain_fraction", &d_subdomain_diagnostics[0],
                                NDIM);
   }

   d_grain_extend_interval.reset(
       new EventInterval(input_db, "GrainExtension", 0.0, "step"));

   d_ncompositions = d_model_parameters.ncompositionFields();

   // Base class does setup of various common things: logfile,
   // restart, basic samrai objects, visit, ALSO VIRTUAL FUNCTIONS for
   // creating and initializing the integrator and registering
   // variables.


   d_model_parameters.readFreeEnergies(model_db);

   EventInterval tmp_interval(input_db, "Visit", 0.0, "step");

   if (tmp_interval.isActive()) {

      std::shared_ptr<tbox::Database> visit_db = input_db->getDatabase("Visit");

      d_model_parameters.readVisitOptions(visit_db);
   }

   d_grains.reset(new Grains(d_qlen,
                             d_model_parameters.with_visit_grain_output(),
                             input_db));

   PFModel::Initialize(input_db, run_name, is_from_restart,
                       restart_read_dirname, restore_num);

   d_grains->initialize(input_db, d_all_periodic);

   // Set up Dirichlet boundary conditions
   if (!d_all_periodic) {
      std::shared_ptr<tbox::Database> bc_db =
          model_db->getDatabase("BoundaryConditions");

      const int phase_id =
          d_model_parameters.with_phase() ? d_phase_scratch_id : -1;
      double factor = d_model_parameters.with_rescaled_temperature()
                          ? 1. / d_model_parameters.rescale_factorT()
                          : -1.;
      d_all_refine_patch_strategy =
          new QuatRefinePatchStrategy("QuatRefinePatchStrategy", bc_db,
                                      phase_id, d_eta_scratch_id,
                                      d_quat_scratch_id, d_conc_scratch_id,
                                      d_temperature_scratch_id, d_nghosts,
                                      factor);

      if (d_model_parameters.needGhosts4PartitionCoeff())
         d_partition_coeff_refine_patch_strategy =
             new PartitionCoeffRefinePatchStrategy(
                 "PartitionCoeffRefinePatchStrategy", bc_db,
                 d_partition_coeff_scratch_id);
   }

   std::shared_ptr<tbox::Database> integrator_db =
       input_db->getDatabase("Integrator");

   initializeTemperature(model_db, integrator_db);

   if (!is_from_restart) {
      d_number_of_grains = INT_MAX;

      // Read initialization database (initial conditions)
      readInitialDatabase(input_db);

      setupHierarchy();

      // Create single level hierarchy for initial data
      setupInitialDataLevel();

      // Fill single level from initial data file
      FieldsInitializer initializer(d_grid_geometry,
                                    d_ratio_of_init_to_coarsest,
                                    d_verbosity->amrLevel());

      int temperature_id =
          d_model_parameters.isTemperatureConstant() ? d_temperature_id : -1;
      int qlen = (d_model_parameters.H_parameter() >= 0.) ? d_qlen : 0;
      int norderp_to_read =
          d_model_parameters.use_FolchPlapp()
              ? 3
              : d_model_parameters.norderpA() + d_model_parameters.norderpB();
      initializer.registerFieldsIds(d_phase_id, d_eta_id, temperature_id,
                                    d_quat_id, qlen, d_conc_id, d_ncompositions,
                                    norderp_to_read);

      if (!d_init_c.empty()) initializer.setCvalue(d_init_c);
      if (!d_init_q.empty()) initializer.setQvalue(d_init_q);
      if (d_init_t >= 0.) initializer.setTvalue(d_init_t);

      initializer.initializeLevelFromData(d_initial_level, d_init_data_filename,
                                          d_slice_index);
      // rescale initial conditions for temperature if we are solving
      // time evolution equation for T since we use reduced units
      math::HierarchyCellDataOpsReal<double> hopscell(d_patch_hierarchy);
      if (d_model_parameters.with_rescaled_temperature()) {
         assert(d_model_parameters.meltingT() == d_model_parameters.meltingT());

         hopscell.scale(d_temperature_id,
                        1. / d_model_parameters.rescale_factorT(),
                        d_temperature_id);
      }

      for (int ll = 0; ll < d_patch_hierarchy->getMaxNumberOfLevels(); ll++) {
         d_tag_buffer_array[ll] = 1;
      }

      if (d_model_parameters.concentrationModelNeedsPhaseConcentrations())
         setAuxilliaryCompositions();

   } else {  // restart case

      const int max_levels = d_patch_hierarchy->getMaxNumberOfLevels();

      d_patch_hierarchy->initializeHierarchy();

      for (int ll = 0; ll < max_levels; ll++) {
         d_tag_buffer_array[ll] = 1;
      }

      for (int ll = 0; ll < max_levels; ll++) {
         std::shared_ptr<hier::PatchLevel> null_level_ptr;

         initializeLevelData(d_patch_hierarchy, ll, -1.0, false, false,
                             null_level_ptr, true);
      }

      d_gridding_algorithm->getTagAndInitializeStrategy()
          ->resetHierarchyConfiguration(
              d_patch_hierarchy, 0, d_patch_hierarchy->getFinestLevelNumber());

      tbox::RestartManager::getManager()->closeRestartFile();

      d_integrator->setTimestep(d_previous_timestep);
   }

   // if temperature is a prescribed field, set it here, for restart case
   // or not
   if (d_model_parameters.isTemperatureUniform() ||
       d_model_parameters.isTemperatureGaussian() ||
       d_model_parameters.isTemperatureGradient()) {
      d_temperature_strategy->setCurrentTemperature(d_patch_hierarchy, 0.0);
   }

   d_quat_grad_strategy = new SimpleQuatGradStrategy(this);
   d_quat_grad_modulus_strategy.reset(new QuatGradModulusStrategy(d_qlen));

   initializeRHSandEnergyStrategies(input_db);

   if (input_db->isDatabase("GrainDiagnostics")) {
      std::shared_ptr<tbox::Database> g_diag_db =
          input_db->getDatabase("GrainDiagnostics");

      if (g_diag_db->keyExists("phase_threshold")) {
         d_phase_threshold = g_diag_db->getDouble("phase_threshold");
      }
   }

   if (d_subdomain_diagnostics[0] > 0.)
      computeVectorWeightsDiagnostics(d_patch_hierarchy);

   InitializeIntegrator();

   copyCurrentToScratch(d_patch_hierarchy, d_time, d_all_refine_patch_strategy);

   d_energy_eval_strategy.reset(new TwoPhasesEnergyEvaluationStrategy(
       d_model_parameters, d_qlen, d_phase_scratch_id, d_quat_scratch_id,
       d_quat_grad_side_id, d_weight_id, d_f_l_id, d_f_a_id, d_temperature_id,
       d_energy_diag_id));

   if (d_model_parameters.withMultipleOrderP()) {
      d_multiorderp_energy.reset(new MultiOrderPEnergyEvaluationStrategy(
          d_model_parameters, d_phase_scratch_id, d_weight_id,
          d_energy_diag_id));

      if (d_model_parameters.withRBmotion()) {

         if (d_model_parameters.rbThreshold() > 0.) {
            d_rigid_body_forces.reset(
                new WangRigidBodyForces(d_model_parameters, d_phase_scratch_id,
                                        d_conc_id, d_weight_id));
         } else {
            d_rigid_body_forces.reset(
                new PhaseRigidBodyForces(d_model_parameters, d_phase_scratch_id,
                                         d_weight_id));
         }
         d_integrator->setRigidBodyForces(d_rigid_body_forces);
      }
   }

   math::HierarchyCellDataOpsReal<double> cellops(d_patch_hierarchy);

   tbox::plog << "Set reference auxilliary compositions..." << std::endl;
   const double cLref = d_model_parameters.concL_ref();
   if (cLref >= 0.) {
      tbox::plog << "Set reference cL to " << cLref << std::endl;
      if (d_conc_l_ref_id >= 0)
         cellops.setToScalar(d_conc_l_ref_id, cLref, false);
      cellops.setToScalar(d_conc_l_id, cLref, false);
   }
   const double cAref = d_model_parameters.concA_ref();
   if (cAref >= 0.) {
      tbox::plog << "Set reference ca to " << cAref << std::endl;
      if (d_conc_a_ref_id >= 0)
         cellops.setToScalar(d_conc_a_ref_id, cAref, false);
      cellops.setToScalar(d_conc_a_id, cAref, false);
   }
   const double cBref = d_model_parameters.concB_ref();
   if (cBref >= 0.) {
      tbox::plog << "Set reference cb to " << cBref << std::endl;
      if (d_conc_b_ref_id >= 0)
         cellops.setToScalar(d_conc_b_ref_id, cBref, false);
      cellops.setToScalar(d_conc_b_id, cBref, false);
   }
}

//=======================================================================

void QuatModel::setAuxilliaryCompositions()
{
   assert(d_conc_l_id >= 0);
   assert(d_conc_a_id >= 0);

   int depth = 0;
   for (auto& init_cl : d_init_cl) {
      std::cout << "Set cl to " << init_cl << std::endl;
      for (int ln = d_patch_hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
         std::shared_ptr<hier::PatchLevel> level =
             d_patch_hierarchy->getPatchLevel(ln);
         for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
              ++ip) {
            std::shared_ptr<hier::Patch> patch = *ip;
            std::shared_ptr<pdat::CellData<double> > cl(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_conc_l_id)));
            assert(cl);
            cl->fill(init_cl, depth);
         }
      }
      depth++;
   }
   depth = 0;
   for (auto& init_ca : d_init_ca) {
      std::cout << "Set ca to " << init_ca << std::endl;
      for (int ln = d_patch_hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
         std::shared_ptr<hier::PatchLevel> level =
             d_patch_hierarchy->getPatchLevel(ln);
         for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
              ++ip) {
            std::shared_ptr<hier::Patch> patch = *ip;
            std::shared_ptr<pdat::CellData<double> > ca(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_conc_a_id)));
            assert(ca);
            ca->fill(init_ca, depth);
         }
      }
      depth++;
   }
   depth = 0;
   for (auto& init_cb : d_init_cb) {
      std::cout << "Set cb to " << init_cb << std::endl;
      for (int ln = d_patch_hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
         std::shared_ptr<hier::PatchLevel> level =
             d_patch_hierarchy->getPatchLevel(ln);
         for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
              ++ip) {
            std::shared_ptr<hier::Patch> patch = *ip;
            std::shared_ptr<pdat::CellData<double> > cb(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_conc_b_id)));
            assert(cb);
            cb->fill(init_cb, depth);
         }
      }
      depth++;
   }
}

//=======================================================================

void QuatModel::InitializeIntegrator(void)
{
   tbox::pout << "QuatModel::InitializeIntegrator()" << std::endl;

   if (d_model_parameters.with_heat_equation()) {
      assert(d_temperature_strategy);
      assert(d_heat_capacity_strategy);
   }

   if (d_model_parameters.with_phase()) {
      d_mobility_strategy =
          MobilityFactory::create(this, d_model_parameters, d_conc_l_id,
                                  d_conc_a_id, d_conc_b_id,
                                  d_temperature_scratch_id, d_ncompositions,
                                  d_model_parameters.use_FolchPlapp(),
                                  d_conc_db);
      d_integrator->setPhaseFluxStrategy(d_phase_flux_strategy);
   }

   d_integrator->setQuatGradStrategy(d_quat_grad_strategy);
   d_integrator->setMobilityStrategy(d_mobility_strategy);
   if (d_model_parameters.with_concentration()) {
      d_integrator->setCompositionDiffusionStrategy(
          d_diffusion_for_conc_in_phase);
      d_integrator->setCompositionRHSStrategy(d_composition_rhs_strategy);
   }
   d_integrator->setFreeEnergyStrategy(d_free_energy_strategy);
   if (d_model_parameters.with_partition_coeff())
      d_integrator->setPartitionCoefficientStrategy(d_partition_coeff_strategy);


   if (d_model_parameters.evolveQuat()) {
      d_integrator_quat_only->setQuatGradStrategy(d_quat_grad_strategy);
      d_integrator_quat_only->setMobilityStrategy(d_mobility_strategy);
      d_integrator_quat_only->setFreeEnergyStrategy(d_free_energy_strategy);
      d_integrator_quat_only->setTemperatureStrategy(
          d_temperature_strategy_quat_only);
   }

   d_integrator->setTemperatureStrategy(d_temperature_strategy);

   if (d_model_parameters.with_heat_equation()) {
      d_integrator->setHeatCapacityStrategy(d_heat_capacity_strategy);
   }

   if (d_model_parameters.with_concentration()) {
      d_integrator->setConcentrationModelParameters(
          d_model_parameters.conc_mobility());
   }

   d_integrator->initialize(d_patch_hierarchy);
   if (d_model_parameters.evolveQuat()) {
      const double atol = 1.e-4;
      const double rtol = atol * 1.e-2;
      tbox::pout << "QuatModel --- "
                 << "set tolerance for Quaternion only integrator:" << atol
                 << " and " << rtol << std::endl;
      d_integrator_quat_only->setAbsTol(atol);
      d_integrator_quat_only->setRelTol(rtol);
      d_integrator_quat_only->initialize(d_patch_hierarchy);
   }
}

void QuatModel::setIntegratorModelParameters()
{
   tbox::plog << "QuatModel::setIntegratorModelParameters()" << std::endl;
   d_integrator->setModelParameters(
       d_time, d_end_time, d_model_parameters.H_parameter(),
       d_model_parameters.epsilon_phase(), d_model_parameters.epsilon_eta(),
       d_model_parameters.epsilon_q(), d_model_parameters.quat_grad_floor(),
       d_model_parameters.quat_grad_floor_type(),
       d_model_parameters.phase_well_scale(),
       d_model_parameters.eta_well_scale(),
       d_model_parameters.orient_interp_func_type1(),
       d_model_parameters.orient_interp_func_type2(),
       d_model_parameters.diffq_avg_func_type(),
       d_model_parameters.phase_well_func_type(),
       d_model_parameters.energy_interp_func_type(),
       d_model_parameters.conc_interp_func_type(),
       d_model_parameters.eta_well_func_type());

   if (d_model_parameters.evolveQuat())
      d_integrator_quat_only->setModelParameters(
          d_time,
          1.e8,  // end_time large enough to never be reached
          d_model_parameters.H_parameter(), d_model_parameters.epsilon_phase(),
          d_model_parameters.epsilon_eta(), d_model_parameters.epsilon_q(),
          d_model_parameters.quat_grad_floor(), "s",
          d_model_parameters.phase_well_scale(),
          d_model_parameters.eta_well_scale(),
          d_model_parameters.orient_interp_func_type1(),
          d_model_parameters.orient_interp_func_type2(),
          "a",  // d_avg_func_type,
          d_model_parameters.phase_well_func_type(),
          d_model_parameters.energy_interp_func_type(),
          d_model_parameters.conc_interp_func_type(),
          d_model_parameters.eta_well_func_type());
}

//=======================================================================

void QuatModel::initializeRefineCoarsenAlgorithms()
{
   if (d_model_parameters.with_phase()) {
      d_phase_refine_op =
          d_grid_geometry->lookupRefineOperator(d_temperature_var,
                                                "LINEAR_REFINE");
   }
   if (d_model_parameters.with_third_phase()) {
      d_eta_refine_op =
          d_grid_geometry->lookupRefineOperator(d_eta_var, "LINEAR_REFINE");
   }

   if (d_model_parameters.with_orientation()) {
      if (d_symmetry_aware) {
         assert(d_quat_symm_rotation_id >= 0);
         if (d_verbosity->notSilent()) {
            tbox::pout << "QuatModel: "
                       << "Using symmetry aware refine/coarsen operators"
                       << std::endl;
         }
         d_quat_refine_op.reset(new QuatLinearRefine(d_quat_symm_rotation_id));
         d_quat_coarsen_op.reset(
             new QuatWeightedAverage(true, d_quat_symm_rotation_id));
      } else {
         d_quat_refine_op = d_grid_geometry->lookupRefineOperator(d_quat_var,
                                                                  "LINEAR_"
                                                                  "REFINE");
         d_quat_coarsen_op.reset(new QuatWeightedAverage(false));
      }
   }

   if (d_model_parameters.with_concentration()) {
      d_conc_refine_op =
          d_grid_geometry->lookupRefineOperator(d_conc_var, "LINEAR_REFINE");
   }

   d_curr_to_curr_refine_alg.reset(new xfer::RefineAlgorithm());

   d_curr_to_scr_refine_alg.reset(new xfer::RefineAlgorithm());


   // curr to curr
   if (d_model_parameters.with_phase()) {
      d_curr_to_curr_refine_alg->registerRefine(
          d_phase_id,          // destination
          d_phase_id,          // source
          d_phase_scratch_id,  // temporary
          d_phase_refine_op);
   }
   if (d_model_parameters.with_third_phase()) {
      d_curr_to_curr_refine_alg->registerRefine(d_eta_id,  // destination
                                                d_eta_id,  // source
                                                d_eta_scratch_id,  // temporary
                                                d_eta_refine_op);
   }

   assert(d_temperature_scratch_id >= 0);
   d_curr_to_curr_refine_alg->registerRefine(
       d_temperature_id,          // destination
       d_temperature_id,          // source
       d_temperature_scratch_id,  // temporary
       d_phase_refine_op);

   if (d_model_parameters.with_orientation()) {
      assert(d_quat_scratch_id >= 0);
      d_curr_to_curr_refine_alg->registerRefine(d_quat_id,  // destination
                                                d_quat_id,  // source
                                                d_quat_scratch_id,  // temporary
                                                d_quat_refine_op);
   }

   if (d_model_parameters.with_concentration()) {
      assert(d_conc_scratch_id >= 0);
      d_curr_to_curr_refine_alg->registerRefine(d_conc_id,  // destination
                                                d_conc_id,  // source
                                                d_conc_scratch_id,  // temporary
                                                d_conc_refine_op);
   }

   // curr to scr
   if (d_model_parameters.with_phase()) {
      d_curr_to_scr_refine_alg->registerRefine(
          d_phase_scratch_id,  // destination
          d_phase_id,          // source
          d_phase_scratch_id,  // temporary work space
          d_phase_refine_op);
   }
   if (d_model_parameters.with_third_phase()) {
      d_curr_to_scr_refine_alg->registerRefine(
          d_eta_scratch_id,  // destination
          d_eta_id,          // source
          d_eta_scratch_id,  // temporary work space
          d_eta_refine_op);
   }

   d_curr_to_scr_refine_alg->registerRefine(
       d_temperature_scratch_id,  // destination
       d_temperature_id,          // source
       d_temperature_scratch_id,  // temporary work space
       d_phase_refine_op);

   if (d_model_parameters.with_orientation()) {
      d_curr_to_scr_refine_alg->registerRefine(
          d_quat_scratch_id,  // destination
          d_quat_id,          // source
          d_quat_scratch_id,  // temporary work space
          d_quat_refine_op);
   }

   if (d_model_parameters.with_concentration()) {
      d_curr_to_scr_refine_alg->registerRefine(
          d_conc_scratch_id,  // destination
          d_conc_id,          // source
          d_conc_scratch_id,  // temporary
          d_conc_refine_op);
   }
}


//=======================================================================

void QuatModel::initializeCoarseRefineOperators()
{
   tbox::pout << "QuatModel::InitializeOperators()" << std::endl;

   assert(d_temperature_id >= 0);
   assert(d_temperature_scratch_id >= 0);
   assert(d_grains);
   assert(d_grid_geometry);

   initializeRefineCoarsenAlgorithms();

   d_grains->initializeRefineCoarsenAlgorithms(d_grid_geometry,
                                               d_quat_coarsen_op);

   d_integrator->initializeCoarseRefineOperators(d_gridding_algorithm,
                                                 d_quat_refine_op,
                                                 d_quat_coarsen_op);
   if (d_model_parameters.evolveQuat())
      d_integrator_quat_only->initializeCoarseRefineOperators(
          d_gridding_algorithm, d_quat_refine_op, d_quat_coarsen_op);
}

//=======================================================================

void QuatModel::copyCurrentToScratch(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time,
    QuatRefinePatchStrategy* patch_strategy)
{
   // tbox::plog<<"QuatModel::copyCurrentToScratch()"<<endl;
   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ln++) {
      copyCurrentToScratch(hierarchy, ln, time, patch_strategy);
   }
}

void QuatModel::copyCurrentToScratch(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const int ln,
    const double time, QuatRefinePatchStrategy* patch_strategy)
{
   if (patch_strategy == d_all_refine_patch_strategy) {
      d_curr_to_scr_refine_sched[ln]->fillData(time);

   } else {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      std::shared_ptr<xfer::RefineSchedule> schedule(
          d_curr_to_scr_refine_alg->createSchedule(level, ln - 1, hierarchy,
                                                   patch_strategy));
      schedule->fillData(time);
   }
}

//=======================================================================

bool QuatModel::isSymmetryAware(void) { return d_symmetry_aware; }

//=======================================================================

void QuatModel::setupInitialDataLevel(void)
{
   assert(d_temperature_id >= 0);

   tbox::plog << "\nsetupInitialDataLevel()..." << std::endl;

   PFModel::setupInitialDataLevel();

   assert(d_initial_level);

   if (d_model_parameters.with_phase()) {
      if (!d_initial_level->checkAllocated(d_phase_id)) {
         d_initial_level->allocatePatchData(d_phase_id);
         d_initial_level->setTime(0.0, d_phase_id);
      }
   }

   if (d_model_parameters.with_third_phase()) {
      if (!d_initial_level->checkAllocated(d_eta_id)) {
         d_initial_level->allocatePatchData(d_eta_id);
         d_initial_level->setTime(0.0, d_eta_id);
      }
   }

   if (!d_initial_level->checkAllocated(d_temperature_id)) {
      d_initial_level->allocatePatchData(d_temperature_id);
      d_initial_level->setTime(0.0, d_temperature_id);
   }

   if (d_model_parameters.with_orientation()) {
      if (!d_initial_level->checkAllocated(d_quat_id)) {
         d_initial_level->allocatePatchData(d_quat_id);
         d_initial_level->setTime(0.0, d_quat_id);
      }
   }

   if (d_model_parameters.with_concentration()) {
      if (!d_initial_level->checkAllocated(d_conc_id)) {
         d_initial_level->allocatePatchData(d_conc_id);
         d_initial_level->setTime(0.0, d_conc_id);
      }
   }
}

//=======================================================================

void QuatModel::setupHierarchy(void) { PFModel::setupHierarchy(); }

//=======================================================================


void QuatModel::CreateIntegrator(std::shared_ptr<tbox::Database> input_db)
{
   tbox::plog << "QuatModel::CreateIntegrator()" << std::endl;

   std::shared_ptr<tbox::Database> model_db =
       input_db->getDatabase("ModelParameters");

   std::string time_integration = "unsplit";

   time_integration =
       model_db->getStringWithDefault("time_integration", "unsplit");

   if (time_integration == "unsplit") {

      d_integrator.reset(QuatIntegratorFactory::create(
          "Integrator", d_model_parameters, this, d_grid_geometry, d_qlen,
          d_ncompositions, input_db, d_use_warm_start, d_symmetry_aware,
          d_all_periodic));
      if (d_model_parameters.evolveQuat())
         d_integrator_quat_only.reset(QuatIntegratorFactory::create(
             "quatonly", d_model_parameters, this, d_grid_geometry, d_qlen, 0,
             input_db, false, d_symmetry_aware, d_all_periodic));

   } else {
      tbox::pout << time_integration << std::endl;
      TBOX_ERROR("Invalid time_integration" << std::endl);
   }

   d_integrator->setVerbosity(d_verbosity->basicLevel());
   if (d_model_parameters.evolveQuat())
      d_integrator_quat_only->setVerbosity(d_verbosity->basicLevel());

   setIntegratorModelParameters();
}

void QuatModel::registerPhaseConcentrationVariables()
{
   assert(d_nghosts_aux_conc >= 0);

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   std::shared_ptr<hier::VariableContext> current =
       variable_db->getContext("CURRENT");

   if (!d_conc_l_var)
      d_conc_l_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                        "conc_l",
                                                        d_ncompositions));
   assert(d_conc_l_var);

   if (!d_conc_a_var)
      d_conc_a_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                        "conc_a",
                                                        d_ncompositions));
   assert(d_conc_a_var);

   if (d_model_parameters.withPhaseB()) {
      if (!d_conc_b_var)
         d_conc_b_var.reset(
             new pdat::CellVariable<double>(tbox::Dimension(NDIM), "conc_b",
                                            d_ncompositions));
      assert(d_conc_b_var);
   }

   if (d_model_parameters.isConcentrationModelCALPHAD()) {
      d_conc_l_ref_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM), "conc_l_ref",
                                         d_ncompositions));
      assert(d_conc_l_ref_var);
      d_conc_a_ref_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM), "conc_a_ref",
                                         d_ncompositions));
      assert(d_conc_a_ref_var);
      d_conc_l_ref_id = variable_db->registerVariableAndContext(
          d_conc_l_ref_var, current,
          hier::IntVector(tbox::Dimension(NDIM), d_nghosts_aux_conc));
      d_conc_a_ref_id = variable_db->registerVariableAndContext(
          d_conc_a_ref_var, current,
          hier::IntVector(tbox::Dimension(NDIM), d_nghosts_aux_conc));
      assert(d_conc_l_ref_id >= 0);
      assert(d_conc_a_ref_id >= 0);
      if (d_model_parameters.withPhaseB()) {
         d_conc_b_ref_var.reset(
             new pdat::CellVariable<double>(tbox::Dimension(NDIM), "conc_b_ref",
                                            d_ncompositions));
         assert(d_conc_b_ref_var);
         d_conc_b_ref_id = variable_db->registerVariableAndContext(
             d_conc_b_ref_var, current,
             hier::IntVector(tbox::Dimension(NDIM), d_nghosts_aux_conc));
         assert(d_conc_b_ref_id >= 0);
      }
   }

   d_conc_l_id = variable_db->registerVariableAndContext(
       d_conc_l_var, current,
       hier::IntVector(tbox::Dimension(NDIM), d_nghosts_aux_conc));
   assert(d_conc_l_id >= 0);

   d_conc_a_id = variable_db->registerVariableAndContext(
       d_conc_a_var, current,
       hier::IntVector(tbox::Dimension(NDIM), d_nghosts_aux_conc));
   assert(d_conc_a_id >= 0);

   if (d_model_parameters.withPhaseB()) {
      assert(d_conc_b_var);
      d_conc_b_id = variable_db->registerVariableAndContext(
          d_conc_b_var, current,
          hier::IntVector(tbox::Dimension(NDIM), d_nghosts_aux_conc));
      assert(d_conc_b_id >= 0);
   }
}

void QuatModel::registerPhaseConcentrationDiffusionVariables()
{
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   std::shared_ptr<hier::VariableContext> current =
       variable_db->getContext("CURRENT");

   // "isotropic" stencil needs some ghost values for diffusion field
   const int nghosts = d_model_parameters.useEBSisotropicStencil() ? 1 : 0;
   d_conc_pfm_diffusion_l_var.reset(
       new pdat::SideVariable<double>(tbox::Dimension(NDIM),
                                      "conc_pfm_diffusion_l",
                                      d_ncompositions * d_ncompositions));
   assert(d_conc_pfm_diffusion_l_var);
   d_conc_pfm_diffusion_l_id = variable_db->registerVariableAndContext(
       d_conc_pfm_diffusion_l_var, current,
       hier::IntVector(tbox::Dimension(NDIM), nghosts));
   assert(d_conc_pfm_diffusion_l_id >= 0);

   d_conc_pfm_diffusion_a_var.reset(
       new pdat::SideVariable<double>(tbox::Dimension(NDIM),
                                      "conc_pfm_diffusion_a",
                                      d_ncompositions * d_ncompositions));
   assert(d_conc_pfm_diffusion_a_var);
   d_conc_pfm_diffusion_a_id = variable_db->registerVariableAndContext(
       d_conc_pfm_diffusion_a_var, current,
       hier::IntVector(tbox::Dimension(NDIM), nghosts));
   assert(d_conc_pfm_diffusion_a_id >= 0);

   if (d_model_parameters.withPhaseB()) {
      d_conc_pfm_diffusion_b_var.reset(
          new pdat::SideVariable<double>(tbox::Dimension(NDIM),
                                         "conc_pfm_diffusion_b",
                                         d_ncompositions * d_ncompositions));
      assert(d_conc_pfm_diffusion_b_var);
      d_conc_pfm_diffusion_b_id = variable_db->registerVariableAndContext(
          d_conc_pfm_diffusion_b_var, current,
          hier::IntVector(tbox::Dimension(NDIM), nghosts));
      assert(d_conc_pfm_diffusion_b_id >= 0);
   }

   d_conc_diffusion_coeff_l_var.reset(
       new pdat::SideVariable<double>(tbox::Dimension(NDIM),
                                      "conc_diffusion_coeff_l",
                                      d_ncompositions * d_ncompositions));
   assert(d_conc_diffusion_coeff_l_var);
   d_conc_diffusion_coeff_l_id = variable_db->registerVariableAndContext(
       d_conc_diffusion_coeff_l_var, current,
       hier::IntVector(tbox::Dimension(NDIM), 0));
   assert(d_conc_diffusion_coeff_l_id >= 0);

   d_conc_diffusion_coeff_a_var.reset(
       new pdat::SideVariable<double>(tbox::Dimension(NDIM),
                                      "conc_diffusion_coeff_a",
                                      d_ncompositions * d_ncompositions));
   assert(d_conc_diffusion_coeff_a_var);
   d_conc_diffusion_coeff_a_id = variable_db->registerVariableAndContext(
       d_conc_diffusion_coeff_a_var, current,
       hier::IntVector(tbox::Dimension(NDIM), 0));
   assert(d_conc_diffusion_coeff_a_id >= 0);


   if (d_model_parameters.with_extra_visit_output()) {
      std::cout << "d_model_parameters.concRHSstrategyIsEBS()" << std::endl;
      d_visit_conc_pfm_diffusion_l_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                         "visit_conc_pfm_diffusion_l",
                                         d_ncompositions * d_ncompositions));
      assert(d_visit_conc_pfm_diffusion_l_var);
      d_visit_conc_pfm_diffusion_l_id = variable_db->registerVariableAndContext(
          d_visit_conc_pfm_diffusion_l_var, current,
          hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(d_visit_conc_pfm_diffusion_l_id >= 0);

      d_visit_conc_pfm_diffusion_a_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                         "visit_conc_pfm_diffusion_a",
                                         d_ncompositions * d_ncompositions));
      assert(d_visit_conc_pfm_diffusion_a_var);
      d_visit_conc_pfm_diffusion_a_id = variable_db->registerVariableAndContext(
          d_visit_conc_pfm_diffusion_a_var, current,
          hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(d_visit_conc_pfm_diffusion_a_id >= 0);

      if (d_model_parameters.withPhaseB()) {
         d_visit_conc_pfm_diffusion_b_var.reset(
             new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                            "visit_conc_pfm_diffusion_b",
                                            d_ncompositions * d_ncompositions));
         assert(d_visit_conc_pfm_diffusion_b_var);
         d_visit_conc_pfm_diffusion_b_id =
             variable_db->registerVariableAndContext(
                 d_visit_conc_pfm_diffusion_b_var, current,
                 hier::IntVector(tbox::Dimension(NDIM), 0));
         assert(d_visit_conc_pfm_diffusion_b_id >= 0);
      }
   }
}

void QuatModel::registerConcentrationVariables(void)
{
   tbox::plog << "QuatModel::registerConcentrationVariables()" << std::endl;
   assert(d_nghosts > 0);

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   std::shared_ptr<hier::VariableContext> current =
       variable_db->getContext("CURRENT");
   std::shared_ptr<hier::VariableContext> scratch =
       variable_db->getContext("SCRATCH");
   d_conc_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                   "concentration",
                                                   d_ncompositions));
   assert(d_conc_var);
   d_conc_id = variable_db->registerVariableAndContext(
       d_conc_var, current, hier::IntVector(tbox::Dimension(NDIM), 0));
   d_conc_scratch_id = variable_db->registerVariableAndContext(
       d_conc_var, scratch, hier::IntVector(tbox::Dimension(NDIM), d_nghosts));
   assert(d_conc_id >= 0);
   assert(d_conc_scratch_id >= 0);

   {
      // no need for D, just need D_0
      d_conc_diffusion_var.reset();
   }

   for (int ic = 0; ic < d_ncompositions; ic++) {
      std::shared_ptr<pdat::SideVariable<double> > conc_pfm_diffusion_var;
      conc_pfm_diffusion_var.reset(new pdat::SideVariable<double>(
          tbox::Dimension(NDIM), "conc_pfm_diffusion" + std::to_string(ic)));
      assert(conc_pfm_diffusion_var);
      d_conc_pfm_diffusion_var.push_back(conc_pfm_diffusion_var);
      d_conc_pfm_diffusion_id.push_back(variable_db->registerVariableAndContext(
          d_conc_pfm_diffusion_var[ic], current,
          hier::IntVector(tbox::Dimension(NDIM), d_nghosts - 1)));
      assert(d_conc_pfm_diffusion_id[ic] >= 0);
   }

   d_model_parameters.checkValidityConcRHSstrategy();

   if (d_model_parameters.concRHSstrategyIsKKS() ||
       d_model_parameters.concRHSstrategyIsCahnHilliard() ||
       d_model_parameters.concRHSstrategyIsWangSintering()) {
      if (d_model_parameters.with_extra_visit_output()) {
         d_visit_conc_pfm_diffusion_var.reset(
             new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                            "visit_conc_pfm_diffusion",
                                            d_ncompositions * d_ncompositions));
         assert(d_visit_conc_pfm_diffusion_var);
         d_visit_conc_pfm_diffusion_id =
             variable_db->registerVariableAndContext(
                 d_visit_conc_pfm_diffusion_var, current,
                 hier::IntVector(tbox::Dimension(NDIM), 0));
         assert(d_visit_conc_pfm_diffusion_id >= 0);
      }
   }


   if (d_model_parameters.concRHSstrategyIsKKS() ||
       d_model_parameters.concRHSstrategyIsBeckermann()) {
      d_conc_phase_coupling_diffusion_var.reset(
          new pdat::SideVariable<double>(tbox::Dimension(NDIM),
                                         "conc_phase_coupling_diffusion"));
      assert(d_conc_phase_coupling_diffusion_var);
      d_conc_phase_coupling_diffusion_id =
          variable_db->registerVariableAndContext(
              d_conc_phase_coupling_diffusion_var, current,
              hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(d_conc_phase_coupling_diffusion_id >= 0);
   } else {
      {
         if (d_model_parameters.concentrationModelNeedsPhaseConcentrations())
            registerPhaseConcentrationDiffusionVariables();
      }

      if (d_model_parameters.with_gradT()) {
         assert(d_ncompositions > 0);
         d_conc_Mq_var.reset(
             new pdat::SideVariable<double>(tbox::Dimension(NDIM), "conc_Mq",
                                            d_ncompositions));
         assert(d_conc_Mq_var);
         d_conc_Mq_id = variable_db->registerVariableAndContext(
             d_conc_Mq_var, current, hier::IntVector(tbox::Dimension(NDIM), 0));
         assert(d_conc_Mq_id >= 0);
      }

      if (d_model_parameters.withPhaseB()) {
         d_conc_diffusion_coeff_b_var.reset(
             new pdat::SideVariable<double>(tbox::Dimension(NDIM),
                                            "conc_diffusion_coeff_b",
                                            d_ncompositions * d_ncompositions));
         assert(d_conc_diffusion_coeff_b_var);
         d_conc_diffusion_coeff_b_id = variable_db->registerVariableAndContext(
             d_conc_diffusion_coeff_b_var, current,
             hier::IntVector(tbox::Dimension(NDIM), 0));
         assert(d_conc_diffusion_coeff_b_id >= 0);
      }
   }

   if (d_model_parameters.concentrationModelNeedsPhaseConcentrations()) {
      assert(d_nghosts_aux_conc >= 0);

      registerPhaseConcentrationVariables();

   }  // if d_model_parameters.concentrationModelNeedsPhaseConcentrations()

   if (d_model_parameters.with_third_phase()) {
      d_conc_eta_coupling_diffusion_var.reset(
          new pdat::SideVariable<double>(tbox::Dimension(NDIM),
                                         "conc_eta_coupling_diffusion"));
      assert(d_conc_eta_coupling_diffusion_var);
      d_conc_eta_coupling_diffusion_id =
          variable_db->registerVariableAndContext(
              d_conc_eta_coupling_diffusion_var, current,
              hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(d_conc_eta_coupling_diffusion_id >= 0);
   }

   if (d_model_parameters.with_partition_coeff()) {
      d_partition_coeff_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                         "partition_coeff", 1));
      d_partition_coeff_id = variable_db->registerVariableAndContext(
          d_partition_coeff_var, current,
          hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(d_partition_coeff_id >= 0);

      if (d_model_parameters.needGhosts4PartitionCoeff()) {
         d_partition_coeff_scratch_id = variable_db->registerVariableAndContext(
             d_partition_coeff_var, scratch,
             hier::IntVector(tbox::Dimension(NDIM), 1));
         assert(d_partition_coeff_scratch_id >= 0);
      }
   }
   if (d_model_parameters.with_velocity()) {
      d_velocity_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM), "velocity", 1));
      d_velocity_id = variable_db->registerVariableAndContext(
          d_velocity_var, current, hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(d_velocity_id >= 0);
   }

   if (d_model_parameters.with_concentration())
      d_integrator->RegisterConcentrationVariables(
          d_conc_var, d_conc_pfm_diffusion_var,
          d_conc_phase_coupling_diffusion_var,
          d_conc_eta_coupling_diffusion_var, d_conc_diffusion_var);
}

//=======================================================================

void QuatModel::registerEtaVariables(void)
{
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   std::shared_ptr<hier::VariableContext> current =
       variable_db->getContext("CURRENT");
   std::shared_ptr<hier::VariableContext> scratch =
       variable_db->getContext("SCRATCH");

   d_eta_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                  "et"
                                                  "a"));
   assert(d_eta_var);
   d_eta_id = variable_db->registerVariableAndContext(
       d_eta_var, current, hier::IntVector(tbox::Dimension(NDIM), 0));
   d_eta_scratch_id = variable_db->registerVariableAndContext(
       d_eta_var, scratch, hier::IntVector(tbox::Dimension(NDIM), d_nghosts));

   d_eta_diffs_var.reset(
       new pdat::SideVariable<double>(tbox::Dimension(NDIM), "eta_diffs"));
   assert(d_eta_diffs_var);
   d_eta_diffs_id = variable_db->registerVariableAndContext(
       d_eta_diffs_var, current, hier::IntVector(tbox::Dimension(NDIM), 1));

   d_eta_grad_cell_var.reset(
       new pdat::CellVariable<double>(tbox::Dimension(NDIM), "eta_grad_cell",
                                      NDIM));
   assert(d_eta_grad_cell_var);
   d_eta_grad_cell_id = variable_db->registerVariableAndContext(
       d_eta_grad_cell_var, current, hier::IntVector(tbox::Dimension(NDIM), 0));

   d_eta_grad_side_var.reset(
       new pdat::SideVariable<double>(tbox::Dimension(NDIM), "eta_grad_side",
                                      NDIM));
   assert(d_eta_grad_side_var);
   d_eta_grad_side_id = variable_db->registerVariableAndContext(
       d_eta_grad_side_var, current, hier::IntVector(tbox::Dimension(NDIM), 0));

   d_eta_mobility_var.reset(
       new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                      "eta_"
                                      "mobility"));
   assert(d_eta_mobility_var);
   d_eta_mobility_id = variable_db->registerVariableAndContext(
       d_eta_mobility_var, current, hier::IntVector(tbox::Dimension(NDIM), 1));
}

//=======================================================================

void QuatModel::registerPhaseVariables(void)
{
   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   std::shared_ptr<hier::VariableContext> current =
       variable_db->getContext("CURRENT");
   std::shared_ptr<hier::VariableContext> scratch =
       variable_db->getContext("SCRATCH");

   const int nphases = d_model_parameters.norderp();
   d_phase_var.reset(
       new pdat::CellVariable<double>(tbox::Dimension(NDIM), "phase", nphases));
   assert(d_phase_var);
   d_phase_id = variable_db->registerVariableAndContext(
       d_phase_var, current, hier::IntVector(tbox::Dimension(NDIM), 0));
   d_phase_scratch_id = variable_db->registerVariableAndContext(
       d_phase_var, scratch, hier::IntVector(tbox::Dimension(NDIM), d_nghosts));

   d_phase_diffs_var.reset(
       new pdat::SideVariable<double>(tbox::Dimension(NDIM), "phase_diffs"));
   assert(d_phase_diffs_var);
   d_phase_diffs_id = variable_db->registerVariableAndContext(
       d_phase_diffs_var, current, hier::IntVector(tbox::Dimension(NDIM), 1));

   if (d_model_parameters.with_extra_visit_output()) {
      d_phase_diffs_cell_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                         "phase_diffs_cell", NDIM));
      assert(d_phase_diffs_cell_var);
      d_phase_diffs_cell_id = variable_db->registerVariableAndContext(
          d_phase_diffs_cell_var, current,
          hier::IntVector(tbox::Dimension(NDIM), 0));
   }

   d_phase_grad_cell_var.reset(
       new pdat::CellVariable<double>(tbox::Dimension(NDIM), "phase_grad_cell",
                                      NDIM));
   assert(d_phase_grad_cell_var);
   d_phase_grad_cell_id = variable_db->registerVariableAndContext(
       d_phase_grad_cell_var, current,
       hier::IntVector(tbox::Dimension(NDIM), 0));

   d_phase_grad_side_var.reset(
       new pdat::SideVariable<double>(tbox::Dimension(NDIM), "phase_grad_side",
                                      NDIM));
   assert(d_phase_grad_side_var);
   d_phase_grad_side_id = variable_db->registerVariableAndContext(
       d_phase_grad_side_var, current,
       hier::IntVector(tbox::Dimension(NDIM), 0));

   d_phase_mobility_var.reset(
       new pdat::CellVariable<double>(tbox::Dimension(NDIM), "phase_mobility",
                                      1));
   assert(d_phase_mobility_var);
   d_phase_mobility_id = variable_db->registerVariableAndContext(
       d_phase_mobility_var, current,
       hier::IntVector(tbox::Dimension(NDIM), 1));
}

//=======================================================================

void QuatModel::registerOrientationVariables(void)
{
   assert(d_qlen > 0);

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   std::shared_ptr<hier::VariableContext> current =
       variable_db->getContext("CURRENT");
   std::shared_ptr<hier::VariableContext> scratch =
       variable_db->getContext("SCRATCH");

   const int symm_depth = isSymmetryAware() ? d_qlen * 2 : d_qlen * 1;

   d_quat_var.reset(
       new pdat::CellVariable<double>(tbox::Dimension(NDIM), "quat", d_qlen));
   assert(d_quat_var);
   d_quat_id = variable_db->registerVariableAndContext(
       d_quat_var, current, hier::IntVector(tbox::Dimension(NDIM), 0));
   d_quat_scratch_id = variable_db->registerVariableAndContext(
       d_quat_var, scratch, hier::IntVector(tbox::Dimension(NDIM), d_nghosts));

   if (d_model_parameters.evolveQuat()) {
      d_quat_relax_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM), "quat_relax",
                                         d_qlen));
      assert(d_quat_relax_var);
      d_quat_relax_id = variable_db->registerVariableAndContext(
          d_quat_relax_var, current, hier::IntVector(tbox::Dimension(NDIM), 0));
      d_quat_relax_scratch_id = variable_db->registerVariableAndContext(
          d_quat_relax_var, scratch,
          hier::IntVector(tbox::Dimension(NDIM), d_nghosts));

      d_quat_diffs_var.reset(
          new pdat::SideVariable<double>(tbox::Dimension(NDIM), "quat_diffs",
                                         symm_depth));
      assert(d_quat_diffs_var);
      d_quat_diffs_id = variable_db->registerVariableAndContext(
          d_quat_diffs_var, current, hier::IntVector(tbox::Dimension(NDIM), 1));

      d_quat_grad_cell_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                         "quat_grad_cell", NDIM * d_qlen));
      assert(d_quat_grad_cell_var);
      d_quat_grad_cell_id = variable_db->registerVariableAndContext(
          d_quat_grad_cell_var, current,
          hier::IntVector(tbox::Dimension(NDIM), 0));

      d_quat_grad_side_var.reset(
          new pdat::SideVariable<double>(tbox::Dimension(NDIM),
                                         "quat_grad_side", NDIM * d_qlen));
      assert(d_quat_grad_side_var);
      d_quat_grad_side_id = variable_db->registerVariableAndContext(
          d_quat_grad_side_var, current,
          hier::IntVector(tbox::Dimension(NDIM), 0));

      if (d_model_parameters.with_extra_visit_output()) {
         d_quat_diffs_cell_var.reset(
             new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                            "quat_diffs_cell", d_qlen * NDIM));
         assert(d_quat_diffs_cell_var);
         d_quat_diffs_cell_id = variable_db->registerVariableAndContext(
             d_quat_diffs_cell_var, current,
             hier::IntVector(tbox::Dimension(NDIM), 0));

         d_quat_nonsymm_diffs_cell_var.reset(new pdat::CellVariable<double>(
             tbox::Dimension(NDIM), "quat_nonsymm_diffs_cell", d_qlen * NDIM));
         assert(d_quat_nonsymm_diffs_cell_var);
         d_quat_nonsymm_diffs_cell_id = variable_db->registerVariableAndContext(
             d_quat_nonsymm_diffs_cell_var, current,
             hier::IntVector(tbox::Dimension(NDIM), 0));

         d_quat_norm_error_var.reset(
             new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                            "quat_norm_error"));
         assert(d_quat_norm_error_var);
         d_quat_norm_error_id = variable_db->registerVariableAndContext(
             d_quat_norm_error_var, current,
             hier::IntVector(tbox::Dimension(NDIM), 0));

         if (d_symmetry_aware) {
            d_quat_symm_rotation_cell_var.reset(
                new pdat::CellVariable<int>(tbox::Dimension(NDIM),
                                            "quat_symm_rotation_cell", NDIM));
            assert(d_quat_symm_rotation_cell_var);
            d_quat_symm_rotation_cell_id =
                variable_db->registerVariableAndContext(
                    d_quat_symm_rotation_cell_var, current,
                    hier::IntVector(tbox::Dimension(NDIM), 0));
         }
      }

      if (d_symmetry_aware) {
         d_quat_symm_rotation_var.reset(
             new pdat::SideVariable<int>(tbox::Dimension(NDIM),
                                         "quat_symm_rotation"));
         assert(d_quat_symm_rotation_var);
         d_quat_symm_rotation_id = variable_db->registerVariableAndContext(
             d_quat_symm_rotation_var, current,
             hier::IntVector(tbox::Dimension(NDIM), 1));
      }

      d_quat_grad_modulus_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                         "quat_grad_modulus"));
      assert(d_quat_grad_modulus_var);
      d_quat_grad_modulus_id = variable_db->registerVariableAndContext(
          d_quat_grad_modulus_var, current,
          hier::IntVector(tbox::Dimension(NDIM), 0));

      // we have ghost values for quat_mobility so that values are defined
      // at physical boundaries when needed to compute sqrt_mobility
      d_quat_mobility_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                         "quat_mobility"));
      assert(d_quat_mobility_var);
      d_quat_mobility_id = variable_db->registerVariableAndContext(
          d_quat_mobility_var, current,
          hier::IntVector(tbox::Dimension(NDIM), 1));
   }
}

//=======================================================================

void QuatModel::registerPatchDataForRestart(void)
{
   if (d_model_parameters.with_phase()) {
      hier::PatchDataRestartManager::getManager()->registerPatchDataForRestart(
          d_phase_id);
   }
   if (d_model_parameters.with_third_phase()) {
      hier::PatchDataRestartManager::getManager()->registerPatchDataForRestart(
          d_eta_id);
   }

   hier::PatchDataRestartManager::getManager()->registerPatchDataForRestart(
       d_temperature_id);

   if (d_model_parameters.with_orientation()) {
      hier::PatchDataRestartManager::getManager()->registerPatchDataForRestart(
          d_quat_id);
   }
   if (d_model_parameters.with_concentration()) {
      hier::PatchDataRestartManager::getManager()->registerPatchDataForRestart(
          d_conc_id);
      if (d_model_parameters.concentrationModelNeedsPhaseConcentrations()) {
         hier::PatchDataRestartManager::getManager()
             ->registerPatchDataForRestart(d_conc_l_id);
         hier::PatchDataRestartManager::getManager()
             ->registerPatchDataForRestart(d_conc_a_id);
         if (d_conc_b_id > -1)
            hier::PatchDataRestartManager::getManager()
                ->registerPatchDataForRestart(d_conc_b_id);
      }
   }
}

//=======================================================================

void QuatModel::RegisterVariables(void)
{
   tbox::pout << "QuatModel::RegisterVariables()" << std::endl;

   assert(d_grains);

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   std::shared_ptr<hier::VariableContext> current =
       variable_db->getContext("CURRENT");
   std::shared_ptr<hier::VariableContext> scratch =
       variable_db->getContext("SCRATCH");

   assert(current);
   assert(scratch);

   if (d_model_parameters.with_phase()) {
      registerPhaseVariables();
   }
   if (d_model_parameters.with_third_phase()) {
      registerEtaVariables();
   }

   d_temperature_var.reset(
       new pdat::CellVariable<double>(tbox::Dimension(NDIM), "temperature"));
   assert(d_temperature_var);
   d_temperature_id = variable_db->registerVariableAndContext(
       d_temperature_var, current, hier::IntVector(tbox::Dimension(NDIM), 0));
   d_temperature_scratch_id = variable_db->registerVariableAndContext(
       d_temperature_var, scratch,
       hier::IntVector(tbox::Dimension(NDIM), d_nghosts));
   if (d_model_parameters.with_steady_temperature()) {
      d_temperature_rhs_steady_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                         "temperature_rhs_", 1));
      d_temperature_rhs_steady_id = variable_db->registerVariableAndContext(
          d_temperature_rhs_steady_var, current,
          hier::IntVector(tbox::Dimension(NDIM), 1));
      assert(d_temperature_rhs_steady_id >= 0);
   }

   if (d_model_parameters.with_bias_well()) {
      d_equilibrium_temperature_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                         "equilibrium_temperature", 1));

      d_equilibrium_temperature_id = variable_db->registerVariableAndContext(
          d_equilibrium_temperature_var, current,
          hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(d_equilibrium_temperature_id >= 0);
   }

   if (d_model_parameters.with_orientation()) {
      registerOrientationVariables();
   }

   d_weight_var.reset(
       new pdat::CellVariable<double>(tbox::Dimension(NDIM), "weight"));
   assert(d_weight_var);
   d_weight_id = variable_db->registerVariableAndContext(
       d_weight_var, current, hier::IntVector(tbox::Dimension(NDIM), 0));

   if (d_subdomain_diagnostics[0] > 0.) {
      d_weight_diagnostics_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                         "weight_diagnostics"));
      assert(d_weight_diagnostics_var);
      d_weight_diagnostics_id = variable_db->registerVariableAndContext(
          d_weight_diagnostics_var, current,
          hier::IntVector(tbox::Dimension(NDIM), 0));
   }

   d_work_var.reset(
       new pdat::CellVariable<double>(tbox::Dimension(NDIM), "work"));
   assert(d_work_var);
   d_work_id = variable_db->registerVariableAndContext(
       d_work_var, current, hier::IntVector(tbox::Dimension(NDIM), 0));

   if (d_model_parameters.with_heat_equation()) {
      d_cp_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM), "cp"));
      assert(d_cp_var);
      // ghost value needed in particular for mobility of off-diagonal
      // phase-temperature block in preconditioner
      d_cp_id = variable_db->registerVariableAndContext(
          d_cp_var, current, hier::IntVector(tbox::Dimension(NDIM), 1));
      assert(d_cp_id >= 0);
   }

   d_integrator->RegisterVariables(d_phase_var, d_eta_var, d_quat_var,
                                   d_quat_grad_cell_var, d_quat_grad_side_var,
                                   d_quat_grad_modulus_var,
                                   d_phase_mobility_var, d_eta_mobility_var,
                                   d_quat_mobility_var, d_quat_diffs_var,
                                   d_quat_symm_rotation_var, d_weight_var,
                                   d_temperature_var, d_cp_var);

   if (d_model_parameters.evolveQuat()) {
      d_integrator_quat_only->RegisterVariables(
          d_phase_var, d_eta_var, d_quat_relax_var, d_quat_grad_cell_var,
          d_quat_grad_side_var, d_quat_grad_modulus_var, d_phase_mobility_var,
          d_eta_mobility_var, d_quat_mobility_var, d_quat_diffs_var,
          d_quat_symm_rotation_var, d_weight_var, d_temperature_var, d_cp_var);
   }

   //
   // free energy variables
   //
   d_f_l_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                  "f_"
                                                  "l"));
   assert(d_f_l_var);
   d_f_l_id = variable_db->registerVariableAndContext(
       d_f_l_var, current, hier::IntVector(tbox::Dimension(NDIM), 0));
   assert(d_f_l_id >= 0);

   d_f_a_var.reset(new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                                  "f_"
                                                  "a"));
   assert(d_f_a_var);
   d_f_a_id = variable_db->registerVariableAndContext(
       d_f_a_var, current, hier::IntVector(tbox::Dimension(NDIM), 0));
   assert(d_f_a_id >= 0);

   if (d_model_parameters.withPhaseB()) {
      d_f_b_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM), "f_b"));
      assert(d_f_b_var);
      d_f_b_id = variable_db->registerVariableAndContext(
          d_f_b_var, current, hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(d_f_b_id >= 0);
   }

   d_integrator->RegisterFreeEnergyVariables(d_f_l_var, d_f_a_var, d_f_b_var);
   if (d_model_parameters.evolveQuat())
      d_integrator_quat_only->RegisterFreeEnergyVariables(d_f_l_var, d_f_a_var,
                                                          d_f_b_var);

   // now create and register concentration variables
   if (d_model_parameters.with_concentration()) {
      registerConcentrationVariables();
   }

   registerPatchDataForRestart();

   if (d_model_parameters.with_visit_energy_output()) {
      d_energy_diag_var.reset(
          new pdat::CellVariable<double>(tbox::Dimension(NDIM), "energy_diag"));
      assert(d_energy_diag_var);
      d_energy_diag_id = variable_db->registerVariableAndContext(
          d_energy_diag_var, current,
          hier::IntVector(tbox::Dimension(NDIM), 0));
   }

   d_grains->registerVariables();

   d_integrator->setupBC();

   if (d_model_parameters.evolveQuat()) {
      d_integrator_quat_only->setupBC();
   }
}

//=======================================================================

void QuatModel::RegisterWithVisit(void)
{
   assert(d_visit_data_writer);

   if (d_model_parameters.with_phase()) {
      const int n = d_model_parameters.norderp();
      // dump at most 32 order parameters
      const int i0 = std::max(0, n - 32);
      for (int i = i0; i < n; i++) {
         std::string visit_name("phase" + std::to_string(i));
         d_visit_data_writer->registerPlotQuantity(visit_name, "SCALAR",
                                                   d_phase_id, i);
      }
   }
   if (d_model_parameters.with_third_phase()) {
      assert(d_eta_id >= 0);
      d_visit_data_writer->registerPlotQuantity("eta", "SCALAR", d_eta_id, 0);
   }

   if (!d_model_parameters
            .isTemperatureUniform()) {  // if SCALAR, then constant in space
      assert(d_temperature_id >= 0);
      d_visit_data_writer->registerPlotQuantity("temperature", "SCALAR",
                                                d_temperature_id, 0);
   }

   if (d_model_parameters.with_orientation() &&
       d_model_parameters.evolveQuat()) {
      assert(d_quat_id >= 0);
      for (int n = 0; n < d_qlen; n++) {
         std::string visit_name("q" + std::to_string(n));

         d_visit_data_writer->registerPlotQuantity(visit_name, "SCALAR",
                                                   d_quat_id, n);
      }
   }

   if (d_model_parameters.with_concentration()) {
      assert(d_conc_id >= 0);
      for (int n = 0; n < d_ncompositions; n++) {
         std::string visit_name("concentration" + std::to_string(n));
         d_visit_data_writer->registerPlotQuantity(visit_name, "SCALAR",
                                                   d_conc_id, n);

         if (d_model_parameters.with_extra_visit_output()) {
            if (d_model_parameters
                    .concentrationModelNeedsPhaseConcentrations()) {
               assert(d_conc_l_id >= 0);
               assert(d_conc_a_id >= 0);
               std::string visit_namel("conc_l" + std::to_string(n));
               d_visit_data_writer->registerPlotQuantity(visit_namel, "SCALAR",
                                                         d_conc_l_id, n);
               std::string visit_namea("conc_a" + std::to_string(n));
               d_visit_data_writer->registerPlotQuantity(visit_namea, "SCALAR",
                                                         d_conc_a_id, n);
               if (d_model_parameters.withPhaseB()) {
                  assert(d_conc_b_id >= 0);
                  std::string visit_nameb("conc_b" + std::to_string(n));
                  d_visit_data_writer->registerPlotQuantity(visit_nameb,
                                                            "SCALAR",
                                                            d_conc_b_id, n);
               }
            }
            // diffusion
            if (d_model_parameters.concRHSstrategyIsEBS()) {
               assert(d_visit_conc_pfm_diffusion_l_id >= 0);
               assert(d_visit_conc_pfm_diffusion_a_id >= 0);
               std::string vnamel("conc_pfm_diffusion_l" + std::to_string(n));
               d_visit_data_writer->registerPlotQuantity(
                   vnamel, "SCALAR", d_visit_conc_pfm_diffusion_l_id, n);
               std::string vnamea("conc_pfm_diffusion_a" + std::to_string(n));
               d_visit_data_writer->registerPlotQuantity(
                   vnamea, "SCALAR", d_visit_conc_pfm_diffusion_a_id, n);
               if (d_model_parameters.withPhaseB()) {
                  assert(d_conc_pfm_diffusion_b_id >= 0);
                  std::string vnameb("conc_pfm_diffusion_b" +
                                     std::to_string(n));
                  d_visit_data_writer->registerPlotQuantity(
                      vnameb, "SCALAR", d_visit_conc_pfm_diffusion_b_id, n);
               }
            } else {
               assert(d_visit_conc_pfm_diffusion_id >= 0);
               std::string vnamel("conc_pfm_diffusion" + std::to_string(n));
               d_visit_data_writer->registerPlotQuantity(
                   vnamel, "SCALAR", d_visit_conc_pfm_diffusion_id, n);
            }
         }
      }
   }

   if (d_model_parameters.with_velocity()) {
      assert(d_velocity_id >= 0);
      d_visit_data_writer->registerPlotQuantity("velocity", "SCALAR",
                                                d_velocity_id, 0);
   }

   if (d_model_parameters.with_extra_visit_output()) {
      if (d_phase_diffs_cell_id >= 0) {
         for (int d = 0; d < NDIM; d++) {
            std::string visit_name("phase_diffs_cell" + std::to_string(d));
            d_visit_data_writer->registerPlotQuantity(visit_name, "SCALAR",
                                                      d_phase_diffs_cell_id, d);
         }
      }

      if (d_phase_grad_cell_id >= 0) {
         for (int d = 0; d < NDIM; d++) {
            std::string visit_name("phase_grad_cell" + std::to_string(d));
            d_visit_data_writer->registerPlotQuantity(visit_name, "SCALAR",
                                                      d_phase_grad_cell_id, d);
         }
      }

      if (d_quat_grad_cell_id >= 0) {
         for (int d = 0; d < NDIM; d++) {
            for (int n = 0; n < d_qlen; n++) {
               std::string visit_name("quat_grad_cell_d" + std::to_string(d) +
                                      "_q" + std::to_string(n));

               d_visit_data_writer->registerPlotQuantity(visit_name, "SCALAR",
                                                         d_quat_grad_cell_id,
                                                         d * d_qlen + n);
            }
         }
      }

      if (d_quat_diffs_cell_id >= 0) {
         for (int d = 0; d < NDIM; d++) {
            for (int n = 0; n < d_qlen; n++) {
               std::string visit_symm_name("quat_diffs_cell_d" +
                                           std::to_string(d) + "_q" +
                                           std::to_string(n));

               d_visit_data_writer->registerPlotQuantity(visit_symm_name,
                                                         "SCALAR",
                                                         d_quat_diffs_cell_id,
                                                         d * d_qlen + n);

               std::string visit_nonsymm_name("quat_nonsymm_diffs_cell_d" +
                                              std::to_string(d) + "_q" +
                                              std::to_string(n));

               d_visit_data_writer->registerPlotQuantity(
                   visit_nonsymm_name, "SCALAR", d_quat_nonsymm_diffs_cell_id,
                   d * d_qlen + n);
            }
         }
      }

      if (d_quat_grad_modulus_id >= 0) {
         d_visit_data_writer->registerPlotQuantity("quat_grad_modulus",
                                                   "SCALAR",
                                                   d_quat_grad_modulus_id, 0);
      }

      if (d_quat_norm_error_id >= 0) {
         d_visit_data_writer->registerPlotQuantity("quat_norm_error", "SCALAR",
                                                   d_quat_norm_error_id, 0);
      }

      if (d_phase_mobility_id >= 0) {
         d_visit_data_writer->registerPlotQuantity("phase_mobility", "SCALAR",
                                                   d_phase_mobility_id, 0);
      }

      if (d_eta_mobility_id >= 0) {
         d_visit_data_writer->registerPlotQuantity("eta_mobility", "SCALAR",
                                                   d_eta_mobility_id, 0);
      }

      if (d_quat_mobility_id >= 0) {
         d_visit_data_writer->registerPlotQuantity("quat_mobility", "SCALAR",
                                                   d_quat_mobility_id, 0);
      }

      if (d_model_parameters.with_orientation() && d_symmetry_aware) {

         for (int d = 0; d < NDIM; d++) {
            std::string visit_name("quat_symm_rotation_cell" +
                                   std::to_string(d));

            d_visit_data_writer->registerPlotQuantity(
                visit_name, "SCALAR", d_quat_symm_rotation_cell_id, d);
         }
      }

      if (d_model_parameters.with_partition_coeff()) {
         d_visit_data_writer->registerPlotQuantity("partition_coeff", "SCALAR",
                                                   d_partition_coeff_id, 0);
      }

      if (d_cp_id >= 0) {
         d_visit_data_writer->registerPlotQuantity("cp", "SCALAR", d_cp_id, 0);
      }
   }  // if ( d_model_parameters.with_extra_visit_output() )

   if (d_model_parameters.with_visit_grain_output() &&
       d_grain_diag_interval->isActive()) {
      d_grains->registerWithVisit(d_visit_data_writer);
   }

   if (d_model_parameters.with_visit_energy_output()) {
      d_visit_data_writer->registerPlotQuantity("energy", "SCALAR",
                                                d_energy_diag_id, 0);
   }

   d_integrator->RegisterWithVisit(d_visit_data_writer);
}

//=======================================================================

void QuatModel::Run(void)
{
   if (d_model_parameters.with_orientation() && d_symmetry_aware)
      computeSymmetryRotations(d_patch_hierarchy, d_time);

   return PFModel::Run();
}

//=======================================================================

void QuatModel::extendGrainOrientation(void)
{
   assert(d_grains);
   assert(d_quat_scratch_id >= 0);
   assert(d_quat_id >= 0);
   assert(d_phase_id >= 0);

   // Fill ghosts of original quat data
   copyCurrentToScratch(d_patch_hierarchy, d_time, d_all_refine_patch_strategy);

   d_grains->extendGrainOrientation(d_patch_hierarchy, d_time,
                                    d_quat_scratch_id, d_phase_id, d_quat_id);

   copyCurrentToScratch(d_patch_hierarchy, d_time, d_all_refine_patch_strategy);
}

//=======================================================================

bool QuatModel::resetGrains(void)
{
   assert(d_grains);
   assert(d_quat_relax_id >= 0);
   assert(d_integrator_quat_only);

   // test if number of grains has changed first...
   int old_number_of_grains = d_number_of_grains;
   d_number_of_grains = d_grains->getNumberOfGrains();

   // if the number of grains has not decreased, don't do anything
   if (old_number_of_grains <= d_number_of_grains &&
       !d_grain_extend_interval->hasIntervalPassed(d_cycle, d_time)) {
      return false;
   }

   t_resetGrains_timer->start();

   tbox::pout << "Old number of grains: " << old_number_of_grains << std::endl;
   tbox::pout << "New number of grains: " << d_number_of_grains << std::endl;

   double total_energy, interface_energy;
   double well_energy, free_energy, eta_energy;
   double original_orient_energy;
   double original_qint_energy;

   evaluateEnergy(d_patch_hierarchy, d_time, total_energy, interface_energy,
                  eta_energy, original_orient_energy, original_qint_energy,
                  well_energy, free_energy);

   if (d_extra_energy_detail) {
      tbox::pout << std::setprecision(8);
      tbox::pout << "  Total energy      = " << total_energy << std::endl;
      tbox::pout << "  interface energy  = " << interface_energy << std::endl;
      tbox::pout << "    orient energy   = " << original_orient_energy
                 << std::endl;
      tbox::pout << "    qint energy    = " << original_qint_energy
                 << std::endl;
      tbox::pout << "    free energy    = " << free_energy << std::endl;
      if (d_model_parameters.with_third_phase()) {
         tbox::pout << "    eta energy     = " << eta_energy << std::endl;
      }
   }

   extendGrainOrientation();

   if (d_symmetry_aware) {
      if (d_fundamental_interval->isActive()) {
         makeQuatFundamental(d_patch_hierarchy, d_time);
      }

      computeSymmetryRotations(d_patch_hierarchy, d_time);
   }

   double orient_energy = 10.e6;
   double qint_energy = 10.e6;
   for (int it = 0; it < 10; it++) {
      smoothQuat(d_patch_hierarchy, d_time);
      evaluateEnergy(d_patch_hierarchy, d_time, total_energy, interface_energy,
                     eta_energy, orient_energy, qint_energy, well_energy,
                     free_energy);

      tbox::pout << std::setprecision(12);
      tbox::pout << "Smooth out quaternions, orient energy  = " << orient_energy
                 << ", qint energy    = " << qint_energy
                 << ", total quat energy    = " << orient_energy + qint_energy
                 << std::endl;
   }

   math::HierarchyCellDataOpsReal<double> cellops(d_patch_hierarchy);
   cellops.copyData(d_quat_relax_id, d_quat_id, false);
   // commented out jlf 10/30/2021
   //   d_integrator_quat_only->initialize(d_patch_hierarchy);

   // double time=d_time;
   double time = 0.;

   tbox::pout << "Relax quaternions..." << std::endl;
   for (int it = 1; it < 200; it++) {

      const double dt = d_integrator_quat_only->Advance(d_patch_hierarchy);
      time += dt;

      double dqe = qint_energy;
      double doe = orient_energy;

      // diagnostics
      cellops.copyData(d_quat_id, d_quat_relax_id, false);
      evaluateEnergy(d_patch_hierarchy, d_time, total_energy, interface_energy,
                     eta_energy, orient_energy, qint_energy, well_energy,
                     free_energy);

      dqe -= qint_energy;
      dqe = fabs(dqe);
      doe -= orient_energy;
      doe = fabs(doe);

      tbox::pout << std::setprecision(12);
      tbox::pout << "Smooth out quaternions with dt=" << dt
                 << ", orient energy  = " << orient_energy
                 << ", qint energy    = " << qint_energy
                 << ", total quat energy    = " << orient_energy + qint_energy
                 << std::endl;

      const double tol = 0.1;
      if (dqe < tol * interface_energy * dt &&
          qint_energy < original_qint_energy + tol) {
         tbox::pout << "Quaternions converged after " << it << " iterations..."
                    << std::endl;
         break;
      }
   }

   cellops.copyData(d_quat_id, d_quat_relax_id, false);

   if (d_extra_energy_detail) {
      tbox::pout << std::setprecision(8);
      tbox::pout << "  Total energy      = " << total_energy << std::endl;
      tbox::pout << "  interface energy  = " << interface_energy << std::endl;
      tbox::pout << "    orient energy   = " << orient_energy << std::endl;
      tbox::pout << "    qint energy     = " << qint_energy << std::endl;
      tbox::pout << "    free energy     = " << free_energy << std::endl;
      if (d_model_parameters.with_third_phase()) {
         tbox::pout << "    eta energy     = " << eta_energy << std::endl;
      }
   }

   t_resetGrains_timer->stop();

   return true;
}

//=======================================================================

double QuatModel::Advance(void)
{
   bool update_solution = false;

   if (d_model_parameters.with_concentration() &&
       d_model_parameters.isConcentrationModelCALPHAD())
      resetRefPhaseConcentrations();

   if (d_model_parameters.with_orientation()) {
      if (d_grain_extend_interval->isActive()) {
         update_solution = resetGrains();
      }

      if (d_symmetry_aware && !update_solution) {
         if (d_fundamental_interval->hasIntervalPassed(d_cycle, d_time)) {
            makeQuatFundamental(d_patch_hierarchy, d_time);
            update_solution = true;

            if (d_symmetry_aware)
               computeSymmetryRotations(d_patch_hierarchy, d_time);
         }
      }
   }

   if (d_test_interval->hasIntervalPassed(d_cycle, d_time)) {
      // do some test, e.g.,
      // update_solution = true;
   }

   // jlf, 10/30/2021: this call is not working, possibly due to some
   // objects not being clean up properly before being reinitialized
   // however, call not necessary if mesh not adapted
   // if (update_solution) {
   //   d_integrator->updateSolution(d_patch_hierarchy, 0,
   //                                d_patch_hierarchy->getFinestLevelNumber());
   //}

   const double dt = d_integrator->Advance(d_patch_hierarchy);

   d_time += dt;

   if (d_model_parameters.with_heat_equation())
      d_heat_capacity_strategy->setCurrentValue(d_patch_hierarchy);


   return dt;
}

//-----------------------------------------------------------------------

void QuatModel::postAdvanceDiagnostics(void)
{
   if (d_scalar_diag_interval->hasIntervalPassed(d_cycle, d_time)) {

      printScalarDiagnostics();
   }
}

//-----------------------------------------------------------------------

void QuatModel::preRunDiagnostics(void)
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   math::HierarchyCellDataOpsReal<double> mathops(d_patch_hierarchy);

   double surface_e = d_model_parameters.surfaceEnergy();
   tbox::pout << "Surface energy (J/m^2): " << surface_e << std::endl;
   double width = d_model_parameters.interfacialWidth();
   tbox::pout << "Interfacial width (um): " << width << std::endl;

   if (d_model_parameters.with_concentration())
      d_composition_rhs_strategy->printDiagnostics(d_patch_hierarchy);

   const double temperature =
       d_temperature_strategy->getCurrentMinTemperature(d_patch_hierarchy,
                                                        d_time);

   if (d_fenergy_diag_filename != "") {

      // pre-run diagnostics
      if (d_model_parameters.with_concentration() &&
          d_model_parameters.isConcentrationModelCALPHAD()) {
         assert(temperature > 0.);
         if (mpi.getRank() == 0) {
            assert(d_free_energy_strategy_for_diffusion);
            d_free_energy_strategy_for_diffusion->preRunDiagnostics(
                temperature);
         }

         // compute equilibrium composition for pair L,A
         double phi_min = mathops.min(d_phase_id);
         // double phi_max = mathops.max( d_phase_id );

         double ceq[4];  // 2 phases x 2 compositions max.
         bool found_ceq = false;
         if (phi_min < 0.1) {
            found_ceq = computeCeq(temperature, Thermo4PFM::PhaseIndex::phaseL,
                                   Thermo4PFM::PhaseIndex::phaseA, &ceq[0]);
         }

         mpi.Barrier();

         if (found_ceq && !d_is_from_restart) {
            setRefPhaseConcentrationsToEquilibrium(ceq);
            // if (d_model_parameters.initPhaseConcAtEq())
            //   setPhaseConcentrationsToEquilibrium(ceq);
         }

      }  // with_concentration
   }

   mpi.Barrier();

   if (d_scalar_diag_interval->includeInitial(d_time)) {
      printScalarDiagnostics();
   }

   if (d_model_parameters.with_concentration() &&
       d_model_parameters.isConcentrationModelCALPHAD())
      preRunDiagnosticsMobilityInPhases(temperature, d_model_parameters,
                                        d_calphad_db);
}

//-----------------------------------------------------------------------

bool QuatModel::computeCeq(const double temperature,
                           const Thermo4PFM::PhaseIndex pi0,
                           const Thermo4PFM::PhaseIndex pi1, double* ceq) const
{
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   math::HierarchyCellDataOpsReal<double> mathops(d_patch_hierarchy);

   double cmin = mathops.min(d_conc_id);
   double cmax = mathops.max(d_conc_id);
   double dc = cmax - cmin;
   cmin = std::max(0.25 * cmin, cmin - 0.2 * dc);
   cmax = std::min(1. - 0.25 * (1. - cmax), cmax + 0.2 * dc);
   if (cmax - cmin <= 1.e-8) {
      cmax = std::min(1., cmin + 0.1);
      cmin = std::max(0., cmin - 0.1);
   }
   tbox::pout << "QuatModel::computeCeq(): " << std::endl
              << "Try to estimate equilibrium concentrations between cmin="
              << cmin << " and cmax=" << cmax << "..." << std::endl;
   tbox::pout << "T=" << temperature << std::endl;

   double ceq_init0 = cmin;
   double ceq_init1 = cmax;
   double lceq[4] = {ceq_init0, ceq_init1, ceq_init0, ceq_init1};

   if (d_model_parameters.knownInitCinPhase()) {
      const unsigned int offset = d_ncompositions;
      for (int ic = 0; ic < d_ncompositions; ic++) {
         lceq[ic] = d_model_parameters.getInitCphaseL(ic);
         lceq[offset + ic] = d_model_parameters.getInitCphaseA(ic);
      }
   }

   // compute equilibrium concentrations
   bool found_ceq = false;
   if (mpi.getRank() ==
       0)  // do it on PE0 only to avoid error message prints from all PEs
   {
      std::shared_ptr<ConcFreeEnergyStrategy> free_energy_strategy =
          std::dynamic_pointer_cast<ConcFreeEnergyStrategy>(
              d_free_energy_strategy);
      assert(free_energy_strategy);

      found_ceq =
          free_energy_strategy->computeCeqT(temperature, pi0, pi1, &lceq[0]);
      if (lceq[0] > 1.) found_ceq = false;
      if (lceq[0] < 0.) found_ceq = false;
      if (lceq[1] > 1.) found_ceq = false;
      if (lceq[1] < 0.) found_ceq = false;

      if (!found_ceq) {
         tbox::pout << "Try again with different initial conditions..."
                    << std::endl;
         lceq[0] = ceq_init1;
         lceq[1] = ceq_init0;
         found_ceq =
             free_energy_strategy->computeCeqT(temperature, pi0, pi1, &lceq[0]);
         if (lceq[0] > 1.) found_ceq = false;
         if (lceq[0] < 0.) found_ceq = false;
         if (lceq[1] > 1.) found_ceq = false;
         if (lceq[1] < 0.) found_ceq = false;
      }

      if (found_ceq) {
         tbox::plog << "Found equilibrium concentrations: " << lceq[0] << ", "
                    << lceq[1] << "..." << std::endl;
         if (d_ncompositions > 1)
            tbox::plog << "                                  " << lceq[2]
                       << ", " << lceq[3] << "..." << std::endl;
      } else {
         tbox::plog << "ERROR: Equilibrium concentrations not found... "
                    << std::endl;
      }

      if (!found_ceq) mpi.abort();
   }

   mpi.Bcast(&lceq[0], 4, MPI_DOUBLE, 0);

   int flag = (int)found_ceq;
   mpi.Bcast(&flag, 1, MPI_INT, 0);
   found_ceq = (bool)flag;

   ceq[0] = lceq[0];
   ceq[1] = lceq[1];
   ceq[2] = lceq[2];
   ceq[3] = lceq[3];

   return found_ceq;
}

//-----------------------------------------------------------------------

void QuatModel::postRunDiagnostics(void)
{
   if (d_scalar_diag_interval->includeFinal(d_time)) {

      printScalarDiagnostics();
   }

   d_integrator->printSolverTotals();
}

//-----------------------------------------------------------------------

void QuatModel::writeRestartFile(void) { PFModel::writeRestartFile(); }

//=======================================================================

void QuatModel::printScalarDiagnostics(void)
{
   double total_energy, interface_energy, orient_energy, qint_energy;
   double well_energy, bulk_energy, eta_energy;

   evaluateEnergy(d_patch_hierarchy, d_time, total_energy, interface_energy,
                  eta_energy, orient_energy, qint_energy, well_energy,
                  bulk_energy, d_model_parameters.grand_potential());

   if (d_model_parameters.with_heat_equation()) {
      double thermal_energy = computeThermalEnergy(d_patch_hierarchy);
      tbox::pout << "Thermal energy [pJ]= " << thermal_energy << std::endl;
   }

   if (!d_time_info_interval->eventOccurredAtTime(d_time)) {
      tbox::pout << "cycle # " << d_cycle << " : t = " << d_time << std::endl;
   }

   if (d_extra_energy_detail) {
      tbox::pout << std::setprecision(8);
      tbox::pout << "  Total energy [pJ]     = " << total_energy << std::endl;
      tbox::pout << "  interface energy [pJ] = " << interface_energy
                 << std::endl;
      tbox::pout << "    orient energy [pJ]  = " << orient_energy << std::endl;
      tbox::pout << "    qint energy [pJ]    = " << qint_energy << std::endl;
      tbox::pout << "    bulk energy [pJ]    = " << bulk_energy << std::endl;
      if (d_model_parameters.grand_potential()) {
         tbox::pout << "  GP [pJ] = " << total_energy << std::endl;
      }
      if (d_model_parameters.with_third_phase()) {
         tbox::pout << "    eta energy [pJ]    = " << eta_energy << std::endl;
      }
   } else {
      tbox::pout << "  Total energy [pJ] = " << total_energy << std::endl;
   }

   if (d_temperature_strategy) {
      double t =
          d_temperature_strategy->getCurrentMinTemperature(d_patch_hierarchy,
                                                           d_time);
      tbox::pout << "  Min. Temperature = " << t << std::endl;
      t = d_temperature_strategy->getCurrentMaxTemperature(d_patch_hierarchy,
                                                           d_time);
      tbox::pout << "  Max. Temperature = " << t << std::endl;
      t = d_temperature_strategy->getCurrentAverageTemperature(
          d_patch_hierarchy, d_time);
      tbox::pout << "  Average Temperature = " << t << std::endl;
   }

   // computes volume of physical domain
   const double* low = d_grid_geometry->getXLower();
   const double* up = d_grid_geometry->getXUpper();
   double vol = 1.;
   for (int d = 0; d < NDIM; d++)
      vol *= (up[d] - low[d]);

   double vphi = vol;
   if (d_model_parameters.with_phase()) {
      if (d_model_parameters.norderp() > 1) {
         for (int i = 0; i < d_model_parameters.norderp(); i++) {
            vphi = evaluatePhaseFraction(d_patch_hierarchy, d_phase_id, i);
            tbox::pout << "  Volume fraction of phase " << i << " = "
                       << vphi / vol << std::endl;
         }

      } else {
         vphi = evaluateVolumeSolid(d_patch_hierarchy, d_phase_id);
         tbox::pout << "  Volume fraction of solid phase = " << vphi / vol
                    << std::endl;
      }
      math::HierarchyCellDataOpsReal<double> mathops(d_patch_hierarchy);
      double m = mathops.max(d_phase_mobility_id);
      tbox::pout << "  Max. Phase mobility = " << m << std::endl;
      m = mathops.min(d_phase_mobility_id);
      tbox::pout << "  Min. Phase mobility = " << m << std::endl;
   }

   if (d_model_parameters.with_third_phase()) {
      const double vphi_eta = evaluateVolumeSolid(d_patch_hierarchy, d_eta_id);
      tbox::pout << "  Volume fraction of eta phase = " << vphi_eta / vol
                 << std::endl;
   }

   if (d_model_parameters.with_concentration()) {
      assert(d_work_id != -1);

      math::HierarchyCellDataOpsReal<double> mathops(d_patch_hierarchy);

      for (int ic = 0; ic < d_ncompositions; ic++) {
         const double c0V0 =
             evaluateIntegralConcentration(d_patch_hierarchy, ic);
         tbox::pout << "  Integral concentration " << ic << "= " << c0V0
                    << std::endl;

         copyDepthCellData(d_patch_hierarchy, d_work_id, 0, d_conc_id, ic);

         double cmax = mathops.max(d_work_id);
         tbox::pout << "  Max. concentration " << ic << "= " << cmax
                    << std::endl;

         if (d_model_parameters.norderp() < 1) {
            // average concentration
            const double c0 = c0V0 / vol;

            // now computes coring factor according to HBSM formula
            const double cphi =
                evaluateIntegralPhaseConcentration(d_patch_hierarchy, ic);

            const double cex = (cphi - c0 * vphi) / c0V0;
            tbox::pout << "  Cex (HBSM) for component " << ic << " = " << cex
                       << std::endl;
         }
      }
   }

   if (d_model_parameters.norderp() > 1 && d_subdomain_diagnostics[0] > 0.) {
      const double density = computeDensityDiagnostics();
      tbox::pout << "  Density = " << density << std::endl;
   }
}

//=======================================================================

void QuatModel::findAndNumberGrains(void)
{
   tbox::pout << "findAndNumberGrains" << std::endl;
   assert(d_grains);

   d_grains->findAndNumberGrains(d_patch_hierarchy, d_phase_id, d_weight_id,
                                 d_time);
}

//=======================================================================

void QuatModel::computeGrainDiagnostics(void)
{
   if (d_model_parameters.with_orientation()) {
      tbox::pout << "Computing grain diagnostics" << std::endl;

      findAndNumberGrains();
      d_grains->computeGrainVolumes(d_patch_hierarchy, d_weight_id);

      if (d_model_parameters.with_concentration()) {
         d_grains->computeGrainConcentrations(d_patch_hierarchy, d_time,
                                              d_conc_id, d_weight_id);
      }
   }
}

double QuatModel::computeDensityDiagnostics(void)
{
   assert(d_weight_diagnostics_id != -1);
   assert(d_work_id != -1);

   math::HierarchyCellDataOpsReal<double> mathops(d_patch_hierarchy);

   // compute volume fraction of porosity
   const int depth_void =
       d_model_parameters.norderpA() + d_model_parameters.norderpB();
   const double volume =
       mathops.sumControlVolumes(d_weight_id, d_weight_diagnostics_id);
   // tbox::pout << "volume = " << volume << std::endl;

   if (depth_void < d_model_parameters.norderp()) {
      // use order parameter associated with vacuum/liquid
      copyDepthCellData(d_patch_hierarchy, d_work_id, 0, d_phase_id,
                        depth_void);
      const double integral =
          mathops.integral(d_work_id, d_weight_diagnostics_id);
      // tbox::pout << "integral = " << integral << std::endl;
      assert(volume > 0.);

      // return relative density
      return 1. - integral / volume;
   } else {
      // use concentration field
      const double integral =
          mathops.integral(d_conc_id, d_weight_diagnostics_id);
      // tbox::pout << "integral = " << integral << std::endl;
      assert(volume > 0.);

      // return relative density
      return integral / volume;
   }
}

//=======================================================================

void QuatModel::Regrid(const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   d_regrid_refine_alg.reset(new xfer::RefineAlgorithm());

   d_regrid_refine_alg->registerRefine(d_phase_id,          // destination
                                       d_phase_id,          // source
                                       d_phase_scratch_id,  // temporary
                                       d_phase_refine_op);

   d_regrid_refine_alg->registerRefine(d_temperature_id,          // destination
                                       d_temperature_id,          // source
                                       d_temperature_scratch_id,  // temporary
                                       d_phase_refine_op);

   if (d_model_parameters.with_third_phase()) {
      d_regrid_refine_alg->registerRefine(d_eta_id,          // destination
                                          d_eta_id,          // source
                                          d_eta_scratch_id,  // temporary
                                          d_eta_refine_op);
   }

   if (d_model_parameters.with_orientation()) {
      assert(d_quat_scratch_id >= 0);
      d_regrid_refine_alg->registerRefine(d_quat_id,          // destination
                                          d_quat_id,          // source
                                          d_quat_scratch_id,  // temporary
                                          d_quat_refine_op);
      if (d_model_parameters.evolveQuat())
         d_regrid_refine_alg->registerRefine(
             d_quat_relax_id,          // destination
             d_quat_relax_id,          // source
             d_quat_relax_scratch_id,  // temporary
             d_quat_refine_op);
   }

   if (d_model_parameters.with_heat_equation()) {
      assert(d_temperature_scratch_id >= 0);
      d_regrid_refine_alg->registerRefine(
          d_temperature_id,          // destination
          d_temperature_id,          // source
          d_temperature_scratch_id,  // temporary
          d_phase_refine_op);
   }

   if (d_model_parameters.with_concentration()) {
      assert(d_conc_scratch_id >= 0);
      d_regrid_refine_alg->registerRefine(d_conc_id,          // destination
                                          d_conc_id,          // source
                                          d_conc_scratch_id,  // temporary
                                          d_conc_refine_op);
   }

#ifdef USE_CPODE
   if (d_use_warm_start) {

      std::set<int> cpodes_id_set;
      std::set<int> phase_id_set;
      std::set<int> eta_id_set;
      std::set<int> orient_id_set;
      std::set<int> conc_id_set;
      std::set<int> temp_id_set;

      d_integrator->getCPODESIdsRequiringRegrid(cpodes_id_set, phase_id_set,
                                                eta_id_set, orient_id_set,
                                                conc_id_set, temp_id_set);

      std::set<int>::iterator it;

      static std::map<int, int> id_map;

      static bool first_time = true;

      if (first_time) {
         hier::VariableDatabase* variable_db =
             hier::VariableDatabase::getDatabase();

         std::shared_ptr<hier::VariableContext> scratch =
             variable_db->getContext("SCRATCH");

         for (it = cpodes_id_set.begin(); it != cpodes_id_set.end(); it++) {

            std::ostringstream name;
            name << "warm_start_tmp_" << *it;

            std::shared_ptr<pdat::CellVariable<double> > var(
                new pdat::CellVariable<double>(tbox::Dimension(NDIM),
                                               name.str()));

            int id = variable_db->registerVariableAndContext(
                var, scratch,
                hier::IntVector(tbox::Dimension(NDIM), d_nghosts));

            id_map[*it] = id;
         }

         first_time = false;
      }

      for (it = phase_id_set.begin(); it != phase_id_set.end(); it++) {
         d_regrid_refine_alg->registerRefine(*it, *it, id_map[*it],
                                             d_phase_refine_op);
      }

      for (it = eta_id_set.begin(); it != eta_id_set.end(); it++) {
         d_regrid_refine_alg->registerRefine(*it, *it, id_map[*it],
                                             d_phase_refine_op);
      }

      for (it = orient_id_set.begin(); it != orient_id_set.end(); it++) {
         d_regrid_refine_alg->registerRefine(*it, *it, id_map[*it],
                                             d_quat_refine_op);
      }

      for (it = conc_id_set.begin(); it != conc_id_set.end(); it++) {
         d_regrid_refine_alg->registerRefine(*it, *it, id_map[*it],
                                             d_conc_refine_op);
      }

      for (it = temp_id_set.begin(); it != temp_id_set.end(); it++) {
         d_regrid_refine_alg->registerRefine(*it, *it, id_map[*it],
                                             d_phase_refine_op);
      }

   }  // if ( d_use_warm_start )
#endif

   PFModel::Regrid(hierarchy);
}

//=======================================================================
//
// Methods inherited from Serializable
//

void QuatModel::putToRestart(const std::shared_ptr<tbox::Database>& db) const
{
   PFModel::putToRestart(db);
}

//=======================================================================
//
// Method inherited from StandardTagAndInitStrategy
//

void QuatModel::initializeLevelData(
    /*! Hierarchy to initialize */
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    /*! Level to initialize */
    const int level_number, const double time, const bool can_be_refined,
    /*! Whether level is being introduced for the first time */
    const bool initial_time,
    /*! Level to copy data from */
    const std::shared_ptr<hier::PatchLevel>& old_level,
    const bool allocate_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((level_number >= 0) &&
               (level_number <= hierarchy->getFinestLevelNumber()));
   if (old_level) {
      TBOX_ASSERT(level_number == old_level->getLevelNumber());
   }
   TBOX_ASSERT((hierarchy->getPatchLevel(level_number)));
#endif

   tbox::pout << "QuatModel::initializeLevelData()" << std::endl;

   assert(d_curr_to_curr_refine_alg);

   // Note that this method is pure virtual in PFModel and MUST be
   // implemented here. However, we might also have some default behavior in
   // PFModel, so it should be called.
   PFModel::initializeLevelData(hierarchy, level_number, time, can_be_refined,
                                initial_time, old_level, allocate_data);

   std::shared_ptr<hier::PatchLevel> level(
       hierarchy->getPatchLevel(level_number));
   assert(level);

   AllocateLocalPatchData(level, time, allocate_data);

   d_integrator->initializeLevelData(hierarchy, level_number, time,
                                     can_be_refined, initial_time, old_level,
                                     allocate_data);
   if (d_model_parameters.evolveQuat())
      d_integrator_quat_only->initializeLevelData(hierarchy, level_number, time,
                                                  can_be_refined, initial_time,
                                                  old_level, allocate_data);

   d_grains->initializeLevelData(hierarchy, level_number, time, can_be_refined,
                                 initial_time, old_level, allocate_data);


   if (initial_time) {
      tbox::pout << "Fill data for initial time..." << std::endl;
      if (level_number == d_level_of_init_data) {

         d_initial_level =
             d_patch_hierarchy->getPatchLevel(d_level_of_init_data);
         assert(d_initial_level);

         std::shared_ptr<xfer::RefineSchedule> schedule(
             d_curr_to_curr_refine_alg->createSchedule(
                 level, d_initial_level, d_all_refine_patch_strategy));

         schedule->fillData(time, false);

      } else if (level_number > d_level_of_init_data) {

         std::shared_ptr<xfer::RefineSchedule> schedule(
             d_curr_to_curr_refine_alg->createSchedule(
                 level, old_level, level_number - 1, hierarchy,
                 d_all_refine_patch_strategy));

         schedule->fillData(time, false);
      } else {
         xfer::CoarsenAlgorithm coarsen_alg(tbox::Dimension(NDIM));

         coarsen_alg.registerCoarsen(d_phase_id, d_phase_id,
                                     d_grid_geometry->lookupCoarsenOperator(
                                         d_phase_var, "CONSERVATIVE_COARSEN"));

         if (d_model_parameters.with_third_phase()) {
            coarsen_alg.registerCoarsen(
                d_eta_id, d_eta_id,
                d_grid_geometry->lookupCoarsenOperator(d_eta_var,
                                                       "CONSERVATIVE_"
                                                       "COARSEN"));
         }

         if (d_model_parameters.with_heat_equation()) {
            coarsen_alg.registerCoarsen(
                d_temperature_id, d_temperature_id,
                d_grid_geometry->lookupCoarsenOperator(d_temperature_var,
                                                       "CONSERVATIVE_"
                                                       "COARSEN"));
         }

         if (d_model_parameters.with_concentration()) {
            coarsen_alg.registerCoarsen(
                d_conc_id, d_conc_id,
                d_grid_geometry->lookupCoarsenOperator(d_conc_var,
                                                       "CONSERVATIVE_"
                                                       "COARSEN"));
         }

         if (d_model_parameters.with_orientation()) {
            assert(d_quat_coarsen_op);

            coarsen_alg.registerCoarsen(d_quat_id, d_quat_id,
                                        d_quat_coarsen_op);
         }
         std::shared_ptr<xfer::CoarsenSchedule> schedule =
             coarsen_alg.createSchedule(level, d_initial_level, nullptr);

         schedule->coarsenData();
      }

   } else {  // if ( initial_time )

      if ((level_number > 0) || old_level) {
         // move solution from old_level to level
         d_regrid_refine_alg
             ->createSchedule(level, old_level, level_number - 1, hierarchy,
                              d_all_refine_patch_strategy)
             ->fillData(time);
      }
   }

   if (!d_model_parameters.with_concentration()) {
      math::HierarchyCellDataOpsReal<double> mathops(hierarchy);

      if (d_model_parameters.free_energy_type()[0] == 's') {
         mathops.setToScalar(d_f_l_id, d_model_parameters.free_energy_liquid());
      }
      mathops.setToScalar(d_f_a_id, d_model_parameters.free_energy_solid_A());
      if (d_model_parameters.withPhaseB()) {
         mathops.setToScalar(d_f_b_id,
                             d_model_parameters.free_energy_solid_B());
      }
   }
}

//=======================================================================

void QuatModel::AllocateQuatLocalPatchData(
    const std::shared_ptr<hier::PatchLevel> level, const double time,
    const bool zero_data)
{
   assert(d_quat_id >= 0);

   if (!level->checkAllocated(d_quat_id))
      level->allocatePatchData(d_quat_id, time);

   AllocateAndZeroData<pdat::CellData<double> >(d_quat_scratch_id, level, time,
                                                zero_data);

   if (d_model_parameters.evolveQuat()) {
      if (!level->checkAllocated(d_quat_relax_id))
         level->allocatePatchData(d_quat_relax_id, time);
      AllocateAndZeroData<pdat::CellData<double> >(d_quat_relax_scratch_id,
                                                   level, time, zero_data);
      AllocateAndZeroData<pdat::CellData<double> >(d_quat_mobility_id, level,
                                                   time, zero_data);

      AllocateAndZeroData<pdat::CellData<double> >(d_quat_grad_cell_id, level,
                                                   time, zero_data);

      AllocateAndZeroData<pdat::SideData<double> >(d_quat_grad_side_id, level,
                                                   time, zero_data);

      AllocateAndZeroData<pdat::CellData<double> >(d_quat_grad_modulus_id,
                                                   level, time, zero_data);

      AllocateAndZeroData<pdat::SideData<double> >(d_quat_diffs_id, level, time,
                                                   zero_data);

      if (d_symmetry_aware) {
         assert(d_quat_symm_rotation_id >= 0);
         AllocateAndZeroData<pdat::SideData<int> >(d_quat_symm_rotation_id,
                                                   level, time, zero_data);
      }

      if (d_model_parameters.with_extra_visit_output()) {
         AllocateAndZeroData<pdat::CellData<double> >(d_quat_diffs_cell_id,
                                                      level, time, zero_data);
         AllocateAndZeroData<pdat::CellData<double> >(
             d_quat_nonsymm_diffs_cell_id, level, time, zero_data);

         AllocateAndZeroData<pdat::CellData<double> >(d_quat_norm_error_id,
                                                      level, time, zero_data);

         if (d_symmetry_aware) {
            assert(d_quat_symm_rotation_cell_id >= 0);
            AllocateAndZeroData<pdat::CellData<int> >(
                d_quat_symm_rotation_cell_id, level, time, zero_data);
         }
      }
   }

   if (d_model_parameters.with_visit_energy_output()) {
      AllocateAndZeroData<pdat::CellData<double> >(d_energy_diag_id, level,
                                                   time, zero_data);
   }
}
//=======================================================================
// Some data may already have been initialized with restart data, so we
// need to check if already allocated. Otherwise that data may be reset
// to 0.
void QuatModel::AllocateLocalPatchData(
    const std::shared_ptr<hier::PatchLevel> level, const double time,
    const bool zero_data)
{
   assert(d_temperature_scratch_id >= 0);

   if (d_phase_id >= 0) {
      if (!level->checkAllocated(d_phase_id)) {
         level->allocatePatchData(d_phase_id, time);
      }

      AllocateAndZeroData<pdat::CellData<double> >(d_phase_scratch_id, level,
                                                   time, zero_data);
   }

   if (d_model_parameters.with_third_phase()) {
      if (!level->checkAllocated(d_eta_id)) {
         level->allocatePatchData(d_eta_id, time);
      }

      AllocateAndZeroData<pdat::CellData<double> >(d_eta_scratch_id, level,
                                                   time, zero_data);
   }

   if (!level->checkAllocated(d_temperature_id)) {
      level->allocatePatchData(d_temperature_id, time);
   }

   AllocateAndZeroData<pdat::CellData<double> >(d_temperature_scratch_id, level,
                                                time, zero_data);
   if (d_model_parameters.with_steady_temperature()) {
      assert(d_temperature_rhs_steady_id >= 0);
      AllocateAndZeroData<pdat::CellData<double> >(d_temperature_rhs_steady_id,
                                                   level, time, zero_data);
   }
   if (d_model_parameters.with_bias_well()) {
      assert(d_equilibrium_temperature_id >= 0);
      AllocateAndZeroData<pdat::CellData<double> >(d_equilibrium_temperature_id,
                                                   level, time, zero_data);
   }

   if (d_model_parameters.with_heat_equation()) {
      assert(d_cp_id >= 0);

      AllocateAndZeroData<pdat::CellData<double> >(d_cp_id, level, time,
                                                   zero_data);
   }

   if (d_model_parameters.with_concentration()) {
      if (!level->checkAllocated(d_conc_id)) {
         AllocateAndZeroData<pdat::CellData<double> >(d_conc_id, level, time,
                                                      zero_data);
      }
      tbox::plog << "AllocateAndZeroData for d_conc_scratch_id" << std::endl;
      AllocateAndZeroData<pdat::CellData<double> >(d_conc_scratch_id, level,
                                                   time, zero_data);
      for (int ic = 0; ic < d_ncompositions; ic++) {
         AllocateAndZeroData<pdat::SideData<double> >(
             d_conc_pfm_diffusion_id[ic], level, time, zero_data);
      }
      if (d_model_parameters.concRHSstrategyIsKKS() ||
          d_model_parameters.concRHSstrategyIsBeckermann()) {
         AllocateAndZeroData<pdat::SideData<double> >(
             d_conc_phase_coupling_diffusion_id, level, time, zero_data);
      }

      if (d_model_parameters.concentrationModelNeedsPhaseConcentrations()) {
         if (d_conc_l_ref_id >= 0)
            if (!level->checkAllocated(d_conc_l_ref_id)) {
               AllocateAndZeroData<pdat::CellData<double> >(d_conc_l_ref_id,
                                                            level, time,
                                                            zero_data);
               AllocateAndZeroData<pdat::CellData<double> >(d_conc_a_ref_id,
                                                            level, time,
                                                            zero_data);
               if (d_model_parameters.withPhaseB())
                  AllocateAndZeroData<pdat::CellData<double> >(d_conc_b_ref_id,
                                                               level, time,
                                                               zero_data);
            }
         if (!level->checkAllocated(d_conc_l_id)) {
            AllocateAndZeroData<pdat::CellData<double> >(d_conc_l_id, level,
                                                         time, zero_data);
         }
         if (!level->checkAllocated(d_conc_a_id)) {
            AllocateAndZeroData<pdat::CellData<double> >(d_conc_a_id, level,
                                                         time, zero_data);
         }
         if (d_model_parameters.withPhaseB() &&
             !level->checkAllocated(d_conc_b_id))
            AllocateAndZeroData<pdat::CellData<double> >(d_conc_b_id, level,
                                                         time, zero_data);

         if (d_model_parameters.concRHSstrategyIsEBS()) {
            AllocateAndZeroData<pdat::SideData<double> >(
                d_conc_pfm_diffusion_l_id, level, time, zero_data);
            AllocateAndZeroData<pdat::SideData<double> >(
                d_conc_pfm_diffusion_a_id, level, time, zero_data);
            if (d_model_parameters.withPhaseB()) {
               AllocateAndZeroData<pdat::SideData<double> >(
                   d_conc_pfm_diffusion_b_id, level, time, zero_data);
            }
            AllocateAndZeroData<pdat::SideData<double> >(
                d_conc_diffusion_coeff_l_id, level, time, zero_data);
            AllocateAndZeroData<pdat::SideData<double> >(
                d_conc_diffusion_coeff_a_id, level, time, zero_data);
            if (d_model_parameters.withPhaseB()) {
               AllocateAndZeroData<pdat::SideData<double> >(
                   d_conc_diffusion_coeff_b_id, level, time, zero_data);
            }
         }
         if (d_model_parameters.with_extra_visit_output()) {
            if (d_model_parameters.concRHSstrategyIsEBS()) {
               AllocateAndZeroData<pdat::CellData<double> >(
                   d_visit_conc_pfm_diffusion_l_id, level, time, zero_data);
               AllocateAndZeroData<pdat::CellData<double> >(
                   d_visit_conc_pfm_diffusion_a_id, level, time, zero_data);
               if (d_model_parameters.withPhaseB())
                  AllocateAndZeroData<pdat::CellData<double> >(
                      d_visit_conc_pfm_diffusion_b_id, level, time, zero_data);
            } else {
               assert(d_visit_conc_pfm_diffusion_id >= 0);
               AllocateAndZeroData<pdat::CellData<double> >(
                   d_visit_conc_pfm_diffusion_id, level, time, zero_data);
            }
         }
      }  // concentrationModelNeedsPhaseConcentrations()

      if (d_model_parameters.with_extra_visit_output() &&
          !d_model_parameters.concRHSstrategyIsEBS()) {
         assert(d_visit_conc_pfm_diffusion_id >= 0);
         AllocateAndZeroData<pdat::CellData<double> >(
             d_visit_conc_pfm_diffusion_id, level, time, zero_data);
      }

      if (d_model_parameters.with_third_phase()) {
         AllocateAndZeroData<pdat::SideData<double> >(
             d_conc_eta_coupling_diffusion_id, level, time, zero_data);
      }

      if (d_model_parameters.with_gradT()) {
         AllocateAndZeroData<pdat::SideData<double> >(d_conc_Mq_id, level, time,
                                                      zero_data);
      }
      if (d_model_parameters.with_partition_coeff()) {
         AllocateAndZeroData<pdat::CellData<double> >(d_partition_coeff_id,
                                                      level, time, zero_data);
         if (d_model_parameters.needGhosts4PartitionCoeff())
            AllocateAndZeroData<pdat::CellData<double> >(
                d_partition_coeff_scratch_id, level, time, zero_data);
      }

   }  // d_model_parameters.with_concentration()

   if (d_model_parameters.with_velocity()) {
      AllocateAndZeroData<pdat::CellData<double> >(d_velocity_id, level, time,
                                                   zero_data);
   }

   if (d_model_parameters.with_phase()) {
      AllocateAndZeroData<pdat::SideData<double> >(d_phase_diffs_id, level,
                                                   time, zero_data);
      if (d_model_parameters.with_extra_visit_output()) {
         AllocateAndZeroData<pdat::CellData<double> >(d_phase_diffs_cell_id,
                                                      level, time, zero_data);
      }

      AllocateAndZeroData<pdat::CellData<double> >(d_phase_grad_cell_id, level,
                                                   time, zero_data);

      AllocateAndZeroData<pdat::CellData<double> >(d_phase_mobility_id, level,
                                                   time, zero_data);
   }

   if (d_model_parameters.with_third_phase()) {
      AllocateAndZeroData<pdat::SideData<double> >(d_eta_diffs_id, level, time,
                                                   zero_data);

      AllocateAndZeroData<pdat::CellData<double> >(d_eta_grad_cell_id, level,
                                                   time, zero_data);

      AllocateAndZeroData<pdat::CellData<double> >(d_eta_mobility_id, level,
                                                   time, zero_data);
   }

   // AllocateAndZeroData< pdat::SideData<double> >(
   //    d_phase_grad_side_id, level, time, zero_data );

   if (d_model_parameters.with_orientation()) {
      AllocateQuatLocalPatchData(level, time, zero_data);
   }

   AllocateAndZeroData<pdat::CellData<double> >(d_weight_id, level, time,
                                                zero_data);

   if (d_subdomain_diagnostics[0] > 0.)
      AllocateAndZeroData<pdat::CellData<double> >(d_weight_diagnostics_id,
                                                   level, time, zero_data);

   AllocateAndZeroData<pdat::CellData<double> >(d_work_id, level, time,
                                                zero_data);

   AllocateAndZeroData<pdat::CellData<double> >(d_f_l_id, level, time,
                                                zero_data);

   AllocateAndZeroData<pdat::CellData<double> >(d_f_a_id, level, time,
                                                zero_data);

   if (d_model_parameters.withPhaseB()) {
      AllocateAndZeroData<pdat::CellData<double> >(d_f_b_id, level, time,
                                                   zero_data);
   }
}
//-----------------------------------------------------------------------

void QuatModel::DeallocateIntermediateLocalPatchData(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   for (int ln0 = 0; ln0 <= d_patch_hierarchy->getFinestLevelNumber(); ln0++) {
      std::shared_ptr<hier::PatchLevel> level(hierarchy->getPatchLevel(ln0));

      DeallocateIntermediateLocalPatchData(level);
   }
}

void QuatModel::DeallocateIntermediateLocalPatchData(
    const std::shared_ptr<hier::PatchLevel> level)
{
   if (d_conc_l_ref_id >= 0) {
      if (level->checkAllocated(d_conc_l_ref_id)) {
         level->deallocatePatchData(d_conc_l_ref_id);
         level->deallocatePatchData(d_conc_a_ref_id);
         if (d_model_parameters.withPhaseB())
            level->deallocatePatchData(d_conc_b_ref_id);
      }
   }

   if (d_model_parameters.evolveQuat()) {
      level->deallocatePatchData(d_quat_grad_cell_id);
      level->deallocatePatchData(d_quat_grad_side_id);
      level->deallocatePatchData(d_quat_grad_modulus_id);
      level->deallocatePatchData(d_quat_diffs_id);
      if (d_model_parameters.evolveQuat()) {
         level->deallocatePatchData(d_quat_relax_id);
         level->deallocatePatchData(d_quat_relax_scratch_id);
      }
      if (d_symmetry_aware) {
         level->deallocatePatchData(d_quat_symm_rotation_id);
         if (d_model_parameters.with_extra_visit_output()) {
            level->deallocatePatchData(d_quat_symm_rotation_cell_id);
         }
      }
   }
}

//-----------------------------------------------------------------------

template <typename T>
void QuatModel::AllocateAndZeroData(
    const int data_id, const std::shared_ptr<hier::PatchLevel> level,
    const double time, const bool zero_data)
{
   assert(data_id >= 0);

   if (!level->checkAllocated(data_id)) level->allocatePatchData(data_id, time);

   if (zero_data) {
      for (hier::PatchLevel::iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<T> data(SAMRAI_SHARED_PTR_CAST<T, hier::PatchData>(
             p->getPatchData(data_id)));
         assert(data);
         data->fillAll(0);
      }
   }
}

//=======================================================================
//
// Method inherited from StandardTagAndInitStrategy

void QuatModel::resetHierarchyConfiguration(
    const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const int coarsest_level, const int finest_level)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(hierarchy);
   TBOX_ASSERT((coarsest_level >= 0) && (coarsest_level <= finest_level) &&
               (finest_level <= hierarchy->getFinestLevelNumber()));
   for (int ln0 = 0; ln0 <= finest_level; ln0++) {
      TBOX_ASSERT((hierarchy->getPatchLevel(ln0)));
   }
#endif
   tbox::pout << "QuatModel::resetHierarchyConfiguration()" << std::endl;

   int nlev = hierarchy->getNumberOfLevels();

   // Note that this method is pure virtual in PFModel and MUST be
   // implemented here. However, we might have some default behavior in
   // PFModel, so it should be called.
   PFModel::resetHierarchyConfiguration(hierarchy, coarsest_level,
                                        finest_level);

   /*
    * Check for overlapping boxes and abort if we find any.  For the time
    * being, we need to avoid this situation since the post-1.13 versions
    * of the Hypre SMG solver do not tolerate overlapping boxes.
    */
   for (int ln = coarsest_level; ln <= finest_level; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
      hier::BoxContainer bl(level->getBoxes());
      if (bl.boxesIntersect()) {
         for (hier::BoxContainer::iterator bli = bl.begin(); bli != bl.end();
              ++bli) {
            tbox::pout << *bli << std::endl;
         }
         std::stringstream message;
         message << "Boxes intersect on level " << ln;
         tbox::Utilities::abort(message.str(), __FILE__, __LINE__);
      }
   }


   d_integrator->resetHierarchyConfiguration(hierarchy, coarsest_level,
                                             finest_level);
   if (d_model_parameters.evolveQuat())
      d_integrator_quat_only->resetHierarchyConfiguration(hierarchy,
                                                          coarsest_level,
                                                          finest_level);

   d_curr_to_scr_refine_sched.resize(nlev);

   for (int ln = coarsest_level; ln <= finest_level; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      std::shared_ptr<xfer::RefineSchedule> schedule(
          d_curr_to_scr_refine_alg->createSchedule(
              level, ln - 1, hierarchy, d_all_refine_patch_strategy));
      d_curr_to_scr_refine_sched[ln] = schedule;
   }

   if (d_grain_diag_interval->isActive() ||
       d_grain_extend_interval->isActive()) {
      d_grains->resetHierarchyConfiguration(hierarchy, coarsest_level,
                                            finest_level);
   }

   computeVectorWeights(d_patch_hierarchy);
}

//=======================================================================
//
// Method inherited from StandardTagAndInitStrategy

void QuatModel::applyGradientDetector(
    std::shared_ptr<hier::PatchHierarchy>& hierarchy, int level_number,
    double time, int tag_index, const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
   std::shared_ptr<hier::PatchLevel> level =
       hierarchy->getPatchLevel(level_number);

   copyCurrentToScratch(hierarchy, level_number, time,
                        d_all_refine_patch_strategy);

   HierarchyStencilOps stencil_ops;

   if (d_tag_phase) {
      stencil_ops.computeDiffs(level, d_phase_scratch_id, d_phase_diffs_id);
      stencil_ops.computeGradCell(level, d_phase_diffs_id,
                                  d_phase_grad_cell_id);
   }

   if (d_tag_eta) {
      stencil_ops.computeDiffs(level, d_eta_scratch_id, d_eta_diffs_id);
      stencil_ops.computeGradCell(level, d_eta_diffs_id, d_eta_grad_cell_id);
   }

   if (d_tag_quat) {
      computeQuatDiffs(level, d_quat_scratch_id, d_quat_diffs_id, time);

      computeQuatGradCell(level, d_quat_diffs_id, d_quat_grad_cell_id, time);
   }

   for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
        ++ip) {
      std::shared_ptr<hier::Patch> patch = *ip;

      tagGradientDetectorCells(*patch, time, initial_time, tag_index,
                               uses_richardson_extrapolation_too);
   }
}

void QuatModel::tagGradientDetectorCells(
    hier::Patch& patch, const double regrid_time, const bool initial_error,
    const int tag_index, const bool uses_richardson_extrapolation_too)
{
   (void)regrid_time;
   (void)initial_error;
   (void)uses_richardson_extrapolation_too;

   std::shared_ptr<pdat::CellData<int> > tags(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
           patch.getPatchData(tag_index)));
   assert(tags);
   tags->fillAll(0);

   const hier::Index& patch_lower = patch.getBox().lower();
   const hier::Index& patch_upper = patch.getBox().upper();

   const hier::Box& tag_ghost_box = tags->getGhostBox();
   const hier::Index& tag_gbox_lower = tag_ghost_box.lower();
   const hier::Index& tag_gbox_upper = tag_ghost_box.upper();

   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));
   TBOX_ASSERT(patch_geom);

   const double* dx = patch_geom->getDx();

   if (d_tag_phase) {
      std::shared_ptr<pdat::CellData<double> > phase_grad(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(d_phase_grad_cell_id)));
      assert(phase_grad);

      const hier::Box& gbox = phase_grad->getGhostBox();
      const hier::Index& gbox_lower = gbox.lower();
      const hier::Index& gbox_upper = gbox.upper();

      TAGFROMGRADS(patch_lower(0), patch_upper(0), patch_lower(1),
                   patch_upper(1),
#if (NDIM == 3)
                   patch_lower(2), patch_upper(2),
#endif
                   phase_grad->getPointer(0), phase_grad->getPointer(1),
#if (NDIM == 3)
                   phase_grad->getPointer(2),
#endif
                   gbox_lower(0), gbox_upper(0), gbox_lower(1), gbox_upper(1),
#if (NDIM == 3)
                   gbox_lower(2), gbox_upper(2),
#endif
                   dx, tags->getPointer(0), tag_gbox_lower(0),
                   tag_gbox_upper(0), tag_gbox_lower(1), tag_gbox_upper(1),
#if (NDIM == 3)
                   tag_gbox_lower(2), tag_gbox_upper(2),
#endif
                   true, d_phase_threshold_untagged);
   }

   if (d_tag_eta) {
      std::shared_ptr<pdat::CellData<double> > eta_grad(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(d_eta_grad_cell_id)));
      assert(eta_grad);

      const hier::Box& gbox = eta_grad->getGhostBox();
      const hier::Index& gbox_lower = gbox.lower();
      const hier::Index& gbox_upper = gbox.upper();

      TAGFROMGRADS(patch_lower(0), patch_upper(0), patch_lower(1),
                   patch_upper(1),
#if (NDIM == 3)
                   patch_lower(2), patch_upper(2),
#endif
                   eta_grad->getPointer(0), eta_grad->getPointer(1),
#if (NDIM == 3)
                   eta_grad->getPointer(2),
#endif
                   gbox_lower(0), gbox_upper(0), gbox_lower(1), gbox_upper(1),
#if (NDIM == 3)
                   gbox_lower(2), gbox_upper(2),
#endif
                   dx, tags->getPointer(0), tag_gbox_lower(0),
                   tag_gbox_upper(0), tag_gbox_lower(1), tag_gbox_upper(1),
#if (NDIM == 3)
                   tag_gbox_lower(2), tag_gbox_upper(2),
#endif
                   true, d_eta_threshold_untagged);
   }

   if (d_tag_quat) {

      std::shared_ptr<pdat::CellData<double> > quat_grad(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(d_quat_grad_cell_id)));
      assert(quat_grad);

      const hier::Box& gbox = quat_grad->getGhostBox();
      const hier::Index& gbox_lower = gbox.lower();
      const hier::Index& gbox_upper = gbox.upper();

      TAGFROMQUATGRADS(patch_lower(0), patch_upper(0), patch_lower(1),
                       patch_upper(1),
#if (NDIM == 3)
                       patch_lower(2), patch_upper(2),
#endif
                       quat_grad->getPointer(0 * d_qlen),
                       quat_grad->getPointer(1 * d_qlen),
#if (NDIM == 3)
                       quat_grad->getPointer(2 * d_qlen),
#endif
                       gbox_lower(0), gbox_upper(0), gbox_lower(1),
                       gbox_upper(1),
#if (NDIM == 3)
                       gbox_lower(2), gbox_upper(2),
#endif
                       d_qlen, dx, tags->getPointer(0), tag_gbox_lower(0),
                       tag_gbox_upper(0), tag_gbox_lower(1), tag_gbox_upper(1),
#if (NDIM == 3)
                       tag_gbox_lower(2), tag_gbox_upper(2),
#endif
                       true, d_quat_threshold_untagged);
   }
}
//=======================================================================

// Read initialization database

void QuatModel::readInitialDatabase(std::shared_ptr<tbox::Database> input_db)
{
   PFModel::readInitialDatabase(input_db);

   hier::BoxContainer boxes = d_grid_geometry->getPhysicalDomain();
   assert(boxes.size() == 1);
}

//=======================================================================

void QuatModel::WriteInitialConditionsFile(std::string filename, int level)
{
   FieldsWriter writer(d_model_parameters, filename, level, d_grid_geometry,
                       d_phase_id, d_phase_scratch_id, d_temperature_id,
                       d_temperature_scratch_id, d_quat_id, d_quat_scratch_id,
                       d_conc_id, d_conc_scratch_id, d_eta_id, d_eta_scratch_id,
                       d_ncompositions, d_qlen, d_all_refine_patch_strategy);

   writer.writeInitialConditionsFile(d_patch_hierarchy, d_time);
}

//=======================================================================

// phase_id is CellData
// phase_diffs_id is SideData

void QuatModel::computePhaseDiffs(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
    int& phase_diffs_id, const double time, const CACHE_TYPE cache)
{
   t_phase_diffs_timer->start();

   static double old_time = tbox::IEEE::getSignalingNaN();

   if (time == old_time && cache == CACHE) return;
   old_time = time;

   if (phase_id < 0) phase_id = d_phase_scratch_id;
   if (phase_diffs_id < 0) phase_diffs_id = d_phase_diffs_id;

   int maxln = hierarchy->getFinestLevelNumber();

   HierarchyStencilOps stencil_ops;

   stencil_ops.computeDiffs(hierarchy, phase_id, phase_diffs_id);

   if (d_model_parameters.with_extra_visit_output()) {

      for (int ln = 0; ln <= maxln; ln++) {
         std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

         for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
              ++p) {
            std::shared_ptr<hier::Patch> patch = *p;
            const hier::Box& box = patch->getBox();

            std::shared_ptr<pdat::SideData<double> > diff_data(
                SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                    patch->getPatchData(phase_diffs_id)));
            assert(diff_data);

            std::shared_ptr<pdat::CellData<double> > cell_diffs_data(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_phase_diffs_cell_id)));
            assert(cell_diffs_data);

            assert(diff_data->getDepth() == cell_diffs_data->getDepth());

            pdat::CellIterator iend(pdat::CellGeometry::end(box));
            for (pdat::CellIterator i(pdat::CellGeometry::begin(box));
                 i != iend; ++i) {
               const pdat::CellIndex ccell = *i;
               const pdat::SideIndex xside(ccell, pdat::SideIndex::X,
                                           pdat::SideIndex::Lower);
               const pdat::SideIndex yside(ccell, pdat::SideIndex::Y,
                                           pdat::SideIndex::Lower);
#if (NDIM == 3)
               const pdat::SideIndex zside(ccell, pdat::SideIndex::Z,
                                           pdat::SideIndex::Lower);
#endif
               (*cell_diffs_data)(ccell, 0) = (*diff_data)(xside);
               (*cell_diffs_data)(ccell, 1) = (*diff_data)(yside);
#if (NDIM == 3)
               (*cell_diffs_data)(ccell, 2) = (*diff_data)(zside);
#endif
            }
         }
      }
   }
   t_phase_diffs_timer->stop();
}

//=======================================================================

// eta_id is CellData
// eta_diffs_id is SideData

void QuatModel::computeEtaDiffs(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& eta_id,
    int& eta_diffs_id, const double time, const CACHE_TYPE cache)
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if (time == old_time && cache == CACHE) return;
   old_time = time;

   if (eta_id < 0) eta_id = d_eta_scratch_id;
   if (eta_diffs_id < 0) eta_diffs_id = d_eta_diffs_id;

   HierarchyStencilOps stencil_ops;

   stencil_ops.computeDiffs(hierarchy, eta_id, eta_diffs_id);
}

//=======================================================================

// var_id is CellData
// diff_id is SideData

void QuatModel::computeVarDiffs(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& var_id,
    int& diffs_id, const double time, const CACHE_TYPE cache)
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if (time == old_time && cache == CACHE) return;
   old_time = time;

   if (var_id < 0) var_id = d_phase_scratch_id;
   if (diffs_id < 0) diffs_id = d_phase_diffs_id;

   HierarchyStencilOps stencil_ops;

   stencil_ops.computeDiffs(hierarchy, var_id, diffs_id);
}

void QuatModel::smoothQuat(const std::shared_ptr<hier::PatchLevel> level)
{
   assert(d_quat_diffs_id >= 0);

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      std::shared_ptr<hier::Patch> patch = *p;
      const hier::Box& box = patch->getBox();
      const hier::Index& ifirst = box.lower();
      const hier::Index& ilast = box.upper();

      std::shared_ptr<pdat::CellData<double> > quat(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_quat_id)));
      assert(quat);
      assert(quat->getGhostCellWidth() ==
             hier::IntVector(tbox::Dimension(NDIM), 0));

      std::shared_ptr<pdat::CellData<double> > quat_scratch(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_quat_scratch_id)));
      assert(quat_scratch);

      std::shared_ptr<pdat::SideData<double> > diff_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(d_quat_diffs_id)));
      assert(diff_data);

      std::shared_ptr<pdat::CellData<double> > phase(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_phase_id)));
      assert(phase);

      SMOOTHQUAT(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                 ifirst(2), ilast(2),
#endif
                 d_qlen, quat_scratch->getPointer(),
                 quat_scratch->getGhostCellWidth()[0], quat->getPointer(), 0,
                 diff_data->getPointer(0, 0), diff_data->getPointer(1, 0),
#if (NDIM == 3)
                 diff_data->getPointer(2, 0),
#endif
                 diff_data->getGhostCellWidth()[0], phase->getPointer(),
                 phase->getGhostCellWidth()[0], d_phase_threshold);
   }
}

void QuatModel::smoothQuat(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time)
{
   // tbox::pout<<"smoothQuat..."<<endl;

   // Fill ghosts of original quat data
   copyCurrentToScratch(hierarchy, time, d_all_refine_patch_strategy);

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> patch_level =
          hierarchy->getPatchLevel(ln);

      computeQuatDiffs(patch_level, d_quat_scratch_id, d_quat_diffs_id, time);

      smoothQuat(patch_level);
   }
}

//=======================================================================

// Computes phase gradients at cell centers.

// phase_diffs_id is SideData
// phase_grad_id is CellData with no ghosts with depth NDIM

void QuatModel::computePhaseGradCell(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_diffs_id,
    int& phase_grad_cell_id, const double time, const CACHE_TYPE cache)
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if (time == old_time && cache == CACHE) return;
   old_time = time;

   if (phase_diffs_id < 0) phase_diffs_id = d_phase_diffs_id;
   if (phase_grad_cell_id < 0) phase_grad_cell_id = d_phase_grad_cell_id;

   HierarchyStencilOps stencil_ops;
   stencil_ops.computeGradCell(hierarchy, phase_diffs_id, phase_grad_cell_id);
}

void QuatModel::computePhaseGradCell(
    const std::shared_ptr<hier::PatchLevel> patch_level, int& phase_diffs_id,
    int& phase_grad_cell_id, const double time)
{
   if (phase_diffs_id < 0) phase_diffs_id = d_phase_diffs_id;
   if (phase_grad_cell_id < 0) phase_grad_cell_id = d_phase_grad_cell_id;

   HierarchyStencilOps stencil_ops;
   stencil_ops.computeGradCell(patch_level, phase_diffs_id, phase_grad_cell_id);
}

//=======================================================================

// Computes gradients at cell centers.

// diff_id is SideData
// grad_id is CellData with no ghosts with depth NDIM

void QuatModel::computeVarGradCell(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& diffs_id,
    int& grad_cell_id, const double time, const CACHE_TYPE cache)
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if (time == old_time && cache == CACHE) return;
   old_time = time;

   HierarchyStencilOps stencil_ops;
   stencil_ops.computeGradCell(hierarchy, diffs_id, grad_cell_id);
}

//-----------------------------------------------------------------------

// Computes gradients at sides.

// diff_id is SideData with at least 1 ghost layer
// grad_id is SideData with no ghosts, depth of NDIM

void QuatModel::computeVarGradSide(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& diffs_id,
    int& grad_side_id, const double time, const CACHE_TYPE cache)
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if (time == old_time && cache == CACHE) return;
   old_time = time;

   if (diffs_id < 0) diffs_id = d_phase_diffs_id;
   if (grad_side_id < 0) grad_side_id = d_phase_grad_side_id;

   HierarchyStencilOps stencil_ops;
   stencil_ops.computeGradSide(hierarchy, diffs_id, grad_side_id);
}

//=======================================================================

// Computes differences at cell sides. This will always compute the
// "normal" diffs, and will also compute the symmetry-aware diffs
// if symmetry is on.

// quat_id is CellData with ghosts
// diff_id is SideData with ghosts

void QuatModel::computeQuatDiffs(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& quat_id,
    int& quat_diffs_id, const double time, const CACHE_TYPE cache)
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if (time == old_time && cache == CACHE) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> patch_level =
          hierarchy->getPatchLevel(ln);

      computeQuatDiffs(patch_level, quat_id, quat_diffs_id, time);
   }

   assert(quat_diffs_id >= 0);
}

void QuatModel::computeQuatDiffs(const std::shared_ptr<hier::PatchLevel> level,
                                 int& quat_id, int& quat_diffs_id,
                                 const double time)
{
   (void)time;

   if (quat_id < 0) quat_id = d_quat_scratch_id;
   if (quat_diffs_id < 0) quat_diffs_id = d_quat_diffs_id;

   assert(quat_id >= 0);
   assert(quat_diffs_id >= 0);

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      std::shared_ptr<hier::Patch> patch = *p;

      std::shared_ptr<pdat::CellData<double> > quat_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(quat_id)));
      std::shared_ptr<pdat::SideData<double> > diff_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(quat_diffs_id)));
      if (d_symmetry_aware) {
         std::shared_ptr<pdat::SideData<int> > rotation_index(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<int>, hier::PatchData>(
                 patch->getPatchData(d_quat_symm_rotation_id)));
         computeQDiffs(quat_data, diff_data, true, rotation_index);
      } else {
         computeQDiffs(quat_data, diff_data, false, nullptr);
      }

      if (d_model_parameters.with_extra_visit_output()) {
         const int symm_depth_offset = 0;
         const int nonsymm_depth_offset = d_symmetry_aware ? d_qlen : 0;

         std::shared_ptr<pdat::CellData<double> > cell_diffs_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_quat_diffs_cell_id)));
         assert(cell_diffs_data);

         std::shared_ptr<pdat::CellData<double> > nonsymm_cell_diffs_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_quat_nonsymm_diffs_cell_id)));
         assert(nonsymm_cell_diffs_data);

         const hier::Box& box = patch->getBox();
         pdat::CellIterator iend(pdat::CellGeometry::end(box));
         for (pdat::CellIterator i(pdat::CellGeometry::begin(box)); i != iend;
              ++i) {
            const pdat::CellIndex ccell = *i;
            const pdat::SideIndex xside(ccell, pdat::SideIndex::X,
                                        pdat::SideIndex::Lower);
            const pdat::SideIndex yside(ccell, pdat::SideIndex::Y,
                                        pdat::SideIndex::Lower);
#if (NDIM == 3)
            const pdat::SideIndex zside(ccell, pdat::SideIndex::Z,
                                        pdat::SideIndex::Lower);
#endif
            for (int q = 0; q < d_qlen; q++) {
               int d = symm_depth_offset + q;
               (*cell_diffs_data)(ccell, q + 0 * d_qlen) =
                   (*diff_data)(xside, d);
               (*cell_diffs_data)(ccell, q + 1 * d_qlen) =
                   (*diff_data)(yside, d);
#if (NDIM == 3)
               (*cell_diffs_data)(ccell, q + 2 * d_qlen) =
                   (*diff_data)(zside, d);
#endif
               d = nonsymm_depth_offset + q;
               (*nonsymm_cell_diffs_data)(ccell, q + 0 * d_qlen) =
                   (*diff_data)(xside, d);
               (*nonsymm_cell_diffs_data)(ccell, q + 1 * d_qlen) =
                   (*diff_data)(yside, d);
#if (NDIM == 3)
               (*nonsymm_cell_diffs_data)(ccell, q + 2 * d_qlen) =
                   (*diff_data)(zside, d);
#endif
            }
         }
      }
   }
}

//-----------------------------------------------------------------------

// Computes gradients at cell centers. This will always compute the
// symmetry-aware gradient if symmetry is on.

// diff_id is SideData with ghosts, depth of NDIM*d_qlen
// grad_id is CellData with no ghosts

void QuatModel::computeQuatGradCell(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& quat_diffs_id,
    int& grad_cell_id, const double time, const CACHE_TYPE cache)
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if (time == old_time && cache == CACHE) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> patch_level =
          hierarchy->getPatchLevel(ln);

      computeQuatGradCell(patch_level, quat_diffs_id, grad_cell_id, time);
   }
}

void QuatModel::computeQuatGradCell(
    const std::shared_ptr<hier::PatchLevel> level, int& quat_diffs_id,
    int& grad_cell_id, const double time)
{
   (void)time;

   if (quat_diffs_id < 0) quat_diffs_id = d_quat_diffs_id;
   if (grad_cell_id < 0) grad_cell_id = d_quat_grad_cell_id;
   assert(quat_diffs_id >= 0);
   assert(grad_cell_id >= 0);

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      std::shared_ptr<hier::Patch> patch = *p;

      std::shared_ptr<pdat::SideData<double> > diff_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(quat_diffs_id)));
      assert(diff_data);
      assert(diff_data->getGhostCellWidth()[0] > 0);

      std::shared_ptr<pdat::CellData<double> > grad_cell_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(grad_cell_id)));
      assert(grad_cell_data);
      assert(grad_cell_data->getGhostCellWidth() ==
             hier::IntVector(tbox::Dimension(NDIM), 0));
      assert(grad_cell_data->getDepth() == NDIM * d_qlen);

      std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
          SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                                 hier::PatchGeometry>(
              patch->getPatchGeometry()));
      TBOX_ASSERT(patch_geom);

      const double* dx = patch_geom->getDx();

      if (d_symmetry_aware) {

         std::shared_ptr<pdat::SideData<int> > rotation_index(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<int>, hier::PatchData>(
                 patch->getPatchData(d_quat_symm_rotation_id)));
         assert(rotation_index);

         computeQGrad(diff_data, grad_cell_data, dx, true, rotation_index);

      } else {

         computeQGrad(diff_data, grad_cell_data, dx, false, nullptr);
      }
   }
}

//-----------------------------------------------------------------------

// Computes gradients at side.

// diff_id is SideData with ghosts, depth of NDIM*d_qlen
// grad_id is SideData with no ghosts

void QuatModel::computeQuatGradSide(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& quat_diffs_id,
    int& grad_side_id, const double time, const CACHE_TYPE cache)
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if (time == old_time && cache == CACHE) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> patch_level =
          hierarchy->getPatchLevel(ln);

      computeQuatGradSide(patch_level, quat_diffs_id, grad_side_id, time);
   }
}

/*
 * Compute gradients, using differences
 */
void QuatModel::computeQuatGradSide(
    const std::shared_ptr<hier::PatchLevel> level, int& quat_diffs_id,
    int& grad_side_id, const double time)
{
   (void)time;

   // tbox::pout<<"computeQuatGradSide()"<<endl;

   if (quat_diffs_id < 0) quat_diffs_id = d_quat_diffs_id;
   if (grad_side_id < 0) grad_side_id = d_quat_grad_side_id;
   assert(grad_side_id >= 0);

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      std::shared_ptr<hier::Patch> patch = *p;

      std::shared_ptr<pdat::SideData<double> > diff_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(quat_diffs_id)));
      assert(diff_data);

      std::shared_ptr<pdat::SideData<double> > grad_side_data(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(grad_side_id)));
      assert(grad_side_data);
      assert(grad_side_data->getDepth() == NDIM * d_qlen);

      std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
          SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                                 hier::PatchGeometry>(
              patch->getPatchGeometry()));
      assert(patch_geom);

      const double* dx = patch_geom->getDx();

      if (d_symmetry_aware) {
         std::shared_ptr<pdat::SideData<int> > rotation_index(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<int>, hier::PatchData>(
                 patch->getPatchData(d_quat_symm_rotation_id)));
         assert(rotation_index);
         computeQGradSide(diff_data, grad_side_data, dx, true,
                          d_model_parameters.useIsotropicStencil(),
                          rotation_index);
      } else {
         computeQGradSide(diff_data, grad_side_data, dx, false,
                          d_model_parameters.useIsotropicStencil(), nullptr);
      }
   }
}

//-----------------------------------------------------------------------

// Computes gradients at cell centers.

// grad_cell_id is CellData with no ghosts, depth of NDIM*d_qlen
// grad_modulus_id is CellData with no ghosts, depth of 1

void QuatModel::computeQuatGradModulus(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& grad_cell_id,
    int& grad_modulus_id, const double time, const CACHE_TYPE cache)
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if (time == old_time && cache != FORCE) return;
   old_time = time;

   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> patch_level =
          hierarchy->getPatchLevel(ln);

      computeQuatGradModulus(patch_level, grad_cell_id, grad_modulus_id, time);
   }
}

void QuatModel::computeQuatGradModulusFromSides(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& grad_side_id,
    int& grad_modulus_id, const double time, const CACHE_TYPE cache)
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if (time == old_time && cache != FORCE) return;
   old_time = time;

   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> patch_level =
          hierarchy->getPatchLevel(ln);

      computeQuatGradModulusFromSides(patch_level, grad_side_id,
                                      grad_modulus_id, time);
   }
}

void QuatModel::computeQuatGradModulus(
    const std::shared_ptr<hier::PatchLevel> level, int& grad_cell_id,
    int& grad_modulus_id, const double time)
{
   (void)time;

   if (grad_cell_id < 0) grad_cell_id = d_quat_grad_cell_id;
   if (grad_modulus_id < 0) grad_modulus_id = d_quat_grad_modulus_id;

   d_quat_grad_modulus_strategy->computeQuatGradModulus(level, grad_cell_id,
                                                        grad_modulus_id);
}

void QuatModel::computeQuatGradModulusFromSides(
    const std::shared_ptr<hier::PatchLevel> level, int& grad_side_id,
    int& grad_modulus_id, const double time)
{
   (void)time;

   if (grad_side_id < 0) grad_side_id = d_quat_grad_side_id;
   assert(grad_side_id >= 0);
   if (grad_modulus_id < 0) grad_modulus_id = d_quat_grad_modulus_id;

   d_quat_grad_modulus_strategy->computeQuatGradModulusFromSides(
       level, grad_side_id, grad_modulus_id);
}

//=======================================================================

void QuatModel::normalizeQuat(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const int quat_id)
{
   if (quat_id < 0) return;
   if (d_qlen == 1) return;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> patch_level =
          hierarchy->getPatchLevel(ln);

      normalizeQuat(patch_level, quat_id);
   }
}

void QuatModel::normalizeQuat(const std::shared_ptr<hier::PatchLevel> level,
                              const int quat_id)
{
   assert(quat_id >= 0);
   if (d_qlen == 1) return;

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      std::shared_ptr<hier::Patch> patch = *p;

      std::shared_ptr<pdat::CellData<double> > quat(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(quat_id)));
      assert(quat);
      const hier::Box& gbox = quat->getGhostBox();

      // issue:mew: Potentially replace this loop with fortran kernel.

      pdat::CellIterator iend(pdat::CellGeometry::end(gbox));
      for (pdat::CellIterator i(pdat::CellGeometry::begin(gbox)); i != iend;
           ++i) {
         pdat::CellIndex cell = *i;
         double qnorm2 = 0.;
         for (int q = 0; q < d_qlen; q++) {
            qnorm2 += (*quat)(cell, q) * (*quat)(cell, q);
         }
         const double invqnorm = 1. / sqrt(qnorm2);
         for (int q = 0; q < d_qlen; q++) {
            (*quat)(cell, q) *= invqnorm;
         }
      }
   }
}

//=======================================================================

// Computes phase mobility at cell centers.

// phase_id is CellData with ghosts
// mobility_id is CellData with no ghosts

void QuatModel::computeUniformPhaseMobility(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
    int& mobility_id, const double time, const CACHE_TYPE cache)
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if (time == old_time && cache == CACHE) return;
   old_time = time;

   if (phase_id < 0) phase_id = d_phase_scratch_id;
   if (mobility_id < 0) mobility_id = d_phase_mobility_id;

   math::HierarchyCellDataOpsReal<double> mathops(hierarchy);
   mathops.setToScalar(mobility_id, d_model_parameters.phase_mobility());
}

//-----------------------------------------------------------------------

void QuatModel::computePhaseMobilityPatch(
    const hier::Box& pbox, std::shared_ptr<pdat::CellData<double> > cd_temp,
    std::shared_ptr<pdat::CellData<double> > cd_mobility)
{
   double* ptr_temp = cd_temp->getPointer();
   double* ptr_m = cd_mobility->getPointer();

   const hier::Box& temp_gbox = cd_temp->getGhostBox();
   int imin_temp = temp_gbox.lower(0);
   int jmin_temp = temp_gbox.lower(1);
   int jp_temp = temp_gbox.numberCells(0);
   int kmin_temp = 0;
   int kp_temp = 0;
#if (NDIM == 3)
   kmin_temp = temp_gbox.lower(2);
   kp_temp = jp_temp * temp_gbox.numberCells(1);
#endif

   const hier::Box& m_gbox = cd_mobility->getGhostBox();
   int imin_m = m_gbox.lower(0);
   int jmin_m = m_gbox.lower(1);
   int jp_m = m_gbox.numberCells(0);
   int kmin_m = 0;
   int kp_m = 0;
#if (NDIM == 3)
   kmin_m = m_gbox.lower(2);
   kp_m = jp_m * m_gbox.numberCells(1);
#endif

   int imin = pbox.lower(0);
   int imax = pbox.upper(0);
   int jmin = pbox.lower(1);
   int jmax = pbox.upper(1);
   int kmin = 0;
   int kmax = 0;
#if (NDIM == 3)
   kmin = pbox.lower(2);
   kmax = pbox.upper(2);
#endif

   for (int kk = kmin; kk <= kmax; kk++) {
      for (int jj = jmin; jj <= jmax; jj++) {
         for (int ii = imin; ii <= imax; ii++) {

            const int idx_temp = (ii - imin_temp) + (jj - jmin_temp) * jp_temp +
                                 (kk - kmin_temp) * kp_temp;

            const int idx_m =
                (ii - imin_m) + (jj - jmin_m) * jp_m + (kk - kmin_m) * kp_m;

            double t = ptr_temp[idx_temp];

            double Rt = t * gas_constant_R_JpKpmol;

            ptr_m[idx_m] = d_model_parameters.phase_mobility() *
                           exp(-d_model_parameters.q0_phase_mobility() / Rt);
         }
      }
   }
}

//=======================================================================

// Computes eta mobility at cell centers.

// eta_id is CellData with ghosts
// mobility_id is CellData with no ghosts

void QuatModel::computeEtaMobility(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
    int& mobility_id, const double time, const CACHE_TYPE cache)
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if (time == old_time && cache == CACHE) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> patch_level =
          hierarchy->getPatchLevel(ln);

      computeEtaMobility(patch_level, phase_id, mobility_id, time);
   }
}

void QuatModel::computeEtaMobility(
    const std::shared_ptr<hier::PatchLevel> level, int& phase_id,
    int& mobility_id, const double time)
{
   (void)time;

   if (phase_id < 0) phase_id = d_phase_scratch_id;
   if (mobility_id < 0) mobility_id = d_eta_mobility_id;

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      std::shared_ptr<hier::Patch> patch = *p;

      const hier::Box& pbox = patch->getBox();

      std::shared_ptr<pdat::CellData<double> > temp_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_temperature_id)));
      assert(temp_data);

      std::shared_ptr<pdat::CellData<double> > phase_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(phase_id)));
      assert(phase_data);

      std::shared_ptr<pdat::CellData<double> > mobility_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(mobility_id)));
      assert(mobility_data);

      if (d_model_parameters.min_eta_mobility() ==
              d_model_parameters.eta_mobility() &&
          d_model_parameters.q0_eta_mobility() == 0.0) {

         mobility_data->fillAll(d_model_parameters.eta_mobility());

      } else {

         computeEtaMobilityPatch(pbox, temp_data, mobility_data, phase_data);
      }
   }
}

//-----------------------------------------------------------------------

void QuatModel::computeEtaMobilityPatch(
    const hier::Box& pbox, std::shared_ptr<pdat::CellData<double> > cd_temp,
    std::shared_ptr<pdat::CellData<double> > cd_mobility,
    std::shared_ptr<pdat::CellData<double> > cd_phi)
{
   double* ptr_temp = cd_temp->getPointer();
   double* ptr_m = cd_mobility->getPointer();
   double* ptr_phi = cd_phi->getPointer();

   const hier::Box& temp_gbox = cd_temp->getGhostBox();
   int imin_temp = temp_gbox.lower(0);
   int jmin_temp = temp_gbox.lower(1);
   int jp_temp = temp_gbox.numberCells(0);
   int kmin_temp = 0;
   int kp_temp = 0;
#if (NDIM == 3)
   kmin_temp = temp_gbox.lower(2);
   kp_temp = jp_temp * temp_gbox.numberCells(1);
#endif

   const hier::Box& m_gbox = cd_mobility->getGhostBox();
   int imin_m = m_gbox.lower(0);
   int jmin_m = m_gbox.lower(1);
   int jp_m = m_gbox.numberCells(0);
   int kmin_m = 0;
   int kp_m = 0;
#if (NDIM == 3)
   kmin_m = m_gbox.lower(2);
   kp_m = jp_m * m_gbox.numberCells(1);
#endif

   const hier::Box& phi_gbox = cd_phi->getGhostBox();
   int imin_phi = phi_gbox.lower(0);
   int jmin_phi = phi_gbox.lower(1);
   int jp_phi = phi_gbox.numberCells(0);
   int kmin_phi = 0;
   int kp_phi = 0;
#if (NDIM == 3)
   kmin_phi = phi_gbox.lower(2);
   kp_phi = jp_phi * phi_gbox.numberCells(1);
#endif

   int imin = pbox.lower(0);
   int imax = pbox.upper(0);
   int jmin = pbox.lower(1);
   int jmax = pbox.upper(1);
   int kmin = 0;
   int kmax = 0;
#if (NDIM == 3)
   kmin = pbox.lower(2);
   kmax = pbox.upper(2);
#endif

   for (int kk = kmin; kk <= kmax; kk++) {
      for (int jj = jmin; jj <= jmax; jj++) {
         for (int ii = imin; ii <= imax; ii++) {

            const int idx_temp = (ii - imin_temp) + (jj - jmin_temp) * jp_temp +
                                 (kk - kmin_temp) * kp_temp;

            const int idx_m =
                (ii - imin_m) + (jj - jmin_m) * jp_m + (kk - kmin_m) * kp_m;

            const int idx_phi = (ii - imin_phi) + (jj - jmin_phi) * jp_phi +
                                (kk - kmin_phi) * kp_phi;

            double t = ptr_temp[idx_temp];
            double phi = ptr_phi[idx_phi];

            double h_phi = INTERP_FUNC(phi, "p");
            double Rt = t * gas_constant_R_JpKpmol;

            double m = d_model_parameters.eta_mobility() *
                       exp(-d_model_parameters.q0_eta_mobility() / Rt);

            ptr_m[idx_m] =
                m - h_phi * (m - d_model_parameters.min_eta_mobility());
         }
      }
   }
}

//=======================================================================

// Computes Quaternion mobility at cell centers.

// phase_id is CellData with ghosts
// mobility_id is CellData with no ghosts

void QuatModel::computeQuatMobility(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
    int& mobility_id, const double time, const CACHE_TYPE cache)
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if (time == old_time && cache == CACHE) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> patch_level =
          hierarchy->getPatchLevel(ln);

      computeQuatMobility(patch_level, phase_id, mobility_id, time);
   }
}

void QuatModel::computeQuatMobility(
    const std::shared_ptr<hier::PatchLevel> level, int& phase_id,
    int& mobility_id, const double time)
{
   (void)time;
   assert(d_model_parameters.evolveQuat());

   if (phase_id < 0) phase_id = d_phase_scratch_id;
   if (mobility_id < 0) mobility_id = d_quat_mobility_id;

   double alt_scale_factor = d_model_parameters.quatMobilityScaleFactor();

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      std::shared_ptr<hier::Patch> patch = *p;

      const hier::Box& pbox = patch->getBox();
      const hier::Index& ifirst = pbox.lower();
      const hier::Index& ilast = pbox.upper();

      std::shared_ptr<pdat::CellData<double> > phase_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(phase_id)));
      assert(phase_data);

      std::shared_ptr<pdat::CellData<double> > mobility_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(mobility_id)));
      assert(mobility_data);

      QUATMOBILITY(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                   ifirst(2), ilast(2),
#endif
                   phase_data->getPointer(), phase_data->getGhostCellWidth()[0],
                   mobility_data->getPointer(),
                   mobility_data->getGhostCellWidth()[0],
                   d_model_parameters.quat_mobility(),
                   d_model_parameters.min_quat_mobility(),
                   d_model_parameters.quat_mobility_func_type().c_str(),
                   alt_scale_factor);
   }
}

//=======================================================================

// Computes Derivative of Quaternion mobility versus Phase at cell centers.

// phase_id is CellData with ghosts
// mobility_deriv_id is CellData with no ghosts

void QuatModel::computeQuatMobilityDeriv(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int& phase_id,
    int& mobility_deriv_id, const double time, const CACHE_TYPE cache)
{
   static double old_time = tbox::IEEE::getSignalingNaN();

   if (time == old_time && cache == CACHE) return;
   old_time = time;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> patch_level =
          hierarchy->getPatchLevel(ln);

      computeQuatMobilityDeriv(patch_level, phase_id, mobility_deriv_id, time);
   }
}

void QuatModel::computeQuatMobilityDeriv(
    const std::shared_ptr<hier::PatchLevel> level, int& phase_id,
    int& mobility_deriv_id, const double time)
{
   (void)time;

   if (phase_id < 0) phase_id = d_phase_scratch_id;
   assert(mobility_deriv_id >= 0);

   double alt_scale_factor = d_model_parameters.quatMobilityScaleFactor();

   for (hier::PatchLevel::Iterator p(level->begin()); p != level->end(); ++p) {
      std::shared_ptr<hier::Patch> patch = *p;
      const hier::Box& pbox = patch->getBox();
      const hier::Index& ifirst = pbox.lower();
      const hier::Index& ilast = pbox.upper();

      std::shared_ptr<pdat::CellData<double> > phase_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(phase_id)));
      assert(phase_data);

      std::shared_ptr<pdat::CellData<double> > mobility_deriv_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(mobility_deriv_id)));
      assert(mobility_deriv_data);
      assert(mobility_deriv_data->getGhostCellWidth() ==
             hier::IntVector(tbox::Dimension(NDIM), 0));

      QUATMOBILITYDERIV(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                        ifirst(2), ilast(2),
#endif
                        phase_data->getPointer(),
                        phase_data->getGhostCellWidth()[0],
                        mobility_deriv_data->getPointer(), 0,
                        d_model_parameters.quat_mobility(),
                        d_model_parameters.min_quat_mobility(),
                        d_model_parameters.quat_mobility_func_type().c_str(),
                        alt_scale_factor);
   }
}

//=======================================================================

void QuatModel::checkQuatNorm(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double tol)
{
   assert(tol >= 0.);

   if (!d_model_parameters.with_orientation()) return;
   if (d_qlen == 1) return;

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {

      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         std::shared_ptr<pdat::CellData<double> > y(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_quat_id)));
         const hier::Box& pbox = patch->getBox();

         pdat::CellIterator icend(pdat::CellGeometry::end(pbox));
         for (pdat::CellIterator ic(pdat::CellGeometry::begin(pbox));
              ic != icend; ++ic) {
            pdat::CellIndex cell = *ic;
            double qnorm2 = 0.;
            for (int q = 0; q < d_qlen; q++) {
               qnorm2 += (*y)(cell, q) * (*y)(cell, q);
            }
            double qnorm = sqrt(qnorm2);
            if (fabs(qnorm - 1.) > tol) {
               std::cerr << std::setprecision(10) << std::scientific;
               std::cerr << "WARNING: q norm=" << qnorm << std::endl;
            }
         }
      }
   }
}

/*
 * Set weight appropriate for computing std::vector norms.
 *
 * If you this function to set the weights used when you
 * SAMRAIVectorReal::addComponent, you can use the
 * std::vector norm functions of SAMRAIVectorReal, and
 * the weights will be used to blank out coarse grid
 * regions under fine grids.
 *
 * The weights computed are specific to the cell-centered
 * discretization used by this class.  The weight is equal
 * to the cell volume if the cell has not been refined,
 * and zero if it has.
 *
 * This function is state-independent.  All inputs are in
 * the argument list.
 *
 * hierarchy:   Hierarchy configuration to compute weights for
 * coarsest_ln: Coarsest level number.  Must be included
 *              in hierarchy.  Must not be greater than finest_ln.
 *              Default to 0.
 * finest_ln:   Finest level number.  Must be included
 *              in hierarchy.  Must not be less than coarsest_ln.
 *              Default to finest level in hierarchy.
 */
void QuatModel::computeVectorWeights(
    std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   assert(d_weight_id != -1);

   int coarsest_ln = 0;
   int finest_ln = hierarchy->getFinestLevelNumber();
   assert(finest_ln >= coarsest_ln);

   for (int ln = finest_ln; ln >= coarsest_ln; --ln) {

      // On every level, first assign cell volume to weight
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;
         std::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
             SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                                    hier::PatchGeometry>(
                 patch->getPatchGeometry()));
         assert(patch_geometry);

         const double* dx = patch_geometry->getDx();
         double cell_vol = dx[0];
         if (NDIM > 1) {
            cell_vol *= dx[1];
         }
         if (NDIM > 2) {
            cell_vol *= dx[2];
         }

         std::shared_ptr<pdat::CellData<double> > w(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_weight_id)));
         assert(w);
         w->fillAll(cell_vol);
      }
   }

   zeroOutLowerLevelWeights(hierarchy);
}

void QuatModel::computeVectorWeightsDiagnostics(
    std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   assert(d_weight_diagnostics_id != -1);

   int coarsest_ln = 0;
   int finest_ln = hierarchy->getFinestLevelNumber();
   assert(finest_ln >= coarsest_ln);

   // get bounds of physical domain
   const double* low = d_grid_geometry->getXLower();
   const double* up = d_grid_geometry->getXUpper();

   // create diagnostics domain defined by dlow, dup by removing
   // a fraction of the original domain
   double dd[NDIM];
   for (int d = 0; d < NDIM; d++)
      dd[d] = (up[d] - low[d]);
   assert(dd[0] > 0.);

   // compute size of part to remove
   for (int d = 0; d < NDIM; d++)
      dd[d] *= (1. - d_subdomain_diagnostics[d]);

   // remove half from each side
   double dlow[NDIM];
   for (int d = 0; d < NDIM; d++)
      dlow[d] = low[d] + 0.5 * dd[d];
   double dup[NDIM];
   for (int d = 0; d < NDIM; d++)
      dup[d] = up[d] - 0.5 * dd[d];
   tbox::plog << "Diagnostics subdomain: [" << dlow[0] << "," << dup[0]
              << "] x [" << dlow[1] << "," << dup[1] << "]"
#if NDIM > 2
              << " x [" << dlow[2] << "," << dup[2] << "]"
#endif
              << std::endl;

   for (int ln = finest_ln; ln >= coarsest_ln; --ln) {

      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;
         std::shared_ptr<geom::CartesianPatchGeometry> patch_geometry(
             SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                                    hier::PatchGeometry>(
                 patch->getPatchGeometry()));
         assert(patch_geometry);

         const hier::Box& pbox = patch->getBox();

         const double* dx = patch_geometry->getDx();
         double cell_vol = dx[0];
         if (NDIM > 1) {
            cell_vol *= dx[1];
         }
         if (NDIM > 2) {
            cell_vol *= dx[2];
         }
         assert(cell_vol > 0.);

         std::shared_ptr<pdat::CellData<double> > w(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_weight_diagnostics_id)));
         assert(w);

         w->fillAll(cell_vol);
         // zero out cells outside diagnostics domain
         pdat::CellIterator iend(pdat::CellGeometry::end(pbox));
         for (pdat::CellIterator i(pdat::CellGeometry::begin(pbox)); i != iend;
              ++i) {
            pdat::CellIndex cell = *i;
            for (int d = 0; d < NDIM; d++) {
               double x = low[d] + cell(d) * dx[d];

               if (x < dlow[d]) (*w)(cell) = 0.;
               if (x > dup[d]) (*w)(cell) = 0.;
            }
         }
      }
   }

   zeroOutLowerLevelWeights(hierarchy);
}

void QuatModel::zeroOutLowerLevelWeights(
    std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   int coarsest_ln = 0;
   int finest_ln = hierarchy->getFinestLevelNumber();

   // On all but the finest level, assign 0 to
   // weight to cells covered by finer cells.
   for (int ln = finest_ln - 1; ln >= coarsest_ln; --ln) {

      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      // First get the boxes that describe index space of the next finer
      // level and coarsen them to describe corresponding index space
      // at this level.
      std::shared_ptr<hier::PatchLevel> next_finer_level =
          hierarchy->getPatchLevel(ln + 1);
      hier::BoxContainer coarsened_boxes = next_finer_level->getBoxes();
      hier::IntVector coarsen_ratio = next_finer_level->getRatioToLevelZero();
      coarsen_ratio /= level->getRatioToLevelZero();
      coarsened_boxes.coarsen(coarsen_ratio);

      // Then set weight to 0 wherever there is
      // a nonempty intersection with the next finer level.
      // Note that all assignments are local.
      for (hier::PatchLevel::iterator p(level->begin()); p != level->end();
           ++p) {

         std::shared_ptr<hier::Patch> patch = *p;
         for (hier::BoxContainer::const_iterator i = coarsened_boxes.begin();
              i != coarsened_boxes.end(); ++i) {

            hier::Box coarse_box = *i;
            hier::Box intersection = coarse_box * (patch->getBox());
            if (!intersection.empty()) {
               std::shared_ptr<pdat::CellData<double> > w(
                   SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>,
                                          hier::PatchData>(
                       patch->getPatchData(d_weight_id)));
               w->fillAll(0.0, intersection);

            }  // assignment only in non-empty intersection
         }     // loop over coarsened boxes from finer level
      }        // loop over patches in level
   }           // loop over levels
}

//=======================================================================

void QuatModel::evaluateEnergy(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time,
    double& total_energy, double& total_interface_e, double& total_eta_e,
    double& total_orient_e, double& total_qint_e, double& total_well_e,
    double& total_bulk_e, const bool gp)
{
   assert(d_weight_id != -1);
   if (d_model_parameters.with_visit_energy_output())
      assert(d_energy_diag_id != -1);

   total_energy = 0.;
   total_interface_e = 0.;
   total_eta_e = 0.;
   total_orient_e = 0.;
   total_qint_e = 0.;
   total_bulk_e = 0.;
   total_well_e = 0.;

   copyCurrentToScratch(hierarchy, time, d_all_refine_patch_strategy);

   if (d_model_parameters.evolveQuat()) {
      int diff_id = -1;
      d_quat_grad_strategy->computeDiffs(hierarchy, d_quat_scratch_id, diff_id,
                                         time, QuatGradStrategy::FORCE);
      assert(diff_id >= 0);

      // Compute gradients on cell faces
      d_quat_grad_strategy->computeGradSide(hierarchy, diff_id,
                                            d_quat_grad_side_id, time,
                                            QuatGradStrategy::FORCE);
   }
   // compute free energies function of phase, concentration and temperature

   if (d_model_parameters.with_concentration()) {
      if (d_model_parameters.partition_phase_concentration()) {
         // computeVelocity(hierarchy,ydot_phase_id);

         assert(d_partition_coeff_strategy != nullptr);

         d_partition_coeff_strategy->evaluate(hierarchy);
         if (d_model_parameters.needGhosts4PartitionCoeff())
            fillPartitionCoeffGhosts();
      }
      if (d_model_parameters.concentrationModelNeedsPhaseConcentrations()) {
         assert(d_phase_conc_strategy != nullptr);
         d_phase_conc_strategy->computePhaseConcentrations(
             hierarchy, d_temperature_scratch_id, d_phase_scratch_id,
             d_conc_scratch_id);
      }
   }

   if (d_free_energy_strategy) {
      d_free_energy_strategy->computeFreeEnergyLiquid(hierarchy,
                                                      d_temperature_id,
                                                      d_f_l_id, gp);

      d_free_energy_strategy->computeFreeEnergySolidA(hierarchy,
                                                      d_temperature_id,
                                                      d_f_a_id, gp);
      if (d_model_parameters.withPhaseB()) {
         d_free_energy_strategy->computeFreeEnergySolidB(hierarchy,
                                                         d_temperature_id,
                                                         d_f_b_id, gp);
      }
   }

   if (!d_model_parameters.use_FolchPlapp() && d_model_parameters.with_phase())
      d_energy_eval_strategy->evaluateEnergy(hierarchy, time, total_energy,
                                             total_interface_e, total_orient_e,
                                             total_qint_e, total_well_e,
                                             total_bulk_e, gp);

   if (d_model_parameters.withRBmotion()) {
      d_multiorderp_energy->evaluatePairEnergy(hierarchy);
      d_multiorderp_energy->printPairEnergy(tbox::plog);

      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
      mpi.Barrier();

      assert(d_rigid_body_forces);

      d_rigid_body_forces->printPairForces(tbox::plog);

      mpi.Barrier();
   }
}

//=======================================================================

void QuatModel::computePhaseConcentrations(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   assert(d_phase_conc_strategy != nullptr);

   t_phase_conc_timer->start();

#ifdef DEBUG_CHECK_ASSERTIONS
   math::HierarchyCellDataOpsReal<double> cellops(hierarchy);
   double maxphi = cellops.max(d_phase_scratch_id);
   assert(maxphi == maxphi);
   double minphi = cellops.min(d_phase_scratch_id);
   assert(maxphi >= 0.);
   assert(maxphi < 1.1);
   assert(minphi >= -0.1);
   assert(minphi <= 1.);
   double maxc = cellops.max(d_conc_scratch_id);
   assert(maxc == maxc);
#endif

   // tbox::pout<<"Evaluate k..."<<endl;
   if (d_model_parameters.with_partition_coeff()) {
      d_partition_coeff_strategy->evaluate(hierarchy);
      if (d_model_parameters.needGhosts4PartitionCoeff())
         fillPartitionCoeffGhosts();
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(cellops.max(d_phase_scratch_id) == cellops.max(d_phase_scratch_id));
#endif

   // phase concentrations are computed for ghost values too,
   // so data with ghosts is needed for phase, conc and temperature
   d_phase_conc_strategy->computePhaseConcentrations(hierarchy,
                                                     d_temperature_scratch_id,
                                                     d_phase_scratch_id,
                                                     d_conc_scratch_id);

   t_phase_conc_timer->stop();
}


//=======================================================================

void QuatModel::computeSymmetryRotations(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time)
{
   (void)time;

   assert(d_quat_scratch_id >= 0);
   assert(d_quat_symm_rotation_id >= 0);

   tbox::pout << "compute symmetry rotations..." << std::endl;

   // Fill ghosts of original quat data
   copyCurrentToScratch(d_patch_hierarchy, d_time, d_all_refine_patch_strategy);

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         const hier::Box& pbox = patch->getBox();
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast = pbox.upper();

         std::shared_ptr<pdat::CellData<double> > quat(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_quat_scratch_id)));
         assert(quat);
         assert(quat->getGhostCellWidth()[0] > 0);

         std::shared_ptr<pdat::SideData<int> > rotation_index(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<int>, hier::PatchData>(
                 patch->getPatchData(d_quat_symm_rotation_id)));
         assert(rotation_index);
         assert(rotation_index->getGhostCellWidth()[0] > 0);

         QUAT_SYMM_ROTATION(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                            ifirst(2), ilast(2),
#endif
                            quat->getPointer(), quat->getGhostCellWidth()[0],
                            d_qlen, rotation_index->getPointer(0),
                            rotation_index->getPointer(1),
#if (NDIM == 3)
                            rotation_index->getPointer(2),
#endif
                            rotation_index->getGhostCellWidth()[0]);

         if (d_model_parameters.with_extra_visit_output()) {
            assert(d_quat_symm_rotation_cell_id >= 0);
            std::shared_ptr<pdat::CellData<int> > cell_rot_data(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<int>, hier::PatchData>(
                    patch->getPatchData(d_quat_symm_rotation_cell_id)));
            assert(cell_rot_data);

            pdat::CellIterator iend(pdat::CellGeometry::end(pbox));
            for (pdat::CellIterator i(pdat::CellGeometry::begin(pbox));
                 i != iend; ++i) {
               const pdat::CellIndex ccell = *i;
               const pdat::SideIndex xside(ccell, pdat::SideIndex::X,
                                           pdat::SideIndex::Lower);
               const pdat::SideIndex yside(ccell, pdat::SideIndex::Y,
                                           pdat::SideIndex::Lower);
#if (NDIM == 3)
               const pdat::SideIndex zside(ccell, pdat::SideIndex::Z,
                                           pdat::SideIndex::Lower);
#endif
               (*cell_rot_data)(ccell, 0) = (*rotation_index)(xside);
               (*cell_rot_data)(ccell, 1) = (*rotation_index)(yside);
#if (NDIM == 3)
               (*cell_rot_data)(ccell, 2) = (*rotation_index)(zside);
#endif
            }
         }
      }
   }
}

//=======================================================================

void QuatModel::makeQuatFundamental(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time)
{
   (void)time;

   assert(d_quat_id >= 0);

   if (d_verbosity->notSilent()) {
      tbox::pout << "Setting fundamental orientation" << std::endl;
   }

   int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         const hier::Box& pbox = patch->getBox();
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast = pbox.upper();

         std::shared_ptr<pdat::CellData<double> > quat(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_quat_id)));
         assert(quat);
         const hier::Box& quat_gbox = quat->getGhostBox();
         const hier::Index& q_lower = quat_gbox.lower();
         const hier::Index& q_upper = quat_gbox.upper();

         QUAT_FUNDAMENTAL(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                          ifirst(2), ilast(2),
#endif
                          quat->getPointer(), q_lower[0], q_upper[0],
                          q_lower[1], q_upper[1],
#if (NDIM == 3)
                          q_lower[2], q_upper[2],
#endif
                          d_qlen);
      }
   }
}

//=======================================================================

double QuatModel::evaluateIntegralConcentration(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const int component)
{
   assert(d_weight_id != -1);

   if (!d_model_parameters.with_concentration()) return 0.;

   math::HierarchyCellDataOpsReal<double> mathops(hierarchy);

   int conc_id = d_conc_id;
   if (d_ncompositions > 1) {
      assert(d_work_id != -1);

      copyDepthCellData(hierarchy, d_work_id, 0, d_conc_id, component);
      conc_id = d_work_id;
   }

   double value = mathops.integral(conc_id, d_weight_id);

   return value;
}

//=======================================================================

double QuatModel::evaluateIntegralPhaseConcentration(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const int component)
{
   assert(d_weight_id != -1);

   if (!d_model_parameters.with_concentration()) return 0.;

   math::HierarchyCellDataOpsReal<double> mathops(hierarchy);

   int conc_id = d_conc_id;
   if (d_ncompositions > 1) {
      copyDepthCellData(hierarchy, d_work_id, 0, d_conc_id, component);
      conc_id = d_work_id;
   }

   // assumes d_phase_id has depth 1
   mathops.multiply(d_phase_scratch_id, conc_id, d_phase_id);

   double value = mathops.L1Norm(d_phase_scratch_id, d_weight_id);

   return value;
}

//=======================================================================

double QuatModel::evaluatePhaseFraction(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const int phase_id,
    const int depth)
{
   assert(d_weight_id != -1);

   math::HierarchyCellDataOpsReal<double> mathops(hierarchy);

   copyDepthCellData(hierarchy, d_work_id, 0, d_phase_id, depth);

   return mathops.integral(d_work_id, d_weight_id);
}

//=======================================================================

double QuatModel::evaluateVolumeSolid(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const int phase_id)
{
   assert(d_weight_id != -1);

   math::HierarchyCellDataOpsReal<double> mathops(hierarchy);

   return mathops.L1Norm(phase_id, d_weight_id);
}

//=======================================================================

void QuatModel::fillPartitionCoeffGhosts(void)
{
   assert(d_partition_coeff_id >= 0);
   assert(d_partition_coeff_scratch_id >= 0);
   if (!d_all_periodic)
      assert(d_partition_coeff_refine_patch_strategy != nullptr);

   // tbox::pout<<"QuatModel::fillPartitionCoeffGhosts"<<endl;

   xfer::RefineAlgorithm copy_to_scratch;

   std::shared_ptr<hier::RefineOperator> refine_op =
       d_grid_geometry->lookupRefineOperator(d_partition_coeff_var,
                                             "LINEAR_REFINE");

   copy_to_scratch.registerRefine(
       d_partition_coeff_scratch_id,  // destination
       d_partition_coeff_id,          // source
       d_partition_coeff_scratch_id,  // temporary work space
       refine_op);

   const int maxl = d_patch_hierarchy->getNumberOfLevels();

   for (int ln = 0; ln < maxl; ln++) {
      std::shared_ptr<hier::PatchLevel> level =
          d_patch_hierarchy->getPatchLevel(ln);

      copy_to_scratch
          .createSchedule(level, ln - 1, d_patch_hierarchy,
                          d_partition_coeff_refine_patch_strategy)
          ->fillData(d_time);
   }
}

//=======================================================================

void QuatModel::resetRefPhaseConcentrations()
{
   assert(d_conc_l_id >= 0);
   assert(d_conc_a_id >= 0);
   assert(d_conc_l_ref_id >= 0);
   assert(d_conc_a_ref_id >= 0);

   // tbox::pout << "QuatModel::resetRefPhaseConcentrations()" << std::endl;

   math::HierarchyCellDataOpsReal<double> cellops(d_patch_hierarchy);

   cellops.copyData(d_conc_l_ref_id, d_conc_l_id, false);
   cellops.copyData(d_conc_a_ref_id, d_conc_a_id, false);
   if (d_model_parameters.withPhaseB())
      cellops.copyData(d_conc_b_ref_id, d_conc_b_id, false);
}

//=======================================================================
/*
void QuatModel::setPhaseConcentrationsToEquilibrium(const double* const ceq)
{
   assert(d_conc_l_id >= 0);
   assert(d_conc_a_id >= 0);

   tbox::pout << "QuatModel::setPhaseConcentrationsToEquilibrium(ceq)"
              << std::endl;
   int offset = d_ncompositions;
   for (int ln = 0; ln <= d_patch_hierarchy->getFinestLevelNumber(); ln++) {
      std::shared_ptr<hier::PatchLevel> level =
          d_patch_hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
           ip++) {
         std::shared_ptr<hier::Patch> patch = *ip;

         std::shared_ptr<pdat::CellData<double> > concl(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_conc_l_id)));
         assert(concl);

         for (int i = 0; i < d_ncompositions; i++)
            concl->fill(ceq[i], i);

         std::shared_ptr<pdat::CellData<double> > conca(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_conc_a_id)));
         assert(conca);

         for (int i = 0; i < d_ncompositions; i++)
            conca->fill(ceq[offset + i], i);
         if (d_model_parameters.with_three_phases()) {
            std::shared_ptr<pdat::CellData<double> > concb(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_conc_b_id)));
            assert(concb);

            for (int i = 0; i < d_ncompositions; i++)
               concb->fill(ceq[2 * offset + i], i);
         }
      }
   }
}
*/
//=======================================================================

void QuatModel::setRefPhaseConcentrationsToEquilibrium(const double* const ceq)
{
   assert(d_conc_l_ref_id >= 0);
   assert(d_conc_a_ref_id >= 0);

   tbox::pout << "QuatModel::setRefPhaseConcentrationsToEquilibrium(ceq)"
              << std::endl;

   math::HierarchyCellDataOpsReal<double> cellops(d_patch_hierarchy);
   cellops.setToScalar(d_conc_l_ref_id, ceq[0], false);
   cellops.setToScalar(d_conc_a_ref_id, ceq[1], false);
   if (d_model_parameters.withPhaseB())
      cellops.setToScalar(d_conc_b_ref_id, ceq[2], false);
}

//=======================================================================

// evaluate velocity field at every cell
void QuatModel::computeVelocity(std::shared_ptr<hier::Patch> patch,
                                int phi_dot_id)
{
   assert(d_phase_grad_cell_id >= 0);
   assert(d_velocity_id >= 0);
   assert(phi_dot_id >= 0);

   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   std::shared_ptr<pdat::CellData<double> > grad_cell_data(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_phase_grad_cell_id)));

   std::shared_ptr<pdat::CellData<double> > phi_dot_data(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(phi_dot_id)));

   std::shared_ptr<pdat::CellData<double> > velocity_data(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_velocity_id)));

   assert(grad_cell_data);
   assert(phi_dot_data);
   assert(velocity_data);

   assert(grad_cell_data->getGhostCellWidth()[0] == 0);
   assert(phi_dot_data->getGhostCellWidth()[0] == 0);
   assert(velocity_data->getGhostCellWidth()[0] == 0);

   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch->getPatchGeometry()));
   TBOX_ASSERT(patch_geom);

   const double* dx = patch_geom->getDx();
   const double threshold = 0.02 / dx[0];
   // tbox::pout<<"QuatModel::computeVelocity() with threshold
   // "<<threshold<<endl;

   VELOCITY(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
            ifirst(2), ilast(2),
#endif
            threshold, grad_cell_data->getPointer(0),
            grad_cell_data->getPointer(1),
#if (NDIM == 3)
            grad_cell_data->getPointer(2),
#endif
            phi_dot_data->getPointer(), velocity_data->getPointer());
}

//=======================================================================

void QuatModel::computeVelocity(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, int phi_dot_id)
{
   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;
         computeVelocity(patch, phi_dot_id);
      }
   }
}

//=======================================================================

double QuatModel::computeThermalEnergy(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   math::HierarchyCellDataOpsReal<double> cellops(hierarchy);

   double lenergy = 0.;
   if (d_model_parameters.with_phase()) {
      lenergy = cellops.integral(d_phase_id, d_weight_id);
      lenergy *= (-1. * d_model_parameters.latent_heat());
   }
   // double
   // refenergy=d_model_parameters.rescale_factorT()*cellops.integral(d_cp_id,
   // d_weight_id ); if ( d_model_parameters.with_rescaled_temperature()
   // )refenergy/d_model_parameters.rescale_factorT();

   // store product cp*T in d_fl_id
   cellops.multiply(
       d_f_l_id, d_cp_id,
       d_temperature_id);  // rescaling of cp compensates rescaling of T

   double cenergy = cellops.integral(d_f_l_id, d_weight_id);

   return lenergy + cenergy;  //-refenergy;
}

//=======================================================================

void QuatModel::setDiffusionVisit(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   if (!d_model_parameters.with_extra_visit_output()) return;

   if (d_model_parameters.concRHSstrategyIsEBS()) {
      assert(d_visit_conc_pfm_diffusion_l_id >= 0);
      assert(d_conc_pfm_diffusion_l_id >= 0);

      sideToCell(hierarchy, d_visit_conc_pfm_diffusion_l_id, 0,
                 d_conc_pfm_diffusion_l_id, 0);
      sideToCell(hierarchy, d_visit_conc_pfm_diffusion_a_id, 0,
                 d_conc_pfm_diffusion_a_id, 0);
      if (d_conc_pfm_diffusion_b_id >= 0)
         sideToCell(hierarchy, d_visit_conc_pfm_diffusion_b_id, 0,
                    d_conc_pfm_diffusion_b_id, 0);
   } else {
      assert(d_visit_conc_pfm_diffusion_id >= 0);
      assert(d_conc_pfm_diffusion_id[0] >= 0);
      sideToCell(hierarchy, d_visit_conc_pfm_diffusion_id, 0,
                 d_conc_pfm_diffusion_id[0], 0);
   }
}
