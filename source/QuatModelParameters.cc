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
#include "QuatModelParameters.h"
#include "PhysicalConstants.h"

const double gas_constant_R_JpKpmol = GASCONSTANT_R_JPKPMOL;

static double def_val = tbox::IEEE::getSignalingNaN();

static void readSpeciesCP(std::shared_ptr<tbox::Database> cp_db,
                          std::map<short, double>& cp)
{
   double tmp = cp_db->getDouble("a");
   cp.insert(std::pair<short, double>(0, tmp));
   if (cp_db->keyExists("b")) {
      tmp = cp_db->getDouble("b");
      cp.insert(std::pair<short, double>(1, tmp));
   }
   if (cp_db->keyExists("dm2")) {
      tmp = cp_db->getDouble("dm2");
      cp.insert(std::pair<short, double>(-2, tmp));
   }
}

QuatModelParameters::QuatModelParameters() : d_moving_frame_velocity(def_val)
{
   d_norderp = 1;
   d_norderp_A = -1;
   d_norderp_B = -1;
   d_sigma = -1.;
   d_delta = -1.;
   d_gamma = 1.5;
   d_H_parameter = def_val;
   d_epsilon_phase = def_val;
   d_epsilon_anisotropy = def_val;
   d_epsilon_eta = def_val;
   d_epsilon_q = def_val;
   d_noise_amplitude = def_val;
   d_phase_mobility = def_val;
   d_eta_mobility = def_val;
   d_quat_mobility = def_val;
   d_min_eta_mobility = def_val;
   d_min_quat_mobility = def_val;
   d_max_quat_mobility = def_val;
   d_exp_scale_quat_mobility = def_val;
   d_q0_phase_mobility = def_val;
   d_q0_eta_mobility = def_val;
   d_phase_well_scale = def_val;
   d_eta_well_scale = def_val;
   d_free_energy_liquid = def_val;
   d_free_energy_solid_A = def_val;
   d_free_energy_solid_B = def_val;
   d_molar_volume_liquid = def_val;
   d_molar_volume_solid_A = def_val;
   d_molar_volume_solid_B = def_val;
   d_D_liquid = def_val;
   d_D_solid_A = def_val;
   d_D_solid_B = def_val;
   d_D0_AB = def_val;
   d_D0_LA = def_val;
   d_D0_LB = def_val;
   d_Q0_liquid = def_val;
   d_Q0_solid_A = def_val;
   d_Q0_solid_B = def_val;
   d_Q0_AB = def_val;
   d_Q0_LA = def_val;
   d_Q0_LB = def_val;

   // defaults values of -1 means they were not set
   d_concL_ref = -1.;
   d_concA_ref = -1.;
   d_concB_ref = -1.;

   d_conc_mobility = def_val;
   d_thermal_diffusivity = def_val;
   d_latent_heat = def_val;
   d_meltingT = def_val;
   d_interface_mobility = def_val;

   d_well_bias_alpha = 0.;
   d_well_bias_gamma = def_val;
   d_liquidus_slope = def_val;
   d_average_concentration = def_val;
   d_stochio_cB = -1.;

   for (short i = 0; i < 3; i++) {
      d_ceq0_solidA[i] = def_val;
      d_ceq0_liquid[i] = def_val;
   }

   d_avg_func_type = "";
   d_diffq_avg_func_type = "";
   d_phase_well_func_type = "";

   d_conc_interp_func_type = Thermo4PFM::ConcInterpolationType::UNDEFINED;
   d_energy_interp_func_type = Thermo4PFM::EnergyInterpolationType::UNDEFINED;
   d_conc_three_args_interp_func_type =
       ConcThreeArgsInterpolationType::UNDEFINED;
   d_energy_three_args_interp_func_type =
       EnergyThreeArgsInterpolationType::UNDEFINED;
   d_diffusion_interp_type = DiffusionInterpolationType::UNDEFINED;
   d_eta_well_func_type = "";
   d_eta_interp_func_type = Thermo4PFM::EnergyInterpolationType::UNDEFINED;
   d_eta_well_func_type = "";
   d_quat_mobility_func_type = "";

   d_conc_avg_func_type = "";

   d_conc_model = ConcModel::UNDEFINED;
   d_conc_rhs_strategy = ConcRHSstrategy::UNKNOWN;
   d_conc_diffusion_type = ConcDiffusionType::UNDEFINED;
   d_temperature_type = TemperatureType::UNDEFINED;

   d_with_concentration = false;
   d_ncompositions = -1;

   d_with_phase = true;
   d_with_third_phase = false;
   d_FolchPlapp = false;
   d_with_heat_equation = false;
   d_with_steady_temperature = false;
   d_with_gradT = false;
   d_with_antitrapping = false;
   d_with_bias_well = false;
   d_with_rescaled_temperature = false;

   d_partition_coeff = "";
   d_phase_concentration_model = "";
   d_vd = -1.;
   d_with_velocity = false;
   d_use_diffs_to_compute_flux = false;
}

//=======================================================================

void QuatModelParameters::readMolarVolumes(std::shared_ptr<tbox::Database> db)
{
   bool data_read = false;

   if (db->keyExists("molar_volume")) {
      d_molar_volume_solid_A = db->getDouble("molar_volume");
      d_molar_volume_solid_B = d_molar_volume_solid_A;
      d_molar_volume_liquid = d_molar_volume_solid_A;
      data_read = true;
   } else if (db->keyExists("molar_volume_solid_A")) {
      d_molar_volume_liquid = db->getDouble("molar_volume_liquid");
      d_molar_volume_solid_A = db->getDouble("molar_volume_solid_A");
      if (db->keyExists("molar_volume_solid_B"))
         d_molar_volume_solid_B = db->getDouble("molar_volume_solid_B");
      else
         d_molar_volume_solid_B = d_molar_volume_solid_A;
      data_read = true;
   } else if (db->keyExists("ConcentrationModel")) {
      readMolarVolumes(db->getDatabase("ConcentrationModel"));
   }

   if (data_read)
      tbox::plog << "Molar volume liquid [m^3/mol]: " << d_molar_volume_liquid
                 << std::endl;
}

//=======================================================================

void QuatModelParameters::readNumberSpecies(
    std::shared_ptr<tbox::Database> conc_db)
{
   int nspecies = conc_db->getIntegerWithDefault("nspecies", 2);
   d_ncompositions = nspecies - 1;
}

void QuatModelParameters::readEquilibriumCompositions(
    std::shared_ptr<tbox::Database> conc_db)
{
   std::string db_name("Equilibrium");
   if (conc_db->keyExists(db_name)) {
      std::shared_ptr<tbox::Database> db = conc_db->getDatabase(db_name);

      d_Tref = db->getDouble("Tref");

      std::shared_ptr<tbox::Database> dbl = db->getDatabase("liquid");
      d_ceq0_liquid[0] = dbl->getDouble("c");
      d_ceq0_liquid[1] = dbl->getDouble("b");
      d_ceq0_liquid[2] = dbl->getDouble("a");

      std::shared_ptr<tbox::Database> dba = db->getDatabase("solidA");
      d_ceq0_solidA[0] = dba->getDouble("c");
      d_ceq0_solidA[1] = dba->getDouble("b");
      d_ceq0_solidA[2] = dba->getDouble("a");
   }
}

//=======================================================================

void QuatModelParameters::readConcDB(std::shared_ptr<tbox::Database> conc_db)
{
   d_with_concentration = true;

   // if concentration is ON, it means we have at least two species
   assert(d_ncompositions > 0);

   std::string conc_model = conc_db->getStringWithDefault("model", "undefined");
   if (conc_model.compare("calphad") == 0) {
      d_conc_model = ConcModel::CALPHAD;
   } else if (conc_model.compare("quadratic") == 0) {
      d_conc_model = ConcModel::QUADRATIC;
   } else if (conc_model.compare("parabolic") == 0) {
      d_conc_model = ConcModel::PARABOLIC;
   } else if (conc_model.compare("linear") == 0) {
      d_conc_model = ConcModel::LINEAR;
   } else if (conc_model.compare("independent") == 0) {
      d_conc_model =
          ConcModel::INDEPENDENT;  // energy independent of composition
   } else if (conc_model.compare("dilute") == 0) {
      d_conc_model = ConcModel::KKSdilute;
   } else if (conc_model.compare("cahn_hilliard") == 0) {
      d_conc_model = ConcModel::CahnHilliard;
   } else if (conc_model.compare("wang_sintering") == 0) {
      d_conc_model = ConcModel::WangSintering;
   } else {
      tbox::plog << "conc_model = " << conc_model << std::endl;
      TBOX_ERROR("Error: unknown concentration model in QuatModelParameters");
   }

   {
      std::string conc_rhs_strategy =
          conc_db->getStringWithDefault("rhs_form", "kks");
      if (d_conc_model == ConcModel::WangSintering) {
         d_conc_rhs_strategy = ConcRHSstrategy::WangSintering;
      } else {
         if (conc_rhs_strategy.compare("kks") == 0) {
            d_conc_rhs_strategy = ConcRHSstrategy::KKS;
         } else if (conc_rhs_strategy.compare("ebs") == 0) {
            d_conc_rhs_strategy = ConcRHSstrategy::EBS;
            // stencil type can be regular, isotropic, or order4
            d_ebs_stencil_type =
                conc_db->getStringWithDefault("ebs_stencil", "regular");
            if (d_ebs_stencil_type != "regular" &&
                d_ebs_stencil_type != "isotropic" &&
                d_ebs_stencil_type != "order4") {
               tbox::plog << "EBS stencil: " << d_ebs_stencil_type << std::endl;
               TBOX_ERROR("Error: unknown stencil type for EBS");
            }
         } else if (conc_rhs_strategy.compare("cahn_hilliard") == 0) {
            d_conc_rhs_strategy = ConcRHSstrategy::CahnHilliard;
            readCahnHilliard(conc_db);
         } else if (conc_rhs_strategy.compare("spinodal") == 0) {
            d_conc_rhs_strategy = ConcRHSstrategy::SPINODAL;
         } else if (conc_rhs_strategy[0] == 'u' ||
                    conc_rhs_strategy[0] == 'B' ||
                    conc_rhs_strategy[0] == 'b') {
            tbox::plog << "Using Beckermann's model" << std::endl;
            d_conc_rhs_strategy = ConcRHSstrategy::Beckermann;
         } else {
            TBOX_ERROR("Error: unknown concentration r.h.s. strategy");
         }
      }
   }

   // default setup so that older inputs files need not to be changed
   std::string default_concdiff_type =
       d_conc_rhs_strategy == ConcRHSstrategy::EBS ? "composition_dependent"
                                                   : "temperature_dependent";
   std::string conc_diffusion_strategy =
       conc_db->getStringWithDefault("diffusion_type", default_concdiff_type);
   if (conc_diffusion_strategy.compare("composition_dependent") == 0) {
      d_conc_diffusion_type = ConcDiffusionType::CTD;
   } else if (conc_diffusion_strategy.compare("mobility") == 0) {
      d_conc_diffusion_type = ConcDiffusionType::MOB;
   } else if (conc_diffusion_strategy.compare("temperature_dependent") == 0) {
      d_conc_diffusion_type = ConcDiffusionType::TD;
   }

   if (d_conc_rhs_strategy == ConcRHSstrategy::Beckermann) {
      tbox::plog << "Read diffusion constants for Beckermann's model"
                 << std::endl;
      d_D_liquid = conc_db->getDouble("D_liquid");
      d_D_solid_A = conc_db->getDouble("D_solid_A");
   }
   if (d_conc_diffusion_type == ConcDiffusionType::TD) {
      tbox::plog << "Read T-dependent diffusions" << std::endl;
      d_D_liquid = conc_db->getDouble("D_liquid");
      d_Q0_liquid = conc_db->getDoubleWithDefault("Q0_liquid", 0.);

      if (!withPhaseB()) {
         tbox::plog << "No phase B to read..." << std::endl;
         if (conc_db->keyExists("D_solid_A"))
            d_D_solid_A = conc_db->getDouble("D_solid_A");
         else
            d_D_solid_A = conc_db->getDouble("D_solid");
         if (conc_db->keyExists("Q0_solid_A"))
            d_Q0_solid_A = conc_db->getDoubleWithDefault("Q0_solid_A", 0.);
         else
            d_Q0_solid_A = conc_db->getDoubleWithDefault("Q0_solid", 0.);
      } else {
         tbox::plog << "Read phases A and B coefficients..." << std::endl;
         d_D_solid_A = conc_db->getDouble("D_solid_A");
         d_D_solid_B = conc_db->getDouble("D_solid_B");
         // optional diffusion values
         d_Q0_solid_A = conc_db->getDoubleWithDefault("Q0_solid_A", 0.);
         d_Q0_solid_B = conc_db->getDoubleWithDefault("Q0_solid_B", 0.);
         // AB interface diffusion 0. by default
         d_D0_LB = conc_db->getDoubleWithDefault("D0_LB", 0.);
         d_Q0_LB = conc_db->getDoubleWithDefault("Q0_LB", 0.);
         d_D0_AB = conc_db->getDoubleWithDefault("D0_AB", 0.);
         d_Q0_AB = conc_db->getDoubleWithDefault("Q0_AB", 0.);
         d_D0_BB = conc_db->getDoubleWithDefault("D0_BB", 0.);
         d_Q0_BB = conc_db->getDoubleWithDefault("Q0_BB", 0.);
      }
      d_D0_LA = conc_db->getDoubleWithDefault("D0_LA", 0.);
      d_Q0_LA = conc_db->getDoubleWithDefault("Q0_LA", 0.);
      d_D0_AA = conc_db->getDoubleWithDefault("D0_AA", 0.);
      d_Q0_AA = conc_db->getDoubleWithDefault("Q0_AA", 0.);
   }
   d_conc_mobility = conc_db->getDoubleWithDefault("mobility", 1.);

   d_conc_avg_func_type =
       conc_db->getStringWithDefault("avg_func_type", d_avg_func_type);
   if (d_conc_avg_func_type[0] != 'a' && d_conc_avg_func_type[0] != 'h') {
      TBOX_ERROR("Error: invalid value for avg_func_type");
   }

   if (conc_db->keyExists("gradT_Q0")) {
      d_with_gradT = true;
      // read heat of transport
      d_Q_heat_transport.resize(2);
      d_Q_heat_transport[0] = conc_db->getDouble("gradT_Q0");
      d_Q_heat_transport[1] = conc_db->getDouble("gradT_Q1");
   } else {
      d_with_gradT = false;
   }

   if (conc_db->keyExists("concL_ref")) {
      d_concL_ref = conc_db->getDouble("concL_ref");
   }
   if (conc_db->keyExists("concA_ref")) {
      d_concA_ref = conc_db->getDouble("concA_ref");
   }
   if (conc_db->keyExists("concB_ref")) {
      d_concB_ref = conc_db->getDouble("concB_ref");
   }

   d_with_antitrapping = conc_db->getBoolWithDefault("antitrapping", false);

   d_grand_potential = conc_db->getBoolWithDefault("gc", false);

   d_partition_coeff = conc_db->getStringWithDefault("partition_coeff", "none");
   tbox::plog << "Partition coefficient type: " << d_partition_coeff
              << std::endl;

   if (d_partition_coeff.compare("Aziz") == 0) {
      d_vd = conc_db->getDouble("vd");

      d_with_velocity = true;

      // a value of -1 for k_eq means it needs to be computed on the fly
      d_keq = conc_db->getDoubleWithDefault("keq", -1.);
      tbox::plog << "Aziz partition coefficient with Keq: " << d_keq
                 << std::endl;
   }
   if (d_partition_coeff.compare("uniform") == 0) {
      d_keq = conc_db->getDouble("keq");
      tbox::plog << "Uniform Keq: " << d_keq << std::endl;
   }

   assert(d_partition_coeff.compare("Aziz") == 0 ||
          d_partition_coeff.compare("uniform") == 0 ||
          d_partition_coeff.compare("none") == 0);

   std::string default_model = "none";
   if (d_conc_model == ConcModel::CALPHAD ||
       d_conc_model == ConcModel::QUADRATIC ||
       d_conc_model == ConcModel::PARABOLIC ||
       d_conc_model == ConcModel::KKSdilute)
      default_model = "kks";
   d_phase_concentration_model =
       conc_db->getStringWithDefault("phase_concentration_model",
                                     default_model);
   assert(d_phase_concentration_model.compare("none") == 0 ||
          d_phase_concentration_model.compare("kks") == 0 ||
          d_phase_concentration_model.compare("partition") == 0);
   tbox::plog << "phase_concentration_model: " << d_phase_concentration_model
              << std::endl;

   if (d_phase_concentration_model.compare("partition") == 0)
      assert(d_partition_coeff.compare("none") != 0);

   if (d_conc_model == ConcModel::LINEAR) {
      assert(d_meltingT == d_meltingT);

      d_liquidus_slope = conc_db->getDoubleWithDefault("liquidus_slope", 0.);
      if (fabs(d_liquidus_slope) > 0.)
         d_average_concentration = conc_db->getDouble("average_concentration");
      d_well_bias_alpha = conc_db->getDouble("alpha");
      if (d_well_bias_alpha > 0.) {
         d_with_bias_well = true;
      }
      d_well_bias_gamma = conc_db->getDoubleWithDefault("gamma", -1.);
      if (d_well_bias_gamma < 0.) {
         d_bias_well_beckermann = true;
      } else {
         d_bias_well_beckermann = false;
      }
      if (d_with_rescaled_temperature) {
         tbox::plog << "Rescale liquidus_slope and gamma..." << std::endl;
         d_liquidus_slope /= d_rescale_factorT;
         d_well_bias_gamma *= d_rescale_factorT;
      }
   }
   if (d_conc_model == ConcModel::KKSdilute) {
      readDiluteAlloy(conc_db);
   }

   if (d_conc_model == ConcModel::WangSintering) {
      readWangSintering(conc_db);
   }

   if (conc_db->keyExists("initc_in_phase")) {
      int nterms = conc_db->getArraySize("initc_in_phase");
      assert(nterms == 2 * d_ncompositions);
      d_initc_in_phase.resize(nterms);
      conc_db->getDoubleArray("initc_in_phase", &d_initc_in_phase[0], nterms);
   }

   d_init_phase_conc_eq =
       conc_db->getBoolWithDefault("init_phase_conc_eq", true);

   d_stochio_cB = conc_db->getDoubleWithDefault("stochio_cB", -1.);

   readEquilibriumCompositions(conc_db);
}

void QuatModelParameters::readDiluteAlloy(
    std::shared_ptr<tbox::Database> conc_db)
{
   d_liquidus_slope = conc_db->getDouble("liquidus_slope");
   d_keq = conc_db->getDouble("keq");
   d_meltingT = conc_db->getDouble("meltingT");  // in [K]
}

//=======================================================================

void QuatModelParameters::readVisitOptions(
    std::shared_ptr<tbox::Database> visit_db)
{
   d_extra_visit_output = visit_db->getBoolWithDefault("extra_output", false);
   d_visit_energy_output = visit_db->getBoolWithDefault("energy_output", false);
   d_visit_grain_output = visit_db->getBoolWithDefault("grain_output", false);
   d_with_velocity =
       visit_db->getBoolWithDefault("velocity_output", d_with_velocity);
   d_rhs_visit_output = visit_db->getBoolWithDefault("rhs_output", false);
}

//=======================================================================

void QuatModelParameters::readCahnHilliard(std::shared_ptr<tbox::Database> db)
{
   std::shared_ptr<tbox::Database> ch_db = db->getDatabase("CahnHilliard");

   d_CH_ca = ch_db->getDouble("ca");
   d_CH_cb = ch_db->getDouble("cb");
   d_CH_well_scale = ch_db->getDouble("well_scale");
   d_CH_kappa = ch_db->getDouble("kappa");
}

void QuatModelParameters::readWangSintering(std::shared_ptr<tbox::Database> db)
{
   if (d_sigma > 0. && d_delta > 0.) {
      tbox::plog << "Estimate WangSintering model parameters from sigma and "
                    "delta"
                 << std::endl;
      // see Biswas, Schwen, Tomar, J.Mater.Sci.(2018)
      d_WangSintering_A = 5. * d_sigma / (4. * d_delta);
      d_WangSintering_B = d_sigma / (4. * d_delta);
      d_WangSintering_beta_rho = 3. * d_delta * d_sigma;
      tbox::plog << "A        = " << d_WangSintering_A << std::endl;
      tbox::plog << "B        = " << d_WangSintering_B << std::endl;
      tbox::plog << "beta_rho = " << d_WangSintering_beta_rho << std::endl;
   } else {
      std::shared_ptr<tbox::Database> ws_db = db->getDatabase("WangSintering");

      d_WangSintering_A = ws_db->getDouble("A");
      d_WangSintering_B = ws_db->getDouble("B");
      d_WangSintering_beta_rho = ws_db->getDouble("beta");
   }
}

//=======================================================================

void QuatModelParameters::readTemperatureModel(
    std::shared_ptr<tbox::Database> model_db)
{
   std::shared_ptr<tbox::Database> temperature_db;
   std::string temperature_type = "";
   if (model_db->keyExists("Temperature")) {
      temperature_db = model_db->getDatabase("Temperature");
      temperature_type = temperature_db->getString("type");
   } else {
      temperature_db = model_db;
      temperature_type =
          temperature_db->getStringWithDefault("temperature_type", "scalar");
   }

   // fix first letter so we don't need to check for capital letters later
   if (temperature_type[0] == 'S') temperature_type[0] = 's';
   if (temperature_type[0] == 'F') temperature_type[0] = 'f';
   if (temperature_type[0] == 'G') temperature_type[0] = 'g';
   if (temperature_type[0] == 'C') temperature_type[0] = 'c';
   if (temperature_type[0] == 'H') temperature_type[0] = 'h';

   // check if option is valid
   if (temperature_type.compare("scalar") != 0 &&    // scalar
       temperature_type.compare("frozen") != 0 &&    // Frozen (gradient)
       temperature_type.compare("gaussian") != 0 &&  // Gaussian
       temperature_type.compare("constant") != 0 &&  // constant
       temperature_type.compare("heat") != 0 &&      // heat equation
       temperature_type.compare("file") != 0         // read from file
   ) {
      TBOX_ERROR("Error: invalid value for temperature_type");
   }

   d_meltingT =
       temperature_db->getDoubleWithDefault("meltingT", -1.);  // in [K]

   if (temperature_type.compare("scalar") == 0) {

      d_temperature_type = TemperatureType::SCALAR;
      d_with_heat_equation = false;

   } else if (temperature_type.compare("gaussian") == 0) {

      d_temperature_type = TemperatureType::GAUSSIAN;
      d_with_heat_equation = false;

   } else if (temperature_type.compare("frozen") == 0)  // frozen approx.
   {

      d_temperature_type = TemperatureType::GRADIENT;
      d_with_heat_equation = false;

   } else if (temperature_type.compare("heat") == 0) {

      assert(d_molar_volume_liquid == d_molar_volume_liquid);
      assert(d_molar_volume_liquid > 0.);
      assert(d_molar_volume_liquid < 1.e15);

      d_temperature_type = TemperatureType::CONSTANT;
      d_with_heat_equation = true;

      tbox::plog << "Read heat equation parameters..." << std::endl;

      std::string method =
          temperature_db->getStringWithDefault("equation_type", "steady");

      // heat capacity value
      std::shared_ptr<tbox::Database> cp_db = temperature_db->getDatabase("cp");
      std::map<short, double> empty_map;

      d_with_rescaled_temperature = ((method != "steady") && (d_meltingT > 0.));
      if (d_with_rescaled_temperature) {
         d_rescale_factorT = d_meltingT;
         tbox::plog << "Solve temperature equation with rescaled Ti, factor "
                    << d_rescale_factorT << std::endl;
      } else {
         d_rescale_factorT = -1;
      }
      d_cp.push_back(empty_map);
      readSpeciesCP(cp_db->getDatabase("SpeciesA"), d_cp[0]);
      if (d_ncompositions > 0) {
         assert(cp_db->keyExists("SpeciesB"));
         d_cp.push_back(empty_map);
         readSpeciesCP(cp_db->getDatabase("SpeciesB"), d_cp[1]);
      }
      if (d_ncompositions > 1) {
         assert(cp_db->keyExists("SpeciesC"));
         d_cp.push_back(empty_map);
         readSpeciesCP(cp_db->getDatabase("SpeciesC"), d_cp[2]);
      }

      tbox::plog << "Cp for each species: " << std::endl;
      for (std::vector<std::map<short, double> >::iterator it = d_cp.begin();
           it != d_cp.end(); ++it) {
         for (std::map<short, double>::iterator itm = it->begin();
              itm != it->end(); ++itm) {
            itm->second *=
                (1.e-6 / d_molar_volume_liquid);  // conversion from [J/mol*K]
                                                  // to [pJ/(mu m)^3*K]
            tbox::plog << "Cp [pJ/(mu m)^3*K]: " << itm->second << std::endl;
         }
      }
      if (d_with_rescaled_temperature) {
         for (std::vector<std::map<short, double> >::iterator it = d_cp.begin();
              it != d_cp.end(); ++it) {
            for (std::map<short, double>::iterator itm = it->begin();
                 itm != it->end(); ++itm) {
               itm->second *= d_rescale_factorT;
               tbox::plog << "rescaled Cp: " << itm->second << std::endl;
            }
         }
      }

      d_thermal_diffusivity = temperature_db->getDouble("thermal_diffusivity");
      tbox::plog << "Thermal diffusivity [cm^2/s]:        "
                 << d_thermal_diffusivity << std::endl;
      d_thermal_diffusivity *= 1.e8;  // cm^2/s -> um^2/s
      tbox::plog << "Thermal diffusivity [um^2/s]:        "
                 << d_thermal_diffusivity << std::endl;

      if (method == "steady") {
         d_with_steady_temperature = true;

         if (temperature_db->keyExists("heat_source_type")) {
            d_heat_source_type = temperature_db->getString("heat_source_type");
            if (d_heat_source_type == "composition") {
               size_t nterms = temperature_db->getArraySize("source");
               d_T_source.resize(nterms);
               temperature_db->getDoubleArray("source", &d_T_source[0], nterms);
               for (size_t i = 0; i < nterms; ++i)
                  d_T_source[i] *=
                      (1.e-6 /
                       d_molar_volume_liquid);  // conversion from [J/mol] to
                                                // [pJ/(mu m)^3]
            } else {
               assert(d_heat_source_type == "gaussian");
            }
         }
      } else {
         d_with_steady_temperature = false;
      }

      if (d_with_rescaled_temperature) {  // rescale units
         d_H_parameter *= d_rescale_factorT;
      }

   } else if (temperature_type.compare("constant") == 0) {
      d_temperature_type = TemperatureType::CONSTANT;
   } else {
      TBOX_ERROR("ERROR: Temperature type needs to be specified!");
   }

   if (temperature_db->keyExists("latent_heat")) {
      d_latent_heat = temperature_db->getDouble("latent_heat");  // in [J/mol]
      assert(d_molar_volume_liquid == d_molar_volume_liquid);
      // conversion from [J/mol] to [pJ/(mu m)^3]
      d_latent_heat *= (1.e-6 / d_molar_volume_liquid);
      tbox::plog << "Latent heat [pJ/(mu m)^3]:           " << d_latent_heat
                 << std::endl;
      assert(d_latent_heat > 0.);
      assert(d_latent_heat < 1.e32);
   } else {
      d_latent_heat = def_val;
   }
}

//=======================================================================

void QuatModelParameters::initializeOrientation(
    std::shared_ptr<tbox::Database> model_db)
{
   if (model_db->keyExists("noise_amplitude")) {
      d_noise_amplitude = model_db->getDouble("noise_amplitude");
   } else {
      d_noise_amplitude = 0.;
   }

   if (d_H_parameter > 0.) {
      if (model_db->keyExists("orient_mobility")) {
         d_quat_mobility = model_db->getDouble("orient_mobility");
      } else if (model_db->keyExists("quat_mobility")) {
         d_quat_mobility = model_db->getDouble("quat_mobility");
         printDeprecated("quat_mobility", "orient_mobility");
      } else if (model_db->keyExists("tau_quat")) {
         double tau = model_db->getDouble("tau_quat");
         d_quat_mobility = 1. / tau;
         printDeprecated("tau_quat", "orient_mobility");
      } else {
         TBOX_ERROR("Error: quaternion mobility not specified");
      }

      d_min_quat_mobility = 1.e-6;
      if (model_db->keyExists("min_orient_mobility")) {
         d_min_quat_mobility = model_db->getDouble("min_orient_mobility");
      } else if (model_db->keyExists("min_quat_mobility")) {
         d_min_quat_mobility = model_db->getDouble("min_quat_mobility");
         printDeprecated("min_quat_mobility", "min_orient_mobility");
      }

      if (model_db->keyExists("epsilon_orient")) {
         d_epsilon_q = model_db->getDouble("epsilon_orient");
      } else if (model_db->keyExists("epsilon_q")) {
         d_epsilon_q = model_db->getDouble("epsilon_q");
         printDeprecated("epsilon_q", "epsilon_orient");
      } else if (model_db->keyExists("epsilon_quat")) {
         d_epsilon_q = model_db->getDouble("epsilon_quat");
         printDeprecated("epsilon_quat", "epsilon_orient");
      } else {
         TBOX_ERROR("Error: epsilon_quat not specified");
      }

      d_quat_grad_floor =
          model_db->getDoubleWithDefault("orient_grad_floor", 1.e-2);
      if (!model_db->keyExists("orient_grad_floor"))
         if (model_db->keyExists("quat_grad_floor")) {
            d_quat_grad_floor = model_db->getDouble("quat_grad_floor");
            printDeprecated("quat_grad_floor", "orient_grad_floor");
         }

      // options for "smooth floor" are: 1./|grad(q)| =...
      // "max": max( 1./|grad(q)|, 1./d_quat_grad_floor )
      // "tanh": tanh(grad_norm/d_quat_grad_floor)/ grad_norm
      // "sqrt": 1./sqrt( |grad(q)|**2+d_quat_grad_floor**2 )
      d_quat_grad_floor_type =
          model_db->getStringWithDefault("orient_grad_floor_type", "max");
      if (model_db->keyExists("quat_smooth_floor")) {
         printDeprecated("quat_smooth_floor", "orient_grad_floor_type");
      }
      if (model_db->keyExists("orient_smooth_floor")) {
         printDeprecated("quat_smooth_floor", "orient_grad_floor_type");
      }

      d_quat_grad_modulus_type =
          model_db->getStringWithDefault("quat_grad_modulus_type", "cells");

      if (model_db->keyExists("diff_interp_func_type")) {
         d_orient_interp_func_type1 =
             model_db->getString("diff_interp_func_type");
         printDeprecated("diff_interp_func_type", "orient_interp_func_type");
      } else if (model_db->keyExists("orient_interp_func_type")) {
         d_orient_interp_func_type1 =
             model_db->getString("orient_interp_func_type");
         printDeprecated("orient_interp_func_type", "orient_interp_func_type1");
      } else {
         d_orient_interp_func_type1 =
             model_db->getStringWithDefault("orient_interp_func_type1",
                                            "quadratic");
      }
      d_orient_interp_func_type2 =
          model_db->getStringWithDefault("orient_interp_func_type2",
                                         "constan"
                                         "t");
      if (d_orient_interp_func_type1[0] != 'q' &&
          d_orient_interp_func_type1[0] != 'w' &&
          d_orient_interp_func_type1[0] != 'p' &&
          d_orient_interp_func_type1[0] != 'l' &&
          d_orient_interp_func_type1[0] != 't' &&
          d_orient_interp_func_type1[0] != 's' &&
          d_orient_interp_func_type1[0] != '3' &&
          d_orient_interp_func_type1[0] != 'c') {
         TBOX_ERROR("Error: invalid value for orient_interp_func_type1");
      }

      if (model_db->keyExists("quat_mobility_func_type")) {
         d_quat_mobility_func_type =
             model_db->getString("quat_mobility_func_type");
         printDeprecated("quat_mobility_func_type",
                         "orient_mobility_func_"
                         "type");
      } else {
         d_quat_mobility_func_type =
             model_db->getStringWithDefault("orient_mobility_func_type", "pbg");
      }
      if (d_quat_mobility_func_type[0] != 'p' &&
          d_quat_mobility_func_type[0] != 'i' &&
          d_quat_mobility_func_type[0] != 'e') {
         TBOX_ERROR("Error: invalid value for orient_mobility_func_type");
      }

      if (d_quat_mobility_func_type[0] == 'i') {
         d_max_quat_mobility = 1.e6;
         if (model_db->keyExists("max_orient_mobility")) {
            d_max_quat_mobility = model_db->getDouble("max_orient_mobility");
         } else if (model_db->keyExists("max_quat_mobility")) {
            d_max_quat_mobility = model_db->getDouble("max_quat_mobility");
            printDeprecated("max_quat_mobility", "max_orient_mobility");
         }
      }

      if (d_quat_mobility_func_type[0] == 'e' ||
          d_quat_mobility_func_type[0] == 'E') {
         d_exp_scale_quat_mobility = 1.e6;
         if (model_db->keyExists("exp_scale_orient_mobility")) {
            d_exp_scale_quat_mobility =
                model_db->getDouble("exp_scale_orient_mobility");
         } else if (model_db->keyExists("exp_scale_quat_mobility")) {
            d_exp_scale_quat_mobility =
                model_db->getDouble("exp_scale_quat_mobility");
            printDeprecated("exp_scale_quat_mobility",
                            "exp_scale_orient_"
                            "mobility");
         }
      }
   }
}

//=======================================================================

void QuatModelParameters::initializeEta(
    std::shared_ptr<tbox::Database> model_db)
{
   d_epsilon_eta = model_db->getDouble("epsilon_eta");

   d_eta_mobility = model_db->getDouble("eta_mobility");
   d_min_eta_mobility =
       model_db->getDoubleWithDefault("min_eta_mobility", d_eta_mobility);

   d_q0_eta_mobility = model_db->getDoubleWithDefault("q0_eta_mobility", 0.0);

   d_eta_well_scale = model_db->getDoubleWithDefault("eta_well_scale", 1.0);

   d_eta_well_func_type =
       model_db->getStringWithDefault("eta_well_func_type", "double");
   if (d_eta_well_func_type[0] != 's' && d_eta_well_func_type[0] != 'd') {
      TBOX_ERROR("Error: invalid value for eta_well_func_type");
   }

   std::string eta_interp_func_type =
       model_db->getStringWithDefault("eta_interp_func_type", "pbg");
   switch (eta_interp_func_type[0]) {
      case 'l':
      case 'L':
         d_eta_interp_func_type = Thermo4PFM::EnergyInterpolationType::LINEAR;
         break;
      case 'p':
      case 'P':
         d_eta_interp_func_type = Thermo4PFM::EnergyInterpolationType::PBG;
         break;
      case 'h':
      case 'H':
         d_eta_interp_func_type = Thermo4PFM::EnergyInterpolationType::HARMONIC;
         break;
      default:
         tbox::plog << "eta_interp_func_type=" << eta_interp_func_type
                    << std::endl;
         TBOX_ERROR("Error: invalid eta_interp_func_type!!!");
   }
}

//=======================================================================

void QuatModelParameters::readModelParameters(
    std::shared_ptr<tbox::Database> model_db)
{
   d_FolchPlapp = model_db->getBoolWithDefault("three_phases", false);

   if (d_FolchPlapp) {
      d_norderp_A = 1;
      d_norderp_B = 1;
      d_norderp = 3;
   } else {
      // unless otherwise specifies, we use one order parameter
      // phi, and define second phase as 1.-phi
      d_norderp = model_db->getIntegerWithDefault("norderp", 1);
      const int def_norderp_A = (d_norderp == 1) ? 1 : d_norderp - 1;
      d_norderp_A = model_db->getIntegerWithDefault("norderp_A", def_norderp_A);
      if (d_norderp_A > 0) {
         d_norderp_B = d_norderp - d_norderp_A - 1;
         if (d_norderp_B < 0) d_norderp_B = 0;
      }
   }
   tbox::plog << "norderp_A = " << d_norderp_A << std::endl;
   tbox::plog << "norderp_B = " << d_norderp_B << std::endl;
   tbox::plog << "norderp = " << d_norderp << std::endl;

   // Set d_H_parameter to negative value, to turn off orientation terms
   d_H_parameter = model_db->getDoubleWithDefault("H_parameter", -1.);

   // Interface energy
   if (model_db->keyExists("Interface")) {
      std::shared_ptr<tbox::Database> interface_db =
          model_db->getDatabase("Interface");
      if (interface_db->keyExists("sigma")) {
         d_sigma = interface_db->getDouble("sigma");
         if (!interface_db->keyExists("delta"))
            TBOX_ERROR("Interface: sigma and delta  needed together!");
         d_delta = interface_db->getDouble("delta");
         d_epsilon_phase = sqrt(6. * d_sigma * d_delta);
         tbox::plog << "Epsilon_phi = " << d_epsilon_phase << std::endl;
         // factor 16 is AMPE convention
         d_phase_well_scale = (3. * d_sigma / d_delta) / 16.;
         tbox::plog << "Double Well scale = " << d_phase_well_scale
                    << std::endl;
      } else {
         d_epsilon_phase = interface_db->getDouble("epsilon_phi");
         d_phase_well_scale = interface_db->getDouble("phi_well_scale");
      }
      d_with_phase = true;
   } else {  // specify directly epsilon and double well scale

      if (model_db->keyExists("epsilon_phi")) {
         d_epsilon_phase = model_db->getDouble("epsilon_phi");
      } else if (model_db->keyExists("epsilon_phase")) {
         d_epsilon_phase = model_db->getDouble("epsilon_phase");
         printDeprecated("epsilon_phase", "epsilon_phi");
      } else if (model_db->keyExists("epsilon_parameter")) {
         d_epsilon_phase = model_db->getDouble("epsilon_parameter");
         printDeprecated("epsilon_parameter", "epsilon_phi");
      } else {
         d_with_phase = false;
         tbox::pout << "No epsilon specified -> run without phase..."
                    << std::endl;
      }

      if (d_with_phase) {
         if (model_db->keyExists("phi_well_scale")) {
            d_phase_well_scale = model_db->getDouble("phi_well_scale");
         } else if (model_db->keyExists("scale_energy_well")) {
            d_phase_well_scale = model_db->getDouble("scale_energy_well");
            printDeprecated("scale_energy_well", "phi_well_scale");
         }
      }
   }

   d_epsilon_anisotropy =
       model_db->getDoubleWithDefault("epsilon_anisotropy", -1.);
   // we need quaternions to define anisotropy
   if (d_epsilon_anisotropy > 0. && d_H_parameter < 0.) d_H_parameter = 0.;

   d_q0_phase_mobility = model_db->getDoubleWithDefault("q0_phi_mobility", 0.0);

   d_phase_well_func_type = "double";
   if (model_db->keyExists("phi_well_func_type")) {
      d_phase_well_func_type = model_db->getString("phi_well_func_type");
   } else if (model_db->keyExists("energy_well_func_type")) {
      d_phase_well_func_type = model_db->getString("energy_well_func_type");
      printDeprecated("energy_well_func_type", "phi_well_func_type");
   }
   if (d_phase_well_func_type[0] != 's' && d_phase_well_func_type[0] != 'd') {
      TBOX_ERROR("Error: invalid value for phi_well_func_type");
   }

   if (!model_db->keyExists("ConcentrationModel")) {
      if (model_db->keyExists("bias_well_alpha")) {
         d_well_bias_alpha = model_db->getDouble("bias_well_alpha");
         if (d_well_bias_alpha > 0.) {
            d_with_bias_well = true;
         }
      }
      if (model_db->keyExists("bias_well_gamma")) {
         d_well_bias_gamma = model_db->getDouble("bias_well_gamma");
      }
   }


   // Molar volumes
   readMolarVolumes(model_db);

   // Interpolation
   std::string energy_interp_func_type =
       model_db->getStringWithDefault("phi_interp_func_type", "pbg");
   {
      if (model_db->keyExists("energy_interp_func_type")) {
         energy_interp_func_type =
             model_db->getString("energy_interp_func_type");
         printDeprecated("energy_interp_func_type", "phi_interp_func_type");
      }
      switch (energy_interp_func_type[0]) {
         case 'l':
         case 'L':
            d_energy_interp_func_type =
                Thermo4PFM::EnergyInterpolationType::LINEAR;
            break;
         case 'p':
         case 'P':
            d_energy_interp_func_type =
                Thermo4PFM::EnergyInterpolationType::PBG;
            break;
         case 'h':
         case 'H':
            d_energy_interp_func_type =
                Thermo4PFM::EnergyInterpolationType::HARMONIC;
            break;
         case 'u':
         case 'U':
            d_energy_interp_func_type =
                Thermo4PFM::EnergyInterpolationType::UNDEFINED;
            break;
         default:
            tbox::plog << "energy_interp_func_type=" << energy_interp_func_type
                       << std::endl;
            TBOX_ERROR("Error: invalid energy_interp_func_type!!!");
      }
   }

   {
      std::string conc_interp_func_type =
          model_db->getStringWithDefault("conc_interp_func_type",
                                         energy_interp_func_type);
      switch (conc_interp_func_type[0]) {
         case 'l':
         case 'L':
            d_conc_interp_func_type = Thermo4PFM::ConcInterpolationType::LINEAR;
            break;
         case 'p':
         case 'P':
            d_conc_interp_func_type = Thermo4PFM::ConcInterpolationType::PBG;
            break;
         case 'h':
         case 'H':
            d_conc_interp_func_type =
                Thermo4PFM::ConcInterpolationType::HARMONIC;
            break;
         default:
            tbox::plog << "conc_interp_func_type=" << conc_interp_func_type
                       << std::endl;
            TBOX_ERROR("Error: invalid conc_interp_func_type!!!");
      }
   }

   std::string energy_three_args_interp_func_type =
       model_db->getStringWithDefault("energy_three_args_interp_func_type",
                                      "F");
   switch (energy_three_args_interp_func_type[0]) {
      case 'm':
      case 'M':
         d_energy_three_args_interp_func_type =
             EnergyThreeArgsInterpolationType::MOELANS2011;
         break;
      case 'f':
      case 'F':
         d_energy_three_args_interp_func_type =
             EnergyThreeArgsInterpolationType::FOLCHPLAPP2005;
         break;
      default:
         tbox::pout << "energy_three_args_interp_func_type = "
                    << energy_three_args_interp_func_type << std::endl;
         TBOX_ERROR("Error: invalid energy_three_args_interp_func_type!!!");
   }

   std::string conc_three_args_interp_func_type =
       model_db->getStringWithDefault("conc_three_args_interp_func_type", "F");
   switch (conc_three_args_interp_func_type[0]) {
      case 'm':
      case 'M':
         d_conc_three_args_interp_func_type =
             ConcThreeArgsInterpolationType::MOELANS2011;
         break;
      case 'f':
      case 'F':
         d_conc_three_args_interp_func_type =
             ConcThreeArgsInterpolationType::FOLCHPLAPP2005;
         break;
      default:
         tbox::pout << "conc_three_args_interp_func_type = "
                    << conc_three_args_interp_func_type << std::endl;
         TBOX_ERROR("Error: invalid conc_three_args_interp_func_type!!!");
   }

   std::string diffusion_interp_type =
       model_db->getStringWithDefault("diffusion_interp_func_type",
                                      "linea"
                                      "r");
   switch (diffusion_interp_type[0]) {
      case 'l':
      case 'L':
         d_diffusion_interp_type = DiffusionInterpolationType::LINEAR;
         break;
      case 'p':
      case 'P':
         d_diffusion_interp_type = DiffusionInterpolationType::PBG;
         break;
      case 'b':
         d_diffusion_interp_type = DiffusionInterpolationType::BIASED;
         break;
      default:
         tbox::plog << "diffusion_interp_type=" << diffusion_interp_type
                    << std::endl;
         TBOX_ERROR("Error: invalid diffusion_interp_type!!!");
   }

   // Currently "arithmetic" or "harmonic"
   // arithmetic: (x1+x2)/2
   // harmonic:   2/(1/x1+1/x2)
   d_avg_func_type = model_db->getStringWithDefault("avg_func_type",
                                                    "harmoni"
                                                    "c");
   if (d_avg_func_type[0] != 'a' && d_avg_func_type[0] != 'h') {
      TBOX_ERROR("Error: invalid value for avg_func_type");
   }

   d_diffq_avg_func_type =
       model_db->getStringWithDefault("diffq_avg_func_type", d_avg_func_type);
   if (d_diffq_avg_func_type[0] != 'a' && d_diffq_avg_func_type[0] != 'h') {
      TBOX_ERROR("Error: invalid value for avg_func_type");
   }

   d_use_diffs_to_compute_flux =
       model_db->getBoolWithDefault("use_diffs_to_compute_flux", false);
   d_stencil_type = model_db->getStringWithDefault("stencil_type", "normal");

   if (model_db->keyExists("MovingFrame")) {
      std::shared_ptr<tbox::Database> db(model_db->getDatabase("MovingFrame"));

      d_moving_frame_velocity = db->getDoubleWithDefault("velocity", 0.);
      d_moving_frame_adapt = db->getBoolWithDefault("adapt", false);
      d_moving_frame_upwind = db->getBoolWithDefault("upwind", false);
   } else {
      d_moving_frame_velocity = 0.;
      d_moving_frame_adapt = false;
      d_moving_frame_upwind = false;
   }

   if (model_db->keyExists("RigidBody")) {
      assert(d_norderp_A > 0);
      std::shared_ptr<tbox::Database> db(model_db->getDatabase("RigidBody"));
      d_with_rb_motion = true;
      if (db->keyExists("external_force"))
         db->getDoubleArray("external_force", &d_rb_external_force[0], NDIM);
      else
         for (int i = 0; i < NDIM; i++)
            d_rb_external_force[i] = 0.;
      d_rb_mobility = db->getDoubleWithDefault("mobility", 1.);
      d_rb_stiffness = db->getDoubleWithDefault("stiffness", 0.);
      d_rb_threshold = db->getDoubleWithDefault("threshold", 0.);
      if (d_rb_threshold > 0.) d_rb_equil_gb = db->getDouble("equil_gb");
   } else {
      for (int i = 0; i < NDIM; i++)
         d_rb_external_force[i] = def_val;
      d_with_rb_motion = false;
   }

   if (with_third_phase()) initializeEta(model_db);

   if (with_orientation()) initializeOrientation(model_db);

   // we need to know how many species we have before reading temperature
   // model which may include species dependent Cp
   if (model_db->keyExists("ConcentrationModel")) {
      readNumberSpecies(model_db->getDatabase("ConcentrationModel"));
   } else {
      d_ncompositions = 0;
   }

   readTemperatureModel(model_db);

   if (model_db->keyExists("ConcentrationModel")) {
      std::shared_ptr<tbox::Database> conc_db(
          model_db->getDatabase("ConcentrationModel"));
      readConcDB(conc_db);
   }

   // Mobility
   if (model_db->keyExists("PhaseMobility")) {
      assert(d_with_phase);
      std::shared_ptr<tbox::Database> mob_db(
          model_db->getDatabase("PhaseMobility"));
      readPhaseMobility(mob_db);
   } else if (d_with_phase) {
      d_phase_mobility = model_db->getDouble("phi_mobility");
      d_phi_mobility_type = "scalar";
   }
}

void QuatModelParameters::readPhaseMobility(
    std::shared_ptr<tbox::Database> model_db)
{
   d_zetaTref = -1.;
   d_phi_mobility_type = model_db->getStringWithDefault("type", "scalar");
   if (isPhaseMobilityScalar()) {
      if (model_db->keyExists("value")) {
         d_phase_mobility = model_db->getDouble("value");
      } else {
         TBOX_ERROR("Error: phi_mobility not specified");
      }
   } else {
      if (d_phi_mobility_type.compare("Kim") == 0) d_phi_mobility_type = "kim";
      if (d_phi_mobility_type.compare("kim") != 0) {
         tbox::pout << "phi_mobility_type=" << d_phi_mobility_type << std::endl;
         TBOX_ERROR("Error: unknown phi_mobility_type");
      }
      double beta = model_db->getDoubleWithDefault("kinetics_coefficient", -1.);
      if (beta < 0.) {
         d_interface_mobility =
             model_db->getDoubleWithDefault("interface_mobility", -1.);
         if (model_db->keyExists("zeta") || model_db->keyExists("zetaLA")) {
            d_zetaTref = model_db->getDouble("Tref");
            double val[3];
            if (model_db->keyExists("zeta"))
               model_db->getDoubleArray("zeta", &val[0], 3);
            else
               model_db->getDoubleArray("zetaLA", &val[0], 3);
            for (int i = 0; i < 3; i++)
               d_zetaKimMobilityLA[i] = val[i];
            if (model_db->keyExists("zetaLB")) {
               model_db->getDoubleArray("zetaLB", &val[0], 3);
               for (int i = 0; i < 3; i++)
                  d_zetaKimMobilityLB[i] = val[i];
            }
         }
      } else {
         // compute interface_mobility from beta value
         assert(beta > 0.);
         assert(d_meltingT == d_meltingT);
         assert(d_liquidus_slope < 0.);
         d_interface_mobility =
             (d_molar_volume_liquid / (gas_constant_R_JpKpmol * d_meltingT));
         // convert to um^3/pJ
         d_interface_mobility *= 1.e6;

         d_interface_mobility *= (-1. * d_liquidus_slope);
         d_interface_mobility /= (1. - d_keq);
         d_interface_mobility /= beta;
         tbox::plog << "Beta : " << beta << std::endl;
         tbox::plog << "Molar volume : " << d_molar_volume_liquid << std::endl;
         tbox::plog << "Melting T : " << d_meltingT << std::endl;
         tbox::plog << "Liquidus slope : " << d_liquidus_slope << std::endl;
         tbox::plog << "k_eq : " << d_keq << std::endl;
         tbox::plog << "Interface mobility computed : " << d_interface_mobility
                    << std::endl;
         assert(d_interface_mobility == d_interface_mobility);
      }
   }
}

void QuatModelParameters::readFreeEnergies(
    std::shared_ptr<tbox::Database> model_db)
{
   std::shared_ptr<tbox::Database> db(model_db->keyExists("FreeEnergyModel")
                                          ? model_db->getDatabase("FreeEnerg"
                                                                  "y"
                                                                  "Model")
                                          : model_db);

   d_free_energy_type = db->getStringWithDefault("type", "none");

   if (d_free_energy_type[0] == 's') {
      d_free_energy_liquid = db->getDouble("free_energy_liquid");

      if (!d_FolchPlapp) {
         d_free_energy_solid_A = db->getDouble("free_energy_solid");
      } else {
         d_free_energy_solid_A = db->getDouble("free_energy_solid_A");
         d_free_energy_solid_B = db->getDouble("free_energy_solid_B");
      }
   } else {
      if (d_free_energy_type[0] == 'l') {
         d_free_energy_solid_A = 0.;
      }
   }
}

double QuatModelParameters::quatMobilityScaleFactor() const
{
   if (d_quat_mobility_func_type[0] == 'e' ||
       d_quat_mobility_func_type[0] == 'E') {
      return d_exp_scale_quat_mobility;
   }
   if (d_quat_mobility_func_type[0] == 'i' ||
       d_quat_mobility_func_type[0] == 'I') {
      return d_max_quat_mobility;
   }

   return tbox::IEEE::getSignalingNaN();
}
