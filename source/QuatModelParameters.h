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
#ifndef included_QuatModelParameters
#define included_QuatModelParameters

#include "tools.h"
#include "CompositionDiffusionStrategy.h"
#include "InterpolationType.h"
#include "ThreeArgsInterpolationType.h"

// Headers for basic SAMRAI objects
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/IEEE.h"

#include <vector>
#include <string>
#include <map>
#include <cmath>

using namespace SAMRAI;

class QuatModelParameters
{
 public:
   QuatModelParameters();

   enum class ConcModel {
      CALPHAD,
      QUADRATIC,
      LINEAR,
      INDEPENDENT,  // energy does not depend on c
      KKSdilute,
      UNDEFINED
   };

   enum class ConcRHSstrategy { KKS, EBS, SPINODAL, Beckermann, UNKNOWN };

   enum class TemperatureType {
      CONSTANT,  // read from a file and left constant in time
      SCALAR,    // constant in space, but possibly not in time
      GAUSSIAN,  // Gaussian in space, varying in time
      GRADIENT,  // linear in space and time
      UNDEFINED
   };

   enum class ConcDiffusionType {
      TD,   // depends on T
      CTD,  // depends on C and T
      UNDEFINED
   };

   void readModelParameters(std::shared_ptr<tbox::Database> quat_db);
   void readPhaseMobility(std::shared_ptr<tbox::Database> model_db);
   void readDiluteAlloy(std::shared_ptr<tbox::Database> conc_db);
   void readTemperatureModel(std::shared_ptr<tbox::Database> model_db);
   void readConcDB(std::shared_ptr<tbox::Database> conc_db);
   void readFreeEnergies(std::shared_ptr<tbox::Database> model_db);
   void readVisitOptions(std::shared_ptr<tbox::Database> visit_db);
   void initializeOrientation(std::shared_ptr<tbox::Database> quat_db);
   void initializeEta(std::shared_ptr<tbox::Database> quat_db);

   // accessors for parameters
   int norderp() const { return (with_three_phases() ? 3 : d_norderp); }
   int norderpA() const { return d_norderp_A; }
   int norderpB() const { return d_norderp_B; }

   bool withPhaseB() const { return (with_three_phases() || d_norderp_B > 0); }
   bool withMultipleOrderP() const { return (d_norderp_A > -1); }

   double H_parameter() const { return d_H_parameter; }
   bool with_orientation() const { return (d_H_parameter >= 0.); }
   bool evolveQuat() const
   {
      assert(d_H_parameter == d_H_parameter);
      return (d_H_parameter > 0.);
   }
   double epsilon_phase() const { return d_epsilon_phase; }
   double epsilon_anisotropy() const { return d_epsilon_anisotropy; }
   double epsilon_eta() const { return d_epsilon_eta; }
   double epsilon_q() const { return d_epsilon_q; }
   double noise_amplitude() const { return d_noise_amplitude; }
   double phase_mobility() const
   {
      assert(d_phase_mobility > 0.);
      return d_phase_mobility;
   }
   double eta_mobility() const { return d_eta_mobility; }
   double quat_mobility() const { return d_quat_mobility; }
   double min_eta_mobility() const { return d_min_eta_mobility; }
   double min_quat_mobility() const { return d_min_quat_mobility; }
   double max_quat_mobility() const { return d_max_quat_mobility; }
   double exp_scale_quat_mobility() const { return d_exp_scale_quat_mobility; }
   double q0_phase_mobility() const { return d_q0_phase_mobility; }
   double q0_eta_mobility() const { return d_q0_eta_mobility; }
   double quat_grad_floor() const { return d_quat_grad_floor; }
   double conc_mobility() const { return d_conc_mobility; }
   double thermal_diffusivity() const { return d_thermal_diffusivity; }
   double latent_heat() const { return d_latent_heat; }
   const std::map<short, double>& cp(const unsigned short isp) const
   {
      assert(isp < d_cp.size());
      return d_cp[isp];
   }
   const std::vector<std::map<short, double> >& cp() const { return d_cp; }

   int ncompositions() const
   {
      assert(d_ncompositions >= 0);
      return d_ncompositions;
   }
   bool knownInitCinPhase() const { return !d_initc_in_phase.empty(); }
   double getInitCphaseL(const int index) const
   {
      return d_initc_in_phase[index];
   }
   double getInitCphaseA(const int index) const
   {
      return d_initc_in_phase[d_ncompositions + index];
   }
   double meltingT() const { return d_meltingT; }
   double interfaceMobility() const { return d_interface_mobility; }
   double rescale_factorT() const { return d_rescale_factorT; }
   double vd() const
   {
      assert(d_vd > 0.);
      return d_vd;
   }
   double keq() const { return d_keq; }
   double liquidus_slope() const { return d_liquidus_slope; }
   double average_concentration() const { return d_average_concentration; }


   const std::vector<double>& T_source() const { return d_T_source; }
   const std::vector<double>& Q_heat_transport() const
   {
      return d_Q_heat_transport;
   }

   double phase_well_scale() const { return d_phase_well_scale; }
   double eta_well_scale() const { return d_eta_well_scale; }
   double free_energy_liquid() const { return d_free_energy_liquid; }
   double free_energy_solid_A() const { return d_free_energy_solid_A; }
   double free_energy_solid_B() const { return d_free_energy_solid_B; }
   std::string free_energy_type() const { return d_free_energy_type; }

   double well_bias_alpha() const { return d_well_bias_alpha; }
   double well_bias_gamma() const { return d_well_bias_gamma; }

   std::string orient_interp_func_type1() const
   {
      return d_orient_interp_func_type1;
   }
   std::string orient_interp_func_type2() const
   {
      return d_orient_interp_func_type2;
   }
   Thermo4PFM::ConcInterpolationType conc_interp_func_type() const
   {
      return d_conc_interp_func_type;
   }
   Thermo4PFM::EnergyInterpolationType energy_interp_func_type() const
   {
      assert(d_energy_interp_func_type !=
             Thermo4PFM::EnergyInterpolationType::UNDEFINED);
      return d_energy_interp_func_type;
   }
   Thermo4PFM::EnergyInterpolationType eta_interp_func_type() const
   {
      return d_eta_interp_func_type;
   }
   DiffusionInterpolationType diffusion_interp_func_type() const
   {
      return d_diffusion_interp_type;
   }
   ConcThreeArgsInterpolationType conc_three_args_interp_func_type() const
   {
      return d_conc_three_args_interp_func_type;
   }
   EnergyThreeArgsInterpolationType energy_three_args_interp_func_type() const
   {
      return d_energy_three_args_interp_func_type;
   }

   std::string avg_func_type() const { return d_avg_func_type; }
   std::string diffq_avg_func_type() const { return d_diffq_avg_func_type; }
   std::string phase_well_func_type() const { return d_phase_well_func_type; }
   std::string eta_well_func_type() const { return d_eta_well_func_type; }
   std::string quat_mobility_func_type() const
   {
      return d_quat_mobility_func_type;
   }
   std::string quat_grad_floor_type() const { return d_quat_grad_floor_type; }

   bool quat_grad_modulus_from_cells() const
   {
      return (d_quat_grad_modulus_type.compare("cells") == 0);
   }

   bool isPhaseMobilityScalar() const
   {
      return (d_phi_mobility_type.compare("scalar") == 0);
   }

   double gamma() const { return d_gamma; }
   double m_moelans2011() const { return 8. * d_phase_well_scale; }

   double molar_volume_liquid() const { return d_molar_volume_liquid; }
   double molar_volume_solid_A() const { return d_molar_volume_solid_A; }
   double molar_volume_solid_B() const { return d_molar_volume_solid_B; }
   double D_liquid() const { return d_D_liquid; }
   double D_solid_A() const { return d_D_solid_A; }
   double D_solid_B() const { return d_D_solid_B; }
   double Q0_liquid() const { return d_Q0_liquid; }
   double Q0_solid_A() const { return d_Q0_solid_A; }
   double Q0_solid_B() const { return d_Q0_solid_B; }
   std::string conc_avg_func_type() const { return d_conc_avg_func_type; }

   bool with_phase() const { return d_with_phase; }
   bool with_concentration() const { return d_with_concentration; }
   bool with_third_phase() const { return d_with_third_phase; }
   bool with_three_phases() const { return d_with_three_phases; }
   bool with_heat_equation() const { return d_with_heat_equation; }
   bool with_unsteady_heat_equation() const
   {
      return d_with_heat_equation && !d_with_steady_temperature;
   }
   bool with_steady_temperature() const { return d_with_steady_temperature; }
   bool with_rescaled_temperature() const
   {
      return d_with_rescaled_temperature;
   }
   bool with_gradT() const { return d_with_gradT; }
   bool with_antitrapping() const { return d_with_antitrapping; }
   bool grand_potential() const { return d_grand_potential; }
   bool with_bias_well() const { return d_with_bias_well; }
   bool with_Aziz_partition_coeff() const
   {
      return (d_partition_coeff.compare("Aziz") == 0);
   }
   bool with_uniform_partition_coeff() const
   {
      return (d_partition_coeff.compare("uniform") == 0);
   }
   bool with_partition_coeff() const
   {
      return (d_partition_coeff.compare("uniform") == 0 ||
              d_partition_coeff.compare("Aziz") == 0);
   }
   bool partition_phase_concentration() const
   {
      return (d_phase_concentration_model.compare("partition") == 0);
   }
   bool kks_phase_concentration() const
   {
      return (d_phase_concentration_model.compare("kks") == 0);
   }
   bool with_velocity() const { return d_with_velocity; }
   bool use_diffs_to_compute_flux() const
   {
      return d_use_diffs_to_compute_flux;
   }
   bool useIsotropicStencil() const
   {
      return d_stencil_type.compare("isotropic") == 0;
   }
   bool useEBSisotropicStencil() const
   {
      return d_ebs_stencil_type.compare("isotropic") == 0;
   }
   bool useEBS4thOrderStencil() const
   {
      return d_ebs_stencil_type.compare("order4") == 0;
   }
   bool useUpwindScheme() const { return d_moving_frame_upwind; }
   int nghosts_required() const
   {
      if (useEBS4thOrderStencil() || useUpwindScheme())
         return 2;
      else
         return 1;
   }
   bool wellBiasBeckermann() const { return d_bias_well_beckermann; }

   double quatMobilityScaleFactor() const;

   bool isConcentrationModelLinear() const
   {
      return (d_conc_model == ConcModel::LINEAR);
   }

   bool isConcentrationModelCALPHADorQuadratic() const
   {
      return (d_conc_model == ConcModel::CALPHAD ||
              d_conc_model == ConcModel::QUADRATIC);
   }

   bool concentrationModelNeedsPhaseConcentrations() const
   {
      return (d_conc_model == ConcModel::CALPHAD ||
              d_conc_model == ConcModel::QUADRATIC ||
              d_conc_model == ConcModel::LINEAR ||
              d_conc_model == ConcModel::KKSdilute ||
              (d_conc_model == ConcModel::INDEPENDENT && d_with_concentration));
   }

   bool isConcentrationModelCALPHAD() const
   {
      assert(d_conc_model != ConcModel::UNDEFINED);

      return (d_conc_model == ConcModel::CALPHAD);
   }

   bool isConcentrationModelQuadratic() const
   {
      assert(d_conc_model != ConcModel::UNDEFINED);

      return (d_conc_model == ConcModel::QUADRATIC);
   }

   bool isConcentrationModelKKSdilute() const
   {
      assert(d_conc_model != ConcModel::UNDEFINED);

      return (d_conc_model == ConcModel::KKSdilute);
   }

   bool isTemperatureUniform() const
   {
      return (d_temperature_type == TemperatureType::SCALAR);
   }
   bool isTemperatureConstant() const
   {
      return (d_temperature_type == TemperatureType::CONSTANT);
   }
   bool isTemperatureGaussian() const
   {
      return (d_temperature_type == TemperatureType::GAUSSIAN);
   }
   bool isTemperatureGradient() const
   {
      return (d_temperature_type == TemperatureType::GRADIENT);
   }

   void checkValidityConcRHSstrategy() const
   {
      assert(d_conc_rhs_strategy == ConcRHSstrategy::KKS ||
             d_conc_rhs_strategy == ConcRHSstrategy::EBS ||
             d_conc_rhs_strategy == ConcRHSstrategy::SPINODAL ||
             d_conc_rhs_strategy == ConcRHSstrategy::Beckermann);
   }

   bool needGhosts4PartitionCoeff() const
   {
      return (d_conc_rhs_strategy == ConcRHSstrategy::Beckermann);
   }

   bool concRHSstrategyIsKKS() const
   {
      return (d_conc_rhs_strategy == ConcRHSstrategy::KKS);
   }
   bool concRHSstrategyIsEBS() const
   {
      return (d_conc_rhs_strategy == ConcRHSstrategy::EBS);
   }
   bool concRHSstrategyIsSPINODAL() const
   {
      return (d_conc_rhs_strategy == ConcRHSstrategy::SPINODAL);
   }
   bool concRHSstrategyIsBeckermann() const
   {
      return (d_conc_rhs_strategy == ConcRHSstrategy::Beckermann);
   }

   bool conDiffusionStrategyIsCTD() const
   {
      return (d_conc_diffusion_type == ConcDiffusionType::CTD);
   }

   bool isHeatSourceCompositionDependent() const
   {
      return d_heat_source_type == "composition";
   }

   bool with_extra_visit_output() const { return d_extra_visit_output; }

   bool with_rhs_visit_output() const { return d_rhs_visit_output; }

   bool with_visit_energy_output() const { return d_visit_energy_output; }

   bool with_visit_grain_output() const { return d_visit_grain_output; }

   unsigned ncompositionFields() const
   {
      assert(d_ncompositions >= 0);
      return d_ncompositions;
   }

   double surfaceEnergy() const
   {
      return d_epsilon_phase * sqrt(16. * d_phase_well_scale) / (3. * sqrt(2.));
   }

   /*
    * Interfacial width based on Boettinget at al. formula
    */
   double interfacialWidth() const
   {
      return d_epsilon_phase / sqrt(32. * d_phase_well_scale);
   }

   bool initPhaseConcAtEq() const { return d_init_phase_conc_eq; }

   bool inMovingFrame() const
   {
      return (std::abs(d_moving_frame_velocity) > 1.e-16);
   }

   double movingVelocity() const { return d_moving_frame_velocity; }

   bool adaptMovingFrame() const { return d_moving_frame_adapt; }

   bool upwindMovingFrame() const { return d_moving_frame_upwind; }

   double zetaFactorLA(const double temp) const
   {
      if (d_zetaTref < 0.) return -1.;

      const double dT = (temp - d_zetaTref);
      return d_zetaKimMobilityLA[0] + d_zetaKimMobilityLA[1] * dT +
             d_zetaKimMobilityLA[2] * dT * dT;
   }

   double zetaFactorLB(const double temp) const
   {
      if (d_zetaTref < 0.) return -1.;

      const double dT = (temp - d_zetaTref);
      return d_zetaKimMobilityLB[0] + d_zetaKimMobilityLB[1] * dT +
             d_zetaKimMobilityLB[2] * dT * dT;
   }

   double zetaFactor(const double temp) const { return zetaFactorLA(temp); }

   double concL_ref() { return d_concL_ref; }
   double concA_ref() { return d_concA_ref; }
   double concB_ref() { return d_concB_ref; }


 private:
   void readNumberSpecies(std::shared_ptr<tbox::Database> conc_db);

   // total number of order parameters
   int d_norderp;

   // number of order parameters associated with phase A
   int d_norderp_A;
   // number of order parameters associated with phase B
   int d_norderp_B;

   // Model parameters
   double d_H_parameter;
   double d_epsilon_phase;
   double d_epsilon_anisotropy;
   double d_epsilon_eta;
   double d_epsilon_q;
   double d_noise_amplitude;
   double d_phase_mobility;
   double d_eta_mobility;
   double d_quat_mobility;
   double d_min_eta_mobility;
   double d_min_quat_mobility;
   double d_max_quat_mobility;
   double d_exp_scale_quat_mobility;
   double d_q0_phase_mobility;
   double d_q0_eta_mobility;
   double d_quat_grad_floor;
   double d_conc_mobility;
   double d_thermal_diffusivity;
   double d_latent_heat;
   // cp for each species
   std::vector<std::map<short, double> > d_cp;
   double d_meltingT;
   double d_interface_mobility;

   // inverse zeta polynomial expansion coefficients in Kim's mobility
   double d_zetaTref;
   std::array<double, 3> d_zetaKimMobilityLA;
   std::array<double, 3> d_zetaKimMobilityLB;

   std::string d_heat_source_type;
   std::vector<double> d_T_source;
   std::vector<double> d_Q_heat_transport;

   /*
    * Initial compositions in each phase:
    * cL0, cL1, ..., cS0, cS1, ...
    */
   std::vector<double> d_initc_in_phase;

   // free energy parameters:
   // f(phi) = d_phase_well_scale * g(phi)
   //        + p(phi)*( d_free_energy_solid - d_free_energy_liquid)
   // where g is a well potential and p an interpolation function
   // s.t. p(0)=0 and p(1)=1
   double d_phase_well_scale;
   double d_eta_well_scale;
   double d_free_energy_liquid;
   double d_free_energy_solid_A;
   double d_free_energy_solid_B;
   std::string d_free_energy_type;

   double d_well_bias_alpha;
   double d_well_bias_gamma;
   double d_liquidus_slope;
   double d_average_concentration;

   // multiple order model parameters
   double d_gamma;

   /*!
    * polynomial types for phi interpolation in quaternion diffusion
    * coefficient (1st and 2nd order grad q terms)
    * valid options are:
    *    "q" for quadratic, phi^2 (default for |grad q|)
    *    "3" for cubic, phi^3
    *    "p" for "pbg"
    *    "c" for constant, h(phi)=1 (default for |grad q|^2)
    */
   std::string d_orient_interp_func_type1;
   std::string d_orient_interp_func_type2;

   /*!
    * form of h_r(phi)
    */
   Thermo4PFM::ConcInterpolationType d_conc_interp_func_type;

   /*!
    * form of h_p(phi)
    */
   Thermo4PFM::EnergyInterpolationType d_energy_interp_func_type;

   Thermo4PFM::EnergyInterpolationType d_eta_interp_func_type;

   /*!
    * form of h_d(phi)
    */
   DiffusionInterpolationType d_diffusion_interp_type;

   EnergyThreeArgsInterpolationType d_energy_three_args_interp_func_type;
   ConcThreeArgsInterpolationType d_conc_three_args_interp_func_type;

   std::string d_avg_func_type;
   std::string d_diffq_avg_func_type;
   std::string d_phase_well_func_type;
   std::string d_eta_well_func_type;
   std::string d_quat_mobility_func_type;
   std::string d_quat_grad_floor_type;

   /*
    * option what grad quat data to use to compute modulus of quaternion
    * gradient possible options: "cells" (default) or "sides"
    */
   std::string d_quat_grad_modulus_type;

   std::string d_phi_mobility_type;

   double d_molar_volume_liquid;
   double d_molar_volume_solid_A;
   double d_molar_volume_solid_B;
   double d_D_liquid;
   double d_D_solid_A;
   double d_D_solid_B;
   double d_Q0_liquid;
   double d_Q0_solid_A;
   double d_Q0_solid_B;
   std::string d_conc_avg_func_type;

   // reference auxilliary compositions to use for initialization
   double d_concL_ref;
   double d_concA_ref;
   double d_concB_ref;

   ConcModel d_conc_model;
   ConcRHSstrategy d_conc_rhs_strategy;

   TemperatureType d_temperature_type;
   ConcDiffusionType d_conc_diffusion_type;

   bool d_with_phase;
   bool d_with_concentration;
   bool d_with_third_phase;
   bool d_with_three_phases;
   bool d_with_heat_equation;
   bool d_with_steady_temperature;
   bool d_with_gradT;
   bool d_with_antitrapping;
   bool d_grand_potential;
   bool d_with_bias_well;
   bool d_bias_well_beckermann;
   bool d_with_rescaled_temperature;

   double d_rescale_factorT;

   int d_ncompositions;

   std::string d_partition_coeff;
   // diffusion speed corresponding to interface used for partition coefficient
   double d_vd;
   double d_keq;
   std::string d_phase_concentration_model;

   bool d_with_velocity;

   /*!
    * Specify use of grad q at sides to compute fluxes in dq/dt equation
    */
   bool d_use_diffs_to_compute_flux;

   /*!
    * Specify stencil type to use for grad q at sides ("normal" or "isotropic")
    */
   std::string d_stencil_type;
   std::string d_ebs_stencil_type;

   bool d_extra_visit_output;
   bool d_rhs_visit_output;
   bool d_visit_energy_output;
   bool d_visit_grain_output;

   bool d_init_phase_conc_eq;

   double d_moving_frame_velocity;
   bool d_moving_frame_adapt;
   bool d_moving_frame_upwind;

   void readMolarVolumes(std::shared_ptr<tbox::Database> db);
};

#endif
