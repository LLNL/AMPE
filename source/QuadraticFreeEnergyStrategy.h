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
#ifndef included_QuadraticFreeEnergyStrategy
#define included_QuadraticFreeEnergyStrategy

#include "FreeEnergyStrategy.h"
#include "InterpolationType.h"

#include <string>

class QuadraticFreeEnergyStrategy : public FreeEnergyStrategy
{
 public:
   QuadraticFreeEnergyStrategy(
       std::shared_ptr<tbox::Database> input_db,
       const EnergyInterpolationType energy_interp_func_type, const double vml,
       const double vma, const double vmb, const int conc_l_id,
       const int conc_a_id, const int conc_b_id, const bool with_third_phase);

   ~QuadraticFreeEnergyStrategy(){};

   void computeFreeEnergyLiquid(hier::Patch& patch, const int temperature_id,
                                const int fl_id, const bool gp = false);

   void computeFreeEnergySolidA(hier::Patch& patch, const int temperature_id,
                                const int fs_id, const bool gp = false);

   void computeFreeEnergySolidB(hier::Patch& patch, const int temperature_id,
                                const int fs_id, const bool gp = false);

   void addDrivingForce(const double time, hier::Patch& patch,
                        const int temperature_id, const int phase_id,
                        const int eta_id, const int conc_id, const int f_l_id,
                        const int f_a_id, const int f_b_id, const int rhs_id);

   void addDrivingForceEta(const double time, hier::Patch& patch,
                           const int temperature_id, const int phase_id,
                           const int eta_id, const int conc_id,
                           const int f_l_id, const int f_a_id, const int f_b_id,
                           const int rhs_id);

   void computeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);
   void computeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);
   void computeSecondDerivativeEnergyPhaseB(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);

   double computeLiquidConcentration(const double hphi, const double heta,
                                     const double c) const;

   double computeSolidAConcentration(const double hphi, const double heta,
                                     const double c) const;

   double computeSolidBConcentration(const double hphi, const double heta,
                                     const double cB) const;

   void preRunDiagnostics(const double temperature){};

 private:
   double computeFreeEnergy(const double temperature, const double conc,
                            const double A, const double Ceq,
                            const double energy_factor,
                            const bool gp = false) const;

   void computeFreeEnergy(hier::Patch& patch, const int temperature_id,
                          const double A, const double Ceq, const int f_id,
                          const int c_i_id, const double energy_factor);

   void computeDerivFreeEnergy(hier::Patch& patch, const int temperature_id,
                               const double A, const double Ceq, const int f_id,
                               const int c_i_id, const double energy_factor);

   double computeLiquidConcentration(const double hphi, const double heta,
                                     const double c, const double Al,
                                     const double Aa, const double Ab,
                                     const double Ceql, const double CeqA,
                                     const double CeqB) const;

   double computeSolidAConcentration(const double hphi, const double heta,
                                     const double c, const double Al,
                                     const double Aa, const double Ab,
                                     const double Ceql, const double CeqA,
                                     const double CeqB) const;

   double computeSolidBConcentration(const double hphi, const double heta,
                                     const double c, const double Al,
                                     const double Aa, const double Ab,
                                     const double Ceql, const double CeqA,
                                     const double CeqB) const;

   void computeFreeEnergy(
       const hier::Box& pbox, std::shared_ptr<pdat::CellData<double> > cd_temp,
       const double A, const double Ceq,
       std::shared_ptr<pdat::CellData<double> > cd_free_energy,
       std::shared_ptr<pdat::CellData<double> > cd_conc_i,
       const double energy_factor);

   void computeDerivFreeEnergy(
       const hier::Box& pbox, std::shared_ptr<pdat::CellData<double> > cd_temp,
       const double A, const double Ceq,
       std::shared_ptr<pdat::CellData<double> > cd_free_energy,
       std::shared_ptr<pdat::CellData<double> > cd_conc_i,
       const double energy_factor);

   void addDrivingForceOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_rhs,
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_eta,
       std::shared_ptr<pdat::CellData<double> > cd_f_l,
       std::shared_ptr<pdat::CellData<double> > cd_f_a,
       std::shared_ptr<pdat::CellData<double> > cd_f_b,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b, const hier::Box& pbox);

   void addDrivingForceEtaOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_rhs,
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_eta,
       std::shared_ptr<pdat::CellData<double> > cd_f_l,
       std::shared_ptr<pdat::CellData<double> > cd_f_a,
       std::shared_ptr<pdat::CellData<double> > cd_f_b,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b, const hier::Box& pbox);

   double computeMu(const double t, const double c);

   double computeDerivFreeEnergy(const double temperature, const double conc,
                                 const double A, const double Ceq,
                                 const double energy_factor) const;

   EnergyInterpolationType d_energy_interp_func_type;

   double d_vm_L;  // molar volume
   double d_vm_A;  // molar volume
   double d_vm_B;  // molar volume
   // double d_vm; // molar volume
   // double d_jpmol2pjpmumcube;

   double d_energy_conv_factor_L;  // molar volume
   double d_energy_conv_factor_A;  // molar volume
   double d_energy_conv_factor_B;  // molar volume

   double d_A_liquid;
   double d_A_solid_A;
   double d_A_solid_B;
   double d_Ceq_liquid;
   double d_Ceq_solid_A;
   double d_Ceq_solid_B;

   bool d_with_third_phase;

   int d_conc_l_id;
   int d_conc_a_id;
   int d_conc_b_id;
};

#endif
