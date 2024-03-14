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

#include "QuadraticFreeEnergyFunctionsBinary.h"

#include <string>

class QuadraticFreeEnergyStrategy : public FreeEnergyStrategy
{
 public:
   QuadraticFreeEnergyStrategy(
       std::shared_ptr<tbox::Database> input_db,
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       const double vml, const double vma, const int conc_l_id,
       const int conc_a_id);

   ~QuadraticFreeEnergyStrategy(){};

   void computeFreeEnergyLiquid(hier::Patch& patch, const int temperature_id,
                                const int fl_id,
                                const bool gp = false) override;

   void computeFreeEnergySolidA(hier::Patch& patch, const int temperature_id,
                                const int fs_id,
                                const bool gp = false) override;

   void computeFreeEnergySolidB(hier::Patch& patch, const int temperature_id,
                                const int fs_id,
                                const bool gp = false) override;

   void addDrivingForce(const double time, hier::Patch& patch,
                        const int temperature_id, const int phase_id,
                        const int eta_id, const int conc_id, const int f_l_id,
                        const int f_a_id, const int f_b_id,
                        const int rhs_id) override;

   void computeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2,
       const bool use_internal_units = true) override;
   void computeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2,
       const bool use_internal_units = true) override;

   double computeLiquidConcentration(const double temp, const double hphi,
                                     const double c) const;

   double computeSolidAConcentration(const double temp, const double hphi,
                                     const double c) const;

   void preRunDiagnostics(const double temperature) override{};

 private:
   double computeFreeEnergy(const double conc, const double A, const double Ceq,
                            const double energy_factor,
                            const bool gp = false) const;
   double computeDerivFreeEnergy(const double conc, const double A,
                                 const double Ceq,
                                 const double energy_factor) const;

   void computeFreeEnergy(hier::Patch& patch, const int temperature_id,
                          const int f_id, const int c_i_id,
                          Thermo4PFM::PhaseIndex pi,
                          const double energy_factor);

   void computeDerivFreeEnergy(hier::Patch& patch, const int temperature_id,
                               const int f_id, const int c_i_id,
                               Thermo4PFM::PhaseIndex pi,
                               const double energy_factor);

   void computeFreeEnergy(
       const hier::Box& pbox, std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::CellData<double> > cd_free_energy,
       std::shared_ptr<pdat::CellData<double> > cd_conc_i,
       Thermo4PFM::PhaseIndex pi, const double energy_factor);

   void computeDerivFreeEnergy(
       const hier::Box& pbox, std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::CellData<double> > cd_free_energy,
       std::shared_ptr<pdat::CellData<double> > cd_conc_i,
       Thermo4PFM::PhaseIndex pi, const double energy_factor);

   void addDrivingForceOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_rhs,
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_f_l,
       std::shared_ptr<pdat::CellData<double> > cd_f_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a, const hier::Box& pbox);

   double computeMu(const double t, const double c);

   std::shared_ptr<Thermo4PFM::QuadraticFreeEnergyFunctionsBinary>
       d_quadratic_fenergy;

   Thermo4PFM::EnergyInterpolationType d_energy_interp_func_type;

   double d_vm_L;  // molar volume
   double d_vm_A;  // molar volume
   // double d_vm; // molar volume
   // double d_jpmol2pjpmumcube;

   double d_energy_conv_factor_L;  // molar volume
   double d_energy_conv_factor_A;  // molar volume

   int d_conc_l_id;
   int d_conc_a_id;
};

#endif
