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
#ifndef included_QuadraticFreeEnergyBinary
#define included_QuadraticFreeEnergyBinary

#include "FreeEnergyStrategyBinary.h"
#include "InterpolationType.h"

#include "QuadraticFreeEnergyFunctionsBinary.h"

#include <string>

class QuadraticFreeEnergyBinary : public FreeEnergyStrategyBinary
{
 public:
   QuadraticFreeEnergyBinary(
       std::shared_ptr<tbox::Database> input_db,
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       const double vml, const double vma, const int conc_l_id,
       const int conc_a_id);

   ~QuadraticFreeEnergyBinary(){};

   void addDrivingForce(const double time, hier::Patch& patch,
                        const int temperature_id, const int phase_id,
                        const int conc_id, const int f_l_id, const int f_a_id,
                        const int f_b_id, const int rhs_id) override;

   void computeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2,
       const bool use_internal_units = true) override;
   void computeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2,
       const bool use_internal_units = true) override;

   void preRunDiagnostics(const double temperature) override{};

 private:
   void addDrivingForceOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_rhs,
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_f_l,
       std::shared_ptr<pdat::CellData<double> > cd_f_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a, const hier::Box& pbox);

   double computeMuL(const double t, const double c);
   double computeMuA(const double t, const double c);

   double computeFreeEnergy(const double temperature, double* c_i,
                            const Thermo4PFM::PhaseIndex pi, const bool gp);

   double computeDerivFreeEnergy(const double temperature, double* c_i,
                                 const Thermo4PFM::PhaseIndex pi);

   double energyFactor(Thermo4PFM::PhaseIndex pi)
   {
      if (pi == Thermo4PFM::PhaseIndex::phaseL)
         return d_energy_conv_factor_L;
      else
         return d_energy_conv_factor_A;
   }

   std::shared_ptr<Thermo4PFM::QuadraticFreeEnergyFunctionsBinary>
       d_quadratic_fenergy;

   Thermo4PFM::EnergyInterpolationType d_energy_interp_func_type;

   double d_vm_L;  // molar volume
   double d_vm_A;  // molar volume

   double d_energy_conv_factor_L;
   double d_energy_conv_factor_A;
};

#endif
