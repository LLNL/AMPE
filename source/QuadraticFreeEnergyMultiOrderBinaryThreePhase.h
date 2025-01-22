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
#ifndef included_QuadraticFreeEnergyMultiOrderBinaryThreePhase
#define included_QuadraticFreeEnergyMultiOrderBinaryThreePhase

#include "QuadraticFreeEnergyFunctionsBinaryThreePhase.h"
#include "MultiOrderBinaryThreePhasesDrivingForce.h"
#include "FreeEnergyStrategyBinary.h"

class QuadraticFreeEnergyMultiOrderBinaryThreePhase
    : public FreeEnergyStrategyBinary
{
 public:
   QuadraticFreeEnergyMultiOrderBinaryThreePhase(
       std::shared_ptr<tbox::Database> input_db,
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       const short norderp_A, const double vml, const double vma,
       const double vmb, const int conc_l_id, const int conc_a_id,
       const int conc_b_id);

   ~QuadraticFreeEnergyMultiOrderBinaryThreePhase();

 private:
   //
   // number of order parameters associated with phase A
   //
   const short d_norderp_A;

   double d_energy_conv_factor_L;
   double d_energy_conv_factor_A;
   double d_energy_conv_factor_B;

   std::shared_ptr<Thermo4PFM::QuadraticFreeEnergyFunctionsBinaryThreePhase>
       d_quadratic_fenergy;

   std::shared_ptr<MultiOrderBinaryThreePhasesDrivingForce> d_driving_force;

   void computeSecondDerivativeEnergyPhaseL(
       const double temp, const std::vector<double>& c,
       std::vector<double>& d2fdc2,
       const bool use_internal_units = true) override;
   void computeSecondDerivativeEnergyPhaseA(
       const double temp, const std::vector<double>& c,
       std::vector<double>& d2fdc2,
       const bool use_internal_units = true) override;
   void computeSecondDerivativeEnergyPhaseB(
       const double temp, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);

   double computeMuL(const double t, const double c0);
   double computeMuA(const double t, const double c0);

   double computeFreeEnergy(const double temperature, double* c_i,
                            const Thermo4PFM::PhaseIndex pi,
                            const bool gp) override;

   double computeDerivFreeEnergy(const double temperature, double* c_i,
                                 const Thermo4PFM::PhaseIndex pi) override;

   double energyFactor(Thermo4PFM::PhaseIndex pi)
   {
      if (pi == Thermo4PFM::PhaseIndex::phaseL)
         return d_energy_conv_factor_L;
      else
         return d_energy_conv_factor_A;
   }

   void addDrivingForce(const double time, hier::Patch& patch,
                        const int temperature_id, const int phase_id,
                        const int conc_id, const int f_l_id, const int f_a_id,
                        const int f_b_id, const int rhs_id) override;
};

#endif
