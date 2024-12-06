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
#ifndef included_QuadraticFreeEnergyMultiOrderTernaryThreePhase
#define included_QuadraticFreeEnergyMultiOrderTernaryThreePhase

#include "QuadraticFreeEnergyFunctionsTernaryThreePhase.h"
#include "FreeEnergyStrategyThreePhase.h"
#include "InterpolationType.h"

class QuadraticFreeEnergyMultiOrderTernaryThreePhase
    : public FreeEnergyStrategyThreePhase
{
 public:
   QuadraticFreeEnergyMultiOrderTernaryThreePhase(
       std::shared_ptr<tbox::Database> input_db,
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       const short norderp_A, const double vml, const double vma,
       const double vmb, const int conc_l_id, const int conc_a_id,
       const int conc_b_id);

   ~QuadraticFreeEnergyMultiOrderTernaryThreePhase();

 private:
   //
   // number of order parameters associated with phase A
   //
   const short d_norderp_A;

   std::shared_ptr<Thermo4PFM::QuadraticFreeEnergyFunctionsTernaryThreePhase>
       d_quadratic_fenergy;

   void computeSecondDerivativeEnergyPhaseL(
       const std::vector<double>& c, std::vector<double>& d2fdc2,
       const bool use_internal_units) override;
   void computeSecondDerivativeEnergyPhaseA(
       const std::vector<double>& c, std::vector<double>& d2fdc2,
       const bool use_internal_units) override;
   void computeSecondDerivativeEnergyPhaseB(
       const std::vector<double>& c, std::vector<double>& d2fdc2,
       const bool use_internal_units) override;

   void computeMuL(const double t, const double c0, const double c1,
                   double* mu);
   void computeMuA(const double t, const double c0, const double c1,
                   double* mu);
   void computeMuB(const double t, const double c0, const double c1,
                   double* mu);

   void addDrivingForceOnPatch(
       std::shared_ptr<pdat::CellData<double> > cd_rhs,
       std::shared_ptr<pdat::CellData<double> > cd_temperature,
       std::shared_ptr<pdat::CellData<double> > cd_phi,
       std::shared_ptr<pdat::CellData<double> > cd_f_l,
       std::shared_ptr<pdat::CellData<double> > cd_f_a,
       std::shared_ptr<pdat::CellData<double> > cd_f_b,
       std::shared_ptr<pdat::CellData<double> > cd_c_l,
       std::shared_ptr<pdat::CellData<double> > cd_c_a,
       std::shared_ptr<pdat::CellData<double> > cd_c_b,
       const hier::Box& pbox) override;

   void computeFreeEnergy(
       const hier::Box& pbox, std::shared_ptr<pdat::CellData<double> > cd_temp,
       std::shared_ptr<pdat::CellData<double> > cd_free_energy,
       std::shared_ptr<pdat::CellData<double> > cd_conc_i,
       Thermo4PFM::PhaseIndex pi, const double energy_factor) override;
};

#endif
