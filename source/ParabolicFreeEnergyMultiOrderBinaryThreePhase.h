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
#ifndef included_ParabolicFreeEnergyMultiOrderBinaryThreePhase
#define included_ParabolicFreeEnergyMultiOrderBinaryThreePhase

#include "ParabolicFreeEnergyFunctionsBinaryThreePhase.h"
#include "FreeEnergyStrategyBinary.h"
#include "MolarVolumeStrategy.h"
#include "MultiOrderBinaryThreePhasesDrivingForce.h"

class ParabolicFreeEnergyMultiOrderBinaryThreePhase
    : public FreeEnergyStrategyBinary
{
 public:
   ParabolicFreeEnergyMultiOrderBinaryThreePhase(
       std::shared_ptr<tbox::Database> input_db,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       const short norderp_A, MolarVolumeStrategy* mvstrategy,
       const int conc_l_id, const int conc_a_id, const int conc_b_id);

   ~ParabolicFreeEnergyMultiOrderBinaryThreePhase();

   void addDrivingForce(const double time, hier::Patch& patch,
                        const int temperature_id, const int phase_id,
                        const int conc_id, const int f_l_id, const int f_a_id,
                        const int f_b_id, const int rhs_id) override;

   void preRunDiagnostics(const double temperature){};

   bool computeCeqT(const double temperature, const Thermo4PFM::PhaseIndex pi0,
                    const Thermo4PFM::PhaseIndex pi1, double* ceq);

 protected:
   double computeFreeEnergy(const double temperature, double* c_i,
                            const Thermo4PFM::PhaseIndex pi,
                            const bool gp) override;

   double computeDerivFreeEnergy(const double temperature, double* c_i,
                                 const Thermo4PFM::PhaseIndex pi) override;

 private:
   MolarVolumeStrategy* d_mv_strategy;

   std::shared_ptr<Thermo4PFM::ParabolicFreeEnergyFunctionsBinaryThreePhase>
       d_parabolic_fenergy;

   std::shared_ptr<MultiOrderBinaryThreePhasesDrivingForce>
       d_multiorder_driving_force;

   void computeSecondDerivativeEnergyPhaseL(
       const double temp, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units) override;
   void computeSecondDerivativeEnergyPhaseA(
       const double temp, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units) override;
   void computeSecondDerivativeEnergyPhaseB(const double temp,
                                            const std::vector<double>& c,
                                            std::vector<double>& d2fdc2,
                                            const bool use_internal_units);

   double computeMuL(const double t, const double c);
   double computeMuA(const double t, const double c);
   double computeMuB(const double t, const double c);
};

#endif
