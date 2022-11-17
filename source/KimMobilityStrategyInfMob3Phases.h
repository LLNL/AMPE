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
#ifdef HAVE_THERMO4PFM

#ifndef included_KimMobilityStrategyInfMob3Phases
#define included_KimMobilityStrategyInfMob3Phases

#include "KimMobilityStrategy.h"
#include "CALPHADFreeEnergyFunctionsBinary2Ph1Sl.h"

template <class FreeEnergyType>
class KimMobilityStrategyInfMob3Phases
    : public KimMobilityStrategy<FreeEnergyType>
{
 public:
   KimMobilityStrategyInfMob3Phases(
       QuatModel* quat_model, const int conc_l_id, const int conc_a_id,
       const int conc_b_id, const int temp_id, const double epsilon,
       const double phase_well_scale,
       const EnergyThreeArgsInterpolationType energy_interp_func_type,
       const ConcInterpolationType conc_interp_func_type,
       std::shared_ptr<tbox::Database> conc_db, const unsigned ncompositions,
       const double DL, const double Q0, const double mv);

 private:
   double compute_zeta(const double* const cl, const double* const cs,
                       const double temp);

   double evaluateMobility(const double temp,
                           const std::vector<double>& phaseconc,
                           const std::vector<double>& phi);

   double d_DL;
   double d_Q0;

   std::vector<double> d_d2fdc2;

   double d_factor;

   std::shared_ptr<CALPHADFreeEnergyFunctionsBinary2Ph1Sl> d_free_energy_LA;
   std::shared_ptr<CALPHADFreeEnergyFunctionsBinary2Ph1Sl> d_free_energy_LB;
};

#endif
#endif
