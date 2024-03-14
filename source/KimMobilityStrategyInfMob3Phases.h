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
       const QuatModelParameters& parameters, QuatModel* quat_model,
       const int conc_l_id, const int conc_a_id, const int conc_b_id,
       const int temp_id, const double epsilon, const double phase_well_scale,
       const EnergyThreeArgsInterpolationType energy_interp_func_type,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       std::shared_ptr<tbox::Database> conc_db, const unsigned ncompositions,
       const double mv);

 private:
   double evaluateMobility(const double temp,
                           const std::vector<double>& phaseconc,
                           const std::vector<double>& phi);

   const QuatModelParameters& d_model_parameters;

   // molar volume
   const double d_mv;
   double d_factor;
};

#endif
