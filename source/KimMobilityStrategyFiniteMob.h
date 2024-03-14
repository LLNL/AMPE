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
#ifndef included_KimMobilityStrategyFiniteMob
#define included_KimMobilityStrategyFiniteMob

#include "KimMobilityStrategy.h"

template <class FreeEnergyType>
class KimMobilityStrategyFiniteMob : public KimMobilityStrategy<FreeEnergyType>
{
 public:
   KimMobilityStrategyFiniteMob(
       QuatModel* quat_model, const int conc_l_id, const int conc_s_id,
       const int temp_id, const double interface_mobility, const double epsilon,
       const double phase_well_scale,
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       std::shared_ptr<tbox::Database> conc_db, const unsigned ncompositions);

 private:
   double evaluateMobility(const double temp,
                           const std::vector<double>& phaseconc,
                           const std::vector<double>& phi)
   {
      (void)temp;
      (void)phaseconc;
      (void)phi;

      return 1. / d_alpha;
   }

   double d_alpha;
};

#endif
