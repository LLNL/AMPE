// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#include "QuatModelParameters.h"

// see Kim. Acta Mat. 55 (2007), 4391
double kks_mobility_factor(
    const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
    const double epsilon, const double phase_well_scale)
{
   double a2 = 0.;
   switch (energy_interp_func_type) {
      case Thermo4PFM::EnergyInterpolationType::PBG: a2 = 47. / 60.; break;
      case Thermo4PFM::EnergyInterpolationType::HARMONIC: a2 = 0.5; break;
      default:
         TBOX_ERROR("Invalid interpolation function in mobility_factor()");
   }

   // factor 16 due to AMPE convention for phase_well_scale
   double mobility_factor =
       16. * phase_well_scale / (3. * epsilon * epsilon * a2);

   return mobility_factor;
}
