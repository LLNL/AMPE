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
#include "QuatModelParameters.h"

double kks_mobility_factor(
    const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
    const double epsilon, const double phase_well_scale);
