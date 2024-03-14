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
#include "ThreeArgsInterpolationType.h"
#include "InterpolationType.h"

// return two-phases interpolation type the three-phases
// interpolation reduces to when one phase is 0
Thermo4PFM::EnergyInterpolationType getTwoPhasesInterpolationType(
    EnergyThreeArgsInterpolationType type)
{
   if (type == EnergyThreeArgsInterpolationType::FOLCHPLAPP2005)
      return Thermo4PFM::EnergyInterpolationType::PBG;
   else
      return Thermo4PFM::EnergyInterpolationType::UNDEFINED;
}
