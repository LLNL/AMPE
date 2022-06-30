// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
#include "InterpolationType.h"

#include "SAMRAI/tbox/Utilities.h"

namespace ampe_thermo
{

char concInterpChar(ConcInterpolationType interp_func_type)
{
   switch (interp_func_type) {
      case ConcInterpolationType::LINEAR: return 'l';
      case ConcInterpolationType::PBG: return 'p';
      case ConcInterpolationType::HARMONIC: return 'h';
      default: return '0';
   }
}

char energyInterpChar(EnergyInterpolationType interp_func_type)
{
   switch (interp_func_type) {
      case EnergyInterpolationType::LINEAR: return 'l';
      case EnergyInterpolationType::PBG: return 'p';
      case EnergyInterpolationType::HARMONIC: return 'h';
      default: return '0';
   }
}

}  // namespace ampe_thermo
