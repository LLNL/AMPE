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
#ifndef ThreeArgsInterpolationType_H
#define ThreeArgsInterpolationType_H

#include "InterpolationType.h"

enum class ConcThreeArgsInterpolationType {
   FOLCHPLAPP2005,
   MOELANS2011,
   UNDEFINED
};

enum class EnergyThreeArgsInterpolationType {
   FOLCHPLAPP2005,
   MOELANS2011,
   UNDEFINED
};

Thermo4PFM::EnergyInterpolationType getTwoPhasesInterpolationType(
    EnergyThreeArgsInterpolationType type);

#endif
