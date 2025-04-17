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
#ifndef MultiOrderBinaryThreePhasesDrivingForceStochioB_H
#define MultiOrderBinaryThreePhasesDrivingForceStochioB_H

#include "FreeEnergyStrategyBinary.h"

#include "SAMRAI/hier/Patch.h"

using namespace SAMRAI;

class MultiOrderBinaryThreePhasesDrivingForceStochioB
{
 public:
   MultiOrderBinaryThreePhasesDrivingForceStochioB(
       FreeEnergyStrategyBinary* fenergy_strategy, const int norderp_A);

   void addDrivingForce(hier::Patch& patch, const int temperature_id,
                        const int phase_id, const int conc_l_id,
                        const int conc_a_id, const int conc_b_id,
                        const int f_l_id, const int f_a_id, const int f_b_id,
                        const int rhs_id);

 private:
   FreeEnergyStrategyBinary* d_fenergy_strategy;

   int d_norderp_A;
};

#endif
