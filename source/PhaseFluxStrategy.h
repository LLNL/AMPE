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
//
#ifndef included_PhaseFluxStrategy
#define included_PhaseFluxStrategy

#include "SAMRAI/hier/PatchLevel.h"

using namespace SAMRAI;

class PhaseFluxStrategy
{
 public:
   PhaseFluxStrategy(){};

   virtual ~PhaseFluxStrategy(){};

   virtual void computeFluxes(const std::shared_ptr<hier::PatchLevel> level,
                              const int phase_id, const int quat_id,
                              const int flux_id) = 0;
};

#endif
