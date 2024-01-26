// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef included_PhaseRHSStrategy
#define included_PhaseRHSStrategy

#include "SAMRAI/hier/PatchHierarchy.h"

#include <memory>

using namespace SAMRAI;

class PhaseRHSStrategy
{
 public:
   virtual void evaluateRHS(const double time,
                            std::shared_ptr<hier::PatchHierarchy> hierarchy,
                            const int ydot_phase_id, const bool eval_flag) = 0;
   virtual void setup(std::shared_ptr<hier::PatchHierarchy> hierarchy) = 0;
};

#endif
