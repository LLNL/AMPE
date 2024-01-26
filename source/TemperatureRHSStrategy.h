// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef included_TemperatureRHSStrategy
#define included_TemperatureRHSStrategy

#include "SAMRAI/hier/PatchHierarchy.h"

#include <memory>

using namespace SAMRAI;

class TemperatureRHSStrategy
{
 public:
   void evaluateRHS(std::shared_ptr<hier::PatchHierarchy> hierarchy,
                    const int ydot_phase_id, const int dphidt_id);
   virtual void evaluateRHS(std::shared_ptr<hier::Patch> hierarchy,
                            const int ydot_phase_id, const int dphidt_id) = 0;
};

#endif
