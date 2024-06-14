// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef included_PhaseRigidBodyForces
#define included_PhaseRigidBodyForces

#include "QuatModelParameters.h"
#include "RigidBodyForces.h"

#include <vector>
#include <array>

class PhaseRigidBodyForces : public RigidBodyForces
{
 public:
   PhaseRigidBodyForces(const QuatModelParameters& model_parameters,
                        const int phase_scratch_id, const int weight_id);

 private:
   void evaluatePairForces(std::shared_ptr<hier::Patch> patch,
                           std::vector<double>& forces);
};

#endif
