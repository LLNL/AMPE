// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef included_WangRigidBodyForces
#define included_WangRigidBodyForces

#include "QuatModelParameters.h"
#include "RigidBodyForces.h"

#include <vector>
#include <array>

class WangRigidBodyForces : public RigidBodyForces
{
 public:
   WangRigidBodyForces(const QuatModelParameters& model_parameters,
                       const int phase_scratch_id, const int conc_id,
                       const int weight_id);

 private:
   const int d_conc_id;

   /*
    * Number of order parameters to include in body forces
    */
   short d_norderp;

   void evaluatePairForces(std::shared_ptr<hier::Patch> patch,
                           std::vector<double>& forces);
};

#endif
