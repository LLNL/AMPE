// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef included_RigidBodyForces
#define included_RigidBodyForces

#include "QuatModelParameters.h"

#include <vector>
#include <array>

class RigidBodyForces
{
 public:
   RigidBodyForces(const QuatModelParameters& model_parameters,
                   const int phase_scratch_id, const int weight_id);

   virtual ~RigidBodyForces(){};

   void evaluatePairForces(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy);

   void printPairForces(std::ostream& os);

   void getPairForce(const short i, const short j, std::array<double, NDIM>& f);

 protected:
   const QuatModelParameters& d_model_parameters;

   const int d_phase_id;
   const int d_weight_id;

   /*
    * force acting on each particle, from all the other particles
    */
   std::vector<std::vector<std::array<double, NDIM>>> d_forces;

   virtual void evaluatePairForces(std::shared_ptr<hier::Patch> patch,
                                   std::vector<double>& forces) = 0;

 private:
   // Timers
   std::shared_ptr<tbox::Timer> t_eval_pair_forces_timer;
};

#endif
