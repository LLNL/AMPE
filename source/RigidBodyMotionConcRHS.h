// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
#ifndef included_RigidBodyMotionConcRHS
#define included_RigidBodyMotionConcRHS

#include "SAMRAI/hier/PatchHierarchy.h"

#include <vector>
#include <array>

using namespace SAMRAI;

class RigidBodyMotionConcRHS
{
 public:
   RigidBodyMotionConcRHS(const int data_scratch_id, int norderp,
                          const int weight_id, const double mobility);

   void addRHS(std::shared_ptr<hier::PatchHierarchy> hierarchy,
               const int ydot_id,
               const std::vector<std::array<double, NDIM>>& forces);

 private:
   const int d_data_id;
   const int d_norderp;

   const int d_weight_id;

   /*!
    * Translational mobilities relating velocities to forces
    */
   const double d_mobility;

   /*!
    * Volumes of each grain
    */
   std::vector<double> d_volumes;

   void computeVolumes(std::shared_ptr<hier::PatchHierarchy> hierarchy);
};

#endif
