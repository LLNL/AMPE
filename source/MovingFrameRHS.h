// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
#ifndef included_MovingFrameRHS
#define included_MovingFrameRHS

#include "SAMRAI/hier/PatchHierarchy.h"

using namespace SAMRAI;

class MovingFrameRHS
{
 public:
   MovingFrameRHS(const int phase_scratch_id);

   void addRHS(std::shared_ptr<hier::PatchHierarchy> hierarchy,
               const int ydot_phase_id, const double frame_velocity);

 private:
   const int d_phase_scratch_id;
};

#endif
