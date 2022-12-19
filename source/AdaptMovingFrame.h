// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef included_AdaptMovingFrame
#define included_AdaptMovingFrame

#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/PatchHierarchy.h"

class QuatModel;

class AdaptMovingFrame
{
 public:
   AdaptMovingFrame(
       std::shared_ptr<SAMRAI::geom::CartesianGridGeometry> grid_geometry,
       QuatModel* quat_model, const double frame_velocity);

   double adaptVelocity(
       const double time,
       const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
       const double current_frame_velocity, const int phase_id);

 private:
   // 1./volume of computational domain
   double d_vol_inv;

   // length of domain in x-direction
   double d_Lx;

   // safety limit on acceleration
   double d_max_accel;

   // internal values to keep track of previous steps
   double d_previous_time;
   double d_previous_fs;

   // model to calculate solid fraction
   QuatModel* d_quat_model;
};

#endif
