// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#include "AdaptMovingFrame.h"
#include "QuatModel.h"

using namespace SAMRAI;

AdaptMovingFrame::AdaptMovingFrame(
    std::shared_ptr<geom::CartesianGridGeometry> grid_geometry,
    QuatModel* quat_model, const double frame_velocity)
    : d_max_accel(std::abs(frame_velocity)), d_quat_model(quat_model)
{
   const double* low = grid_geometry->getXLower();
   const double* up = grid_geometry->getXUpper();
   double vol = 1.;
   for (int d = 0; d < NDIM; d++)
      vol *= (up[d] - low[d]);
   d_vol_inv = 1. / vol;
   d_Lx = up[0] - low[0];
}

double AdaptMovingFrame::adaptVelocity(
    const double time, const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
    const double current_frame_velocity, const int phase_id)
{
   // compute current solid fraction
   double fs =
       d_quat_model->evaluateVolumeSolid(hierarchy, phase_id) * d_vol_inv;
   tbox::pout << "Moving frame: fs = " << fs << std::endl;

   static bool first_time = true;
   double acceleration = 0.;
   if (!first_time) {
      double dfsk = (fs - d_previous_fs);
      double dtinv = 1. / (time - d_previous_time);

      // difference in velocity between current frame velocity
      // and moving front velocity
      // Properties:
      // * would be 0 if fs did not change from previous value
      // * larger than 0 if fs growing (current frame velocity smaller
      //   than moving front velocity)
      double deltavk = d_Lx * dfsk * dtinv;

      acceleration = deltavk;
      if (std::abs(acceleration) > d_max_accel)
         acceleration *= (d_max_accel / std::abs(acceleration));
   }
   tbox::pout << "Moving frame: acceleration = " << acceleration << std::endl;
   double frame_velocity = current_frame_velocity + acceleration;
   tbox::pout << "Moving frame: new velocity = " << frame_velocity << std::endl;

   first_time = false;

   // save current values for next call
   d_previous_time = time;
   d_previous_fs = fs;

   return frame_velocity;
}
