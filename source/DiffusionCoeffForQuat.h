// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef included_DiffusionCoeffForQuat
#define included_DiffusionCoeffForQuat

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

#include <cstring>

using namespace SAMRAI;

class DiffusionCoeffForQuat
{
 public:
   DiffusionCoeffForQuat(
       std::shared_ptr<geom::CartesianGridGeometry> grid_geometry,
       const double H_parameter, const double quat_grad_floor, const int qlen,
       std::string orient_interp_func_type, std::string avg_func_type,
       std::string quat_smooth_floor_type, const bool precond_has_dquatdphi,
       std::shared_ptr<pdat::SideVariable<double> > quat_diffusion_var,
       const int quat_diffusion_id, const int quat_diffusion_deriv_id);

   // set variables that may change over time
   void setup(const std::shared_ptr<hier::PatchHierarchy> hierarchy);

   // calculate diffusion coefficient
   // output in d_quat_diffusion_id
   void setDiffusion(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                     const int phase_scratch_id,
                     const int temperature_scratch_id);

   void setDerivDiffusion(const std::shared_ptr<hier::PatchHierarchy> hierarchy,
                          const int phase_scratch_id,
                          const int temperature_scratch_id,
                          const int quat_grad_side_id);

 private:
   double d_H_parameter;
   double d_quat_grad_floor;
   int d_qlen;

   std::string d_orient_interp_func_type;
   std::string d_avg_func_type;
   std::string d_quat_smooth_floor_type;

   bool d_precond_has_dquatdphi;

   int d_quat_diffusion_id;
   int d_quat_diffusion_deriv_id;

   xfer::CoarsenAlgorithm d_quat_diffusion_coarsen;
   xfer::CoarsenAlgorithm d_quat_diffusion_deriv_coarsen;

   std::vector<std::shared_ptr<xfer::CoarsenSchedule> >
       d_quat_diffusion_coarsen_schedule;
   std::vector<std::shared_ptr<xfer::CoarsenSchedule> >
       d_quat_diffusion_deriv_coarsen_schedule;

   void setDiffusionPatch(hier::Patch& patch, const int phase_scratch_id,
                          const int temperature_scratch_id);
   void setDerivDiffusionPatch(hier::Patch& patch, const int phase_scratch_id,
                               const int temperature_scratch_id,
                               const int quat_grad_side_id);
};

#endif
