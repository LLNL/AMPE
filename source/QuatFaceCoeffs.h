// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef included_SAMRAI_config
#include "SAMRAI/SAMRAI_config.h"
#endif

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/pdat/SideData.h"

#include <string>
#include <memory>

class QuatFaceCoeffs
{
 public:
   QuatFaceCoeffs(const int qlen, const double epsilon_q,
                  const double gradient_floor,
                  const std::string grad_floor_type);

   // Evaluate coefficient eps_q^2+D_q(phi)/|nabla q|
   void computeCoeffs(
       const std::shared_ptr<SAMRAI::hier::PatchHierarchy> hierarchy,
       const int diffusion_coef_id, const int grad_q_id,
       const int face_coef_id);

 private:
   const int d_qlen;
   const double d_epsilon_q;
   const double d_gradient_floor;
   const std::string d_grad_floor_type;

   void computeCoeffs(const SAMRAI::hier::Patch& patch,
                      SAMRAI::pdat::SideData<double>& diffusion_coef_data,
                      SAMRAI::pdat::SideData<double>& grad_q_data,
                      SAMRAI::pdat::SideData<double>& face_coef_data);
};
