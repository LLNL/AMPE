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
//
#ifndef included_QuatFaceCoeff_h
#define included_QuatFaceCoeff_h

#ifndef included_SAMRAI_config
#include "SAMRAI/SAMRAI_config.h"
#endif

#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/PatchHierarchy.h"


#include <vector>
#include <string>

using namespace SAMRAI;

class QuatFaceCoeff
{

 public:
   QuatFaceCoeff(const int qlen, const double epsilon_q,
                 const double gradient_floor,
                 const std::string d_grad_floor_type);

   ~QuatFaceCoeff();

   void computeFaceCoefs(std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
                         const int diffusion_coef_id, const int grad_q_id,
                         const int face_coef_id);  // output

 private:
   const int d_qlen;
   // norm(grad q)^2 coefficient
   const double d_epsilon_q;
   const double d_gradient_floor;
   const std::string d_grad_floor_type;

   void computeFaceCoefsOnPatch(
       const hier::Patch& patch, pdat::SideData<double>& diffusion_coef_data,
       pdat::SideData<double>& grad_q_data,
       pdat::SideData<double>& face_coef_data) const;  // output
};
#endif
