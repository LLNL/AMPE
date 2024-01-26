// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
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
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/hier/PatchHierarchy.h"


#include <vector>
#include <string>

using namespace SAMRAI;

class QuatFaceCoeff
{

 public:
   QuatFaceCoeff(const int qlen, const double epsilon_q,
                 const double gradient_floor, const std::string grad_floor_type,
                 const double Hparameter, const std::string interp_type1,
                 const std::string interp_type2, const std::string avg_type);

   ~QuatFaceCoeff();

   void computeFaceCoefs(std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
                         const int phase_id, const int temp_id,
                         const int grad_q_id,
                         const int face_coef_id);  // output

 private:
   const int d_qlen;
   // norm(grad q)^2 coefficient
   const double d_epsilon_q;
   const double d_gradient_floor;
   const std::string d_grad_floor_type;
   const double d_Hparameter;
   // interpolation polynomial for phase fraction in |\nabla q| term
   const std::string d_interp_type1;
   // interpolation polynomial for phase fraction in |\nabla q|^2 term
   const std::string d_interp_type2;
   const std::string d_avg_type;

   void computeFaceCoefsOnPatch(
       const hier::Patch& patch, pdat::CellData<double>& phase_data,
       pdat::CellData<double>& temp_data, pdat::SideData<double>& grad_q_data,
       pdat::SideData<double>& face_coef_data) const;  // output
};
#endif
