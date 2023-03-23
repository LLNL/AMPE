// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#include "QuatFaceCoeffs.h"
#include "QuatFort.h"

#include "SAMRAI/math/HierarchySideDataOpsReal.h"

using namespace SAMRAI;

QuatFaceCoeffs::QuatFaceCoeffs(const int qlen, const double epsilon_q,
                               const double gradient_floor,
                               const std::string grad_floor_type)
    : d_qlen(qlen),
      d_epsilon_q(epsilon_q),
      d_gradient_floor(gradient_floor),
      d_grad_floor_type(grad_floor_type)
{
}

// Evaluate coefficient eps_q^2+D_q(phi)/|nabla q|
void QuatFaceCoeffs::computeCoeffs(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int diffusion_coef_id, const int grad_q_id, const int face_coef_id)
{
   // Check for negative diffusion coefficients.  We don't like them.
   // Zero is ok since epsilon^2 is added later
   math::HierarchySideDataOpsReal<double> sideops(
       hierarchy, 0, -hierarchy->getFinestLevelNumber());
   double diffusion_coef_min = sideops.min(diffusion_coef_id);
   if ((diffusion_coef_min + d_epsilon_q * d_epsilon_q) < 0.) {
      TBOX_ERROR(
          "QuatFaceCoeffs::computeCoeffs(), Negative diffusion coefficient!");
   }

   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator pi(level->begin()); pi != level->end();
           pi++) {
         std::shared_ptr<hier::Patch> patch = *pi;

         std::shared_ptr<pdat::SideData<double> > diffusion_coef_data(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(diffusion_coef_id)));
         std::shared_ptr<pdat::SideData<double> > grad_q_data(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(grad_q_id)));
         std::shared_ptr<pdat::SideData<double> > face_coef_data(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(face_coef_id)));

         assert(diffusion_coef_data->getDepth() == 1);
         computeCoeffs(*patch, *diffusion_coef_data, *grad_q_data,
                       *face_coef_data);
      }
   }
}

// face_coef_data: output
void QuatFaceCoeffs::computeCoeffs(const hier::Patch& patch,
                                   pdat::SideData<double>& diffusion_coef_data,
                                   pdat::SideData<double>& grad_q_data,
                                   pdat::SideData<double>& face_coef_data)
{
   assert(patch.inHierarchy());
   assert(diffusion_coef_data.getDepth() == 1);
   assert(grad_q_data.getDepth() == NDIM * d_qlen);
   assert(face_coef_data.getDepth() == d_qlen);

   const hier::Box& box = patch.getBox();
   const hier::Index& lower = box.lower();
   const hier::Index& upper = box.upper();

   const hier::Box& dc_gbox = diffusion_coef_data.getGhostBox();
   const hier::Index& dcglower = dc_gbox.lower();
   const hier::Index& dcgupper = dc_gbox.upper();

   const hier::Box& gq_gbox = grad_q_data.getGhostBox();
   const hier::Index& gqlower = gq_gbox.lower();
   const hier::Index& gqupper = gq_gbox.upper();

   const hier::Box& d_gbox = face_coef_data.getGhostBox();
   const hier::Index& dlower = d_gbox.lower();
   const hier::Index& dupper = d_gbox.upper();

#if NDIM == 2
   COMPUTE_FACE_COEF2D(lower[0], upper[0], lower[1], upper[1], d_qlen,
                       d_epsilon_q, diffusion_coef_data.getPointer(0),
                       dcglower[0], dcgupper[0] + 1, dcglower[1], dcgupper[1],
                       diffusion_coef_data.getPointer(1), dcglower[0],
                       dcgupper[0], dcglower[1], dcgupper[1] + 1,
                       grad_q_data.getPointer(0), gqlower[0], gqupper[0] + 1,
                       gqlower[1], gqupper[1], grad_q_data.getPointer(1),
                       gqlower[0], gqupper[0], gqlower[1], gqupper[1] + 1,
                       face_coef_data.getPointer(0), dlower[0], dupper[0] + 1,
                       dlower[1], dupper[1],  // output
                       face_coef_data.getPointer(1), dlower[0], dupper[0],
                       dlower[1], dupper[1] + 1,  // output
                       d_gradient_floor, d_grad_floor_type.c_str());
#endif
#if NDIM == 3
   COMPUTE_FACE_COEF3D(
       lower[0], upper[0], lower[1], upper[1], lower[2], upper[2], d_qlen,
       d_epsilon_q, diffusion_coef_data.getPointer(0), dcglower[0],
       dcgupper[0] + 1, dcglower[1], dcgupper[1], dcglower[2], dcgupper[2],
       diffusion_coef_data.getPointer(1), dcglower[0], dcgupper[0], dcglower[1],
       dcgupper[1] + 1, dcglower[2], dcgupper[2],
       diffusion_coef_data.getPointer(2), dcglower[0], dcgupper[0], dcglower[1],
       dcgupper[1], dcglower[2], dcgupper[2] + 1, grad_q_data.getPointer(0),
       gqlower[0], gqupper[0] + 1, gqlower[1], gqupper[1], gqlower[2],
       gqupper[2], grad_q_data.getPointer(1), gqlower[0], gqupper[0],
       gqlower[1], gqupper[1] + 1, gqlower[2], gqupper[2],
       grad_q_data.getPointer(2), gqlower[0], gqupper[0], gqlower[1],
       gqupper[1], gqlower[2], gqupper[2] + 1, face_coef_data.getPointer(0),
       dlower[0], dupper[0] + 1, dlower[1], dupper[1], dlower[2], dupper[2],
       face_coef_data.getPointer(1), dlower[0], dupper[0], dlower[1],
       dupper[1] + 1, dlower[2], dupper[2], face_coef_data.getPointer(2),
       dlower[0], dupper[0], dlower[1], dupper[1], dlower[2], dupper[2] + 1,
       d_gradient_floor, d_grad_floor_type.c_str());
#endif
}
