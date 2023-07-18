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
#include "QuatFaceCoeff.h"
#include "QuatFort.h"

#include "SAMRAI/math/HierarchySideDataOpsReal.h"

#include <cassert>

QuatFaceCoeff::QuatFaceCoeff(const int qlen, const double epsilon_q,
                             const double gradient_floor,
                             const std::string grad_floor_type,
                             const double Hparameter,
                             const std::string interp_type1,
                             const std::string interp_type2,
                             const std::string avg_type)
    : d_qlen(qlen),
      d_epsilon_q(epsilon_q),
      d_gradient_floor(gradient_floor),
      d_grad_floor_type(grad_floor_type),
      d_Hparameter(Hparameter),
      d_interp_type1(interp_type1),
      d_interp_type2(interp_type2),
      d_avg_type(avg_type)
{
   assert(d_epsilon_q > 0.);
}

QuatFaceCoeff::~QuatFaceCoeff() {}

void QuatFaceCoeff::computeFaceCoefs(
    std::shared_ptr<hier::PatchHierarchy> hierarchy, const int phase_id,
    const int temp_id, const int grad_q_id,
    const int face_coef_id)  // output
{
   for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator pi(level->begin()); pi != level->end();
           pi++) {
         std::shared_ptr<hier::Patch> patch = *pi;

         std::shared_ptr<pdat::CellData<double> > phase_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(phase_id)));
         std::shared_ptr<pdat::CellData<double> > temp_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(temp_id)));
         std::shared_ptr<pdat::SideData<double> > grad_q_data(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(grad_q_id)));
         std::shared_ptr<pdat::SideData<double> > face_coef_data(
             SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
                 patch->getPatchData(face_coef_id)));

         computeFaceCoefsOnPatch(*patch, *phase_data, *temp_data, *grad_q_data,
                                 *face_coef_data);
      }
   }
#if 0
   double norm = sideops.L1Norm( face_coef_id );
   tbox::plog<<"L1 Norm face_coef_id="<<norm<<endl;
#endif
}


/*
*******************************************************************
*                                                                 *
* AMR-unaware patch-centered computational kernels.               *
*                                                                 *
*******************************************************************
*/
void QuatFaceCoeff::computeFaceCoefsOnPatch(
    const hier::Patch& patch, pdat::CellData<double>& phase_data,
    pdat::CellData<double>& temp_data, pdat::SideData<double>& grad_q_data,
    pdat::SideData<double>& face_coef_data) const  // output
{
   assert(d_Hparameter > 0.);

   assert(patch.inHierarchy());
   assert(grad_q_data.getDepth() == NDIM * d_qlen);
   assert(face_coef_data.getDepth() == 1);

   const hier::Box& box = patch.getBox();
   const hier::Index& lower = box.lower();
   const hier::Index& upper = box.upper();

#if NDIM == 2
   COMPUTE_FACE_COEF2D(
       lower[0], upper[0], lower[1], upper[1], d_qlen, d_epsilon_q,
       phase_data.getPointer(), phase_data.getGhostCellWidth()[0],
       temp_data.getPointer(), temp_data.getGhostCellWidth()[0],
       2. * d_Hparameter, grad_q_data.getPointer(0), grad_q_data.getPointer(1),
       grad_q_data.getGhostCellWidth()[0], face_coef_data.getPointer(0),
       face_coef_data.getPointer(1),
       face_coef_data.getGhostCellWidth()[0],  // output
       d_gradient_floor, d_grad_floor_type.c_str(), d_interp_type1.c_str(),
       d_interp_type2.c_str(), d_avg_type.c_str());
#endif
#if NDIM == 3
   COMPUTE_FACE_COEF3D(
       lower[0], upper[0], lower[1], upper[1], lower[2], upper[2], d_qlen,
       d_epsilon_q, phase_data.getPointer(), phase_data.getGhostCellWidth()[0],
       temp_data.getPointer(), temp_data.getGhostCellWidth()[0],
       2. * d_Hparameter, grad_q_data.getPointer(0), grad_q_data.getPointer(1),
       grad_q_data.getPointer(2), grad_q_data.getGhostCellWidth()[0],
       face_coef_data.getPointer(0), face_coef_data.getPointer(1),
       face_coef_data.getPointer(2), face_coef_data.getGhostCellWidth()[0],
       d_gradient_floor, d_grad_floor_type.c_str(), d_interp_type1.c_str(),
       d_interp_type2.c_str(), d_avg_type.c_str());
#endif
}
