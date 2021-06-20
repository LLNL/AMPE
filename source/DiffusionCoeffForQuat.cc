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
#include "DiffusionCoeffForQuat.h"
#include "QuatFort.h"
#include "QuatParams.h"

#include "SAMRAI/pdat/CellData.h"

DiffusionCoeffForQuat::DiffusionCoeffForQuat(
    std::shared_ptr<geom::CartesianGridGeometry> grid_geometry,
    const double H_parameter, const double quat_grad_floor, const int qlen,
    std::string orient_interp_func_type, std::string avg_func_type,
    std::string quat_smooth_floor_type, const bool precond_has_dquatdphi,
    std::shared_ptr<pdat::SideVariable<double> > quat_diffusion_var,
    const int quat_diffusion_id, const int quat_diffusion_deriv_id)
    : d_H_parameter(H_parameter),
      d_quat_grad_floor(quat_grad_floor),
      d_qlen(qlen),
      d_orient_interp_func_type(orient_interp_func_type),
      d_avg_func_type(avg_func_type),
      d_quat_smooth_floor_type(quat_smooth_floor_type),
      d_precond_has_dquatdphi(precond_has_dquatdphi),
      d_quat_diffusion_id(quat_diffusion_id),
      d_quat_diffusion_deriv_id(quat_diffusion_deriv_id),
      d_quat_diffusion_coarsen(tbox::Dimension(NDIM)),
      d_quat_diffusion_deriv_coarsen(tbox::Dimension(NDIM))
{
   assert(d_quat_diffusion_id >= 0);

   std::shared_ptr<hier::CoarsenOperator> diff_coarsen_op =
       grid_geometry->lookupCoarsenOperator(quat_diffusion_var,
                                            "CONSERVATIVE_COARSEN");

   d_quat_diffusion_coarsen.registerCoarsen(d_quat_diffusion_id,
                                            d_quat_diffusion_id,
                                            diff_coarsen_op);
   if (d_precond_has_dquatdphi)
      d_quat_diffusion_deriv_coarsen.registerCoarsen(d_quat_diffusion_deriv_id,
                                                     d_quat_diffusion_deriv_id,
                                                     diff_coarsen_op);
}

void DiffusionCoeffForQuat::setup(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   d_quat_diffusion_coarsen_schedule.resize(hierarchy->getNumberOfLevels());
   if (d_precond_has_dquatdphi)
      d_quat_diffusion_deriv_coarsen_schedule.resize(
          hierarchy->getNumberOfLevels());

   for (int ln = 0; ln < hierarchy->getNumberOfLevels() - 1; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
      std::shared_ptr<hier::PatchLevel> finer_level =
          hierarchy->getPatchLevel(ln + 1);

      d_quat_diffusion_coarsen_schedule[ln] =
          d_quat_diffusion_coarsen.createSchedule(level, finer_level);
      if (d_precond_has_dquatdphi) {
         d_quat_diffusion_deriv_coarsen_schedule[ln] =
             d_quat_diffusion_deriv_coarsen.createSchedule(level, finer_level);
      }
   }
}

void DiffusionCoeffForQuat::setDiffusion(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int phase_scratch_id, const int temperature_scratch_id)
{
   assert(hierarchy);

   // set diffusion coefficients
   const int maxl = hierarchy->getNumberOfLevels();
   for (int amr_level = 0; amr_level < maxl; amr_level++) {
      std::shared_ptr<hier::PatchLevel> level =
          hierarchy->getPatchLevel(amr_level);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         setDiffusionPatch(**p, phase_scratch_id, temperature_scratch_id);
      }
   }

   // coarsen diffusion coefficients
   for (int amr_level = hierarchy->getFinestLevelNumber() - 1; amr_level >= 0;
        amr_level--) {
      d_quat_diffusion_coarsen_schedule[amr_level]->coarsenData();
   }
}

//-----------------------------------------------------------------------
// set diffusion coefficient to 2*H*T*[1-p(phi)] as in Pusztai et al.
//
void DiffusionCoeffForQuat::setDiffusionPatch(hier::Patch& patch,
                                              const int phase_scratch_id,
                                              const int temperature_scratch_id)
{
   // tbox::pout<<"QuatIntegrator::setDiffusionCoeffForQuatPatch()"<<endl;
   assert(phase_scratch_id >= 0);
   assert(temperature_scratch_id >= 0);
   assert(d_quat_diffusion_id >= 0);

   const hier::Box& pbox = patch.getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   std::shared_ptr<pdat::CellData<double> > phi(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(phase_scratch_id)));

   std::shared_ptr<pdat::CellData<double> > temperature(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_scratch_id)));

   std::shared_ptr<pdat::SideData<double> > diffusion(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(d_quat_diffusion_id)));

   assert(diffusion->getDepth() == 1);
   assert(phi);
   assert(diffusion);

#ifdef DEBUG_CHECK_ASSERTIONS
   const hier::Box& phi_gbox = phi->getGhostBox();
   const hier::Index& philower = phi_gbox.lower();
   const hier::Index& phiupper = phi_gbox.upper();
   for (int i = 0; i < NDIM; i++) {
      assert(philower[i] <= ifirst(i));
      assert(phiupper[i] >= ilast(i));
   }

   const hier::Box& diff_gbox = diffusion->getGhostBox();
   const hier::Index& difflower = diff_gbox.lower();
   const hier::Index& diffupper = diff_gbox.upper();
   for (int i = 0; i < NDIM; i++) {
      assert(difflower[i] <= ifirst(i));
      assert(diffupper[i] >= ilast(i));
   }

   assert(diffusion->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));
   assert(phi->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), NGHOSTS));
   assert(temperature->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), NGHOSTS));

   for (int d = 0; d < NDIM; d++) {
      assert(diffusion->getPointer(d) != nullptr);
   }

#endif

   QUATDIFFUSION(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                 ifirst(2), ilast(2),
#endif
                 2. * d_H_parameter, temperature->getPointer(),
                 temperature->getGhostCellWidth()[0], phi->getPointer(),
                 NGHOSTS, diffusion->getPointer(0), diffusion->getPointer(1),
#if (NDIM == 3)
                 diffusion->getPointer(2),
#endif
                 0, d_orient_interp_func_type.c_str(), d_avg_func_type.c_str());
}

void DiffusionCoeffForQuat::setDerivDiffusion(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int phase_scratch_id, const int temperature_scratch_id,
    const int quat_grad_side_id)
{
   assert(hierarchy);
   assert(d_precond_has_dquatdphi);

   const int maxl = hierarchy->getNumberOfLevels();

   // set derivative of diffusion coefficients
   for (int amr_level = 0; amr_level < maxl; amr_level++) {
      std::shared_ptr<hier::PatchLevel> level =
          hierarchy->getPatchLevel(amr_level);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {
         std::shared_ptr<hier::Patch> patch = *p;

         setDerivDiffusionPatch(*patch, phase_scratch_id,
                                temperature_scratch_id, quat_grad_side_id);
      }
   }

   for (int amr_level = hierarchy->getFinestLevelNumber() - 1; amr_level >= 0;
        amr_level--) {
      d_quat_diffusion_deriv_coarsen_schedule[amr_level]->coarsenData();
   }
}

void DiffusionCoeffForQuat::setDerivDiffusionPatch(
    hier::Patch& patch, const int phase_scratch_id,
    const int temperature_scratch_id, const int quat_grad_side_id)
{
   assert(phase_scratch_id >= 0);
   assert(temperature_scratch_id >= 0);
   assert(d_quat_diffusion_deriv_id >= 0);
   assert(quat_grad_side_id >= 0);

   const hier::Box& pbox = patch.getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   std::shared_ptr<pdat::CellData<double> > phi(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(phase_scratch_id)));

   std::shared_ptr<pdat::CellData<double> > temperature(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch.getPatchData(temperature_scratch_id)));

   std::shared_ptr<pdat::SideData<double> > diffusion_deriv(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(d_quat_diffusion_deriv_id)));

   std::shared_ptr<pdat::SideData<double> > grad_q(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch.getPatchData(quat_grad_side_id)));

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(phi);
   const hier::Box& phi_gbox = phi->getGhostBox();
   const hier::Index& philower = phi_gbox.lower();
   const hier::Index& phiupper = phi_gbox.upper();
   for (int i = 0; i < NDIM; i++) {
      assert(philower[i] <= ifirst(i));
      assert(phiupper[i] >= ilast(i));
   }
   assert(phi->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), NGHOSTS));
   assert(temperature->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), NGHOSTS));

   assert(diffusion_deriv);
   assert(diffusion_deriv->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));
   const hier::Box& diffderiv_gbox = diffusion_deriv->getGhostBox();
   const hier::Index& diffderivlower = diffderiv_gbox.lower();
   const hier::Index& diffderivupper = diffderiv_gbox.upper();
   for (int i = 0; i < NDIM; i++) {
      assert(diffderivlower[i] <= ifirst(i));
      assert(diffderivupper[i] >= ilast(i));
   }

   for (int d = 0; d < NDIM; d++) {
      assert(diffusion_deriv->getPointer(d) != nullptr);
   }

   assert(grad_q);
#endif

   QUATDIFFUSIONDERIV(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                      ifirst(2), ilast(2),
#endif
                      2. * d_H_parameter, temperature->getPointer(),
                      temperature->getGhostCellWidth()[0], phi->getPointer(),
                      NGHOSTS, d_qlen, grad_q->getPointer(0),
                      grad_q->getPointer(1),
#if (NDIM == 3)
                      grad_q->getPointer(2),
#endif
                      0, diffusion_deriv->getPointer(0),
                      diffusion_deriv->getPointer(1),
#if (NDIM == 3)
                      diffusion_deriv->getPointer(2),
#endif
                      0, d_quat_grad_floor, d_quat_smooth_floor_type.c_str(),
                      d_orient_interp_func_type.c_str(),
                      d_avg_func_type.c_str());
}
