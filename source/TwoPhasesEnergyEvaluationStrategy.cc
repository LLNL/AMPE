// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.

#include "TwoPhasesEnergyEvaluationStrategy.h"
#include "QuatFort.h"

#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/math/PatchSideDataNormOpsReal.h"
#include "SAMRAI/math/PatchCellDataNormOpsReal.h"

using namespace SAMRAI;

TwoPhasesEnergyEvaluationStrategy::TwoPhasesEnergyEvaluationStrategy(
    const QuatModelParameters& model_parameters, const int qlen,
    const int phase_scratch_id, const int quat_scratch_id,
    const int quat_grad_side_id, const int weight_id, const int f_l_id,
    const int f_a_id, const int temperature_id, const int energy_diag_id)
    : EnergyEvaluationStrategy(model_parameters, energy_diag_id),
      d_model_parameters(model_parameters),
      d_qlen(qlen),
      d_phase_scratch_id(phase_scratch_id),
      d_quat_scratch_id(quat_scratch_id),
      d_quat_grad_side_id(quat_grad_side_id),
      d_weight_id(weight_id),
      d_f_l_id(f_l_id),
      d_f_a_id(f_a_id),
      d_temperature_id(temperature_id),
      d_energy_diag_id(energy_diag_id)
{
}

void TwoPhasesEnergyEvaluationStrategy::evaluateEnergy(
    std::shared_ptr<hier::Patch> patch, const double time, double& total_energy,
    double& total_interface_e, double& total_orient_e, double& total_qint_e,
    double& total_well_e, double& bulk_energy, const bool gp)
{
   std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch->getPatchGeometry()));
   TBOX_ASSERT(patch_geom);

   const double* dx = patch_geom->getDx();

   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   double* pgrad_quat[NDIM];
   if (d_model_parameters.evolveQuat()) {
      std::shared_ptr<pdat::SideData<double> > grad_quat(
          SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
              patch->getPatchData(d_quat_grad_side_id)));
      assert(grad_quat);
      assert(grad_quat->getGhostCellWidth() ==
             hier::IntVector(tbox::Dimension(NDIM), 0));
      for (int d = 0; d < NDIM; d++) {
         pgrad_quat[d] = grad_quat->getPointer(d);
      }
#ifdef DEBUG_CHECK_ASSERTIONS
      SAMRAI::math::PatchSideDataNormOpsReal<double> sops;
      double l2gq = sops.L2Norm(grad_quat, pbox);
      assert(l2gq == l2gq);
#endif
   } else {
      for (int d = 0; d < NDIM; d++) {
         pgrad_quat[d] = nullptr;
      }
   }
   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_phase_scratch_id)));
   std::shared_ptr<pdat::CellData<double> > weight(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_weight_id)));
   std::shared_ptr<pdat::CellData<double> > fl(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_f_l_id)));
   std::shared_ptr<pdat::CellData<double> > fa(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_f_a_id)));
   std::shared_ptr<pdat::CellData<double> > temperature(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_temperature_id)));

   double* quat_ptr = nullptr;
   int ngquat = 0;
   const double epsilon_anisotropy = d_model_parameters.epsilon_anisotropy();
   if (epsilon_anisotropy >= 0.) {
      std::shared_ptr<pdat::CellData<double> > quat(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_quat_scratch_id)));
      quat_ptr = quat->getPointer();
      ngquat = quat->getGhostCellWidth()[0];
      assert(quat_ptr != nullptr);
   }

   assert(phase);
   assert(weight);
   assert(fl);
   assert(fa);

#ifdef DEBUG_CHECK_ASSERTIONS
   math::PatchCellDataNormOpsReal<double> ops;
   double l2phi = ops.L2Norm(phase, pbox);
   assert(l2phi == l2phi);

   double l2t = ops.L2Norm(temperature, pbox);
   assert(l2t == l2t);

   double l2fl = ops.L2Norm(fl, pbox);
   assert(l2fl == l2fl);

   double l2fa = ops.L2Norm(fa, pbox);
   assert(l2fa == l2fa);
#endif
   int per_cell = 0;
   double* ptr_energy = nullptr;
   if (d_model_parameters.with_visit_energy_output()) {
      per_cell = 1;
      assert(d_energy_diag_id >= 0);
      std::shared_ptr<pdat::CellData<double> > energy(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_energy_diag_id)));
      ptr_energy = energy->getPointer();
   }

   assert(weight->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));

   if (d_model_parameters.evolveQuat())
      for (int i = 0; i < NDIM; i++)
         assert(pgrad_quat[i] != nullptr);

   const char interpf =
       energyInterpChar(d_model_parameters.energy_interp_func_type());
   // printf("orient_interp_func_type2 =%s\n",
   // d_model_parameters.orient_interp_func_type2().c_str()[0]);
   QUATENERGY(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
              ifirst(2), ilast(2),
#endif
              d_qlen, dx, pgrad_quat[0], pgrad_quat[1],
#if (NDIM == 3)
              pgrad_quat[2],
#endif
              0, phase->getPointer(), phase->getGhostCellWidth()[0], 1,
              quat_ptr, ngquat, d_model_parameters.epsilon_phase(),
              d_model_parameters.epsilon_q(), epsilon_anisotropy, 4,
              2. * d_model_parameters.H_parameter(), temperature->getPointer(),
              temperature->getGhostCellWidth()[0],
              d_model_parameters.phase_well_scale(), weight->getPointer(),
              total_energy, total_interface_e, total_orient_e, total_qint_e,
              ptr_energy, per_cell, &interpf,
              d_model_parameters.phase_well_func_type().c_str(),
              d_model_parameters.orient_interp_func_type1().c_str(),
              d_model_parameters.orient_interp_func_type2().c_str(),
              d_model_parameters.avg_func_type().c_str(),
              d_model_parameters.quat_grad_floor_type().c_str(),
              d_model_parameters.quat_grad_floor());

   BULKENERGY(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
              ifirst(2), ilast(2),
#endif
              phase->getPointer(), phase->getGhostCellWidth()[0],
              fl->getPointer(), fa->getPointer(), weight->getPointer(),
              total_energy, bulk_energy, ptr_energy, per_cell, &interpf);
}
