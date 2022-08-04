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
    const double H_parameter, const int qlen,
    std::string orient_interp_func_type, std::string avg_func_type,
    const int quat_diffusion_id)
    : d_H_parameter(H_parameter),
      d_qlen(qlen),
      d_orient_interp_func_type(orient_interp_func_type),
      d_avg_func_type(avg_func_type),
      d_quat_diffusion_id(quat_diffusion_id),
      d_quat_diffusion_coarsen(tbox::Dimension(NDIM))
{
   assert(d_quat_diffusion_id >= 0);

   hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();

   std::shared_ptr<hier::Variable> variable;
   vdb->mapIndexToVariable(d_quat_diffusion_id, variable);
   std::shared_ptr<hier::CoarsenOperator> diff_coarsen_op =
       grid_geometry->lookupCoarsenOperator(variable, "CONSERVATIVE_COARSEN");

   d_quat_diffusion_coarsen.registerCoarsen(d_quat_diffusion_id,
                                            d_quat_diffusion_id,
                                            diff_coarsen_op);
}

void DiffusionCoeffForQuat::setup(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   d_quat_diffusion_coarsen_schedule.resize(hierarchy->getNumberOfLevels());

   for (int ln = 0; ln < hierarchy->getNumberOfLevels() - 1; ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
      std::shared_ptr<hier::PatchLevel> finer_level =
          hierarchy->getPatchLevel(ln + 1);

      d_quat_diffusion_coarsen_schedule[ln] =
          d_quat_diffusion_coarsen.createSchedule(level, finer_level);
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
// set diffusion coefficient to 2*H*T*p(phi) as in Pusztai et al.
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

   assert(diffusion->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));
   for (int d = 0; d < NDIM; d++) {
      assert(diffusion->getPointer(d) != nullptr);
   }

   QUATDIFFUSION(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                 ifirst(2), ilast(2),
#endif
                 2. * d_H_parameter, temperature->getPointer(),
                 temperature->getGhostCellWidth()[0], phi->getPointer(),
                 phi->getGhostCellWidth()[0], diffusion->getPointer(0),
                 diffusion->getPointer(1),
#if (NDIM == 3)
                 diffusion->getPointer(2),
#endif
                 0, d_orient_interp_func_type.c_str(), d_avg_func_type.c_str());
}
