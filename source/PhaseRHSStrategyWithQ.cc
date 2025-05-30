// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE

#include "PhaseRHSStrategyWithQ.h"
#include "UniformNoise.h"
#include "QuatFort.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"
#include "SAMRAI/math/PatchSideDataOpsReal.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

PhaseRHSStrategyWithQ::PhaseRHSStrategyWithQ(
    const QuatModelParameters& model_parameters, const int phase_scratch_id,
    const int conc_scratch_id, const int quat_scratch_id,
    const int temperature_scratch_id, const int eta_scratch_id,
    const int f_l_id, const int f_a_id, const int f_b_id,
    const int phase_mobility_id, const int flux_id,
    const int quat_grad_modulus_id, const int noise_id,
    const int phase_rhs_visit_id, QuatIntegrator* time_integrator,
    std::shared_ptr<FreeEnergyStrategy> free_energy_strategy,
    std::shared_ptr<geom::CartesianGridGeometry> grid_geom,
    std::shared_ptr<PhaseFluxStrategy> phase_flux_strategy)
    : d_model_parameters(model_parameters),
      d_phase_well_scale(model_parameters.phase_well_scale()),
      d_eta_well_scale(model_parameters.eta_well_scale()),
      d_H_parameter(model_parameters.H_parameter()),
      d_energy_interp_func_type(model_parameters.energy_interp_func_type()),
      d_orient_interp_func_type1(model_parameters.orient_interp_func_type1()),
      d_orient_interp_func_type2(model_parameters.orient_interp_func_type2()),
      d_phase_scratch_id(phase_scratch_id),
      d_conc_scratch_id(conc_scratch_id),
      d_quat_scratch_id(quat_scratch_id),
      d_temperature_scratch_id(temperature_scratch_id),
      d_eta_scratch_id(eta_scratch_id),
      d_f_l_id(f_l_id),
      d_f_a_id(f_a_id),
      d_f_b_id(f_b_id),
      d_phase_mobility_id(phase_mobility_id),
      d_flux_id(flux_id),
      d_quat_grad_modulus_id(quat_grad_modulus_id),
      d_noise_id(noise_id),
      d_phase_rhs_visit_id(phase_rhs_visit_id),
      d_time_integrator(time_integrator),
      d_free_energy_strategy(free_energy_strategy),
      d_grid_geometry(grid_geom),
      d_flux_coarsen_algorithm(tbox::Dimension(NDIM)),
      d_phase_flux_strategy(phase_flux_strategy)
{
   hier::VariableDatabase* vdb = hier::VariableDatabase::getDatabase();

   std::shared_ptr<hier::Variable> variable;
   vdb->mapIndexToVariable(d_flux_id, variable);
   d_flux_coarsen_op = d_grid_geometry->lookupCoarsenOperator(variable,
                                                              "CONSERVATIVE_"
                                                              "COARSEN");

   d_flux_coarsen_algorithm.registerCoarsen(d_flux_id, d_flux_id,
                                            d_flux_coarsen_op);
   d_deltat = 1.e9;

   tbox::TimerManager* tman = tbox::TimerManager::getManager();
   t_phase_rhs_timer =
       tman->getTimer("AMPE::PhaseRHSStrategyWithQ::evaluatePhaseRHS()");
}

void PhaseRHSStrategyWithQ::setup(
    std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   d_flux_coarsen_schedule.resize(hierarchy->getNumberOfLevels());

   for (int ln = 0; ln < hierarchy->getFinestLevelNumber(); ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
      std::shared_ptr<hier::PatchLevel> finer_level =
          hierarchy->getPatchLevel(ln + 1);

      d_flux_coarsen_schedule[ln] =
          d_flux_coarsen_algorithm.createSchedule(level, finer_level);
   }
}

void PhaseRHSStrategyWithQ::evaluateRHS(
    const double time, std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int ydot_phase_id, const bool eval_flag)
{
   assert(ydot_phase_id >= 0);
   assert(d_phase_scratch_id >= 0);
   assert(d_temperature_scratch_id >= 0);
   assert(d_phase_mobility_id >= 0);
   assert(d_phase_flux_strategy);

   t_phase_rhs_timer->start();

   static double old_time = -1.;

   math::HierarchyCellDataOpsReal<double> cellops(hierarchy);
#ifdef DEBUG_CHECK_ASSERTIONS
   const double norm_y = cellops.L2Norm(d_phase_scratch_id);
   assert(norm_y == norm_y);
#endif

   // get time of last accepted step
   double last_time = d_time_integrator->lastAcceptedTime();

   // if this is not a FD operation and the time has been updated
   // turn flag ON to recompute random noise
   d_newtime = false;
   if (time != old_time && eval_flag) {
      d_deltat = time - last_time;
      // tbox::pout<<"deltat="<<d_deltat<<endl;
      d_newtime = true;
   }

   // Loop from finest coarsest levels.  We assume that ghost cells
   // on all levels have already been filled by a prior call to
   // setCoefficients (via fillScratch).
   for (int ln = hierarchy->getFinestLevelNumber(); ln >= 0; --ln) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      d_phase_flux_strategy->computeFluxes(level, d_phase_scratch_id,
                                           d_quat_scratch_id, d_flux_id);

      // Coarsen flux data from next finer level so that
      // the computed flux becomes the composite grid flux.
      if (ln < hierarchy->getFinestLevelNumber()) {
         d_flux_coarsen_schedule[ln]->coarsenData();
      }

      for (hier::PatchLevel::Iterator ip(level->begin()); ip != level->end();
           ++ip) {
         std::shared_ptr<hier::Patch> patch = *ip;

         evaluateRHS(time, patch, ydot_phase_id, eval_flag);
      }
   }

   old_time = time;

   t_phase_rhs_timer->stop();
};

void PhaseRHSStrategyWithQ::evaluateRHS(const double time,
                                        std::shared_ptr<hier::Patch> patch,
                                        const int ydot_phase_id,
                                        const bool eval_flag)
{
   math::PatchCellDataOpsReal<double> mathops;

   if (d_free_energy_strategy) {
      d_free_energy_strategy->computeFreeEnergyLiquid(*patch,
                                                      d_temperature_scratch_id,
                                                      d_f_l_id, false);

      d_free_energy_strategy->computeFreeEnergySolidA(*patch,
                                                      d_temperature_scratch_id,
                                                      d_f_a_id, false);
   }

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch->getPatchGeometry()));
   const double* dx = patch_geom->getDx();

   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_phase_scratch_id)));
   assert(phase);

   std::shared_ptr<pdat::CellData<double> > phase_rhs(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(ydot_phase_id)));
   assert(phase_rhs);

   if (d_free_energy_strategy) {
      std::shared_ptr<pdat::CellData<double> > fl(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_f_l_id)));
      assert(fl);
      std::shared_ptr<pdat::CellData<double> > fa(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_f_a_id)));
      assert(fa);
   }

   std::shared_ptr<pdat::SideData<double> > phase_flux(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch->getPatchData(d_flux_id)));
   assert(phase_flux);
   assert(phase_flux->getDepth() == 1);

   std::shared_ptr<pdat::CellData<double> > temperature(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_temperature_scratch_id)));
   assert(temperature);

   int with_orient = 0;
   double* ptr_quat_grad_modulus = nullptr;
   if (d_model_parameters.evolveQuat()) {
      with_orient = 1;
      assert(d_quat_grad_modulus_id >= 0);
      std::shared_ptr<pdat::CellData<double> > qgm(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_quat_grad_modulus_id)));
      ptr_quat_grad_modulus = qgm->getPointer();
   }

   int with_third_phase = 0;
   double* ptr_eta = nullptr;
   int ngeta = 0;
   if (d_eta_scratch_id >= 0) {
      with_third_phase = 1;
      std::shared_ptr<pdat::CellData<double> > eta(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_eta_scratch_id)));
      ptr_eta = eta->getPointer();
      ngeta = eta->getGhostCellWidth()[0];
   }

   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   assert(phase_rhs->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));
#ifdef DEBUG_CHECK_ASSERTIONS
   SAMRAI::math::PatchCellDataNormOpsReal<double> opc;
   double l2t = opc.L2Norm(temperature, pbox);
   assert(l2t == l2t);

   SAMRAI::math::PatchSideDataNormOpsReal<double> ops;
   double l2f = ops.L2Norm(phase_flux, pbox);
   assert(l2f == l2f);
#endif
   char well_func_type = 'd';
   const char interpf = energyInterpChar(d_energy_interp_func_type);

   // first compute component from interfacial energy
   COMPUTERHSPBG(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                 ifirst(2), ilast(2),
#endif
                 dx, 2.0 * d_H_parameter, d_model_parameters.epsilon_q(),
                 phase_flux->getPointer(0), phase_flux->getPointer(1),
#if (NDIM == 3)
                 phase_flux->getPointer(2),
#endif
                 phase_flux->getGhostCellWidth()[0], temperature->getPointer(),
                 temperature->getGhostCellWidth()[0], d_phase_well_scale,
                 d_eta_well_scale, phase->getPointer(),
                 phase->getGhostCellWidth()[0], ptr_eta, ngeta,
                 ptr_quat_grad_modulus, 0, phase_rhs->getPointer(), 0,
                 &well_func_type, &well_func_type, &interpf,
                 d_orient_interp_func_type1.c_str(),
                 d_orient_interp_func_type2.c_str(), with_orient,
                 with_third_phase);

#ifdef DEBUG_CHECK_ASSERTIONS
   double l2rhs = opc.L2Norm(phase_rhs, pbox);
   assert(l2rhs == l2rhs);
#endif

   // then add component from chemical energy
   if (d_free_energy_strategy) {
      d_free_energy_strategy->addDrivingForce(
          time, *patch, d_temperature_scratch_id, d_phase_scratch_id,
          d_conc_scratch_id, d_f_l_id, d_f_a_id, d_f_b_id, ydot_phase_id);
   }

#ifdef DEBUG_CHECK_ASSERTIONS
   l2rhs = opc.L2Norm(phase_rhs, pbox);
   assert(l2rhs == l2rhs);
#endif

   std::shared_ptr<pdat::CellData<double> > phase_mobility(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_phase_mobility_id)));
   assert(phase_mobility);

   if (d_model_parameters.with_rhs_visit_output() && eval_flag) {
      std::shared_ptr<pdat::CellData<double> > phase_rhs_visit(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_phase_rhs_visit_id)));
      assert(phase_rhs_visit);
      mathops.copyData(phase_rhs_visit, phase_rhs, pbox);
   }

   // multiply by mobility
   mathops.multiply(phase_rhs, phase_mobility, phase_rhs, pbox);

   // add noise
   if (d_model_parameters.noise_amplitude() > 0. && d_deltat > 0.) {
      const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
      UniformNoise& noise(*(UniformNoise::instance(mpi.getRank())));
      if (d_newtime) {
         noise.setField(patch, d_noise_id, d_phase_scratch_id);
      }
      std::shared_ptr<pdat::CellData<double> > noise_field(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch->getPatchData(d_noise_id)));
      double alpha = d_model_parameters.noise_amplitude() / sqrt(d_deltat);
      mathops.axpy(phase_rhs, alpha, noise_field, phase_rhs, patch->getBox());
   }
}
