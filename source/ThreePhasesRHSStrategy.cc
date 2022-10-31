// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE

#include "ThreePhasesRHSStrategy.h"
#include "UniformNoise.h"
#include "QuatFort.h"
#include "QuatParams.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"
#include "SAMRAI/math/PatchSideDataOpsReal.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"

ThreePhasesRHSStrategy::ThreePhasesRHSStrategy(
    const QuatModelParameters& model_parameters, const int phase_scratch_id,
    const int conc_scratch_id, const int temperature_scratch_id,
    const int f_l_id, const int f_a_id, const int f_b_id,
    const int phase_mobility_id, const int flux_id,
#ifdef USE_CPODE
    CPODESSolver* sundials_solver,
#else
    CVODESolver* sundials_solver,
#endif
    std::shared_ptr<FreeEnergyStrategy> free_energy_strategy,
    std::shared_ptr<geom::CartesianGridGeometry> grid_geom,
    std::shared_ptr<PhaseFluxStrategy> phase_flux_strategy)
    : d_model_parameters(model_parameters),
      d_phase_well_scale(model_parameters.phase_well_scale()),
      d_phase_scratch_id(phase_scratch_id),
      d_conc_scratch_id(conc_scratch_id),
      d_temperature_scratch_id(temperature_scratch_id),
      d_f_l_id(f_l_id),
      d_f_a_id(f_a_id),
      d_f_b_id(f_b_id),
      d_phase_mobility_id(phase_mobility_id),
      d_flux_id(flux_id),
      d_sundials_solver(sundials_solver),
      d_free_energy_strategy(free_energy_strategy),
      d_grid_geometry(grid_geom),
      d_flux_coarsen_algorithm(tbox::Dimension(NDIM)),
      d_phase_flux_strategy(phase_flux_strategy)
{
   assert(d_f_b_id >= 0);
   assert(d_flux_id >= 0);

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
       tman->getTimer("AMPE::ThreePhasesRHSStrategy::evaluatePhaseRHS()");
}

void ThreePhasesRHSStrategy::setup(
    std::shared_ptr<hier::PatchHierarchy> hierarchy)
{
   d_patch_hierarchy = hierarchy;

   d_flux_coarsen_schedule.resize(hierarchy->getNumberOfLevels());

   for (int ln = 0; ln < hierarchy->getFinestLevelNumber(); ln++) {
      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);
      std::shared_ptr<hier::PatchLevel> finer_level =
          hierarchy->getPatchLevel(ln + 1);

      d_flux_coarsen_schedule[ln] =
          d_flux_coarsen_algorithm.createSchedule(level, finer_level);
   }
}

void ThreePhasesRHSStrategy::evaluateRHS(
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
   double last_time =
       d_sundials_solver->getActualFinalValueOfIndependentVariable();

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

      d_phase_flux_strategy->computeFluxes(level, d_phase_scratch_id, -1,
                                           d_flux_id);

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

void ThreePhasesRHSStrategy::evaluateRHS(const double time,
                                         std::shared_ptr<hier::Patch> patch,
                                         const int ydot_phase_id,
                                         const bool eval_flag)
{
   math::PatchCellDataOpsReal<double> mathops;

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch->getPatchGeometry()));
   const double* dx = patch_geom->getDx();

   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_phase_scratch_id)));
   assert(phase);
   assert(phase->getDepth() == 3);

   std::shared_ptr<pdat::CellData<double> > phase_rhs(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(ydot_phase_id)));
   assert(phase_rhs);
   assert(phase_rhs->getDepth() == 3);

   std::shared_ptr<pdat::SideData<double> > phase_flux(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch->getPatchData(d_flux_id)));
   assert(phase_flux);
   assert(phase_flux->getDepth() == 3);

   std::shared_ptr<pdat::CellData<double> > temperature(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_temperature_scratch_id)));
   assert(temperature);


   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   assert(phase->getDepth() == 3);
   assert(phase->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), NGHOSTS));
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

   // first compute component from interfacial energy
   // for all 3 phases
   COMPUTERHSTHREEPHASES(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                         ifirst(2), ilast(2),
#endif
                         dx, phase_flux->getPointer(0),
                         phase_flux->getPointer(1),
#if (NDIM == 3)
                         phase_flux->getPointer(2),
#endif
                         phase_flux->getGhostCellWidth()[0], d_phase_well_scale,
                         phase->getPointer(), NGHOSTS, phase_rhs->getPointer(),
                         0);

#ifdef DEBUG_CHECK_ASSERTIONS
   double l2rhs = opc.L2Norm(phase_rhs, pbox);
   assert(l2rhs == l2rhs);
#endif

   if (d_free_energy_strategy) {
      d_free_energy_strategy->computeFreeEnergyLiquid(*patch,
                                                      d_temperature_scratch_id,
                                                      d_f_l_id, false);

      d_free_energy_strategy->computeFreeEnergySolidA(*patch,
                                                      d_temperature_scratch_id,
                                                      d_f_a_id, false);

      d_free_energy_strategy->computeFreeEnergySolidB(*patch,
                                                      d_temperature_scratch_id,
                                                      d_f_b_id, false);

      // then add component from chemical energy
      d_free_energy_strategy->addDrivingForce(
          time, *patch, d_temperature_scratch_id, d_phase_scratch_id, -1,
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

   // multiply by mobility
   assert(phase_mobility->getDepth() == phase_rhs->getDepth());
   mathops.multiply(phase_rhs, phase_mobility, phase_rhs, pbox);
}

void ThreePhasesRHSStrategy::projectPhases(const int phase_id,
                                           const int corr_id, const int err_id)
{
   std::shared_ptr<hier::PatchLevel> level =
       d_patch_hierarchy->getPatchLevel(0);
   hier::PatchLevel::Iterator pi(level->begin());
   for (; pi != level->end(); pi++) {
      hier::Patch& patch = **pi;

      std::shared_ptr<pdat::CellData<double> > p_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(phase_id)));
      std::shared_ptr<pdat::CellData<double> > corr_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(corr_id)));
      std::shared_ptr<pdat::CellData<double> > err_data(
          SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
              patch.getPatchData(err_id)));

      const hier::Box& box = patch.getBox();
      const hier::Index& lower = box.lower();
      const hier::Index& upper = box.upper();

      const hier::Box& p_gbox = p_data->getGhostBox();
      const hier::Index& plower = p_gbox.lower();
      const hier::Index& pupper = p_gbox.upper();

      const hier::Box& c_gbox = corr_data->getGhostBox();
      const hier::Index& clower = c_gbox.lower();
      const hier::Index& cupper = c_gbox.upper();

      const hier::Box& e_gbox = err_data->getGhostBox();
      const hier::Index& elower = e_gbox.lower();
      const hier::Index& eupper = e_gbox.upper();

#if NDIM == 2
      PROJECTPHI2D(lower[0], upper[0], lower[1], upper[1], 3,
                   p_data->getPointer(), plower[0], pupper[0], plower[1],
                   pupper[1], corr_data->getPointer(), clower[0], cupper[0],
                   clower[1], cupper[1], err_data->getPointer(), elower[0],
                   eupper[0], elower[1], eupper[1]);
#endif
#if NDIM == 3
      PROJECTPHI3D(lower[0], upper[0], lower[1], upper[1], lower[2], upper[2],
                   3, p_data->getPointer(), plower[0], pupper[0], plower[1],
                   pupper[1], plower[2], pupper[2], corr_data->getPointer(),
                   clower[0], cupper[0], clower[1], cupper[1], clower[2],
                   cupper[2], err_data->getPointer(), elower[0], eupper[0],
                   elower[1], eupper[1], elower[2], eupper[2]);
#endif
   }
}
