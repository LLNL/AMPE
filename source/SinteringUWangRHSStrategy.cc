// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE

#include "SinteringUWangRHSStrategy.h"
#include "QuatFort.h"
#include "ArrayOperation.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"
#include "SAMRAI/math/PatchSideDataOpsReal.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/ArrayDataOperationUtilities.h"

SinteringUWangRHSStrategy::SinteringUWangRHSStrategy(
    const QuatModelParameters& model_parameters, const int phase_id,
    const int conc_id, const int temperature_id, const int phase_mobility_id,
    const int flux_id, QuatIntegrator* time_integrator,
    std::shared_ptr<geom::CartesianGridGeometry> grid_geom,
    std::shared_ptr<PhaseFluxStrategy> phase_flux_strategy)
    : d_model_parameters(model_parameters),
      d_B(model_parameters.WangSintering_B()),
      d_phase_id(phase_id),
      d_conc_id(conc_id),
      d_temperature_id(temperature_id),
      d_phase_mobility_id(phase_mobility_id),
      d_flux_id(flux_id),
      d_time_integrator(time_integrator),
      d_grid_geometry(grid_geom),
      d_flux_coarsen_algorithm(tbox::Dimension(NDIM)),
      d_phase_flux_strategy(phase_flux_strategy)
{
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
       tman->getTimer("AMPE::SinteringUWangRHSStrategy::evaluatePhaseRHS()");
}

void SinteringUWangRHSStrategy::setup(
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

void SinteringUWangRHSStrategy::evaluateRHS(
    const double time, std::shared_ptr<hier::PatchHierarchy> hierarchy,
    const int ydot_phase_id, const bool eval_flag)
{
   assert(ydot_phase_id >= 0);
   assert(d_phase_id >= 0);
   assert(d_temperature_id >= 0);
   assert(d_phase_flux_strategy);

   t_phase_rhs_timer->start();

   static double old_time = -1.;

   math::HierarchyCellDataOpsReal<double> cellops(hierarchy);
#ifndef NDEBUG
   const double norm_y = cellops.L2Norm(d_phase_id);
   assert(!std::isnan(norm_y));
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

      d_phase_flux_strategy->computeFluxes(level, d_phase_id, -1, d_flux_id);

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

void SinteringUWangRHSStrategy::evaluateRHS(const double time,
                                            std::shared_ptr<hier::Patch> patch,
                                            const int ydot_phase_id,
                                            const bool eval_flag)
{
   // tbox::plog << "SinteringUWangRHSStrategy::evaluateRHS()..." << std::endl;
   math::PatchCellDataOpsReal<double> mathops;

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch->getPatchGeometry()));
   const double* dx = patch_geom->getDx();

   std::shared_ptr<pdat::CellData<double> > phase(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_phase_id)));
   assert(phase);
   const int norder = phase->getDepth();
   assert(norder > 1);

   std::shared_ptr<pdat::CellData<double> > conc(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_conc_id)));
   assert(conc);

   std::shared_ptr<pdat::CellData<double> > phase_rhs(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(ydot_phase_id)));
   assert(phase_rhs);
   assert(phase_rhs->getDepth() == norder);

   std::shared_ptr<pdat::SideData<double> > phase_flux(
       SAMRAI_SHARED_PTR_CAST<pdat::SideData<double>, hier::PatchData>(
           patch->getPatchData(d_flux_id)));
   assert(phase_flux);
   assert(phase_flux->getDepth() == norder);

   const hier::Box& pbox = patch->getBox();
   const hier::Index& ifirst = pbox.lower();
   const hier::Index& ilast = pbox.upper();

   assert(phase_rhs->getGhostCellWidth() ==
          hier::IntVector(tbox::Dimension(NDIM), 0));
#ifndef NDEBUG
   SAMRAI::math::PatchSideDataNormOpsReal<double> ops;
   double l2f = ops.L2Norm(phase_flux, pbox);
   assert(!std::isnan(l2f));
#endif

   std::vector<double> tmp(pbox.size());

   assert(phase_flux->getDepth() == phase_flux->getDepth());
   assert(d_B > 0.);
   COMPUTERHS_WANG_SINTERING(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                             ifirst(2), ilast(2),
#endif
                             norder, dx, phase_flux->getPointer(0),
                             phase_flux->getPointer(1),
#if (NDIM == 3)
                             phase_flux->getPointer(2),
#endif
                             phase_flux->getGhostCellWidth()[0], d_B,
                             phase->getPointer(), phase->getGhostCellWidth()[0],
                             conc->getPointer(), conc->getGhostCellWidth()[0],
                             tmp.data(), phase_rhs->getPointer(),
                             phase_rhs->getGhostCellWidth()[0]);

#ifndef NDEBUG
   SAMRAI::math::PatchCellDataNormOpsReal<double> opc;
   double l2rhs = opc.L2Norm(phase_rhs, pbox);
   assert(!std::isnan(l2rhs));
#endif

#ifndef NDEBUG
   l2rhs = opc.L2Norm(phase_rhs, pbox);
   assert(!std::isnan(l2rhs));
#endif

   std::shared_ptr<pdat::CellData<double> > phase_mobility(
       SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
           patch->getPatchData(d_phase_mobility_id)));
   assert(phase_mobility);
#ifndef NDEBUG
   l2rhs = opc.L1Norm(phase_mobility, pbox);
   assert(!std::isnan(l2rhs));
#endif

   // multiply by mobility
   assert(phase_mobility->getDepth() == 1);
   const hier::IntVector src_shift(pbox.getDim(), 0);
   MultiplyOperation<double> multop;

   const unsigned int src_depth = 0;
   const unsigned int num_depth = 1;
   for (int dst_depth = 0; dst_depth < norder; dst_depth++) {
      pdat::ArrayDataOperationUtilities<double, MultiplyOperation<double> >::
          doArrayDataOperationOnBox(phase_rhs->getArrayData(),
                                    phase_mobility->getArrayData(), pbox,
                                    src_shift, dst_depth, src_depth, num_depth,
                                    multop);
   }
}
