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
//
#include "SteadyStateTemperatureCompositionSource.h"

#include "TemperatureFACSolver.h"
#include "QuatFort.h"
#include "TemperatureFACOps.h"

#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/pdat/CellData.h"

SteadyStateTemperatureCompositionSource::
    SteadyStateTemperatureCompositionSource(
        const int temperature_scratch_id, const int composition_id,
        const int rhs_id, const int weight_id, const double thermal_diffusivity,
        const int cp_id, const std::vector<double>& T_source,
        std::shared_ptr<tbox::Database> temperature_sys_solver_database,
        HeatCapacityStrategy* heat_capacity_strategy,
        solv::LocationIndexRobinBcCoefs* bc_coefs)
    : SteadyStateTemperatureStrategy(temperature_scratch_id, rhs_id, weight_id,
                                     temperature_sys_solver_database, bc_coefs),
      d_thermal_diffusivity(thermal_diffusivity),
      d_T_source(T_source),
      d_heat_capacity_strategy(heat_capacity_strategy)
{
   assert(temperature_scratch_id >= 0);
   assert(thermal_diffusivity > 0.);
   assert(thermal_diffusivity < 1.e15);
   if (T_source.size() > 0) assert(composition_id >= 0);

   d_composition_id = composition_id;
   d_cp_id = cp_id;
}

void SteadyStateTemperatureCompositionSource::setCurrentTemperature(
    std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
{
   (void)time;

   // tbox::pout<<"SteadyStateTemperatureCompositionSource: solve for steady
   // state T"<<endl;

   assert(d_temperature_sys_solver);
   assert(d_heat_capacity_strategy);

   d_temperature_sys_solver->setOperatorCoefficients(1., 0.,
                                                     -1. *
                                                         d_thermal_diffusivity);

   d_heat_capacity_strategy->setCurrentValue(patch_hierarchy);

   int maxln = patch_hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {

      std::shared_ptr<hier::PatchLevel> level =
          patch_hierarchy->getPatchLevel(ln);
      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           p++) {
         std::shared_ptr<hier::Patch> patch = *p;

         const hier::Box& pbox = patch->getBox();
         const hier::Index& ifirst = pbox.lower();
         const hier::Index& ilast = pbox.upper();

         std::shared_ptr<pdat::CellData<double> > rhs(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_rhs_id)));
         std::shared_ptr<pdat::CellData<double> > conc(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_composition_id)));
         std::shared_ptr<pdat::CellData<double> > cp(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_cp_id)));

         SOURCE_TEMPERATURE(ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                            ifirst(2), ilast(2),
#endif
                            conc->getPointer(), conc->getGhostCellWidth()[0],
                            rhs->getPointer(), rhs->getGhostCellWidth()[0],
                            cp->getPointer(), cp->getGhostCellWidth()[0],
                            &d_T_source[0],
                            static_cast<int>(d_T_source.size()));
      }
   }

   this->solveSystem();

   // math::HierarchyCellDataOpsReal<double> mathops( patch_hierarchy );
   // tbox::pout<<"max. T after solve heat
   // equation="<<mathops.max(d_temperature_scratch_id)<<endl; tbox::pout<<"min.
   // T after solve heat
   // equation="<<mathops.min(d_temperature_scratch_id)<<endl;
}
