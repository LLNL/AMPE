// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#include "EnergyEvaluationStrategy.h"
#include "tools.h"

#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

void EnergyEvaluationStrategy::evaluateEnergy(
    const std::shared_ptr<hier::PatchHierarchy> hierarchy, const double time,
    double& total_energy, double& total_phase_e, double& total_orient_e,
    double& total_qint_e, double& total_well_e, double& total_bulk_e,
    const bool gp)
{
   total_energy = 0.;
   total_phase_e = 0.;
   total_orient_e = 0.;
   total_qint_e = 0.;
   total_bulk_e = 0.;
   total_well_e = 0.;

   const int maxln = hierarchy->getFinestLevelNumber();
   for (int ln = 0; ln <= maxln; ln++) {

      std::shared_ptr<hier::PatchLevel> level = hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           ++p) {

         std::shared_ptr<hier::Patch> patch = *p;

         evaluateEnergy(patch, time, total_energy, total_phase_e,
                        total_orient_e, total_qint_e, total_well_e,
                        total_bulk_e, gp);
      }
   }

   total_energy = sumReduction(total_energy);
   total_phase_e = sumReduction(total_phase_e);
   total_orient_e = sumReduction(total_orient_e);
   total_qint_e = sumReduction(total_qint_e);
   total_well_e = sumReduction(total_well_e);
   total_bulk_e = sumReduction(total_bulk_e);

   math::HierarchyCellDataOpsReal<double> mathops(hierarchy);

   if (d_model_parameters.with_visit_energy_output()) {
      double emin = mathops.min(d_energy_diag_id);
      double emax = mathops.max(d_energy_diag_id);
      tbox::plog << "Min. energy density = " << emin << std::endl;
      tbox::plog << "Max. energy density = " << emax << std::endl;
   }
}
