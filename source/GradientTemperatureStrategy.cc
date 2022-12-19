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
#include "GradientTemperatureStrategy.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

#include "ConcFort.h"

GradientTemperatureStrategy::GradientTemperatureStrategy(
    const int temperature_id, const int temperature_scratch_id,
    const double temperature0, const double velocity_x,
    std::shared_ptr<tbox::Database> temperature_db)
{
   assert(temperature_id >= 0);
   assert(temperature_scratch_id >= 0);
   assert(temperature0 > 0.);

   d_temperature_id = temperature_id;
   d_temperature_scratch_id = temperature_scratch_id;
   d_temperature0 = temperature0;
   d_velocity_x = velocity_x;
   d_ref_time = 0.;

   d_dtemperaturedt =
       temperature_db->getDoubleWithDefault("dtemperaturedt", 0.0);
   d_target_temperature =
       temperature_db->getDoubleWithDefault("target_temperature", 0.0);

   size_t nterms = temperature_db->getArraySize("gradient");
   assert(nterms == NDIM);

   temperature_db->getDoubleArray("gradient", &d_gradient[0], NDIM);

   if (temperature_db->isDouble("center"))
      temperature_db->getDoubleArray("center", &d_ref_center[0], NDIM);
   else
      for (short i = 0; i < NDIM; i++)
         d_ref_center[i] = 0.;

   for (short i = 0; i < NDIM; i++)
      d_center[i] = d_ref_center[i];

   assert(d_ref_center[0] < 1.e10);
   assert(d_center[0] < 1.e10);
}

double GradientTemperatureStrategy::getCurrentMaxTemperature(
    std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
{
   (void)time;

   math::HierarchyCellDataOpsReal<double> ops(patch_hierarchy);

   return ops.max(d_temperature_id);
}

double GradientTemperatureStrategy::getCurrentMinTemperature(
    std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
{
   (void)time;

   math::HierarchyCellDataOpsReal<double> ops(patch_hierarchy);

   return ops.min(d_temperature_id);
}

double GradientTemperatureStrategy::getCurrentAverageTemperature(
    std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
{
   (void)time;

   math::HierarchyCellDataOpsReal<double> ops(patch_hierarchy);

   return 0.5 * (ops.max(d_temperature_id) + ops.min(d_temperature_id));
}

// set temperature field according to current temperature profile
void GradientTemperatureStrategy::setCurrentTemperature(
    std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
{
   assert(d_temperature_id >= 0);
   assert(d_temperature_scratch_id >= 0);

   // update temperature according to cooling rate
   double temperature0 = d_temperature0 + d_dtemperaturedt * time;

   // from a moving frame perspective, the "center"
   // is moving backward
   d_center[0] = d_ref_center[0] - (time - d_ref_time) * d_velocity_x;
   // tbox::pout << "GradientTemperatureStrategy:" << std::endl;
   // tbox::pout << "d_center[0] = " << d_center[0] << std::endl;
   // tbox::pout << "d_ref_center[0] = " << d_ref_center[0] << std::endl;
   // tbox::pout << "d_velocity_x = " << d_velocity_x << std::endl;

   int maxln = patch_hierarchy->getFinestLevelNumber();

   for (int ln = 0; ln <= maxln; ln++) {

      std::shared_ptr<hier::PatchLevel> level =
          patch_hierarchy->getPatchLevel(ln);

      for (hier::PatchLevel::Iterator p(level->begin()); p != level->end();
           p++) {

         std::shared_ptr<hier::Patch> patch = *p;

         std::shared_ptr<pdat::CellData<double> > t_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_temperature_id)));

         setCurrentTemperature(temperature0, *patch, t_data);

         std::shared_ptr<pdat::CellData<double> > ts_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_temperature_scratch_id)));

         setCurrentTemperature(temperature0, *patch, ts_data);
      }
   }
}

//=======================================================================

void GradientTemperatureStrategy::setCurrentTemperature(
    const double temperature0, hier::Patch& patch,
    std::shared_ptr<pdat::CellData<double> > cd_temp)
{
   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));

   const double* dx = patch_geom->getDx();
   const double* xlo = patch_geom->getXLower();
   const double* xhi = patch_geom->getXUpper();

   hier::IntVector ghost_cells = cd_temp->getGhostCellWidth();
   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   INITGRADIENT(dx, xlo, xhi, ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                ifirst(2), ilast(2),
#endif
                ghost_cells(0), ghost_cells(1),
#if (NDIM == 3)
                ghost_cells(2),
#endif
                cd_temp->getPointer(), &d_center[0], temperature0,
                &d_gradient[0]);
}
