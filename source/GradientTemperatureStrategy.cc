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
#include "SAMRAI/tbox/IEEE.h"

#include "ConcFort.h"

GradientTemperatureStrategy::GradientTemperatureStrategy(
    const int temperature_id, const int temperature_scratch_id,
    const double velocity_x, std::shared_ptr<tbox::Database> temperature_db)
{
   assert(temperature_id >= 0);
   assert(temperature_scratch_id >= 0);

   d_temperature_id = temperature_id;
   d_temperature_scratch_id = temperature_scratch_id;
   d_velocity_x = velocity_x;
   d_ref_time = 0.;

   d_temperature0 = temperature_db->getDoubleWithDefault("temperature", -1.);
   d_dtemperaturedt =
       temperature_db->getDoubleWithDefault("dtemperaturedt", 0.0);

   std::string temperature_type = temperature_db->getString("type");
   assert(temperature_type.compare("frozen") == 0);
   d_gradient[0] = tbox::IEEE::getSignalingNaN();

   if (temperature_db->keyExists("filename"))  // from file
   {
      d_temperature_history.reset(new TemperatureHistory(tbox::pout));
      std::string filename = temperature_db->getString("filename");
      tbox::plog << "Read temperature data from file " << filename << std::endl;
      d_temperature_history->readCSV(filename);
   } else {
      assert(temperature_db->getArraySize("gradient") == NDIM);

      temperature_db->getDoubleArray("gradient", &d_gradient[0], NDIM);
      tbox::plog << "GradientTemperatureStrategy with grad T = "
                 << d_gradient[0] << std::endl;
   }

   if (temperature_db->isDouble("center"))
      temperature_db->getDoubleArray("center", &d_ref_center[0], NDIM);
   else
      for (short i = 0; i < NDIM; i++)
         d_ref_center[i] = 0.;

   for (short i = 0; i < NDIM; i++)
      d_center[i] = d_ref_center[i];
   for (short i = 0; i < NDIM; i++)
      d_init_center[i] = d_ref_center[i];

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

   double temperature0;

   if (d_temperature_history)  // from file
   {
      // we track how we move so we can pick T at right position
      d_center[0] = d_ref_center[0] + (time - d_ref_time) * d_velocity_x;
      // update temperature and gradient according to prescribed temperature
      // history
      // tbox::pout << "Read T from file..." << std::endl;
      int ret =
          d_temperature_history->getTandGradT(time, d_center[0], temperature0,
                                              d_gradient[0]);
      d_gradient[1] = 0.;
#if (NDIM == 3)
      d_gradient[2] = 0.;
#endif
      if (ret == 1) TBOX_ERROR("Failed to set T to a frozen gradient");
   } else {
      // from a moving frame perspective, the "center"
      // is moving backward
      d_center[0] = d_ref_center[0] - (time - d_ref_time) * d_velocity_x;
      assert(d_temperature0 > 0.);
      // update temperature according to cooling rate
      temperature0 = d_temperature0 + d_dtemperaturedt * time;
   }
   // tbox::pout << "GradientTemperatureStrategy:" << std::endl;
   // tbox::pout << "d_center[0] = " << d_center[0] << std::endl;
   // tbox::pout << "d_ref_center[0] = " << d_ref_center[0] << std::endl;
   // tbox::pout << "d_velocity_x = " << d_velocity_x << std::endl;
   // tbox::pout << "temperature0 = " << temperature0 << std::endl;
   // tbox::pout << "gradient = " << d_gradient[0] << std::endl;
   assert(!std::isnan(d_gradient[0]));

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
         assert(t_data);

         setCurrentTemperature(temperature0, *patch, t_data);

         std::shared_ptr<pdat::CellData<double> > ts_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_temperature_scratch_id)));
         assert(ts_data);

         setCurrentTemperature(temperature0, *patch, ts_data);
      }
   }

   assert(!std::isnan(getCurrentMaxTemperature(patch_hierarchy, time)));
   assert(!std::isnan(getCurrentMinTemperature(patch_hierarchy, time)));
}

//=======================================================================

void GradientTemperatureStrategy::setCurrentTemperature(
    const double temperature0, hier::Patch& patch,
    std::shared_ptr<pdat::CellData<double> > cd_temp)
{
   assert(!std::isnan(d_gradient[0]));
   assert(!std::isnan(d_gradient[1]));
#if (NDIM == 3)
   assert(!std::isnan(d_gradient[2]));
#endif

   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry,
                              hier::PatchGeometry>(patch.getPatchGeometry()));

   const double* dx = patch_geom->getDx();
   double xlo[NDIM] = {
      patch_geom->getXLower()[0],
      patch_geom->getXLower()[1]
#if (NDIM == 3)
      ,
      patch_geom->getXLower()[2]
#endif
   };

   // when center moves, so does the edge of the domain
   if (d_temperature_history)
      for (int i = 0; i < NDIM; i++)
         xlo[i] += (d_center[i] - d_init_center[i]);

   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   INITGRADIENT(dx, xlo, ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                ifirst(2), ilast(2),
#endif
                cd_temp->getPointer(), cd_temp->getGhostCellWidth()[0],
                &d_center[0], temperature0, &d_gradient[0]);
}
