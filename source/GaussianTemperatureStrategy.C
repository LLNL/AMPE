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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, UT BATTELLE, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
#include "GaussianTemperatureStrategy.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/tbox/IEEE.h"

#include "ConcFort.h"

GaussianTemperatureStrategy::GaussianTemperatureStrategy(
    const int temperature_id, const int temperature_scratch_id,
    const int weight_id, std::shared_ptr<tbox::Database> temperature_db,
    std::shared_ptr<geom::CartesianGridGeometry> grid_geometry)
    : d_grid_geometry(grid_geometry)
{
   assert(temperature_id >= 0);
   assert(temperature_scratch_id >= 0);

   d_temperature_id = temperature_id;
   d_temperature_scratch_id = temperature_scratch_id;
   d_weight_id = weight_id;

   tbox::plog << "Read Gaussian Temperature profile parameters..." << std::endl;
   d_temperature_base = temperature_db->getDouble("base");
   // default is temperature initially uniform
   d_initial_temperature =
       temperature_db->getDoubleWithDefault("initial", d_temperature_base);
   d_target_temperature = temperature_db->getDouble("target");
   d_dtemperaturedt =
       temperature_db->getDoubleWithDefault("dtemperaturedt", 0.0);
   d_standard_dev = temperature_db->getDouble("standard_deviation");

   if (temperature_db->isDouble("center"))
      temperature_db->getDoubleArray("center", &d_center0[0], NDIM);
   else
      for (short i = 0; i < NDIM; i++)
         d_center0[i] = 0.;

   if (temperature_db->isDouble("velocity"))
      temperature_db->getDoubleArray("velocity", &d_velocity[0], NDIM);
   else
      for (short i = 0; i < NDIM; i++)
         d_velocity[i] = 0.;

   for (short i = 0; i < NDIM; i++)
      d_center[i] = tbox::IEEE::getSignalingNaN();

   tbox::Dimension dim(NDIM);
   const hier::IntVector& one_vec = hier::IntVector::getOne(dim);
   const hier::IntVector periodic = grid_geometry->getPeriodicShift(one_vec);

   for (int dd = 0; dd < NDIM; dd++) {
      d_periodic_flag[dd] = (periodic[dd] != 0);
   }
}

double GaussianTemperatureStrategy::getCurrentMaxTemperature(
    std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
{
   (void)time;

   math::HierarchyCellDataOpsReal<double> mathops(patch_hierarchy);
   return mathops.max(d_temperature_scratch_id);
}

double GaussianTemperatureStrategy::getCurrentMinTemperature(
    std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
{
   (void)time;

   math::HierarchyCellDataOpsReal<double> mathops(patch_hierarchy);
   return mathops.min(d_temperature_scratch_id);
}
double GaussianTemperatureStrategy::getCurrentAverageTemperature(
    std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
{
   (void)time;
   math::HierarchyCellDataOpsReal<double> cellops(patch_hierarchy);

   return cellops.integral(d_temperature_id, d_weight_id) /
          cellops.sumControlVolumes(d_temperature_id, d_weight_id);
}
double GaussianTemperatureStrategy::getCurrentTemperature(const double time)
{
   double t_peak = d_initial_temperature + d_dtemperaturedt * time;

   // saturate temperature to target value if reached.

   if (d_dtemperaturedt < 0. && t_peak < d_target_temperature) {

      t_peak = d_target_temperature;

   } else if (d_dtemperaturedt > 0. && t_peak > d_target_temperature) {

      t_peak = d_target_temperature;
   }

   return t_peak;
}

void GaussianTemperatureStrategy::setCurrentTemperature(
    std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
{
   const double* low = d_grid_geometry->getXLower();
   const double* up = d_grid_geometry->getXUpper();

   for (short d = 0; d < NDIM; d++) {
      d_center[d] = d_center0[d] + time * d_velocity[d];
      if (d_periodic_flag[d]) {
         double ll = (up[d] - low[d]);
         d_periodic_length[d] = ll;
         while (d_center[d] < low[d])
            d_center[d] += ll;
         while (d_center[d] > up[d])
            d_center[d] -= ll;
      } else {
         d_periodic_length[d] = -1.;
      }
   }

   if (d_dtemperaturedt != 0.0 || time == 0.0) {

      double tgaussian = getCurrentTemperature(time);

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
            setCurrentTemperaturePrivatePatch(*patch, t_data, tgaussian);

            std::shared_ptr<pdat::CellData<double> > ts_data(
                SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                    patch->getPatchData(d_temperature_scratch_id)));
            setCurrentTemperaturePrivatePatch(*patch, ts_data, tgaussian);
         }
      }
   }

   // double t = getCurrentMinTemperature( patch_hierarchy, time );
   // tbox::pout << "setCurrentTemperature():  Min. Temperature = " << t <<
   // std::endl; t = getCurrentMaxTemperature( patch_hierarchy, time );
   // tbox::pout << "setCurrentTemperature():  Max. Temperature = " << t <<
   // std::endl;
}

//=======================================================================

void GaussianTemperatureStrategy::setCurrentTemperaturePrivatePatch(
    hier::Patch& patch, std::shared_ptr<pdat::CellData<double> > cd_temp,
    const double tgaussian)
{
   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
           patch.getPatchGeometry()));

   const double* dx = patch_geom->getDx();
   const double* xlo = patch_geom->getXLower();
   const double* xhi = patch_geom->getXUpper();

   const hier::Index ifirst(patch.getBox().lower());
   const hier::Index ilast(patch.getBox().upper());

   FORT_INITGAUSSIAN(dx, xlo, xhi, ifirst(0), ilast(0), ifirst(1), ilast(1),
#if (NDIM == 3)
                     ifirst(2), ilast(2),
#endif
                     cd_temp->getPointer(), cd_temp->getGhostCellWidth()[0],
                     d_center, d_periodic_length, d_standard_dev,
                     d_temperature_base, tgaussian);
}
