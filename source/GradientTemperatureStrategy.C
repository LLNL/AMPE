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
#include "GradientTemperatureStrategy.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"

#include "ConcFort.h"

GradientTemperatureStrategy::GradientTemperatureStrategy(
    const int temperature_id, const int temperature_scratch_id,
    const double temperature0, std::shared_ptr<tbox::Database> temperature_db)
{
   assert(temperature_id >= 0);
   assert(temperature_scratch_id >= 0);
   assert(temperature0 > 0.);

   d_temperature_id = temperature_id;
   d_temperature_scratch_id = temperature_scratch_id;
   d_temperature0 = temperature0;

   d_dtemperaturedt =
       temperature_db->getDoubleWithDefault("dtemperaturedt", 0.0);
   d_target_temperature =
       temperature_db->getDoubleWithDefault("target_temperature", 0.0);

   size_t nterms = temperature_db->getArraySize("gradient");
   assert(nterms == NDIM);

   temperature_db->getDoubleArray("gradient", &d_gradient[0], NDIM);

   if (temperature_db->isDouble("center"))
      temperature_db->getDoubleArray("center", &d_center[0], NDIM);
   else
      for (short i = 0; i < NDIM; i++)
         d_center[i] = 0.;
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

double GradientTemperatureStrategy::getCurrentTemperature(const double time)
{
   double t = d_temperature0 + d_dtemperaturedt * time;

   // limit temperature to target if set.

   if (d_dtemperaturedt < 0. && t < d_target_temperature) {

      t = d_target_temperature;

   } else if (d_target_temperature > 0.0 && d_dtemperaturedt > 0. &&
              t > d_target_temperature) {

      t = d_target_temperature;
   }

   return t;
}

void GradientTemperatureStrategy::setCurrentTemperature(
    std::shared_ptr<hier::PatchHierarchy> patch_hierarchy, const double time)
{
   assert(d_temperature_id >= 0);
   assert(d_temperature_scratch_id >= 0);

   double temperature0 = d_temperature0 + d_dtemperaturedt * time;
   // tbox::pout<<"GradientTemperatureStrategy: set center T to "
   //          <<temperature0<<std::endl;
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

         setCurrentTemperaturePrivatePatch(temperature0, *patch, t_data);

         std::shared_ptr<pdat::CellData<double> > ts_data(
             SAMRAI_SHARED_PTR_CAST<pdat::CellData<double>, hier::PatchData>(
                 patch->getPatchData(d_temperature_scratch_id)));

         setCurrentTemperaturePrivatePatch(temperature0, *patch, ts_data);
      }
   }
}

//=======================================================================

void GradientTemperatureStrategy::setCurrentTemperaturePrivatePatch(
    const double temperature0, hier::Patch& patch,
    std::shared_ptr<pdat::CellData<double> > cd_temp)
{
   const std::shared_ptr<geom::CartesianPatchGeometry> patch_geom(
       SAMRAI_SHARED_PTR_CAST<geom::CartesianPatchGeometry, hier::PatchGeometry>(
           patch.getPatchGeometry()));

   const double* dx = patch_geom->getDx();
   const double* xlo = patch_geom->getXLower();
   const double* xhi = patch_geom->getXUpper();

   hier::IntVector ghost_cells = cd_temp->getGhostCellWidth();
   const hier::Index ifirst = patch.getBox().lower();
   const hier::Index ilast = patch.getBox().upper();

   FORT_INITGRADIENT(dx, xlo, xhi, ifirst(0), ilast(0), ifirst(1), ilast(1),
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
