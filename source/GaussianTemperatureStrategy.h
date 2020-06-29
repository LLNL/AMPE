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
#ifndef included_GaussianTemperatureStrategy
#define included_GaussianTemperatureStrategy

#include "TemperatureStrategy.h"

#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

#include <cassert>

// This strategy fills the temperature field from a Gaussian distribution
// which may depend on time.

class GaussianTemperatureStrategy : public TemperatureStrategy
{
 public:
   GaussianTemperatureStrategy(
       const int temperature_id, const int temperature_scratch_id,
       const int weight_id, std::shared_ptr<tbox::Database> temperature_db,
       std::shared_ptr<geom::CartesianGridGeometry> grid_geometry);

   ~GaussianTemperatureStrategy(){};

   virtual double getCurrentMaxTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time);
   virtual double getCurrentMinTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time);
   virtual double getCurrentAverageTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time);

   virtual void setCurrentTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time);

 private:
   int d_temperature_id;
   int d_temperature_scratch_id;
   int d_weight_id;

   // temperature at tails of Gaussian
   double d_temperature_base;

   double d_dtemperaturedt;

   // target temperature for peak of Gaussian
   double d_target_temperature;

   // initail temperature for peak of Gaussian
   double d_initial_temperature;

   // Gaussian parameters
   double d_center[3];

   // Gaussian center at t=0
   double d_center0[3];
   double d_periodic_length[3];

   // velocity of Gaussian center
   double d_velocity[3];

   double d_standard_dev;

   std::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;
   int d_periodic_flag[3];

   void setCurrentTemperaturePrivatePatch(
       hier::Patch& patch, std::shared_ptr<pdat::CellData<double> > cd_temp,
       const double tgaussian);

   double getCurrentTemperature(const double time);
};

#endif
