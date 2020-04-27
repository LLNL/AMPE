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
#ifndef included_TemperatureStrategyFactory
#define included_TemperatureStrategyFactory

#include "QuatModelParameters.h"
#include "HeatCapacityStrategy.h"

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"

class TemperatureStrategy;
using namespace SAMRAI;

class TemperatureStrategyFactory
{
 public:
   TemperatureStrategyFactory(
       const int temperature_id, const int temperature_scratch_id,
       const int conc_id, const int weight_id, const int temperature_rhs_id,
       const int cp_id, const double molar_volume,
       const bool with_concentration,
       std::shared_ptr<geom::CartesianGridGeometry> grid_geometry,
       HeatCapacityStrategy* heat_capacity_strategy);

   TemperatureStrategy* create(std::shared_ptr<tbox::Database> model_db,
                               std::shared_ptr<tbox::Database> integrator_db,
                               const QuatModelParameters&);

 private:
   const int d_temperature_id;
   const int d_temperature_scratch_id;
   const int d_conc_id;
   const int d_weight_id;
   const int d_temperature_rhs_id;
   const int d_cp_id;

   double d_molar_volume;
   bool d_with_concentration;
   solv::LocationIndexRobinBcCoefs* d_temperature_bc_coefs;
   std::shared_ptr<geom::CartesianGridGeometry> d_grid_geometry;

   TemperatureStrategy* d_strategy;
   HeatCapacityStrategy* d_heat_capacity_strategy;

   double readTemperature0(
       std::shared_ptr<tbox::Database> temperature_db,
       QuatModelParameters::TemperatureType temperature_type);
};

#endif
