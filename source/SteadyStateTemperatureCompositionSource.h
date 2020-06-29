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
#ifndef included_SteadyStateTemperatureCompositionSource
#define included_SteadyStateTemperatureCompositionSource

#include "SteadyStateTemperatureStrategy.h"
#include "HeatCapacityStrategy.h"

#include "SAMRAI/solv/LocationIndexRobinBcCoefs.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/tbox/Database.h"

#include <cassert>
#include <vector>

// This strategy fills the temperature field from teh solution of
// diffusion equation
class TemperatureFACSolver;

class SteadyStateTemperatureCompositionSource
    : public SteadyStateTemperatureStrategy
{
 public:
   SteadyStateTemperatureCompositionSource(
       const int temperature_scratch_id,
       const int composition_id,  // source depends on composition
       const int rhs_id,  // used internally only, but allocated outside class
       const int weight_id, const double thermal_diffusivity, const int cp_id,
       const std::vector<double>& T_source,
       std::shared_ptr<tbox::Database> temperature_sys_solver_database,
       HeatCapacityStrategy* heat_capacity_strategy,
       solv::LocationIndexRobinBcCoefs* bc_coefs);

   ~SteadyStateTemperatureCompositionSource(){};

   void setCurrentTemperature(
       std::shared_ptr<hier::PatchHierarchy> patch_hierarchy,
       const double time);

 private:
   const double d_thermal_diffusivity;
   const std::vector<double> d_T_source;

   int d_composition_id;
   int d_cp_id;

   HeatCapacityStrategy* d_heat_capacity_strategy;
};


#endif
