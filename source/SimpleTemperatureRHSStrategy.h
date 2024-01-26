// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
#ifndef included_SimpleTemperatureRHSStrategy
#define included_SimpleTemperatureRHSStrategy

#include "TemperatureRHSStrategy.h"

class SimpleTemperatureRHSStrategy : public TemperatureRHSStrategy
{
 public:
   SimpleTemperatureRHSStrategy(const double thermal_diffusivity,
                                const double latent_heat,
                                const int temperature_scratch_id,
                                const int cp_id);

   void evaluateRHS(std::shared_ptr<hier::Patch> patch,
                    const int temperature_rhs_id, const int dphidt_id);

 private:
   double d_thermal_diffusivity;
   double d_latent_heat;

   int d_temperature_scratch_id;
   int d_cp_id;
};
#endif
