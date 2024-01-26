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
#ifndef LinearMeltingTemperatureStrategy_H
#define LinearMeltingTemperatureStrategy_H

#include "MeltingTemperatureStrategy.h"

using namespace SAMRAI;

// class to compute melting temperature as a function of concentration
// based on linearized phase-diagram
class LinearMeltingTemperatureStrategy : public MeltingTemperatureStrategy
{
 public:
   LinearMeltingTemperatureStrategy(const double Tref, const double c0,
                                    const double liquidus_slope,
                                    const int concentration_id,
                                    const int equilibrium_temperature_id);

   void evaluate(hier::Patch& patch);

   int equilibrium_temperature_id() { return d_equilibrium_temperature_id; }

 private:
   double d_Tref;

   // nominal composition
   double d_c0;
   double d_liquidus_slope;

   // concentration variable to use in linear model
   int d_concentration_id;

   int d_equilibrium_temperature_id;
};

#endif
