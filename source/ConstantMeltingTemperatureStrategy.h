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
#ifndef ConstantMeltingTemperatureStrategy_H
#define ConstantMeltingTemperatureStrategy_H

#include "MeltingTemperatureStrategy.h"

using namespace SAMRAI;

// class to compute melting temperature as a function of concentration
// based on linearized phase-diagram
class ConstantMeltingTemperatureStrategy : public MeltingTemperatureStrategy
{
 public:
   ConstantMeltingTemperatureStrategy(const double Tref,
                                      const int equilibrium_temperature_id)
       : d_Tref(Tref), d_equilibrium_temperature_id(equilibrium_temperature_id)
   {
      assert(d_equilibrium_temperature_id >= 0);
   }

   void evaluate(hier::Patch& patch);

   int equilibrium_temperature_id() { return d_equilibrium_temperature_id; }

 private:
   double d_Tref;
   int d_equilibrium_temperature_id;
};

#endif
