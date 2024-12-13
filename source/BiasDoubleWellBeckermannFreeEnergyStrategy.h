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
#ifndef included_BiasDoubleWellBeckermannFreeEnergyStrategy
#define included_BiasDoubleWellBeckermannFreeEnergyStrategy

#include "BiasDoubleWellFreeEnergyStrategy.h"
#include "MeltingTemperatureStrategy.h"

/*!
 * @brief Class BiasDoubleWellFreeEnergyStrategy implements the free energy
 * used by Beckermann et al., J. Comput. Phys. 154 (1999)
 * with
 *    m(T)=alpha*(Te-T)
 *
 */

class BiasDoubleWellBeckermannFreeEnergyStrategy
    : public BiasDoubleWellFreeEnergyStrategy
{
 public:
   BiasDoubleWellBeckermannFreeEnergyStrategy(
       const double alpha, MeltingTemperatureStrategy* meltingTstrat);

   ~BiasDoubleWellBeckermannFreeEnergyStrategy(){};

   void addDrivingForce(const double time, hier::Patch& patch,
                        const int temperature_id, const int phase_id,
                        const int conc_id, const int fl_id, const int fa_id,
                        const int fb_id, const int rhs_id);

   void preRunDiagnostics(const double temperature){};


 private:
   double d_alpha;
   MeltingTemperatureStrategy* d_meltingTstrat;
};

#endif
