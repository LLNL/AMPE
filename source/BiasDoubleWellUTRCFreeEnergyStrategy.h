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
#ifndef included_BiasDoubleWellUTRCFreeEnergyStrategy
#define included_BiasDoubleWellUTRCFreeEnergyStrategy

#include "BiasDoubleWellFreeEnergyStrategy.h"
#include "MeltingTemperatureStrategy.h"

/*!
 * @brief Class BiasDoubleWellFreeEnergyStrategy implements the free energy
 * used by Acharya et al., Acta Mat. 124 (2017)
 * with
 *    m(T)=(alpha/pi)tan^-1(gamma(Te-T))
 *
 */

class BiasDoubleWellUTRCFreeEnergyStrategy
    : public BiasDoubleWellFreeEnergyStrategy
{
 public:
   BiasDoubleWellUTRCFreeEnergyStrategy(
       const double alpha, const double gamma,
       MeltingTemperatureStrategy* meltingTstrat);

   ~BiasDoubleWellUTRCFreeEnergyStrategy(){};

   void addDrivingForce(const double time, hier::Patch& patch,
                        const int temperature_id, const int phase_id,
                        const int conc_id, const int fl_id, const int fa_id,
                        const int fb_id, const int rhs_id);

   void preRunDiagnostics(const double temperature){};

 private:
   double d_alpha;
   double d_gamma;
   MeltingTemperatureStrategy* d_meltingTstrat;
};

#endif
