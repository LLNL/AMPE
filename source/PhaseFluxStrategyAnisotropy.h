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
#ifndef included_PhaseFluxStrategyAnisotropy
#define included_PhaseFluxStrategyAnisotropy

#include "PhaseFluxStrategy.h"

class PhaseFluxStrategyAnisotropy : public PhaseFluxStrategy
{
 public:
   PhaseFluxStrategyAnisotropy(const double epsilon_phase, const double nu,
                               const int knumber)
       : d_epsilon_phase(epsilon_phase), d_nu(nu), d_knumber(knumber)
   {
   }

   void computeFluxes(const std::shared_ptr<hier::PatchLevel> level,
                      const int phase_id, const int quat_id, const int flux_id);

 private:
   const double d_epsilon_phase;
   const double d_nu;
   const int d_knumber;
};

#endif
