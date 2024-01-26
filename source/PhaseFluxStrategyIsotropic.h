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
#ifndef included_PhaseFluxStrategyIsotropic
#define included_PhaseFluxStrategyIsotropic

#include "PhaseFluxStrategy.h"

class PhaseFluxStrategyIsotropic : public PhaseFluxStrategy
{
 public:
   PhaseFluxStrategyIsotropic(const double epsilon_phase)
       : d_epsilon_phase(epsilon_phase)
   {
      tbox::plog << "Uses PhaseFluxStrategyIsotropic class..." << std::endl;
   }

   void computeFluxes(const std::shared_ptr<hier::PatchLevel> level,
                      const int phase_id, const int quat_id, const int flux_id);

 private:
   const double d_epsilon_phase;
};

#endif
