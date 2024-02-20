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
#ifndef included_PhaseFluxStrategySimple
#define included_PhaseFluxStrategySimple

#include "PhaseFluxStrategy.h"

#include "SAMRAI/tbox/PIO.h"

class PhaseFluxStrategySimple : public PhaseFluxStrategy
{
 public:
   PhaseFluxStrategySimple(const double epsilon_phase)
       : d_epsilon_phase(epsilon_phase)
   {
      assert(d_epsilon_phase > 0.);
      tbox::plog << "PhaseFluxStrategySimple with epsilon = " << d_epsilon_phase
                 << std::endl;
   }

   void computeFluxes(const std::shared_ptr<hier::PatchLevel> level,
                      const int phase_id, const int quat_id, const int flux_id);

 private:
   const double d_epsilon_phase;
};

#endif
