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
#ifndef included_PhaseFluxStrategyAnisotropyMultiOrder
#define included_PhaseFluxStrategyAnisotropyMultiOrder

#include "PhaseFluxStrategy.h"

#include <vector>
#include <array>

class PhaseFluxStrategyAnisotropyMultiOrder : public PhaseFluxStrategy
{
 public:
   PhaseFluxStrategyAnisotropyMultiOrder(
       const double epsilon_phase, const double nu, const int knumber,
       std::vector<std::array<double, 4>> quat)
       : d_epsilon_phase(epsilon_phase),
         d_nu(nu),
         d_knumber(knumber),
         d_quat(quat)
   {
      tbox::plog << "PhaseFluxStrategyAnisotropyMultiOrder..." << std::endl;
   }

   void computeFluxes(const std::shared_ptr<hier::PatchLevel> level,
                      const int phase_id, const int quat_id, const int flux_id);

 private:
   const double d_epsilon_phase;
   const double d_nu;
   const int d_knumber;

   std::vector<std::array<double, 4>> d_quat;
};

#endif
