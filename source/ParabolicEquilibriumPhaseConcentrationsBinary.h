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
#ifndef included_ParabolicEquilibriumPhaseConcentrationsBinary
#define included_ParabolicEquilibriumPhaseConcentrationsBinary

#include "EquilibriumPhaseConcentrationsBinary.h"
#include "ParabolicFreeEnergyFunctionsBinary.h"

#include <string>

class ParabolicEquilibriumPhaseConcentrationsBinary
    : public EquilibriumPhaseConcentrationsBinary
{
 public:
   ParabolicEquilibriumPhaseConcentrationsBinary(
       const int conc_l_id, const int conc_a_id,
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       std::shared_ptr<tbox::Database> conc_db);

   ~ParabolicEquilibriumPhaseConcentrationsBinary() {}

 private:
   std::shared_ptr<Thermo4PFM::ParabolicFreeEnergyFunctionsBinary> d_fenergy;

   int computePhaseConcentrations(const double temperature, const double conc,
                                  const double hphi, double* sol) override
   {
      double c = conc;
      double phi[2] = {1. - hphi, hphi};
      return d_fenergy->computePhaseConcentrations(temperature, &c, phi, sol);
   }
};

#endif
