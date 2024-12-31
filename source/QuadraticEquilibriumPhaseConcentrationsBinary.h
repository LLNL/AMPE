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
#ifndef included_QuadraticEquilibriumPhaseConcentrationsBinary
#define included_QuadraticEquilibriumPhaseConcentrationsBinary

#include "EquilibriumPhaseConcentrationsBinary.h"
#include "QuadraticFreeEnergyFunctionsBinary.h"

#include <string>

class QuadraticEquilibriumPhaseConcentrationsBinary
    : public EquilibriumPhaseConcentrationsBinary
{
 public:
   QuadraticEquilibriumPhaseConcentrationsBinary(
       const int conc_l_id, const int conc_a_id,
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       const Thermo4PFM::ConcInterpolationType conc_interp_func_type,
       std::shared_ptr<tbox::Database> conc_db);

   ~QuadraticEquilibriumPhaseConcentrationsBinary() {}

 private:
   std::shared_ptr<Thermo4PFM::QuadraticFreeEnergyFunctionsBinary> d_fenergy;

   int computePhaseConcentrations(const double temperature, const double conc,
                                  const double hphi, double* sol) override
   {
      double c = conc;
      double phi = hphi;
      return d_fenergy->computePhaseConcentrations(temperature, &c, &phi, sol);
   }
};

#endif
