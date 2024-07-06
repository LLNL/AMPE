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
#ifndef included_QuadraticEquilibriumPhaseConcentrationsStrategyMultiOrder
#define included_QuadraticEquilibriumPhaseConcentrationsStrategyMultiOrder

#include "PhaseConcentrationsStrategy.h"
#include "QuatModelParameters.h"
#include "QuadraticFreeEnergyFunctionsBinary.h"
#include "EquilibriumPhaseConcentrationsBinaryMultiOrder.h"

#include <string>

class QuadraticEquilibriumPhaseConcentrationsStrategyMultiOrder
    : public EquilibriumPhaseConcentrationsBinaryMultiOrder
{
 public:
   QuadraticEquilibriumPhaseConcentrationsStrategyMultiOrder(
       const int conc_l_id, const int conc_a_id,
       const QuatModelParameters& model_parameters,
       std::shared_ptr<tbox::Database> conc_db);

   ~QuadraticEquilibriumPhaseConcentrationsStrategyMultiOrder() {}

 protected:
   virtual int computePhaseConcentrations(const double t, double* c,
                                          double* hphi, double* x)
   {
      return d_fenergy->computePhaseConcentrations(t, c, hphi, x);
   }

 private:
   std::shared_ptr<Thermo4PFM::QuadraticFreeEnergyFunctionsBinary> d_fenergy;
};

#endif
