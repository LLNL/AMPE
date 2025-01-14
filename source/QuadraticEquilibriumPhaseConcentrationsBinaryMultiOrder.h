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
#ifndef included_QuadraticEquilibriumPhaseConcentrationsBinaryMultiOrder
#define included_QuadraticEquilibriumPhaseConcentrationsBinaryMultiOrder

#include "PhaseConcentrationsStrategy.h"
#include "QuatModelParameters.h"
#include "QuadraticFreeEnergyFunctionsBinary.h"
#include "EquilibriumPhaseConcentrationsBinaryMultiOrder.h"

#include <string>

class QuadraticEquilibriumPhaseConcentrationsBinaryMultiOrder
    : public EquilibriumPhaseConcentrationsBinaryMultiOrder
{
 public:
   QuadraticEquilibriumPhaseConcentrationsBinaryMultiOrder(
       const int conc_l_id, const int conc_a_id,
       const QuatModelParameters& model_parameters,
       std::shared_ptr<tbox::Database> conc_db);

   ~QuadraticEquilibriumPhaseConcentrationsBinaryMultiOrder() {}

 protected:
   virtual int computePhaseConcentrations(const double t, double* c,
                                          double* hphi, double* x)
   {
      // Thermo4PFM uses solid fraction only for two phases
      double phi = hphi[1];
      return d_fenergy->computePhaseConcentrations(t, c, &phi, x);
   }

 private:
   std::shared_ptr<Thermo4PFM::QuadraticFreeEnergyFunctionsBinary> d_fenergy;
};

#endif
