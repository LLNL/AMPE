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
#ifndef included_CALPHADequilibriumPhaseConcentrationsStrategyMultiOrder
#define included_CALPHADequilibriumPhaseConcentrationsStrategyMultiOrder

#include "CALPHADFreeEnergyFunctionsBinary.h"
#include "EquilibriumPhaseConcentrationsBinaryMultiOrder.h"

#include <string>

class CALPHADequilibriumPhaseConcentrationsStrategyMultiOrder
    : public EquilibriumPhaseConcentrationsBinaryMultiOrder
{
 public:
   CALPHADequilibriumPhaseConcentrationsStrategyMultiOrder(
       const int conc_l_id, const int conc_a_id,
       std::shared_ptr<tbox::Database> conc_db,
       std::shared_ptr<tbox::Database> newton_db);

   ~CALPHADequilibriumPhaseConcentrationsStrategyMultiOrder() {}

 protected:
   virtual int computePhaseConcentrations(const double t, double* c,
                                          double* hphi, double* x)
   {
      double phi = hphi[1];
      return d_fenergy->computePhaseConcentrations(t, c, &phi, x);
   }

 private:
   std::shared_ptr<Thermo4PFM::CALPHADFreeEnergyFunctionsBinary> d_fenergy;
};

#endif
