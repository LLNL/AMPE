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
#ifndef included_ParabolicEquilibriumThreePhasesBinaryMultiOrder
#define included_ParabolicEquilibriumThreePhasesBinaryMultiOrder

#include "EquilibriumPhaseConcentrationsBinaryMultiOrderThreePhases.h"
#include "ParabolicFreeEnergyFunctionsBinaryThreePhase.h"

class ParabolicEquilibriumThreePhasesBinaryMultiOrder
    : public EquilibriumPhaseConcentrationsBinaryMultiOrderThreePhases
{
 public:
   ParabolicEquilibriumThreePhasesBinaryMultiOrder(
       const short norderp_A, const int conc_l_id, const int conc_a_id,
       const int conc_b_id, const QuatModelParameters& model_parameters,
       std::shared_ptr<tbox::Database> conc_db);

   ~ParabolicEquilibriumThreePhasesBinaryMultiOrder();

 protected:
   virtual int computePhaseConcentrations(const double t, double* c,
                                          double* hphi, double* x)
   {
      return d_fenergy->computePhaseConcentrations(t, c, hphi, x);
   }

 private:
   std::shared_ptr<Thermo4PFM::ParabolicFreeEnergyFunctionsBinaryThreePhase>
       d_fenergy;
};

#endif
