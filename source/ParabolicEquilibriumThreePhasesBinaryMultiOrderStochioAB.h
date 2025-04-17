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
#ifndef included_ParabolicEquilibriumThreePhasesBinaryMultiOrderStochioAB
#define included_ParabolicEquilibriumThreePhasesBinaryMultiOrderStochioAB

#include "EquilibriumPhaseConcentrationsBinaryMultiOrderThreePhases.h"
#include "QuatModelParameters.h"

class ParabolicEquilibriumThreePhasesBinaryMultiOrderStochioAB
    : public EquilibriumPhaseConcentrationsBinaryMultiOrderThreePhases
{
 public:
   ParabolicEquilibriumThreePhasesBinaryMultiOrderStochioAB(
       const short norderp_A, const int conc_l_id, const int conc_a_id,
       const int conc_b_id, const QuatModelParameters& model_parameters,
       std::shared_ptr<tbox::Database> conc_db);

   ~ParabolicEquilibriumThreePhasesBinaryMultiOrderStochioAB();

 protected:
   virtual int computePhaseConcentrations(const double t, double* c,
                                          double* hphi, double* x);

 private:
   const QuatModelParameters d_model_parameters;
};

#endif
