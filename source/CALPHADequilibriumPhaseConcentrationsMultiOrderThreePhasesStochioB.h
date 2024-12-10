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
#ifndef included_CALPHADequilibriumPhaseConcentrationsMultiOrderThreePhasesStochioB
#define included_CALPHADequilibriumPhaseConcentrationsMultiOrderThreePhasesStochioB

#include "EquilibriumPhaseConcentrationsBinaryMultiOrderThreePhases.h"
#include "InterpolationType.h"
#include "QuatModelParameters.h"

// Thermo4PFM
#include "CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB.h"

#include "SAMRAI/tbox/InputManager.h"

class CALPHADequilibriumPhaseConcentrationsMultiOrderThreePhasesStochioB
    : public EquilibriumPhaseConcentrationsBinaryMultiOrderThreePhases
{
 public:
   CALPHADequilibriumPhaseConcentrationsMultiOrderThreePhasesStochioB(
       const short norderp_A, const int conc_l_id, const int conc_a_id,
       const int conc_b_id, const QuatModelParameters& model_parameters,
       std::shared_ptr<tbox::Database> conc_db,
       std::shared_ptr<tbox::Database> newton_db);

   ~CALPHADequilibriumPhaseConcentrationsMultiOrderThreePhasesStochioB() {}

 private:
   QuatModelParameters d_model_parameters;

   std::unique_ptr<
       Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB>
       d_calphad_fenergy;

   /*
    * Implement specific method to compute auxilliary compositions when
    * phase is nearly 100% the sochiometric phase B
    */
   int computePhaseConcentrations(const double temp, double* c, double* hphi,
                                  double* x) override;
};

#endif
