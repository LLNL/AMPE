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
#ifndef included_CALPHADequilibriumPhaseConcentrationsThreePhasesStochioB
#define included_CALPHADequilibriumPhaseConcentrationsThreePhasesStochioB

#include "EquilibriumPhaseConcentrationsThreePhases.h"
#include "InterpolationType.h"
#include "QuatModelParameters.h"

// Thermo4PFM
#include "CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB.h"

#include "SAMRAI/tbox/InputManager.h"

#include <boost/property_tree/ptree.hpp>

class CALPHADequilibriumPhaseConcentrationsThreePhasesStochioB
    : public EquilibriumPhaseConcentrationsThreePhases<
          Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB>
{
 public:
   CALPHADequilibriumPhaseConcentrationsThreePhasesStochioB(
       QuatModelParameters& model_parameters, const int conc_l_id,
       const int conc_a_id, const int conc_b_id, const int conc_l_ref_id,
       const int conc_a_ref_id, const int conc_b_ref_id,
       const Thermo4PFM::EnergyInterpolationType energy_interp_func_type,
       boost::property_tree::ptree calphad_pt,
       std::shared_ptr<tbox::Database> newton_db, const unsigned ncompositions);

   ~CALPHADequilibriumPhaseConcentrationsThreePhasesStochioB() {}

 private:
   QuatModelParameters d_model_parameters;

   /*
    * Implement specific method to compute auxilliary compositions when
    * phase is nearly 100% the sochiometric phase B
    */
   int computeAuxilliaryConcentrations(const double temp, double* c,
                                       double* hphi, double* x) override;
};

#endif
