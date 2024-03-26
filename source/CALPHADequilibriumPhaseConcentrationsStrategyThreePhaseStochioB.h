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
#ifndef included_CALPHADequilibriumPhaseConcentrationsStrategyThreePhaseStochioB
#define included_CALPHADequilibriumPhaseConcentrationsStrategyThreePhaseStochioB

#include "CALPHADequilibriumPhaseConcentrationsStrategy.h"

#include "CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB.h"

#include "SAMRAI/tbox/InputManager.h"

#include <boost/property_tree/ptree.hpp>

class CALPHADequilibriumPhaseConcentrationsStrategyThreePhaseStochioB
    : public CALPHADequilibriumPhaseConcentrationsStrategy<
          Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB>
{
 public:
   CALPHADequilibriumPhaseConcentrationsStrategyThreePhaseStochioB(
       const double cB, const int conc_l_id, const int conc_a_id,
       const int conc_b_id, const int conc_l_ref_id, const int conc_a_ref_id,
       const int conc_b_ref_id, boost::property_tree::ptree calphad_pt,
       std::shared_ptr<tbox::Database> newton_db);

   ~CALPHADequilibriumPhaseConcentrationsStrategyThreePhaseStochioB() {}
};

#endif
