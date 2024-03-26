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
#include "CALPHADequilibriumPhaseConcentrationsStrategyThreePhaseStochioB.h"
#include "Database2JSON.h"

namespace pt = boost::property_tree;

CALPHADequilibriumPhaseConcentrationsStrategyThreePhaseStochioB::
    CALPHADequilibriumPhaseConcentrationsStrategyThreePhaseStochioB(
        const double cB, const int conc_l_scratch_id,
        const int conc_a_scratch_id, const int conc_b_scratch_id,
        const int conc_l_ref_id, const int conc_a_ref_id,
        const int conc_b_ref_id, pt::ptree calphad_pt,
        std::shared_ptr<tbox::Database> newton_db)
    : CALPHADequilibriumPhaseConcentrationsStrategy<
          Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB>(
          conc_l_scratch_id, conc_a_scratch_id, conc_b_scratch_id,
          conc_l_ref_id, conc_a_ref_id, conc_b_ref_id,
          Thermo4PFM::EnergyInterpolationType::LINEAR,
          Thermo4PFM::ConcInterpolationType::LINEAR, false, calphad_pt,
          newton_db, 1)
{
   pt::ptree newton_pt;
   if (newton_db) copyDatabase(newton_db, newton_pt);

   d_calphad_fenergy = std::unique_ptr<
       Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB>(
       new Thermo4PFM::CALPHADFreeEnergyFunctionsBinaryThreePhaseStochioB(
           cB, calphad_pt, newton_pt,
           Thermo4PFM::EnergyInterpolationType::LINEAR,
           Thermo4PFM::ConcInterpolationType::LINEAR));
}
